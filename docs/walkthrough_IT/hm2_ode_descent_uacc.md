# HM2_powered_descent/ode_descent_uacc.m

## Ruolo del file nel progetto

Questo file e' **la stessa fisica di `ode_descent.m` riscritta con un'altra
parametrizzazione del controllo**: invece della spinta T = [Tx; Ty], il controllo
e' l'**accelerazione comandata** u = T/m. Non e' un dettaglio cosmetico: e' il
passo che apre la strada alla convessificazione. Con u come controllo, le righe
di velocita' perdono il fattore 1/m (che era la sorgente della bilinearita' fra
controllo e stato) e diventano un doppio integratore puro; e la riga della massa
diventa **omogenea in m**, il che significa che nella variabile z = ln m diventa
lineare.

E' il RHS "nativo" della **variante (d) di Task 2 -- GFOLD log-mass**
(Acikmese & Blackmore). La trascrizione GFOLD non lavora in coordinate (m, T) ma
in coordinate (z = ln m, u = T/m, sigma >= ||u||), dove la dinamica e'
**esattamente LTI** e viene discretizzata con un unico esponenziale di matrice in
`lti_zoh.m`. Questo file esiste per due ragioni precise:

1. **Definire il plant fisico corrispondente al ZOH sull'accelerazione.** La
   convenzione ZOH di GFOLD tiene costante u, non T: sull'intervallo la spinta
   fisica T(t) = m(t)*u fluttua seguendo la massa che cala. Se si replayasse una
   soluzione GFOLD con `ode_descent` (T-hold) si otterrebbe una traiettoria
   diversa e la validazione sarebbe scorretta.
2. **Fare da ground-truth non lineare** per la trascrizione LTI. Il replay
   `fwd_integrate_uacc` (`main_task2.m` righe 1285-1303) integra questo RHS con
   `ode45` a tolleranze strette *nelle coordinate fisiche originali*
   [x; y; vx; vy; **m**], e confronta il risultato con la predizione lineare.
   Il confronto avviene su due metriche distinte, e vale la pena tenerle
   separate: `node_err` (`main_task2.m` righe 859-868) misura **solo posizione e
   velocita'** (`X(:,1:4)` -- la docstring dice esplicitamente *"mass
   excluded"*), ed e' li' che il README riporta 7.3e-12 nondim, cioe' il
   pavimento dell'integratore; il canale di massa e' invece controllato dal
   drift di m_f calcolato da `land` (righe 142-144, che legge `X(end,5)`). La
   distinzione non e' pedanteria: le righe di posizione e velocita' non
   contengono sigma, quindi uno slack non aderente (sigma > ||u||) sarebbe
   **invisibile** all'errore di nodo e comparirebbe solo nel drift di massa.

Chi lo chiama: `main_task2.m` riga 1298 (dentro `fwd_integrate_uacc`, usato sia
per la validazione finale sia -- riga 1251 -- dentro il **ratio test** della SCvx
GFOLD, dove J_act e' misurato sul modello non lineare), `proto_gfold_logmass.m`
riga 249, e i test `tests/gfoldLogMassTest.m` (righe 20-52 e 92).

---

## `ode_descent_uacc` (righe 1-21)

```matlab
function dx = ode_descent_uacc(x, uacc, Vc)
...
umag = sqrt(uacc(1)^2 + uacc(2)^2);
dx = [ x(3); x(4); uacc(1); uacc(2) - 1; -Vc * x(5) * umag ];
```

- **Riga 1**: firma identica a `ode_descent` per forma (stato, controllo, Vc) ma
  con `uacc` al posto di `u`, a rimarcare che il secondo argomento e'
  un'**accelerazione**, non una forza. Lo stato resta
  `x = [x; y; vx; vy; m]` -- **in massa, non in log-massa**. Questa e' una scelta
  di design importante: il RHS non lavora in z. Il motivo e' che il replay deve
  produrre una traiettoria confrontabile con quelle delle altre tre varianti:
  `land` (`main_task2.m` righe 142-144) legge `X(end,5)` come **massa fisica**
  per il drift di m_f, e cosi' fanno il campo `sol.m` e i plot; `node_err`
  (righe 859-868) si limita invece alle colonne 1:4 (posizione e velocita').

- **Righe 2-16**: la docstring, che e' insolitamente ricca e contiene gia' la
  derivazione. Riga 11-13: *"With u held constant the thrust T = m(t)*u floats
  with the (depleting) mass, so vx_dot = ux and vy_dot = uy - 1 are exact, while
  the mass row reads m_dot = -Vc*||T|| = -Vc*m*||u||"*. Righe 13-14 mettono nero
  su bianco il contrasto con `ode_descent.m`. Riga 16 dichiara l'assenza del
  blocco `arguments` (hot loop dentro `ode45`).

- **Riga 18**: `umag = sqrt(uacc(1)^2 + uacc(2)^2)` -- il modulo
  dell'**accelerazione** comandata, non della spinta. Adimensionalmente, poiche'
  T_ref = m_ref*g e a_ref = g, la relazione u = T/m in variabili adimensionali
  resta letteralmente u' = T'/m' senza fattori di conversione: e' uno dei regali
  della scelta T_ref = m_ref*g_ref.

- **Riga 19**: le cinque righe della dinamica.

### Derivazione riga per riga

Partiamo dalla dinamica dimensionale (identica a quella di `ode_descent`):

    vx_dot = Tx/m ,   vy_dot = Ty/m - g ,   m_dot = -||T||/c

e sostituiamo la definizione del nuovo controllo, T = m*u:

- **vx_dot = ux**, **vy_dot = uy - 1** (adimensionale). La massa **sparisce
  completamente** dalle righe di velocita': erano bilineari in (T, 1/m), ora sono
  lineari in u e indipendenti dallo stato. Questo e' il guadagno strutturale
  principale. Il test `testAccelerationIsDirect` (`gfoldLogMassTest.m` riga 29)
  verifica proprio questo: con la stessa `uacc` e masse 0.9 e 0.3, dx(3:4) e'
  identico.

- **m_dot = -Vc*||T|| = -Vc*||m*u|| = -Vc*m*||u||** (perche' m > 0, quindi il
  modulo si porta fuori). E' l'unica riga che ancora accoppia stato e controllo,
  ed e' **bilineare** -- apparentemente siamo tornati da capo. Il colpo di scena
  e' che il termine e' *omogeneo di grado 1 in m*, quindi:

      d(ln m)/dt = m_dot / m = -Vc * ||u||

  cioe' **con z = ln m la riga della massa non dipende piu' dallo stato affatto**
  ed e' funzione del solo controllo. La non linearita' residua e' la norma ||u||,
  che si elimina introducendo lo slack sigma >= ||u|| (vincolo di **cono di
  secondo ordine**, convesso) e scrivendo z_dot = -Vc*sigma, ora **lineare**.
  Questo e' l'intero trucco di Acikmese-Blackmore: lo si vede istanziato in
  `lti_zoh.m` riga 24, dove `B(5,3) = -Vc` e' il coefficiente di sigma nella riga
  di z.

  Il test `testMassFlowScalesWithMass` (riga 39) verifica l'omogeneita': con
  ||u|| = 1, m_dot vale -Vc*1.0 a massa 1.0 e -Vc*0.5 a massa 0.5, cioe'
  d(ln m)/dt e' la stessa.

- **Terzo/quarto stato**: la gravita' resta un -1 costante, come in
  `ode_descent`, e finisce nel vettore di offset `c = [0;0;0;-1;0]` della forma
  LTI (`lti_zoh.m` riga 25). Non e' un controllo: e' un termine affine noto, ed
  e' per questo che la discretizzazione ZOH ha bisogno del terzo blocco `cbar`
  oltre ad `Abar` e `Bbar`.

### Il punto sottile: sigma e' davvero uguale a ||u||?

La dinamica LTI usa sigma, il RHS non lineare usa ||u||. **Coincidono solo se il
vincolo conico e' attivo all'ottimo.** Se il solver restituisse sigma > ||u||, la
predizione LTI consumerebbe **piu'** massa di quanta ne consumi il razzo vero, e
il replay con `ode_descent_uacc` darebbe una m finale piu' alta della predizione.

Perche' allora e' attivo? Perche' il costo e' `-XI(5,N)` (`main_task2.m` riga
1183), cioe' si **massimizza z_N**, e integrando z_dot = -Vc*sigma si ha
z_N = z_0 - Vc*sum(sigma_k*dt). Massimizzare z_N equivale a minimizzare la somma
degli sigma, quindi l'ottimizzatore li spinge sul loro lower bound, che e'
esattamente ||u_k||. Il cono si chiude **da solo**, non serve imporlo. E'
l'argomento di *lossless convexification*, ed e' esattamente cio' che il test
`testMassRowConsistency` (`gfoldLogMassTest.m` riga 85) verifica numericamente:
integrando questo RHS con ode45 e ponendo sigma = ||u||, la variazione di ln(m)
non lineare e la predizione LTI coincidono a 1e-9.

> **Possibile domanda d'esame** -- Se la dinamica in coordinate (z, u, sigma) e'
> esattamente LTI, perche' la variante (d) usa comunque un ciclo SCvx e non un
> singolo SOCP?
> *Risposta:* Perche' e' **la dinamica** a essere esatta, non tutti i vincoli. Il
> bound superiore di spinta ||T|| <= Tmax diventa, in coordinate GFOLD,
> sigma <= Tmax*exp(-z): un vincolo *non convesso* (esponenziale come maggiorante
> di una variabile). Il codice lo linearizza attorno a z_ref con lo sviluppo di
> Taylor `sigma <= Tmax*exp(-z_ref)*(1 - (z - z_ref))` (`main_task2.m` righe
> 1158-1159) e itera nella SCvx per far convergere z_ref. Il vero single-shot
> SOCP (lossless convexification completa, con la maggiorazione quadratica del
> lower bound e il bound superiore trattato una volta sola su un profilo di massa
> a priori) e' ancora aperto nel backlog come ticket T006.

---

## Perche' esistono entrambi i RHS

La domanda "perche' non un solo file?" ha una risposta precisa: **descrivono due
plant discreti diversi**, non due implementazioni della stessa cosa. Il ZOH e'
parte della definizione del problema discretizzato, non un dettaglio di
integrazione. Su un intervallo [t_k, t_{k+1}]:

- con `ode_descent` (T costante) l'accelerazione a(t) = T/m(t) **cresce** perche'
  m cala. La massa e' esattamente affine in t (m_dot = -Vc*||T|| e' costante).
- con `ode_descent_uacc` (u costante) l'accelerazione e' piatta e la **spinta
  fisica cala**, T(t) = m(t)*u. Il log-massa e' esattamente affine in t
  (z_dot = -Vc*||u|| e' costante), non la massa.

Sono due modelli di attuatore diversi. Non c'e' un "giusto" in assoluto: il
computer di bordo puo' benissimo comandare una spinta costante (T-hold) o
un'accelerazione costante compensando la massa (u-hold). Cio' che conta e' che la
trascrizione e la validazione usino **lo stesso** ZOH -- ed e' esattamente per
questo che il codice tiene due replay separati (`fwd_integrate` con `ode_descent`
per le varianti a/b/c, `fwd_integrate_uacc` con questo file per la variante d,
scelta esplicitamente alle righe 124-127 di `main_task2.m`).

Osservazione onesta: nel **limite dt -> 0** le due convenzioni convergono alla
stessa traiettoria continua, quindi la differenza e' un effetto di
discretizzazione, non di fisica. Con N = 50 e dt_nd ~ 0.044 la differenza fra i
due plant e' piccola ma non nulla, ed e' sufficiente a rovinare una validazione
fatta con il RHS sbagliato.

---

## Nessun blocco `arguments`: e' voluto

Riga 16: *"No arguments validation by design: hot-loop RHS inside ode45."* Qui il
carico e' meno estremo che per `ode_descent` (non c'e' `fmincon` con differenze
finite sopra), ma questo RHS viene comunque chiamato dentro `ode45` **a ogni
iterazione SCvx**, su N-1 = 49 intervalli, per calcolare J_act nel ratio test
(riga 1251). Con RelTol 1e-10 / AbsTol 1e-12 `ode45` fa decine di passi per
intervallo, ciascuno con 6 valutazioni del RHS. La validazione sta nel chiamante:
`fwd_integrate_uacc` **ha** il suo blocco `arguments` (righe 1289-1292).

---

## Possibili domande d'esame

**D: Qual e' esattamente la differenza fra `ode_descent` e `ode_descent_uacc`?**
R: La parametrizzazione del controllo. Il primo prende come controllo la
**spinta** T = [Tx; Ty] (forza, adimensionalizzata con T_ref = m0*g) e produce
vx_dot = Tx/m, m_dot = -Vc*||T||. Il secondo prende come controllo
l'**accelerazione comandata** u = T/m e produce vx_dot = ux (nessun 1/m),
m_dot = -Vc*m*||u||. Nel ZOH questo si traduce in due plant discreti diversi: nel
primo caso e' la spinta a essere piecewise-constant e l'accelerazione a crescere
mentre il razzo si alleggerisce; nel secondo e' l'accelerazione a essere piatta e
la spinta a calare con la massa.

**D: E' la versione "uacc" quella che rende la dinamica affine nel controllo?**
R: Quasi, ma va detto con precisione. Le righe di velocita' diventano
**esattamente affini** in u (vx_dot = ux, vy_dot = uy - 1: e' un doppio
integratore con offset costante), e questo elimina la bilinearita' T/m. Ma la
riga della massa, **scritta cosi' come e' nel file** (in variabile m), e' ancora
-Vc*m*||u||: bilineare in m e non lineare in u. La convessificazione richiede due
passi ulteriori, che *non* stanno in questo file ma in `lti_zoh.m` e nel SOCP:
(i) z = ln m, che rende la riga indipendente dallo stato: z_dot = -Vc*||u||;
(ii) lo slack sigma >= ||u||, che sostituisce la norma con una variabile lineare:
z_dot = -Vc*sigma, e sposta la non linearita' in un vincolo conico convesso.
Quindi: `ode_descent_uacc` e' il *ponte* verso la forma convessa, ma di per se'
e' ancora un modello non lineare -- ed e' voluto, perche' e' il ground truth con
cui si verifica la trascrizione lineare.

**D: Da dove viene esattamente `m_dot = -Vc*m*||u||`?**
R: Dal fatto che il consumo dipende dalla **spinta fisica**, non
dall'accelerazione. La legge del razzo e' m_dot = -||T||/c. Sostituendo la
definizione del controllo, T = m*u, e usando m > 0 per portare la massa fuori dal
modulo: m_dot = -||m*u||/c = -m*||u||/c. Adimensionalizzando (vedi la pagina di
`ode_descent.m`) il fattore 1/c diventa Vc = V_ref/c, da cui
m_dot = -Vc*m*||u||. Il fattore m e' cio' che dice "un razzo pesante deve
bruciare piu' propellente per ottenere la stessa accelerazione".

**D: Cosa succederebbe se si validasse la soluzione GFOLD con `fwd_integrate`
(cioe' con il RHS T-hold) invece che con `fwd_integrate_uacc`?**
R: Si misurerebbe un errore di nodo spurio, perche' si starebbe integrando un
plant diverso da quello che la trascrizione ha imposto. Il solver GFOLD ha
prodotto (ux, uy) piecewise-constant; ricostruendo T_k = m_k*u_k e tenendola
costante sull'intervallo, l'accelerazione crescerebbe (m cala) mentre la
trascrizione la assumeva piatta. L'errore sarebbe di ordine dt volte il tasso di
variazione della massa, cioe' O(Vc*||u||*dt) relativo: piccolo ma ordini di
grandezza sopra il 7.3e-12 riportato, e verrebbe scambiato per un difetto della
trascrizione LTI che invece e' esatta.

**D: La dinamica in z e' lineare. Vuol dire che il problema di minimo consumo e'
convesso?**
R: La *dinamica* si', ma il problema completo no, per due vincoli: il bound
superiore di spinta sigma <= Tmax*exp(-z) (non convesso, va linearizzato o
maggiorato) e -- in un caso realistico con Tmin > 0 -- il bound inferiore
sigma >= Tmin*exp(-z), che nella forma originale ||T|| >= Tmin e' un vincolo
"fuori da una palla", palesemente non convesso. In HM2 Tmin = 0
(`main_task2.m` riga 32), quindi il secondo problema non si pone e resta solo il
primo, che il codice affronta con SCvx. Il vincolo di glide-slope e il cono
sigma >= ||u|| sono invece gia' convessi.

**D: Perche' il file integra la massa m e non direttamente z = ln m, visto che in
z tutto e' lineare?**
R: Perche' il suo compito e' fare da **verifica indipendente**. Se integrasse z
con la stessa legge lineare usata dalla trascrizione, il confronto sarebbe
circolare e non proverebbe nulla. Integrando invece m con la legge non lineare
m_dot = -Vc*m*||u|| e confrontando exp(z_k) con m_k, si testa davvero che il
cambio di variabili e la ricostruzione T = m*u siano corretti. In piu', il layout
[x; y; vx; vy; m] e' quello atteso dal post-processing: `land` e i plot leggono
la quinta colonna come massa fisica, `node_err` legge le colonne 1:4 (posizione
e velocita', massa esclusa).
