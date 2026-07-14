# HM0_falcon9_ascent/main2.m

## Ruolo del file nel progetto

`main2.m` risolve **lo stesso identico problema fisico** di `main.m` -- ascesa del
primo stadio del Falcon 9, stesso veicolo, stesso sito, stesse tre fasi, stesse
forze -- ma con due cambiamenti strutturali che non toccano la fisica:

1. **Adimensionalizzazione completa** dello stato e dei parametri, con scale di
   riferimento `L_ref = R_E`, `V_ref = 7800` m/s (prima velocita' cosmica),
   `T_ref = L_ref/V_ref`, `m_ref = m0`.
2. **Riparametrizzazione del tempo**: al posto del tempo fisico si integra su una
   variabile `tau` in [0, 3], in cui **ciascuno dei tre archi di volo occupa
   esattamente un'unita'** di tau. La durata dimensionale di ciascun arco entra
   nelle equazioni come un fattore moltiplicativo `Delta_k` davanti al RHS.

Il file **non e' una versione migliore** di `main.m`, ne' un doppione. E' un
**file preparatorio**: entrambe queste scelte sono esattamente cio' che serve nei
problemi di ottimizzazione di traiettoria che seguono (HM1: shooting indiretto,
HM2: collocation diretta). L'adimensionalizzazione condiziona bene il problema
numerico (tutti i termini O(1), cosa che conta enormemente quando i costati o il
Jacobiano di un NLP mescolano metri, secondi e chilogrammi); la struttura a tre
archi con `Delta_k` in ingresso e' la forma **canonica dei problemi a durata di
arco libera**, dove i `Delta_k` diventano *variabili di decisione* e la griglia di
integrazione resta fissa.

Nella pipeline di HM0, `main2.m` ha un ruolo preciso: e' la **seconda
implementazione indipendente** contro cui `tests/falcon9AscentTest.m`
cross-valida `main.m`. Due codici scritti con scalature e parametrizzazioni
diverse che convergono agli stessi numeri sono una verifica molto piu' forte di
qualunque singolo test di regressione.

Differenze operative rispetto a `main.m`: `main2.m` **non esporta le figure in
PNG** (non ha il blocco di export), aggiunge una nona figura con lo stato
adimensionale contro tau, stampa a console le scale di riferimento e lo stato
finale adimensionale, e usa una event function (che, come si vedra' piu' sotto,
**non fa quello che i commenti dichiarano**).

---

## Intestazione: il contratto della riparametrizzazione (righe 1-35)

```matlab
%  Integration variable:
%    A single parameter tau in [0, 3] covers all three arcs.
%      tau in [0,1]  -> Arc 1: t* =  tau      * Delta1
%      tau in [1,2]  -> Arc 2: t* = (tau - 1) * Delta2 + t1*
%      tau in [2,3]  -> Arc 3: t* = (tau - 2) * Delta3 + t2*
%    The EOM is scaled accordingly:  dy*/dtau = Delta_k * (dy*/dt*)
```

- Righe 6-16: le scale di riferimento e quelle **derivate** (accelerazione, forza,
  pressione). Il punto chiave e' che si scelgono liberamente **quattro** scale
  indipendenti (L, V, m -- e T che ne discende), e tutte le altre sono *obbligate*
  per consistenza dimensionale:

      a_ref = V_ref^2 / L_ref          (= 9.539 m/s^2)
      F_ref = m_ref * a_ref            (= 5.43 MN)
      p_ref = F_ref / L_ref^2 = m_ref*V_ref^2 / L_ref^3

  Non c'e' liberta' di scelta qui: se si sbaglia una sola di queste derivate, il
  bilancio delle forze non torna e la traiettoria e' sbagliata di un fattore.
- Riga 32: **la chain rule**, il cuore matematico del file:

      dy*/dtau = (dy*/dt*) * (dt*/dtau) = Delta_k * (dy*/dt*)

  perche' entro l'arco k si ha `t* = t_start* + (tau - (k-1)) * Delta_k`, e quindi
  `dt*/dtau = Delta_k` costante a tratti. E' l'unica modifica alle equazioni: si
  moltiplica tutto il RHS per `Delta_k`.
- Riga 35: `clear`, con la stessa conseguenza sulla testabilita' vista in
  `main.m`.

**Perche' fare tutto questo?** La risposta onesta e' che **per HM0 non serve a
niente**: il problema e' ben condizionato anche in forma dimensionale, e la
riparametrizzazione (come si vedra') *peggiora* il RHS introducendo due
discontinuita' artificiali. Il guadagno arriva dopo:

- In un problema a **durate di arco libere** (tipico dell'ottimizzazione: quando
  finisce il pitchover? quando avviene lo staging?), i `Delta_k` sono incogniti.
  Con la parametrizzazione tau, l'intervallo di integrazione **resta [0, 3]
  qualunque siano i Delta_k**: la griglia non si muove, i nodi di collocation non
  si muovono, e i `Delta_k` entrano come semplici parametri moltiplicativi nel
  RHS -- quindi si possono differenziare rispetto a essi in modo banale. Se invece
  si integrasse in t con estremi variabili, ogni valutazione dell'obiettivo
  cambierebbe il dominio di integrazione: un incubo per un NLP.
- L'adimensionalizzazione con tutte le variabili O(1) e' cio' che rende
  utilizzabile un `fsolve`/`fmincon`: con r ~ 6.4e6, m ~ 5.7e5 e u ~ 1e2 nello
  stesso vettore, il Jacobiano ha numeri di condizionamento assurdi e le
  tolleranze relative/assolute perdono significato.

---

## `CONSTANTS AND VEHICLE PARAMETERS` (righe 37-83) e `LAUNCH SITE` (righe 85-90)

Identici a `main.m` righe 12-54 e 60-61: stessi valori numerici, stessi commenti.
Valgono quindi **le stesse osservazioni oneste** fatte in `hm0_main.md`:

- Riga 52: `Tsid = 86136` s, contro il giorno sidereo standard di 86164.09 s
  (errore relativo 3.2e-4, trascurabile ma reale).
- Riga 48: `Hscale = 8000` m non e' consistente con `Tamb = 288.15` K sotto
  l'ipotesi di gas ideale isotermo, che darebbe H = R*T/g = 8435 m.
- Riga 58: `a_sound` = 340.30 m/s **costante a tutte le quote** -- il Mach in
  quota e' sottostimato.
- Riga 67: `Qdot1 = Tvac1/cvac1` = 2536.07 kg/s; `Qdot1*tb1` = 410 843 kg contro
  `mp1` = 410 900 kg dichiarati.

Il fatto che le costanti siano **duplicate** anziche' condivise (per esempio in un
`params_falcon9.m` chiamato da entrambi) e' un debito tecnico reale: qualunque
correzione (per esempio `Tsid`) va fatta in due posti, e se se ne dimentica uno il
test di cross-validazione `main.m` vs `main2.m` fallirebbe segnalando una
"differenza fisica" che invece e' un errore di copia-incolla.

---

## `REFERENCE QUANTITIES` (righe 92-108)

```matlab
L_ref = RE;            % [m]
V_ref = 7800;          % [m/s]
T_ref = L_ref / V_ref; % [s]
m_ref = m0;            % [kg]
```

- Riga 96: `L_ref = RE`. Scelta con una conseguenza elegante: `r0* = 1`
  **esattamente** e la quota adimensionale e' semplicemente `r* - 1` (riga 507
  nella `eom_nd`). Nessuna sottrazione fra numeri grandi quasi uguali -- che in
  forma dimensionale e' un vero rischio di cancellazione catastrofica: al lift-off
  `h = r - RE` e' una differenza fra due numeri da 6.4e6 che vale 0. In doppia
  precisione non e' un problema serio, ma in singola lo sarebbe.
- Riga 97: `V_ref = 7800` m/s, la "prima velocita' cosmica" (velocita' orbitale
  circolare rasente). E' una scelta **fisicamente motivata**, non arbitraria: e' la
  velocita' che il veicolo dovra' *comunque* raggiungere. Cosi' la velocita'
  adimensionale si legge come "frazione di orbita gia' guadagnata": con la
  velocita' relativa a MECO di 2773.6 m/s (README), il primo stadio consegna circa
  il **36%** di V_ref.
- Riga 98: `T_ref = L_ref/V_ref` = **817.71 s**. Non e' libera: discende dalle
  altre due. Notare che il tempo di combustione, 162 s, corrisponde a
  `t*_b = 0.198`: tutta la missione HM0 avviene in un quinto del tempo di
  riferimento.
- Riga 101: `m_ref = m0`, quindi `m0* = 1` e la massa adimensionale e' la frazione
  di massa iniziale rimasta. A MECO vale ~0.278.

---

## `NON-DIMENSIONAL PARAMETERS` (righe 110-123)

```matlab
mu_nd     = mu     / (V_ref^2 * L_ref);
omegaE_nd = omegaE * T_ref;
rho0_nd   = rho0   * L_ref^3 / m_ref;
p0_nd     = p0     * L_ref^3 / (m_ref * V_ref^2);
H_nd      = Hscale / L_ref;
Tvac_nd   = Tvac1  * L_ref  / (m_ref * V_ref^2);
Aex_nd    = Aex1   / L_ref^2;
Qdot_nd   = Qdot1  * T_ref  / m_ref;
Sref_nd   = Sref   / L_ref^2;
```

Ogni riga e' una divisione per la scala di riferimento della grandezza
corrispondente. Vale la pena verificarle una per una, perche' un solo esponente
sbagliato produce una traiettoria sbagliata e silenziosamente plausibile.

- Riga 114: `[mu] = L^3/T^2`. La scala e' `L_ref^3/T_ref^2 = L_ref^3 *
  V_ref^2/L_ref^2 = L_ref * V_ref^2`. Quindi `mu_nd = mu/(V_ref^2 * L_ref)` =
  **1.0272**. E' un numero bellissimo: dice che `V_ref` e' *quasi esattamente* la
  velocita' circolare rasente, perche' quella vera sarebbe
  `sqrt(mu/RE)` = 7905 m/s (contro i 7800 assunti). Il fatto che `mu_nd ~ 1`
  conferma che la scala di velocita' e' stata scelta bene.
- Riga 115: `omegaE_nd = omegaE*T_ref` = **0.05965**. Una frequenza si
  adimensionalizza moltiplicando (non dividendo) per il tempo di riferimento.
- Riga 116: `[rho] = m/L^3` -> `rho_nd = rho * L_ref^3/m_ref`. Verifica: 1.225 *
  (6.378e6)^3 / 5.691e5 e' un numero enorme (~5.6e14). **Non tutti i parametri
  adimensionali sono O(1)** -- e non e' un errore: cio' che deve essere O(1) sono i
  *termini delle equazioni*, non i singoli parametri. Nel termine di drag
  `rho_nd * Sref_nd * CD * Vrel_nd^2`, il fattore `L_ref^3` di rho_nd si semplifica
  contro l'`1/L_ref^2` di Sref_nd, e il risultato e' O(1). Se all'orale si fa
  notare che rho_nd e' 5.6e14, la risposta e': guarda il *prodotto*, non il
  fattore.
- Riga 117: la pressione. `p_ref = F_ref/L_ref^2 = m_ref*V_ref^2/L_ref^3`, quindi
  `p0_nd = p0 * L_ref^3/(m_ref*V_ref^2)`. Corretto.
- Riga 119: `Tvac_nd = Tvac * L_ref/(m_ref*V_ref^2)` = **1.5155**. Una forza si
  divide per `F_ref = m_ref*V_ref^2/L_ref`. Il valore ~1.5 e' significativo:
  poiche' `mu_nd ~ 1` (cioe' l'accelerazione di gravita' adimensionale a r* = 1 e'
  ~1), **`Tvac_nd/m0_nd ~ 1.5` e' direttamente il rapporto spinta/peso in vuoto**.
  L'adimensionalizzazione, quando le scale sono ben scelte, fa emergere i numeri
  di merito del problema come valori dei parametri.
- Riga 121: `Qdot_nd = Qdot*T_ref/m_ref` = **3.644**. Significa che, se il motore
  potesse funzionare per un intero `T_ref`, brucerebbe 3.64 volte la massa
  iniziale del veicolo. Confrontato con `tb*` = 0.198, spiega perche' viene
  consumato il 72% della massa: 3.644 * 0.198 = 0.722.
- Riga 123: il commento "CD is already dimensionless" e' corretto e non banale: e'
  l'unico parametro che non richiede conversione.

> **Possibile domanda d'esame** -- L'adimensionalizzazione dovrebbe rendere tutti
> i termini O(1), ma `rho0_nd` vale circa 5.6e14. Non e' una contraddizione?
> *Risposta:* No, perche' cio' che deve essere O(1) sono i **termini delle
> equazioni**, non i parametri presi singolarmente. `rho0_nd` compare solo dentro
> il prodotto `rho_nd * Sref_nd * CD * Vrel_nd^2 / m_nd`, dove il fattore
> `L_ref^3` di `rho_nd` si cancella contro l'`1/L_ref^2` di `Sref_nd` lasciando un
> singolo `L_ref`, che a sua volta e' assorbito dalla scala di forza. Il risultato
> e' un'accelerazione adimensionale O(1). Il criterio corretto e' guardare il RHS
> assemblato, non i coefficienti.

---

## `NON-DIMENSIONAL INITIAL CONDITIONS` (righe 125-135)

```matlab
r0_nd = 1;
u0_nd = 0;
v0_nd = omegaE * RE * cos(lat0) / V_ref;
w0_nd = 0;
m0_nd = 1;
y0_nd = [r0_nd; lon0; lat0; u0_nd; v0_nd; w0_nd; m0_nd];
```

- Riga 129 e riga 133: `r0* = 1` e `m0* = 1` per costruzione, come voluto.
- Riga 131: `v0*` = 408.6/7800 = **0.0524**. La rotazione terrestre regala poco
  piu' del 5% della velocita' orbitale.
- Riga 135: **`theta` e `phi` restano dimensionali** (radianti). Corretto: gli
  angoli sono gia' adimensionali, non c'e' nulla da scalare. E' facile sbagliare
  qui e dividere anche loro per qualcosa. Da notare che questo rende il vettore di
  stato **eterogeneo**: 5 componenti adimensionalizzate e 2 angoli in radianti.
  Non e' un problema (i radianti sono O(1)), ma va detto.

---

## `ARC DURATIONS` (righe 137-147)

```matlab
Delta1 = t1_nd;              % arc 1 nd duration
Delta2 = t2_nd - t1_nd;      % arc 2 nd duration
Delta3 = tb_nd - t2_nd;      % arc 3 nd duration
```

Numericamente (calcolati con `T_ref = 817.71` s):

    Delta1 = 5   / 817.71 = 0.006115     (ascesa verticale, 5 s)
    Delta2 = 10  / 817.71 = 0.012229     (pitchover, 10 s)
    Delta3 = 147 / 817.71 = 0.179770     (gravity turn, 147 s)

- Notare `Delta3/Delta2 = 14.7` e `Delta2/Delta1 = 2`. Questi due rapporti sono
  **esattamente l'ampiezza dei salti che il RHS subisce** attraversando tau = 1 e
  tau = 2 (vedi sotto): l'arco 3 e' quasi 15 volte piu' lungo dell'arco 2, e il
  fattore Delta davanti al RHS salta di altrettanto.
- **Qui il codice mostra il suo limite come "file preparatorio"**: i `Delta_k` sono
  *calcolati* da tempi fissi (5, 15, 162 s), non sono variabili. La struttura e'
  pronta a riceverli come incognite (basterebbe metterli in `par` dall'esterno) ma
  in HM0 non lo sono. Il beneficio della parametrizzazione tau e' quindi
  **potenziale, non realizzato** in questo file.

---

## `PARAMETER STRUCTURE` (righe 149-169)

- Righe 153-162: i parametri adimensionali.
- Righe 168-169: qui c'e' un dettaglio che merita attenzione. La struct porta
  `par.T_ref` e `par.t_end1` **in unita' dimensionali**, in mezzo a parametri
  altrimenti tutti adimensionali. Il motivo e' alle righe 545-546: la legge di
  pitchover e' definita in **gradi al secondo dimensionali** (0.05 deg/s), quindi
  dentro l'EOM adimensionale bisogna **ri-dimensionalizzare il tempo** per
  valutarla.

  E' una piccola inconsistenza di design: la legge di guida non e' stata
  adimensionalizzata insieme al resto. La scelta alternativa "pulita" sarebbe stata
  esprimere il rateo di pitchover in `rad` per unita' di `t*` (cioe'
  0.05 * pi/180 * T_ref = 0.7135 rad per unita' di t*) e togliere `T_ref` dalla
  struct. Funzionalmente e' identico; ma la forma attuale garantisce che
  `main.m` e `main2.m` usino **esattamente la stessa legge di pitchover**, il che
  e' importante per il test di cross-validazione.

---

## `SINGLE ode45 CALL` (righe 171-184)

```matlab
opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12, ...
              'Events', @arc_boundary_events);

[tau_sol, Y_nd] = ode45(@(tau, y) eom_nd(tau, y, par), [0 3], y0_nd, opts);
```

- Riga 184: **una sola chiamata** a `ode45` copre tutta la missione, su `[0 3]`.
- Rispetto a `main.m`: stesse tolleranze (`RelTol 1e-10`, `AbsTol 1e-12`), ma
  **nessun `MaxStep`**. Non serve, perche' l'unita' di tau vale al massimo 147 s
  dimensionali: `ode45` sceglie passi in tau che si traducono automaticamente in
  passi dimensionali di grana fine.

### La event function NON fa quello che i commenti dicono

Il commento alla riga 178 afferma:

> `%  Arc boundaries are exact (Events function stops the step there).`

e `documentazione.txt` (sezione 10) rincara: *"il solver non si ferma ma inserisce
un punto esatto della soluzione in tau = 1 e tau = 2 [...] senza eventi, un
singolo step di Runge-Kutta potrebbe attraversare il cambio di scala producendo un
errore di integrazione."*

**Entrambe le affermazioni sono false, e questo e' stato verificato
sperimentalmente.** Con `isterminal = [0; 0]`, `ode45`:

- rileva l'evento e ne calcola l'istante esatto per interpolazione;
- lo restituisce nelle uscite **`te`, `ye`, `ie`** -- cioe' nella *quarta, quinta e
  sesta* uscita di `ode45`;
- **non modifica ne' la sequenza dei passi ne' la griglia di output `t`/`y`**;
- **non impedisce a un passo di attraversare il confine**.

Riproducendo lo schema esatto di `arc_boundary_events` su un problema di prova, la
griglia di output con e senza event function risulta **bit-per-bit identica**
(369 punti in entrambi i casi, `isequal(t, t_senza_eventi) = 1`), e **nessun punto
di output cade esattamente su tau = 1 o tau = 2** (il piu' vicino a tau = 1 e'
1.00241...).

In piu' -- colpo di grazia -- la riga 184 chiama `ode45` con **due sole uscite**:

```matlab
[tau_sol, Y_nd] = ode45(...);
```

quindi `te`, `ye`, `ie` **non vengono nemmeno raccolte**. La event function e'
quindi completamente **inerte**: rimuoverla non cambierebbe di un bit il
risultato. E' codice morto documentato come se fosse essenziale.

**Cosa protegge davvero la soluzione ai confini degli archi**, allora? Solo il
controllo di passo adattivo di `ode45`: attraversando tau = 1 e tau = 2 il RHS
salta (di un fattore 2 e di un fattore 14.7 rispettivamente), lo stimatore
d'errore locale esplode, il passo viene rifiutato e ridotto ripetutamente finche'
non rientra in `RelTol = 1e-10`. Funziona -- i risultati sono corretti, come
conferma il test di cross-validazione contro `main.m` -- ma e' costoso, non e'
"esatto", e non e' merito della event function.

**La correzione corretta** sarebbe `isterminal = [1; 1]` con **tre chiamate
separate** a `ode45` (una per arco), riavviando ogni volta dallo stato finale del
precedente. Cosi' nessun passo attraverserebbe mai una discontinuita'. E' esatta-
mente la struttura *multi-arco* che serve comunque negli homework di ottimizzazione.

> **Possibile domanda d'esame** -- La parametrizzazione in tau e' presentata come
> un miglioramento numerico. Lo e' davvero, per HM0?
> *Risposta:* No, per HM0 e' un **peggioramento**. In `main.m` il RHS e' continuo a
> t = 5 s (la legge di pitchover a t = 5 restituisce esattamente 90 deg,
> riproducendo la fase 1) e ha una piccola discontinuita' di direzione a t = 15 s.
> In `main2.m`, invece, il fattore `Delta_k` moltiplica l'**intero** RHS e cambia a
> gradino: a tau = 1 il RHS **raddoppia** di colpo (Delta2/Delta1 = 2) e a tau = 2
> viene **moltiplicato per 14.7** (Delta3/Delta2 = 14.7). La riparametrizzazione
> *introduce* due discontinuita' grandi che nella formulazione dimensionale non
> esistevano, costringendo `ode45` a rifiutare passi ai confini. Il guadagno
> arriva solo quando i `Delta_k` diventano **variabili di decisione** in un
> problema di ottimizzazione a durate di arco libere: li' il dominio di
> integrazione fisso [0, 3] vale ampiamente il prezzo.

---

## `RECOVER GLOBAL ND TIME AND DIMENSIONAL QUANTITIES` (righe 186-268)

### Ricostruzione del tempo (righe 190-200)

```matlab
mask1 =              tau_sol <= 1;
mask2 = tau_sol > 1 & tau_sol <= 2;
mask3 = tau_sol > 2;
t_nd(mask1) =  tau_sol(mask1)       * Delta1;
t_nd(mask2) = (tau_sol(mask2) - 1)  * Delta2  + t1_nd;
t_nd(mask3) = (tau_sol(mask3) - 2)  * Delta3  + t2_nd;
t = t_nd * T_ref;
```

- Righe 195-197: e' l'inversa **esatta** della mappa tau -> t*. Il tempo
  dimensionale `t(tau)` risulta **lineare a tratti con tre pendenze diverse**,
  proporzionali a Delta1, Delta2, Delta3. La Figura 9, pannello (3,2), lo mostra
  esplicitamente ed e' il modo piu' rapido per verificare che la mappa sia giusta:
  se le tre pendenze non stanno nel rapporto 5 : 10 : 147, c'e' un errore.
- **Nota sulla coerenza**: la stessa mappa e' implementata **due volte** -- qui in
  post-processing (righe 195-197) e dentro `eom_nd` (righe 480-496). Sono
  duplicate a mano. Se una delle due venisse modificata senza l'altra, il codice
  produrrebbe silenziosamente risultati sbagliati (la dinamica userebbe un tempo,
  il post-processing un altro). E' un rischio reale di manutenzione.

### Ri-dimensionalizzazione e post-processing (righe 202-268)

- Righe 203-209: si moltiplica per le scale. `theta` e `phi` (righe 204-205)
  **non** vengono scalati, coerentemente con il fatto che erano gia' radianti.
- Righe 212-268: da qui in poi il post-processing e' **riga per riga identico a
  quello di `main.m`** (righe 114-188): quota, Vrel, qdyn, Mach, gammaT/gammaA,
  eventi, ENU, ground track. Valgono le stesse osservazioni:
  - `qdyn` con la velocita' **relativa** (riga 224), corretto;
  - `gammaT(k) = gammaA(k)` in fase 3 (riga 242) e' **assegnato**, quindi la
    coincidenza delle due curve nella Figura 7 e' una tautologia, non una verifica;
  - `im1` e `imQ` (righe 249-250) sono estratti **sui campioni**, senza
    interpolazione. Qui pero' **non c'e' `MaxStep`**: la griglia e' scelta
    liberamente da `ode45`. E' plausibile che la griglia in tau sia piu' rada, in
    termini dimensionali, di quella di `main.m` -- motivo per cui il test tollera
    `RelTol = 1e-2` sul confronto di `qmax` fra i due script (vedi
    `hm0_test_falcon9Ascent.md`).
- **Riga 223**: `p_traj = p0 * exp(-h/Hscale)` viene calcolata e **mai usata** nel
  seguito. E' una variabile morta (in `main.m` non c'e' affatto). Innocua, ma e'
  rumore.

---

## `CONSOLE SUMMARY` (righe 270-300)

Identico a `main.m` (righe 194-212) fino alla riga 292, poi aggiunge (righe
293-299) lo **stato finale adimensionale**: `r*`, `u*`, `v*`, `w*`, `m*`, `t*`.
E' l'output piu' utile per la verifica: `m*(end)` deve valere
`1 - Qdot_nd * tb_nd` = 1 - 3.644*0.198 = 0.278 e `r*(end)` deve valere
`1 + h_end/RE ~ 1.013`. Sono controlli che si possono fare a mente.

---

## `PLOTS` (righe 302-462)

- Figure 1-8 (righe 311-415): **le stesse otto figure di `main.m`**, con gli stessi
  eventi marcati e la stessa palette, ma con `[ND]` nel `Name` e "(ND simulation)"
  nel titolo. Servono al confronto visivo fra i due script.
- **`main2.m` non ha il blocco di export PNG.** Le figure restano a schermo e non
  vengono scritte su disco. E' una scelta deliberata e corretta: `figures/` e'
  di proprieta' della pipeline di `main.m` (come dice il README), e due script che
  scrivono gli stessi nomi file si sovrascriverebbero a vicenda in modo
  imprevedibile. Il commento della classe di test (`falcon9AscentTest.m`, righe
  7-8) dice "running the scripts regenerates the PNGs": e' vero solo per `main.m`.
- **Figura 9** (righe 417-462): l'unica figura che `main.m` non ha. Un
  `tiledlayout(3,2)` con lo **stato adimensionale contro tau**: `r*`, le componenti
  `u*,v*,w*`, `m*`, `h* = r*-1`, `|V*|`, e `t_dim(tau)`.
  I pannelli hanno `xline(1)` e `xline(2)` a marcare i confini d'arco. E' il
  pannello **diagnostico** della riparametrizzazione: guardando (3,2) si legge
  immediatamente la mappa lineare a tratti tau -> t; guardando (1,2) si vede che
  in tau i tre archi occupano lo stesso spazio orizzontale pur avendo durate
  fisiche di 5, 10 e 147 s -- e' proprio l'effetto che si voleva.

---

## `eom_nd` -- equazioni del moto adimensionali (righe 468-578)

### Identificazione dell'arco e ricostruzione del tempo (righe 479-496)

```matlab
if tau <= 1
    Delta = par.Delta1;  t_arc_start = 0;         tau_k = tau;      phase = 1;
elseif tau <= 2
    Delta = par.Delta2;  t_arc_start = par.t1_nd; tau_k = tau - 1;  phase = 2;
else
    Delta = par.Delta3;  t_arc_start = par.t2_nd; tau_k = tau - 2;  phase = 3;
end
t_nd = t_arc_start + tau_k * Delta;
```

- E' la stessa cascata di `if` di `main.m` (righe 425-450), ma con una funzione in
  piu': oltre a scegliere la fase, deve **ricostruire il tempo adimensionale
  globale `t*`** dalla variabile di integrazione `tau`, perche' la legge di
  pitchover ne ha bisogno.
- **Confini**: i test sono `tau <= 1` e `tau <= 2`, quindi il valore esatto
  `tau = 1` e' assegnato all'**arco 1** e `tau = 2` all'**arco 2**. E' arbitrario ma
  innocuo: nei punti di confine le due definizioni coincidono comunque (a tau = 1,
  la legge di pitchover con t_dim = 5 s restituisce 90 deg, che e' esattamente la
  fase 1). A tau = 2, invece, le due fasi **non** coincidono, quindi il valore del
  RHS esattamente in tau = 2 e' quello della fase 2. Poiche' e' un insieme di
  misura nulla, non ha effetto sull'integrazione.

### Corpo delle equazioni (righe 498-574)

Da riga 498 a riga 574, la struttura e' **identica** a quella di `main.m` (righe
387-463). Questa e' la proprieta' piu' importante e vale la pena enunciarla
esplicitamente:

> **Le equazioni adimensionali hanno la stessa forma algebrica di quelle
> dimensionali.**

Non e' un caso: e' la conseguenza di aver scelto le scale in modo
**auto-consistente** (a_ref = V_ref^2/L_ref, F_ref = m_ref*a_ref, ecc.).
Concretamente:

- riga 507: `alt = r - 1` invece di `alt = r - RE` (perche' `R_E* = 1`);
- riga 523: `g_u = -par.mu/r^2` -- stessa formula, `mu` adimensionale;
- riga 528: stesso coefficiente di drag vettoriale;
- riga 535: `Tmag = par.Tvac - patm*par.Aex` -- **la correzione di contropressione
  ha la stessa identica forma**, con tutti i simboli adimensionali;
- righe 564-571: cinematica e dinamica **letteralmente identiche** (stessi termini
  di trasporto `(v^2+w^2)/r`, `-u*v/r`, `v*w*tan(phi)/r`, ...).

Se una qualunque delle scale derivate fosse sbagliata, **questa proprieta' si
romperebbe** e comparirebbero fattori spuri nelle equazioni. Il fatto che il corpo
di `eom_nd` sia un copia-incolla di `eom` e' quindi la miglior prova che le
adimensionalizzazioni delle righe 114-122 sono corrette.

La derivazione dei termini di trasporto (perche' `(v^2+w^2)/r`, perche' non ci
sono termini di Coriolis, ecc.) e' identica a quella dimensionale ed e' svolta per
esteso in **`hm0_main.md`**, sezione `eom`.

### La legge di pitchover (righe 543-550)

```matlab
case 2
    t_dim = t_nd * par.T_ref;
    gT    = deg2rad(90 - 0.05 * (t_dim - par.t_end1));
    Tu    = Tmag * sin(gT);
    Tv    = Tmag * cos(gT);
    Tw    = 0;
```

- Riga 545: **si ri-dimensionalizza il tempo** per valutare la legge di guida. E'
  la conseguenza di quella scelta di design discussa sopra: la legge (0.05 deg/s)
  e' definita in unita' dimensionali. Funziona, ma significa che `eom_nd` non e'
  "puramente adimensionale": porta dentro `T_ref`.
- Righe 548-550: qui c'e' la **differenza numerica reale** rispetto a `main.m`.
  `main.m` (righe 437-439) scrive:

      Tv = Tmag * cos(gT) * sin(psi)     con psi = deg2rad(90)
      Tw = Tmag * cos(gT) * cos(psi)

  e poiche' `cos(pi/2)` in floating point vale 6.1e-17 e non 0, `main.m` ha una
  componente Nord di spinta **non esattamente nulla** (dell'ordine di 1e-12 N:
  Tmag ~ 7e6 N per cos(gT) <= 8.7e-3 per cos(psi) = 6.1e-17).
  `main2.m` scrive direttamente `Tv = Tmag*cos(gT)` e `Tw = 0`, sostituendo a mano
  `sin(90 deg) = 1` e `cos(90 deg) = 0`. E' piu' pulito e piu' veloce, e spiega
  perche' i due script non possono coincidere bit-per-bit. L'effetto e' comunque
  del tutto trascurabile rispetto alla differenza dovuta alle diverse griglie
  temporali.
- Nota che `main2.m` **perde la generalita'**: se domani si volesse un azimuth
  diverso da 90 deg, `main.m` basterebbe cambiare `psi`, mentre `main2.m` andrebbe
  riscritto. E' un trade-off consapevole (o forse no -- il codice non lo dichiara).

### Chain rule (righe 576-577)

```matlab
% Chain rule: dy*/dtau = Delta_k * (dy*/dt*)
dydt = Delta * [dr; dtheta; dphi; du; dv; dw; dm];
```

- Riga 577: **l'unica riga che distingue davvero questo EOM da quello
  dimensionale.** Tutto il RHS, **tutte e sette le componenti massa inclusa**,
  viene moltiplicato per `Delta`.
- E' cruciale che sia *tutto* il vettore: se si dimenticasse di scalare `dm`, la
  massa evolverebbe con un tempo diverso da quello della dinamica e il veicolo
  arriverebbe a MECO con la massa sbagliata. Il test `testPropellantBookkeeping`
  (che verifica `m(t_b) = m0 - Qdot*t_b` a `RelTol 1e-8` **anche per `main2.m`**) e'
  precisamente il controllo che smaschererebbe un errore in questa riga: e' un
  oracolo analitico esatto sulla componente piu' semplice del sistema, e passa
  solo se la mappa tau -> t e la scalatura `Delta_k` sono entrambe giuste.

---

## `arc_boundary_events` (righe 582-589)

```matlab
function [value, isterminal, direction] = arc_boundary_events(tau, ~, ~)
    value      = [tau - 1;   tau - 2];
    isterminal = [0;          0];
    direction  = [+1;         +1];
end
```

- Riga 586: due funzioni evento che si annullano a tau = 1 e tau = 2.
- Riga 587: `isterminal = [0; 0]` -- **il solver non si ferma**.
- Riga 588: `direction = [+1; +1]` -- evento rilevato solo attraversando in
  direzione crescente. Corretto (tau cresce sempre), ma irrilevante.
- **Come documentato sopra, questa funzione e' inerte**: con `isterminal = 0`
  `ode45` non tocca ne' il passo ne' la griglia di output, e le uscite `te`/`ye`
  non vengono nemmeno richieste alla riga 184. Il commento di intestazione (righe
  583-585, "so ode45 never straddles an arc boundary") **descrive un comportamento
  che non si verifica**.
- La firma ignora sia lo stato sia i parametri (`~, ~`): l'evento dipende solo
  dalla variabile di integrazione. Questo di per se' e' un indizio: un evento che
  dipende solo da `tau` e' un evento a *istante noto*, e gli eventi a istante noto
  non hanno bisogno di essere "trovati" -- vanno gestiti spezzando l'integrazione.

---

## Possibili domande d'esame

**D: Qual e' la differenza sostanziale fra `main.m` e `main2.m`?**
R: **Nessuna, sul piano fisico.** Stesso veicolo, stesse forze (gravita' 1/r^2,
drag con CD costante, spinta corretta per contropressione), stesse tre fasi, stesse
equazioni in coordinate sferiche con velocita' UEN, stessi valori numerici di tutte
le costanti. Le differenze sono **due, entrambe di formulazione**: (1)
adimensionalizzazione con L_ref = R_E, V_ref = 7800 m/s, m_ref = m0; (2)
riparametrizzazione del tempo su tau in [0,3] con un arco per unita' di tau e il
fattore `Delta_k` davanti al RHS. Differenze operative minori: `main2.m` non
esporta i PNG, aggiunge una nona figura diagnostica sullo stato adimensionale, e
scrive `Tw = 0` esattamente in fase 2 dove `main.m` calcola `cos(90 deg)` in
floating point (6.1e-17). Il test di cross-validazione conferma che i due
convergono agli stessi risultati entro lo 0.5%.

**D: A cosa serve davvero la parametrizzazione in tau, visto che in HM0 le durate
degli archi sono fisse?**
R: Non serve a niente in HM0 -- anzi, peggiora il condizionamento del RHS
introducendo due discontinuita' (fattori 2 e 14.7) ai confini d'arco. Serve **per
gli homework di ottimizzazione**. Quando le durate degli archi diventano variabili
di decisione (durata del pitchover, istante di staging, tempo finale libero), la
parametrizzazione tau garantisce che il **dominio di integrazione resti [0, 3]
qualunque siano i Delta_k**: la griglia di collocation non si muove, i nodi non si
muovono, e i `Delta_k` entrano nel problema come semplici parametri moltiplicativi
nel RHS -- quindi con derivate banali rispetto alle variabili di decisione.
Integrando invece nel tempo fisico con estremi mobili, ogni valutazione
dell'obiettivo cambierebbe il dominio, cosa che un NLP gestisce malissimo. E' la
forma canonica dei problemi multi-arco a durata libera.

**D: La event function `arc_boundary_events` serve a evitare che `ode45`
attraversi le discontinuita' ai confini d'arco. Vero o falso?**
R: **Falso**, e i commenti nel codice (riga 178) e in `documentazione.txt` sono
sbagliati su questo punto. Con `isterminal = [0; 0]`, `ode45` rileva l'evento e
lo riporta nelle uscite `te`/`ye`/`ie`, ma **non modifica ne' la sequenza dei
passi ne' la griglia di output**: verificato sperimentalmente, la griglia con e
senza event function e' identica e nessun punto di output cade su tau = 1 o
tau = 2. In piu' la riga 184 chiama `ode45` con sole **due** uscite, quindi
`te`/`ye`/`ie` sono scartate: la funzione e' completamente inerte. Cio' che
protegge la soluzione ai confini e' semplicemente il controllo di passo adattivo,
che rifiuta e riduce i passi finche' l'errore locale rientra in tolleranza. Per
ottenere davvero l'effetto dichiarato servirebbe `isterminal = [1; 1]` e **tre
chiamate `ode45` separate** con restart a ogni confine.

**D: `mu_nd` vale 1.027 e `Tvac_nd` vale 1.516. Sono numeri "casuali" o dicono
qualcosa?**
R: Dicono molto, ed e' il segno che le scale sono state scelte bene. `mu_nd ~ 1`
significa che `V_ref = 7800` m/s e' quasi esattamente la velocita' orbitale
circolare rasente (quella vera e' sqrt(mu/R_E) = 7905 m/s): l'accelerazione di
gravita' adimensionale alla superficie vale quindi ~1. Di conseguenza
`Tvac_nd / m0_nd = 1.516 / 1 = 1.516` **e' direttamente il rapporto spinta/peso in
vuoto** del veicolo. Un'adimensionalizzazione ben fatta trasforma i parametri nei
numeri di merito del problema: si legge il regime di volo direttamente dai valori
dei coefficienti, senza fare conti.

**D: Perche' `theta` e `phi` non vengono adimensionalizzati?**
R: Perche' sono **gia'** adimensionali: un angolo in radianti e' un rapporto fra
lunghezze (arco / raggio). Non esiste una "scala di riferimento per gli angoli" da
applicare. Il vettore di stato risulta quindi eterogeneo -- 5 componenti scalate e 2
angoli in radianti -- ma questo non e' un problema, perche' i radianti sono gia'
O(1) e quindi l'obiettivo dell'adimensionalizzazione (avere tutte le componenti di
grandezza confrontabile) e' comunque soddisfatto.

**D: Perche' `eom_nd` porta dentro `par.T_ref`, che e' un'unita' dimensionale?**
R: Perche' la legge di pitchover e' definita in gradi al secondo **dimensionali**
(0.05 deg/s), e la riga 546 deve valutarla: la riga 545 ricostruisce percio'
`t_dim = t_nd * T_ref`. E' una piccola inconsistenza di design: la legge di guida
non e' stata adimensionalizzata insieme al resto del modello. Si poteva evitare
esprimendo il rateo come 0.05 * (pi/180) * T_ref = 0.7135 rad per unita' di `t*` e
togliendo `T_ref` dalla struct. Il vantaggio della forma attuale e' che garantisce
che `main.m` e `main2.m` usino **la stessa identica legge**, il che e' essenziale
perche' il test di cross-validazione abbia senso.
