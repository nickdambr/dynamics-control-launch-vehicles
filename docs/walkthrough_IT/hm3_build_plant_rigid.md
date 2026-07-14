# HM3/build_plant_rigid.m

## Ruolo del file nel progetto

Costruisce il **modello di corpo rigido** del lanciatore nel piano di beccheggio
al max-q: 4 stati `[z, zdot, theta, thetadot]`, 2 ingressi `[delta, alpha_w]`,
7 uscite. E' la **Eq. (1) della traccia privata delle due righe di bending**, ed
e' il plant su cui gira tutto Task 1 (progetto PD, carta di Nichols, risposta alla
raffica) e su cui `design_controller.m` ricava i guadagni iniziali.

Il file e' volutamente banale come codice -- tre matrici e un `ss` -- ma **e' il
punto in cui tutta la fisica del problema viene messa per iscritto**. All'orale e'
qui che si gioca la partita: ogni elemento delle matrici va saputo derivare.

Chi lo chiama: `main_task1.m` (riga 18), `main_task2.m` (riga 23, come termine di
paragone rigido), `main_montecarlo.m`, `make_pz_figures.m`, `init_simulink_hm3.m`,
`run_simulink_closed_loop.m`, e in ciclo su `t_ref`
`LTV_FULL_ASCENT/main_full_ascent.m`. Dipende solo dalla struct `p` di
`load_hw3_params`. Il suo output finisce in `assemble_loop(G, K, Wact)`, che chiude
il loop PD e restituisce `L` (per Nichols) e `T` (per le simulazioni).

Il risultato chiave, verificabile in una riga (`pole(build_plant_rigid(p))`):

    poli = { 0 ,  +0.0291 ,  +1.8165 ,  -1.8610 }   [rad/s]

cioe' **due poli a parte reale positiva**: il lanciatore rigido e' instabile in
open loop, e la retroazione non e' un'opzione ma una condizione di esistenza.

---

## Firma e contratto (righe 1-9)

```matlab
function G = build_plant_rigid(p)
%   G - ss, 4 states [z zdot theta thetadot], in [delta alpha_w],
%       out [theta_m thetadot_m z_m zdot_m theta z zdot]
```

- Riga 1: unico argomento, la struct `p`. **Nessun blocco `arguments`**: coerente
  con lo stile del repo, che tiene la validazione al confine di ingresso -- qui
  sta nel blocco `arguments` di `load_hw3_params`, che e' l'unica sorgente di `p`.
  (Gli unici usi in ciclo di questa funzione sono gli sweep su `t_ref` di
  `LTV_FULL_ASCENT/main_full_ascent.m` riga 43 e `main_q_scheduling.m` riga 56;
  `main_montecarlo.m` la chiama **una sola volta**, riga 41, fuori dal `parfor` --
  dentro al ciclo Monte Carlo gira `build_plant_full`.)
- Righe 6-9: il contratto delle uscite. Le **prime 4** sono le *misure* che
  alimentano il controllore (`theta_m, thetadot_m, z_m, zdot_m`); le **ultime 3**
  sono segnali di plotting (`theta, z, zdot` veri). Nel modello rigido misure e
  stati veri coincidono -- la distinzione esiste solo perche' in
  `build_plant_full.m` **non** coincideranno piu' (contaminazione INS del bending),
  e le due funzioni devono esporre la **stessa interfaccia** per essere
  intercambiabili in `assemble_loop`.

---

## La matrice A: le equazioni del moto (righe 11-14)

```matlab
A = [0   1        0               0;
     0   p.a1     p.a1*p.V+p.a4   0;
     0   0        0               1;
     0   p.A6/p.V p.A6            0];
```

Le due equazioni fisiche dietro queste quattro righe sono l'equazione dei
**momenti** (beccheggio) e quella delle **forze laterali** (drift), entrambe
linearizzate attorno alla traiettoria nominale e congelate a `t = 72 s`.

### Equazione dei momenti (riga 4 di A)

Il momento attorno al baricentro ha due contributi: quello aerodinamico
(`N_alpha * alpha` applicata a distanza `l_alpha` davanti al CG) e quello di
controllo (`Tc * delta` applicata a distanza `l_c` dietro il CG):

    Iyy*theta_ddot = N_alpha*l_alpha*alpha + Tc*l_c*delta

    theta_ddot = A6*alpha + K1*delta,     A6 = N_alpha*l_alpha/Iyy
                                          K1 = Tc*l_c/Iyy

L'angolo d'attacco **non e' uno stato**: e' una combinazione degli stati piu' il
disturbo. Con `zdot` la velocita' laterale del veicolo e `V` la velocita' assiale,

    alpha = theta + zdot/V - alpha_w

(sul segno di `alpha_w` vedi la sezione dedicata piu' avanti). Sostituendo:

    theta_ddot = A6*theta + (A6/V)*zdot - A6*alpha_w + K1*delta

che e' **esattamente** la riga 14 del codice (`[0, A6/V, A6, 0]`) piu' le colonne
di ingresso `K1` e `-A6`. Il termine `(A6/V)*zdot` e' il **coupling drift ->
assetto**: e' quello che rende il problema un vero sistema 4x4 e non due
sottosistemi separati, ed e' la ragione per cui i margini vanno letti sul loop
accoppiato e non sulla sola dinamica rotazionale.

- Riga 13 (`[0 0 0 1]`): semplicemente `d(theta)/dt = thetadot`.

### Equazione delle forze laterali (riga 2 di A)

    m*z_ddot = -N_alpha*alpha - (T - D)*theta + Tc*delta

Dividendo per `m` e usando le definizioni `a1 = -N_alpha/(m*V)`,
`a4 = -(T-D)/m`, `a3 = Tc/m`:

    z_ddot = a1*V*alpha + a4*theta + a3*delta
           = a1*V*(theta + zdot/V - alpha_w) + a4*theta + a3*delta
           = a1*zdot + (a1*V + a4)*theta + a3*delta - a1*V*alpha_w

che riproduce **carattere per carattere** la riga 12 del codice
(`[0, a1, a1*V+a4, 0]`) e le colonne `a3` e `-a1*V`. Il significato dei tre
coefficienti:

| coeff. | formula | segno | significato fisico |
|---|---|---|---|
| `a1` | `-N_alpha/(m*V)` | `-0.0154` (< 0) | smorzamento aerodinamico del drift: la forza normale si oppone alla velocita' laterale |
| `a3` | `Tc/m` | `+20.61` (> 0) | componente **laterale della spinta** quando l'ugello e' deflesso: e' il canale con cui il TVC muove il veicolo di lato |
| `a1*V + a4` | `-[(T-D) + N_alpha]/m` | `-41.73` (< 0) | accelerazione laterale per unita' di assetto: la spinta assiale **inclinata** di `theta` piu' la forza normale generata da `theta` |

> **Le formule sono definizioni, i numeri vengono dalla tabella.** Attenzione a
> non farsi prendere in castagna su `a4`: la definizione e' `a4 = -(T-D)/m`, ma
> con i valori di Tabella 1 (`T-D = 1.71e6 N`, `m = 7.38e4 kg`) darebbe `-23.17`,
> mentre la tabella -- e il set LPV, che concorda -- elenca `a4 = -27.2710`.
> L'incoerenza e' nota e documentata in `load_hw3_params.m` righe 36-38: nella
> dinamica entra `-27.2710`, ed e' quel valore (non il rapporto ricalcolato) a
> produrre tutti i numeri di questa pagina.

- Riga 11 (`[0 1 0 0]`): `d(z)/dt = zdot`.
- **La prima colonna di A e' identicamente nulla**: `z` non retroagisce su nulla.
  Il drift e' un puro integratore della velocita' laterale, e infatti uno dei
  quattro poli e' esattamente `0`. Questo ha una conseguenza pratica: la
  retroazione su `z_m` (guadagni `Kp_z = Kd_z = -1e-3`, fissati in
  `design_controller.m`) non puo' *stabilizzare* nulla -- serve a **contenere la
  deriva accumulata**. Attenzione a non chiamarla "load relief": disattivandola
  il picco di incidenza cambia di appena l'1% (`0.584` contro `0.577` deg).
  Compra *drift control*, non alleggerimento del carico.

> **Possibile domanda d'esame** -- Perche' `a3` (`= Tc/m`) e' positivo mentre `K1`
> (`= Tc*l_c/Iyy`) genera un momento che si oppone? Non e' contraddittorio?
> *Risposta:* no, e' la caratteristica **non-minimum-phase** del TVC. Deflettere
> l'ugello di `delta > 0` produce simultaneamente (i) una forza laterale
> istantanea sul veicolo nella direzione della componente di spinta (`a3*delta`) e
> (ii) un momento che ruota il veicolo nel verso opposto (`K1*delta`). Nel breve
> termine il veicolo si sposta "dalla parte sbagliata", nel lungo termine ruota e
> la spinta assiale lo riporta. E' l'analogo del "wrong-way" dell'elevatore
> aeronautico e si manifesta come uno zero a destra nella funzione di
> trasferimento `delta -> z`.

---

## Le colonne di ingresso (righe 16-17)

```matlab
Bd = [0; p.a3; 0; p.K1];          % delta column
Bw = [0; -p.a1*p.V; 0; -p.A6];    % alpha_w column
```

- Riga 16: il TVC entra in **due** righe: forza laterale `a3` e momento `K1`.
  Sono i due effetti simultanei discussi sopra.
- Riga 17: il vento entra con `-a1*V` e `-A6`, cioe' **con il segno opposto** a
  come `theta` entra tramite `alpha`. Ed e' l'unica informazione che serve per
  ricostruire la convenzione di segno adottata: se il disturbo entrasse come
  `alpha = ... + alpha_w`, la colonna sarebbe `[0; +a1*V; 0; +A6]`.

### Il segno di `alpha_w`: perche' il meno e' l'unico coerente con l'Eq. (1)

Oggi tutta la catena usa lo **stesso** segno, il meno:

- **il plant** (`build_plant_rigid.m` riga 17, e identicamente
  `build_plant_full.m` riga 28) implementa

      alpha = theta + zdot/V - alpha_w        ->   Bw = [0; -a1*V; 0; -A6]

- **il post-processing** (`simulate_gust_response.m` riga 29, usato da
  `main_task1.m` e da `main_task2.m` riga 188) ricostruisce

      r.alpha = r.theta + r.zdot/w.V - r.alphaw     % total angle of attack

**Perche' il meno non e' una convenzione ma un vincolo.** `alpha_w = v_w/V` con
`v_w` la velocita' laterale **dell'aria** (definizione di `load_wind_profile.m`,
riga 2). L'incidenza e' l'angolo fra l'asse del veicolo e il vento **relativo**,
quindi dipende dalla velocita' laterale del veicolo **rispetto all'aria**:

    alpha = theta + (zdot - v_w)/V = theta + zdot/V - alpha_w

Il plant **lo conferma da solo**, e basta guardare `Bw`: le due componenti non
nulle sono `-a1*V` e `-A6`, cioe' **esattamente** i coefficienti con cui `theta`
entra nella dinamica attraverso `alpha`, **cambiati di segno**. Se l'incidenza
fosse `theta + zdot/V + alpha_w`, la colonna del vento dovrebbe essere
`[0; +a1*V; 0; +A6]`. Il segno di `Bw` **e'** la definizione di `alpha`: non
sono due informazioni indipendenti. E' anche la convenzione degli appunti del
corso (Lez. 17, `-alpha_w`) e quella derivata nel report
(`Introduction.tex`, righe 61-64).

> **Nota storica -- il tipo di domanda che un esaminatore adora.** Fino a poco fa
> `simulate_gust_response.m` ricostruiva l'incidenza con il segno **piu'**,
> mentre il plant aveva gia' il meno. **Il modello non e' mai stato sbagliato**:
> il difetto era solo nella formula a valle, cioe' nella ricostruzione di `alpha`
> usata per l'indicatore di carico. Ma le conseguenze non erano cosmetiche:
>
> | formula | picco \|alpha\| (Task 1) | `qbar*alpha` |
> |---|---|---|
> | `+alpha_w` (vecchio, sbagliato) | 0.255 deg | 20.7 kPa*deg |
> | **`-alpha_w` (attuale, corretto)** | **0.577 deg** | **46.8 kPa*deg** |
>
> Piu' del doppio -- e con la **conclusione fisica rovesciata**. Con il `+`, il
> picco di incidenza cadeva *sotto* il contributo del solo vento (0.390 deg) e
> sembrava che l'anello facesse **load relief**. Con il segno giusto il picco lo
> **supera**: per tenere l'assetto il loop becca il muso **dentro** il vento
> relativo e il suo contributo si **somma** a quello del vento. Un puro
> attitude-hold e' **load-aggravating**. Coerentemente, con `A6 > 0` il centro di
> pressione sta davanti al baricentro: il momento aerodinamico e' divergente e
> **non c'e' nessuna stabilita' a banderuola** che allevi il carico da sola.

Cosa la correzione **non** ha toccato (ed e' importante saperlo dire):

- **Nulla** su stabilita', margini, Nichols, progetto del notch, corner di
  Task 3, Monte Carlo. Motivo: `alpha_w` e' un **ingresso di disturbo** e non
  entra nel loop di retroazione. `L(s)` e' calcolata da `assemble_loop` come
  `getLoopTransfer(T,'delta',-1)`, cioe' sul percorso `delta -> misure -> delta`:
  la colonna `Bw` **non compare**. Tutti i numeri di margine sono immuni.
- **Nulla** sui picchi di `theta`, `z` e `delta`: quelle sono uscite del plant,
  che era gia' giusto.
- **Solo** la quantita' ricostruita a valle -- il picco di incidenza totale e
  quindi `qbar*alpha` -- e con essa la lettura del trade attitude-hold /
  load-relief.

Da dire all'orale cosi' com'e': il modello e' sempre stato corretto, il segno lo
detta la colonna `Bw`, e la metrica di carico ora la rispetta.

---

## Le uscite: misure e segnali di plot (righe 19-28)

```matlab
C = [0 0 1 0;    % theta_m   = theta
     0 0 0 1;    % thetadot_m= thetadot
     1 0 0 0;    % z_m       = z
     0 1 0 0;    % zdot_m    = zdot
     0 0 1 0;    % theta
     1 0 0 0;    % z
     0 1 0 0];   % zdot
D = zeros(7, 2);
```

- Righe 21-27: le prime 4 righe (misure) sono **identiche** alle ultime 3 (verita')
  a meno della permutazione. E' esattamente cio' che il test
  `testRigidMeasurementsEqualTrueStates` (`hm3PlantTest.m` righe 65-70) verifica:
  `G.C(1,:) == G.C(5,:)` e `G.C(3,:) == G.C(6,:)`. La ridondanza e' **voluta**:
  serve a mantenere l'interfaccia identica a quella del plant flessibile, dove le
  due famiglie divergono.
- Riga 28: `D = 0`, il plant e' **strettamente proprio**: nessun percorso diretto
  `delta -> misure`. E' questo che permette a `getLoopTransfer` di rompere il loop
  su `delta` senza cappi algebrici.

---

## L'oggetto `ss` e i suoi nomi (righe 30-34)

```matlab
G = ss(A, [Bd Bw], C, D);
G.StateName  = {'z','zdot','theta','thetadot'};
G.InputName  = {'delta','alpha_w'};
G.OutputName = {'theta_m','thetadot_m','z_m','zdot_m','theta','z','zdot'};
```

- Righe 31-33: **i nomi non sono cosmetici, sono load-bearing.** Due usi concreti:
  (i) `assemble_loop.m` (riga 35) usa `connect(G, Kc, Wa, {'alpha_w','theta_ref'},
  {'theta','z','zdot','delta'}, {'delta'})`, che collega i blocchi **per nome**;
  (ii) `design_controller.m` (righe 41-45) ri-estrae i coefficienti dal plant per
  nome:

      iTh = strcmp(G.StateName,'theta');  iTd = strcmp(G.StateName,'thetadot');
      A6  = G.A(iTd, iTh);   K1 = G.B(iTd, iDe);

  cioe' legge `A6` e `K1` **dalla matrice**, non dalla struct `p`. Rinominare uno
  stato rompe silenziosamente il seed dei guadagni.

---

## Il polo instabile: da dove viene esattamente

La domanda "perche' il polo e' `+sqrt(A6)`" ha una risposta esatta e una piu'
onesta.

**Risposta esatta (dinamica rotazionale disaccoppiata).** Se si ignora il drift
(`zdot -> 0`), resta

    theta_ddot = A6*theta + K1*delta   ->   s^2 - A6 = 0   ->   s = +/- sqrt(A6)

con `A6 = 3.3818` -> `sqrt(A6) = 1.839 rad/s`. E' una **coppia** di poli reali
simmetrici, uno stabile e uno instabile: la firma della classica instabilita'
statica a sella (nessuna oscillazione, divergenza esponenziale con costante di
tempo `1/1.84 = 0.54 s`). E' il modello su cui D'Antuono ricava i guadagni
canonici `Kp = 2*A6/K1`, `Kd = sqrt(A6)/K1` usati come punto di partenza in
`design_controller.m` (righe 55).

**Risposta onesta (plant 4x4 completo).** La prima colonna di `A` e' nulla, quindi
un polo e' in `0` (l'integratore su `z`). I restanti tre autovalori sono quelli
del blocco `[zdot, theta, thetadot]`, il cui polinomio caratteristico e'

    s^3 - a1*s^2 - A6*s - A6*a4/V = 0

Sostituendo i valori a `t = 72 s`:

    s^3 + 0.01542*s^2 - 3.3818*s + 0.09835 = 0
    ->  s = -1.8610 ,  +1.8165 ,  +0.02908

- il **polo veloce instabile** `+1.8165` e' lo `sqrt(A6) = 1.839` spostato di
  ~1.2 % dal coupling del drift (e' proprio la tolleranza del 2 % del test
  `testRigidAirframeUnstablePole`, riga 61-62 di `hm3PlantTest.m`);
- il **polo lento instabile** `+0.0291` **non esiste nel modello rotazionale
  puro**: nasce dal termine noto `-A6*a4/V`, e vale approssimativamente

      s_lento ~= -a4/V = 27.271/937.71 = 0.02908 rad/s

  (coincide a 4 cifre). E' la **divergenza di drift**: con l'assetto tenuto fermo,
  la componente laterale della spinta accelera lentamente il veicolo di lato e il
  moto laterale a sua volta modifica l'angolo d'attacco. Costante di tempo ~34 s,
  quindi non e' un problema di stabilita' pratica in una finestra di volo di 10 s,
  ma **e' la ragione per cui il loop e' condizionalmente stabile** e per cui un
  singolo numero di `margin()` non ha senso: con due poli instabili la curva di
  Nichols deve tenere il punto critico (-180 deg, 0 dB) **incastrato** fra i suoi
  attraversamenti della fase critica, e servono margini classificati
  per banda di frequenza (`classify_margins.m`).

> **Possibile domanda d'esame** -- Quanti poli instabili ha il lanciatore rigido, e
> cosa implica per la lettura dei margini?
> *Risposta:* due (`+1.82` e `+0.029 rad/s`), piu' un integratore in `0`. Per il
> criterio di Nyquist con `P = 2` poli instabili in anello aperto, la stabilita'
> richiede **N = 2 encerchiamenti antiorari** del punto critico: il loop e'
> *conditionally stable*, cioe' diventa instabile sia aumentando **sia riducendo**
> il guadagno. Per questo il codice non si fida di `margin()` e classifica i
> margini per banda: il **gain margin aerodinamico** all'attraversamento della
> fase critica `-180 deg` a bassa frequenza (~0.6 rad/s, protegge dal *ridurre*
> il guadagno), i margini di fase/guadagno **rigidi** al crossover di controllo
> (~2.5 rad/s) e il delay margin. (Sulla carta di Nichols il punto critico sta a
> (-180 deg, 0 dB), convenzione degli appunti del corso da 1 + L = 0; D'Antuono
> etichetta lo stesso punto +180 -- fase mod 360 -- ed e' la rietichettatura che
> il codice usava fino a poco fa.)

---

## Possibili domande d'esame

**D: Scrivi le equazioni del moto del pitch-plane e mostra come diventano la
matrice A del codice.**
R: Momenti: `theta_ddot = A6*alpha + K1*delta` con `A6 = N_alpha*l_alpha/Iyy`,
`K1 = Tc*l_c/Iyy`. Forze: `z_ddot = a1*V*alpha + a4*theta + a3*delta` con
`a1 = -N_alpha/(m*V)`, `a4 = -(T-D)/m`, `a3 = Tc/m`. Il coupling e'
`alpha = theta + zdot/V - alpha_w`: sostituendolo, `theta_ddot` acquista il
termine `(A6/V)*zdot` (riga 4 di A) e `z_ddot` acquista `a1*zdot` e il termine
combinato `(a1*V + a4)*theta` (riga 2 di A). Le colonne di ingresso escono da
`delta` (`a3` sulla forza, `K1` sul momento) e da `alpha_w` (`-a1*V`, `-A6`).

**D: `A6` e' un momento aerodinamico. Perche' e' *destabilizzante* e non
*stabilizzante* come su un aeroplano?**
R: Perche' su un lanciatore snello il centro di pressione sta **davanti** al
baricentro (il corpo e' quasi un cilindro con un'ogiva, non ci sono superfici di
coda a spostare il cp all'indietro). Con `l_alpha > 0` misurato dal CG verso
prua, `A6 = N_alpha*l_alpha/Iyy > 0` e la dinamica di beccheggio isolata e'
`theta_ddot = +A6*theta`: un angolo d'attacco positivo genera un momento che lo
fa crescere. Su un aeroplano il cp e' dietro il CG, il segno si inverte e si
ottiene una coppia di poli immaginari (short period stabile).

**D: Perche' `z` compare come stato se la prima colonna di A e' nulla?**
R: Perche' e' la **variabile di uscita** che interessa (il drift laterale
accumulato, che entra nel budget di traiettoria) e perche' serve la sua misura
`z_m` per chiudere il debole loop di deriva. Dinamicamente e' un puro
integratore: non retroaziona su nulla, e infatti contribuisce l'autovalore `0`.
La retroazione su `z` con guadagni `-1e-3` non stabilizza (non puo'): riduce la
deriva accumulata. Non e' un dispositivo di load relief -- togliendola il picco
di incidenza si muove dell'1% -- ed e' bene dirlo prima che lo chiedano.

**D: Cosa succede se togli il termine `(A6/V)*zdot` dalla matrice A?**
R: Si ottiene il modello rotazionale disaccoppiato: i poli diventano
`{0, a1, +sqrt(A6), -sqrt(A6)} = {0, -0.0154, +1.839, -1.839}`. Il blocco
`[zdot, theta, thetadot]` si triangolarizza -- il drift non retroaziona piu' sul
beccheggio -- quindi `zdot` resta con il proprio autovalore `a1` e il beccheggio
con la coppia simmetrica `+/-sqrt(A6)`. Sparisce cosi' il polo lento **instabile**
a `+0.029` (al suo posto compare lo smorzamento aerodinamico **stabile** `a1`) e
con esso la parte a bassa frequenza del problema. Con i guadagni canonici il gain
margin aerodinamico "sulla carta" sarebbe ~6 dB; sul loop accoppiato reale scende
a ~4 dB, ed e' proprio la ragione per cui `design_controller.m` **ri-tuna** i
guadagni sul loop completo invece di fermarsi alla formula chiusa (lo dice
esplicitamente nel suo header, righe 6-10).

**D: Il modello e' valido solo a `t = 72 s`. Perche' progettare li'?**
R: Perche' e' il punto di **massima pressione dinamica** (`qbar = 81 kPa`), quindi
quello di massimo momento aerodinamico destabilizzante (`A6` scala con `qbar`) e
di massimo carico strutturale (`qbar*alpha`). E' il caso peggiore del volo
atmosferico: un controllore che chiude i margini richiesti li' e' tipicamente
adeguato altrove. La verifica di questa affermazione e' l'estensione LPV in
`LTV_FULL_ASCENT/`, che confronta il design frozen-time con uno gain-scheduled
sull'intera ascesa.
