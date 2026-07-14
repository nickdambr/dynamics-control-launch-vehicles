# HM3/main_task1.m

## Ruolo del file nel progetto

E' il **primo dei tre entry point** di HM3, quello che il professore fa girare per
il Task 1 della traccia ("Rigid LV"). Lo script prende il lanciatore fittizio di
Greensite congelato all'istante di **max-q_bar** (t = 72 s), lo modella come
**corpo rigido con attuatore ideale** (niente dinamica TVC, niente ritardo,
niente bending), progetta un controllore **PD di assetto** con una debole
retroazione negativa di deriva laterale, ne verifica la stabilita' sul **piano di
Nichols** e la valida con una **risposta temporale a raffica di vento**. Sono
esattamente i due deliverable che la traccia chiede: "verify closed-loop
stability using classical control tools in the frequency domain (e.g., the
Nichols plot)" e "simulate the system's time-domain response to a wind gust.
Plot the time histories of theta, z, zdot, delta".

Lo script e' un **orchestratore**: non contiene matematica propria se non due
righe di cross-check analitico (righe 25-26). Tutto il lavoro sta nei moduli, che
chiama in cascata:

    load_hw3_params  ->  build_plant_rigid  ->  design_controller
                                                    |
                                       (assemble_loop + classify_margins)
                                                    |
                        load_wind_profile  ->  simulate_gust_response
                                                    |
                                             plot_nichols_lv

Il punto concettuale che governa tutto il file e' che il velivolo a max-q_bar e'
**aerodinamicamente instabile**: il centro di pressione sta davanti al centro di
massa, quindi il momento aerodinamico e' *divergente* (A_6 > 0). Il polo di corpo
rigido sta a +sqrt(A_6) = +1.84 rad/s. La retroazione non e' un miglioramento
opzionale, e' l'unica cosa che tiene insieme il lanciatore. La conseguenza
tecnica, che ricorre in tutti e tre i task, e' che l'anello e' **condizionalmente
stabile**: si destabilizza sia *alzando* troppo il guadagno sia *abbassandolo*
troppo. Per questo lo script non legge mai i margini con un `margin()` secco, ma
li fa **classificare per banda di frequenza** da `classify_margins`.

I gain di deriva `Kp_z = Kd_z = -1e-3` non vengono toccati dallo script: sono i
default di `design_controller`, e vengono dalla guideline della traccia ("The
values of KP,z and KD,z must be negative in the order of 10^-3").

---

## `%% Model and parameters (Table 1 @ t = 72 s)` (righe 12-18)

```matlab
p = load_hw3_params();
fprintf('Parameters source: %s\n', p.src);
G = build_plant_rigid(p);
```

- Riga 13: `load_hw3_params()` senza argomenti -> scaling di incertezza a 1
  (nominale). La funzione prova a leggere `General/hw3-v3/GreensiteLPV_DATA.mat`
  e a interpolare i coefficienti tempo-varianti a t = 72 s; se il file manca,
  ricade sui letterali di Table 1. Il campo `p.src` dice quale delle due strade
  e' stata presa, ed e' stampato apposta a riga 14: e' l'unico modo, a run-time,
  per sapere se i numeri vengono dal dataset o dai letterali.
- Righe 15-16: stampa i tre coefficienti che contano e, soprattutto,
  `sqrt(p.A6)`, etichettato come "unstable airframe pole". Nella run reale:

      A6 = 3.3818 1/s^2,  K1 = 4.5647 1/s^2,  V = 937.7 m/s
      polo instabile a +sqrt(A6) = +1.839 rad/s

  Da dove viene: nella dinamica rotazionale disaccoppiata (trascurando lo
  zdot-coupling) si ha theta_ddot = A_6*theta, cioe' s^2 - A_6 = 0, quindi i due
  poli +/- sqrt(A_6). A_6 = N_alpha*l_alpha/I_yy > 0 perche' il braccio l_alpha
  del centro di pressione e' misurato *davanti* al CG: momento destabilizzante.
- Riga 18: `build_plant_rigid(p)` costruisce lo state-space a 4 stati
  [z, zdot, theta, thetadot], 2 ingressi [delta, alpha_w], 7 uscite (4
  misurate + 3 di plotting). E' Eq. (1) della traccia privata delle righe di
  bending.

> **Possibile domanda d'esame** -- I poli in anello aperto del modello rigido
> sono davvero solo +/- sqrt(A_6)?
> *Risposta:* No. Quelli sono i poli della **dinamica rotazionale disaccoppiata**.
> Il plant a 4 stati completo, che include l'accoppiamento con la deriva laterale
> (termine A_6/V su zdot e il termine a1*V + a4 su theta), ha poli
> [0, -1.861, +0.0291, +1.8165]. Ci sono quindi **due poli a parte reale
> positiva** (+1.8165, il modo aerodinamico, e +0.0291, un modo lento di deriva)
> piu' un **integratore in zero** (la posizione laterale z e' l'integrale di
> zdot). Lo script stampa +1.839 perche' e' il valore analitico di riferimento,
> non l'autovalore esatto del plant accoppiato.

---

## `%% Controller design (PD pitch + weak negative drift feedback)` (righe 20-31)

```matlab
[K, m] = design_controller(G, []);    % [] => ideal actuator

% Pole-placement cross-check: CL is s^2 + K1*Kd*s + (K1*Kp - A6).
wc_eq = sqrt(p.K1*K.Kp_th - p.A6);
ze_eq = p.K1*K.Kd_th/(2*wc_eq);
```

- Riga 22: il secondo argomento `[]` e' l'alias di "attuatore ideale"
  (`assemble_loop` lo converte in `tf(1)`). E' questo che rende il Task 1 il caso
  rigido puro della traccia.

### La coppia canonica di partenza (dentro `design_controller`, riga 55)

`design_controller` non parte da guadagni a caso: parte dalla **coppia in forma
chiusa** ottenuta col pole placement sulla dinamica rotazionale disaccoppiata.
Vale la pena rifare la derivazione, perche' e' la domanda d'orale piu' probabile
su questo task.

Sulla sola rotazione, con theta_ref = 0 e delta = Kp*(0 - theta) - Kd*thetadot:

    theta_ddot = A_6*theta + K_1*delta
               = A_6*theta - K_1*Kp*theta - K_1*Kd*thetadot

cioe' il polinomio caratteristico in anello chiuso

    s^2 + K_1*Kd*s + (K_1*Kp - A_6) = 0

Confrontandolo con la forma standard s^2 + 2*zeta*wn*s + wn^2:

    wn   = sqrt(K_1*Kp - A_6)
    zeta = K_1*Kd / (2*wn)

Ora si scelgono Kp e Kd per ottenere una coppia "bella". La scelta canonica
(D'Antuono Eq. 3.6-3.7) e':

    Kp0 = 2*A_6/K_1      ->  wn = sqrt(2*A_6 - A_6) = sqrt(A_6)
    Kd0 = sqrt(A_6)/K_1  ->  2*zeta*wn = sqrt(A_6)  ->  zeta = 0.5

**Interpretazione**: Kp0 e' il guadagno che **specchia il polo instabile
+sqrt(A_6) nel semipiano sinistro alla stessa frequenza naturale**, e Kd0 e'
quello che gli da' smorzamento 0.5. E' la scelta minimale: non si chiede al
controllo di essere piu' veloce dell'instabilita' che deve domare, solo di
invertirne il segno. Numericamente:

    Kp0 = 2*3.3818/4.5647 = 1.4817
    Kd0 = sqrt(3.3818)/4.5647 = 0.4029

- Righe 25-26: **cross-check a posteriori**. Dopo che il tuner ha finito, lo
  script ricalcola wn e zeta *con le stesse due formule di sopra* ma usando i
  guadagni tunati. Serve a controllare che l'auto-tuner non sia scappato in una
  regione assurda. Risultato reale:

      Kp = 1.7845, Kd = 0.4433
      wc_eq = sqrt(4.5647*1.7845 - 3.3818) = 2.18 rad/s
      ze_eq = 4.5647*0.4433/(2*2.18)        = 0.46

  Il commento a riga 27 dice "course band 1-4": la banda di crossover tipica di
  un lanciatore sta fra 1 e 4 rad/s (sopra si eccitano i modi flessibili, sotto
  non si domina l'instabilita' aerodinamica). 2.18 rad/s ci sta comodamente.

**Attenzione al limite**: wc_eq e ze_eq sono grandezze della **dinamica
disaccoppiata**. Non sono i poli veri dell'anello chiuso a 4 stati (che includono
la deriva). Sono un sanity check, non un risultato. Lo script non lo dice
esplicitamente, ma il commento a riga 24 ("Pole-placement cross-check") lo lascia
intendere.

Per confronto, i **poli veri** dell'anello chiuso a 4 stati sono due coppie:

    -0.9533 +/- 1.9047i    (wn = 2.13 rad/s, zeta = 0.45)  <- la coppia piazzata dal PD
    -0.0559 +/- 0.2329i    (wn = 0.24 rad/s, zeta = 0.23)  <- il modo lento di deriva

La coppia veloce e' ben approssimata dal cross-check disaccoppiato (2.18 / 0.46
contro 2.13 / 0.45): il seed analitico regge. La coppia lenta, invece, il modello
disaccoppiato **non la vede proprio** -- e' il modo che la retroazione di deriva
crea chiudendo l'integratore di posizione, ed e' quello che detta l'orizzonte di
simulazione (tau ~ 18 s, da cui Tend = 80 s).

### Perche' il tuner sposta i guadagni dai valori canonici

Da 1.4817/0.4029 (canonici) a 1.7845/0.4433 (tunati) c'e' un +20% su Kp e un +10%
su Kd. Il motivo e' spiegato nell'header di `design_controller`: la coppia
canonica da' 6 dB di margine sull'anello **disaccoppiato**, ma sull'anello
**completo** la retroazione di deriva erode il margine aerodinamico (scende a
circa 4 dB). Il tuner `fminsearch` rialza il guadagno finche' il margine
*classificato sull'anello vero* torna a 6 dB.

Il costo minimizzato (design_controller, riga 88) e':

    c = (|AeroGM_dB| - 6)^2 + (RigidPM_deg - 30)^2      [+ 1e4 se instabile]

- La ricerca avviene in **coordinate logaritmiche**, x = log([Kp Kd]) (riga 55),
  cosi' i guadagni restano positivi per costruzione senza vincoli espliciti.
- **Limite onesto**: il costo somma dB^2 e deg^2, cioe' pesa implicitamente
  1 dB = 1 grado. E' una scelta arbitraria (per quanto convenzionale). Il codice
  non la giustifica e non espone un peso relativo. In pratica funziona perche'
  entrambi i target vengono centrati esattamente (6.00 dB / 30.0 deg), quindi il
  minimo e' praticamente zero e la pesatura non discrimina.

- Righe 30-31: `L = m.L` (anello aperto SISO, rotto sul segnale `delta`) e
  `T = m.T` (anello chiuso completo a 4 stati). Da qui in poi lo script usa solo
  questi due oggetti: L per la frequenza, T per il tempo.

> **Possibile domanda d'esame** -- Perche' non usare `place()` o `lqr()` invece
> di un `fminsearch` su due margini?
> *Risposta:* Perche' la specifica della traccia **non e' sui poli, e' sui
> margini** ("values of KP,theta and KD,theta that yield a GM about 6 dB and PM
> about 30 deg"). Un pole placement ti da' i poli che vuoi ma non ti garantisce i
> margini, e viceversa: sull'anello *completo* (con deriva e attuatore) la
> relazione fra poli e margini non e' piu' quella della coppia disaccoppiata. La
> struttura scelta e' la piu' onesta: usa il pole placement come **seed
> analitico** e poi lascia che l'ottimizzatore centri la specifica vera sul loop
> vero.

---

## `%% Rigid-body stability margins` (righe 33-41)

- Righe 34-41: solo stampe, i margini sono gia' dentro `m` (calcolati da
  `classify_margins` dentro `design_controller`). Output reale:

      Aero gain margin  : |GM| = 6.00 dB @ 0.59 rad/s (low-freq gain-reduction)
      Rigid phase margin: PM   = 30.0 deg @ 2.45 rad/s
      Delay margin      : DM   = 213 ms
      Rigid GM / Flex margins: none (ideal actuator, no bending)
      Full 4-state closed loop stable (isstable): 1

### Perche' i margini vanno classificati e non letti con `margin()`

E' il cuore concettuale di tutto HM3. Il plant ha **due poli instabili e un
integratore**. Per il criterio di Nyquist, con P poli instabili in anello aperto
servono P *encirclement* antiorari del punto critico per essere stabili in anello
chiuso. Questo significa che il diagramma **deve** passare "dalla parte giusta"
del punto critico, e quindi:

- **abbassando** il guadagno il diagramma si abbassa e perde l'encirclement ->
  instabile;
- **alzandolo** troppo lo si perde dall'altra parte -> instabile.

L'anello e' **condizionalmente stabile**: esiste una *banda* di guadagni ammessi,
non una semiretta. Ne segue che il margine di guadagno a bassa frequenza e' un
margine di **riduzione** (gmdb < 0): e' quanto guadagno puoi *togliere* prima di
cadere. `margin()` di default restituisce un solo numero (il crossover piu'
"vicino") e su un anello del genere quel numero e' privo di significato fisico:
tipicamente becca uno degli attraversamenti di deriva a bassissima frequenza, che
non sono margini ma artefatti dell'integratore di posizione.

`classify_margins` risolve questo prendendo **tutti** gli attraversamenti da
`allmargin()` e assegnandoli a bande fisiche:

- **Aero GM** (riga 44): fra gli attraversamenti della fase critica (-180 deg,
  mod 360) con `gmdb < 0` prende
  quello a **frequenza piu' bassa** (mode `'minf'`). La maschera e'
  `gf > 0 & gf < w_flex & gmdb < 0`: **esclude solo la voce di DC** (`gf == 0`,
  l'integratore), **non** gli attraversamenti sotto `w_drift`. Qui: -6.00 dB @
  0.59 rad/s. (Che 0.59 cada *sopra* `w_drift = 0.55` e' un fatto di questa run,
  non una regola imposta dal codice -- comodo saperlo se un giorno il tuner
  spostasse il crossover.)
- **Rigid PM** (riga 49): margine di fase al crossover di corpo rigido. Qui si',
  la maschera filtra esplicitamente `pf > w_drift` (e prende il *massimo* valore,
  mode `'maxv'`). Qui: 30.0 deg @ 2.45 rad/s.
- **Rigid GM** (riga 45): margine di *aumento* di guadagno (`gmdb > 0`). Qui
  **assente** (NaN).
- **drift_w** (riga 62): gli attraversamenti di fase sotto
  `w_drift = 0.3*sqrt(A6) = 0.55 rad/s`, raccolti a parte e marcati esplicitamente
  come "non margini".

### Perche' il Rigid GM e' assente nel Task 1

Con attuatore ideale l'anello aperto e' L = Kc * G, dove il PD contribuisce uno
zero (Kp + Kd*s) e il plant rigido ha 4 poli. La fase ad alta frequenza tende
all'asintoto di -90 deg (relativo grado 1 sulla catena di assetto: il derivativo
cancella un ordine). **Non riesce mai a riavvolgersi fino a un secondo
attraversamento della fase critica -180 deg**, quindi non esiste un margine di aumento di
guadagno: puoi alzare il guadagno quanto vuoi e l'anello resta stabile. E' la
riga 40 che lo dichiara. Nel Task 2 l'attuatore TVC (2 poli), il ritardo e il
notch aggiungono ritardo di fase sufficiente a creare quel crossover, e infatti
li' spunta un Rigid GM = 7.56 dB @ 11.11 rad/s.

> **Possibile domanda d'esame** -- Il delay margin e' 213 ms. Ma nel Task 1
> l'attuatore e' ideale e non c'e' nessun ritardo. A cosa serve quel numero?
> *Risposta:* E' una **previsione**: dice quanto ritardo puro l'anello rigido
> tollererebbe prima di perdere stabilita'. Serve come budget da spendere nel
> Task 2, dove entrano un ritardo reale di 20 ms e la dinamica TVC (che a bassa
> frequenza si comporta come ritardo equivalente). 213 ms e' molto sopra il
> requisito tipico di 100 ms citato a riga 39, quindi il Task 1 lascia margine
> abbondante. Nel Task 2, dopo il notch, scende a 165 ms; e con i guadagni del
> Task 1 non ri-tunati scenderebbe a 98 ms, cioe' sotto il requisito.

---

## `%% Wind-gust time response` (righe 43-55)

```matlab
w = load_wind_profile(p, Tend=80);
r = simulate_gust_response(T, w);
```

- Riga 44: la raffica di default e' un **1-cosine discrete gust** di durata
  Tg = 3 s con onset a t0 = 1 s. L'ampiezza Vg non e' hard-coded: viene
  interpolata dalla dispersione `drywind.mat` (severita' `severe`, il default)
  alla quota di volo Alt = 15143 m. Risultato reale: **Vg = 6.38 m/s**, quindi

      alpha_w_peak = Vg/V = 6.38/937.7 = 6.80e-3 rad = 0.39 deg

  Il profilo e' vw(t) = 0.5*Vg*(1 - cos(2*pi*(t - t0)/Tg)) sulla finestra della
  raffica, zero altrove: parte da zero con derivata nulla e ci torna, quindi non
  inietta discontinuita' nell'anello.
- Il commento a riga 44 giustifica `Tend=80`: il modo dominante in anello chiuso
  ha tau ~ 18 s (wn = 0.24, zeta = 0.23), e servono ~5 costanti di tempo per
  vedere l'assetto tornare a zero. E' il **modo di deriva lento**, non la coppia
  di beccheggio a 2.18 rad/s: la deriva laterale, chiusa con guadagni da 1e-3, e'
  lentissima. Se si simulasse a 12 s (il default di `load_wind_profile`) si
  vedrebbe solo il transitorio di assetto e si perderebbe il ritorno a zero.
- Riga 45: `simulate_gust_response` fa un `lsim` di T con ingresso
  [alpha_w(t), theta_ref = 0]: e' una simulazione di **pura reiezione di
  disturbo**, senza comando di assetto.

Risultati reali:

      peak |theta| = 0.261 deg
      peak |z|     = 2.27 m
      peak |delta| = 0.528 deg
      peak |zdot|  = 0.68 m/s

- Righe 51-52: il **budget di incidenza** e il carico aerodinamico. A max-q_bar
  il parametro dimensionante non e' theta ne' delta, e' il prodotto
  **q_bar * alpha** (pressione dinamica per incidenza), che e' proporzionale al
  momento flettente alla base del lanciatore. Con q_bar = 81.1 kPa e un picco di
  incidenza di 0.577 deg, lo script riporta **46.8 kPa*deg**.
- Righe 53-55: dichiara la deriva laterale come "canale di load relief" e non
  come margine di Nichols, e la confronta con due soglie (500 m, 15 m/s) che pero'
  **non sono nel codice come costanti verificate**: sono numeri scritti a mano
  nella stringa di formato. Non c'e' nessun `assert`, nessun confronto: sono
  puramente decorativi. **E l'etichetta "load relief" e' un residuo storico
  fuorviante**: la retroazione di deriva contiene la deriva e stabilizza
  l'integratore di posizione, ma sull'incidenza incide per l'1% (vedi sotto).

### Il segno di alpha: perche' e' un MENO (e perche' ribalta la narrazione)

`simulate_gust_response.m` (riga 29) calcola:

    r.alpha = r.theta + r.zdot/w.V - r.alphaw      (segno MENO su alpha_w)

ed e' **l'unico segno coerente col plant**. La **Eq. (1) della traccia** ha la
colonna del disturbo pari a [0; -a1*V; 0; -A_6], e il codice la riproduce
identica (`build_plant_rigid.m` riga 17). Sviluppando la riga di theta_ddot:

    theta_ddot = A_6*theta + (A_6/V)*zdot + K_1*delta - A_6*alpha_w
               = A_6*(theta + zdot/V - alpha_w) + K_1*delta

cioe' l'incidenza aerodinamica che il **momento vede davvero** e'

    alpha = theta + zdot/V - alpha_w

Idem per la riga di zdot_dot, che si raggruppa in a1*V*(theta + zdot/V -
alpha_w). Il meno ha una lettura fisica immediata: `alpha` e' l'incidenza rispetto
alla velocita' **relativa all'aria**, quindi un vento laterale si **sottrae** allo
zdot del veicolo. Il caso limite lo conferma: con un gradino di vento sostenuto e
il feedback di deriva spento, la simulazione converge a `zdot_ss = +V*alpha_w`
(rapporto verificato 1.0000) con `theta_ss = 0` e dunque **alpha_ss = 0**: il
veicolo viene spinto sottovento finche' non vola *insieme* all'aria, e l'incidenza
relativa si annulla.

**Nota storica (buona da sapere all'orale).** Fino a poco fa questa riga aveva un
`+ r.alphaw`. Era un **bug di post-processing**, non di modello: il plant e' sempre
stato giusto, sbagliava solo la formula con cui si ri-costruiva `alpha` a valle di
`lsim`. Corretto. Cosa cambia e cosa no:

    alpha col vecchio + (buggato):  picco 0.255 deg  ->  q_bar*alpha = 20.7 kPa*deg
    alpha col meno (corretto):      picco 0.577 deg  ->  q_bar*alpha = 46.8 kPa*deg

Un **fattore 2.3** sul carico. **Non cambia nulla** di margini, Nichols,
stabilita', time history di theta/z/zdot/delta e Task 3: `alpha_w` non entra in
L(s), e la riga e' pura post-elaborazione. Cambiano solo i due numeri di carico
(righe 51-52) e la figura f3 -- e con essi la **narrazione fisica**.

### La conseguenza: un puro attitude hold e' LOAD-AGGRAVATING

Il picco di incidenza totale (**0.577 deg**) **supera** il contributo del solo
vento (**0.390 deg**). Il vento genera un momento -A_6*alpha_w che spinge il muso
in negativo; il PD, che ha come unico obiettivo theta -> 0, reagisce **beccheggiando
il muso dentro il vento relativo**, e quel theta (-0.178 deg al picco di raffica)
entra in alpha = theta + zdot/V - alpha_w con lo **stesso segno** di -alpha_w: si
**somma** al vento invece di cancellarlo.

Una legge di puro attitude-hold e' quindi **load-aggravating**, non
load-relieving -- ed e' esattamente il motivo per cui i lanciatori reali
aggiungono un **termine esplicito di load relief** (accelerometro laterale o alpha
stimata). Il velivolo non offre aiuto: con A_6 > 0 il centro di pressione sta
davanti al baricentro, il momento aerodinamico e' divergente, non c'e' nessuna
stabilita' a banderuola che riallinei spontaneamente il muso col vento.

> **Possibile domanda d'esame** -- Il picco di q_bar*alpha e' 46.8 kPa*deg. E'
> tanto o poco?
> *Risposta:* Poco: circa un ordine di grandezza sotto i limiti tipici di un
> lanciatore snello (alcune centinaia di kPa*deg). Il motivo pero' non e' merito
> del controllore -- che anzi **aggrava** il carico rispetto al solo vento
> (0.577 vs 0.390 deg di incidenza) -- ma della raffica, che e' modesta:
> Vg = 6.4 m/s, la sigma della dispersione `drywind` a 15 km, non un caso di vento
> estremo. Il caso di carico vero si ottiene con
> `load_wind_profile(..., 'profile','strongwind')`, che usa il generatore di vento
> del professore.

---

## `%% Figures` (righe 57-99)

Tre figure, tutte con `'Color','w'` esplicito.

### f1 -- Nichols (righe 62-65)

Delegata a `plot_nichols_lv(L, m, ...)`, con `wrange = [1e-2 1e2]` e
`xlim = [-720 0]`. La convenzione e' quella degli appunti del corso (lo dichiara
il commento alle righe 58-61): il punto critico sta a **(-180 deg, 0 dB)**,
derivato da 1 + L = 0 <=> L = -1. La fase e' **phase-matched** in modo che il
crossover rigido cada sul ramo che contiene -180: lo shift e' un **multiplo
esatto di 360 deg** (`sh = ph + 360*round(...)`, riga 49 di `plot_nichols_lv`),
quindi e' una pura **rietichettatura del ramo di fase**, non una rotazione: la
curva non viene distorta (nel Task 1 lo shift risulta addirittura 0: il ramo
naturale di MATLAB e' gia' quello del corso). Verificato sulla figura:
l'attraversamento aerodinamico casca esattamente a -180 deg, 6 dB sopra il punto
critico, e il crossover rigido a -150 deg, cioe' 30 deg a destra del punto
critico -- il PM si legge direttamente dal grafico.

(D'Antuono Fig. 3.2 mostra la **stessa carta rietichettata di +360 deg**, punto
critico a +180: +180 e -180 sono lo stesso punto, fase mod 360. Il codice usava
quella rietichettatura fino a poco fa ed e' stato riallineato alla convenzione
del corso; nessun margine cambia, e' solo l'etichetta dell'asse.)

Il senso della figura: l'anello **viene dall'alto** (a DC il guadagno tende
a infinito per via dell'integratore di posizione laterale z) e passa **sopra**
il punto critico con 6 dB di scarto. E' la firma grafica della stabilita'
condizionale (a un solo lato, nel Task 1: manca l'attraversamento superiore).
I marker sovrapposti sono i margini classificati
(quadrato = Aero GM, rombo = Rigid PM) e le croci nere sono gli attraversamenti
di deriva, etichettati letteralmente "not a margin".

### f2 -- Gust response (righe 68-82)

Tiled 2x2: theta, z, zdot, delta. Sono **esattamente le quattro variabili che la
traccia elenca**. Nessuna elaborazione, plot diretti di `r`.

### f3 -- Alpha budget e carico (righe 86-99)

Il primo pannello scompone alpha nei suoi tre contributi (theta, zdot/V,
alpha_w), il secondo mostra q_bar*alpha (picco **46.8 kPa*deg**). Il commento a
riga 84 e il titolo del pannello a riga 96 riportano entrambi la formula col segno
giusto, `alpha = theta + zdot/V - alpha_w`.

E' **la figura piu' interessante di tutto il Task 1**, ed e' quella da mostrare
all'orale: si vede a occhio che la curva di `theta` e quella di `alpha_w` vanno
nella **stessa direzione** nel bilancio (theta negativo, -alpha_w negativo), e che
la `alpha` risultante e' **piu' grande in modulo** del solo contributo del vento.
E' la prova grafica del load aggravation discusso sopra.

### Export (righe 101-113)

- Riga 102: `fig_dir` risolto via `mfilename('fullpath')`, quindi lo script
  scrive sempre in `HM3/figures/` **indipendentemente dalla working directory**.
  Buona pratica, e la ragione per cui il professore puo' lanciarlo da qualunque
  cartella.
- Righe 105-109: `theme(f,'light')` dentro un `try/catch` con fallback su
  `f.Color = 'w'`. Serve a forzare figure chiare anche se il desktop MATLAB e' in
  dark mode; `theme()` esiste solo da R2025a, da cui il catch.
- Righe 110-111: `exportgraphics` a 200 dpi, nomi `task1_<Name>.png` dove `<Name>` e'
  la proprieta' `'Name'` della figura.

---

## Possibili domande d'esame

**D: Cosa vuol dire che l'anello e' "condizionalmente stabile" e come lo si vede
sul Nichols?**
R: Vuol dire che la stabilita' richiede il guadagno dentro una **banda**, non
sopra una soglia: sia riducendolo sia aumentandolo troppo si perde. Nasce dal
fatto che il plant ha poli instabili in anello aperto (+1.82 e +0.03 rad/s, piu'
un integratore), quindi per Nyquist il diagramma deve produrre un numero preciso
di encirclement del punto critico. Sul Nichols si riconosce perche' la curva
**arriva dall'alto** (guadagno infinito a DC) e attraversa la fase critica
-180 deg **sopra** il punto critico (-180, 0 dB), senza toccarlo: c'e' margine
di guadagno verso il basso e margine di fase al crossover. Il margine
aerodinamico e' quello di **riduzione** (gmdb < 0), qui 6.00 dB a 0.59 rad/s.

**D: Perche' i guadagni di deriva Kp_z e Kd_z sono negativi e cosi' piccoli, e
cosa fanno davvero?**
R: Il segno negativo lo impone la traccia ("must be negative in the order of
10^-3"). La loro funzione primaria e' un **requisito di stabilita'**, non di load
relief: chiudono l'**integratore libero della posizione laterale** (il plant ha un
polo esatto in s = 0, perche' z e' l'integrale di zdot). Verifica diretta,
spegnendoli (Kp_z = Kd_z = 0):

| | poli di anello chiuso | picco alpha | picco z |
|---|---|---|---|
| con drift feedback | -0.056 +/- 0.233i, -0.953 +/- 1.905i | 0.577 deg | 2.27 m |
| senza (Kp_z = Kd_z = 0) | **0**, -0.076, -0.98 +/- 1.93i | 0.584 deg | **9.54 m** |

Senza di loro resta un polo **esattamente nell'origine**: il sistema e' solo
**marginalmente stabile** e z non torna mai a zero. Sull'incidenza incidono per
l'**1%** (0.584 -> 0.577 deg: irrilevante), ma tagliano la deriva **da 9.5 m a
2.3 m**. Chiamarli "load relief" (come fanno i commenti del codice e il README) e'
fuorviante: un vero load relief richiederebbe un termine dedicato in alpha o
accelerazione laterale. Il valore piccolo (1e-3) serve perche' il canale di deriva
non interferisca con la banda di assetto.

**D: L'auto-tuner centra 6.00 dB e 30.0 deg esattamente. Non e' sospetto?**
R: No, e' atteso: il costo e' una somma di due quadrati con due gradi di liberta'
(Kp, Kd) e due target. Il sistema e' quadrato, quindi esiste generalmente una
soluzione a costo (quasi) nullo e `fminsearch` la trova. La conseguenza pero' e'
che **la pesatura relativa dB^2 vs deg^2 non conta**: se i target fossero
incompatibili (piu' target che gradi di liberta') quella pesatura arbitraria
diventerebbe determinante e andrebbe giustificata.

**D: Che differenza c'e' fra i guadagni canonici (1.48, 0.40) e quelli tunati
(1.78, 0.44), e perche' il tuner li alza?**
R: I canonici Kp0 = 2*A_6/K_1 e Kd0 = sqrt(A_6)/K_1 sono esatti sulla dinamica
rotazionale **disaccoppiata** e li' danno 6 dB di margine con zeta = 0.5. Ma
sull'anello **completo** l'aggiunta della retroazione di deriva (e del suo
integratore) abbassa il margine aerodinamico a circa 4 dB. Il tuner alza Kp del
20% per recuperare i 6 dB sul loop vero. La lezione e': la formula chiusa e' un
ottimo **seed**, non il progetto finale.

**D: Il codice riporta un picco di alpha di 0.577 deg, ma il vento da solo ne fa
0.390. Come fa l'incidenza totale a essere PIU' GRANDE del disturbo?**
R: Perche' il controllore la peggiora. Il bilancio e' alpha = theta + zdot/V -
alpha_w (il meno viene dalla colonna di disturbo [0; -a1*V; 0; -A_6] di Eq. (1):
alpha e' l'incidenza rispetto all'aria, quindi il vento si sottrae allo zdot del
veicolo). Il momento -A_6*alpha_w spinge il muso in negativo; il PD di assetto,
che vuole solo theta -> 0, reagisce ruotando il muso **dentro** il vento relativo,
e quel theta (-0.178 deg al picco) entra nel bilancio con lo **stesso segno** di
-alpha_w. I due contributi si **sommano**. Morale: un puro attitude hold e'
**load-aggravating**; il load relief vero richiede un termine esplicito
(accelerometro laterale / alpha stimata) che questo controllore non ha. Nota
storica: fino a poco fa `simulate_gust_response.m` aveva un `+` su alpha_w e
riportava 0.255 deg / 20.7 kPa*deg, facendo sembrare l'anello load-relieving. Era
un bug di **post-processing** -- il plant e' sempre stato giusto -- e infatti la
correzione non tocca ne' i margini ne' le time history di theta, z, zdot, delta.

**D: Perche' simulare 80 s se la raffica dura 3 s?**
R: Perche' il modo dominante in anello chiuso non e' quello di beccheggio
(2.18 rad/s, tau < 1 s) ma il **modo lento di deriva laterale**, con wn = 0.24
rad/s e zeta = 0.23, cioe' tau ~ 18 s. Servono ~5 costanti di tempo per vedere
l'assetto e la deriva rientrare. Con il default di 12 s si vedrebbe solo il
transitorio veloce e si perderebbe la parte piu' interessante (il ritorno a zero
e il picco di deriva).
