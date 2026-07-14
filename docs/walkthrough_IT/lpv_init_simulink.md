# HM3/LTV_FULL_ASCENT/init_simulink_lpv.m

## Ruolo del file nel progetto

Questo file e' il **setup unico** dello showcase LPV "full ascent" (la parte di
HM3 che va oltre la traccia: il design a max-q sollevato a tutta l'ascesa,
0-140 s). Fa una cosa sola ma la fa per due consumatori diversi:

1. **carica e pre-elabora tutti i dati che dipendono dal tempo** -- i
   coefficienti del plant LPV, il guadagno schedulato, il vento generato dal
   modello del professore, i dati del modo flessibile, il notch e il TVC;
2. li **restituisce due volte**, in due formati:
   - come **struct `S` di `griddedInterpolant`** -> e' quello che consumano i
     RHS `ode_lpv_ascent` / `ode_lpv_flex` dentro `ode45` (la *baseline*, la
     "fonte di verita'" dichiarata dal README);
   - come **variabili nel base workspace** (righe 119-138) -> e' quello che
     consumano i due modelli Simulink `hm3_full_ascent.slx` e
     `hm3_full_ascent_flex.slx`, che referenziano quelle variabili per nome
     dentro i blocchi Lookup Table, Gain, Constant e Transfer Fcn.

Il punto chiave, e la ragione per cui questo file esiste come modulo separato,
e' proprio il **numero 2**: script e modello Simulink devono partire *dagli
stessi identici numeri*, altrimenti il confronto Simulink-vs-ode45 non
misurerebbe l'errore di integrazione ma una differenza di dati. `init_simulink_lpv`
e' il collo di bottiglia comune, il punto in cui i due binari si toccano.

E' il corrispettivo LPV di `HM3/init_simulink_hm3.m` (quello congelato a
`t_ref = 72 s`): la' le matrici di `build_plant_rigid` sono **costanti**, qui
diventano **storie temporali** campionate da `General/hw3-v3/GreensiteLPV_DATA.mat`
(151 campioni, passo 1 s, t da 0 a 150 s -- verificato caricando il `.mat`).

Chiamato da: `main_full_ascent.m`, `main_flex.m`, `main_q_scheduling.m`,
`run_full_ascent_simulink.m`, `run_flex_simulink.m`. Dipende da (cartella
padre `HM3/`): `load_hw3_params`, `build_plant_rigid`, `design_controller`,
`build_tvc`.

---

## Firma, contratto e opzioni (righe 1-37)

```matlab
function S = init_simulink_lpv(o)
arguments
    o.tsched_step (1,1) ... = 5      % passo griglia schedule [s]
    o.t0          (1,1) ... = 5      % inizio schedule / sim [s]
    o.Tstop       (1,1) ... = 140    % orizzonte [s]
    o.push        (1,1) logical = true
end
```

- Riga 1: unica uscita `S`, tutto il resto (le variabili base) e' un
  **side effect** controllato da `o.push`. Chi vuole solo la baseline ode45 puo'
  chiamare `init_simulink_lpv('push', false)`.
- Righe 2-30: l'header documenta il contratto piu' importante del file, quello
  dei **coefficienti efficaci** (righe 7-12):

      c1 = a1        (* zdot)      c5 = A6/V   (* zdot)
      c2 = a1*V + a4 (* theta)     c6 = A6     (* theta e * (-alpha_w))
      c3 = a3        (* delta)     c7 = K1     (* delta)
      c4 = a1*V      (* alpha_w)   invV = 1/V  (alpha_w = v_w * invV)

  cioe' **ogni termine dell'equazione e' UN lookup moltiplicato per UN segnale**.
  Non e' un dettaglio estetico: e' la scelta che rende possibile l'accordo a
  1e-7 rad (vedi la sezione sui coefficienti, sotto).
- Righe 39-42: `addpath` della cartella padre `HM3/` -- il file riusa gli helper
  dell'homework "vero" invece di duplicarli.

---

## Caricamento dati di riferimento (righe 44-47)

```matlab
D = load(fullfile(gdir,'GreensiteLPV_DATA.mat'));  L = D.GreensiteLPV;
W = load(fullfile(gdir,'drywind.mat'));            drywind = W.drywind;
GreensiteLPV = L;                     % lo vuole anche il generatore
```

- Riga 45: `GreensiteLPV` e' una struct di `timeseries` (campi `A6, K1, a1, a3,
  a4, sigma_ins, phi_tvc, phi_ins, omega, Q, V, Tc, h, aqk, Mach`).
- Riga 47: e' un **alias di comodo** (`GreensiteLPV` e `L` sono lo stesso
  dataset), non un vincolo funzionale. Il generatore riceve il dataset come
  **argomento** (riga 88: `run_wind_generator(gdir, o.Tstop, drywind, L)`), e il
  nome viene imposto come **stringa letterale** sia da
  `in.setVariable('GreensiteLPV', ...)` (riga 164) sia dal fieldname del push in
  base workspace (riga 127): il nome della variabile *chiamante* e' quindi
  irrilevante. Cio' che conta davvero e' che nel **base workspace** esista una
  variabile di nome `GreensiteLPV`, perche' la copia del generatore dentro i
  `.slx` la risolve per nome.

---

## Coefficienti efficaci e breakpoints (righe 49-68)

```matlab
tg = L.V.Time(:);  tg = tg(tg <= o.Tstop + eps);
at = @(f) interp1(L.(f).Time, squeeze(L.(f).Data), tg);
...
Vsafe = max(V, 1);          % V(0)=0: guard A6/V e 1/V al lift-off
c1 = a1;  c2 = a1.*V + a4;  c3 = a3;  c4 = a1.*Vsafe;
c5 = A6./Vsafe;  c6 = A6;   c7 = K1;  invV = 1./Vsafe;
```

- Righe 50-51: i **breakpoints** sono i tempi del dataset stesso, tagliati a
  `Tstop`. Con i default: 141 punti, `t = 0, 1, ..., 140 s`. Non c'e'
  raffinamento: la griglia della lookup e' la griglia dei dati.
- Riga 52: `at` ricampiona ogni serie su `tg` con `interp1` lineare. Poiche' le
  serie condividono gia' la stessa base tempi, qui e' di fatto un'identita'.
- Righe 53-58: campionamento di `V, A6, K1, a1, a3, a4, Q, h` e -- per il
  modello flessibile (T008) -- `omega` (frequenza del primo modo di bending),
  `sigma_ins`, `phi_ins` (accoppiamento INS-bending) e `aqk` (forzante del TVC
  sul bending).
- Righe 61-68: le equazioni del plant rigido sono

      zddot     = a1*zdot + (a1*V + a4)*theta + a3*delta - a1*V*alpha_w
      thetaddot = (A6/V)*zdot + A6*theta      + K1*delta - A6*alpha_w

  e i coefficienti `c1..c7` sono esattamente i **fattori davanti a ciascun
  segnale**, gia' combinati. Nota che `c6` compare *due volte* nella seconda
  equazione (davanti a `theta` e davanti a `-alpha_w`): infatti nel modello
  Simulink un solo lookup `c6` alimenta due Product (`P6` e `P8`).

**Perche' non tabulare `a1`, `V`, `a4` separatamente e moltiplicare a runtime?**
Perche' **l'interpolazione non commuta con il prodotto**:

    interp(a1*V)(t)  !=  interp(a1)(t) * interp(V)(t)     per t fra due breakpoint

Su questo dataset (breakpoint a 1 s) ho misurato lo scarto:

| quantita' | max scarto assoluto | scarto relativo |
|---|---|---|
| `c2 = a1*V + a4` | 2.24e-3 | 5.5e-5 |
| `invV = 1/V` | 1.23e-7 | 1.1e-4 |
| `c5 = A6/V` | 4.9e-7 | ~4e-4 |
| `omega^2` | 6.8e-3 | 1.9e-5 |

Se il modello Simulink usasse `c2` (tabella del prodotto) e l'ODE usasse
`fa1(t)*fV(t) + fa4(t)` (prodotto delle tabelle), i due starebbero integrando
**due sistemi LTV diversi**, e il residuo del confronto non potrebbe mai
scendere sotto ~5e-5 in relativo. Con i coefficienti efficaci, invece,
`ode_lpv_ascent` (righe 26-28) chiama `M.fc1..M.fc7`, cioe' i
`griddedInterpolant` costruiti sulle **stesse identiche tabelle** che finiscono
nei blocchi Lookup -- e i due sistemi coincidono punto per punto.

> **Possibile domanda d'esame** -- Perche' hai tabulato `a1*V + a4` come un
> unico coefficiente invece di tabulare `a1`, `V`, `a4` e fare il prodotto nel
> modello?
> *Risposta:* Perche' l'interpolazione lineare non commuta con la
> moltiplicazione: la tabella del prodotto e il prodotto delle tabelle
> coincidono solo sui breakpoint, e in mezzo differiscono (qui ~5e-5 in
> relativo). Tabulando il coefficiente *efficace* il modello Simulink e il RHS
> ode45 vedono la stessa funzione del tempo, quindi integrano lo stesso sistema
> LTV: il residuo del confronto misura solo l'errore numerico, non una
> differenza di modello.

> **Attenzione -- il guard `Vsafe` e' codice morto e il commento e' sbagliato.**
> La riga 59 dice `V(0)=0: guard A6/V and 1/V at lift-off`. Ho caricato il
> `.mat` in repo: `min(V) = 410.4 m/s` (a t = 7 s), `V(0) = 413.0 m/s`. Quindi
> `max(V,1) == V` identicamente e il guard **non entra mai in funzione**. E'
> innocuo, ma il commento descrive un dataset diverso da quello presente.
> (Nota collaterale: i primi campioni del dataset sono comunque strani --
> `Q(0) = -16.4 Pa`, negativa, e `Mach(0) = 0.03` mentre `V(0) = 413 m/s`:
> incoerenze del dato di partenza, non del codice.)

---

## Controllore congelato di riferimento (righe 70-72)

```matlab
p0 = load_hw3_params();                       % t_ref = 72 s
K0 = design_controller(build_plant_rigid(p0), [], 'verbose', false);
```

- Riga 71: `load_hw3_params` senza argomenti -> `t_ref = 72 s`, il punto di
  max-q su cui e' costruito tutto HM3.
- Riga 72: `Wact = []` -> attuatore ideale (Task 1 di HM3). `K0` e' la coppia PD
  di pitch *piu'* i guadagni di drift `Kp_z = Kd_z = -1e-3` (default di
  `design_controller`, righe 28-29 di quel file). Questi ultimi **non vengono
  mai schedulati** in nessuno dei due modelli: restano blocchi Gain costanti.

---

## Gain schedule (righe 74-85)

```matlab
tsched   = (o.t0:o.tsched_step:o.Tstop).';    % 5:5:140 -> 28 punti
Kprev = [2.0 1.4];                            % first warm start
for i = 1:numel(tsched)
    pk = load_hw3_params('t_ref', tsched(i));
    Kk = design_controller(build_plant_rigid(pk), [], 'K0', Kprev, ...);
    Kp_sched(i) = Kk.Kp_th;  Kd_sched(i) = Kk.Kd_th;
    Kprev = [Kk.Kp_th Kk.Kd_th];              % continuation
end
```

- Riga 75: griglia dello schedule, 5 s di passo (default), da `t0 = 5` a 140 s
  -> 28 punti. Molto piu' rada dei breakpoint del plant (1 s): i guadagni
  variano lentamente, i coefficienti aerodinamici no.
- Righe 77-78: `margin()` emette un warning su ogni loop condizionatamente
  stabile; viene silenziato e ripristinato con `onCleanup` (pattern corretto:
  il warning torna al suo stato anche se la funzione esce per errore).
- Righe 80-85: la logica e' **"design at frozen plants"**, la ricetta classica
  del gain scheduling: si congela il plant a ogni istante `tsched(i)`, si
  ri-progetta il PD con lo stesso tuner di HM3 (`design_controller`, che
  ritocca `Kp_th, Kd_th` con `fminsearch` finche' GM e PM classificati non
  colpiscono 6 dB / 30 gradi), e si tabula il risultato contro il tempo.

> **Attenzione -- la "continuation" (warm start) non esiste.**
> Le righe 79 e 84 costruiscono `Kprev`, che la riga 82 passa come `'K0'`, e il commento
> (e il README della cartella) lo chiamano *continuation*. Ma
> `design_controller.m` dichiara `o.K0 (1,2) ... = [0 0]   % accepted, unused`
> (riga 32) e nell'header `K0  ignored (kept for call compatibility)` (riga 19):
> **il valore viene ignorato**. Il `fminsearch` riparte ogni volta dal punto
> analitico di D'Antuono `x0 = log([2*A6/K1, sqrt(A6)/K1])` (riga 55 di
> `design_controller`). Ogni punto dello schedule e' quindi progettato in modo
> **indipendente**, non a continuazione. Funziona lo stesso (il punto di
> partenza analitico e' gia' un'ottima stima), ma commento e README sono
> fuorvianti e vanno corretti.

---

## Vento: il generatore del professore, una volta sola (righe 87-90 e 141-178)

```matlab
wg     = run_wind_generator(gdir, o.Tstop, drywind, L);
Vwg    = max(interp1(L.V.Time, squeeze(L.V.Data), wg.t), 1);
alphaw = wg.vw ./ Vwg;                 % angolo d'attacco da vento
```

- Riga 88: simula `strong_wind.slx` sull'intero orizzonte. In HM3 "vero", invece,
  `load_wind_profile` ritaglia solo 12 s di vento attorno a max-q: qui il vento
  e' full-ascent perche' anche il plant lo e'.
- Righe 89-90: `alpha_w = v_w / V(t)` -- piccolo angolo: la componente di vento
  normale alla traiettoria divisa per la velocita' relativa produce un angolo
  d'attacco apparente.

### `run_wind_generator` (righe 141-178)

- Righe 153-154: `load_system` + `onCleanup(@() close_system('strong_wind', 0))`.
  Lo `0` significa **chiudi senza salvare**: il modello del professore viene
  modificato in memoria (si accende il logging sulle porte) ma **mai riscritto
  su disco**. E' il modo pulito di rispettare il vincolo "non toccare il
  materiale del corso".
- Righe 156-160: il logging viene messo sulle **porte sorgente** del sottosistema
  (`get_param(..., 'PortHandles')` -> `set_param(ph.Outport(k), 'DataLogging',
  'on')`), non su blocchi To Workspace aggiunti. Motivo pratico: aggiungere
  blocchi al modello altrui e' invasivo, mentre marcare una porta e' una
  proprieta' del segnale che si perde alla chiusura senza salvataggio.
- Righe 162-167: `Simulink.SimulationInput` con `setVariable` per iniettare
  `drywind` e `GreensiteLPV` (isolamento: non si sporca il base workspace del
  chiamante) e `setModelParameter` per `StopTime` e `SignalLogging`.
- Righe 169-177: le due uscite sono la **raffica media** `v_wp` (profilo di vento
  medio) e la **turbolenza** (Dryden, schedulata in quota). I log a passo
  variabile possono contenere **tempi ripetuti** (riga 171: `unique`), quindi si
  costruisce la griglia unione e si somma:

      v_w(t) = v_wp(t) + turbolenza(t)

  Il commento (riga 145) parla di *fixed seeds*: i seed stanno dentro i blocchi
  del generatore del professore, il codice qui **non li imposta** -- la
  riproducibilita' e' ereditata, non garantita da questo file.

> **Possibile domanda d'esame** -- Il vento e' generato una volta qui e anche,
> di nuovo, dentro il modello Simulink (che contiene una copia del generatore).
> Le due realizzazioni sono identiche?
> *Risposta:* No, non necessariamente. I seed sono fissi, quindi la *sequenza di
> rumore* e' la stessa, ma i filtri di Dryden sono a stati continui e vengono
> integrati con solver e passi diversi nei due modelli (`strong_wind.slx`
> standalone contro la copia dentro `hm3_full_ascent.slx`, che gira con
> `MaxStep = 0.02`). Le uscite differiscono quindi di poco. E' esattamente per
> questo che `run_full_ascent_simulink` **non** usa `S.windfun`, ma ripiega
> l'`alpha_w` *loggato dal modello stesso*.

---

## Struct di ritorno e interpolanti (righe 92-117)

```matlab
gi = @(y) griddedInterpolant(tg, y, 'linear', 'nearest');
S.fc1 = gi(c1); ... S.fc7 = gi(c7);
S.fKp = griddedInterpolant(tsched, Kp_sched, 'linear', 'nearest');
S.windfun = griddedInterpolant(wg.t, alphaw, 'linear', 'nearest');
```

- Riga 101: `'linear'` = metodo di **interpolazione**, `'nearest'` = metodo di
  **estrapolazione** (tiene il valore dell'estremo fuori dal supporto). Questa
  seconda scelta e' quella che permette allo schedule (che parte da `t0 = 5 s`)
  di essere valutato anche per `t < 5 s`: tiene `Kp_sched(1)`. Nel modello
  Simulink lo stesso comportamento e' ottenuto con `'ExtrapMethod','Clip'` sui
  lookup dei guadagni (`build_hm3_full_ascent.m`, righe 108-109) -- e il
  commento la' lo dichiara esplicitamente: e' cio' che rende *esatto* l'overlay
  del caso scheduled anche nei primi 5 s.
- Righe 110-111: per il modello flessibile vengono esposti anche i coefficienti
  **grezzi** `fa1, fa3, fa4, fA6, fK1, fomega, faqk, fsig, fphi`. **Qui si
  rompe la simmetria**: `ode_lpv_flex` (righe 19, 38-39) usa questi grezzi e
  ricostruisce `(a1*V + a4)` e `A6/V` a runtime, mentre
  `hm3_full_ascent_flex.slx` continua a usare i lookup `c1..c7` gia' combinati.
  Fra i breakpoint i due **non sono lo stesso sistema** (vedi la tabella degli
  scarti sopra). E' la spiegazione piu' probabile del perche' il residuo
  flessibile dichiarato dal README (5e-7 rad) sia ~5 volte quello rigido
  (1.1e-7 / 2.2e-7): non e' solo tolleranza del solver, e' anche una piccola
  incoerenza di modello.
- Righe 112-115: il TVC (`build_tvc(p0, 3)` = attuatore del 2o ordine per un
  Pade del 3o ordine sul ritardo di 20 ms -> **5 stati**) viene esportato in
  **due forme**: `ssdata` (per il RHS ode45, che integra `x_tvc`) e `tfdata`
  (per il blocco Transfer Fcn del modello). E' **LTI**, quindi le due
  realizzazioni hanno lo stesso comportamento ingresso-uscita e la
  duplicazione e' sicura (vedi la domanda d'esame in fondo).
- Righe 116-117: le costanti del notch profondo di HM3: `zN = 0.002` (numeratore
  poco smorzato -> notch stretto e profondo), `zD = 0.7` (denominatore ben
  smorzato), piu' `wn72 = omega(72) = 18.9 rad/s`, la frequenza di bending al
  punto di design -- il riferimento contro cui si misura il detuning.

---

## Push al base workspace (righe 119-138)

```matlab
base = struct('lpv_t', tg, 'lpv_c1', c1, ..., 'sched', 0, ...
              'drywind', drywind, 'GreensiteLPV', GreensiteLPV, ...);
fn = fieldnames(base);
for i = 1:numel(fn), assignin('base', fn{i}, base.(fn{i})); end
```

- Righe 121-133: **il contratto con i due `.slx`**. Ogni nome qui e' referenziato
  letteralmente in un blocco:
  - `lpv_t` -> `BreakpointsForDimension1` di tutte le lookup;
  - `lpv_c1..lpv_c7`, `lpv_invV` -> `Table` delle lookup dei coefficienti;
  - `Kp_th0, Kd_th0, Kp_z0, Kd_z0` -> blocchi Constant / Gain (guadagni congelati);
  - `tsched, Kp_sched, Kd_sched` -> breakpoints e tabelle dei lookup dei guadagni;
  - `sched` -> il blocco Constant che pilota i due Switch (0 = frozen, 1 = scheduled);
  - `drywind`, `GreensiteLPV` -> li richiede la copia del generatore;
  - `Tstart, Tstop` -> `StopTime` del modello e' letteralmente la stringa `'Tstop'`;
  - `lpv_omega, lpv_omega2, lpv_2zBMw, lpv_aqk, lpv_sig, lpv_phi`,
    `notch_zN, notch_zD`, `tvc_num, tvc_den` -> solo modello flessibile.
- Riga 130: nota che vengono spinte **sia** `lpv_omega` **sia** `lpv_omega2 =
  omega.^2`. Il modello flessibile usa la tabella `omega2` per il termine di
  rigidezza `-omega^2*eta` e per il notch; `ode_lpv_flex` invece calcola
  `w = M.fomega(t)` e poi `w^2` (riga 40). Di nuovo: `interp(omega^2) !=
  interp(omega)^2` (scarto relativo misurato 1.9e-5). `lpv_2zBMw = 2*zBM*omega`
  e' invece **lineare** in `omega`, quindi tabella e prodotto coincidono
  esattamente: la discrepanza c'e' solo sui termini quadratici.
- Riga 135: `assignin('base', ...)` in un ciclo. E' l'accoppiamento
  script-modello piu' fragile che c'e' (variabili globali per nome, nessun
  controllo a compile-time), ma e' anche il modo standard di parametrizzare un
  modello Simulink senza Data Dictionary. Conseguenza pratica: **nessuno dei due
  `.slx` puo' essere simulato se prima non si e' chiamato `init_simulink_lpv`**.

> **Possibile domanda d'esame** -- Come si rappresenta un plant *LPV* in
> Simulink, visto che il blocco State-Space accetta solo matrici costanti?
> *Risposta:* Qui si scompone la dinamica in blocchi elementari: **Clock ->
> Lookup Table 1-D (breakpoints = tempo di volo, tabella = coefficiente) ->
> Product -> Sum -> Integrator**. Ogni coefficiente variabile diventa un segnale
> generato da una lookup e moltiplicato per lo stato corrispondente da un blocco
> Product; gli integratori chiudono la catena. Le alternative (il blocco *LPV
> System*, che interpola un array di `ss`, oppure un *MATLAB Function block* che
> valuta il RHS) funzionerebbero ma nascondono i coefficienti dentro un oggetto
> o dentro del codice: la scelta fatta li lascia tutti ispezionabili sul
> canvas -- ed e' anche lo stile del generatore di vento del professore. Il
> codice non discute esplicitamente le alternative: questa e' la lettura della
> scelta fatta.

---

## Possibili domande d'esame

**D: Che cos'e' un plant LPV e perche' HM3 "vero" non ne aveva bisogno?**
R: LPV = *Linear Parameter-Varying*: la struttura del modello resta lineare ma
le matrici dipendono da un parametro che varia (qui il tempo di volo, che fa da
proxy per Mach, quota e pressione dinamica). HM3 congela tutto a `t = 72 s`
(max-q) e ottiene un LTI, perche' la traccia chiede un design di punto: e'
legittimo perche' max-q e' l'istante peggiore (`A6`, cioe' l'instabilita'
aerodinamica, e' massima li'). Questo showcase toglie il congelamento per
vedere quanto lontano si estende quel design.

**D: Su cosa si schedulano i guadagni, e come e' realizzata la schedulazione nel
modello?**
R: Su **tempo di volo**, con una Lookup Table 1-D per `Kp_theta` e una per
`Kd_theta`: breakpoints `tsched = 5:5:140` (28 punti), tabelle `Kp_sched` /
`Kd_sched`, interpolazione **lineare** ed estrapolazione **Clip** sotto i 5 s.
Un blocco Switch pilotato dalla costante `sched` sceglie fra il ramo schedulato
e la coppia congelata di max-q (`Kp_th0`, `Kd_th0`). I guadagni di drift
`Kp_z, Kd_z` restano fissi in entrambi i casi. Il file spinge anche `lpv_Q`, per
cui una lookup indicizzata sulla pressione dinamica sarebbe possibile -- ma i
due modelli documentati qui non la usano (lo studio su `q` sta in
`main_q_scheduling.m`, e conclude che `q` e' una cattiva variabile di
scheduling perche' non e' monotona e genera isteresi).

**D: Come vengono calcolati i guadagni schedulati?**
R: Con la ricetta "design at frozen plants": per ogni `t` della griglia si
ricostruisce il plant congelato (`build_plant_rigid(load_hw3_params('t_ref',
t))`) e si ri-esegue il tuner di HM3 (`design_controller`), che cerca con
`fminsearch` la coppia `(Kp_th, Kd_th)` che porta il gain margin aerodinamico a
6 dB e il phase margin rigido a 30 gradi. Il risultato viene tabulato contro il
tempo. Il codice *dice* di usare una continuation (warm start dalla soluzione
precedente), ma `design_controller` ignora l'opzione `K0`: ogni punto riparte
dalla formula chiusa di D'Antuono `Kp = 2*A6/K1`, `Kd = sqrt(A6)/K1`.

**D: Perche' il TVC viene esportato sia come `ss` sia come `tf`? Non e' una
duplicazione pericolosa?**
R: E' una duplicazione, ma sicura, **perche' il TVC e' LTI**: due realizzazioni
diverse (stato-spazio per `ode45`, forma canonica interna del blocco Transfer
Fcn) della stessa funzione di trasferimento hanno lo stesso comportamento
ingresso-uscita, quindi lo stato interno puo' benissimo essere diverso. La
stessa liberta' **non** vale per il notch variabile: se i coefficienti dipendono
dal tempo, realizzazioni diverse della stessa funzione di trasferimento
"congelata" danno sistemi LTV *diversi*. Ed e' infatti per questo che sia
`ode_lpv_flex` sia il modello Simulink usano la **stessa** forma canonica di
controllabilita' per il notch.

**D: Se dovessi fidarti di un solo numero prodotto da questo file, quale
guarderesti con sospetto?**
R: `Vsafe = max(V,1)` e il commento associato: dichiara di proteggere da
`V(0) = 0`, ma nel dataset in repo `V` non scende mai sotto 410 m/s, quindi il
guard non entra mai in funzione e il commento e' semplicemente sbagliato. Non
cambia i risultati, ma e' il tipo di commento che porta un lettore a credere che
i primi secondi siano trattati in modo speciale, quando non lo sono. Da
segnalare anche il warm start della schedule, che il codice descrive ma non
esegue.
