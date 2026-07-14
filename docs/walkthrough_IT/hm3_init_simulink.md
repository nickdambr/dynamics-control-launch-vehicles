# HM3/init_simulink_hm3.m

## Ruolo del file nel progetto

`init_simulink_hm3.m` e' il **ponte fra il track script e il track Simulink** di
HM3. Non contiene ne' dinamica ne' progetto di controllo: ricalcola, chiamando
gli stessi moduli usati dai `main_task*.m` (`load_hw3_params`,
`build_plant_rigid`, `build_plant_full`, `build_tvc`, `build_notch_filter`,
`design_controller`, `load_wind_profile`), **ogni singolo parametro che i
blocchi di `models/hm3_closed_loop.slx` devono leggere**, e li spara nel base
workspace con `assignin`. Il modello Simulink non contiene nemmeno un numero
scritto a mano: nelle maschere dei blocchi ci sono solo *nomi di variabile*
(`A_full`, `tvc_num`, `Kp_th`, `wind_ts`, `Tstop`, `task`), che Simulink
risolve nel base workspace alla compilazione. Questo e' verificabile aprendo lo
`.slx` come zip: nel `system_root.xml` il blocco Gain `Controller_PD` ha
`Gain = [Kp_th, -Kp_th, -Kd_th, -Kp_z, -Kd_z]` e il blocco `From Workspace` ha
`VariableName = wind_ts`; i blocchi State-Space stanno invece nei file dei due
rami del variant -- `system_22.xml` (`Plant_Full`: `A = A_full`,
`B = [Bdelta_full, Bwind_full]`, `C = [C_meas_full; C_plot_full]`,
`D = zeros(7,2)`) e `system_7.xml` (`Plant_Rigid`, con gli omologhi `*_rigid`).

**Perche' esiste il track Simulink**, se gli script bastano gia' a produrre
Nichols, margini e risposte al gust? Due ragioni, entrambe difendibili
all'orale:

1. **Validazione incrociata.** Il track script vive interamente nel dominio
   delle funzioni di trasferimento: `assemble_loop` costruisce `T(s)` con
   `connect` e `simulate_gust_response` la integra con `lsim`, cioe' con una
   discretizzazione esatta del sistema LTI. Il track Simulink e' invece una
   integrazione numerica a passo variabile (`ode45`) di un diagramma a blocchi
   assemblato in modo indipendente: realizzazioni di stato diverse, ordine dei
   blocchi diverso, aritmetica diversa. Se le due risposte coincidono, gli
   errori possibili si riducono drasticamente: un segno sbagliato in una C, una
   riga della B scambiata, un guadagno con il segno invertito nel PD
   *romperebbero* la sovrapposizione. E' una prova indipendente, non una
   ripetizione.
2. **Formato di consegna industriale.** Il controllo d'assetto, nell'industria,
   si consegna come modello (Simulink -> autocodifica -> HIL), non come script
   di analisi. Il track Simulink e' la traduzione del progetto classico nel
   formato in cui verrebbe effettivamente integrato, con il vantaggio che il
   modello e' *scriptato dai dati* e non ri-digitato.

La filosofia dichiarata (README HM3 e `models/SIMULINK_GUIDE.md`) e'
**script-first**: gli script sono la fonte di verita', `init_simulink_hm3` e' il
compilatore che li traduce in parametri di blocco, `run_simulink_closed_loop` e'
il test di regressione. Il file e' chiamato da `run_simulink_closed_loop.m`
(riga 33) e a mano dall'utente prima di aprire il modello.

**Attenzione -- c'e' una discrepanza reale fra questo file e i `main_task2/3`,
verificata numericamente e documentata in fondo alla pagina** (sezione
"Discrepanze verificate"): i guadagni PD che questo file esporta *non sono*
quelli del progetto finale di Task 2/3.

---

## Firma e `arguments` block (righe 1-26)

```matlab
function S = init_simulink_hm3(task, o)
arguments
    task (1,1) {mustBeMember(task, [1 2 3])} = 2
    o.mu_alpha_scale (1,1) {mustBeNumeric, mustBeReal} = 1.0
    o.mu_c_scale     (1,1) {mustBeNumeric, mustBeReal} = 1.0
    o.severity ... = 'severe'
    o.profile  {mustBeTextScalar} = 'gust'
    o.push (1,1) logical = true
end
```

- Riga 1: la funzione restituisce una **struct `S`** con tutto cio' che ha
  calcolato. Il valore di ritorno esiste perche' `run_simulink_closed_loop` ha
  bisogno degli stessi dati per costruire la baseline analitica (non li rilegge
  dal base workspace, se li fa restituire): e' il modo pulito per garantire che
  script e modello vedano *esattamente* gli stessi numeri.
- Riga 20: `task` seleziona il **variant** del modello, non cosa viene
  calcolato. Il file calcola comunque **sia** le matrici rigide **sia** quelle
  flessibili, ed esporta entrambe. Il default e' `task = 2`.
- Righe 21-22: `mu_alpha_scale` e `mu_c_scale` sono i fattori moltiplicativi
  sui coefficienti mu_alpha = A_6 (momento aerodinamico destabilizzante) e
  mu_c = K_1 (efficacia del controllo) per i vertici +/- 30 % di Task 3. Vengono
  solo inoltrati a `load_hw3_params`.
- Righe 23-24: `severity` e `profile` scelgono il vento; `'strongwind'` fa
  girare il generatore Simulink del professore (`General/hw3-v3/strong_wind.slx`)
  invece del gust analitico 1-cosine.
- Riga 25: `push` (default `true`) e' l'unica opzione che ha un effetto
  collaterale: senza di essa il modello Simulink non troverebbe le variabili.

> **Possibile domanda d'esame** -- Perche' `task` non cambia cosa viene calcolato?
> *Risposta:* perche' il modello e' un unico `.slx` con un **Variant Subsystem**
> (`Plant_and_Actuator`) a due scelte (`Rigid_Ideal` con espressione di controllo
> `task == 1`, `Full_TVC_Notch` con `task ~= 1`). Il blocco non imposta
> `VariantActivationTime`, quindi vale il default *update diagram*: Simulink
> valuta l'espressione e **pota** la scelta inattiva prima della compilazione --
> le variabili del ramo non attivo non vengono nemmeno risolte, e in linea di
> principio potrebbero mancare. Esportare comunque **entrambi** i set (22
> variabili, righe 67-71) e' quindi una scelta di comodita', non un obbligo del
> compilatore: rende il base workspace valido per tutti e tre i task, cosi' che
> cambiare `task` e ri-simulare non richieda di rilanciare l'init.

---

## Parametri e controllore (righe 28-38)

```matlab
p  = load_hw3_params('mu_alpha_scale',o.mu_alpha_scale, ...
                     'mu_c_scale',o.mu_c_scale);
p0 = load_hw3_params();                    % nominal design point
K  = design_controller(build_plant_rigid(p0), [], 'verbose', false);
```

- Riga 32: `p` = parametri **perturbati** (con gli scale factor del corner). Da
  qui escono le matrici del plant che il modello simulera'.
- Riga 33: `p0` = parametri **nominali**. Servono per il controllore e per il
  notch, che devono restare *ignari* della perturbazione -- questo e' il cuore
  concettuale di uno studio di robustezza: il controllore vola con la sua
  taratura nominale, e' il veicolo a essere diverso da come lo si era immaginato.
- Riga 34: qui si progetta il PD. **Ed e' qui il problema.** La chiamata e'
  `design_controller(build_plant_rigid(p0), [], ...)`: plant **rigido**
  (4 stati, niente bending) e `Wact = []`, che dentro `design_controller`
  (riga 38 di quel file) diventa `tf(1)`, cioe' **attuatore ideale**. Questa e'
  esattamente la configurazione di **Task 1**. I `main_task2.m` (riga 153) e
  `main_task3.m` (righe 23-24) chiamano invece
  `design_controller(Gfull, Wtvc*Hn, ...)`, cioe' ri-tarano il PD **sul loop
  completo** (plant flessibile + TVC + ritardo + notch). I due tuner
  convergono a guadagni diversi:

      init_simulink_hm3 :  Kp_th = 1.7845   Kd_th = 0.4433   (Task 1)
      main_task2 / 3    :  Kp_th = 1.7318   Kd_th = 0.6867   (Task 2 re-tuned)

  Numeri ottenuti eseguendo i due tuner in MATLAB, non stimati.
- Righe 36-38: la struct `S` viene inizializzata e i quattro guadagni
  (`Kp_th`, `Kd_th`, `Kp_z`, `Kd_z`) copiati dentro. `Kp_z = Kd_z = -1e-3` sono i
  guadagni di drift, fissi per default in `design_controller` (riga 28-29): il
  segno negativo e' il *load relief* -- il lanciatore si lascia derivare un po'
  sottovento invece di inseguire z = 0, riducendo l'incidenza aerodinamica.

Il commento alle righe 29-31 dice: *"Controller and notch FROZEN at the nominal
point, as in main_task3: robustness = fixed nominal gains on the perturbed
plant"*. Lo **spirito** e' giusto (non si ri-tara sul corner). Il **valore
congelato**, pero', e' il controllore di Task 1, non quello di Task 2. Il
commento della testata (righe 16-17) -- *"Gains/filters are exactly those
designed by the scripts; the model only mirrors them"* -- e' quindi **falso per
task = 2 e task = 3**, ed e' precisamente il tipo di contraddizione
commento/codice che va segnalata.

> **Possibile domanda d'esame** -- Cosa cambia, in pratica, far volare il loop
> completo con i guadagni di Task 1 invece di quelli ri-tarati?
> *Risposta:* il loop resta **stabile**, ma con margini molto peggiori. Valutato
> sul loop di Task 2 (plant flessibile INS + TVC + Pade + notch profondo):
> con Kp = 1.784 / Kd = 0.443 -> Aero |GM| = 6.08 dB, **Rigid PM = 14.6 deg**,
> delay margin 98 ms; con Kp = 1.732 / Kd = 0.687 -> Aero |GM| = 6.00 dB,
> **Rigid PM = 30.0 deg**, delay margin 165 ms. E' esattamente il confronto
> "BEFORE / AFTER re-tuning" che `main_task2.m` stampa alle righe 146-161. Il
> derivativo piu' alto e' cio' che ricompra la fase mangiata da servo + ritardo
> + notch.

---

## Matrici del plant (righe 40-51)

```matlab
Gr = build_plant_rigid(p);
S.A_rigid  = Gr.A;
S.Bdelta_rigid = Gr.B(:,1);  S.Bwind_rigid = Gr.B(:,2);
S.C_meas_rigid = Gr.C(1:4,:);   % [theta_m thetadot_m z_m zdot_m]
S.C_plot_rigid = Gr.C(5:7,:);   % [theta z zdot]
```

- Righe 41-45: plant rigido, 4 stati `[z, zdot, theta, thetadot]`, 2 ingressi
  `[delta, alpha_w]`, 7 uscite.
- Riga 43: **la B viene spaccata in due colonne**. Questo e' l'unico vero
  "adattamento" a Simulink presente nel file. In MATLAB `Gr.B` e' 4x2 e `lsim`
  accetta un ingresso vettoriale a 2 componenti. In Simulink i due ingressi
  arrivano da posti fisici diversi -- `u_pd` (o l'uscita del TVC) da un lato,
  `alpha_w` dal blocco `From Workspace` dall'altro -- e vanno ricomposti in un
  segnale 2-wide con un **Mux** prima del blocco State-Space. Esportare le due
  colonne separatamente permette di scrivere nella maschera
  `B = [Bdelta_rigid, Bwind_rigid]`, che documenta l'ordine (delta = porta 1
  del Mux, alpha_w = porta 2) invece di nasconderlo.
- Righe 44-45: la C viene spaccata in **righe di misura** (1-4, cio' che il
  controllore legge) e **righe di plot** (5-7, gli stati veri). La distinzione
  e' fisica, non estetica: nel plant flessibile le righe 1-4 sono contaminate
  dal bending (Eq. 2, coefficienti `sigma_ins`, `phi_ins`), le righe 5-7 no.
  Nel modello questa separazione sopravvive nel `Demux` a 7 uscite: le prime 4
  tornano al Mux del controllore (chiudono l'anello), le ultime 3 vanno ai
  blocchi `To Workspace`.
- Righe 47-51: idem per il plant **completo**, 6 stati (aggiunge `eta`,
  `etadot` del primo modo di bending), con `meas = 'ins'`, cioe' la variante in
  cui il bending **entra nelle misure**. E' questa contaminazione -- non il modo
  di bending in se' -- che destabilizza il loop e motiva il notch.

> **Possibile domanda d'esame** -- Perche' non si esporta direttamente l'oggetto
> `ss` e si usa un blocco LTI System?
> *Risposta:* si potrebbe (esiste il blocco `LTI System` della Control System
> Toolbox). Il codice sceglie il blocco **State-Space** base, che usa solo
> Simulink core: rende il modello leggibile senza toolbox aggiuntive e, cosa
> piu' importante didatticamente, **espone le matrici** nel diagramma, dove si
> vede che B ha due colonne e C sette righe. Con un blocco LTI il modello
> sarebbe una scatola nera con dentro un oggetto MATLAB.

---

## Attuatore TVC + notch (righe 53-57)

```matlab
Wtvc = build_tvc(p,3);
[S.tvc_num, S.tvc_den] = tfdata(tf(Wtvc),'v');
Hx = build_notch_filter(p0.wBM, 0.002, 0.7, +1);
[S.notch_num, S.notch_den] = tfdata(Hx,'v');
```

- Riga 54: `build_tvc(p, 3)` = servo TVC di 2o ordine
  (omega_TVC = 70 rad/s, zeta_TVC = 0.7) **in serie con l'approssimante di Pade
  di ordine 3** del ritardo puro tau = 20 ms. Nota: `p` e' il set perturbato,
  ma `wTVC`, `zTVC`, `tau` non sono toccati dagli scale factor, quindi
  `build_tvc(p,3)` e `build_tvc(p0,3)` danno la stessa `tf` -- la scelta di `p`
  qui e' innocua.
- Riga 55: `tfdata(tf(Wtvc),'v')` estrae numeratore e denominatore come **vettori
  riga** (`'v'` = vector, non cell array). Il risultato:
  denominatore di grado 5 (2 poli del servo + 3 poli di Pade), numeratore di
  grado 3 (gli zeri di Pade), restituito con **due zeri iniziali di padding** per
  avere la stessa lunghezza del denominatore. Il blocco `Transfer Fcn` di
  Simulink accetta esattamente questo formato. Poli e zeri di Pade(3) per
  tau = 0.02 s (calcolati):

      poli:  -183.9 +/- 175.4i,  -232.2      (semipiano sinistro)
      zeri:  +183.9 +/- 175.4i,  +232.2      (semipiano DESTRO, speculari)

  La simmetria e' il punto: modulo unitario a tutte le frequenze (all-pass),
  solo fase. Pade **non e'** il ritardo vero, e' la sua approssimazione
  razionale -- l'unica cosa che `margin`, `pole`, `isstable`, `minreal` e il
  tuner `fminsearch` sanno maneggiare.
- Riga 56: `build_notch_filter(p0.wBM, 0.002, 0.7, +1)`, cioe'

      H_x(s) = (s^2 + 2*zeta_N*w_x*s + w_x^2) / (s^2 + 2*zeta_D*w_x*s + w_x^2)

  con w_x = omega_BM = 18.9 rad/s, zeta_N = 0.002, zeta_D = 0.7, `numSign = +1`
  (numeratore a fase minima -> notch simmetrico, non lead-lag). I coefficienti
  esportati sono (calcolati):

      notch_num = [1,  0.0756,  357.2]      (2*zeta_N*w_x = 0.0756)
      notch_den = [1, 26.46  ,  357.2]      (2*zeta_D*w_x = 26.46, w_x^2 = 357.2)

  Profondita' del notch = 20*log10(zeta_N/zeta_D) = **-50.9 dB** alla frequenza
  centrale. Da segnalare onestamente: `build_notch_filter` documenta nel suo
  proprio header la *guideline* `zN` in 0.1-0.3 -- qui si usa 0.002, cioe' **due
  ordini di grandezza fuori dalla linea guida della traccia**. E' voluto (e'
  il "deep notch" scelto dal trade di `main_task2.m`), ma va saputo difendere:
  il prezzo e' che il notch e' **strettissimo** -- la coppia di zeri e' larga
  2*zeta_N*w_x = 0.076 rad/s su una centrale di 18.9 rad/s -- ed e' proprio questa
  strettezza che obbliga a conoscere omega_BM con precisione (Step D di
  `main_task2` mostra che +5 % di errore su omega_BM destabilizza).
- Riga 56, dettaglio importante: il notch e' centrato su **`p0.wBM`** (nominale)
  e non su `p.wBM`. Concettualmente e' la scelta giusta (il filtro non deve
  conoscere la perturbazione). Numericamente qui non cambia nulla, perche' gli
  scale factor toccano solo `A6` e `K1`, mai `wBM`.

> **Possibile domanda d'esame** -- Perche' precalcolare `tvc_num`/`tvc_den` in
> MATLAB invece di scrivere l'espressione nella maschera del blocco?
> *Risposta:* perche' i coefficienti non sono scrivibili a mano. Sono il
> risultato del prodotto polinomiale servo x Pade(3): sei coefficienti a
> denominatore che nessuno vuole ricalcolare a mano ogni volta che cambia tau o
> l'ordine di Pade. Metterli in maschera come espressione MATLAB (`pade(0.02,3)`
> dentro il blocco) sarebbe possibile ma **duplicherebbe la definizione**: da
> quel momento in poi il modello e lo script potrebbero divergere senza che
> nessuno se ne accorga. Con il precalcolo esiste **una sola** definizione
> (`build_tvc.m`) e il modello e' un consumatore passivo.

---

## Vento come `timeseries` (righe 59-64)

```matlab
w = load_wind_profile(p,'severity',o.severity,'profile',o.profile);
S.wind_ts = timeseries(w.alphaw(:), w.t(:), 'Name','alpha_w');
S.Tstop = w.t(end);
S.task = task;
```

- Riga 60: il vento e' generato **una volta sola**, dallo stesso
  `load_wind_profile` che usano gli script. Il gust di default e' un 1-cosine
  di 3 s con V_g preso dalla dispersione `drywind.mat` alla quota di volo
  (V_g = 6.38 m/s -> picco alpha_w = 0.39 deg).
- Riga 61: viene impacchettato come oggetto `timeseries` perche' e' il formato
  che il blocco **From Workspace** legge nativamente (`VariableName = wind_ts`
  nel `.slx`; `OutputAfterFinalValue = 'Holding final value'`). Dentro il
  blocco Simulink **interpola linearmente** fra i campioni: la griglia del vento
  ha dt = 0.005 s, quella di `ode45` e' variabile, quindi l'interpolazione c'e'
  ed e' una delle sorgenti (piccole) di scostamento fra i due track.
- Riga 62: `Tstop` = fine del profilo di vento = **12 s** con i default
  (`load_wind_profile` ha `Tend = 12` di default). Attenzione: `main_task2.m` e
  `main_task3.m` chiamano invece `load_wind_profile(p, Tend=80)` -- 80 s, per
  vedere il modo lento di drift (tau ~ 20 s). **Le figure Simulink coprono
  quindi un orizzonte diverso (12 s) da quelle degli script (80 s).** Per il
  corner peggiore V3 i picchi cadono comunque entro i primi 12 s (verificato:
  stessi picchi a 12 e a 80 s), quindi non e' un errore, ma e' una differenza da
  conoscere se si confrontano le figure.
- Riga 64: `S.task = task` -- questa e' la variabile che l'espressione di
  controllo del Variant Subsystem (`task == 1` / `task ~= 1`) legge.

---

## Push nel base workspace (righe 66-75)

```matlab
if o.push
    fn = fieldnames(S);
    for i = 1:numel(fn)
        assignin('base', fn{i}, S.(fn{i}));
    end
    fprintf('init_simulink_hm3: pushed %d variables to base ...');
end
```

- Righe 67-71: **22 variabili** finiscono nel base workspace (`p`, i 4 guadagni,
  5 matrici rigide, 5 matrici full, 4 vettori di coefficienti, `wind_ts`,
  `Tstop`, `task`). Il loop su `fieldnames` fa si' che aggiungere un campo a `S`
  lo esporti automaticamente: nessuna lista da tenere aggiornata a mano.
- **Perche' il base workspace e non il model workspace?** Simulink risolve i nomi
  nelle maschere risalendo la gerarchia: mask workspace -> model workspace ->
  data dictionary -> **base workspace**. Il modello HM3 non ha model workspace
  ne' data dictionary, quindi il base workspace e' l'unico posto dove le
  variabili possono stare. Il prezzo e' che il base workspace e' **globale**:
  se qualcuno ha una variabile `task` o `p` gia' definita, viene sovrascritta.
  Una versione piu' industriale userebbe `Simulink.SimulationInput.setVariable`
  (che infatti `load_wind_profile` usa alle righe 92-96 per il modello del
  professore) -- qui non lo si fa, e vale la pena saperlo dire.
- Riga 72: il `fprintf` di conferma stampa quante variabili e con quali scale
  factor. E' l'unica traccia visibile che il push e' avvenuto.

---

## Discrepanze verificate fra i due track

Tutto quanto segue e' stato **eseguito**, non dedotto:

1. **I guadagni PD del modello Simulink, per task = 2 e task = 3, non sono
   quelli del progetto consegnato.** `init_simulink_hm3` (riga 34) tara il PD sul
   plant rigido con attuatore ideale (Kp = 1.7845, Kd = 0.4433, i guadagni di
   Task 1), mentre `main_task2.m` (riga 153) e `main_task3.m` (riga 23) lo
   ri-tarano sul loop completo (Kp = 1.7318, Kd = 0.6867). Il loop risultante e'
   stabile ma con Rigid PM = 14.6 deg invece di 30 deg e delay margin 98 ms
   invece di 165 ms. Il commento in testata del file ("Gains/filters are exactly
   those designed by the scripts") **contraddice il codice**.
2. **Conseguenza sulla validazione:** la sovrapposizione script-vs-Simulink
   *continua a funzionare*, perche' `run_simulink_closed_loop` costruisce la sua
   baseline analitica con la **stessa** `K` restituita da `init_simulink_hm3`
   (riga 43 di quel file). I due track coincidono -- ma coincidono su un
   controllore che non e' quello del report. La validazione e' internamente
   coerente e **esternamente sfasata** rispetto a Task 2/3.
3. **I numeri di Task 3 in `models/SIMULINK_GUIDE.md` (righe 222-224) non
   tornano.** La guida dice: corner V3, picco theta ~ 0.35 deg e picco z ~ 4.2 m,
   *"matching main_task3.m"*. In realta' `main_task3.m` (guadagni ri-tarati)
   da' **0.878 deg / 4.92 m** (come riporta il README), e il modello Simulink
   cosi' com'e' inizializzato (guadagni di Task 1) da' **1.402 deg / 4.88 m**.
   La guida non coincide con nessuno dei due: e' testo stantio.
4. **Nessun conflitto sul Pade.** Verificando lo `.slx` (`system_22.xml`), il
   ramo `Full_TVC_Notch` contiene **solo due blocchi `Transfer Fcn`**
   (`notch_num/notch_den` e `tvc_num/tvc_den`): **non c'e' nessun blocco
   Transport Delay**. Il ritardo e' approssimato con Pade in *entrambi* i track,
   ed e' esattamente la stessa `tf`. Su questo punto i due track sono coerenti.

Come si aggiusterebbe la (1), se lo si volesse: sostituire la riga 34 con una
progettazione condizionata a `task` -- per `task == 1` la chiamata attuale, per
`task >= 2` `design_controller(build_plant_full(p0,'ins'), build_tvc(p0,3)*Hx,
'w_flex',0.6*p0.wBM, 'w_flex_hi',1.5*p0.wBM, 'w_bending',p0.wBM)`. Non lo faccio
qui perche' la consegna e' documentare, non modificare.

---

## Possibili domande d'esame

**D: Perche' costruire un modello Simulink se gli script producono gia' tutti i
risultati richiesti dalla traccia?**
R: Per due motivi. Il primo e' la **validazione incrociata**: due
implementazioni indipendenti dello stesso anello chiuso (funzioni di
trasferimento + `lsim` da un lato, diagramma a blocchi + `ode45` dall'altro) che
devono dare la stessa risposta temporale. Un errore di segno, una C sbagliata,
un guadagno invertito farebbero divergere le due curve: la coincidenza e' quindi
una prova, non una ripetizione. Il secondo e' che il **controllo d'assetto si
consegna come modello**, non come script: Simulink e' il formato da cui si fa
autocodifica e HIL, e il modello e' anche il documento che si mostra a chi
integra.

**D: Perche' precalcolare ogni parametro in MATLAB invece di mettere espressioni
nelle maschere dei blocchi?**
R: Perche' i parametri non sono costanti: i guadagni PD escono da un
`fminsearch` su margini classificati per banda (`design_controller`), i
coefficienti del TVC dal prodotto servo x Pade(3), le matrici del plant da
`build_plant_full`. Scrivere queste espressioni in maschera significherebbe
**duplicare la definizione** e permettere a modello e script di divergere in
silenzio; inoltre le espressioni in maschera vengono rivalutate a ogni compile
del modello, il che vorrebbe dire rieseguire il tuner ogni volta. Con il
precalcolo esiste una sola fonte di verita' (gli `.m`) e il modello e' un
consumatore passivo di nomi di variabile.

**D: Come fa un unico `.slx` a coprire tre task con plant e attuatori diversi?**
R: Con un **Variant Subsystem** (`Plant_and_Actuator`) a due scelte:
`Rigid_Ideal` (espressione di controllo `task == 1`: plant rigido a 4 stati,
attuatore ideale, delta = u_pd) e `Full_TVC_Notch` (`task ~= 1`: notch -> TVC con
Pade -> plant flessibile a 6 stati). La variabile `task`, esportata da
`init_simulink_hm3` (riga 64), decide quale ramo e' attivo a inizio simulazione.
Task 2 e Task 3 usano lo *stesso* ramo: differiscono solo per i valori dentro
`A_full`/`Bdelta_full`, che cambiano quando si passano `mu_alpha_scale` e
`mu_c_scale`. Questa e' proprio la definizione operativa di studio di
robustezza: **cambia il veicolo, non il controllore**.

**D: Perche' il notch e' centrato su `p0.wBM` e i guadagni sono calcolati su
`p0`, mentre le matrici del plant usano `p`?**
R: Perche' e' la differenza fra **quello che il controllore sa** e **quello che
il veicolo e'**. Il controllore vola con la sua taratura nominale: non conosce
il vero A_6, non conosce il vero K_1. Il plant, invece, e' perturbato. Se
ri-tarassimo il PD su ogni corner, non staremmo piu' misurando la robustezza
del progetto ma la bravura del tuner: e infatti il commento alle righe 29-31
avverte che il re-tuning sul corner destabilizza il loop con il bending. Nel
caso specifico gli scale factor toccano solo A_6 e K_1, quindi `p.wBM == p0.wBM`
e il notch e' identico nei due casi -- ma la struttura del codice e' quella
giusta e resta corretta se domani si perturbasse anche omega_BM (cosa che
`main_montecarlo.m` fa).

**D: Il modello Simulink implementa davvero il controllore che hai consegnato
nel report?**
R: **Per Task 1 si', per Task 2 e 3 no** -- ed e' onesto dirlo. La riga 34 di
`init_simulink_hm3.m` tara il PD sul plant *rigido* con attuatore *ideale*
(Kp = 1.784, Kd = 0.443), che e' il progetto di Task 1; i `main_task2/3`
ri-tarano invece sul loop completo (Kp = 1.732, Kd = 0.687). Il loop simulato
resta stabile, ma con margine di fase rigido 14.6 deg invece di 30 deg. La
validazione script-vs-Simulink resta valida (entrambi i track usano la stessa
`K`), ma valida la configurazione "before re-tuning". La correzione e' di una
riga: condizionare il progetto del controllore al valore di `task`.

**D: Cosa succede se si dimentica di rilanciare `init_simulink_hm3` dopo aver
cambiato un guadagno negli script?**
R: Il modello simula tranquillamente **con i vecchi valori**, perche' i blocchi
leggono le variabili gia' presenti nel base workspace e non hanno modo di
sapere che l'`.m` e' cambiato. Non c'e' nessun meccanismo di invalidazione:
questa e' la fragilita' strutturale della strategia "push nel base workspace", e
il motivo per cui l'header del file (riga 17) raccomanda "*Re-run whenever a
gain changes*". Una versione robusta userebbe `Simulink.SimulationInput` con
`setVariable`, cosi' che i parametri viaggino con la singola simulazione invece
di vivere in uno stato globale.
