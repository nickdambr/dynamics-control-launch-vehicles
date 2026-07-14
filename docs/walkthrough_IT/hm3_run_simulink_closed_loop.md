# HM3/run_simulink_closed_loop.m

## Ruolo del file nel progetto

`run_simulink_closed_loop.m` e' il **driver di validazione** del track Simulink.
Fa tre cose in sequenza: (1) chiama `init_simulink_hm3` per popolare il base
workspace con guadagni, matrici, filtri e vento; (2) ricostruisce con gli stessi
dati la **baseline analitica** degli script (`assemble_loop` +
`simulate_gust_response`); (3) simula `models/hm3_closed_loop.slx` con `sim()`,
estrae i segnali loggati e **sovrappone** le due risposte in un'unica figura,
`figures/task<N>_simulink_vs_script.png`.

E' l'unico punto della repo in cui i due track si toccano. La logica e' quella
di un **test di equivalenza**: due implementazioni indipendenti dello stesso
anello chiuso -- una in funzioni di trasferimento integrata con `lsim`, l'altra
come diagramma a blocchi integrato con `ode45` -- devono produrre la stessa
storia temporale a fronte dello **stesso identico vento**. La riga 46 e' la
chiave del rigore: il vento non viene rigenerato, viene **rimesso in gioco**
leggendolo dal `timeseries` che il modello stessa usa
(`S.wind_ts`), cosi' che l'unica differenza fra i due track sia numerica e non
di ingresso.

Ordine di chiamata: `run_simulink_closed_loop(1|2|3, ...)` a mano dall'utente.
Dipende da `init_simulink_hm3`, `build_plant_rigid`, `build_plant_full`,
`build_tvc`, `build_notch_filter`, `assemble_loop`, `simulate_gust_response`, e
dal file `models/hm3_closed_loop.slx`, che **non genera** (se manca, stampa un
messaggio e ritorna `[]`).

Due cose vanno dette subito, perche' sono il punto debole del file:
**non calcola nessuna metrica di scostamento** (il confronto e' puramente
visivo: produce una figura, non un numero) e **non esiste nessun test** in
`HM3/tests/` che copra il track Simulink (verificato con grep: nessun match su
`simulink` o `hm3_closed_loop`). I numeri di accordo citati nel report sono
stati calcolati fuori da questa funzione.

---

## Firma, `arguments` e guardia sul modello (righe 1-29)

```matlab
function out = run_simulink_closed_loop(task, o)
here  = fileparts(mfilename('fullpath'));
model = 'hm3_closed_loop';
mdlfile = fullfile(here,'models',[model '.slx']);
if ~isfile(mdlfile)
    fprintf(['[run_simulink_closed_loop] Model not found:\n  %s\n' ...
             'Build it first by following models/SIMULINK_GUIDE.md ...']);
    out = [];  return;
end
```

- Righe 13-19: l'`arguments` block replica quello di `init_simulink_hm3` meno
  l'opzione `push` (che qui deve restare `true`, altrimenti il modello non
  troverebbe le variabili). Le opzioni vengono poi inoltrate in blocco alla riga
  32 con `namedargs2cell(o)`, che trasforma la struct name-value in una cell
  array `{'mu_alpha_scale',1.3,'mu_c_scale',0.7,...}` da splattare nella
  chiamata: e' il modo canonico per fare *pass-through* di name-value senza
  ricopiarli uno a uno.
- Righe 21-23: `mfilename('fullpath')` rende la funzione **indipendente dalla
  current directory**: tutti i percorsi (modello, figure) sono ancorati alla
  cartella `HM3/`, non a dove l'utente si trova.
- Righe 24-29: **guardia**. Lo `.slx` e' un binario che gli agenti/script non
  generano: se manca, la funzione non esplode ma stampa l'istruzione ("costruisci
  il modello seguendo `models/SIMULINK_GUIDE.md`") e ritorna `[]`. Il commento in
  testata (righe 9-10) dice *"The .slx is built by hand ..., not auto-generated"* --
  ed e' coerente con lo stato della repo, dove il `.slx` esiste ed e' committato.

---

## Init del workspace e baseline degli script (righe 31-47)

```matlab
optArgs = namedargs2cell(o);
S = init_simulink_hm3(task, optArgs{:});
p = S.p;
if task == 1
    G = build_plant_rigid(p);  Wact = [];
else
    G = build_plant_full(p,'ins');
    Wact = build_tvc(p,3) * build_notch_filter(p.wBM,0.002,0.7,+1);
end
K.Kp_th=S.Kp_th; ... 
[~,T] = assemble_loop(G,K,Wact);
w = struct('t', S.wind_ts.Time(:), 'alphaw', squeeze(S.wind_ts.Data), 'V', p.V);
rs = simulate_gust_response(T,w);
```

- Riga 33: `init_simulink_hm3` viene chiamato **prima** di tutto. E' lui che
  spara le 22 variabili nel base workspace e restituisce `S`. Il fatto che la
  funzione poi legga i guadagni **da `S`** (riga 43) e non li ricalcoli e'
  proprio cio' che garantisce che baseline e modello condividano lo stesso `K`.
- Righe 37-42: la baseline analitica viene ricostruita **replicando a mano** cio'
  che il modello contiene. Per `task == 1`: plant rigido, `Wact = []` che
  `assemble_loop` interpreta come `tf(1)` (attuatore ideale) -- nel modello e'
  il ramo `Rigid_Ideal`, dove `delta` e' semplicemente `u_pd` derivato dalla
  linea. Per `task >= 2`: plant flessibile `'ins'` (bending nelle misure) e
  `Wact = TVC * notch`.
- Riga 41, **dettaglio sull'ordine**: nello script la catena e'
  `build_tvc(...) * build_notch_filter(...)`, cioe' TVC *poi* notch; nel modello
  (verificato nell'XML dello `.slx`, `system_22.xml`) le linee sono
  `u_pd -> Notch_Hx -> TVC -> plant`, cioe' notch *poi* TVC. **Non e' un bug**:
  sono due blocchi SISO LTI in serie e il prodotto e' commutativo, la funzione di
  trasferimento risultante e' identica. Cambiano solo le realizzazioni di stato
  interne (e quindi, marginalmente, gli errori di arrotondamento) -- che e'
  esattamente il tipo di indipendenza che si vuole in una validazione incrociata.
- Riga 41, **insidia latente**: qui il notch e' centrato su `p.wBM`
  (parametri **perturbati**), mentre `init_simulink_hm3` (riga 56) lo centra su
  `p0.wBM` (**nominale**). Oggi non fa differenza, perche' `mu_alpha_scale` e
  `mu_c_scale` toccano solo A_6 e K_1 e mai omega_BM. Ma se domani si perturbasse
  omega_BM (cosa che `main_montecarlo.m` fa), i due track userebbero **notch
  diversi** e la sovrapposizione fallirebbe -- non per un errore di Simulink, ma
  per questa asimmetria. Vale la pena saperlo.
- Riga 44: `assemble_loop(G,K,Wact)` restituisce `T`, il closed loop
  `{alpha_w, theta_ref} -> {theta, z, zdot, delta}`.
- Righe 45-46: **il punto piu' importante del file.** Il vento non viene
  ri-generato con `load_wind_profile`: viene ricostruito **dal `timeseries` che
  il modello legge**, con `S.wind_ts.Time` e `S.wind_ts.Data`. Se lo si
  rigenerasse, per il gust analitico si otterrebbe lo stesso risultato -- ma per
  `profile = 'strongwind'` (che gira il generatore stocastico del professore) si
  potrebbe ottenere una realizzazione diversa, e la sovrapposizione perderebbe
  ogni significato. Cosi' com'e', **e' garantito per costruzione** che i due
  track vedano lo stesso alpha_w(t) campione per campione.
- Riga 47: `simulate_gust_response(T,w)` fa `lsim(T, [alphaw, 0], t)` sulla
  griglia uniforme del vento (dt = 0.005 s).

> **Possibile domanda d'esame** -- Se `lsim` e' esatto per un sistema LTI, perche'
> le due curve non coincidono a machine precision?
> *Risposta:* `lsim` e' esatto **dato un modello di ricostruzione dell'ingresso**
> (mantiene l'ingresso costante o lineare fra due campioni della griglia) e
> discretizza esattamente la dinamica su quel passo. Simulink invece integra con
> `ode45` a passo variabile e ricostruisce l'ingresso interpolando linearmente il
> `timeseries` del blocco `From Workspace` sui suoi istanti. Le due sorgenti di
> differenza sono quindi: (a) la **ricostruzione dell'ingresso** su griglie
> diverse, e (b) l'**errore di integrazione** di `ode45` (controllato da
> `RelTol = 1e-6`). Entrambe sono O(1e-6) relative, e infatti e' quello che si
> misura.

---

## Simulazione del modello (righe 49-57)

```matlab
addpath(fullfile(here,'models'));
so = sim(model, 'StopTime', num2str(S.Tstop));
get_ts = @(nm) get_logged_signal(so, nm);
[sl.t, sl.theta] = get_ts('theta_sl');
[~,    sl.z]     = get_ts('z_sl');
[~,    sl.zdot]  = get_ts('zdot_sl');
[~,    sl.delta] = get_ts('delta_sl');
```

- Riga 50: `addpath` sulla cartella `models/` perche' `sim('hm3_closed_loop')`
  cerca il modello **sul path**, non nella cartella corrente. Effetto
  collaterale: il path resta sporcato dopo l'uscita (nessun `onCleanup` che lo
  rimuova). Difetto minore ma reale.
- Riga 51: `sim(model, 'StopTime', num2str(S.Tstop))` -- forma comando classica,
  con override esplicito dello stop time. Nota che il modello ha gia'
  `StopTime = 'Tstop'` nella sua config (letto dallo `.slx`), quindi l'override
  e' **ridondante ma esplicito**: rende il codice leggibile senza dover aprire
  il modello. `Tstop` vale 12 s con il vento di default.
- Riga 51, cosa restituisce: nel modello `ReturnWorkspaceOutputs = 'on'` (letto
  in `configSet0.xml`), quindi `so` e' un oggetto `Simulink.SimulationOutput` e
  le variabili dei blocchi `To Workspace` (`theta_sl`, `z_sl`, `zdot_sl`,
  `delta_sl`, salvate in formato `Timeseries`) ne diventano **proprieta'**,
  invece di essere sparate nel base workspace.
- Righe 53-57: estrazione dei quattro segnali con un handle anonimo. **Il
  contratto e' sui nomi delle variabili, non su quelli dei blocchi**: nel root
  del modello i quattro `To Workspace` si chiamano `log_theta`, `log_z`,
  `log_zdot`, `log_delta` (verificato nell'XML), ma cio' che questa funzione
  cerca sono i loro `VariableName`: `theta_sl`, `z_sl`, `zdot_sl`, `delta_sl`.
  Rinominare un blocco e' quindi innocuo; cambiare un `VariableName` rompe la
  validazione -- e, come si vede sotto, la rompe *in silenzio*.

### Configurazione del solver (dalla config del modello, non dallo script)

Estratta direttamente da `hm3_closed_loop.slx` (`simulink/configSet0.xml`):

    SolverName = ode45        (Dormand-Prince, esplicito, PASSO VARIABILE)
    RelTol     = 1e-6
    AbsTol     = auto
    MaxStep    = auto,   FixedStep = auto
    StopTime   = Tstop   (variabile di workspace)

Vale la pena **derivare** perche' questa e' la scelta giusta, perche' l'istinto
("c'e' un modo poco smorzato -> sistema stiff -> usa `ode15s`") porta fuori
strada.

Poli del loop chiuso di Task 2 (calcolati): il piu' lento ha |lambda| = 0.233
rad/s (il modo di drift laterale), il piu' veloce |lambda| = 256 rad/s (i poli di
Pade a -183.9 +/- 175.4i e -232.2). **Rapporto di stiffness ~ 1.1e3**: modesto.
La condizione di stabilita' di `ode45` impone h ~< 2.8/256 ~ 11 ms; con
`RelTol = 1e-6` l'accuratezza impone comunque passi molto piu' piccoli, quindi
la stabilita' non e' il vincolo. Su un orizzonte di 12 s il costo e' irrisorio.

Il modo di bending non e' "stiff": **e' un problema di accuratezza, non di
stabilita'**. Con omega_BM = 18.9 rad/s e zeta_BM = 0.005 si ha
Q = 1/(2*zeta) = 100, periodo 0.332 s, tempo di decadimento
1/(zeta*omega_BM) = 10.6 s. Su 12 s di simulazione il modo compie ~36 periodi e
**non fa in tempo a smorzarsi**. Un solver che introduce smorzamento numerico
(le formule BDF di `ode15s`, o `ode23tb`) attenuerebbe artificialmente proprio
il modo che il notch deve tenere sotto controllo, e la sovrapposizione con la
baseline `lsim` peggiorerebbe di ordini di grandezza. `ode45` esplicito con
tolleranza stretta e' la scelta conservativa: **paga passi piccoli, ma non
falsifica l'ampiezza della risonanza**.

### Algebraic loop: non c'e', e si puo' dimostrare

Percorso di retroazione nel modello: `Controller_PD` (Gain -> feedthrough
diretto) -> `Plant_and_Actuator` -> `Demux` -> `Meas_Mux` -> `Controller_PD`.
Il loop si spezza dentro il variant:

- ramo `Full_TVC_Notch`: `Notch_Hx` e' un `Transfer Fcn` **biproprio** (num e den
  entrambi di grado 2) e quindi **ha feedthrough**; ma subito dopo c'e' il
  blocco `TVC`, con numeratore di grado 3 e denominatore di grado 5, cioe'
  **strettamente proprio** (grado relativo 2, nessun feedthrough). Il plant
  State-Space ha `D = zeros(7,2)`. Il feedthrough e' interrotto **due volte**.
- ramo `Rigid_Ideal`: `u_pd -> Mux -> State-Space` con `D = zeros(7,2)`: nessun
  feedthrough. (Il ramo `delta = u_pd` e' feedthrough, ma va solo a un
  `To Workspace`, non chiude nessun anello.)

Quindi Simulink non ha nulla da risolvere iterativamente: nessun warning di
algebraic loop, nessun solutore algebrico. **Ma e' fragile**: il notch e'
biproprio, e basterebbe che il plant avesse un `D` non nullo (per esempio se la
misura di theta contenesse un termine proporzionale a delta) perche' l'anello
diventasse algebrico e Simulink dovesse risolvere il punto fisso a ogni passo.
Il fatto che `D = 0` non e' un caso: e' la struttura fisica del problema (la
deflessione dell'ugello genera un momento, quindi accelerazione, quindi arriva
sulle uscite solo dopo due integrazioni).

---

## Overlay e figura (righe 59-83)

```matlab
nexttile; plot(rs.t,rs.theta*180/pi,'b-', ...
               sl.t,sl.theta*180/pi,'r--','LineWidth',1.3);
```

- Righe 60-68: tre riquadri -- theta [deg], z [m], delta [deg] -- con lo script in
  **blu continuo** e Simulink in **rosso tratteggiato**. Convenzione grafica
  volutamente esplicita: se le due curve sono distinguibili a occhio, qualcosa
  non va.
- Righe 72-76: `theme(f,'light')` con fallback su `f.Color = 'w'` -- forza il
  tema chiaro anche se il desktop MATLAB e' in dark mode (altrimenti le PNG
  esportate uscirebbero con sfondo scuro nel report).
- Righe 77-81: la figura di default (`profile = 'gust'`) si chiama
  `task<N>_simulink_vs_script.png`; qualunque altro profilo (es. `strongwind`)
  aggiunge un suffisso, cosi' che la figura del gust nominale non venga
  sovrascritta. Nella repo esistono infatti
  `task1_/task2_/task3_simulink_vs_script.png` e
  `task2_simulink_vs_script_strongwind.png`.
- Riga 83: `out = struct('script',rs,'simulink',sl)` -- le due strutture vengono
  restituite, cosi' che il chiamante *possa* calcolarsi lo scostamento. **Ma la
  funzione non lo calcola.**

### Che tolleranza si accetta, e cosa significherebbe uno scostamento

Il codice **non impone nessuna soglia**. La quantificazione esiste solo nel
report (`HM3/report/chapters/SimulinkReproduction.tex`, tabella alle righe
106-120), calcolata fuori da questa funzione:

    caso                    max|d theta| [rad]   max|d z| [m]   max|d delta| [rad]
    Task 1 - gust                3.2e-6            8.8e-4            5.4e-6
    Task 2 - gust                1.7e-7            1.0e-5            6.2e-7
    Task 2 - strong wind         8.7e-7            7.7e-5            1.8e-5
    Task 3 - corner V3           3.7e-7            2.3e-5            9.4e-7

Da leggere in **termini relativi**, con il normalizzatore giusto: il picco a cui
rapportare l'errore e' quello del **run Simulink**, che gira con i guadagni
pre-retune (punto 2 piu' sotto), cioe' ~0.29 deg = 5e-3 rad -- **non** gli
0.231 deg del progetto ri-tarato di `main_task2`. Quindi 1.7e-7 rad e' un errore
relativo ~3e-5, coerente con `RelTol = 1e-6` accumulato su 12 s piu' la
differente ricostruzione dell'ingresso. Su z (picco ~2.2 m) l'errore
sub-millimetrico e' dello stesso ordine relativo.

Cosa significherebbe **uno scostamento vero**:

- **1e-3 relativo o piu', ma con la stessa forma d'onda** -> problema numerico:
  tolleranza troppo lasca, `MaxStep` troppo grande, o interpolazione del vento su
  una griglia troppo rada. Si stringe `RelTol` e si ricontrolla.
- **Forme d'onda diverse, o segno opposto, o ampiezza doppia** -> errore
  **strutturale**: una riga di C scambiata, la B con le colonne invertite (delta
  e alpha_w scambiati nel Mux!), un segno sbagliato nel vettore dei guadagni, il
  demux ricablato nell'ordine sbagliato. E' *esattamente* la classe di errori che
  questa validazione esiste per pescare, e la ragione per cui il track Simulink
  non e' un esercizio di stile.
- **Divergenza a fine transitorio, con oscillazione a ~19 rad/s** -> il modo di
  bending sta venendo smorzato numericamente da una parte e non dall'altra:
  solver o tolleranza sbagliati.

> **Possibile domanda d'esame** -- La sovrapposizione perfetta dimostra che il
> progetto e' corretto?
> *Risposta:* no. Dimostra che le **due implementazioni** dello stesso progetto
> coincidono. Se il modello matematico e' sbagliato (segno di A_6 invertito,
> coefficiente della traccia letto male), entrambi i track sbagliano allo stesso
> modo e la sovrapposizione e' perfetta lo stesso. La validazione incrociata
> pesca gli errori di *trascrizione*, non quelli di *modellazione*: per quelli
> servono altri controlli (il polo instabile a +sqrt(A_6) = +1.84 rad/s che
> corrisponde alla fisica, i margini che tornano con i valori attesi, la
> risposta al gust dell'ordine di grandezza giusto).

---

## `get_logged_signal` (righe 86-102)

```matlab
function [t,y] = get_logged_signal(so, name)
t = []; y = [];
try
    if isprop(so,'logsout') && ~isempty(so.logsout) && ...
            any(strcmp(so.logsout.getElementNames, name))
        e = so.logsout.getElement(name);
        t = e.Values.Time;  y = e.Values.Data;
    elseif isprop(so,name) || isfield(so,name)
        v = so.(name);  t = v.Time;  y = v.Data;
    end
catch
    warning('run_simulink_closed_loop:signal', ...
            'Could not read signal "%s".',name);
end
```

- Riga 86: helper locale. Serve perche' in Simulink **ci sono due modi** di
  tirar fuori un segnale da una simulazione, e la funzione li prova entrambi.
- Righe 92-95: prima strada -- `logsout`, il *signal logging dataset*. Si popola
  quando una **linea** del modello e' marcata con `DataLogging = 'on'` (e' quello
  che fa `load_wind_profile` alle righe 87-89 per il modello del professore).
  Nel modello HM3 la config ha `SignalLogging = 'on'`, ma **nessuna linea e'
  marcata**, quindi `logsout` e' vuoto e questo ramo non scatta mai.
- Righe 96-97: seconda strada -- i blocchi `To Workspace`. Con
  `ReturnWorkspaceOutputs = 'on'`, le loro variabili diventano proprieta'
  dell'oggetto `SimulationOutput`, quindi `so.theta_sl` e' un `timeseries` con
  `.Time` e `.Data`. **E' il ramo effettivamente usato.** (Nota: `isfield` su un
  oggetto `SimulationOutput` ritorna sempre `false`; l'`isprop` e' cio' che
  regge.)
- Righe 90 e 99-101: **degradazione silenziosa.** Se il segnale non esiste sotto
  nessuna delle due forme, la funzione non solleva errore: ritorna `t = []`,
  `y = []`. Il `plot` a valle disegnera' semplicemente la curva blu e nessuna
  curva rossa, e la figura verra' salvata come se nulla fosse. Il `catch` avvisa
  solo se l'accesso *lancia* un errore, non se il nome semplicemente non c'e'.
  Questo e' un difetto vero: **il fallimento della validazione puo' passare
  inosservato**. Un `assert(~isempty(sl.theta), ...)` dopo la riga 57 costerebbe
  una riga.

---

## Discrepanze e limiti verificati

Tutto quanto segue e' stato ispezionato sul codice e sullo `.slx` (aperto come
archivio zip), non ipotizzato:

1. **Il ritardo puro e' Pade in ENTRAMBI i track.** Il sospetto ricorrente
   ("script con Pade, Simulink con blocco `Transport Delay`") **non si verifica**:
   nel ramo `Full_TVC_Notch` del modello ci sono solo due blocchi `Transfer Fcn`,
   `notch_num/notch_den` e `tvc_num/tvc_den`, e quest'ultimo *contiene gia'* il
   Pade di ordine 3 calcolato da `build_tvc`. Nessun blocco `Transport Delay`
   esiste nel modello. I due track condividono percio' la stessa
   approssimazione. Se cosi' non fosse, la differenza sarebbe sostanziale:
   `exp(-tau*s)` e' trascendente (fase esattamente -omega*tau a ogni omega, modulo
   1); Pade(3) e' una razionale all-pass che riproduce quella fase solo in banda
   (a omega_BM: omega*tau = 0.378 rad = 21.7 deg, dove Pade(3) e' praticamente
   esatto; l'errore cresce a decine di rad/s piu' in alto). Inoltre il blocco
   `Transport Delay` in Simulink e' implementato con un **buffer di stati
   passati interpolati**, non con un'equazione di stato: introduce un errore di
   interpolazione proprio, forza `ode45` a passi ridotti vicino alle
   discontinuita', e -- cosa decisiva -- **non e' linearizzabile**, quindi il
   progetto in frequenza (`margin`, Nichols, `allmargin`, il tuner `fminsearch`
   di `design_controller`) non potrebbe girarci sopra. Il Pade non e' una scorciatoia:
   e' il prezzo per poter fare controllo classico. Usare il ritardo esatto in
   Simulink e il Pade nello script farebbe divergere i due track **per
   costruzione**, e la validazione perderebbe di significato.
2. **I guadagni PD del track Simulink non sono quelli del progetto di Task 2/3.**
   `init_simulink_hm3` (riga 34) tara il PD sul plant *rigido* con attuatore
   *ideale* (Kp = 1.7845, Kd = 0.4433 = progetto Task 1); `main_task2.m` (riga
   153) e `main_task3.m` (riga 23) lo ri-tarano sul loop completo
   (Kp = 1.7318, Kd = 0.6867). Sul loop completo la differenza e' seria:
   **Rigid PM 14.6 deg contro 30.0 deg**, delay margin 98 ms contro 165 ms
   (entrambi stabili). La sovrapposizione script-vs-Simulink *funziona lo stesso*,
   perche' la riga 43 di questo file prende `K` da `S` -- quindi anche la baseline
   analitica e' costruita con i guadagni "sbagliati". La validazione e'
   **internamente coerente ma valida il loop pre-retune**, non quello consegnato.
3. **Orizzonte diverso.** Qui si simula fino a `Tstop` = 12 s (default di
   `load_wind_profile`), mentre `main_task2/3` usano `Tend = 80 s` per catturare
   il modo di drift lento (tau ~ 20 s). Confrontando le figure Simulink con quelle
   del report, gli assi temporali non coincidono. (Per i picchi non cambia nulla:
   verificato che sul corner V3 i massimi cadono comunque entro i primi 12 s.)
4. **Nessuna asserzione, nessun test.** La funzione produce una figura, non un
   verdetto. Non c'e' un `max(abs(rs.theta - interp1(sl.t,sl.theta,rs.t)))`
   confrontato con una soglia, e in `HM3/tests/` non c'e' nessun test che tocchi
   il track Simulink. La tolleranza "accettata" e' quella riportata a mano nel
   report.
5. **Numeri stantii nella guida.** `models/SIMULINK_GUIDE.md` (righe 11-13)
   dichiara accordo "~1e-7 rad per Task 1 e ~1e-8 rad per Task 2", mentre la
   tabella del report misura 3.2e-6 (Task 1) e 1.7e-7 (Task 2). E le righe
   222-224 della guida danno per il corner V3 picchi di 0.35 deg / 4.2 m
   "matching main_task3.m": `main_task3.m` da' 0.878 deg / 4.92 m, e il modello
   cosi' inizializzato da' 1.402 deg / 4.88 m. La guida non coincide con nessuno
   dei due.

---

## Possibili domande d'esame

**D: A cosa serve davvero il track Simulink, se gli script producono gia' tutti i
risultati?**
R: A due cose. **Validazione incrociata**: due implementazioni indipendenti dello
stesso anello (trasferimento + `lsim` da un lato; blocchi + `ode45` dall'altro,
con realizzazioni di stato e ordine dei blocchi diversi) che devono dare la
stessa risposta temporale al medesimo vento. Un segno invertito nel PD, le
colonne di B scambiate, il demux ricablato male romperebbero la sovrapposizione.
E' un test di *trascrizione*. Seconda ragione: il controllo d'assetto **si
consegna come modello** (Simulink -> autocodifica -> HIL), non come script di
analisi; il diagramma e' il formato con cui il progetto entra in un flusso
industriale.

**D: Quale solver usa il modello e perche' non uno stiff, visto che c'e' un modo
a zeta = 0.005?**
R: `ode45` a **passo variabile**, `RelTol = 1e-6` (letto dalla configurazione
del modello). Un modo poco smorzato **non e' un modo stiff**: la stiffness nasce
da scale temporali *decadenti* molto separate, e qui il rapporto fra il polo
piu' veloce (Pade, |lambda| = 256 rad/s) e il piu' lento (drift, 0.233 rad/s) e'
solo ~1.1e3 -- modesto, e su 12 s di simulazione un solver esplicito lo digerisce
senza problemi. Il bending pone invece un problema di **accuratezza**: con
Q = 100 e tempo di decadimento 10.6 s, il modo compie ~36 periodi senza
smorzarsi, e un solver implicito con smorzamento numerico (BDF di `ode15s`,
`ode23tb`) lo attenuerebbe artificialmente -- falsificando proprio la grandezza
che il notch deve controllare. `ode45` con tolleranza stretta paga passi piccoli
ma preserva l'ampiezza della risonanza.

**D: Nel modello c'e' un algebraic loop?**
R: No, e si dimostra guardando i feedthrough. Il `Gain` del PD ha feedthrough; il
`Transfer Fcn` del notch e' biproprio (grado 2 su grado 2) e quindi anche lui ha
feedthrough; ma il `Transfer Fcn` del TVC e' strettamente proprio (numeratore di
grado 3, denominatore di grado 5) e il blocco State-Space ha `D = zeros(7,2)`.
La catena di feedthrough e' interrotta due volte, quindi il ciclo di retroazione
non e' algebrico. E' fragile pero': se il plant avesse un `D` non nullo, il notch
biproprio chiuderebbe l'anello algebrico e Simulink dovrebbe risolvere un punto
fisso a ogni passo. Che `D = 0` non e' un caso: la deflessione dell'ugello agisce
sull'accelerazione, e arriva alle uscite solo dopo due integrazioni.

**D: Che tolleranza di accordo si accetta fra i due track, e cosa vorrebbe dire
uno scostamento?**
R: Il codice non ne impone nessuna -- produce una figura, non un numero (ed e' un
limite). Le misure riportate nel report stanno a max|d theta| = 1.7e-7 rad
(Task 2) e 3.2e-6 rad (Task 1), con drift sub-millimetrico: rapportate al picco
del run Simulink (~5e-3 rad) sono errori relativi ~3e-5 per il Task 2, cioe'
**errore numerico puro** (RelTol = 1e-6 +
diversa ricostruzione dell'ingresso: `lsim` mantiene l'ingresso fra i campioni,
`From Workspace` lo interpola linearmente sugli istanti di `ode45`). Uno
scostamento di ordine 1e-3 relativo, ma con la stessa forma d'onda, indicherebbe
un problema di tolleranza/passo; forme d'onda diverse, segno opposto o ampiezza
raddoppiata indicherebbero un errore **strutturale** di cablaggio -- ed e' quella
la classe di errori che questa validazione esiste per catturare.

**D: Come fai a garantire che i due track vedano lo stesso vento?**
R: Non lo rigenero. La riga 46 ricostruisce il vento **dal `timeseries` che il
modello stesso legge** (`S.wind_ts.Time`, `S.wind_ts.Data`), quindi lo stesso
oggetto alimenta il blocco `From Workspace` di Simulink e il `lsim` della
baseline. E' essenziale con `profile = 'strongwind'`, dove il vento viene dal
generatore stocastico del professore: rigenerarlo darebbe una realizzazione
diversa (anche a seme fisso, con orizzonti diversi) e la sovrapposizione perderebbe
significato.

**D: Se domani cambi il `VariableName` di un blocco `To Workspace`, cosa succede?**
R: **Niente di visibile**, ed e' un problema. (Rinominare il *blocco* -- oggi
`log_theta`, `log_z`, `log_zdot`, `log_delta` -- non fa invece nulla: il
contratto e' sul `VariableName`, non sul nome del blocco.) `get_logged_signal`
(righe 86-102) cerca il nome prima in `logsout`, poi fra le proprieta' del
`SimulationOutput`; se non lo trova ritorna `t = []`, `y = []` senza errore ne'
warning (il `catch` scatta solo se l'accesso *lancia*). La figura viene disegnata
con la sola curva blu dello script e salvata come se tutto fosse a posto. Il
contratto sui nomi (`theta_sl`, `z_sl`, `zdot_sl`, `delta_sl`) e' documentato ma
non verificato: basterebbe un `assert(~isempty(sl.theta))` dopo la riga 57.

**D: Il modello Simulink e' il progetto che hai consegnato in Task 2?**
R: Onestamente no. La baseline e il modello girano entrambi con i guadagni
tarati da `init_simulink_hm3` sul plant rigido con attuatore ideale
(Kp = 1.784, Kd = 0.443, cioe' il progetto di Task 1), non con quelli ri-tarati
sul loop completo dai `main_task2/3` (Kp = 1.732, Kd = 0.687). Il loop resta
stabile ma con margine di fase rigido 14.6 deg invece di 30 deg -- cioe' la
configurazione che `main_task2.m` stampa come "BEFORE re-tuning". La validazione
resta valida come test di equivalenza (i due track condividono la stessa `K`),
ma valida il loop pre-retune. E' un bug da una riga in `init_simulink_hm3.m`
(riga 34), non un problema del modello.
