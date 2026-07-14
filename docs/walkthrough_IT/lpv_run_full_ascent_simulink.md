# HM3/LTV_FULL_ASCENT/run_full_ascent_simulink.m

## Ruolo del file nel progetto

Questo e' il file di **validazione incrociata** del track LPV rigido: simula il
modello `hm3_full_ascent.slx` (nelle due varianti *frozen* e *scheduled*) e
sovrappone il risultato alla **baseline MATLAB** integrata con `ode45` sul RHS
`ode_lpv_ascent`. Produce la figura `figures/fullascent_simulink_vs_script.png`
e stampa i residui massimi su `theta`, `z` e `delta`.

La domanda a cui risponde e' semplice ma non banale: *il modello Simulink,
costruito da script con lookup table e blocchi elementari, e' davvero lo stesso
sistema che integro in MATLAB?* La risposta dichiarata dal README della cartella
e' "si', a ~1e-7 rad su theta" (1.1e-7 frozen, 2.2e-7 scheduled). Questa pagina
spiega **come** si ottiene un numero cosi' stretto e **cosa** lo romperebbe.

Il rapporto gerarchico fra i due track e' esplicito nel codice e nel README:
**`ode45` e' la fonte di verita'**, Simulink e' la replica. Non e' una scelta
scontata (in molti progetti industriali e' il contrario); qui ha senso perche'
il design -- plant, PD, margini -- e' tutto scritto in MATLAB, e il `.slx` serve
a dimostrare che la stessa legge di controllo gira in un ambiente di simulazione
"di volo", con il generatore di vento del professore **dentro l'anello**.

Chiamato a mano dall'utente (vedi il *How to run* del README). Dipende da:
`init_simulink_lpv` (dati + variabili base), `build_hm3_full_ascent` (autore del
`.slx`), `ode_lpv_ascent` (RHS della baseline).

---

## Firma e opzioni (righe 1-26)

```matlab
function out = run_full_ascent_simulink(o)
arguments
    o.rebuild (1,1) logical = false
end
...
if o.rebuild || ~isfile(mdlfile)
    build_hm3_full_ascent();
end
```

- Riga 1: ritorna `out`, una struct con un campo per variante (`frozen`,
  `scheduled`), ciascuno con `t` e `err`.
- Righe 24-26: **il `.slx` e' un artefatto derivato**. Se manca, viene
  ri-generato da `build_hm3_full_ascent` (che *e'* la definizione del modello:
  ~150 righe di `add_block` / `add_line`). Con `'rebuild', true` lo si rigenera
  sempre. Nota pero' che il `.slx` **e' comunque committato in repo** -- quindi
  in pratica coesistono la sorgente (lo script) e il binario (il modello); il
  meccanismo garantisce che il binario sia riproducibile, non che sia sempre
  aggiornato. Se si modifica `build_hm3_full_ascent.m` e non si ri-esegue con
  `'rebuild', true`, si simula un modello vecchio senza accorgersene.

> **Possibile domanda d'esame** -- Perche' costruire il modello Simulink da
> script invece che disegnarlo a mano?
> *Risposta:* Tre motivi concreti. (1) **Riproducibilita'**: il `.slx` e' un
> binario che il version control non sa diffare; lo script `.m` invece si legge,
> si diffa e si rivede in code review -- lo script *e'* il modello. (2) **La
> sorgente resta unica**: i parametri dei blocchi non sono numeri incollati nel
> canvas ma nomi (`lpv_c1`, `Kp_th0`, `tsched`) risolti dal workspace, quindi il
> modello eredita automaticamente qualunque cambio di design fatto in MATLAB.
> (3) **Rigenerabilita'**: se cambia il design (nuovo passo dello schedule,
> nuovo orizzonte, nuovo notch) basta ri-eseguire il build. Il contrasto e' con
> il modello a mano di HM3 (`hm3_closed_loop.slx`, costruito seguendo
> `models/SIMULINK_GUIDE.md`), che va rifatto a mano a ogni cambio.

---

## Setup: dati, tolleranze, modello (righe 28-37)

```matlab
S    = init_simulink_lpv();
odeo = odeset('RelTol', 1e-9, 'AbsTol', 1e-11);
load_system(mdlfile);
variants = {'frozen', 0; 'scheduled', 1};
```

- Riga 28: **`init_simulink_lpv()` con i default** -> `push = true`: spinge nel
  base workspace tutte le variabili che il modello risolve per nome, e
  restituisce `S` (gli interpolanti) per la baseline. Le due strade partono
  quindi *dagli stessi numeri*: e' la precondizione numero uno dell'accordo a
  1e-7.
  **Costo nascosto:** questa chiamata ri-simula il generatore di vento
  (`strong_wind.slx`) e **ri-progetta i 28 punti dello schedule** con
  `fminsearch` -- decine di secondi -- e il vento cosi' generato (`S.wind`,
  `S.windfun`) **non viene poi usato** (vedi `pack_model`). Lavoro sprecato: non
  esiste un'opzione per saltare il generatore.
- Riga 29: tolleranze della baseline: `RelTol = 1e-9`, `AbsTol = 1e-11`. Molto
  strette, coerenti con la convenzione della repo per il lavoro di precisione.
  **Nota importante:** il *modello Simulink* NON usa queste tolleranze -- la
  `set_param` di `build_hm3_full_ascent` (righe 35-37) fissa `Solver = ode45`,
  `SolverType = Variable-step`, `MaxStep = 0.02`, piu' `StartTime`/`StopTime` e il
  logging, ma **non** tocca ne' `RelTol` ne' `AbsTol`: il modello gira quindi con le
  tolleranze **di default** di Simulink (`RelTol = 1e-3`). I due integratori
  hanno tolleranze diverse di sei ordini di grandezza: l'accordo a 1e-7 **non**
  viene da tolleranze appaiate (vedi sotto).
- Riga 32: le due varianti sono solo due valori dello scalare `sched`: lo stesso
  canvas, un blocco Switch. Non ci sono due modelli.

---

## Ciclo sulle due varianti (righe 39-70)

### Simulazione e lettura dei log (righe 40-45)

```matlab
assignin('base', 'sched', sc);
so = sim(mdl, 'StopTime', num2str(S.Tstop));
th = so.theta_sl;  z = so.z_sl;  de = so.delta_sl;  aw = so.alpha_w_sl;
tt = th.Time;
```

- Riga 41: la variante si seleziona **scrivendo nel base workspace**. Il blocco
  Constant `sched` legge quella variabile, il blocco Switch confronta con la
  soglia 0.5 e passa il ramo schedulato (`u2 >= 0.5`) o quello congelato.
  Accoppiamento fragile (nessun controllo statico), ma e' il pattern Simulink
  standard senza Data Dictionary.
- Riga 42: `so` e' un `Simulink.SimulationOutput`; i cinque segnali arrivano dai
  blocchi To Workspace del modello, salvati in formato `Timeseries` con
  `SampleTime = -1` e `MaxDataPoints = inf` (`build_hm3_full_ascent`, riga 136):
  **ogni major step del solver, nessun tetto a 1000 punti**. Il tetto di default
  avrebbe decimato il log e reso impossibile un confronto a 1e-7.
- Riga 45: `tt` sono i tempi dei passi Simulink. Con `MaxStep = 0.02` su 140 s
  sono circa 7000 punti, di fatto uniformi.
  **Fragilita':** `tt` (cioe' `th.Time`) viene passato tale e quale come `tspan` a
  `ode45` (riga 49), mentre la griglia del `griddedInterpolant` del vento e'
  `aw.Time` (`pack_model`, riga 92) -- **non** `tt`: sono due log distinti prodotti
  dagli stessi passi del solver, quindi in pratica gli stessi istanti, ma sono due
  vettori diversi. Entrambi gli usi pretendono tempi **strettamente crescenti**, e
  i log a passo variabile possono ripetere un istante (e infatti
  `init_simulink_lpv`, righe 171-172, si difende con `unique`) -- qui quella
  difesa **non c'e'**, ne' su `tt` ne' su `aw.Time`. In pratica non capita, ma e'
  una fragilita' reale.

### Replay ode45 sullo stesso vento (righe 47-52)

```matlab
M = pack_model(S, sc, aw);
[~, x] = ode45(@(t, x) ode_lpv_ascent(t, x, M), tt, zeros(4,1), odeo);
if sc, Kp = S.fKp(tt); Kd = S.fKd(tt);
else,  Kp = S.K0.Kp_th*ones(size(tt)); Kd = S.K0.Kd_th*ones(size(tt)); end
del_ode = -(Kp.*x(:,3) + Kd.*x(:,4) + S.K0.Kp_z*x(:,1) + S.K0.Kd_z*x(:,2));
```

- Riga 48: `pack_model` (righe 80-93) e' il cuore metodologico del file:
  costruisce la struct per `ode_lpv_ascent` usando **il vento che il modello ha
  appena prodotto**, non `S.windfun`:

      'windfun', griddedInterpolant(aw.Time, squeeze(aw.Data), ...)

  Questo elimina d'un colpo la piu' grande sorgente di scostamento. Se si usasse
  `S.windfun`, si confronterebbero **due realizzazioni diverse dello stesso
  processo stocastico**: il generatore standalone (simulato dentro
  `init_simulink_lpv`) e la sua copia dentro `hm3_full_ascent.slx` girano con
  solver e passi diversi, e i filtri di Dryden sono a stati continui -- stesso
  seed, ma uscite leggermente diverse. In piu' il modello calcola
  `alpha_w = v_w * interp(1/V)` mentre `init_simulink_lpv` (riga 90) calcola
  `alpha_w = v_w / interp(V)`: **due funzioni diverse fra i breakpoint** (scarto
  relativo misurato ~1.1e-4). Nessuno dei due errori sarebbe compatibile con un
  residuo di 1e-7 rad.
- Riga 49: `ode45` con `tspan = tt` (vettore lungo): **non** forza i passi
  interni a coincidere con quelli di Simulink -- restituisce solo l'uscita
  *densa* valutata su `tt`. I due integratori seguono percorsi diversi; il fatto
  che arrivino allo stesso posto e' proprio la verifica.
- Righe 50-52: `delta` **non e' uno stato** (attuatore ideale nel modello
  rigido): va ricostruita algebricamente dagli stati con gli stessi guadagni.
  Nel caso scheduled si usa `S.fKp(tt)`, lo stesso interpolante che alimenta il
  lookup del modello; sotto `tsched(1) = 5 s` l'interpolante estrapola
  `'nearest'` (tiene l'estremo) e il lookup del modello ha `ExtrapMethod =
  'Clip'` (tiene l'estremo): **stessa funzione**. E' esattamente cio' che dichiara
  il commento in `build_hm3_full_ascent`, righe 106-107. Se il lookup avesse
  l'estrapolazione lineare di default, l'overlay scheduled avrebbe un errore
  visibile nei primi 5 s.
- Nota di coerenza: nel ramo `else` i guadagni di drift restano `S.K0.Kp_z` e
  `S.K0.Kd_z` **in entrambe le varianti** -- lo schedule tocca solo il PD di
  pitch, e cosi' fa il modello (i blocchi `Gain_Kpz` / `Gain_Kdz` sono costanti).

### Residui e figura (righe 54-69)

```matlab
err = struct('theta', max(abs(x(:,3) - squeeze(th.Data))), ...
             'z',     max(abs(x(:,1) - squeeze(z.Data))), ...
             'delta', max(abs(del_ode - squeeze(de.Data))));
```

- Righe 54-56: **norma infinito** dell'errore, non RMS: la metrica piu' severa,
  giusta per una validazione.
- Righe 62-69: overlay con baseline `ode45` in linea continua e Simulink
  tratteggiato -- la convenzione classica: se il tratteggiato copre il continuo,
  i due coincidono. Solo il tile di `theta` riporta il residuo nel titolo
  (riga 64, `\Delta=%.1e rad`); i tile `z` e `delta` (righe 67 e 69) hanno il
  solo nome della variante -- `err.z` e `err.delta` finiscono unicamente nella
  `fprintf` (righe 58-59).

---

## Perche' l'accordo e' 1e-7 rad -- cosa deve essere identico

Un accordo di **1e-7 rad su theta**, con picco di `theta ~ 0.95 gradi = 1.7e-2
rad`, e' un errore **relativo di ~1e-5**. Perche' scenda cosi' in basso devono
valere *tutte* queste condizioni, e ognuna e' verificabile sul codice:

1. **Stessi coefficienti fra i breakpoint.** Il modello interpola linearmente
   le tabelle `lpv_c1..lpv_c7`; `ode_lpv_ascent` (righe 26-28) chiama
   `M.fc1..M.fc7`, cioe' `griddedInterpolant('linear')` costruiti **sulle stesse
   tabelle e sugli stessi breakpoint** (`init_simulink_lpv`, righe 101-103).
   Interpolazione lineare in entrambi i casi -> **stesso sistema LTV**, non solo
   sui nodi ma ovunque. E' possibile solo perche' `init_simulink_lpv` tabula i
   *coefficienti efficaci* (`c2 = a1*V + a4` in un'unica tabella): se il modello
   tabulasse `a1` e `V` separatamente e li moltiplicasse, sarebbe un sistema
   diverso (scarto relativo ~5e-5, misurato).
2. **Stesso vento.** Garantito per costruzione: la baseline riceve l'`alpha_w`
   *loggato dal modello*. Resta un residuo, perche' il replay interpola
   linearmente fra i campioni loggati mentre il generatore, dentro Simulink, e'
   un filtro continuo eccitato da rumore a passo 0.1 s: fra due campioni le due
   funzioni non coincidono. Con `MaxStep = 0.02` i campioni sono fitti e questo
   errore e' piccolo -- l'header del file (righe 3-5) lo dichiara come **la**
   sorgente residua.
3. **Stessa famiglia di solver:** ode45 (Dormand-Prince) da entrambe le parti.
4. **Stesse condizioni iniziali:** tutti gli integratori del modello hanno
   `InitialCondition = 0`, la baseline parte da `zeros(4,1)`.
5. **Stessa interpolazione dei guadagni**, estrapolazione inclusa (Clip contro
   'nearest'): vedi sopra.

**Cosa NON e' identico, e va detto:**

- **Le tolleranze del solver.** Il modello gira con i default di Simulink
  (`RelTol = 1e-3`), la baseline con `1e-9`. L'accordo non viene dalle
  tolleranze appaiate: viene da **`MaxStep = 0.02`**, imposto nel build (riga
  36). Con passo forzato a 0.02 s su una dinamica di scala ~1 s, l'errore locale
  di ode45 (ordine 5) e' dell'ordine di `h^5 ~ 3e-9` per passo: il controllo di
  tolleranza a 1e-3 non morde mai, e l'accuratezza effettiva e' fissata dal
  passo. Il commento nel build (righe 33-34) motiva `MaxStep` con la necessita'
  di risolvere le innovazioni di rumore a 0.1 s del generatore -- l'accuratezza
  e' un effetto collaterale (fortunato). Se si alzasse `MaxStep`, il residuo
  peggiorerebbe, e non per colpa del modello.
- **I passi.** Simulink cammina a 0.02 s; `ode45` in MATLAB sceglie i suoi passi
  (grandi, vista la tolleranza) e produce l'uscita su `tt` per **dense output**
  (interpolante di ordine 4 interno a `ode45`). Le due traiettorie di
  integrazione sono diverse: coincidono perche' entrambe convergono.

**Cosa lo romperebbe (in ordine di gravita'):**

- usare `S.windfun` invece dell'`alpha_w` loggato (due realizzazioni diverse del
  vento + asimmetria `1/interp(V)` contro `interp(1/V)`) -> errore di ordini di
  grandezza superiore;
- tabulare i coefficienti in modo diverso fra i due track (e' *esattamente* cio'
  che succede nel modello flessibile, vedi `lpv_run_flex_simulink.md`);
- cambiare il metodo di interpolazione di un lookup (es. `Flat` / `Nearest`
  invece di `Linear`), o l'`ExtrapMethod` dei guadagni;
- alzare `MaxStep` o togliere `MaxDataPoints = inf` dai To Workspace (il log
  verrebbe decimato e il confronto perderebbe risoluzione);
- introdurre nel modello uno stato che la baseline non ha (per esempio un
  attuatore reale al posto di quello ideale) senza aggiornare `ode_lpv_ascent`.

> **Possibile domanda d'esame** -- Il residuo e' 1e-7 rad: e' "tolleranza del
> solver" come dice il README?
> *Risposta:* In parte. Le due integrazioni **non** hanno la stessa tolleranza
> (Simulink gira ai default 1e-3, la baseline a 1e-9): a rendere accurato il
> lato Simulink e' il vincolo `MaxStep = 0.02`, non la tolleranza. Il residuo e'
> quindi la somma di due contributi: l'errore di troncamento accumulato del lato
> Simulink a passo 0.02 s, e l'errore di interpolazione lineare del vento fra i
> campioni loggati. Dire "tolleranza del solver" e' un'abbreviazione accettabile
> ma non e' letteralmente cio' che il codice fa.

---

## `pack_model` (righe 80-93)

```matlab
M = struct('fc1', S.fc1, ..., 'fc7', S.fc7, 'fKp', S.fKp, 'fKd', S.fKd, ...
           'Kp_th0', S.K0.Kp_th, 'Kd_th0', S.K0.Kd_th, ...
           'Kp_z', S.K0.Kp_z, 'Kd_z', S.K0.Kd_z, ...
           'sched', logical(sched), ...
           'windfun', griddedInterpolant(aw.Time, squeeze(aw.Data), ...));
```

- Righe 88-92: e' la stessa struct che costruisce `make_model` dentro
  `main_full_ascent.m` (righe 130-133), con **una sola differenza**: `windfun`.
  La' e' `S.windfun` (il generatore standalone); qui e' il vento loggato dal
  modello. La duplicazione della struct in due file e' una piccola violazione
  DRY, ma e' voluta: separa il "run scientifico" (main) dal "run di validazione"
  (questo).
- Riga 91: `logical(sched)` -- `ode_lpv_ascent` fa `if M.sched` (riga 17), quindi
  gli basterebbe lo 0/1 numerico; il cast e' difensivo.

---

## Possibili domande d'esame

**D: Perche' la baseline `ode45` e' la "fonte di verita'" e non il modello
Simulink?**
R: Perche' tutto il design (plant `build_plant_rigid`, tuner `design_controller`,
margini, notch) e' scritto in MATLAB e il RHS `ode_lpv_ascent` e' la
trascrizione diretta delle equazioni; il `.slx` e' invece una *ricostruzione*
delle stesse equazioni con blocchi. La direzione della verifica va dalla forma
piu' vicina alla matematica a quella piu' vicina all'implementazione. Il valore
del modello Simulink non e' l'accuratezza (che e' identica) ma il fatto di poter
mettere il **generatore di vento del professore dentro l'anello**, senza
ricampionarlo o ritagliarlo.

**D: Il vento della baseline e' quello generato dal modello. Non e' un modo di
"barare" sul confronto?**
R: No: e' l'unico modo di isolare cio' che si vuole misurare. L'ingresso
esogeno (il vento) e' un processo stocastico generato da un modello Simulink a
stati continui; due simulazioni con solver diversi danno realizzazioni
leggermente diverse anche con lo stesso seed. Se si lasciassero divergere gli
ingressi, il residuo misurerebbe la differenza fra i venti, non fra i due
modelli di veicolo. Fissando l'ingresso, il residuo misura esattamente cio' che
interessa: la fedelta' della trascrizione della dinamica e del controllore. Il
prezzo e' che la validazione **non** copre il generatore stesso (quello resta
non verificato, ma e' materiale del professore preso as-is).

**D: Il modello ha un solo canvas per due controllori. Come fa?**
R: Un blocco Constant legge lo scalare `sched` dal base workspace; due blocchi
Switch (`Kp_sw`, `Kd_sw`, criterio `u2 >= 0.5`) selezionano fra i lookup dello
schedule (`Kp_s`, `Kd_s`, tabelle `Kp_sched`/`Kd_sched` su breakpoint `tsched`)
e le costanti congelate (`Kp_th0`, `Kd_th0`). Il resto del canvas -- plant,
generatore, drift feedback -- e' condiviso. Cosi' le due varianti vedono
**esattamente lo stesso vento e lo stesso plant**, e la differenza fra le due
risposte e' attribuibile solo al controllore.

**D: Quali sono i limiti onesti di questa validazione?**
R: (1) Valida solo il **plant rigido con attuatore ideale**: nessun TVC, nessun
bending. (2) Non valida il generatore di vento (usato come ingresso comune).
(3) Copre un solo scenario di vento -- una sola realizzazione, nessun Monte
Carlo. (4) `delta` non e' confrontata come segnale del modello contro segnale del
modello, ma contro una ricostruzione algebrica dagli stati della baseline: se
il modello sbagliasse *la stessa* formula in modo consistente, il confronto non
lo vedrebbe (in realta' non e' cosi', perche' nel modello `delta` viene dalla
somma dei quattro rami di feedback e non da una formula riscritta). (5) Il
residuo dipende da `MaxStep = 0.02`, non da tolleranze dichiarate: non c'e'
alcuna garanzia formale, e' un accordo empirico.

**D: `run_full_ascent_simulink` chiama `init_simulink_lpv()` a ogni esecuzione.
Costo?**
R: Alto e in parte inutile: ogni chiamata ri-simula `strong_wind.slx` sull'intero
orizzonte e ri-progetta i 28 punti dello schedule con `fminsearch` (ognuno con
decine di valutazioni di `margin`). Il vento cosi' generato viene poi
**scartato** (il replay usa quello loggato dal modello). Sarebbe naturale un
`o.skip_wind` o una cache su disco; il codice non ce l'ha.
