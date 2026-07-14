# HM3/tests/hm3LoopPerformanceTest.m

## Ruolo del file nel progetto

E' l'unica suite di HM3 che **non verifica la correttezza ma il costo**. Estende
`matlab.perftest.TestCase` (non `matlab.unittest.TestCase`) e misura il tempo di
esecuzione dei blocchi che stanno **dentro i cicli caldi** dell'homework. Si lancia
con `runperf('hm3LoopPerformanceTest')`: e' l'unico modo per ottenere i **tempi**.
Attenzione pero' -- `matlab.perftest.TestCase` e' una sottoclasse di
`matlab.unittest.TestCase`, quindi anche `runtests('HM3/tests')` **scopre ed esegue**
questa classe (eseguendo una volta i corpi dei cicli e valutando le `verify*`);
semplicemente non produce misure.

**Perche' un perf-test in un homework di controlli.** Perche' HM3, a differenza di
HM1 e HM2, non risolve *un* problema: ne risolve **migliaia**. Il progetto non e'
un singolo anello, e' una **ricerca** su uno spazio di parametri, e la ricerca e'
costruita interamente attorno a una funzione, `assemble_loop.m`, che viene
richiamata da tre punti diversi:

1. **L'auto-tuner PD** (`design_controller.m`, righe 56-57): `fminsearch` valuta la
   funzione di costo fino a `MaxFunEvals = 400` volte, e **ogni valutazione**
   chiama `assemble_loop` + `minreal` + `classify_margins` + `isstable` (righe
   82-89).
2. **Il trade dei filtri del Task 2**: lo sweep sulle parametrizzazioni del
   lead-lag (il README parla di "best of 75") e sulle varianti di notch richiede
   un riassemblaggio dell'anello per ogni candidato.
3. **Il Monte Carlo** (`main_montecarlo.m`, riga 27: `N = 1500`; riga 87:
   `[L, T] = assemble_loop(Gf, K, Wf)` **dentro il `parfor`**): 1500
   riassemblaggi dell'anello completo, uno per estrazione.

`assemble_loop` **non e' una funzione economica**. Fa `connect()` (interconnessione
per nome, che costruisce una realizzazione in spazio di stato allargata e poi
elimina i segnali interni), `getLoopTransfer`, una conversione `ss -> tf` e un
`minreal`. Misurato su questa macchina:

| Operazione | Costo misurato |
|---|---|
| `assemble_loop` (rigido, 4 stati) | ~10.5 ms/chiamata |
| `assemble_loop` (completo + TVC + notch, 13 stati) | ~9.1 ms/chiamata |
| `classify_margins` (via `allmargin`) | ~3.7 ms/chiamata |
| costo del tuner (`assemble_loop` + `margin` + `isstable`) | ~15.0 ms/chiamata |
| `simulate_gust_response` (`lsim`, 2401 campioni) | ~5.8 ms/chiamata |
| `design_controller` (una ricerca completa) | ~830 ms |

Dieci millisecondi sembrano nulla. Moltiplicati per il Monte Carlo diventano
**~14 secondi solo di assemblaggio** (1500 x 9.1 ms), a cui vanno aggiunti i
margini: la stima realistica e' 20-25 s in seriale. E' esattamente la ragione per
cui `main_montecarlo.m` usa `parfor` (riga 75). Se `assemble_loop` rallentasse di
un fattore 3 -- per esempio perche' qualcuno aggiunge un `minreal` in piu' o
sostituisce `connect` con una costruzione meno efficiente -- il Monte Carlo
passerebbe da "si fa mentre prendi un caffe'" a "si fa stanotte", e la cosa
verrebbe scoperta **al momento sbagliato**. Questa suite serve a scoprirlo prima.

**Questa suite ha ospitato, fino a poco fa, un golden value obsoleto che la
rompeva.** Ora e' allineato e la suite HM3 passa **35/35**, ma la storia e' troppo
istruttiva per cancellarla: vedi la sezione su `Kref` qui sotto.

---

## `Kref` -- un golden value che era diventato STALE (righe 8-13)

```matlab
properties (Constant)
    % Task-1 PD design (pinned), same as hm3LoopTest. Keep the two in sync:
    % testDesignControllerSearch asserts the tuner converges back to Kp_th.
    Kref = struct('Kp_th', 1.7845, 'Kd_th', 0.4433, ...
                  'Kp_z', -1e-3,   'Kd_z', -1e-3)
end
```

Oggi (righe 11-12) `Kref` vale `Kp_th = 1.7845, Kd_th = 0.4433`, **esattamente**
come `hm3LoopTest.m` (riga 11), ed e' il progetto che `design_controller` produce
davvero: `Kp_th = 1.784455`. Il commento delle righe 9-10 e' stato riscritto perche'
dica una cosa in piu' del semplice "same as hm3LoopTest": **avverte di tenere i due
file in sync** e spiega *perche'* (`testDesignControllerSearch` asserisce che il
tuner riconverga su quel numero).

Quell'avvertimento e' li' perche' e' stato pagato.

### Com'era, e perche' era rotto

Fino a poco fa questa property pinnava

    Kref = struct('Kp_th', 1.9800, 'Kd_th', 1.3997, ...)

cioe' i guadagni **PRE-ritaratura**, e il commento sosteneva -- falsamente -- di
essere "same as hm3LoopTest". `hm3LoopTest.m` pinnava gia' `1.7845 / 0.4433`. **I
due `Kref` non coincidevano**, e il secondo non corrispondeva a nessun progetto
prodotto dal codice.

Conseguenza, **osservata eseguendo `runperf`**: il test `testDesignControllerSearch`
(righe 91-98) **falliva**.

    verifyEqual failed.
        Actual Value:    1.784454942663017
        Expected Value:  1.980000000000000
        AbsoluteTolerance: 0.005
    Totals: 0 Valid, 1 Invalid.

Il `1.9800 / 1.3997` era il residuo di **una versione precedente del progetto**,
antecedente al re-tuning sull'anello completo -- quello che ha portato i guadagni
canonici di D'Antuono a essere raffinati contro i **margini classificati**. Quando
la metodologia dei margini e' stata riscritta, il numero pinnato in `hm3LoopTest` e'
stato aggiornato e **questo e' rimasto indietro**.

**Perche' era peggio di un semplice test rosso.** In `matlab.perftest`, una verifica
fallita **invalida la misura**: `runperf` riportava *0 Valid, 1 Invalid*. Il test
non era solo "rosso", **non produceva alcun dato di tempo** per il livello workflow.
Cioe' il benchmark piu' importante della suite -- quanto costa una ricerca completa
del tuner -- non misurava niente.

Gli altri quattro test **passavano lo stesso**, ed e' proprio questo che ha reso il
bug longevo: usano `Kref` solo come *insieme di guadagni con cui chiudere l'anello
per misurare i tempi*, e il tempo di `assemble_loop` **non dipende dal valore
numerico dei guadagni** (la struttura e l'ordine dell'anello sono gli stessi). Con
`1.98 / 1.3997` l'anello rigido era comunque stabile e `margin()` restituiva valori
finiti, quindi anche le asserzioni di sanita' delle righe 78-79 reggevano. Il valore
stale era **invisibile** a tutto tranne che all'unico test che lo confrontava con
l'output del tuner.

### La correzione, e la lezione

La correzione e' stata **una riga** -- allineare `Kref` a `1.7845 / 0.4433` -- piu'
il commento di avvertimento. La correzione *strutturale*, quella che il codice non
ha ancora fatto, sarebbe far derivare entrambe le suite da un'**unica sorgente di
verita'** (una funzione, un `.mat`, una property condivisa) invece di duplicare i
numeri in due file. **La duplicazione e' la causa prima del bug**, non il numero
sbagliato: finche' il valore vive in due posti, prima o poi i due posti divergono.

> **Possibile domanda d'esame** -- Come e' potuto succedere che un test pinni dei
> guadagni che il codice non produce piu'?
> *Risposta:* E' un classico *stale golden value*, e la meccanica e' istruttiva. Il
> valore era **duplicato in due file** (`hm3LoopTest.m` e
> `hm3LoopPerformanceTest.m`); quando ho ri-tarato il progetto sull'anello completo
> -- riscrivendo la metodologia dei margini da `margin()` a `classify_margins` --
> ho aggiornato il test unitario e **non** il perf-test, che e' rimasto sui guadagni
> pre-ritaratura (1.98 / 1.40 invece di 1.78 / 0.44). Il fallimento e' rimasto
> nascosto perche' 4 test su 5 usano `Kref` solo come guadagni con cui chiudere
> l'anello per **cronometrarlo**, e il tempo non dipende dal loro valore: solo
> `testDesignControllerSearch`, che confronta l'output del tuner con il numero
> pinnato, se ne accorgeva -- e in `matlab.perftest` una verify fallita **invalida
> la misura**, quindi il benchmark di livello workflow non produceva piu' alcun
> dato (*0 Valid, 1 Invalid*). **Attenzione a non giustificarlo dicendo che i
> perf-test "sfuggono a `runtests`"**: `matlab.perftest.TestCase` e' una
> sottoclasse di `matlab.unittest.TestCase`, quindi `runtests('HM3/tests')` -- il
> comando prescritto dalle convenzioni del repo -- scopre ed esegue anche questa
> classe, e il `verifyEqual` falliva li' esattamente come sotto `runperf`. La causa
> vera e' che non ho rilanciato la suite per intero dopo il re-tuning. Ora e'
> allineato (35/35), ma **la lezione e' che un golden value non va duplicato**: va
> derivato da un'unica sorgente. Un golden value non e' una verita' eterna, e' una
> **fotografia di un progetto** -- e invecchia esattamente quando il progetto
> cambia, cioe' nel momento in cui sei meno propenso a rileggere i test.

---

## `TestClassSetup` e `TestMethodSetup` (righe 24-45)

```matlab
methods (TestMethodSetup)
    function buildModels(testCase)
        % Model construction outside the measurement boundary; the
        % conditionally stable loop makes margin() warn on every call
        ws = warning('off', 'Control:analysis:MarginUnstable');
        testCase.addTeardown(@() warning(ws));
        testCase.p      = load_hw3_params();
        testCase.Grigid = build_plant_rigid(testCase.p);
        testCase.Gfull  = build_plant_full(testCase.p, 'ins');
        testCase.Wchain = build_tvc(testCase.p) * ...
            build_notch_filter(testCase.p.wBM, 0.002, 0.7, +1);
        [~, testCase.Trigid] = assemble_loop(testCase.Grigid, testCase.Kref);
        testCase.wind = load_wind_profile(testCase.p);
    end
end
```

- Righe 32-44: **tutta la costruzione dei modelli avviene nel setup**, cioe'
  **fuori dal confine di misura**. E' la decisione metodologicamente piu'
  importante del file, ed e' dichiarata nel commento (riga 33).

  Il motivo: `matlab.perftest` misura **solo** cio' che sta dentro il ciclo
  `while testCase.keepMeasuring` (o fra `startMeasuring` e `stopMeasuring`). Se
  `load_hw3_params()`, `build_plant_full()` e `build_tvc()` stessero dentro il
  ciclo, il tempo misurato includerebbe il caricamento del `.mat` LPV e la
  costruzione delle `tf`, che **non fanno parte del ciclo caldo reale**: nel
  Monte Carlo e nel tuner il plant e la catena attuatore sono costruiti **una
  volta sola**, fuori dal loop, e solo `assemble_loop` viene richiamata. Un
  benchmark che includesse il setup misurerebbe il ciclo sbagliato e le
  ottimizzazioni verrebbero fatte nel posto sbagliato.
- Righe 35-36: la solita muta della warning `Control:analysis:MarginUnstable`, che
  su questo anello condizionatamente stabile scatterebbe **a ogni singola
  chiamata** di `margin`. Qui la motivazione e' anche **prestazionale**: emettere
  una warning ha un costo non nullo, e in un ciclo `keepMeasuring` da migliaia di
  iterazioni contaminerebbe la misura con il tempo di formattazione del messaggio.
  Il ripristino e' via `addTeardown`, come altrove.
- `TestMethodSetup` (non `TestClassSetup`): i modelli sono ricostruiti **prima di
  ogni test**. Costa qualche decina di millisecondi in piu' ma garantisce che ogni
  benchmark parta da oggetti freschi, senza cache interne della Control System
  Toolbox eventualmente scaldate dal test precedente.

---

## `testAssembleLoopRigid` (righe 48-57) -- livello *unit*

```matlab
% fminsearch cost kernel: close the rigid loop (connect +
% getLoopTransfer + minreal)
G = testCase.Grigid;  K = testCase.Kref;
while testCase.keepMeasuring
    [L, T] = assemble_loop(G, K);
end
testCase.verifyEqual(order(T), 4);
testCase.verifyNotEmpty(L);
```

- Righe 52-54: **il ciclo `keepMeasuring` e' l'unita' di misura**. Il framework
  esegue il corpo del `while` ripetutamente (con dei *warmup* iniziali scartati),
  raccoglie campioni finche' non raggiunge una confidenza statistica prefissata, e
  riporta media e dispersione. E' il pattern giusto per operazioni **brevi**
  (~10 ms): una singola esecuzione sarebbe dominata dal rumore dello scheduler,
  dalla cache, dal JIT.
- Riga 51: le variabili sono **copiate in locali** prima del ciclo
  (`G = testCase.Grigid`) invece di essere lette dalle property dentro il ciclo.
  Non e' pignoleria: l'accesso a una property di un oggetto MATLAB passa per il
  meccanismo di dispatch della classe e costa qualcosa. Metterlo dentro il loop
  misurerebbe anche quello. E' un dettaglio idiomatico dei perf-test scritti bene.
- Righe 55-56: **le asserzioni stanno FUORI dal ciclo di misura**. Anche questo e'
  deliberato: `verifyEqual` e la macchina di diagnostica del framework hanno un
  costo, e se stessero dentro il `while` verrebbero cronometrate insieme al codice
  sotto test. Le due verifiche servono solo come **sanity check** (il benchmark ha
  effettivamente costruito l'oggetto giusto), non come test.
- `order(T) == 4`: l'anello chiuso rigido ha 4 stati. Il controllore e' un
  **guadagno statico** (`Kc = ss([...])` in `assemble_loop.m` riga 25, senza stati)
  e l'attuatore e' `tf(1)` (ideale): l'unica dinamica e' quella del plant.
- **Perche' questo e' il kernel giusto da misurare.** E' *esattamente* l'operazione
  che `fminsearch` ripete a ogni valutazione della funzione di costo
  (`design_controller.m` riga 82). Misurando questa, si misura il collo di
  bottiglia del tuner.

---

## `testAssembleLoopFullChain` (righe 59-67) -- livello *unit*, caso pesante

```matlab
G = testCase.Gfull;  K = testCase.Kref;  Wa = testCase.Wchain;
while testCase.keepMeasuring
    [L, T] = assemble_loop(G, K, Wa);
end
testCase.verifyEqual(order(T), 13);   % 6 + 5 (TVC+Pade) + 2 (notch)
```

- Riga 65: il conteggio degli stati e' una **verifica di composizione della
  catena**, e il commento esplicita la decomposizione:

      6  stati del plant completo (z, zdot, theta, thetadot, eta, etadot)
    + 5  stati della catena TVC (2 attuatore + 3 Pade)
    + 2  stati del notch (biquad)
    = 13

  Verificato: `order(T) = 13`. Se un `minreal` interno cancellasse qualcosa, o se
  l'ordine del Pade cambiasse, il conteggio non tornerebbe. E' lo stesso invariante
  di `hm3FilterTest.testTvcOrderIsActuatorPlusPade`, ma **a valle**, sull'anello
  assemblato: verifica che `connect` non abbia perso pezzi per strada.
- **Il risultato sorprendente**: l'anello completo a 13 stati costa **~9.1 ms**,
  cioe' *meno* dell'anello rigido a 4 stati (~10.5 ms). E' controintuitivo e vale
  la pena averlo notato. La spiegazione plausibile e' che il costo di
  `assemble_loop` sia dominato non dalla dimensione del sistema (13 contro 4 stati
  sono entrambi minuscoli per un solver denso) ma dall'**overhead fisso** di
  `connect` (parsing dei nomi, costruzione della matrice di interconnessione) e di
  `minreal`, e che la tolleranza `1e-6` di `minreal` faccia lavori diversi nei due
  casi. Comunque sia, e' proprio questo il genere di conclusione che un perf-test
  serve a produrre: **l'intuizione "piu' stati = piu' lento" e' sbagliata qui**, e
  ottimizzare riducendo l'ordine del modello non porterebbe da nessuna parte. Se
  si volesse davvero accelerare il Monte Carlo, il bersaglio sarebbe l'overhead di
  `connect`/`minreal` (per esempio precostruendo l'interconnessione una volta e
  aggiornando solo i guadagni), non il numero di stati.

---

## `testTunerCostEvaluation` (righe 69-80) -- livello *system*

```matlab
% One full cost evaluation as design_controller performs it:
% assemble_loop + margin + isstable on the rigid loop
G = testCase.Grigid;  K = testCase.Kref;
while testCase.keepMeasuring
    [L, T] = assemble_loop(G, K);
    [Gm, Pm] = margin(L);
    stable = isstable(T);
end
testCase.verifyTrue(stable);
testCase.verifyTrue(isfinite(Gm) && isfinite(Pm));
```

- Righe 73-77: misura **una valutazione completa della funzione di costo**, non solo
  il suo pezzo piu' costoso. E' il salto dal livello *unit* al livello *system*:
  quantifica il costo dell'unita' atomica che `fminsearch` ripete.
- Misurato: **~15.0 ms**, contro i ~10.5 ms del solo `assemble_loop`. Quindi
  l'assemblaggio pesa circa il 70 % del costo, e l'analisi dei margini il restante
  30 %. **Questa e' l'informazione azionabile**: se si volesse accelerare il tuner,
  bisogna guardare `assemble_loop`, non `margin`.
- Righe 78-79: sanity check fuori dal ciclo. `isfinite(Gm) && isfinite(Pm)` verifica
  che l'anello abbia effettivamente degli attraversamenti da misurare: se `margin`
  restituisse `Inf` (nessun crossing), la funzione di costo di `design_controller`
  sarebbe degenere e la misura non sarebbe rappresentativa.
- **Nota di onesta'.** Il test replica la funzione di costo usando `margin(L)`, ma
  `design_controller` usa in realta' `classify_margins(Lt, bands{:})` (riga 84 del
  sorgente), che a sua volta chiama `allmargin`. Non e' la stessa cosa: misurato,
  `classify_margins` costa ~3.7 ms, un valore simile a `margin`, quindi la stima
  complessiva regge. Ma il test **non sta misurando esattamente il codice che gira
  nel tuner**, e chiamarlo "one full cost evaluation as design_controller performs
  it" (riga 70) e' una lieve sovra-affermazione. Manca inoltre il `minreal(Lt, 1e-6)`
  che `design_controller` esegue alla riga 83.
- Coerenza dei numeri: ~15 ms per valutazione x ~55 valutazioni = ~830 ms, che e'
  esattamente il costo misurato di una ricerca completa. I conti tornano, il che
  conferma che il modello di costo e' corretto.

---

## `testGustResponseLsim` (righe 82-89) -- livello *system*

```matlab
% Time-domain replay: lsim over the 2401-point severe gust
T = testCase.Trigid;  w = testCase.wind;
while testCase.keepMeasuring
    r = simulate_gust_response(T, w);
end
testCase.verifySize(r.theta, [numel(w.t) 1]);
```

- Righe 85-87: misura la simulazione temporale. `simulate_gust_response` e' un
  wrapper attorno a un solo `lsim` (riga 21 del sorgente) su una griglia di **2401
  punti** (12 s a dt = 5 ms, piu' l'istante iniziale -- verificato:
  `numel(w.t) = 2401`, coerente con il commento della riga 83).
- Misurato: **~5.8 ms**. E' l'operazione **piu' economica** del gruppo, il che e'
  un risultato in se': la simulazione temporale, che intuitivamente sembra la cosa
  costosa (2401 passi di integrazione!), costa **meno della meta'** di un
  assemblaggio d'anello. Il motivo e' che `lsim` su un sistema lineare tempo
  invariante non integra numericamente: **discretizza una volta** (matrice
  esponenziale) e poi itera una ricorsione matrice-vettore, che su 4 stati e' quasi
  gratis. Il collo di bottiglia di HM3 e' l'**algebra dei sistemi**, non
  l'**integrazione**.
- Riga 88: `verifySize(r.theta, [numel(w.t) 1])` -- vettore colonna della lunghezza
  giusta. Sanity check, fuori dal ciclo.
- Perche' vale la pena misurarlo comunque: il Monte Carlo (`main_montecarlo.m`) e i
  corner del Task 3 producono anche **risposte temporali** per ogni caso, non solo
  margini. Con N = 1500 anche 5.8 ms diventano ~9 s.

---

## `testDesignControllerSearch` (righe 91-98) -- livello *workflow*

```matlab
% Workflow level: one complete PD margin-matching search
G = testCase.Grigid;
testCase.startMeasuring();
K = design_controller(G, [], 'verbose', false);
testCase.stopMeasuring();
testCase.verifyEqual(K.Kp_th, testCase.Kref.Kp_th, 'AbsTol', 5e-3);
```

- Righe 94-96: **cambia la API di misura**. Qui non c'e' `keepMeasuring` ma la
  coppia esplicita `startMeasuring` / `stopMeasuring`. La differenza e'
  sostanziale:
  - `keepMeasuring` ripete il corpo **molte volte** finche' la statistica non
    converge. Va bene per operazioni da millisecondi.
  - `startMeasuring`/`stopMeasuring` cronometra **una singola esecuzione** del
    codice fra le due chiamate. E' il pattern per operazioni **costose**, dove
    ripetere centinaia di volte sarebbe proibitivo: qui una ricerca completa costa
    **~830 ms**, quindi anche solo i warmup del framework portano il test a durare
    una quindicina di secondi (misurato: `runperf` su questo singolo test ha
    impiegato 15.7 s).
- **Perche' esiste il livello workflow.** I tre test precedenti misurano i
  *mattoni*; questo misura il *risultato*. Serve a rispondere a una domanda che i
  micro-benchmark non possono rispondere: "quante volte posso permettermi di
  ri-tarare il controllore?". La risposta -- 830 ms -- e' quella che rende
  praticabile il Task 3 (ri-tarare a ogni vertice del box di incertezza) e il Monte
  Carlo (dove pero' il controllore e' **fisso**, riga 42 di `main_montecarlo.m`:
  `K = design_controller(...)` viene chiamato **una volta sola**, fuori dal
  `parfor`, e poi congelato -- che e' anche il senso fisico dello studio di
  robustezza: il controllore non si adatta, subisce).
- **Riga 97: e' la riga che era rotta, e ora passa.** `testCase.Kref.Kp_th` vale
  `1.7845` e `design_controller` restituisce `1.784455`: lo scarto e' 4.5e-5 contro
  una `AbsTol` di 5e-3, comodamente dentro. Fino a poco fa `Kref` pinnava pero'
  `1.9800`, e il test falliva con uno scarto di 0.196 -- **invalidando la misura**
  (`0 Valid, 1 Invalid`), cioe' azzerando il benchmark di livello workflow. Vedi la
  sezione su `Kref` per la storia completa.
- Nota metodologica: usare una `verifyEqual` sui guadagni **dentro un perf-test** e'
  discutibile in se', e questo episodio lo dimostra. Un perf-test dovrebbe misurare;
  la correttezza e' compito di `hm3LoopTest`, che gia' pinna gli stessi guadagni
  (`testDesignControllerMeetsTargets`, righe 84-85). La duplicazione qui non aggiunge
  copertura e ha aggiunto un modo di rompersi -- **quello che poi si e' rotto davvero**.
  Il valore residuo dell'asserzione e' garantire che il benchmark non stia
  cronometrando una ricerca **divergente** (un `fminsearch` che sbatte contro
  `MaxFunEvals` senza convergere costerebbe di piu' e falserebbe la misura): un
  controllo piu' robusto e meno fragile sarebbe verificare i **margini**
  (`m.stable`, `|GM| ~ 6`, `PM ~ 30`) invece dei guadagni, perche' i margini sono i
  **requisiti della traccia** e non si spostano quando il progetto si ri-tara.

---

## Cosa questa suite NON fa (limiti noti)

E' il punto piu' importante da saper dire all'orale, perche' un perf-test mal
capito da' un falso senso di sicurezza.

- **Non c'e' nessuna soglia.** Nessun test asserisce "questa operazione deve durare
  meno di X ms". `runperf` **misura e riporta**, non boccia. Cioe' se domani
  `assemble_loop` diventasse 10 volte piu' lenta, **tutti i test passerebbero
  comunque**: il rallentamento comparirebbe solo nella
  tabella dei tempi, e solo se qualcuno la guarda. Un perf-test senza soglia e' una
  **misura**, non un **test di regressione**. Per farne un gate servirebbe salvare
  una baseline e confrontarsi con quella (`matlab.perftest` lo supporta, ma qui non
  e' usato).
- **Il Monte Carlo (N = 1500) non e' benchmarkato**, pur essendo il ciclo caldo
  citato nel commento di testa. Il perf-test misura il **mattone** (`assemble_loop`)
  ma non il **ciclo** che lo contiene, quindi non cattura l'overhead del `parfor`,
  del trasferimento dei dati ai worker, ne' la degradazione a seriale in assenza di
  Parallel Computing Toolbox.
- **Il trade dei filtri del Task 2 non e' benchmarkato**, idem.
- **`classify_margins` non e' misurata direttamente**, pur essendo dentro la
  funzione di costo del tuner (`testTunerCostEvaluation` usa `margin` al suo posto).
- **Non si misura la memoria**, solo il tempo.
- Piccola tensione con la convenzione del repo: il `CLAUDE.md` prescrive che le
  funzioni di ciclo caldo **non** portino validazione `arguments`, perche' stanno
  dentro loop da milioni di chiamate. `assemble_loop.m` **ha** un blocco `arguments`
  (righe 15-19) ed **e'** una funzione di ciclo caldo. In questo caso la deroga e'
  innocua -- il costo di `mustBeA` e' del tutto trascurabile rispetto ai ~10 ms di
  `connect` + `minreal` -- ma vale la pena notare che **e' proprio questo perf-test
  lo strumento che permette di affermarlo con dati invece che per fede**.

---

## Il commento "~400x" (righe 3-4)

```matlab
%  Unit: one assemble_loop call (fminsearch hits it ~400x, the Task-2
%  sweeps ~150x). System: full tuner cost (assemble_loop + margin +
%  isstable) and the lsim gust replay. Workflow: one design_controller
%  search. Run: runperf('hm3LoopPerformanceTest')
```

Il "~400x" e' il **limite superiore**, non il valore tipico: `400` e' il
`MaxFunEvals` passato a `optimset` in `design_controller.m` (riga 57). Misurando,
una ricerca completa costa ~830 ms a ~15 ms per valutazione, cioe' **~55
valutazioni effettive**: `fminsearch` converge molto prima del tetto, grazie al
buon seed di D'Antuono e alle tolleranze (`TolX = 1e-4`, `TolFun = 1e-3`). Il
commento non e' sbagliato come *worst case*, ma va letto per quello che e'.

Questo, incidentalmente, e' un buon argomento a favore del seed analitico: partire
dalla coppia canonica `Kp0 = 2*A_6/K_1`, `Kd0 = sqrt(A_6)/K_1` invece che da un
punto arbitrario riduce di un ordine di grandezza il numero di valutazioni. **Il
seed teorico non e' eleganza, e' performance.**

---

## Possibili domande d'esame

**D: Perche' hai scritto dei performance test per un homework di controlli? Non
basta che il risultato sia giusto?**
R: Perche' HM3 non calcola un anello, ne calcola migliaia. `assemble_loop` viene
chiamata fino a 400 volte per ogni ricerca dell'auto-tuner, ~150 volte nel trade
dei filtri del Task 2, e **1500 volte** nel Monte Carlo (una per estrazione, dentro
il `parfor` di `main_montecarlo.m`). Misurata, costa ~10 ms a chiamata: solo il
Monte Carlo sono ~14 s di assemblaggio, piu' i margini. Se qualcuno la rallentasse
di un fattore 3 -- aggiungendo un `minreal` di troppo, per dire -- il Monte Carlo
diventerebbe impraticabile e me ne accorgerei nel momento sbagliato. Il perf-test
quantifica il costo dei mattoni e mi dice **dove** ottimizzare se serve: ho
misurato che l'assemblaggio pesa il 70 % di una valutazione della funzione di costo
e l'analisi dei margini il 30 %, quindi il bersaglio e' `connect`/`minreal`, non
`margin`.

**D: Perche' un test usa `keepMeasuring` e un altro `startMeasuring`/`stopMeasuring`?**
R: Dipende dal costo dell'operazione. `keepMeasuring` ripete il corpo del ciclo
molte volte, scartando dei warmup, finche' la misura non e' statisticamente
affidabile: e' il pattern giusto per operazioni brevi (~10 ms), dove una singola
esecuzione sarebbe sepolta dal rumore di scheduler e cache. `startMeasuring` /
`stopMeasuring` cronometra invece **una singola esecuzione**, ed e' il pattern per
operazioni costose: una ricerca completa di `design_controller` dura ~830 ms, e
ripeterla centinaia di volte per convergenza statistica renderebbe il test
inutilizzabile (gia' cosi', con i soli warmup, il test impiega ~16 s). In entrambi i
casi le `verify*` stanno **fuori** dal confine di misura, altrimenti cronometrerei
anche il framework di test.

**D: Che cosa hai scoperto misurando, che non avresti indovinato?**
R: Due cose. Primo, che l'anello **completo** a 13 stati (plant + TVC + Pade +
notch) si assembla in ~9 ms, cioe' **piu' velocemente** dell'anello rigido a 4
stati (~10.5 ms). L'intuizione "piu' stati = piu' lento" e' falsa qui: il costo e'
dominato dall'overhead fisso di `connect` e `minreal`, non dalla dimensione del
sistema. Ridurre l'ordine del modello non accelererebbe nulla. Secondo, che
`lsim` sulla raffica da 2401 punti costa ~5.8 ms, **meno della meta'** di un
assemblaggio d'anello: la simulazione temporale non e' il collo di bottiglia, perche'
su un sistema LTI `lsim` discretizza una volta e poi itera una ricorsione
matrice-vettore, non integra numericamente. Il collo di bottiglia di HM3 e'
l'algebra dei sistemi, non l'integrazione.

**D: I tuoi perf-test possono fallire se il codice rallenta?**
R: No, ed e' un limite che riconosco. Nessuno dei test asserisce una **soglia** di
tempo: `runperf` misura e riporta, non boccia. Se `assemble_loop` diventasse dieci
volte piu' lenta, i test passerebbero comunque e il rallentamento comparirebbe solo
nella tabella dei tempi. Per farne un vero gate di regressione servirebbe salvare
una **baseline** e confrontarsi con quella (`matlab.perftest` lo supporta). Allo
stato, questa suite e' uno strumento di **profiling riproducibile**, non un test di
non-regressione prestazionale.

**D: C'e' stato qualcosa che non andava in questa suite?**
R: Si', ed e' l'aneddoto che racconto piu' volentieri perche' e' un errore
*strutturale*, non di distrazione. La property `Kref` di `hm3LoopPerformanceTest`
pinnava `Kp_th = 1.9800, Kd_th = 1.3997` -- i guadagni **pre-ritaratura** -- mentre
il commento sosteneva di essere "same as hm3LoopTest", che invece pinnava
`1.7845 / 0.4433`. I due non coincidevano, e `design_controller` produce
`1.784455`. Risultato: `testDesignControllerSearch` **falliva**, e siccome in
`matlab.perftest` una verifica fallita **invalida la misura**, il benchmark di
livello workflow non produceva piu' alcun dato (`0 Valid, 1 Invalid`) -- cioe' il
numero piu' importante della suite era silenziosamente sparito. Era uno *stale
golden value* nato dalla **duplicazione dello stesso numero in due file**: quando ho
ri-tarato il progetto sull'anello completo (passando da `margin()` a
`classify_margins`) ho aggiornato il test unitario e non il perf-test, e non ho
rilanciato la suite per intero. Ed e' rimasto nascosto a lungo perche' gli altri
quattro test usano `Kref` solo per **cronometrare** l'anello, e il tempo non dipende
dal valore dei guadagni. Non e' che il perf-test sfugga a `runtests`:
`matlab.perftest.TestCase` deriva da `matlab.unittest.TestCase`, quindi
`runtests('HM3/tests')` -- il comando prescritto dalle convenzioni del repo --
scopre anche questa classe e il fallimento sarebbe emerso li'. **Ora e' allineato e
la suite HM3 passa 35/35.** La correzione vera pero' non e' il numero: e' non
duplicarlo affatto, derivandolo da un'unica sorgente. La lezione generale e' che un
golden value e' la **fotografia di un progetto**, e invecchia esattamente quando il
progetto cambia -- cioe' quando sei meno propenso a rileggere i test. Per questo,
in un perf-test, avrei fatto meglio ad asserire i **margini** (che sono i requisiti
della traccia e non si spostano) invece dei guadagni (che sono un risultato e si
spostano).
