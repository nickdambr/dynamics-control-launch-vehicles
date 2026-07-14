# HM0_falcon9_ascent/tests/falcon9AscentTest.m

## Ruolo del file nel progetto

E' la suite di test di HM0: una classe `matlab.unittest.TestCase` che esegue
**entrambi** gli script (`main.m` dimensionale e `main2.m` adimensionale) una sola
volta, ne raccoglie i risultati e verifica cinque proprieta'. Si lancia con
`runtests('HM0_falcon9_ascent/tests')`.

La sua strategia e' interessante e vale la pena capirla bene, perche' e'
condizionata da un vincolo che il codice sotto test le impone: **gli script HM0
non sono funzioni**. Non hanno una firma, non restituiscono nulla, non accettano
parametri, e cominciano con `clear` (riga 6 di `main.m`, riga 35 di `main2.m`).
Non c'e' quindi nessun modo "pulito" di chiamarli e leggerne l'output: bisogna
eseguirli e poi **raccogliere le variabili dal workspace**. Tutto il design della
classe discende da qui.

La suite e' quindi un **test di integrazione end-to-end**, non un insieme di unit
test. Nessuna delle funzioni interne (`eom`, `eom_nd`, `arc_boundary_events`) e'
testabile singolarmente: sono *local function* dentro gli script, quindi non
raggiungibili dall'esterno. Il test puo' solo osservare il *risultato finale* della
propagazione e chiedersi se e' plausibile.

Il fulcro concettuale e' il **test di cross-validazione** (righe 54-61): due
implementazioni indipendenti, con scalature e parametrizzazioni temporali diverse,
devono produrre gli stessi numeri. E' l'oracolo piu' forte disponibile in assenza
di una soluzione analitica.

---

## Dichiarazione della classe e docstring (righe 1-8)

```matlab
classdef falcon9AscentTest < matlab.unittest.TestCase
    %  Note: running the scripts regenerates the PNGs in figures/ (the
    %  repository pipeline owns those files, so this is intended).
```

- Riga 1: classe di test class-based (non function-based), coerente con la
  convenzione della repo (`CLAUDE.md`).
- Righe 7-8: la docstring avvisa di un **effetto collaterale**: eseguire i test
  **riscrive i PNG in `figures/`**. E' un'onesta' apprezzabile, ma va precisata:
  solo `main.m` esporta i PNG (righe 352-369 di quel file); `main2.m` **non ha il
  blocco di export**. Quindi l'avviso e' corretto solo a meta'.
  Il fatto che una suite di test *scriva su disco nel repository* e' comunque una
  scelta discutibile: idealmente un test dovrebbe essere osservativo. Qui e' stato
  accettato consapevolmente perche' i PNG sono comunque il prodotto della pipeline.

---

## `properties` (righe 10-13)

```matlab
properties
    dim   % results harvested from main.m
    nd    % results harvested from main2.m
end
```

- Due proprieta' che ospitano le due struct di risultati. Sono popolate una volta
  sola nel `TestClassSetup` e poi lette (mai scritte) dai test.
- I nomi `dim` / `nd` **non sono cosmetici**: coincidono con i valori del
  `TestParameter` `impl` alla riga 16, e questo e' cio' che permette
  l'indicizzazione dinamica `testCase.(impl)` alla riga 31. Se si rinominasse una
  proprieta' senza rinominare il parametro, i test parametrizzati fallirebbero con
  un errore oscuro di "campo inesistente".

---

## `properties (TestParameter)` (righe 15-17)

```matlab
properties (TestParameter)
    impl = {'dim', 'nd'}
end
```

- Questo e' il meccanismo di **parametrizzazione** di `matlab.unittest`. Ogni
  metodo di test che dichiara `impl` fra i suoi argomenti viene **eseguito una
  volta per ciascun valore** del cell array.
- Effetto concreto sul conteggio: quattro metodi parametrizzati (righe 29, 36, 42,
  47) x due implementazioni = **8 punti di test**, piu' il test di
  cross-validazione non parametrizzato (riga 54) = **9 test totali**. Ciascuno
  appare separatamente nel report, con nomi tipo
  `testPropellantBookkeeping(impl=nd)`.
- Il guadagno e' che **le stesse quattro proprieta' fisiche sono verificate su
  entrambe le implementazioni senza duplicare una riga di codice**. Se domani
  arrivasse un `main3.m`, basterebbe aggiungere `'x'` al cell array e una
  proprieta' `x`.

---

## `TestClassSetup` -- `runScripts` (righe 19-26)

```matlab
methods (TestClassSetup)
    function runScripts(testCase)
        hm0 = fileparts(fileparts(mfilename('fullpath')));
        testCase.dim = runAscentScript(fullfile(hm0, 'main.m'));
        testCase.nd  = runAscentScript(fullfile(hm0, 'main2.m'));
        testCase.addTeardown(@() close('all'));
    end
end
```

- Riga 19: `TestClassSetup` (non `TestMethodSetup`) significa **una sola
  esecuzione per l'intera classe**, non una per test. E' una scelta obbligata dal
  costo: ogni script fa un'integrazione `ode45` con `RelTol = 1e-10` e genera 8-9
  figure. Con `TestMethodSetup` il setup -- che esegue **entrambi** gli script --
  girerebbe una volta per ciascuno dei 9 punti di test: ogni script girerebbe
  **9 volte** invece di 1, cioe' **18 esecuzioni** di script invece di 2.
  Il prezzo di questa scelta e' che i test **condividono lo stato**: se un test
  modificasse `testCase.dim`, contaminerebbe i successivi. Qui non succede (i test
  sono tutti in sola lettura), ma e' una fragilita' latente.
- Riga 21: `fileparts(fileparts(mfilename('fullpath')))` risale di **due** livelli:
  da `.../HM0_falcon9_ascent/tests/falcon9AscentTest.m` al file, poi a `tests/`,
  poi a `HM0_falcon9_ascent/`. Il doppio `fileparts` e' il modo idiomatico di
  ottenere la cartella genitore. Cosi' il test **non dipende dal working
  directory**: si puo' lanciare da qualunque parte della repo.
- Righe 22-23: le due chiamate alla funzione locale `runAscentScript` -- il cuore
  del pattern, discusso sotto.
- Riga 24: `addTeardown(@() close('all'))` chiude le figure **dopo che tutti i
  test della classe sono finiti**. Senza questo, una `runtests` lascerebbe aperte
  le **9 finestre di `main2.m`** (non 8 + 9 = 17: `main2.m` comincia alla riga 35
  con `clear; close all; clc;`, quindi il suo `close all` chiude le 8 figure
  appena create da `main.m` prima di aprire le proprie 9 -- lo stesso `close all`
  di `main.m` e' cio' che garantisce che il suo export loop, basato su
  `findobj(groot,'Type','figure')` alla riga 359, esporti esattamente 8 PNG).
  Registrato in `TestClassSetup`, il teardown ha durata di classe, coerentemente
  con il setup.

---

## `runAscentScript` -- il pattern "run + harvest" (righe 65-74)

E' la funzione piu' importante del file, ed e' l'unico modo per testare degli
script MATLAB che cominciano con `clear`.

```matlab
function S = runAscentScript(scriptPath)
%  The script's leading `clear` wipes this workspace first, then
%  the script repopulates it; results are collected after the run.
    run(scriptPath);
    S = struct('t', t, 'h', h, 'mass', mass, 'Vmag', Vmag, ...
               'qdyn', qdyn, 'Mach', Mach, 'qmax', qmax, ...
               'm0', m0, 'Qdot1', Qdot1, 'tb1', tb1, ...
               'im1', im1, 'imQ', imQ);
end
```

### Perche' non si puo' chiamare `run()` direttamente in un metodo di test

Questo e' **la domanda d'esame** su questo file, e la risposta va capita a fondo.

`run(scriptPath)` non e' una funzione "normale": esegue lo script **nel workspace
del chiamante**. E' proprio cio' che serve -- altrimenti le variabili prodotte dallo
script sarebbero irraggiungibili -- ma ha una conseguenza brutale.

Immaginiamo di scrivere:

```matlab
function testQualcosa(testCase)   % NON funziona
    run('main.m');                % <- la riga 6 di main.m fa `clear`
    testCase.verifyGreaterThan(h(end), 50e3);   % <- ERRORE
end
```

La prima istruzione eseguita da `main.m` e' `clear` (riga 6). Quel `clear` viene
eseguito **nel workspace del metodo di test**, e cancella **tutto** cio' che c'e'
dentro -- **`testCase` compreso**. Alla riga successiva, `testCase` non esiste piu':
si ottiene un errore di variabile non definita. Non c'e' modo di aggirarlo dal
metodo: qualunque cosa si salvi in una variabile locale prima della `run` viene
spazzata via.

**La soluzione e' il sacrificio di un workspace.** `runAscentScript` e' una
funzione la cui *unica* variabile locale e' `scriptPath` -- che serve solo a `run` e
che, dopo, non serve piu' a nessuno. Il `clear` dello script distrugge quel
workspace usa-e-getta (`scriptPath` incluso!), lo script lo ripopola con le sue
variabili, e la riga 70 le raccoglie in una struct che viene **restituita per
valore** al chiamante. Il workspace del metodo di test, con dentro `testCase`,
non e' mai stato esposto.

In sintesi: **la funzione locale fa da "camera di decontaminazione"** fra il
`clear` dello script e lo stato del test framework. E' esattamente il pattern che
`CLAUDE.md` documenta come intenzionale per HM0 ("the scripts' leading `clear`
makes assertions on `testCase` impossible inside the same workspace, hence the
helper pattern").

### Il prezzo del pattern

- **Accoppiamento per nome, non per contratto.** La riga 70 pesca 12 variabili
  (`t`, `h`, `mass`, `Vmag`, `qdyn`, `Mach`, `qmax`, `m0`, `Qdot1`, `tb1`, `im1`,
  `imQ`) **per nome** dal workspace dello script. Non c'e' nessuna interfaccia,
  nessun contratto, nessun controllo statico. Se domani si rinominasse `qmax` in
  `q_max` dentro `main.m`, il test esploderebbe con un errore di variabile non
  definita, non con un fallimento leggibile. Gli script sono di fatto
  **un'API implicita** definita dai nomi delle loro variabili interne.
- **Fragilita' rispetto ai due script.** Le stesse 12 variabili devono esistere
  con gli **stessi nomi** in `main.m` **e** in `main2.m`. Cosa che al momento e'
  vera (verificato: `main2.m` produce `t` alla riga 200, `h` alla 212, `mass` alla
  209, `Vmag` alla 215, `qdyn` alla 224, `Mach` alla 225, `qmax`/`imQ` alla 250,
  `im1` alla 249, `m0` alla 100, `Qdot1` alla 67, `tb1` alla 64), ma e' un
  invariante non dichiarato e non verificato.
- **Lentezza ed effetti collaterali.** Ogni chiamata esegue l'integrazione
  completa, disegna 8-9 figure e (per `main.m`) scrive 8 PNG su disco. Nessuna di
  queste cose serve ai test.
- `MException` / `warning` degli script si propagano al test come errori di setup.

**L'alternativa vera**, se un giorno si volesse una suite pulita: rifattorizzare
la fisica in una funzione (`[t, Y] = simulate_ascent(par)`), lasciare agli script
solo il ruolo di driver + plotting, e testare la funzione. In quel caso `eom`
diventerebbe testabile isolatamente (per esempio verificando che a velocita' nulla
il RHS dia esattamente `du = -mu/r^2 + T_SL/m0`). Il pattern attuale e' una
**risposta corretta al vincolo dato**, non un design ideale.

> **Possibile domanda d'esame** -- Perche' non basta salvare `testCase` in una
> variabile temporanea prima di `run()` per proteggerla dal `clear`?
> *Risposta:* Perche' `clear` senza argomenti cancella **l'intero workspace**, non
> una lista selettiva: qualunque copia, con qualunque nome, viene distrutta.
> L'unica difesa e' che il `clear` non venga mai eseguito **nello stesso
> workspace** in cui vive `testCase`. Da qui la funzione locale: si sacrifica un
> workspace vuoto e monouso, e si passano i risultati indietro **per valore**
> tramite il return della funzione -- che il `clear` non puo' toccare, perche'
> avviene prima.

---

## `testPropellantBookkeeping` (righe 29-34)

```matlab
function testPropellantBookkeeping(testCase, impl)
    % dm/dt = -Qdot is constant: m(tb) = m0 - Qdot*tb exactly
    S = testCase.(impl);
    testCase.verifyEqual(S.mass(end), S.m0 - S.Qdot1*S.tb1, 'RelTol', 1e-8);
end
```

- Riga 29: la firma dichiara `impl` come secondo argomento -- e' cio' che attiva la
  parametrizzazione. Riga 31: `testCase.(impl)` e' indicizzazione dinamica di
  proprieta' (`testCase.dim` oppure `testCase.nd`).
- Righe 32-33: **questo e' l'unico test con un oracolo analitico esatto**, ed e' il
  piu' potente della suite nonostante l'apparenza banale. L'equazione di massa e'
  `dm/dt = -Qdot` con `Qdot` **costante** (riga 463 di `main.m`, riga 574 di
  `main2.m`): la portata non dipende ne' dalla quota ne' dallo stato, perche'
  l'ugello e' in condizioni critiche. La soluzione esatta e' quindi

      m(t) = m0 - Qdot * t     =>     m(t_b) = m0 - Qdot * t_b

  con `t_b = 162` s. Non c'e' nessuna approssimazione: e' la soluzione vera, non
  una stima.
- Cosa verifica davvero questo test, componente per componente:
  - **Per `main.m`**: che `ode45` propaghi correttamente per 162 s. La tolleranza
    `RelTol 1e-8` e' un controllo diretto sul fatto che le tolleranze
    dell'integratore (`RelTol 1e-10`) stiano facendo il loro lavoro. Se qualcuno
    allentasse `RelTol` a 1e-4, questo test lo prenderebbe.
  - **Per `main2.m`**: molto di piu'. Questo test e' l'unico controllo automatico
    sulla **correttezza dell'intera riparametrizzazione in tau**. Perche' `m*(t*)`
    arrivi al valore giusto devono essere corrette **tutte** queste cose insieme:
    (a) `Qdot_nd = Qdot*T_ref/m_ref` (riga 121); (b) la scalatura chain-rule
    `dydt = Delta * [...]` applicata **anche alla componente di massa** (riga 577);
    (c) le tre `Delta_k` (righe 145-147); (d) la mappa `tau -> t*` dentro `eom_nd`
    (righe 480-496); (e) la mappa inversa in post-processing (righe 195-197); (f)
    la ri-dimensionalizzazione `mass = Y_nd(:,7)*m_ref` (riga 209).
    Sbagliare **una qualunque** di queste fa fallire il test. E' un ottimo
    "canarino" per la parametrizzazione a tre archi, e vale molto piu' di quanto la
    sua semplicita' suggerisca.

> **Possibile domanda d'esame** -- Verificare che `m(t_b) = m0 - Qdot*t_b` sembra
> tautologico: e' un'equazione lineare, cosa c'e' da testare?
> *Risposta:* Per `main.m` e' effettivamente poco piu' di un controllo
> sull'accuratezza dell'integratore. Per `main2.m` invece e' il test **piu'
> informativo** della suite: e' l'unico oracolo **esatto** disponibile, e per
> passare richiede che siano simultaneamente corretti l'adimensionalizzazione
> della portata, la scalatura chain-rule `Delta_k` applicata anche alla componente
> di massa, le tre durate d'arco, la mappa tau -> t dentro l'EOM, la mappa inversa
> in post-processing e la ri-dimensionalizzazione finale. Un errore in una sola di
> queste sei cose lo fa fallire. Testare la componente **piu' semplice** del
> sistema e' spesso il modo migliore per smascherare errori nella **struttura**.

---

## `testMachOneBeforeMaxQ` (righe 36-40)

```matlab
S = testCase.(impl);
testCase.verifyNotEmpty(S.im1);
testCase.verifyLessThan(S.t(S.im1), S.t(S.imQ));
```

- Riga 38: `verifyNotEmpty(S.im1)` non e' pignoleria. `im1` e' prodotto da
  `find(Mach >= 1, 1, 'first')`, che restituisce `[]` se il veicolo non supera
  Mach 1. Senza questo controllo, la riga 39 farebbe `S.t([])` -> array vuoto, e
  `verifyLessThan([], ...)` **passerebbe silenziosamente** (un confronto su array
  vuoto e' vacuamente vero). Il test darebbe verde su una simulazione in cui il
  razzo non decolla. La riga 38 chiude questo buco: e' esattamente il tipo di
  dettaglio che distingue un test scritto bene da uno scritto in fretta.
- Riga 39: verifica l'**ordinamento** degli eventi di missione, Mach 1 prima di
  max-Q. E' un controllo **qualitativo di fisica**, non un confronto numerico: non
  fissa i valori (61.8 s e 74.9 s nella run nominale), fissa la *relazione*.
  Questo lo rende **robusto ai cambiamenti di modello**: se domani si introducesse
  un CD dipendente da Mach, i due istanti cambierebbero ma l'ordine no, e il test
  continuerebbe a essere valido.
- **Limite**: e' un test debole. Un bug che spostasse max-Q da 74.9 s a 150 s lo
  supererebbe tranquillamente. Verifica una condizione necessaria, non
  sufficiente.

---

## `testAltitudeStaysAboveGround` (righe 42-45)

```matlab
testCase.verifyGreaterThanOrEqual(min(S.h), -1e-6);
```

- Riga 44: la quota non deve mai essere negativa (il razzo non scava). La soglia
  `-1e-6` m (un micrometro) e' una **tolleranza di floating point**, non una
  concessione fisica: a t = 0 vale `h = r0 - RE = RE - RE = 0` esattamente, ma i
  passi successivi dell'integratore possono produrre, per errore di
  arrotondamento, valori dell'ordine di -1e-9 m. Senza la tolleranza, il test
  fallirebbe su un artefatto numerico.
- E' un **sanity check strutturale**: cattura un errore di segno nella gravita',
  nella spinta, o nella componente `du` -- cioe' proprio le classi di bug piu'
  probabili in un file di EOM. Se `g_u` avesse segno sbagliato, o se il termine
  centrifugo fosse `-(v^2+w^2)/r`, il razzo sprofonderebbe e il test lo direbbe
  subito.
- **Copre anche un caso limite reale**: il rapporto spinta/peso al decollo e'
  1.275, non enorme. Se qualcuno sbagliasse la correzione di contropressione (per
  esempio **sommando** invece di sottrarre `p*A_ex`, o dimenticandola del tutto e
  usando T_vac al livello del mare), il T/W cambierebbe e questo test in generale
  **non** lo prenderebbe (il razzo decollerebbe comunque, anzi meglio). E' un buco:
  la correzione di contropressione, che e' il pezzo di propulsione piu'
  interessante del codice, **non e' testata da nessuno dei cinque test**.

---

## `testFinalAltitudeInPlausibleBand` (righe 47-52)

```matlab
% Loose sanity band for the Falcon 9 first-stage MECO altitude
testCase.verifyGreaterThan(S.h(end), 50e3);
testCase.verifyLessThan(S.h(end), 120e3);
```

- Righe 50-51: banda 50-120 km sulla quota a MECO. Il valore nominale (README) e'
  **82.84 km**, quindi il test ha ~33 km di margine sotto e ~37 km sopra.
- E' un **test di regressione volutamente lasco**. Il commento alla riga 48 lo
  ammette apertamente ("Loose sanity band"). La banda e' larga abbastanza da
  sopravvivere a un cambio di modello ragionevole (un CD dipendente da Mach, un
  profilo di atmosfera diverso) ma stretta abbastanza da catturare un disastro (un
  errore di unita', una spinta sbagliata di un fattore, un `Qdot` errato).
- **Limite dichiarato**: non e' un confronto con un valore di riferimento della
  traccia, e non stringe abbastanza da catturare regressioni sottili. Un bug che
  spostasse la quota MECO da 82.8 a 95 km passerebbe indisturbato. La banda e'
  onesta su cio' che e': una rete di sicurezza, non un test di accuratezza.

---

## `testDimensionalVsNonDimensionalAgreement` (righe 54-61)

```matlab
% The two scripts integrate the same physics with different
% scalings: results must agree to well below 1%
testCase.verifyEqual(testCase.nd.h(end),    testCase.dim.h(end),    'RelTol', 5e-3);
testCase.verifyEqual(testCase.nd.Vmag(end), testCase.dim.Vmag(end), 'RelTol', 5e-3);
testCase.verifyEqual(testCase.nd.mass(end), testCase.dim.mass(end), 'RelTol', 1e-6);
testCase.verifyEqual(testCase.nd.qmax,      testCase.dim.qmax,      'RelTol', 1e-2);
```

- Riga 54: **l'unico test non parametrizzato** -- per costruzione, perche' li
  confronta entrambi.
- E' il test **concettualmente piu' forte** della suite. In assenza di soluzione
  analitica, l'oracolo migliore e' una **seconda implementazione indipendente**.
  `main.m` e `main2.m` non sono banalmente lo stesso codice: differiscono per
  scalatura di tutte le variabili, parametrizzazione del tempo (t contro tau),
  griglia dell'integratore, e persino per come e' scritta la spinta in fase 2
  (`cos(deg2rad(90))` contro `0` esatto). Che convergano agli stessi numeri e' una
  verifica sostanziale.

**La struttura delle quattro tolleranze e' il punto piu' interessante del file** --
non sono scelte a caso, e ognuna dice qualcosa:

- **`mass(end)` a `RelTol 1e-6`** (riga 59) -- la piu' **stretta**. Perche'? Perche'
  la massa e' l'unica componente con **soluzione esatta** (lineare in t). Entrambi
  gli script devono arrivare allo stesso `m0 - Qdot*tb` indipendentemente da come
  hanno parametrizzato il tempo. Non c'e' errore di modello, solo errore di
  integrazione, che con `RelTol 1e-10` e' minuscolo. Una tolleranza lasca qui
  sarebbe uno spreco: si perderebbe potere diagnostico.
- **`h(end)` e `Vmag(end)` a `RelTol 5e-3`** (righe 57-58) -- 0.5%. Queste sono
  **quantita' integrate**: dipendono da tutta la storia della traiettoria. I due
  script hanno **griglie temporali diverse** (`main.m` ha `MaxStep = 1` s e integra
  in t; `main2.m` non ha `MaxStep` e integra in tau con due discontinuita' del RHS
  ai confini d'arco), quindi accumulano errori di troncamento locale diversi.
  Inoltre c'e' la differenza reale della componente `Tw` in fase 2. Lo 0.5% e' la
  banda che assorbe tutto questo. Su 82.84 km significa ~400 m di tolleranza.
- **`qmax` a `RelTol 1e-2`** (riga 60) -- la **piu' lasca**, 1%. E la ragione e'
  importante e sottile: `qmax` non e' una quantita' integrata, e' un **massimo
  campionato**. Entrambi gli script lo calcolano con `[qmax, imQ] = max(qdyn)`,
  cioe' prendendo il valore piu' alto **fra i campioni della griglia di output**,
  senza interpolare. Il vero massimo cade fra due campioni; il valore riportato e'
  quindi sistematicamente **sottostimato**, e di quanto dipende **da dove capitano
  i campioni**. Poiche' i due script hanno griglie diverse (in t contro in tau), i
  loro campioni cadono in punti diversi attorno al picco e "tagliano" il massimo in
  modo diverso. La tolleranza dell'1% non e' un errore di fisica: e' l'errore di
  **quantizzazione della griglia**. Questa e' esattamente la ragione per cui
  `main.m` ha `MaxStep = 1` s.

> **Possibile domanda d'esame** -- Perche' la tolleranza su `qmax` (1%) e' dieci
> volte piu' lasca di quella su `h(end)` (0.5%) e diecimila volte piu' lasca di
> quella su `mass(end)` (1e-6)?
> *Risposta:* Perche' le tre quantita' hanno **fonti di errore di natura diversa**.
> `mass(end)` ha soluzione analitica esatta: l'unico errore e' quello
> dell'integratore, minuscolo con `RelTol 1e-10`. `h(end)` e `Vmag(end)` sono
> integrate su tutta la traiettoria e risentono delle diverse griglie e
> discontinuita' dei due script: 0.5% e' il margine per l'errore accumulato.
> `qmax` invece e' un **massimo campionato senza interpolazione**: il suo valore
> dipende da *dove cadono i campioni* rispetto al picco, e i due script hanno
> griglie diverse. L'1% assorbe l'errore di quantizzazione della griglia, non un
> errore fisico. Se si volesse stringere, bisognerebbe interpolare il picco (per
> esempio con una parabola sui tre campioni attorno al massimo) invece di prendere
> il campione piu' alto.

---

## Cosa la suite NON verifica (limiti dichiarati)

Vale la pena avere in tasca anche i buchi, perche' all'orale la domanda "e cosa
non hai testato?" e' sempre in agguato.

1. **La correzione di contropressione `T = T_vac - p*A_ex`** non e' testata da
   nessuna parte. E' il pezzo di modellazione propulsiva piu' interessante del
   codice e passerebbe indenne anche se avesse il segno sbagliato (il razzo
   volerebbe *meglio* e la banda 50-120 km probabilmente reggerebbe).
2. **La legge di pitchover** non e' testata: nessun controllo che `gammaT` scenda
   da 90 a 89.5 deg fra t = 5 e t = 15 s.
3. **Il gravity turn** non e' testato: nessun controllo che la spinta sia
   effettivamente allineata a `Vrel` in fase 3. (E il grafico che sembrerebbe
   mostrarlo, Figura 7, e' una tautologia -- vedi `hm0_main.md`.)
4. **`eom` / `eom_nd` non sono unit-testate** e non lo *possono* essere: sono
   *local function* dentro gli script, non raggiungibili dall'esterno. Non si puo'
   verificare, per esempio, che a stato iniziale il RHS restituisca esattamente
   `du = 2.719 m/s^2`.
5. **Nessun valore di riferimento della traccia** e' verificato: max-Q = 29.5 kPa,
   Mach 1 a 61.8 s, quota MECO 82.84 km compaiono nel README ma nessun test li
   fissa. La suite verifica **relazioni** e **bande**, mai numeri.
6. **La event function inerte di `main2.m`** (vedi `hm0_main2.md`) non viene
   smascherata da nessun test -- e non potrebbe esserlo, dato che rimuoverla non
   cambia i risultati.
7. **Nessun `matlab.perftest`**: `CLAUDE.md` prevede benchmark
   `<feature>PerformanceTest.m` per ogni homework; HM0 non ne ha.

Questo non significa che la suite sia mal fatta. Fa esattamente cio' che si puo'
fare **dato il vincolo di testare degli script**: verifica invarianti globali,
ordinamenti fisici, bande di plausibilita' e la cross-validazione fra due
implementazioni. Per andare oltre bisognerebbe rifattorizzare la fisica in
funzioni pure -- che e' un cambio di design, non un test in piu'.

---

## Possibili domande d'esame

**D: Perche' i test non chiamano `run('main.m')` direttamente in un metodo, ma
passano per una funzione locale?**
R: Perche' `run()` esegue lo script **nel workspace del chiamante**, e la prima
istruzione di `main.m` (riga 6) e' `clear`. Chiamato dentro un metodo di test,
quel `clear` cancellerebbe l'intero workspace del metodo, **`testCase` compreso**:
la riga successiva darebbe errore di variabile non definita. Nessuna copia difensiva
salva, perche' `clear` senza argomenti cancella tutto. La soluzione e' isolare la
`run()` in una funzione locale il cui workspace e' **sacrificabile** (contiene solo
`scriptPath`): il `clear` distrugge quel workspace, lo script lo ripopola, la
funzione raccoglie le variabili in una struct e la **restituisce per valore**. Il
workspace del test, con `testCase`, non e' mai stato esposto.

**D: Il `TestClassSetup` esegue gli script una volta sola per tutta la classe.
Quali sono i pro e i contro?**
R: Pro: costo. Ogni script fa un'integrazione con `RelTol = 1e-10` e disegna 8-9
figure; con `TestMethodSetup` il setup girerebbe una volta per ciascuno dei 9 punti
di test, quindi ogni script girerebbe **9 volte** invece di 1 (**18 esecuzioni**
totali invece di 2), e la suite diventerebbe intollerabilmente lenta. Contro: i test **condividono lo stato**. Se un
test modificasse `testCase.dim`, contaminerebbe tutti quelli successivi, e l'ordine
di esecuzione diventerebbe significativo -- cosa che un test framework non
garantisce. Qui il rischio non si concretizza perche' tutti i test sono in sola
lettura, ma e' una fragilita' latente. E' il classico trade-off velocita' contro
isolamento, risolto a favore della velocita' con una scelta consapevole.

**D: Il test di accordo fra `main.m` e `main2.m` usa quattro tolleranze diverse
(1e-6, 5e-3, 5e-3, 1e-2). Come si giustificano?**
R: Rispecchiano tre **meccanismi d'errore diversi**. `mass(end)` a 1e-6: la massa
ha soluzione analitica esatta (lineare, `dm/dt = -Qdot` costante), quindi l'unico
errore possibile e' quello dell'integratore, minimo con `RelTol 1e-10`; stringere
qui massimizza il potere diagnostico. `h(end)` e `Vmag(end)` a 5e-3: sono quantita'
**integrate** su 162 s, e i due script hanno griglie temporali diverse (`main.m`
integra in t con `MaxStep = 1`; `main2.m` integra in tau, senza `MaxStep`, con due
discontinuita' del RHS ai confini d'arco), quindi accumulano errori di troncamento
diversi. `qmax` a 1e-2: e' un **massimo campionato senza interpolazione**, quindi il
suo valore dipende da dove cadono i campioni rispetto al picco -- un errore di
quantizzazione della griglia, non di fisica. Ogni tolleranza e' calibrata sul
meccanismo d'errore che le compete.

**D: `verifyNotEmpty(S.im1)` alla riga 38 sembra ridondante, visto che alla riga
successiva `im1` viene comunque usato. Serve davvero?**
R: Si', ed e' un dettaglio importante. `im1 = find(Mach >= 1, 1, 'first')`
restituisce `[]` se il veicolo non supera mai Mach 1. Senza la riga 38, la riga 39
valuterebbe `S.t([])`, che e' un **array vuoto**, e `verifyLessThan([], [])`
**passa** -- perche' un'asserzione su un array vuoto e' vacuamente vera. Il test
darebbe quindi verde su una simulazione in cui il razzo non decolla nemmeno. La
riga 38 trasforma un falso negativo silenzioso in un fallimento esplicito.

**D: Qual e' il test piu' informativo della suite, e perche'?**
R: Non e' quello che sembra. `testPropellantBookkeeping` sembra banale -- verifica
un'equazione lineare -- ma applicato a `main2.m` e' il **canarino della
riparametrizzazione**: per far tornare `m(t_b) = m0 - Qdot*t_b` devono essere
simultaneamente corretti l'adimensionalizzazione della portata, le tre durate
d'arco, la scalatura chain-rule `Delta_k` applicata **anche** alla componente di
massa, la mappa tau -> t dentro l'EOM, la mappa inversa in post-processing e la
ri-dimensionalizzazione finale. Un errore in una sola di queste sei cose lo fa
fallire. A pari merito c'e' `testDimensionalVsNonDimensionalAgreement`, che e'
l'unico oracolo **indipendente** disponibile in assenza di soluzione analitica: due
implementazioni con scalature e parametrizzazioni diverse che convergono agli
stessi numeri.
