# HM1/tests/odeBurnPerformanceTest.m

## Ruolo del file nel progetto

E' il **benchmark** di `HM1/ode_burn.m`. Non verifica la correttezza (quello lo
fa `odeBurnTest.m`): **misura il tempo**. Esiste per una ragione precisa e molto
concreta, dichiarata nel commento di intestazione: `ode_burn` e' l'*hot loop* di
HM1, "called ~1e6+ times by the fsolve sweeps of main_task1..4" (righe 3-4).

Il file e' l'evidenza sperimentale a supporto di una **decisione di design**
controintuitiva presa in `ode_burn.m` e codificata nel `CLAUDE.md` del repo: la
funzione **non ha un blocco `arguments` di validazione**. Il commento in
`ode_burn.m` (righe 13-14) lo dice esplicitamente:

    % No arguments block by design: called ~1e6+ times inside ode45/fsolve;
    % validate at the call site.

La convenzione del repo e' che la validazione stia nei **boundary helper**
(le funzioni chiamate una volta per run: `nondim`, `solve_*`, plotting), non
nelle funzioni chiamate un milione di volte. Questo benchmark serve a poter
**difendere quella scelta con dei numeri** invece che con un'opinione.

La classe eredita da `matlab.perftest.TestCase` -- che a sua volta eredita da
`matlab.unittest.TestCase`. Conseguenza pratica: i suoi metodi sono **sia**
benchmark (con `runperf`) **sia** test veri (con `runtests`), perche' contengono
delle `verify*`. Si lancia con:

    results = runperf('odeBurnPerformanceTest')

come indicato alla riga 6.

Misura **due unita' diverse**:
1. `testRhsEvaluation` -- il costo di **una singola valutazione** della RHS: la
   grandezza microscopica, quella che viene moltiplicata per ~1e6.
2. `testBurnArcIntegration` -- il costo di **un'intera integrazione** dell'arco
   propulso con `ode45`, parametrizzato su tolleranza lasca vs stretta: la
   grandezza macroscopica, quella che l'utente percepisce come "quanto ci mette
   `main_task1`".

---

## Perche' `ode_burn` e' un hot loop -- il conto degli ordini di grandezza

Vale la pena ricostruire la stima `1e6+` invece di ripeterla a memoria, perche'
e' esattamente la domanda che ti fanno all'orale. (Sono stime di ordine di
grandezza, non misure: il codice non le verifica.)

La catena di chiamata di `main_task1.m` e':

    sweep su Q e yf (continuation)
      -> fsolve  (per ogni punto dello sweep)
           -> residual  (per ogni iterazione + per ogni colonna del Jacobiano)
                -> ode45  (una integrazione completa dell'arco)
                     -> ode_burn  (per ogni stage di ogni passo)

Moltiplicando:

- **`ode45` -> `ode_burn`**: Dormand-Prince e' un Runge-Kutta a 7 stadi con
  proprieta' FSAL, quindi ~6 valutazioni della RHS per passo **accettato**, piu'
  quelle dei passi rifiutati. A `RelTol = 1e-10` (la tolleranza del repo) i passi
  sono molti: ordine delle centinaia/migliaia per arco. Diciamo ~1e3-1e4
  valutazioni per integrazione.
- **`residual` -> `ode45`**: una integrazione completa per ogni valutazione del
  residuo.
- **`fsolve` -> `residual`**: `fsolve` non ha il Jacobiano analitico, quindi lo
  approssima per **differenze finite**: con 4 incognite servono `4 + 1 = 5`
  valutazioni del residuo per ogni Jacobiano, piu' quelle della line search /
  trust region. Con ~10-30 iterazioni per solve, siamo sull'ordine di 1e2
  integrazioni per solve.
- **sweep -> `fsolve`**: `main_task1` fa continuation su `Q` (in avanti e
  all'indietro dal punto di partenza) per **tre** valori di `yf`: la griglia e'
  `linspace(1.8, 7, 80)`, quindi fino a 80 x 3 = 240 solve -- centinaia, non decine.

Il prodotto -- 1e3 valutazioni/integrazione x 1e2 integrazioni/solve x centinaia di
solve -- arriva comodamente a **1e7 chiamate di `ode_burn` e oltre** in un singolo run
di `main_task1` (la stima `~1e6+` del commento e' quindi conservativa). Ecco perche' il
**costo unitario** di quella funzione conta davvero, ed ecco perche' esiste questo file.

> **Possibile domanda d'esame** -- Perche' non fornire il Jacobiano analitico a
> `fsolve` invece di ottimizzare la RHS?
> *Risposta:* Sarebbe la mossa piu' efficace in assoluto: eliminerebbe il fattore
> `n+1 = 5` delle differenze finite e migliorerebbe anche la robustezza della
> convergenza. Ma il Jacobiano del residuo di shooting richiede di propagare le
> **equazioni variazionali** (la matrice di transizione di stato) insieme alla
> traiettoria -- cioe' integrare 6 + 6x4 = 30 stati invece di 6. HM1 non lo fa: e'
> una semplificazione consapevole, che scarica il costo sul numero di
> integrazioni. Il benchmark documenta il prezzo pagato.

---

## Properties e fixture (righe 8-31)

```matlab
properties
    p    % costate/thrust parameter struct
    z0   % initial state [x; y; vx; vy; m; lam_m]
end

properties (TestParameter)
    RelTol = struct('Loose', 1e-6, 'Tight', 1e-10)
end
```

- **Righe 8-11 -- properties normali.** `p` e `z0` sono i dati del benchmark,
  riempiti dal `TestMethodSetup`. Sono property e non variabili locali perche'
  devono essere costruite **fuori** dalla regione misurata (vedi sotto): il
  `setup` non entra nel tempo.

- **Righe 13-15 -- `TestParameter`.** Questo e' il meccanismo di
  **parametrizzazione** di `matlab.unittest`: una struct i cui campi diventano
  altrettanti casi di test. Qui `RelTol` ha due valori -- `Loose = 1e-6` e
  `Tight = 1e-10` -- quindi il metodo `testBurnArcIntegration(testCase, RelTol)`
  viene eseguito **due volte**, e i risultati compaiono come due voci separate
  (`testBurnArcIntegration(RelTol=Loose)` e `(RelTol=Tight)`). I nomi dei campi
  (`Loose`, `Tight`) diventano le etichette nel report: e' il motivo per cui si usa
  una struct e non un cell array.

  Il **valore diagnostico** e' proprio nel rapporto fra i due tempi: quantifica
  quanto costa la tolleranza stretta che la convenzione del repo impone
  (`RelTol = 1e-10, AbsTol = 1e-12` per il lavoro indiretto/shooting). Non e' una
  tolleranza gratuita, e questo benchmark ne mette il prezzo su un numero.

- **Righe 17-22 -- `TestClassSetup`.** Identico a quello di `odeBurnTest.m`: risale
  di due livelli da `mfilename('fullpath')` per trovare `HM1/` e lo aggiunge al
  path con una `PathFixture` (rimozione automatica alla fine). Gira **una volta
  sola** per la classe.

- **Righe 24-31 -- `TestMethodSetup`.** Gira **prima di ogni metodo di test**.

  ```matlab
  % Representative Task 1 solution neighbourhood (Q = 2.5, yf = 0.04)
  testCase.p  = struct('T', 1.5, 'Q', 2.5, 'c', 0.6, ...
                       'lam_vx0', 0.6, 'lam_vy0', 3.8, 'lam_y', 14);
  testCase.z0 = [0; 0; 0; 0; 1; 1];
  ```

  I numeri **non sono arbitrari**: sono presi dall'intorno della soluzione vera di
  Task 1. Il README di HM1 riporta `Q* = 2.52` come portata ottima per `yf = 0.04`,
  e il commento alla riga 26 lo dichiara. La coerenza interna torna:
  `T/c = 1.5/0.6 = 2.5 = Q`, come deve essere visto che `dm/dt = -T/c = -Q`.

  **Perche' e' importante**: un benchmark ha senso solo se misura il **regime
  operativo reale**. Con costati piccoli e `lam_y` nullo, l'integratore farebbe
  pochi passi e il numero misurato sarebbe ottimistico e inutile. Con `lam_y = 14`
  l'angolo di spinta ruota rapidamente e `ode45` deve lavorare come lavora davvero
  dentro `fsolve`. `z0 = [0;0;0;0;1;1]` e' letteralmente la condizione iniziale
  usata da `main_task1` (`ic = [0;0;0;0;1;1]`, con `lam_m0 = 1` per la
  normalizzazione).

---

## `testRhsEvaluation` (righe 34-41)

```matlab
function testRhsEvaluation(testCase)
    pp = testCase.p;
    z  = [0.05; 0.01; 0.4; 0.2; 0.7; 1.2];
    while testCase.keepMeasuring
        dz = ode_burn(0.15, z, pp);
    end
    testCase.verifySize(dz, [6 1]);
end
```

- **Riga 35 -- `pp = testCase.p`.** Sembra pedanteria, non lo e'. L'accesso a una
  **property** di un oggetto MATLAB passa dal meccanismo delle classi ed e'
  sensibilmente piu' lento dell'accesso a una variabile locale. Copiando `p` in
  `pp` **fuori** dal ciclo, la misura riguarda `ode_burn` e non l'overhead di
  `testCase.p`. E' la regola generale del micro-benchmarking: tutto cio' che non
  vuoi misurare va fuori dalla regione misurata.

- **Righe 37-39 -- `while testCase.keepMeasuring`.** E' il costrutto centrale di
  `matlab.perftest`. Delimita la **regione cronometrata**: tutto cio' che sta
  prima (`pp = ...`, `z = ...`) e dopo (`verifySize`) e' **escluso** dal tempo. Il
  framework si occupa di ripetere il corpo quanto serve per superare la
  granularita' del timer (una singola chiamata a `ode_burn` dura frazioni di
  microsecondo, ben sotto la risoluzione di `tic/toc`) e di riportare il tempo
  **per iterazione**.

  Senza `keepMeasuring` (cioe' misurando l'intero metodo) il tempo sarebbe
  dominato dall'overhead del framework e la misura sarebbe priva di significato.

- **Riga 40 -- `verifySize` FUORI dal ciclo.** Due motivi: (a) l'asserzione non
  deve entrare nel tempo misurato; (b) mettendola comunque, il metodo resta un
  **test valido** eseguibile con `runtests`, non solo un benchmark. Questo e' il
  vantaggio di `matlab.perftest.TestCase < matlab.unittest.TestCase`.

  **Nota di stile**: `dz` viene usata dopo il ciclo, quindi il codice assume che
  il corpo del `while` giri **almeno una volta**. In pratica `keepMeasuring`
  garantisce sempre almeno un'iterazione, quindi funziona -- ma e' una dipendenza
  implicita dal comportamento del framework, non un invariante scritto.

- **Cosa misura davvero.** Il costo di: 7 accessi a campi della struct
  (`lam_vx0`, `lam_vy0`, `lam_y`, `Q`, piu' `T` letto tre volte) + 1 `sqrt` + 1
  `atan2` + 1 `sin` + 1 `cos` + 1 `zeros(6,1)` + poche moltiplicazioni. Le
  funzioni trascendenti (`atan2`, `sin`, `cos`) sono la parte dominante. E' il
  numero da moltiplicare per ~1e6 per capire quanto pesa `ode_burn` nel tempo
  totale di `main_task1`.

> **Possibile domanda d'esame** -- Come useresti questo benchmark per giustificare
> l'assenza del blocco `arguments` in `ode_burn.m`?
> *Risposta:* Misurerei `testRhsEvaluation` cosi' com'e' (baseline), poi
> aggiungerei a `ode_burn` un blocco `arguments` con i validatori tipici
> (`mustBeNumeric`, `mustBeReal`, controllo di dimensione su `z`, `mustBeField` su
> `p`) e rimisurerei. Il corpo della funzione costa poche centinaia di nanosecondi
> -- dominato da `atan2`, `sin`, `cos` -- mentre un blocco `arguments` con
> validatori costa tipicamente qualche microsecondo per chiamata. Il rapporto
> puo' quindi essere di un ordine di grandezza, e moltiplicato per 1e6+ chiamate
> diventa la differenza fra uno sweep che gira in decine di secondi e uno che gira
> in minuti. **Onesta': questo A/B il file non lo fa** -- fornisce solo il termine
> di paragone. La misura del ramo "con `arguments`" andrebbe fatta apposta.

---

## `testBurnArcIntegration` (righe 43-52)

```matlab
function testBurnArcIntegration(testCase, RelTol)
    pp   = testCase.p;
    ic   = testCase.z0;
    opts = odeset('RelTol', RelTol, 'AbsTol', RelTol*1e-2);
    while testCase.keepMeasuring
        [~, Z] = ode45(@(t,z) ode_burn(t, z, pp), [0 0.3], ic, opts);
    end
    % Mass equation is linear: exact propellant bookkeeping
    testCase.verifyEqual(Z(end,5), 1 - pp.Q*0.3, 'AbsTol', 1e-8);
end
```

- **Riga 43 -- la firma parametrizzata.** Il secondo argomento `RelTol` viene
  iniettato dal framework, che chiama il metodo una volta per ogni campo della
  `TestParameter`. Non e' un argomento che si passa a mano.

- **Riga 46 -- le tolleranze.** `AbsTol = RelTol * 1e-2`, quindi:

      Loose:  RelTol = 1e-6,   AbsTol = 1e-8
      Tight:  RelTol = 1e-10,  AbsTol = 1e-12

  Il caso `Tight` riproduce **esattamente** le tolleranze usate dai `main_task*.m`
  (`RelTol = 1e-10, AbsTol = 1e-12`, la convenzione del repo per il lavoro
  indiretto/shooting). Il caso `Loose` e' il termine di paragone.

  La regola `AbsTol = RelTol/100` e' un'euristica sensata: `AbsTol` deve dominare
  vicino allo zero (dove `RelTol` perde senso perche' moltiplica una quantita'
  nulla), e le componenti dello stato di HM1 sono O(0.1-1) in nondimensionale,
  quindi due ordini sotto `RelTol` e' un buffer ragionevole.

- **Righe 47-49 -- la regione misurata.** Una **intera integrazione** dell'arco
  propulso su `[0, 0.3]`. E' l'unita' di lavoro che `fsolve` chiama a ogni
  valutazione del residuo, quindi il numero misurato e' quello che, moltiplicato
  per il numero di valutazioni del residuo, da' il tempo di un solve.

  Anche qui `pp` e `ic` sono copiati fuori dal ciclo, e `opts` e' costruita fuori
  (`odeset` non e' gratis).

- **Riga 51 -- l'asserzione di sanita'.** `m(tf) = 1 - Q * 0.3 = 1 - 0.75 = 0.25`,
  esatto. Il commento (riga 50) spiega perche' funziona: **l'equazione di massa e'
  lineare** (`dm/dt = -Q` costante), quindi ogni metodo Runge-Kutta la integra in
  modo **esatto** a meno dell'arrotondamento macchina -- indipendentemente da
  `RelTol`. Ecco perche' la stessa `AbsTol = 1e-8` sull'asserzione passa sia nel
  caso `Loose` che nel caso `Tight`.

  E' un'asserzione ben scelta: e' l'**unica** quantita' del problema che si puo'
  asserire con la stessa tolleranza in entrambi i rami parametrizzati. Asserire
  `vy(tf)` con `1e-8` fallirebbe nel caso `Loose`. E' anche il motivo per cui il
  commento la chiama "exact propellant bookkeeping": non e' un test di
  accuratezza, e' un test che l'integrazione sia **avvenuta** e sull'arco giusto.

- **Il vero output del test non e' l'assert, e' il tempo.** Confrontando i due
  campioni si ottiene il **fattore di costo della tolleranza stretta**. E' un
  numero che vale la pena conoscere prima di dire "tanto le tolleranze si mettono
  strette e basta": in un metodo indiretto `RelTol = 1e-10` non e' un vezzo -- i
  costati sono estremamente sensibili alle condizioni iniziali e un'integrazione
  imprecisa produce un residuo di shooting rumoroso, che manda in confusione il
  Jacobiano alle differenze finite di `fsolve` e uccide la convergenza. Il
  benchmark dice **quanto** si paga quella necessita'.

> **Possibile domanda d'esame** -- Perche' un metodo indiretto ha bisogno di
> `RelTol = 1e-10` mentre un metodo diretto (HM2) sopravvive con tolleranze molto
> piu' larghe?
> *Risposta:* Perche' nello shooting il residuo e' una funzione dell'**intera
> propagazione**: l'errore dell'integratore entra direttamente nel residuo, e
> `fsolve` costruisce il Jacobiano per **differenze finite** perturbando le
> incognite di quantita' piccole. Se il rumore numerico dell'integratore e'
> dello stesso ordine della perturbazione, il Jacobiano stimato e' spazzatura e
> Newton non converge. Inoltre le dinamiche dei costati sono spesso instabili in
> avanti: piccoli errori si amplificano lungo l'arco. Nel collocation invece la
> traiettoria e' una variabile di decisione discretizzata e i vincoli sono
> **locali** (un residuo per intervallo): non c'e' propagazione da amplificare, e
> l'accuratezza si controlla infittendo la griglia, non stringendo un integratore.

---

## Limiti del file (onesta' richiesta)

- **Non e' un regression gate.** Non c'e' nessuna asserzione sul **tempo** (non
  c'e' nessun `verifyLessThan(t, soglia)`). `runperf` produce dei numeri, ma
  nessuno li confronta automaticamente con un baseline salvato: se domani una
  modifica rendesse `ode_burn` 10 volte piu' lenta, il test **passerebbe
  comunque**. Per farne un gate servirebbe salvare un baseline e confrontarlo (per
  esempio con `matlab.perftest.TimeExperiment` e una soglia esplicita).

- **Non misura l'alternativa.** Come detto sopra, non c'e' un ramo "con blocco
  `arguments`" da confrontare: il file fornisce la baseline che *renderebbe
  possibile* quel confronto, ma il confronto non e' fatto. La tesi
  "niente `arguments` perche' hot loop" resta quindi **argomentata ma non
  misurata** dentro il repo.

- **Non misura l'intero `main_task1`.** L'unita' massima e' una singola
  integrazione, non un solve di `fsolve` ne' uno sweep di continuation. La stima
  `1e6+` chiamate resta una stima analitica, non un dato raccolto (per contarle
  davvero basterebbe un contatore `persistent` in `ode_burn`, ma sarebbe intrusivo
  proprio nell'hot loop).

- **`c` e' nel parameter struct ma `ode_burn` non lo usa.** `ode_burn.m` legge
  solo `T`, `Q`, `lam_vx0`, `lam_vy0`, `lam_y`: il campo `c` viene passato ma mai
  letto (serve altrove, per esempio nel calcolo di `H(0)` dentro il residuo). Non
  e' un bug, ma e' un pelo di grasso nella struct che viene copiata a ogni
  chiamata.

---

## Possibili domande d'esame

**D: Che cos'e' `matlab.perftest` e cosa fa in piu' rispetto a `tic/toc`?**
R: E' il framework di performance testing di MATLAB, costruito sopra
`matlab.unittest`. Rispetto a un `tic/toc` a mano aggiunge tre cose. (1)
**Statistica**: `runperf` non fa una misura sola, ne fa molte -- di default esegue
alcune run di warm-up (per assorbire l'effetto della JIT compilation e delle
cache fredde) e poi raccoglie campioni finche' il margine di errore relativo non
scende sotto una soglia (default 5% al 95% di confidenza), fino a un massimo di
campioni. (2) **Isolamento della regione misurata**: `while testCase.keepMeasuring`
esclude setup, teardown e asserzioni dal cronometro, e ripete il corpo quanto
basta a superare la granularita' del timer. (3) **Integrazione con il framework di
test**: parametrizzazione (`TestParameter`), fixture, e il fatto che gli stessi
metodi restano test veri eseguibili con `runtests`.

**D: Perche' `ode_burn.m` non ha un blocco `arguments`, se il resto del repo lo
raccomanda?**
R: Perche' e' l'**hot loop**. `ode_burn` viene chiamata circa 1e6+ volte in un run
di `main_task1` (6 valutazioni per passo di `ode45`, x centinaia/migliaia di passi
per integrazione a `RelTol = 1e-10`, x ~5 integrazioni per Jacobiano alle
differenze finite su 4 incognite, x decine di iterazioni di `fsolve`, x centinaia di
solve nello sweep di continuation -- 80 nodi in `Q` x 3 valori di `yf`). Il corpo
della funzione costa poche centinaia
di nanosecondi; un blocco `arguments` con validatori costa tipicamente qualche
microsecondo -- potenzialmente **piu' della funzione stessa**. La convenzione del
repo (scritta in `CLAUDE.md`) e' quindi: la validazione vive nei **boundary
helper** chiamati una volta per run (`nondim`, `solve_*`, plotting), non nelle RHS
di ODE, nei residui di shooting o nelle callback di `nonlcon`. E' un
**trade-off consapevole**: si accetta un messaggio di errore piu' criptico in
cambio di un ordine di grandezza sul tempo di esecuzione.

**D: Perche' parametrizzare su `RelTol` invece di misurare solo la tolleranza
vera (`1e-10`)?**
R: Perche' un numero isolato non dice niente. Sapere che un'integrazione costa
`X` millisecondi e' inutile senza un termine di paragone; sapere che a `1e-10`
costa `k` volte quello che costa a `1e-6` **quantifica il prezzo della scelta di
tolleranza**, che e' una decisione di progetto documentata nel `CLAUDE.md`. Il
benchmark trasforma una convenzione ("tolleranze strette per il lavoro indiretto")
in un costo misurabile, e permette di rispondere a "e se le allentassi?" con un
numero invece che con un'opinione.

**D: L'asserzione sulla massa (`m(tf) = 1 - Q*tf`) passa con la stessa tolleranza
sia a `RelTol = 1e-6` che a `1e-10`. Come e' possibile?**
R: Perche' `dm/dt = -Q` e' **costante**: la soluzione e' un polinomio di grado 1 in
`t`, e ogni metodo Runge-Kutta di ordine >= 1 integra esattamente i polinomi di
grado 1, a meno dell'arrotondamento macchina. Il controllo di passo adattivo non
c'entra: l'errore locale su quella componente e' zero per costruzione. Ecco perche'
e' l'asserzione giusta da mettere in un test **parametrizzato sulla tolleranza** --
e' l'unica quantita' che si puo' asserire con la stessa `AbsTol` in entrambi i
rami. Se avessi asserito `vy(tf)` con `1e-8`, il ramo `Loose` sarebbe fallito.

**D: Il benchmark protegge da regressioni di performance?**
R: **No, non da solo**, e va detto. Non c'e' nessuna asserzione sul tempo: `runperf`
stampa i tempi ma non li confronta con un baseline salvato. Se `ode_burn` diventasse
10 volte piu' lenta, la suite passerebbe. Per farne un vero gate servirebbe salvare i
risultati di `runperf` come baseline (per esempio serializzando i
`MeasurementResult`) e aggiungere un confronto esplicito, oppure inserire una soglia
assoluta. Oggi il file e' uno **strumento di misura**, non un guardrail.
