# HM2_powered_descent/tests/descentDynamicsPerformanceTest.m

## Ruolo del file nel progetto

Questo file non verifica *correttezza*: verifica *velocita'*. E' l'unico
benchmark della repo HM2 e misura le tre funzioni che stanno sul cammino
critico dell'NLP: `ode_descent` (la RHS continua del punto materiale),
`rk4_zoh` (il propagatore a passo fisso usato dai defect ZOH) e una
replica `ode45` di un singolo intervallo ZOH alle tolleranze usate dai
forward-integrate di verifica.

Il motivo per cui esiste e' una scelta di design che attraversa tutta la
repo. In HM2 il problema di powered descent e' trascritto come NLP e
risolto con `fmincon` (SQP). Il vincolo non lineare `trap_nonlcon`
(`main_task2.m` righe 759-781) chiama `ode_descent` una volta per nodo, e
`zoh_nonlcon` (righe 783-803) chiama `rk4_zoh` una volta per intervallo.
`fmincon` non riceve gradienti analitici (`fmincon_opts`, righe 731-735,
non imposta `SpecifyConstraintGradient`), quindi ricostruisce lo Jacobiano
dei vincoli per differenze finite: con `N = 50` nodi il vettore delle
variabili ha 7*N = 350 componenti, e ogni Jacobiano costa ~351 valutazioni
di `nonlcon`. Moltiplicando, `ode_descent` viene chiamata dell'ordine di
1e7 volte in una singola soluzione. E' un hot loop nel senso stretto del
termine.

Da qui la regola dichiarata in `CLAUDE.md` e ripetuta come commento nei
sorgenti: le RHS e i propagatori **non hanno blocchi `arguments`** di
validazione (`ode_descent.m` righe 10-11, `rk4_zoh.m` righe 12-13), perche'
il costo della validazione verrebbe pagato milioni di volte contro un corpo
funzione che fa poche operazioni in virgola mobile. La validazione vive nei
*boundary helper* chiamati una volta per run (`nondim`, `solve_*`,
`fwd_integrate_*`, `lti_zoh`). Questo benchmark e' il presidio di quella
scelta: se qualcuno "ripulisce" `ode_descent` aggiungendo validazione, o
sostituisce l'indicizzazione esplicita con qualcosa di piu' elegante ma piu'
lento, `runperf` deve mostrare il degrado.

Onesta': il file **non** misura `ode_descent_uacc.m` ne' `lti_zoh.m`, che
pure fanno parte della famiglia. Vedremo piu' avanti perche' l'esclusione e'
difendibile, ma va detto subito che la copertura di performance e' parziale.

---

## Il framework: `matlab.unittest.TestCase` vs `matlab.perftest.TestCase`

Prima del codice, il pezzo di teoria che l'orale puo' chiedere.

**`matlab.unittest.TestCase`** e' la classe base dei test funzionali: ogni
metodo `Test` esegue del codice e fa asserzioni (`verifyEqual`,
`verifyLessThan`, ...). Il framework riporta pass/fail. Il tempo di
esecuzione viene misurato, ma solo come informazione accessoria: e' il tempo
totale del metodo, setup incluso, misurato una volta sola.

**`matlab.perftest.TestCase`** (riga 1) eredita da `matlab.unittest.TestCase`
e aggiunge la nozione di **confine di misura** (*measurement boundary*).
Il punto e' che il tempo di un metodo di test *nella sua interezza* e'
inutile come metrica: comprende la costruzione dell'oggetto test, il
`TestMethodSetup`, l'allocazione delle variabili, le asserzioni finali. Se
il codice che ti interessa dura 60 nanosecondi e il setup ne dura 200
microsecondi, stai misurando il setup. Serve un modo per dire al framework
"cronometra **solo** questa regione".

Ci sono due API per farlo:

1. **`startMeasuring` / `stopMeasuring`.** Racchiudono esplicitamente la
   regione di interesse:

       function testSomething(testCase)
           payload = buildBigThing();   % NON misurato
           testCase.startMeasuring();
           result = functionUnderTest(payload);
           testCase.stopMeasuring();
           testCase.verifyEqual(result, expected);  % NON misurato
       end

   Il framework esegue il metodo per intero, ma registra solo il tempo fra
   `start` e `stop`. Accettano un `label` opzionale per definire piu' confini
   distinti nello stesso metodo.

2. **`keepMeasuring`** (introdotta in R2018b), che e' quella usata da questo
   file. Si usa come condizione di un `while`:

       while testCase.keepMeasuring
           <codice da misurare>
       end

   La differenza non e' cosmetica. `startMeasuring`/`stopMeasuring`
   cronometrano **una** esecuzione della regione: se quella regione dura meno
   della risoluzione del timer, il numero e' rumore puro. `keepMeasuring`
   invece dice al framework: *itera questo ciclo quante volte ti servono per
   ottenere una misura accurata*. Il framework decide autonomamente il numero
   di ripetizioni del `while` e normalizza il risultato. E' l'API giusta per
   codice velocissimo, ed e' esattamente il caso di `ode_descent`.

### Perche' misurare solo la regione di interesse e' l'unico modo

Detto in modo operativo: qualsiasi tempo misurato e' `t_setup + t_utile`. Se
vuoi rilevare una regressione del 20% su `t_utile`, ma `t_setup >> t_utile`,
il segnale che cerchi e' sepolto sotto una costante additiva enorme. Il
rapporto segnale/rumore di un benchmark e'

    SNR ~ t_utile / (t_setup + jitter_del_sistema_operativo)

e l'unica leva che hai sul numeratore e' escludere il setup dal confine di
misura. In questo file l'esclusione avviene su due livelli: il
`TestMethodSetup` (righe 29-35) sta fuori dal `while` per costruzione del
framework, e le variabili locali (`x`, `u`, `vc`, `h`, `opts`) sono
preparate prima del `while` dentro ciascun metodo.

### Come `runperf` fa statistica (e perche' un `timeit` non basta)

`runperf('descentDynamicsPerformanceTest')` costruisce un *time experiment*
di tipo `limitingSamplingError`. I default documentati sono:

| Parametro                | Default | Significato                                  |
| ------------------------ | ------- | -------------------------------------------- |
| `NumWarmups`             | 5       | esecuzioni di riscaldamento, scartate         |
| `MinSamples`             | 4       | campioni minimi raccolti                      |
| `MaxSamples`             | 256     | tetto ai campioni                             |
| `RelativeMarginOfError`  | 0.05    | obiettivo: 5% di margine di errore sulla media |
| `ConfidenceLevel`        | 0.95    | livello di confidenza                         |

La logica e':

- **Warm-up (5 esecuzioni scartate).** La prima esecuzione di un file `.m` in
  MATLAB paga il parsing e la compilazione JIT, e i dati non sono ancora in
  cache. Includere quelle esecuzioni nella statistica falserebbe la media
  verso l'alto. Le si esegue e le si butta.
- **Campionamento adattivo.** Il framework continua a raccogliere campioni
  finche' la media campionaria non e' nota entro il 5% di margine relativo al
  95% di confidenza, oppure finche' non arriva a 256 campioni. Se la macchina
  e' rumorosa (altri processi, throttling termico), servono piu' campioni; se
  e' quieta, il framework si ferma prima. E' un criterio di **arresto
  statistico**, non un conteggio fisso.
- **Riepilogo.** `sampleSummary(results)` produce una tabella con media,
  mediana, deviazione standard, min e max dei campioni. La **mediana** e' la
  statistica che conviene leggere per confrontare due run: e' robusta agli
  outlier, e gli outlier in un benchmark su un sistema operativo general
  purpose sono la norma (un context switch, il garbage collector di un altro
  processo, un interrupt).

Perche' un singolo `timeit` non basta: `timeit` restituisce **un** numero
(esso stesso una mediana di ripetizioni interne, per inciso), senza margine
di errore associato e senza criterio di arresto. Non ti dice se la differenza
fra 61 ns e 66 ns e' una regressione o rumore. Un time experiment ti da' una
distribuzione con un intervallo di confidenza dichiarato: e' l'unica base su
cui puoi affermare che due misure differiscono davvero. Inoltre `timeit` non
si integra con il resto della suite `matlab.unittest` (fixture, parametri,
selezione dei test), che qui serve.

---

## `classdef` e proprieta' (righe 1-20)

```matlab
classdef descentDynamicsPerformanceTest < ...
        matlab.perftest.TestCase
    properties (Constant)
        Vc = 0.0777    % V_ref/c (Table 1 data)
        dt = 0.0444    % one ZOH interval, non-dim
    end
    properties
        x0
        u0
    end
    properties (TestParameter)
        nSub = struct('one', 1, 'two', 2, 'eight', 8)
    end
```

- Riga 1: eredita da `matlab.perftest.TestCase`, non da
  `matlab.unittest.TestCase`. E' questa riga a rendere disponibili
  `keepMeasuring` / `startMeasuring` e a far si' che `runperf` accetti la
  classe. Nota che la classe resta comunque un `matlab.unittest.TestCase`
  per ereditarieta', quindi `runtests` la esegue lo stesso -- eseguirebbe i
  metodi una volta sola, trattando i `verifySize` come normali asserzioni.
- Righe 8-11: le due costanti fisiche del benchmark, **hard-coded**.
  Verifichiamole contro la vera non-dimensionalizzazione
  (`main_task2.m`, `nondim`, righe 203-228):

      L_ref = y0 = 3000 m
      g_ref = 9.81 m/s^2
      t_ref = sqrt(L_ref/g_ref) = sqrt(3000/9.81) = 17.4874 s
      V_ref = sqrt(g_ref*L_ref) = sqrt(9.81*3000) = 171.55 m/s
      c     = Isp*g0 = 225*9.80665 = 2206.50 m/s
      Vc    = V_ref/c = 171.55/2206.50 = 0.077749

  che arrotonda a 0.0777: il valore torna. E per `dt` (con `tf = 38 s`,
  `main_task2.m` riga 35, e `N = 50`, riga 36):

      tf_nd = 38/17.4874 = 2.17299
      dt    = tf_nd/(N-1) = 2.17299/49 = 0.044347

  Attenzione: 0.044347 arrotondato a quattro decimali fa **0.0443**, non
  0.0444. La costante hard-coded e' il risultato di un doppio arrotondamento
  (0.044347 -> 0.04435 -> 0.0444) ed e' alta di circa lo 0.1%. Per un
  benchmark e' innocuo -- `dt` fissa solo il punto di lavoro, non un risultato
  fisico -- ma va detto: il valore *non* e' l'arrotondamento corretto di
  `tf_nd/(N-1)`. Il commento "N = 50" e' invece esatto.
  **Limite onesto:** sono duplicati. La sorgente di verita' e' `nondim`, che
  qui non viene chiamata. Se qualcuno cambiasse lo schema di scale di
  riferimento (per esempio `L_ref = x0` invece di `y0`), il benchmark
  continuerebbe a girare felicemente su costanti stantie, misurando un punto
  di lavoro che non esiste piu' nel problema. Non e' un errore -- e' un
  rischio di manutenzione, e va dichiarato.
- Righe 13-16: `x0` e `u0` sono proprieta' normali (non `Constant`) perche'
  vengono riscritte dal `TestMethodSetup` prima di ogni metodo.
- Righe 18-20: `nSub` e' una `TestParameter`. Il framework **moltiplica** il
  metodo di test che la dichiara come argomento, generando tre casi
  indipendenti (`one`, `two`, `eight`). Non e' zucchero sintattico: ogni caso
  ottiene la sua serie di campioni e la sua statistica, quindi si legge il
  *costo in funzione del numero di substep*. Vedi la sezione su
  `testRk4ZohPropagation` per cosa ci si compra con questo.

---

## `addHm2ToPath` (righe 22-27)

```matlab
methods (TestClassSetup)
    function addHm2ToPath(testCase)
        hm2 = fileparts(fileparts(mfilename('fullpath')));
        testCase.applyFixture( ...
            matlab.unittest.fixtures.PathFixture(hm2));
    end
end
```

- Riga 22: `TestClassSetup` gira **una volta per classe**, non per metodo.
  Giusto: mettere il path e' un'operazione idempotente e costosa, non ha
  senso ripeterla.
- Riga 24: `mfilename('fullpath')` da' il path di questo file
  (`.../HM2_powered_descent/tests/descentDynamicsPerformanceTest`); il primo
  `fileparts` sale a `.../HM2_powered_descent/tests`, il secondo a
  `.../HM2_powered_descent`. Da li' `ode_descent.m` e `rk4_zoh.m` sono
  visibili. Il doppio `fileparts` e' l'idioma standard per "la cartella
  padre della cartella dei test".
- Riga 25: `PathFixture` non e' un semplice `addpath`. E' una *fixture*: il
  framework la installa prima e la **disinstalla dopo**, ripristinando il
  path originale anche se un test fallisce o solleva un'eccezione. E' la
  ragione per cui non si usa `addpath` a mano -- quello lascerebbe la
  sessione MATLAB sporca.

---

## `setupState` (righe 29-35)

```matlab
methods (TestMethodSetup)
    function setupState(testCase)
        % Mid-descent state, near-hover control
        testCase.x0 = [0.1; 0.5; -0.2; -0.4; 0.8];
        testCase.u0 = [0.1; 0.9];
    end
end
```

- Righe 29-30: `TestMethodSetup` gira prima di **ogni** metodo di test. Sta
  fuori dal confine di misura per costruzione, quindi il suo costo non entra
  nel benchmark. Anche se qui e' banale (due assegnazioni), la struttura e'
  quella corretta: se domani lo stato iniziale richiedesse un calcolo, non si
  dovrebbe toccare nulla.
- Riga 32: lo stato e' `[x; y; vx; vy; m]` non-dim. Ridimensionalizzando con
  le scale sopra:

      x  = 0.1  * 3000   =  300 m
      y  = 0.5  * 3000   = 1500 m
      vx = -0.2 * 171.55 = -34.3 m/s
      vy = -0.4 * 171.55 = -68.6 m/s
      m  = 0.8  * 2000   = 1600 kg

  cioe' effettivamente un punto a meta' discesa, in caduta, con circa il 20%
  di propellente gia' bruciato. Il commento non mente.
- Riga 33: `u0 = [0.1; 0.9]` e' la spinta non-dim `[Tx; Ty]`, con
  `T_ref = m0*g = 2000*9.81 = 19620 N`. Verifichiamo il commento "near-hover":
  in unita' non-dim la gravita' vale 1 (vedi `ode_descent.m` riga 14, il `-1`
  nella riga di `vy`), quindi la spinta di hover e' `Ty = m = 0.8`. Qui
  `Ty = 0.9`, appena sopra hover, e `|T| = sqrt(0.1^2+0.9^2) = 0.905`, ben
  dentro il bound `Tmax = 70000/19620 = 3.57`. Il commento e' accurato.

**Perche' questo punto di lavoro e non un altro.** Non e' indifferente:
`ode_descent` contiene una divisione `u/x(5)` e una `sqrt`. Se si scegliesse
`u = [0; 0]`, `Tmag = 0` e la `sqrt(0)` potrebbe imboccare un percorso
denormale/degenere; se `m` fosse vicino a zero, la divisione esploderebbe.
Un punto "generico e sano" evita di misurare accidentalmente un caso
patologico della FPU. Detto questo, il codice **non** documenta questa
motivazione: e' una lettura, non una citazione.

> **Possibile domanda d'esame** -- Perche' il setup dello stato sta in un
> `TestMethodSetup` e non dentro il metodo di test, prima del `while`?
> *Risposta:* Funzionalmente sarebbe equivalente, perche' in entrambi i casi
> resta fuori dal confine di misura. La ragione e' di condivisione: i tre
> metodi di test usano lo **stesso** punto di lavoro, e centralizzarlo
> garantisce che confrontino la stessa cosa. Se un domani si misurasse un
> altro stato, si cambierebbe in un posto solo.

---

## `testOdeDescentEvaluation` (righe 38-48)

```matlab
function testOdeDescentEvaluation(testCase)
    % One call (~60 ns) is below framework precision;
    % measure a batch of 1000 per sample
    x = testCase.x0;  u = testCase.u0;  vc = testCase.Vc;
    while testCase.keepMeasuring
        for k = 1:1000
            dx = ode_descent(x, u, vc);
        end
    end
    testCase.verifySize(dx, [5 1]);
end
```

- Riga 41: le tre variabili sono estratte dalle proprieta' **prima** del
  `while`. Non e' pedanteria: l'accesso a una proprieta' di un oggetto
  (`testCase.x0`) in MATLAB passa dal dispatch delle property, che non e'
  gratuito. Se `testCase.x0` fosse dentro il ciclo, si misurerebbe in parte
  il costo di accesso alla property invece della RHS. Questo e' il tipo di
  contaminazione che rende inutile un microbenchmark.
- Righe 42-46: il confine di misura. Il corpo del `while` **non** e' una
  singola chiamata: e' un batch di 1000. La ragione e' nel commento (righe
  39-40): una chiamata sola sta sotto la risoluzione del framework. Il corpo
  di `ode_descent` (righe 13-14 del suo file) e' una `sqrt`, una radice di
  somma di quadrati, due divisioni e un assemblaggio di vettore 5x1: e'
  dell'ordine delle decine di nanosecondi, mentre l'overhead di una singola
  misura del framework e' ordini di grandezza sopra. Con il batch, il tempo
  misurato per campione e' ~1000 volte quello di una chiamata, e il rapporto
  segnale/overhead diventa favorevole.

  **Attenzione al numero riportato.** `keepMeasuring` normalizza per
  iterazione del `while`, non per iterazione del `for` interno. Il valore che
  `runperf` stampa per questo test e' quindi il tempo di **1000 chiamate**.
  Per avere il costo unitario bisogna dividere a mano per 1000. Questo e' un
  trabocchetto reale, e peggiora nel test successivo (batch da 100): i numeri
  dei due test **non sono confrontabili direttamente**.

  Nota sul `~60 ns` del commento: e' un'affermazione dell'autore, non un
  valore misurato o verificato da questo test. Il test non contiene nessuna
  soglia che lo controlli.
- Riga 44: `dx` viene sovrascritta a ogni iterazione. Gli argomenti `x`, `u`,
  `vc` sono **identici** per tutte e 1000 le chiamate. Questo garantisce che
  i dati siano in cache L1 e che la predizione dei branch sia perfetta: e'
  una misura di **best case**. Dentro `fmincon`, invece, ogni valutazione a
  differenze finite perturba una componente diversa del vettore, e i dati
  scorrono. Il benchmark misura il costo aritmetico della RHS, non il costo
  che la RHS ha realmente nel contesto di uso. Vedi la sezione sui limiti.
- Riga 47: `verifySize(dx, [5 1])` sta **dopo** il `while`, quindi fuori dal
  confine di misura. Non e' un test di correttezza (quella vive in
  `tests/odeDescentTest.m`): e' un guardrail contro il benchmark che misura
  una funzione che nel frattempo ha smesso di restituire un vettore 5x1, e
  contro l'ipotesi remota che l'engine elimini una chiamata il cui risultato
  non viene mai usato. Usare il risultato e' buona igiene da microbenchmark.

---

## `testRk4ZohPropagation` (righe 50-61)

```matlab
function testRk4ZohPropagation(testCase, nSub)
    % Batch of 100 per sample: the 1-substep case is
    % sub-microsecond and noise-dominated measured singly
    x  = testCase.x0;  u = testCase.u0;
    vc = testCase.Vc;  h = testCase.dt;
    while testCase.keepMeasuring
        for k = 1:100
            xn = rk4_zoh(x, u, h, vc, nSub);
        end
    end
    testCase.verifySize(xn, [5 1]);
end
```

- Riga 50: la firma prende `nSub` come secondo argomento. E' il meccanismo
  di parametrizzazione di `matlab.unittest`: il framework riconosce che
  `nSub` e' una `TestParameter` (righe 18-20) e istanzia tre test distinti
  con `nSub = 1, 2, 8`. Nei risultati compaiono come
  `testRk4ZohPropagation(nSub=one)`, `(nSub=two)`, `(nSub=eight)`.
- Righe 55-59: batch da 100 (non 1000). Un `rk4_zoh` con `n_sub = 1` fa 4
  chiamate a `ode_descent` (righe 17-20 di `rk4_zoh.m`), quindi costa circa
  4 volte una RHS piu' l'overhead di chiamata; con `n_sub = 2` (il valore di
  produzione) le valutazioni sono 8, con `n_sub = 8` sono 32. Resta sotto il
  microsecondo, come dice il commento, ma essendo da 4 a 32 volte piu' caro di
  una RHS singola basta un batch dieci volte piu' piccolo per superare la
  risoluzione del timer. Il costo di ogni caso e', in prima approssimazione:

      costo(n_sub) = a + b * n_sub

  dove `b` e' il costo di un substep (4 valutazioni di `ode_descent` piu' la
  combinazione lineare `x + (h/6)*(k1 + 2*k2 + 2*k3 + k4)`, riga 21) e `a` e'
  l'overhead fisso della chiamata di funzione MATLAB piu' la divisione
  `h = dt/n_sub` (riga 15). **E' questo che compra la parametrizzazione**:
  con tre punti (1, 2, 8) si puo' separare `a` da `b`. Se un giorno il tempo
  del caso `nSub=8` non fosse circa `a + 8b` coerente con gli altri due, si
  saprebbe che qualcosa e' cambiato *dentro* il ciclo di substep e non
  nell'overhead di chiamata. E' diagnostica, non solo misura.
- Riga 57: nota che `n_sub = 2` e' il valore effettivamente usato in
  produzione (`main_task2.m` riga 37, `n_sub = 2`). I casi 1 e 8 sono li' per
  la diagnostica di cui sopra, non perche' il solver li usi.

> **Possibile domanda d'esame** -- Perche' il batch di questo test e' 100 e
> quello di `ode_descent` e' 1000? E' un problema?
> *Risposta:* Perche' `rk4_zoh` costa `4*n_sub` valutazioni di RHS piu'
> l'overhead di chiamata -- da 4 volte una singola RHS (`n_sub = 1`) a 32
> volte (`n_sub = 8`) -- quindi bastano meno ripetizioni per superare la
> risoluzione del timer. E' pero' un problema di
> leggibilita': `runperf` riporta il tempo del corpo del `while`, cioe' del
> batch intero, quindi i numeri dei due test sono in "unita'" diverse (tempo
> di 1000 RHS contro tempo di 100 intervalli RK4). Per confrontarli bisogna
> normalizzare a mano. Un design piu' pulito userebbe la stessa dimensione di
> batch, o normalizzerebbe esplicitamente.

---

## `testOde45ZohReplayInterval` (righe 63-71)

```matlab
function testOde45ZohReplayInterval(testCase)
    x  = testCase.x0;  u = testCase.u0;
    vc = testCase.Vc;  h = testCase.dt;
    opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
    while testCase.keepMeasuring
        [~, Y] = ode45(@(t, xx) ode_descent(xx, u, vc), ...
                       [0 h], x, opts);
    end
    testCase.verifySize(Y(end,:), [1 5]);
end
```

Questo e' il test di **livello sistema** annunciato nel commento di classe
(riga 5). Non misura una funzione, misura un'operazione dell'utente finale:
integrare un intervallo ZOH con `ode45` alle tolleranze di verifica.

- Riga 66: `odeset` sta **fuori** dal `while`. Corretto: costruire una struct
  di opzioni e' relativamente costoso e non fa parte di cio' che si vuole
  misurare. Le tolleranze `RelTol = 1e-10`, `AbsTol = 1e-12` non sono
  arbitrarie: sono **esattamente** quelle di `fwd_integrate` (`main_task2.m`
  riga 839) e `fwd_integrate_uacc` (riga 1295), le funzioni che replicano la
  soluzione ottimizzata per misurarne la fedelta'. Il benchmark quindi
  riproduce il carico reale della fase di verifica, non un carico inventato.
  Sono anche coerenti con la convenzione della repo per il lavoro
  indiretto/shooting (`RelTol=1e-10, AbsTol=1e-12`).
- Riga 67: nessun batch. Un `ode45` a tolleranza 1e-10 su un intervallo fa
  molti passi adattivi, ognuno con 6 valutazioni della RHS (Dormand-Prince):
  e' abbastanza costoso da essere misurabile singolarmente. Quindi qui il
  numero riportato e' il tempo di **una** replica di intervallo -- di nuovo,
  un'unita' diversa dai due test precedenti.
- Riga 68: **divergenza dal codice di produzione.** Qui la RHS e'
  `@(t, xx) ode_descent(xx, u, vc)`, che cattura `u` direttamente. In
  `fwd_integrate` (righe 844-852) la RHS e' invece
  `@(tt, x) ode_descent(x, u_fcn(tt), d.Vc)`, dove `u_fcn` e' a sua volta una
  anonymous function. La produzione paga quindi **una dispatch di anonymous
  function in piu' a ogni valutazione della RHS**, che il benchmark non paga.
  Il benchmark e' cioe' sistematicamente ottimista rispetto al percorso reale
  `'zoh'` di `fwd_integrate`. La differenza e' probabilmente piccola in
  rapporto al costo di `ode45`, ma non e' zero e non e' dichiarata da nessuna
  parte nel codice.
- Riga 68 (secondo punto): la anonymous function viene **costruita dentro** il
  `while`, quindi il suo costo di costruzione entra nella misura. E'
  un'inconsistenza rispetto a `opts`, che e' stato accuratamente issato fuori.
  Piccola, ma e' una contaminazione del confine di misura.
- Riga 70: `Y(end,:)` e' `1x5` (`ode45` restituisce le righe come istanti di
  tempo), da cui `verifySize(..., [1 5])` e non `[5 1]`. Consuma il risultato
  per evitare che venga considerato morto.

> **Possibile domanda d'esame** -- Perche' misurare un `ode45` in un test di
> performance, se `ode45` non compare mai dentro l'NLP?
> *Risposta:* Perche' non e' sul cammino critico del solver, ma e' sul cammino
> critico della **verifica**: `fwd_integrate` replica la soluzione intervallo
> per intervallo (`main_task2.m` righe 840-855) per calcolare il node error,
> ed e' l'operazione che produce le metriche di fedelta' riportate nel README
> (errore 1.4e-8 non-dim per la variante ZOH, 7.3e-12 per GFOLD). Con N = 50
> nodi si fanno 49 di queste integrazioni per ogni soluzione da validare, e ce
> ne sono quattro varianti piu' uno sweep. Se il costo per intervallo
> raddoppiasse, raddoppierebbe il tempo di tutta la fase di post-processing.

---

## Il nesso con il design: perche' le RHS non hanno `arguments` -- e perche' `lti_zoh` invece ce l'ha

E' il punto piu' probabile all'orale, perche' sembra un'incoerenza. Non lo e'.

### Le funzioni hot loop: nessuna validazione

`ode_descent.m` (righe 10-11) e `rk4_zoh.m` (righe 12-13) portano lo stesso
commento:

    % No arguments validation by design: hot-loop RHS inside
    % ode45/fmincon; validate at the call site.

Lo stesso commento compare su `trap_nonlcon` (riga 769), `zoh_nonlcon`
(riga 794), `path_ineq` (riga 812) e sulla RHS aumentata dell'Appendix A
(riga 288). Il conteggio delle chiamate, ricavato dal codice:

**Percorso trapezoidale** (`solve_trap` -> `trap_nonlcon`, righe 759-781):

    variabili di decisione   = 7*N = 7*50 = 350
    Jacobiano a diff. finite = 350 + 1 = 351 valutazioni di nonlcon
    ode_descent per nonlcon  = N = 50          (righe 772-774)
    ------------------------------------------------------------
    ode_descent per Jacobiano ~ 351 * 50 = 1.8e4
    MaxIterations = 1000     (riga 733)
    ------------------------------------------------------------
    totale ~ 1.8e7 chiamate (piu' quelle della line search SQP)

**Percorso ZOH** (`solve_zoh` -> `zoh_nonlcon`, righe 783-803):

    rk4_zoh per nonlcon      = N-1 = 49        (righe 797-800)
    ode_descent per rk4_zoh  = 4*n_sub = 4*2 = 8
    ode_descent per nonlcon  = 49*8 = 392
    ------------------------------------------------------------
    ode_descent per Jacobiano ~ 351 * 392 = 1.4e5
    totale al cap di 1000 iter ~ 1.4e8 chiamate

Contro un corpo funzione che fa una `sqrt` e due divisioni, un blocco
`arguments` con validatori (`mustBeFinite`, controllo di size, ecc.) e' un
costo dello stesso ordine o maggiore del calcolo stesso, pagato 1e7-1e8
volte. **Attenzione all'onesta':** questa e' l'argomentazione di design, ma
la repo **non contiene** un benchmark A/B che confronti `ode_descent` con e
senza `arguments`. Questo file misura la funzione **com'e'**; non dimostra
che la scelta fosse corretta, la presidia soltanto contro regressioni future.
Se all'orale ti chiedessero "e di quanto sarebbe piu' lenta con la
validazione?", la risposta onesta e' "non l'ho misurato".

### `lti_zoh` ha un blocco `arguments` (righe 19-22): non e' un'incoerenza

```matlab
function [Abar, Bbar, cbar] = lti_zoh(dt, Vc)
    arguments
        dt (1,1) double {mustBePositive, mustBeFinite}
        Vc (1,1) double {mustBeFinite}
    end
    ...
    E = expm([A, B, c; zeros(4,9)] * dt);
```

La discriminante non e' "che tipo di funzione e'", ma **quante volte viene
chiamata**. E il codice risponde in modo netto: in `solve_gfold_scvx`,
`lti_zoh` e' invocata alla **riga 1218**, mentre il ciclo SCvx comincia alla
**riga 1239** (`for iter = 1:max_iter`). Cioe' e' **fuori dal ciclo**, una
volta sola per risoluzione, con tanto di commento in linea:
`% exact LTI ZOH, once`.

La ragione fisica per cui puo' stare fuori dal ciclo e' il cuore del metodo
GFOLD: il cambio di variabili `z = ln(m)`, `u = T/m` (Acikmese & Blackmore)
rende la dinamica **esattamente LTI** (vedi l'intestazione di `lti_zoh.m`,
righe 3-5). Se il sistema e' tempo-invariante, le matrici discrete
`Abar, Bbar, cbar` sono **costanti su tutta la griglia** e non dipendono
dalla traiettoria di riferimento. Quindi non c'e' nulla da ricalcolare a ogni
iterazione SCvx. Confronta con `compute_ltv_zoh` (riga 305) del percorso LTV,
che invece deve reintegrare l'ODE aumentata dell'Appendix A a **ogni**
iterazione, perche' li' la linearizzazione dipende dal riferimento corrente.

Conseguenza: `lti_zoh` viene chiamata O(1) volte per risoluzione, e il suo
costo e' comunque dominato da una `expm` di una matrice 9x9. La validazione
degli argomenti e' rumore in confronto, e in cambio protegge da un `dt <= 0`
o non finito, che produrrebbe una `expm` silenziosamente sbagliata. **E'
esattamente la regola "la validazione vive nei boundary helper chiamati una
volta per run"**, applicata correttamente. Nessuna incoerenza.

Coerentemente, anche `solve_gfold_scvx` (righe 1211-1217), `nondim`
(righe 210-212), `fmincon_opts` (righe 728-730) e `fwd_integrate_uacc`
(righe 1289-1292) hanno blocchi `arguments`: sono tutti boundary.

### E `ode_descent_uacc`?

Il suo commento (riga 16) dice `No arguments validation by design: hot-loop
RHS inside ode45`, e il codice conferma: e' chiamata solo da
`fwd_integrate_uacc` (`main_task2.m` riga 1298) e da
`proto_gfold_logmass.m` (riga 249), entrambe dentro un `ode45`. Quindi e'
hot **rispetto a `ode45`** (migliaia di valutazioni per intervallo), ma
**non** e' sul cammino critico dell'NLP: nel percorso GFOLD la dinamica
dentro il ciclo SCvx e' gestita da `Abar/Bbar/cbar`, non da valutazioni della
RHS. E' un "warm loop": chiamato molto, ma una volta sola per soluzione, in
fase di verifica. Che non sia benchmarkato in questo file e' difendibile su
questa base; resta il fatto che **non e' benchmarkato**, e il commento di
classe (righe 2-5) non spiega l'esclusione.

---

## Cosa significherebbe una regressione qui

Il valore di questo file sta tutto nel fattore di amplificazione. Usando i
conteggi ricavati sopra, una regressione di **100 ns per chiamata** di
`ode_descent` (cioe' circa un raddoppio, se il commento dei ~60 ns e' nel
giusto ordine di grandezza) si propaga cosi':

| Percorso                | Chiamate `ode_descent` | Costo aggiunto @100 ns |
| ----------------------- | ---------------------- | ---------------------- |
| `solve_trap` (Task 1/2) | ~1.8e7                 | ~1.8 s per solve       |
| `solve_zoh` (variante a)| ~1.4e8                 | ~14 s per solve        |

E se qualcuno aggiungesse un blocco `arguments` che costasse, poniamo, 1
microsecondo per chiamata (ordine di grandezza **plausibile ma non misurato
in questa repo** -- lo dico esplicitamente per non spacciare un numero
inventato per un risultato), gli stessi conteggi darebbero ~18 s e ~140 s di
overhead puro per singola risoluzione. Il README dichiara ~1-1.5 min per
valore di `tf` per il Task 1: significa che una regressione del genere
**raddoppierebbe o triplicherebbe** il tempo di soluzione. E lo sweep di
sensitivita' risolve tre valori di `tf`, piu' uno studio di convergenza di
griglia a N = 25/50/100 -- il fattore si moltiplica ancora.

C'e' un secondo effetto, piu' sottile. Le corse a `tf` corto terminano al cap
`MaxIterations = 1000` con optimality di 1e-3/1e-4 (README, e riga 733 per il
cap). Un solver che gira piu' lentamente non converge peggio, ma rende piu'
costoso *provare* a farlo convergere: alzare `MaxIterations` per curare quel
problema noto diventa proibitivo se ogni iterazione costa il doppio. La
performance qui non e' comodita', e' cio' che rende praticabile l'iterazione
sul metodo.

---

## Limiti onesti del benchmark

Elenco senza sconti.

1. **Non misura una chiamata isolata, misura un batch.** 1000 chiamate per
   `ode_descent`, 100 per `rk4_zoh`, 1 per `ode45`. Tre unita' diverse. I
   numeri stampati da `runperf` **non sono confrontabili fra i tre test**
   senza normalizzare a mano. Nulla nel codice fa questa normalizzazione e
   nulla la documenta oltre ai due commenti alle righe 39-40 e 51-52.

2. **Non e' rappresentativo del carico dentro `fmincon`.** Il batch ripete la
   stessa chiamata con **gli stessi identici argomenti**. Dati caldi in
   cache, branch predetti, nessuna variazione. Dentro `fmincon` le
   valutazioni a differenze finite perturbano ogni volta una componente
   diversa di un vettore da 350 elementi, e la RHS e' chiamata da dentro
   `trap_nonlcon` che fa `reshape`, alloca `f = zeros(5,N)` e cicla. Il costo
   *reale* per chiamata nel contesto d'uso include effetti di memoria e
   overhead del chiamante che questo benchmark **non** cattura. E' un
   microbenchmark di best case, e va letto come tale.

3. **Non misura il vero hot loop.** Cio' che `fmincon` chiama non e'
   `ode_descent`: e' `trap_nonlcon` / `zoh_nonlcon`. Quelle funzioni fanno
   `reshape`, allocano matrici, costruiscono i defect e chiamano
   `path_ineq`. Un benchmark di `trap_nonlcon(z, N, dt, d)` sarebbe molto
   piu' vicino alla realta' e catturerebbe regressioni (per esempio una
   allocazione inutile nel ciclo dei defect) che questo file lascia passare
   del tutto. La scelta di stare al livello della RHS e' difendibile (isola
   la primitiva), ma lascia scoperto lo strato dove sta la maggior parte del
   codice.

4. **Copertura parziale.** Non sono benchmarkati: `lti_zoh` (giustificabile:
   O(1) per solve), `ode_descent_uacc` (giustificabile: fuori dal cammino
   critico dell'NLP), ma soprattutto **`compute_ltv_zoh`** e la RHS aumentata
   `ltv_aug_rhs` (riga 279), che integrano un sistema a **70 stati** a ogni
   iterazione SCvx e sono plausibilmente la parte piu' costosa delle varianti
   b/c. Quella e' un'omissione senza giustificazione ovvia.

5. **Nessuna baseline, nessuna soglia, nessun gate.** `runperf` restituisce
   un array di `TimeResult`, ma la repo non memorizza nessun valore di
   riferimento e non contiene nessun confronto automatico (`comparisonPlot`
   o assert su soglia). In pratica una regressione viene "rilevata" solo se un
   umano lancia `runperf` a mano e si ricorda i numeri di prima. Chiamarlo
   *regression test* e' generoso: e' uno **strumento di misura**, non un
   allarme.

6. **Costanti duplicate.** `Vc` e `dt` (righe 9-10) sono ricopiate a mano da
   `nondim`. Verificate oggi: `Vc = 0.0777` torna, `dt = 0.0444` e' invece un
   doppio arrotondamento di 0.044347 (alto dello 0.1% circa). Non c'e' nulla
   che le tenga sincronizzate con `nondim`.

7. **Contaminazione minore del confine di misura.** La anonymous function
   della riga 68 e' costruita a ogni iterazione misurata, mentre `opts` (riga
   66) e' stato issato fuori. E la RHS del benchmark non replica la doppia
   anonymous function (`u_fcn(tt)`) del vero `fwd_integrate`, quindi
   sottostima quel percorso.

8. **Dipendenza dalla macchina.** Nessun benchmark e' portabile fra CPU
   diverse, frequenze di turbo, stati termici. I numeri hanno senso solo come
   confronto *prima/dopo sulla stessa macchina, nella stessa sessione*. Il
   file non lo dice.

---

## Possibili domande d'esame

**D: Che differenza c'e' fra `matlab.unittest.TestCase` e
`matlab.perftest.TestCase`, e perche' questo file eredita dalla seconda?**
R: `matlab.perftest.TestCase` estende `matlab.unittest.TestCase` aggiungendo i
*confini di misura*: `startMeasuring`/`stopMeasuring` per delimitare
esplicitamente una regione, e `keepMeasuring` per farla ripetere
automaticamente dal framework finche' la misura non e' accurata. Serve perche'
il tempo di un metodo di test *intero* e' inutile: include setup, allocazioni
e asserzioni, che possono dominare il codice di interesse di ordini di
grandezza. Solo delimitando la regione si ottiene un numero che significa
qualcosa. `runperf` inoltre accetta solo classi `perftest` per fare
l'esperimento statistico.

**D: Perche' il test misura un ciclo `for k = 1:1000` invece di una singola
chiamata a `ode_descent`?**
R: Perche' una chiamata singola (dell'ordine delle decine di nanosecondi,
il commento dice ~60 ns) sta sotto la risoluzione temporale del framework:
misurarla darebbe rumore. Ripetendola 1000 volte dentro il confine di misura,
il tempo totale diventa misurabile con precisione e l'overhead della misura
diventa trascurabile. Il prezzo e' che il numero riportato e' il tempo di 1000
chiamate, non di una, e che il benchmark diventa un best case: 1000 chiamate
identiche con dati caldi in cache non riproducono il pattern d'accesso reale
dentro `fmincon`.

**D: Perche' `ode_descent` e `rk4_zoh` non hanno un blocco `arguments`, mentre
`lti_zoh` ce l'ha? Non e' un'incoerenza?**
R: No, il criterio e' la frequenza di chiamata. `ode_descent` e' invocata
dentro `trap_nonlcon`/`zoh_nonlcon`, che `fmincon` chiama a ogni valutazione a
differenze finite: con 350 variabili e fino a 1000 iterazioni si arriva a
1e7-1e8 chiamate, e un validatore costerebbe quanto o piu' del corpo funzione
stesso. `lti_zoh` invece e' chiamata **una volta per risoluzione**, alla riga
1218 di `main_task2.m`, cioe' **prima** del ciclo SCvx che parte alla riga
1239: la dinamica GFOLD in variabili `z = ln(m)`, `u = T/m` e' esattamente LTI,
quindi le matrici discrete sono costanti sulla griglia e non vanno
ricalcolate. Essendo un boundary helper O(1), paga volentieri la validazione.
La regola della repo e' proprio "hot loop senza validazione, boundary con
validazione", e i due casi la rispettano entrambi.

**D: Come fa `runperf` a darti un numero affidabile, e cosa aggiunge rispetto a
un `timeit`?**
R: `runperf` costruisce un time experiment `limitingSamplingError`: scarta 5
esecuzioni di warm-up (che pagano JIT e cache fredde), poi raccoglie campioni
in modo **adattivo** finche' la media campionaria non e' nota entro un margine
di errore relativo del 5% al 95% di confidenza, fra un minimo di 4 e un massimo
di 256 campioni. `sampleSummary` riepiloga media, mediana, deviazione standard,
min e max. `timeit` da' invece un singolo numero senza intervallo di confidenza
e senza criterio d'arresto statistico: non ti permette di distinguere una
regressione vera dal rumore del sistema operativo, che e' esattamente la
domanda a cui un performance test deve rispondere.

**D: Cosa comporterebbe, in concreto, una regressione di prestazioni su
`ode_descent`?**
R: Un'amplificazione enorme, per via del conteggio delle chiamate. Sul percorso
trapezoidale servono ~1.8e7 valutazioni per risoluzione, su quello ZOH (con
`n_sub = 2`, quindi 8 RHS per intervallo) ~1.4e8. Un rallentamento di 100 ns per
chiamata costerebbe rispettivamente ~2 s e ~14 s per singola soluzione; e lo
script risolve tre valori di `tf` piu' uno studio di convergenza a
N = 25/50/100. Su un runtime nominale di ~1-1.5 min per `tf`, un'aggiunta di
overhead per chiamata dell'ordine del microsecondo raddoppierebbe o
triplicherebbe il tempo totale, e renderebbe proibitivo alzare
`MaxIterations = 1000` per curare la convergenza stallata a `tf` corto.

**D: Questo benchmark e' rappresentativo del carico reale dentro `fmincon`?**
R: Solo parzialmente, e va detto. Misura la primitiva in condizioni ideali:
stessi argomenti a ogni chiamata, dati in cache, nessuna variazione. Il vero
hot loop chiamato da `fmincon` e' `trap_nonlcon`/`zoh_nonlcon`, che oltre alla
RHS fanno `reshape`, allocano matrici e valutano `path_ineq`, e che ricevono
ogni volta un vettore perturbato diverso. Un benchmark di `nonlcon` sarebbe piu'
fedele. Inoltre il test `ode45` usa una RHS piu' semplice di quella vera di
`fwd_integrate` (che passa per una anonymous function `u_fcn(tt)` in piu' a ogni
valutazione), quindi sottostima quel percorso. Il file misura la primitiva, non
il sistema.

**D: Perche' c'e' un `verifySize` dentro un test di performance?**
R: Sta **fuori** dal confine di misura (dopo il `while`), quindi non inquina il
numero. Serve a due cose: assicurare che il risultato venga effettivamente
consumato (igiene da microbenchmark, evita che il calcolo possa essere
considerato codice morto) e fare da guardrail minimo -- se `ode_descent`
smettesse di restituire un `5x1`, staresti misurando la velocita' di una
funzione sbagliata. Non e' un test di correttezza: quella vive in
`tests/odeDescentTest.m` e `tests/rk4ZohTest.m`.
