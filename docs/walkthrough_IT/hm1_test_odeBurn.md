# HM1/tests/odeBurnTest.m

## Ruolo del file nel progetto

E' la **suite di unit test** di `HM1/ode_burn.m`, cioe' del cuore numerico di
tutto l'homework: la RHS che integra insieme stato e legge di controllo ottima
(*linear tangent*) per l'arco propulso. Tutti e quattro i `main_task*.m`
chiamano `ode_burn` dentro `ode45`, dentro il residuo di shooting, dentro
`fsolve`. Se `ode_burn` ha un errore di segno o una formula sbagliata, **tutti**
i risultati di HM1 sono sbagliati, e lo sono in modo silenzioso: `fsolve`
converge lo stesso, semplicemente al problema sbagliato.

La classe eredita da `matlab.unittest.TestCase` e contiene **sei test**. Si
lancia con `runtests('HM1/tests')` secondo la convenzione del repo. Il commento
di intestazione dichiara la copertura: "Covers the linear-tangent steering law,
the costate equations, and two analytic limits (ballistic flight, vertical
Tsiolkovsky burn). Non-dimensional (g = 1)" (righe 2-5).

La strategia dei test e' su **due livelli**, ed e' la parte interessante da
difendere all'orale:

1. **Test puntuali sulla RHS** (test 1-4): una sola valutazione di `ode_burn`,
   con parametri scelti in modo che il valore atteso sia calcolabile a mano.
   Verificano *forma* del vettore, kinematica, legge di steering, equazione del
   costato di massa.
2. **Test di limite analitico** (test 5-6): si **integra davvero** con `ode45` e
   si confronta la traiettoria con una soluzione in forma chiusa ottenuta
   degenerando il problema (spinta nulla -> caduta balistica; spinta verticale
   -> Tsiolkovsky con perdita gravitazionale). Sono i test con vero potere
   diagnostico, perche' l'oracolo e' indipendente dal codice.

Per riferimento, la RHS sotto test (`ode_burn.m`, righe 16-33) e':

    lam_vx = p.lam_vx0
    lam_vy = p.lam_vy0 - p.lam_y * t
    phi    = atan2(lam_vy, lam_vx)

    dz = [ vx ;
           vy ;
           (T/m)*cos(phi) ;
           (T/m)*sin(phi) - 1 ;     (g = 1 nondim)
          -Q ;
           (T/m^2)*|lam_v| ]        (dlam_m/dt)

---

## `addHm1ToPath` -- TestClassSetup (righe 7-12)

```matlab
methods (TestClassSetup)
    function addHm1ToPath(testCase)
        hm1 = fileparts(fileparts(mfilename('fullpath')));
        testCase.applyFixture( ...
            matlab.unittest.fixtures.PathFixture(hm1));
    end
end
```

- **Riga 9 -- risoluzione del percorso.** `mfilename('fullpath')` da'
  `.../HM1/tests/odeBurnTest`; il primo `fileparts` toglie il nome file
  (`.../HM1/tests`), il secondo sale di un livello (`.../HM1`). E' un modo
  **relativo al file, non alla working directory**: la suite funziona ovunque sia
  lanciata, il che e' precisamente cio' che serve per `runtests('HM1/tests')`
  eseguito dalla root della repo.

- **Riga 10 -- `PathFixture`.** Aggiunge `HM1/` al path MATLAB *per la durata
  della classe di test* e lo **rimuove automaticamente** alla fine (teardown
  garantito, anche se un test lancia un'eccezione). Senza questo, `ode_burn`
  non sarebbe visibile dalla cartella `tests/`. Usare la fixture invece di un
  `addpath` nudo evita di inquinare il path dell'utente: e' il motivo per cui il
  framework fornisce le fixture.

- `TestClassSetup` (non `TestMethodSetup`): gira **una volta sola** per l'intera
  classe. Corretto -- il path non cambia fra un test e l'altro.

---

## `testKinematicsAndMassFlow` (righe 15-25)

```matlab
p = struct('T', 1.2, 'Q', 2, 'c', 0.6, ...
           'lam_vx0', 1, 'lam_vy0', 0.5, 'lam_y', 0.2);
z = [0.1; 0.2; 0.3; 0.4; 0.9; 1.1];
dz = ode_burn(0.5, z, p);
testCase.verifySize(dz, [6 1]);
testCase.verifyEqual(dz(1), z(3), 'AbsTol', 1e-15);
testCase.verifyEqual(dz(2), z(4), 'AbsTol', 1e-15);
testCase.verifyEqual(dz(5), -p.Q, 'AbsTol', 1e-15);
```

- **Riga 21 -- `verifySize(dz, [6 1])`.** Il test piu' banale e uno dei piu'
  utili. `ode45` richiede un **vettore colonna**: se `ode_burn` restituisse un
  `1x6` (riga), `ode45` fallirebbe con un errore criptico sulle dimensioni. Ma
  soprattutto: le 6 componenti sono `[x; y; vx; vy; m; lam_m]` -- se qualcuno
  aggiungesse o togliesse uno stato, tutti i `main_task*.m` (che costruiscono a
  mano un vettore iniziale a 6 componenti: `ic = [0;0;0;0;1;1]` in Task 1 e 4,
  `ic = [p.x0; p.y0; p.vx0; p.vy0; p.m0; 1]` in Task 2, con l'ultima componente
  lasciata incognita in Task 3) si romperebbero. Il test blocca il **contratto
  d'interfaccia**.

- **Righe 22-23 -- kinematica.** `dx/dt = vx`, `dy/dt = vy`: le prime due righe
  della RHS devono essere **copie esatte** di due componenti dello stato. La
  tolleranza `1e-15` (non `0`) e' formalmente corretta, ma di fatto qui il
  confronto e' bit-esatto. E' un test *anti-refuso*: protegge contro un
  `dz(1) = z(4)` da copia-incolla.

- **Riga 24 -- portata di massa.** `dm/dt = -Q`, **costante**. E' un fatto
  modellistico non ovvio e vale la pena difenderlo: nel modello di HM1 la portata
  e' costante (motore a spinta fissa), quindi la massa decresce **linearmente**
  in `t`, e vale la relazione `Q = T/c` (infatti nei dati del test
  `T/c = 1.2/0.6 = 2 = Q`, coerente). La conseguenza numerica e' importante:
  l'equazione di massa e' esattamente integrabile, e infatti la forma integrata
  `m(t) = m0 - Q*t` viene usata come **check di conservazione** nei test 5-6 di
  questa suite (righe 73 e 88) e nel performance test
  (`odeBurnPerformanceTest.m`, riga 51). Qui, invece, si asserisce solo la
  derivata `dz(5) = -Q`, senza integrare.

- **Il punto della scelta dei parametri**: `lam_vx0`, `lam_vy0`, `lam_y` sono
  tutti **non nulli** (`1`, `0.5`, `0.2`) e `t = 0.5` e' non nullo. Quindi il
  test afferma che queste tre componenti della RHS sono **indipendenti dai
  costati e dal tempo**: qualunque cosa faccia la legge di steering, `dx/dt` resta
  `vx`. E' un'affermazione strutturale, non un caso particolare.

> **Possibile domanda d'esame** -- Perche' testare `verifySize` invece di
> lasciare che sia `ode45` a lamentarsi?
> *Risposta:* Perche' l'errore di `ode45` arriva in fondo a una catena
> `fsolve -> residual -> ode45` e dice qualcosa come "the size of the derivative
> vector must be equal to the size of the initial condition" -- non dice *dove* e'
> il bug. Il test lo localizza in `ode_burn` in un millisecondo e senza dipendere
> dal solver. Piu' in generale, la dimensione dello stato e' un contratto pubblico
> di `ode_burn`: tutti i `main_task*.m` lo assumono quando costruiscono a mano un
> `ic` a 6 componenti (`[0;0;0;0;1;1]` in Task 1 e 4,
> `[p.x0; p.y0; p.vx0; p.vy0; p.m0; 1]` in Task 2).

---

## `testThrustAlongConstantCostate` (righe 27-36)

```matlab
% lam = (1, 0) constant -> phi = 0: horizontal thrust, dvy = -g
p = struct('T', 1.2, 'Q', 2, 'c', 0.6, ...
           'lam_vx0', 1, 'lam_vy0', 0, 'lam_y', 0);
m = 0.8;
z = [0; 0; 0.5; 0.1; m; 1];
dz = ode_burn(0.3, z, p);
testCase.verifyEqual(dz(3), p.T/m, 'AbsTol', 1e-14);
testCase.verifyEqual(dz(4), -1,    'AbsTol', 1e-14);
```

- **Righe 29-30 -- la costruzione.** `lam_vy0 = 0` **e** `lam_y = 0`, quindi
  `lam_vy(t) = 0 - 0*t = 0` per ogni `t`. Con `lam_vx0 = 1 > 0`, l'angolo di
  spinta e'

      phi = atan2(0, 1) = 0

  cioe' **spinta perfettamente orizzontale**.

- **Righe 34-35 -- cosa si verifica.** Con `phi = 0`: `cos(phi) = 1`,
  `sin(phi) = 0`, quindi

      dvx/dt = T/m           = 1.2/0.8 = 1.5
      dvy/dt = 0 - 1         = -1      (sola gravita')

  Il primo assert isola la **normalizzazione dell'accelerazione di spinta**
  (`T/m`, non `T`, non `T*m`): verifica che la massa sia a denominatore. Il
  secondo isola il **termine gravitazionale**: nel modello nondimensionale
  `g = 1`, hard-coded come `- 1` in `ode_burn.m` riga 31. Se qualcuno rimettesse
  `9.81` (dimensionale), o cambiasse segno, questo test si accorgerebbe subito.

  **Perche' e' l'invariante giusto**: e' il caso in cui la legge di controllo
  degenera a un valore banale (`phi = 0`) noto senza calcoli, il che permette di
  testare *separatamente* i due termini di `dvy/dt` -- il contributo di spinta
  (che si annulla) e quello di gravita' (che resta). Con un `phi` generico i due
  termini sarebbero mescolati e il test perderebbe potere diagnostico.

- Il tempo di valutazione e' `t = 0.3` (non zero) ma con `lam_y = 0` non ha
  effetto: e' proprio il test successivo a occuparsi della dipendenza dal tempo.

---

## `testLinearTangentSwitch` (righe 38-49)

```matlab
% lam_vy(t) = lam_vy0 - lam_y*t vanishes at t* = lam_vy0/lam_y,
% where the thrust must be exactly horizontal (phi = 0)
p = struct('T', 1.2, 'Q', 2, 'c', 0.6, ...
           'lam_vx0', 1, 'lam_vy0', 2, 'lam_y', 4);
tStar = p.lam_vy0 / p.lam_y;
...
dz = ode_burn(tStar, z, p);
testCase.verifyEqual(dz(3), p.T/m, 'AbsTol', 1e-14);
testCase.verifyEqual(dz(4), -1,    'AbsTol', 1e-14);
```

- **Questo e' il test piu' intelligente della suite** e va capito bene, perche'
  isola l'unica formula davvero "teorica" della RHS.

- **Il contesto teorico.** Le equazioni dei costati per il problema planare di
  HM1 si ottengono da `lam_dot = -dH/dx`. Siccome `x` non compare nella dinamica
  e `x(tf)` e' libero, `lam_x = 0` identicamente. Siccome `y` non compare nella
  dinamica (Terra piatta, `g` costante, no drag), `lam_y` e' **costante**. Da
  `lam_vx_dot = -lam_x = 0` segue che `lam_vx` e' **costante**; da
  `lam_vy_dot = -lam_y` segue che `lam_vy` e' **lineare nel tempo**:

      lam_vy(t) = lam_vy0 - lam_y * t

  L'angolo ottimo (che massimizza `lam_v . u_hat`) e' allineato al primer vector:

      tan(phi(t)) = lam_vy(t)/lam_vx = (lam_vy0 - lam_y*t)/lam_vx0

  cioe' la **tangente dell'angolo di spinta e' lineare nel tempo**: e' la
  *linear tangent law*, il risultato classico per ascesa a Terra piatta senza
  drag. (E' il caso particolare della *bilinear tangent law* quando `lam_x = 0`.)

- **Righe 42-43 -- la costruzione del test.** `lam_vy0 = 2`, `lam_y = 4`, quindi
  `lam_vy(t)` si annulla in

      t* = lam_vy0/lam_y = 2/4 = 0.5

  Il test valuta `ode_burn` **esattamente in `t*`** e verifica che li' la spinta
  sia orizzontale (`dvx = T/m`, `dvy = -1`), esattamente le stesse asserzioni del
  test precedente.

- **Perche' questo e' l'invariante giusto.** L'assert e' identico a quello di
  `testThrustAlongConstantCostate`, ma il *modo* di ottenerlo e' opposto: qui
  `lam_vy0 = 2 != 0` e sono la **dipendenza dal tempo** e il valore di `lam_y` a
  produrre `lam_vy(t*) = 0`. Quindi il test fallisce **se e solo se** il termine
  `- p.lam_y * t` e' sbagliato: se lo si cancellasse, in `t*` si avrebbe
  `lam_vy = 2`, `phi = atan2(2, 1) ~= 63.4 gradi`, e `dvx` scenderebbe da `1.2` a
  circa `0.54`. Se il segno fosse `+`, `lam_vy(t*) = 4` e sarebbe ancora peggio.
  E' un test **chirurgico** su un singolo termine, costruito sfruttando uno **zero
  noto** della funzione. Questa e' la tecnica generale: cercare i punti dove il
  modello degenera a un valore calcolabile a mano.

- Il "switch" nel nome e' un po' fuorviante: **non** e' lo switch di
  spegnimento motore di Task 3 (funzione di switching `S = |lam_v|/m - lam_m/c`).
  E' semplicemente il passaggio a zero di `lam_vy`, cioe' l'istante in cui il
  vettore spinta attraversa l'orizzontale ruotando verso il basso -- il tipico
  *pitch-over* dell'ascesa ottima.

> **Possibile domanda d'esame** -- Da dove viene la *linear tangent law* e cosa la
> distrugge?
> *Risposta:* Viene dal fatto che, con Terra piatta, gravita' costante e nessun
> drag, la dinamica non dipende ne' da `x` ne' da `y`. Quindi `lam_x = 0` e
> `lam_y` costante, da cui `lam_vx` costante e `lam_vy` lineare, e
> `tan(phi) = lam_vy/lam_vx` lineare in `t`. Basta introdurre il drag (che
> dipende da `y` tramite la densita', e da `v`) o la gravita' inversa-quadrato
> (che dipende da `r`) perche' `lam_y` smetta di essere costante: le equazioni dei
> costati si accoppiano allo stato e vanno integrate numericamente insieme ad
> esso. La legge lineare-tangente e' quindi una **proprieta' del modello
> semplificato**, non una legge generale -- ed e' esattamente il motivo per cui
> `ode_burn` puo' permettersi di **non** integrare `lam_vx` e `lam_vy` (li
> ricostruisce in forma chiusa da `p` e `t`), integrando solo 6 stati invece di 9.

---

## `testCostateMassEquation` (righe 51-59)

```matlab
% dlam_m/dt = (T/m^2) * |lam_v|, with |lam_v| = 5 for lam = (3,4)
p = struct('T', 1.2, 'Q', 2, 'c', 0.6, ...
           'lam_vx0', 3, 'lam_vy0', 4, 'lam_y', 0);
m = 0.5;
z = [0; 0; 0; 0; m; 1];
dz = ode_burn(0, z, p);
testCase.verifyEqual(dz(6), p.T/m^2 * 5, 'AbsTol', 1e-12);
```

- **Righe 53-54 -- la terna pitagorica.** `lam_vx0 = 3`, `lam_vy0 = 4` danno
  `|lam_v| = sqrt(9 + 16) = 5` **esatto**, senza errori di arrotondamento nel
  valore atteso. Trucco deliberato: il valore di riferimento (`1.2/0.25 * 5 = 24`)
  e' rappresentabile esattamente, quindi la tolleranza `1e-12` misura solo
  l'errore del codice, non quello dell'oracolo.

- **La derivazione dell'equazione.** L'Hamiltoniana del problema di ascesa e'

      H = lam_x*vx + lam_y*vy + lam_vx*(T/m)*cos(phi)
          + lam_vy*((T/m)*sin(phi) - 1) - lam_m*Q

  Con `phi` scelto ottimamente, il termine di spinta vale
  `(T/m)*(lam_vx*cos + lam_vy*sin) = (T/m)*|lam_v|`. Allora

      dH/dm = -(T/m^2)*|lam_v|

  e per l'equazione dei costati `lam_m_dot = -dH/dm`:

      lam_m_dot = +(T/m^2)*|lam_v|

  che e' esattamente la riga 33 di `ode_burn.m`. **Segno positivo**: `lam_m`
  **cresce** durante il burn. Interpretazione: `lam_m` e' il valore marginale di
  un chilo di massa; man mano che la massa cala, `T/m` cresce e la stessa massa
  residua "vale di piu'". A `tf`, per trasversalita' con costo `-m(tf)`, si ha
  `lam_m(tf) = 1`. Nei Task 1, 2 e 4 la normalizzazione del repo pone invece
  `lam_m0 = 1`, sfruttando l'omogeneita' di grado 1 di `H` nei costati; in
  Task 3 (arco di coast) `lam_m0` resta un'incognita del BVP e la normalizzazione
  e' imposta al cutoff, `lam_m(tc) = 1` (`main_task3.m`, righe 34 e 39).

- **Onesta' sul potere del test.** Questo assert e' un **mirror della formula**:
  ricopia in MATLAB la stessa espressione che il codice implementa. Protegge da
  regressioni (se domani qualcuno scrive `m^3` invece di `m^2`, il test cade), ma
  **non verifica la derivazione**: se il segno fosse sbagliato *nella teoria*,
  test e codice sarebbero sbagliati insieme e il test passerebbe. E' quello che si
  chiama un *change detector*. I test 5 e 6, che confrontano contro soluzioni
  analitiche indipendenti, sono gli unici a poter smentire il codice.

- **Nota**: `lam_y = 0` fa si' che `|lam_v|` sia costante in `t`, e il test
  valuta comunque in `t = 0`. Quindi il test **non** copre l'interazione fra la
  variazione temporale di `lam_vy` e `|lam_v|` -- cioe' non verifica che
  `|lam_v(t)|` sia ricalcolata ad ogni passo. Un test aggiuntivo con `lam_y != 0`
  e `t != 0` chiuderebbe questa lacuna.

---

## `testBallisticLimitMatchesAnalytic` (righe 61-74)

```matlab
% T = Q = 0: pure ballistic flight under unit gravity
p = struct('T', 0, 'Q', 0, 'c', 0.6, ...
           'lam_vx0', 1, 'lam_vy0', 1, 'lam_y', 0);
z0 = [0; 0; 0.3; 0.5; 1; 1];  tf = 0.8;
opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
[~, Z] = ode45(@(t,z) ode_burn(t, z, p), [0 tf], z0, opts);
testCase.verifyEqual(Z(end,1), z0(3)*tf,            'AbsTol', 1e-9);
testCase.verifyEqual(Z(end,2), z0(4)*tf - 0.5*tf^2, 'AbsTol', 1e-9);
testCase.verifyEqual(Z(end,3), z0(3),               'AbsTol', 1e-10);
testCase.verifyEqual(Z(end,4), z0(4) - tf,          'AbsTol', 1e-9);
testCase.verifyEqual(Z(end,5), z0(5),               'AbsTol', 1e-12);
```

- **Il primo test che integra davvero.** Fin qui si e' valutata la RHS in un
  punto; qui si chiama `ode45` sull'intervallo `[0, 0.8]` con **le tolleranze
  vere del repo** (`RelTol = 1e-10`, `AbsTol = 1e-12`, riga 67). Quindi si sta
  testando anche il **cablaggio** `@(t,z) ode_burn(t,z,p)`, cioe' esattamente il
  modo in cui i `main_task*.m` lo usano.

- **Righe 63-64 -- la degenerazione.** Ponendo `T = 0` e `Q = 0` la RHS collassa a

      dx/dt = vx,  dy/dt = vy,  dvx/dt = 0,  dvy/dt = -1,  dm/dt = 0

  cioe' **moto parabolico sotto gravita' unitaria**. Nota: con `T = 0` anche
  `dlam_m/dt = 0`, ma `lam_m` non viene asserito.

- **Righe 69-73 -- le cinque asserzioni** sono la soluzione esatta del moto
  balistico, verificabile a mano:

      x(tf)  = vx0*tf                = 0.3 * 0.8       = 0.24
      y(tf)  = vy0*tf - 0.5*tf^2     = 0.4 - 0.32      = 0.08
      vx(tf) = vx0                   = 0.3     (costante)
      vy(tf) = vy0 - tf              = 0.5 - 0.8       = -0.3
      m(tf)  = m0                    = 1       (nessun consumo)

  Notare che `vy(tf) < 0`: il proiettile e' gia' nella fase discendente. Non e'
  un problema, e' solo aritmetica -- ma dimostra che il test non e' stato scelto
  per far "sembrare" fisica la traiettoria, e' scelto per avere numeri esatti.

- **Perche' e' l'invariante giusto.** E' il **check di consistenza gravitazionale
  e cinematica** in un caso dove la spinta non c'e': isola il termine `-1` nella
  riga di `dvy/dt` e la doppia integrazione `y -> vy -> -g`. Un errore di fattore
  2 sulla gravita', o un `+1` invece di `-1`, produce un errore di ordine `0.1`
  su `y(tf)`, enormemente sopra la tolleranza `1e-9`.

- **La scelta delle tolleranze di asserzione**: `1e-9` sugli stati integrati,
  `1e-10` su `vx` (che non e' integrato, e' costante), `1e-12` su `m` (idem). Sono
  coerenti col fatto che `ode45` a `RelTol = 1e-10` accumula un errore di
  troncamento di quell'ordine. Le tolleranze **non** sono uniformi: chi le ha
  scritte ha ragionato su quanto errore ciascuna quantita' puo' accumulare. E'
  esattamente cio' che va fatto -- mettere `1e-12` su `y(tf)` renderebbe il test
  fragile senza guadagno diagnostico.

> **Possibile domanda d'esame** -- Se `T = 0` la traiettoria non e' piu' ottima:
> che senso ha testare un caso non ottimo?
> *Risposta:* Il test non verifica l'ottimalita', verifica la **dinamica**.
> `ode_burn` implementa due cose insieme: la dinamica del veicolo e la legge di
> controllo ottima. Spegnendo la spinta si mette a nudo la sola dinamica, che ha
> soluzione chiusa e puo' quindi essere validata contro un oracolo esatto. La
> legge di controllo e' validata separatamente dai test 2-3. E' la strategia
> corretta: **isolare un sottosistema alla volta**, invece di testare tutto
> insieme dove nessun oracolo e' disponibile.

---

## `testVerticalBurnTsiolkovskyWithGravity` (righe 76-89)

```matlab
% lam = (0, 1) constant -> phi = 90 deg: vertical burn.
% Analytic: vy(t) = (T/Q)*ln(m0/(m0 - Q*t)) - t  with m0 = 1
p = struct('T', 1.2, 'Q', 2, 'c', 0.6, ...
           'lam_vx0', 0, 'lam_vy0', 1, 'lam_y', 0);
z0 = [0; 0; 0; 0; 1; 1];  tf = 0.3;
opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
[~, Z] = ode45(@(t,z) ode_burn(t, z, p), [0 tf], z0, opts);
vyAnalytic = (p.T/p.Q) * log(1/(1 - p.Q*tf)) - tf;
testCase.verifyEqual(Z(end,4), vyAnalytic, 'AbsTol', 1e-9);
testCase.verifyEqual(Z(end,3), 0,          'AbsTol', 1e-12);
testCase.verifyEqual(Z(end,5), 1 - p.Q*tf, 'AbsTol', 1e-12);
```

- **Il test piu' forte della suite.** E' l'unico che valida la RHS **completa**
  (spinta accesa + consumo di massa + gravita') contro un oracolo analitico
  indipendente.

- **Righe 79-80 -- la degenerazione.** `lam_vx0 = 0`, `lam_vy0 = 1`, `lam_y = 0`
  danno `phi = atan2(1, 0) = pi/2`: **spinta puramente verticale e costante**.
  La RHS collassa a

      dvx/dt = (T/m)*cos(pi/2) = 0
      dvy/dt = (T/m)*sin(pi/2) - 1 = T/m - 1
      dm/dt  = -Q      ->  m(t) = m0 - Q*t = 1 - Q*t

- **La derivazione dell'equazione di Tsiolkovsky con gravita'.** Sostituendo
  `m(t) = 1 - Q*t` in `dvy/dt`:

      dvy/dt = T/(1 - Q*t) - 1

  Integrando da `0` a `t` con `vy(0) = 0`:

      vy(t) = (T/Q)*ln(1/(1 - Q*t)) - t

  Il primo termine e' il **delta-v ideale di Tsiolkovsky** -- infatti `T/Q = c`
  (velocita' efficace di scarico; nei dati del test `1.2/2 = 0.6 = p.c`, coerente)
  e `1/(1 - Q*t) = m0/m(t)`, quindi il termine e' `c * ln(m0/m)`. Il secondo
  termine, `-t`, e' la **perdita gravitazionale** (`g*t` con `g = 1`): e'
  letteralmente quanto delta-v la gravita' si mangia durante un burn verticale di
  durata `t`.

  Riga 85: `(p.T/p.Q) * log(1/(1 - p.Q*tf)) - tf` -- corrisponde esattamente.

- **Riga 88 -- la conservazione della massa.** `m(tf) = 1 - Q*tf = 1 - 0.6 = 0.4`.
  Positivo, quindi il burn e' fisicamente ammissibile (non si esaurisce la massa
  prima di `tf`). Con `Q = 2`, `tf` non puo' superare `0.5` -- un vincolo implicito
  del test che vale la pena conoscere: `tf = 0.3` e' scelto con margine.

- **Riga 83 -- tolleranze piu' strette.** Qui `RelTol = 1e-12`, `AbsTol = 1e-14`,
  contro `1e-10 / 1e-12` del test balistico. Motivo: `1/(1 - Q*t)` diverge quando
  `Q*t -> 1`, quindi la RHS e' molto piu' "ripida" e l'errore locale
  dell'integratore cresce; per poter asserire `1e-9` su `vy(tf)` serve un margine
  numerico maggiore. E' una scelta consapevole, non un vezzo.

- **Perche' e' l'invariante giusto**: verifica **contemporaneamente** che (a) `phi`
  sia calcolato con `atan2` nel verso giusto (un `atan2(lam_vx, lam_vy)` invertito
  darebbe `phi = 0` e spinta orizzontale, e il test esploderebbe), (b) la massa
  entri a denominatore come `T/m` e non altrimenti, (c) `dm/dt = -Q`, (d) `g = 1`
  sia sottratta. Se **una qualsiasi** di queste e' sbagliata, `vy(tf)` non torna.
  E' l'unico test che accoppia tutti i pezzi.

> **Possibile domanda d'esame** -- Il test si chiama "Tsiolkovsky with gravity":
> perche' non basta l'equazione classica `dv = c*ln(m0/mf)`?
> *Risposta:* Perche' Tsiolkovsky vale in assenza di forze esterne. In salita
> verticale la gravita' agisce per tutta la durata del burn e sottrae `g*tb` al
> delta-v ottenuto, con `tb` durata del burn: questa e' la *gravity loss*. Il test
> ne fa un oracolo esatto perche' in salita verticale (e senza drag) la perdita e'
> integrabile in forma chiusa: `-g*t`, cioe' `-t` in nondimensionale. E' anche il
> motivo fisico per cui l'ascesa ottima **non** e' verticale se non per un breve
> tratto iniziale: minimizzare la gravity loss spinge a ruotare presto verso
> l'orizzontale, ed e' esattamente cio' che la linear-tangent law fa.

---

## Lacune di copertura (onesta' richiesta)

Cose che la suite **non** verifica, e che vale la pena saper elencare all'orale:

- **Conservazione dell'Hamiltoniana.** Il sistema aumentato e' **autonomo** (la
  dipendenza esplicita da `t` in `ode_burn` e' solo la ricostruzione in forma
  chiusa di `lam_vy(t)`, non una vera non autonomia), quindi lungo la traiettoria
  ottima `H` deve essere **costante**, e per il caso a tempo finale libero deve
  valere `H = 0`. Ricostruendo i costati da `p` e `t` si potrebbe calcolare `H` a
  ogni punto della soluzione di `ode45` e asserire `max|H(t) - H(0)| < tol`.
  **Questo test non esiste nella suite.** Sarebbe l'invariante piu' potente
  disponibile, e coprirebbe simultaneamente steering law, costati e dinamica.

- **Confronto derivata analitica vs differenze finite.** Nessun test confronta
  `ode_burn` con una derivata numerica dello stato (`(z(t+h) - z(t))/h`). Sarebbe
  ridondante rispetto ai test 5-6, ma non e' presente.

- **Casi limite / robustezza numerica.** Non c'e' nessun test su `m -> 0`
  (dove `T/m` e `T/m^2` divergono), ne' su `lam_vx = lam_vy = 0` (dove
  `atan2(0,0) = 0` per convenzione MATLAB, ma il primer vector e' indefinito e la
  legge di controllo perde senso), ne' su input di dimensione sbagliata (che con
  la scelta di **non** avere un blocco `arguments` produrrebbero errori oscuri).

- **Test 1-4 sono `change detector`**: ricopiano le formule del codice. Solo i
  test 5-6 hanno un oracolo indipendente.

---

## Possibili domande d'esame

**D: Che cosa rende un test "buono" per una RHS di ODE, e come si vede in questa
suite?**
R: Un test buono confronta il codice con un **oracolo indipendente**, non con se
stesso. Nella suite ci sono due famiglie: i test 1-4 valutano la RHS in un punto e
confrontano con la formula ricopiata a mano (proteggono da regressioni, ma se la
teoria fosse sbagliata passerebbero comunque); i test 5-6 **degenerano** il
problema a un caso con soluzione chiusa (balistica, Tsiolkovsky) e integrano
davvero con `ode45`. Solo questi ultimi possono smentire il codice. La tecnica
generale e' proprio quella: **spegnere termini** (`T = 0`, `Q = 0`) o **allineare i
costati** (`lam = (1,0)`, `lam = (0,1)`) per far collassare il modello su qualcosa
che si sa risolvere a mano.

**D: Perche' `testLinearTangentSwitch` valuta la RHS proprio in
`t* = lam_vy0/lam_y`?**
R: Perche' li' `lam_vy(t*) = 0` per costruzione, quindi `phi = atan2(0, lam_vx0) = 0`
e la spinta e' esattamente orizzontale: l'atteso e' calcolabile a mano
(`dvx = T/m`, `dvy = -1`). E' il modo per testare **chirurgicamente** il termine
`- lam_y*t` -- l'unica parte della RHS che dipende dal tempo. Se quel termine
sparisse o cambiasse segno, in `t*` la spinta non sarebbe piu' orizzontale e il
test cadrebbe con un errore grosso. Testare in un `t` generico non permetterebbe
di scrivere l'atteso senza rifare gli stessi conti del codice.

**D: Perche' `lam_vx0 = 3`, `lam_vy0 = 4` e non due numeri qualsiasi?**
R: Terna pitagorica: `|lam_v| = 5` esatto. Il valore atteso del test
(`T/m^2 * 5`) e' quindi rappresentabile esattamente in floating point, e la
tolleranza `1e-12` misura solo l'errore del codice sotto test, non quello
dell'oracolo. E' una piccola accortezza che rende il test piu' pulito da leggere e
piu' stretto da asserire.

**D: L'equazione `dlam_m/dt = +(T/m^2)*|lam_v|` ha segno positivo. E' corretto? Da
dove viene?**
R: Si'. Viene da `lam_m_dot = -dH/dm` con
`H = ... + (T/m)*|lam_v| - lam_m*Q` (avendo gia' sostituito il `phi` ottimo, che
rende il prodotto `lam_v . u_hat` pari a `|lam_v|`). Derivando, `dH/dm = -(T/m^2)*|lam_v|`,
e il segno meno dell'equazione dei costati lo ribalta in `+`. Quindi `lam_m` **cresce**
durante il burn. Interpretazione fisica: `lam_m` e' il valore marginale della massa;
man mano che il veicolo si alleggerisce l'accelerazione `T/m` cresce, e ogni chilo
residuo diventa piu' prezioso. La condizione di trasversalita' del problema a massa
finale massima da' `lam_m(tf) = 1`, e nei Task 1, 2 e 4 HM1 sfrutta l'omogeneita' di
grado 1 di `H` nei costati per normalizzare invece `lam_m0 = 1`. In Task 3 `lam_m0`
e' una delle incognite del BVP e la normalizzazione `lam_m = 1` e' imposta al cutoff
(`lam_m(tc) = 1`, dove `lam_m` e' costante durante il coast).

**D: Come testeresti che la traiettoria trovata e' davvero ottima, non solo che
l'integrazione e' corretta?**
R: Con l'**invariante Hamiltoniano**, che oggi manca nella suite. Il sistema
aumentato e' autonomo, quindi lungo una traiettoria che soddisfa le condizioni di
Pontryagin `H` deve essere costante; e siccome in HM1 il tempo finale e' libero, la
condizione di trasversalita' impone `H = 0` per tutto l'arco. Ricostruendo `lam_vx`,
`lam_vy(t)` e `lam_y` da `p`, e prendendo `lam_m` dalla sesta componente dello stato
integrato, si puo' valutare `H` su tutta la griglia di `ode45` e asserire
`max|H(t)| < 1e-8` sulla soluzione convergente di `fsolve`. Sarebbe l'unico test che
verifica **l'ottimalita'**, non solo la dinamica -- tutti gli altri passerebbero anche
con costati sbagliati, purche' la RHS li propaghi coerentemente.

**D: Perche' `PathFixture` invece di un semplice `addpath`?**
R: Perche' la fixture garantisce il **teardown**: MATLAB rimuove `HM1/` dal path alla
fine della classe di test anche se un test lancia un'eccezione. Un `addpath` nudo
lascerebbe il path dell'utente sporco, e in una suite piu' grande potrebbe far
"funzionare" test che invece dovrebbero fallire per file non trovati (contaminazione
fra test). E' la stessa filosofia del RAII: l'acquisizione della risorsa e' legata a
uno scope.
