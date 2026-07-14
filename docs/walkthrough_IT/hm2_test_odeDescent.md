# HM2_powered_descent/tests/odeDescentTest.m

## Ruolo del file nel progetto

Questa e' la suite di unit test (class-based `matlab.unittest`) del cuore
dinamico di HM2: `ode_descent.m`, il right-hand side (RHS) adimensionale del
problema di powered descent. La funzione sotto test e' minuscola -- due righe di
codice eseguibile -- ma e' il pezzo su cui poggia TUTTO il resto dell'homework:
la collocazione trapezoidale la usa per costruire i defect, il propagatore
`rk4_zoh.m` la chiama quattro volte per substep, la linearizzazione LTV di SCvx
ne calcola le Jacobiane analitiche, e la replay open-loop con `ode45` la integra
per misurare la fedelta' della soluzione. Un errore di segno o di scala qui non
farebbe fallire nessun solver: produrrebbe semplicemente una traiettoria
ottimale... del problema sbagliato. Da qui l'esigenza di un test che inchiodi la
definizione della derivata prima ancora di parlare di ottimizzazione.

La funzione testata (`HM2_powered_descent/ode_descent.m`, righe 13-14) e':

    Tmag = sqrt(u(1)^2 + u(2)^2)
    dx   = [ x(3); x(4); u(1)/x(5); u(2)/x(5) - 1; -Vc*Tmag ]

cioe', in forma matematica, il punto materiale 2D a Terra piatta, senza
aerodinamica, in variabili adimensionali:

    x_dot  = vx
    y_dot  = vy
    vx_dot = Tx/m
    vy_dot = Ty/m - 1
    m_dot  = -Vc * ||T||

Il "-1" nella riga di vy e' la gravita' adimensionalizzata (g = 1 per costruzione
dello scaling: l'accelerazione di riferimento e' g stessa). `Vc = V_ref/c` e' il
numero di Tsiolkovsky, con `c = Isp*g0` velocita' efficace di scarico; nel codice
di produzione vale 0.0777 (vedi `main_task2.m`).

Attenzione a una convenzione che ritorna in tutta HM2: qui il **controllo e' la
spinta T**, non l'accelerazione. Il file gemello `ode_descent_uacc.m` tiene
costante `u = T/m` ed e' quello usato dalla variante GFOLD log-massa. Le due RHS
NON sono la stessa cosa sotto ZOH (ci torniamo nella pagina di `rk4_zoh`), e
questa suite documenta solo la prima.

I quattro test scelgono quattro stati/controlli in cui la risposta esatta e'
scrivibile a mano. Non c'e' nessun confronto con differenze finite, e vedremo
perche' e' la scelta giusta per una RHS (mentre sarebbe la scelta giusta -- e
manca -- per le Jacobiane analitiche di SCvx).

---

## `TestClassSetup` -- `addHm2ToPath` (righe 5-10)

```matlab
hm2 = fileparts(fileparts(mfilename('fullpath')));
testCase.applyFixture( ...
    matlab.unittest.fixtures.PathFixture(hm2));
```

- Righe 6-9: `mfilename('fullpath')` restituisce il path assoluto di QUESTO file
  di test (`.../HM2_powered_descent/tests/odeDescentTest`); due `fileparts`
  annidati risalgono di due livelli e danno la cartella `HM2_powered_descent/`,
  dove vivono `ode_descent.m` e `rk4_zoh.m`.
- Riga 8: `PathFixture` aggiunge quella cartella al MATLAB path **e la rimuove
  automaticamente al teardown** della classe. E' la differenza fra un test
  pulito e un `addpath` che sporca l'ambiente dell'utente: dopo `runtests` il
  path e' esattamente quello di prima.
- Il blocco e' `TestClassSetup` e non `TestMethodSetup`: la fixture serve una
  volta sola per l'intera classe (i quattro test condividono lo stesso path).
  Metterla per-metodo la applicherebbe e rimuoverebbe quattro volte, senza alcun
  guadagno.
- Effetto pratico: si puo' lanciare `runtests('HM2_powered_descent/tests')` dalla
  root della repo senza `cd` preventivo. E' anche la ragione per cui i test
  girano identici in CI headless (`matlab -batch`).

> **Possibile domanda d'esame** -- Perche' usare `PathFixture` invece di un
> semplice `addpath` in cima al test?
> *Risposta:* Perche' la fixture e' transazionale: il framework garantisce il
> ripristino del path anche se un test fallisce o solleva un'eccezione. Un
> `addpath` "nudo" lascerebbe la cartella nel path del MATLAB dell'utente,
> creando dipendenze invisibili fra test successivi (un test potrebbe passare
> solo perche' un test precedente aveva aggiunto una cartella). L'isolamento fra
> test e' una proprieta' che si progetta, non che si spera.

---

## `testDerivativeDefinition` (righe 13-21)

```matlab
x  = [1; 2; 0.3; -0.4; 0.8];
u  = [0.6; 0.8];                 % |u| = 1
Vc = 0.0777;
dx = ode_descent(x, u, Vc);
expected = [0.3; -0.4; 0.6/0.8; 0.8/0.8 - 1; -Vc];
testCase.verifyEqual(dx, expected, 'AbsTol', 1e-15);
```

- Righe 15-17: lo stato e' `[x; y; vx; vy; m] = [1; 2; 0.3; -0.4; 0.8]` e il
  controllo `u = [0.6; 0.8]`. La scelta della terna pitagorica 3-4-5 riscalata
  non e' casuale: `||u|| = sqrt(0.36 + 0.64) = 1` esattamente, il che permette di
  scrivere la riga di massa attesa come `-Vc` puro, senza dover ricalcolare una
  radice a mano nel test. Un test il cui valore atteso richiede una formula
  complicata e' un test che puo' sbagliare quanto il codice.
- Riga 19: `expected` e' la derivata calcolata **a mano** dalla definizione:
  `[vx; vy; Tx/m; Ty/m - 1; -Vc*||T||] = [0.3; -0.4; 0.75; 0; -0.0777]`. Nota che
  `0.6/0.8` e `0.8/0.8` sono lasciati come espressioni: il test dichiara la
  formula, non il numero, cosi' un lettore vede immediatamente quale riga della
  dinamica sta verificando.
- Riga 20: `AbsTol 1e-15` invece di `verifyEqual` esatto. Serve: `0.6^2 = 0.36`
  non e' rappresentabile esattamente in binario, quindi `sqrt(0.36 + 0.64)` puo'
  differire da 1.0 di un ulp (circa 2.2e-16), e la riga di massa `-Vc*Tmag` eredita
  quell'errore. La tolleranza e' scelta appena sopra la macchina epsilon: e'
  ancora un test di identita' esatta, non un test "circa uguale".
- Questo e' il test di **correttezza di definizione**: se qualcuno invertisse
  `x(3)` e `x(4)`, dimenticasse la divisione per la massa, o mettesse il segno
  sbagliato alla gravita', qui si vede subito.

Un dettaglio non intenzionale ma vero: con `Ty = 0.8` e `m = 0.8`, `dx(4)` viene
esattamente 0, cioe' questo stato e' anche in equilibrio verticale. Il test
funziona lo stesso, ma un lettore distratto potrebbe non accorgersi che la riga
di vy sta verificando `1 - 1 = 0` e non un numero "vivo".

> **Possibile domanda d'esame** -- Perche' non confrontare la RHS con una
> derivata numerica per differenze finite? Non e' quello il test canonico?
> *Risposta:* Le differenze finite sono il test canonico per una **Jacobiana**
> (una derivata di una funzione gia' implementata), non per una RHS. `ode_descent`
> non e' la derivata di qualcosa che il codice conosce: e' il modello fisico
> stesso, un'espressione algebrica chiusa. Confrontarla con differenze finite di
> cosa? Bisognerebbe integrare la traiettoria e derivarla: si introdurrebbero
> errori di troncamento (O(h) o O(h^2)) e si perderebbe la capacita' di
> distinguere un bug da un errore di discretizzazione. Il confronto con il valore
> esatto calcolato a mano e' piu' forte: tolleranza 1e-15 invece di 1e-6. Il test
> a differenze finite serve invece per la funzione `jacobians` in `main_task2.m`
> (righe 255-277), che calcola A e B analitiche per SCvx -- e quel test, ad oggi,
> NON esiste nella suite: e' un buco di copertura reale.

---

## `testBallisticCoast` (righe 23-28)

```matlab
x  = [0.5; 1; 0.2; -0.1; 0.9];
dx = ode_descent(x, [0; 0], 0.0777);
testCase.verifyEqual(dx, [0.2; -0.1; 0; -1; 0], 'AbsTol', 1e-15);
```

- Riga 26: motore spento, `u = [0; 0]`. La derivata attesa (riga 27) e' caduta
  libera pura: le velocita' si propagano nelle posizioni, `vx_dot = 0`,
  `vy_dot = -1` (solo gravita'), `m_dot = 0`.
- **Perche' e' l'invariante giusto**: l'arco di coast e' un pezzo *reale* della
  soluzione ottima. La traccia produce un profilo bang-off-bang (max-coast-max):
  nella soluzione nominale il motore e' spento fra circa t = 14.0 s e t = 33.1 s
  (vedi il README di HM2). Con `tf = 38 s` sono 19.1 s di coast: il solver passera'
  **circa meta'** del tempo di volo a valutare la RHS esattamente in `u = 0`. Se la
  RHS fosse mal definita li', l'ottimizzazione lavorerebbe su spazzatura proprio nel
  tratto piu' lungo.
- Il punto delicato che questo test tocca senza dirlo: `u = 0` e' il **punto di
  non differenziabilita'** della norma `||u|| = sqrt(u1^2 + u2^2)`. Il *valore*
  della RHS li' e' perfettamente definito (e il test lo prova: nessun NaN,
  `m_dot = 0`), ma il *gradiente* rispetto a u non lo e': la norma euclidea ha un
  cono in zero. Il test verifica la continuita' del valore, non la
  differenziabilita' -- ed e' bene sapere che quest'ultima *manca*. E' esattamente
  la causa radice dello stallo di convergenza documentato nel README (fmincon con
  gradienti a differenze finite si blocca a ottimalita' del primo ordine 1e-3 /
  1e-4 sull'arco di coast). La convessificazione lossless con la slack Gamma >= ||T||
  esiste proprio per togliere la norma dalla dinamica.
- Nota di onesta': `u = [0;0]` produce `0/x(5) = 0` solo perche' `x(5) = 0.9 > 0`.
  La funzione **non ha alcuna guardia** su massa nulla: `x(5) = 0` darebbe
  `0/0 = NaN` o `Inf`. Nessun test copre quel caso, e la funzione non ha `arguments`
  block. Il contratto e' delegato al call site: in `main_task2.m` i box bounds
  impongono `m >= 1e-3`.

---

## `testHoverEquilibrium` (righe 30-37)

```matlab
m  = 0.73;
x  = [0; 1; 0; 0; m];
dx = ode_descent(x, [0; m], 0.0777);
testCase.verifyEqual(dx(3), 0, 'AbsTol', 1e-15);
testCase.verifyEqual(dx(4), 0, 'AbsTol', 1e-15);
```

- Righe 32-34: spinta puramente verticale di modulo pari alla massa,
  `u = [0; m]`. Allora `vy_dot = m/m - 1 = 0` e `vx_dot = 0/m = 0`.
- **Perche' e' l'invariante giusto**: e' l'unico test che verifica la *scala*
  della gravita', non solo il suo segno. In variabili dimensionali la condizione
  di hover e' `T = m*g`; nel codice adimensionalizzato diventa `Ty = m` perche'
  `g = 1`. Se lo scaling fosse sbagliato (per esempio se qualcuno avesse lasciato
  `- 9.80665` o un `- g0/a_ref` diverso da 1 nella riga 14 di `ode_descent.m`),
  `testDerivativeDefinition` fallirebbe ma in modo poco parlante; qui il fallimento
  dice esattamente "la tua condizione di hover non e' hover", cioe' punta al
  bug fisico.
- C'e' un legame diretto col solver: `init_guess` in `main_task2.m` (riga 673)
  inizializza il controllo proprio a `T = m0` (hover, gravita' = 1 nondim). Questo
  test e' quindi anche la validazione della *guess iniziale* usata da tutte le
  varianti.
- Sottigliezza importante: il test asserisce solo `dx(3)` e `dx(4)`, **non** l'intero
  vettore. Corretto: l'hover NON e' un punto di equilibrio del sistema a 5 stati,
  perche' `m_dot = -Vc*m < 0`. Il veicolo tiene la quota ma consuma propellente,
  quindi la spinta di hover richiesta cambia istante per istante. Asserire
  `dx == 0` sarebbe stato un errore concettuale.

> **Possibile domanda d'esame** -- L'hover e' un equilibrio del sistema?
> *Risposta:* No. E' un equilibrio del *sottosistema traslazionale* (le derivate
> di vx e vy si annullano istantaneamente), ma non del sistema completo: la riga
> di massa resta `m_dot = -Vc*||T|| < 0`. Il sistema non ha punti fissi con motore
> acceso; l'unico equilibrio genuino sarebbe motore spento e veicolo a terra. Per
> questo il test verifica solo le componenti 3 e 4 del vettore derivata.

---

## `testMassFlowDependsOnlyOnThrustMagnitude` (righe 39-47)

```matlab
x   = [0; 1; 0; 0; 1];
dx1 = ode_descent(x, [1; 0],  Vc);
dx2 = ode_descent(x, [0; -1], Vc);
testCase.verifyEqual(dx1(5), -Vc,    'AbsTol', 1e-15);
testCase.verifyEqual(dx2(5), dx1(5), 'AbsTol', 1e-15);
```

- Righe 43-44: due controlli di modulo unitario ma direzione completamente
  diversa: uno orizzontale (`[1; 0]`), l'altro verticale **verso il basso**
  (`[0; -1]`).
- Righe 45-46: entrambi devono dare `m_dot = -Vc`. E' l'**isotropia del consumo**:
  il razzo brucia propellente in proporzione al modulo della spinta, non alla sua
  direzione. La legge e' `m_dot = -||T||/c` che, adimensionalizzata, diventa
  `m_dot = -Vc*||T||`.
- **Perche' e' l'invariante giusto**: e' un test di *falsificazione mirata*. La
  scelta di `[0; -1]` e' astuta: se qualcuno implementasse la riga di massa con una
  componente con segno invece che con la norma (per esempio `-Vc*u(2)`, errore
  plausibile in un modello 1D verticale poi esteso a 2D), quel controllo darebbe
  `m_dot = +Vc`, cioe' massa **crescente**: il razzo si riempirebbe di propellente
  spingendo verso il basso, e il solver di minimo consumo sfrutterebbe
  immediatamente quel bug per "guadagnare" massa finale. La funzione obiettivo di
  HM2 e' `max m(tf)`: qualsiasi errore sulla riga 5 e' un errore *sull'obiettivo*,
  non solo sulla dinamica. Da qui la dignita' di test dedicato.
- Nota fisica: `u = [0; -1]` e' assurdo per un lander (spinta verso il suolo), ma
  la RHS non deve saperlo -- i vincoli fisici stanno nei bound e nei path
  constraint, non nel modello. Il test sfrutta correttamente questa separazione di
  responsabilita'.

---

## Cosa questa suite NON copre (onesta' sulla copertura)

- **Nessun test sulla singolarita' di massa.** `x(5) -> 0` fa esplodere `u/x(5)`.
  La protezione vive nei box bounds del NLP (`m >= 1e-3`), non nella funzione: e'
  una scelta deliberata (nessun `arguments` block nelle hot loop), ma va detta.
- **Nessun test sulle Jacobiane analitiche.** `jacobians` in `main_task2.m` (righe
  255-277) calcola A e B a mano, con una regolarizzazione `Tmag_reg` per evitare la
  divisione per zero della derivata della norma. Quella si' che meriterebbe un
  confronto con differenze finite -- ed e' assente.
- **Nessun test sulla forma degli input.** La funzione non ha `arguments` block; il
  contratto (colonne 5x1 e 2x1) e' solo nel commento di header.
- **Nessun test di conservazione lungo una traiettoria.** Qui si testa solo la RHS
  puntuale; la monotonia di m(t) e la positivita' della quota sono verificate
  altrove (o non verificate affatto).

---

## Possibili domande d'esame

**D: Quali invarianti verifica questa suite, e perche' proprio quelli?**
R: Quattro. (1) La *definizione* della derivata in un punto generico, contro un
valore calcolato a mano -- inchioda ordine delle righe, divisione per la massa,
segno della gravita'. (2) Il *coast balistico* (u = 0), perche' e' un arco reale
della soluzione bang-off-bang ed e' anche il punto di non differenziabilita' della
norma. (3) L'*hover* (Ty = m), che e' l'unico test che verifica la scala
dell'adimensionalizzazione (g = 1). (4) L'*isotropia del consumo*
(m_dot = -Vc*||T||), che protegge la funzione obiettivo `max m(tf)` da un errore di
segno. Insieme coprono le cinque righe della RHS con test che falliscono in modo
diagnostico, cioe' dicono *quale* pezzo di fisica e' rotto.

**D: Perche' la tolleranza e' 1e-15 e non zero?**
R: Perche' `0.6^2 = 0.36` non e' rappresentabile esattamente in virgola mobile
binaria, quindi `sqrt(0.6^2 + 0.8^2)` puo' differire da 1 di circa un ulp
(2.2e-16). Con `verifyEqual` esatto il test sarebbe fragile rispetto a dettagli di
arrotondamento della libreria matematica. 1e-15 e' appena sopra la macchina
epsilon: resta un test di identita' esatta, non un test di "vicinanza numerica".

**D: L'RHS divide per la massa. Cosa succede se la massa va a zero, e perche' non
c'e' un controllo?**
R: Si otterrebbe Inf o NaN e il solver divergerebbe. Non c'e' controllo perche'
`ode_descent` e' una hot loop: viene chiamata dell'ordine di 1e8 volte dentro
fmincon (49 intervalli x 4 valutazioni RK4 x n_sub substep x ~350 valutazioni per
gradiente a differenze finite x fino a 1000 iterazioni). Un `arguments` block
costerebbe circa un microsecondo per chiamata contro i ~60 ns della funzione
stessa: renderebbe la validazione piu' costosa della fisica. La protezione e'
spostata al call site (box bound `m >= 1e-3` nel NLP). E' la stessa filosofia
documentata nell'header del file e nel CLAUDE.md della repo.

**D: Perche' l'arco di coast e' il punto piu' delicato di tutta HM2?**
R: Perche' e' il punto in cui la dinamica e' perfettamente sana ma il suo
*gradiente* non esiste. `||T||` ha un cono in T = 0; fmincon con differenze finite
vi calcola un gradiente arbitrario, e la convergenza al primo ordine si ferma a
1e-3/1e-4 (documentato nel README). Il test `testBallisticCoast` prova che il
valore e' corretto, non che il problema sia liscio. La cura strutturale e' la
convessificazione lossless: si introduce la slack Gamma >= ||T|| e si sostituisce
la norma nella dinamica con Gamma, che e' lineare -- e si dimostra che
all'ottimo il vincolo e' attivo (Gamma = ||T||), quindi non si perde nulla.

**D: Perche' esistono due RHS diverse (`ode_descent` e `ode_descent_uacc`) e questa
suite ne testa una sola?**
R: Perche' realizzano due convenzioni ZOH diverse. `ode_descent` tiene costante la
**spinta T** sull'intervallo (convenzione della trascrizione nonlineare + RK4);
`ode_descent_uacc` tiene costante l'**accelerazione u = T/m** (convenzione nativa
della trascrizione GFOLD log-massa, dove il cambio di variabile z = ln m rende la
dinamica esattamente LTI). Sotto ZOH le due producono traiettorie *diverse* sullo
stesso intervallo (coincidono solo nel limite dt -> 0). `odeDescentTest` copre la
prima; la seconda, insieme a `lti_zoh.m`, e' coperta da `gfoldLogMassTest.m`.
