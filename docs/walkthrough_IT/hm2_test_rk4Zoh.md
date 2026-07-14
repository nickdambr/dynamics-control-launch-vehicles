# HM2_powered_descent/tests/rk4ZohTest.m

## Ruolo del file nel progetto

Questa suite (`matlab.unittest`, class-based) testa `rk4_zoh.m`, il propagatore
Runge-Kutta a 4 stadi che HM2 usa per la trascrizione ZOH (Task 2, variante a:
"Nonlinear ZOH + RK4", multiple shooting). E' il pezzo che trasforma la dinamica
continua in una mappa discreta:

    x_{k+1} = RK4(x_k, u_k, dt)

e i vincoli di uguaglianza del NLP (i *defect*) sono esattamente
`Z(1:5,k+1) - rk4_zoh(...)` (vedi `zoh_nonlcon` in `main_task2.m`, righe 797-800).
Quindi la mappa discreta **e' il modello** su cui fmincon ottimizza: se il
propagatore ha un bug, il solver converge alla soluzione esatta di una dinamica
che non esiste, e la replay open-loop con `ode45` la smaschera con un errore di
atterraggio. La suite serve a evitare che si arrivi a quel punto.

`rk4_zoh.m` (righe 15-23) e' il classico RK4 con `n_sub` substep interni:

    h = dt / n_sub
    k1 = f(x),  k2 = f(x + h/2*k1),  k3 = f(x + h/2*k2),  k4 = f(x + h*k3)
    x  = x + (h/6)*(k1 + 2*k2 + 2*k3 + k4)

con `f = ode_descent` e il controllo `u` **tenuto costante** su tutto
l'intervallo -- questo e' il significato di *zero-order hold*. In produzione
`main_task2.m` usa `n_sub = 2` (riga 37) e `dt = tf_nd/(N-1)` con N = 50.

Nota di convenzione, cruciale e ricorrente: qui e' la **spinta T** a essere tenuta
costante, non l'accelerazione. La variante GFOLD (`lti_zoh.m` + `ode_descent_uacc.m`)
tiene costante `u = T/m`. Sono due ZOH *diversi*: producono traiettorie diverse
sullo stesso intervallo finito e coincidono solo nel limite dt -> 0. Ci torniamo
alla fine, perche' e' il motivo per cui in questa suite **non** esiste un test
"rk4_zoh riproduce lti_zoh".

I tre test coprono, nell'ordine: accuratezza contro un riferimento fidato,
esattezza strutturale sulla riga di massa, e -- il piu' importante -- **ordine di
convergenza**.

---

## Proprieta' costanti (righe 5-10)

```matlab
properties (Constant)
    Vc = 0.0777                    % V_ref/c
    x0 = [0.3; 1; -0.2; -0.6; 1]   % non-dim initial state
    u  = [0.5; 1.0]                % non-dim ZOH control
    dt = 0.4                       % non-dim ZOH interval
end
```

- Riga 6: `Vc = 0.0777` e' il valore di produzione (numero di Tsiolkovsky,
  `V_ref/c`, dai dati di Tabella 1 della traccia).
- Riga 7: stato iniziale con **tutte le componenti non banali** (posizione fuori
  asse, velocita' orizzontale e verticale entrambe non nulle, massa 1). Serve per
  eccitare tutte e cinque le righe della RHS: un test con `vx = 0` non
  distinguerebbe un bug che azzera la prima riga.
- Riga 8: `u = [0.5; 1.0]`, quindi `||u|| = 1.118` e con `m = 1` l'accelerazione
  supera la gravita' (il veicolo decelera). Anche qui: entrambe le componenti non
  nulle, cosi' la norma nella riga di massa e' "viva".
- Riga 9: **`dt = 0.4`**. Attenzione: NON e' il dt di produzione, che vale
  `tf_nd/(N-1) = 0.0444` (lo si legge, quel valore, in `descentDynamicsPerformanceTest.m`
  e in `gfoldLogMassTest.m`). E' circa 9 volte piu' grande, ed e' una scelta
  deliberata: serve a mantenere l'errore di troncamento di RK4 **ben sopra il
  rumore del riferimento** per tutti gli `n_sub` testati. Se si usasse dt = 0.0444,
  con `n_sub = 8` l'errore RK4 (che scala come h^4 = (dt/n_sub)^4) scenderebbe sotto
  la soglia di accuratezza di ode45 e il test di ordine misurerebbe il rumore del
  riferimento, non l'integratore. E' il principio base di ogni convergence study:
  bisogna restare nel **regime asintotico** ma **sopra il floor numerico**.

> **Possibile domanda d'esame** -- Perche' il test usa dt = 0.4 invece del dt di
> produzione (0.0444)?
> *Risposta:* Perche' il test di ordine ha bisogno che l'errore sia dominato dal
> troncamento di RK4, non dagli errori del riferimento. Con dt piccolo e n_sub = 8
> il passo effettivo h = dt/n_sub e' cosi' piccolo che l'errore RK4 (~h^4) cade sotto
> la tolleranza di ode45 (1e-12 relativa) e il rapporto err(n)/err(2n) perde
> significato: si misurerebbe rumore. Ingrandendo dt si tiene la successione degli
> errori in una finestra dove il termine C*h^4 domina, e la stima log2 del rapporto
> restituisce davvero 4. Il test di accuratezza e il test di ordine hanno esigenze
> opposte sul passo, e questo file risolve il conflitto scegliendo il passo per il
> test piu' esigente.

---

## `TestClassSetup` -- `addHm2ToPath` (righe 12-17)

Identico a quello di `odeDescentTest.m`: due `fileparts` su `mfilename('fullpath')`
risalgono da `tests/` alla cartella `HM2_powered_descent/`, e una `PathFixture`
la aggiunge al path con teardown automatico a fine classe. Rende la suite
lanciabile da qualunque working directory e non lascia residui nel path
dell'utente.

---

## `ode45Reference` -- metodo statico (righe 44-54)

```matlab
opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
[~, Y] = ode45(@(t, x) ode_descent(x, rk4ZohTest.u, rk4ZohTest.Vc), ...
               [0 rk4ZohTest.dt], rk4ZohTest.x0, opts);
xRef = Y(end, :).';
```

- Riga 45: e' un metodo `Static`, quindi puo' essere invocato come
  `rk4ZohTest.ode45Reference()` senza istanza; per la stessa ragione dentro il
  corpo le costanti si leggono col nome della classe (`rk4ZohTest.u`, riga 50) e
  non con `testCase.u`.
- Riga 49: tolleranze **strette** (RelTol 1e-12, AbsTol 1e-14). Non e' pignoleria:
  il riferimento deve essere piu' accurato di cio' che si sta misurando di almeno
  qualche ordine di grandezza, altrimenti il test misura il riferimento.
- Righe 50-51: la closure `@(t, x) ode_descent(x, u, Vc)` **ignora t** e passa
  sempre lo stesso `u`. Questa e' la definizione operativa di ZOH: il controllo e'
  costante sull'intervallo, quindi la traiettoria "vera" contro cui confrontarsi
  e' il flusso della ODE autonoma con quel controllo congelato. Se il riferimento
  interpolasse il controllo (come fa la trascrizione trapezoidale, che e' PWL), il
  confronto sarebbe con un'altra dinamica e il test fallirebbe pur essendo RK4
  corretto.
- Riga 52: si prende solo l'ultimo campione, `Y(end,:).'`, cioe' lo stato a t = dt:
  quello che RK4 deve riprodurre.
- Nota: `ode45Reference` viene ricalcolata da capo in ogni test che ne ha bisogno
  (righe 21 e 36), quindi due volte per run. E' una decina di millisecondi, del
  tutto irrilevante -- ma tecnicamente sarebbe cachabile in una proprieta' di
  classe.
- **Onesta' sul riferimento**: per questa dinamica una soluzione analitica *esiste*.
  Con T costante, `m(t) = m0 - Vc*||T||*t` e' lineare, quindi
  `vx_dot = Tx/m(t)` si integra in forma chiusa (`vx = vx0 - (Tx/(Vc*||T||)) * ln(m(t)/m0)`,
  una Tsiolkovsky per componente), e le posizioni seguono da `integral(ln)`, anch'esso
  elementare. Un riferimento analitico eliminerebbe del tutto il floor a 1e-12.
  Il codice ha scelto ode45: piu' semplice da scrivere e da leggere, ma pone un
  limite inferiore artificiale all'accuratezza verificabile.

---

## `testMatchesOde45Reference` (righe 20-24)

```matlab
xRef = rk4ZohTest.ode45Reference();
xRk4 = rk4_zoh(testCase.x0, testCase.u, testCase.dt, testCase.Vc, 8);
testCase.verifyEqual(xRk4, xRef, 'AbsTol', 1e-8);
```

- Riga 22: RK4 con **8 substep** su dt = 0.4, cioe' h = 0.05.
- Riga 23: tolleranza 1e-8. E' una soglia che vive nel mezzo fra due scale:
  l'errore del riferimento (ordine 1e-12, fissato dalle tolleranze di ode45) e
  l'errore di troncamento globale di RK4 con h = 0.05 (empiricamente sotto 1e-8).
  Passare significa: "l'RK4 a 8 substep e' accurato almeno a 1e-8, e la differenza
  osservata e' dominata dall'RK4, non dal riferimento".
- **Perche' e' l'invariante giusto**: e' il test di *accuratezza* (il propagatore
  approssima davvero il flusso della ODE), ed e' la prima linea di difesa contro i
  bug grossolani -- un `h` invece di `h/2` in `k2`, un pesi sbagliato, un segno.
  Ma da solo NON basta: un integratore rotto in modo sottile (per esempio con i
  pesi `(1,1,1,1)/4` invece di `(1,2,2,1)/6`) puo' comunque restare entro 1e-8 su
  un intervallo cosi' corto, semplicemente perche' converge -- solo piu'
  lentamente. Ecco perche' esiste il terzo test.

---

## `testMassRowIsExact` (righe 26-32)

```matlab
xRk4 = rk4_zoh(testCase.x0, testCase.u, testCase.dt, testCase.Vc, 1);
mExpected = testCase.x0(5) ...
    - testCase.Vc*norm(testCase.u)*testCase.dt;
testCase.verifyEqual(xRk4(5), mExpected, 'AbsTol', 1e-14);
```

- Riga 29: **un solo substep** (`n_sub = 1`), cioe' il caso piu' brutale possibile:
  un unico passo RK4 su tutto dt = 0.4.
- Riga 30: il valore atteso e' la retta `m(dt) = m0 - Vc*||u||*dt`.
- Riga 31: tolleranza 1e-14, cioe' *esattezza a meno di arrotondamenti*, e non un
  errore di troncamento.

**La matematica dietro.** Sotto ZOH sulla spinta, la riga di massa e'

    m_dot = -Vc * ||u||

con u costante: il secondo membro **non dipende dallo stato**, e' una costante.
Applicando RK4 a `x_dot = a` (a costante) si ha `k1 = k2 = k3 = k4 = a`, quindi

    x + (h/6)*(a + 2a + 2a + a) = x + (h/6)*6a = x + h*a

cioe' RK4 degenera nel metodo di Eulero, che per una derivata costante e' la
soluzione **esatta**. Non c'e' errore di troncamento sulla quinta riga, per
qualunque `n_sub`, incluso `n_sub = 1`. Il test verifica proprio questo.

**Perche' e' l'invariante giusto.** Perche' la massa finale *e' la funzione
obiettivo*: `fmincon` massimizza `m(tf)` (in `solve_zoh`, riga 646, si minimizza
`-z(iN_m)`). Se il propagatore introducesse anche solo un errore di troncamento
sulla riga di massa, l'ottimizzatore lo sfrutterebbe: guadagnerebbe massa finale
"fittizia" scegliendo profili di spinta che ingannano il discretizzatore invece di
risparmiare propellente. Un test dedicato che pretende esattezza *strutturale*
(non "vicinanza") sulla riga della funzione obiettivo e' quindi molto piu' di una
verifica di accuratezza: e' una garanzia che il modello discreto non regali
carburante.

**Il contrappunto che spiega tutta HM2.** Questa esattezza vale SOLO con la
convenzione T-ZOH. Con la convenzione accelerazione-ZOH (`ode_descent_uacc.m`) la
riga di massa diventa `m_dot = -Vc * m * ||u||`: dipende dallo stato, la soluzione
e' un esponenziale, e RK4 **non** sarebbe piu' esatto. La mossa di GFOLD e' proprio
questa: passare a `z = ln m`, per cui `z_dot = -Vc*||u||` torna costante -- e infatti
`lti_zoh.m` discretizza il sistema con un singolo `expm`, esattamente. Le tre
varianti di HM2 sono tre modi diversi di preservare l'esattezza del canale di massa.

> **Possibile domanda d'esame** -- Perche' testare la riga di massa separatamente,
> se il test contro ode45 gia' confronta tutte e cinque le componenti?
> *Risposta:* Perche' il test contro ode45 verifica una *approssimazione* (entro
> 1e-8), mentre qui si verifica una *identita' esatta* (entro 1e-14, cioe'
> arrotondamento macchina) valida per qualunque n_sub, incluso n_sub = 1. Sono
> affermazioni di forza diversa: la seconda e' una proprieta' strutturale del
> metodo applicato a quella specifica dinamica, e proteggerla e' importante perche'
> quella riga e' la funzione obiettivo del NLP. Un errore di troncamento sulla massa
> si tradurrebbe in propellente immaginario, non in un piccolo errore numerico.

---

## `testFourthOrderConvergence` (righe 34-41)

```matlab
xRef = rk4ZohTest.ode45Reference();
err = arrayfun(@(n) norm(rk4_zoh(testCase.x0, testCase.u, ...
    testCase.dt, testCase.Vc, n) - xRef), [1 2 4 8]);
order = log2(err(1:end-1) ./ err(2:end));
testCase.verifyGreaterThan(min(order), 3.5);
```

- Righe 37-38: si calcola l'errore (norma euclidea sul vettore di stato finale)
  per `n_sub = 1, 2, 4, 8`. Raddoppiare `n_sub` significa **dimezzare il passo**
  `h = dt/n_sub`.
- Riga 39: `order = log2(err(k)/err(k+1))` per le tre coppie consecutive.
- Riga 40: si pretende che il **minimo** dei tre esponenti stimati superi 3.5.

**La matematica dietro il test dell'ordine.** Un metodo a un passo di ordine p ha
errore globale (sull'intervallo fissato dt) della forma

    e(h) = C * h^p + O(h^(p+1))

con C costante che dipende dalle derivate della soluzione ma **non** da h.
Dimezzando il passo:

    e(h/2) = C * (h/2)^p = C * h^p / 2^p = e(h) / 2^p

quindi il rapporto fra errori consecutivi e'

    e(h) / e(h/2) = 2^p        =>        p = log2( e(h) / e(h/2) )

Per RK4 classico p = 4, quindi il rapporto atteso e' 2^4 = **16**: dimezzando il
passo, l'errore deve calare di un fattore 16. E' esattamente cio' che la riga 39
misura -- il logaritmo in base 2 del rapporto restituisce direttamente l'ordine
osservato.

**Perche' e' l'invariante giusto -- anzi, il piu' importante dei tre.** Un test di
accuratezza a passo singolo ("l'errore e' sotto 1e-8") verifica un *numero*. Il
test di ordine verifica l'*algoritmo*. La differenza e' sostanziale:

- Se si sbagliano i pesi della quadratura (per esempio `(k1+k2+k3+k4)/4` invece di
  `(k1+2k2+2k3+k4)/6`), il metodo resta **consistente** e converge -- ma cade a
  ordine 2. Su un intervallo corto con h piccolo l'errore assoluto puo' ancora
  passare sotto 1e-8 e il test di accuratezza **passerebbe lo stesso**. Il test di
  ordine no: misurerebbe un esponente vicino a 2 e fallirebbe.
- Analogamente, valutare `k4` in `x + 0.5*h*k3` invece che in `x + h*k3` (un errore
  di battitura del tutto plausibile, righe 19-20 del sorgente sono quasi identiche)
  degrada l'ordine senza rompere nulla di visibile.

In altre parole: l'ordine di convergenza e' la **firma** di un integratore. Se la
firma e' 4, allora i coefficienti di Butcher sono giusti, punto. E' il test che
distingue "il codice da' un numero plausibile" da "il codice implementa RK4".

**Perche' la soglia e' 3.5 e non 4.** Ci sono quattro fonti di degrado che rendono
l'esponente misurato leggermente diverso da 4, tutte legittime:

1. il riferimento non e' esatto: `xRef` ha il suo errore (~1e-12), che contamina
   soprattutto `err(4)`, il piu' piccolo -- e quindi *abbassa* l'ultimo rapporto;
2. a `n_sub = 1` con h = 0.4 non si e' pienamente nel regime asintotico: i termini
   `O(h^5)` non sono ancora trascurabili, quindi il primo rapporto puo' scostarsi;
3. la norma mescola cinque canali con costanti C diverse (e uno, la massa, con
   errore identicamente zero -- vedi il test precedente);
4. rumore di arrotondamento in virgola mobile.

3.5 e' una soglia che accetta questo margine ma **rifiuta ordine 3** (e a maggior
ragione ordine 2 o 1): e' abbastanza lasca da non essere fragile, abbastanza
stretta da essere diagnostica.

**Un dettaglio spesso frainteso**: `verifyGreaterThan(min(order), 3.5)` chiede che
**tutti e tre** i rapporti superino 3.5, non che la media lo faccia. E' una
richiesta piu' forte: non basta che l'errore cali complessivamente, deve calare col
tasso giusto a *ogni* raffinamento.

> **Possibile domanda d'esame** -- Se dimezzo il passo e l'errore cala di 16, che
> ordine ha il metodo? E se calasse di 4?
> *Risposta:* 16 = 2^4, quindi ordine 4 (RK4 corretto). Un fattore 4 = 2^2
> significherebbe ordine 2, cioe' un metodo consistente ma con i coefficienti
> sbagliati -- tipicamente pesi errati nella combinazione finale o uno stadio
> valutato nel punto sbagliato. Un fattore 2 sarebbe Eulero esplicito. La relazione
> generale e' e(h)/e(h/2) = 2^p, da cui p = log2 del rapporto, ed e' esattamente
> quello che calcola la riga 39.

---

## Cosa questa suite NON copre (onesta' sulla copertura)

- **Nessun confronto diretto fra `rk4_zoh` e `lti_zoh`.** Il brief lo cita come
  invariante desiderabile ("le due discretizzazioni devono essere consistenti"), ma
  **non esiste** in questo file -- e va detto perche' la versione ingenua di quel
  test sarebbe *sbagliata*. `rk4_zoh` tiene costante la **spinta T**; `lti_zoh`
  discretizza il sistema log-massa in cui e' costante l'**accelerazione u = T/m**.
  Sono ZOH diversi: dentro l'intervallo, con T costante l'accelerazione *cresce*
  (la massa cala), mentre con u costante e' la spinta a *calare*. Le due traiettorie
  divergono a O(dt^2) e coincidono solo nel limite dt -> 0. Un test che pretendesse
  `rk4_zoh(x,u,dt) == [Abar Bbar cbar]*[...]` fallirebbe correttamente. La
  consistenza di `lti_zoh` e' invece verificata dove ha senso, in
  `gfoldLogMassTest.m`: contro la forma chiusa analitica delle matrici discrete
  (`testZohClosedForm`), contro un `ode45` della *stessa* dinamica LTI
  (`testZohMatchesOde45`), e contro la replay nonlineare del canale di massa con
  `ode_descent_uacc` (`testMassRowConsistency`).
- **Nessun test su n_sub non intero o nullo.** `rk4_zoh` non ha `arguments` block
  (scelta deliberata: hot loop). `n_sub = 0` darebbe `h = Inf` e il ciclo non
  girerebbe, restituendo `x` invariato: un fallimento silenzioso. La validazione e'
  al call site (`solve_zoh` in `main_task2.m`, riga 639, valida `n_sub` con
  `mustBeInteger, mustBePositive`).
- **Nessun test di stabilita'.** Non si verifica che dt stia dentro la regione di
  stabilita' assoluta di RK4 per questa dinamica. Per il problema in esame non e'
  un rischio (la dinamica non e' stiff), ma non e' testato.
- **Un solo punto di lavoro.** Tutti i test usano lo stesso `x0`/`u`. Un test
  parametrizzato (o `u = 0`, il caso di coast, dove la norma ha il cono) darebbe
  copertura piu' larga.

---

## Possibili domande d'esame

**D: Quali invarianti verifica questa suite, e in che ordine di forza?**
R: Tre, di forza crescente. (1) *Accuratezza*: RK4 a 8 substep riproduce una
integrazione ode45 a tolleranze strette entro 1e-8 -- prova che il propagatore
approssima il flusso della ODE giusta (con u congelato, cioe' proprio la semantica
ZOH). (2) *Esattezza strutturale sulla massa*: sotto T-ZOH `m_dot = -Vc*||u||` e'
costante, quindi RK4 e' esatto sulla quinta riga per qualunque n_sub -- e questo
protegge la funzione obiettivo `max m(tf)`. (3) *Ordine di convergenza*: dimezzando
il passo l'errore deve calare di 2^4 = 16. Il terzo e' il piu' forte perche' verifica
l'algoritmo (i coefficienti di Butcher) e non solo il numero prodotto.

**D: Spiega la matematica dietro il test di ordine.**
R: L'errore globale di un metodo a un passo di ordine p sull'intervallo fissato dt
si scrive e(h) = C*h^p + O(h^(p+1)), con C indipendente da h. Dimezzando h:
e(h/2) = C*(h/2)^p = e(h)/2^p, quindi log2(e(h)/e(h/2)) = p. Il test valuta l'errore
per n_sub = 1, 2, 4, 8 (h = dt/n_sub, quindi ogni passaggio dimezza h), calcola i tre
rapporti consecutivi e ne prende il log2: sono tre stime dell'ordine. Si pretende che
il minimo sia > 3.5, cioe' che *ogni* raffinamento mostri il tasso quartico. La
soglia non e' 4 esatto perche' il riferimento ode45 ha il suo errore (~1e-12) che
contamina l'errore piu' piccolo, perche' a h = 0.4 i termini di ordine superiore non
sono ancora trascurabili, e per il rumore di arrotondamento.

**D: Perche' RK4 e' esatto sulla riga di massa? Vale anche per la variante GFOLD?**
R: Perche' sotto ZOH sulla *spinta* la derivata della massa e' -Vc*||u||, una
costante che non dipende dallo stato. Applicando RK4 a x_dot = a costante si ottiene
k1 = k2 = k3 = k4 = a e l'aggiornamento diventa x + (h/6)*6a = x + h*a, che e' la
soluzione esatta. Nella variante GFOLD la convenzione e' diversa: si tiene costante
l'accelerazione u = T/m, quindi m_dot = -Vc*m*||u|| dipende dallo stato e la soluzione
e' esponenziale -- RK4 non sarebbe piu' esatto. La contromossa e' il cambio di
variabile z = ln m, che riporta z_dot = -Vc*||u|| a essere costante: la dinamica
diventa esattamente LTI e `lti_zoh.m` la discretizza con un solo `expm`, senza alcun
errore di troncamento.

**D: Perche' non c'e' un test che confronta `rk4_zoh` con `lti_zoh`?**
R: Perche' non sarebbero d'accordo, e giustamente. Le due funzioni implementano due
convenzioni ZOH *diverse*: `rk4_zoh` congela la spinta T (quindi dentro l'intervallo
l'accelerazione T/m cresce man mano che la massa cala), `lti_zoh` congela
l'accelerazione u = T/m (quindi e' la spinta a calare). Le traiettorie sull'intervallo
differiscono a ordine dt^2 e convergono solo per dt -> 0. La consistenza sensata da
verificare -- ed e' quella che `gfoldLogMassTest.m` verifica -- e' che `lti_zoh` coincida
con l'integrazione ode45 della *sua* dinamica e che il suo canale di massa coincida
con la replay nonlineare di `ode_descent_uacc`.

**D: Cosa fallirebbe se sbagliassi un peso di RK4, e quale test lo prenderebbe?**
R: Il metodo resterebbe consistente (convergerebbe comunque a dt -> 0) ma
crollerebbe di ordine, tipicamente a 2 o 1. Il test contro ode45 potrebbe ancora
passare, perche' su un intervallo corto e con 8 substep anche un metodo di ordine 2
puo' finire sotto 1e-8. Il test di ordine invece misurerebbe log2 del rapporto degli
errori e troverebbe circa 2 invece di 4, fallendo la soglia 3.5. E' precisamente la
ragione per cui in una suite di test di un integratore il convergence test non e'
ridondante rispetto al test di accuratezza: sono i due test *diversi* che servono.
