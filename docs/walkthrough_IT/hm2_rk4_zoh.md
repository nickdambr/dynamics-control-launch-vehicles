# HM2_powered_descent/rk4_zoh.m

## Ruolo del file nel progetto

`rk4_zoh` e' il **propagatore discreto** della variante (a) di Task 2 ("Nonlinear
ZOH + RK4"). Prende lo stato al nodo k, il controllo tenuto costante
sull'intervallo, e restituisce lo stato al nodo k+1 integrando `ode_descent` con
un Runge-Kutta classico a 4 stadi e passo fisso.

Serve a costruire il vincolo di difetto della trascrizione ZOH. Mentre Task 1
usa la regola trapezoidale

    x_{k+1} - x_k - (dt/2)*(f_k + f_{k+1}) = 0

la variante ZOH usa un **multiple-shooting** su un solo intervallo:

    x_{k+1} - RK4(x_k, u_k, dt) = 0

(vedi `zoh_nonlcon`, `main_task2.m` righe 797-800). La differenza concettuale e'
profonda: nel trapezoidale il difetto e' una *approssimazione* dell'integrale
esatto (errore O(dt^2)), mentre nel ZOH il difetto e' una *simulazione* del plant
sotto un controllo che e' esattamente quello che verra' comandato. Il residuo di
errore non e' piu' un errore di modello, ma solo l'errore di integrazione di RK4,
che si puo' spingere a piacere aumentando `n_sub`. Il risultato si vede nel
README di HM2: il replay `ode45` della soluzione trapezoidale atterra a 4.3 m dal
pad, quello della soluzione ZOH+RK4 riproduce i nodi a 1.4e-8 nondim.

Chiamato da: `main_task2.m` riga 798 (`zoh_nonlcon`, dentro `fmincon`),
`tests/rk4ZohTest.m`, `tests/descentDynamicsPerformanceTest.m` riga 57.
Chiama: `ode_descent.m` (quattro volte per sotto-passo).

---

## `rk4_zoh` (righe 1-25)

```matlab
function x_next = rk4_zoh(x, u, dt, Vc, n_sub)
...
h = dt / n_sub;
for ii = 1:n_sub
    k1 = ode_descent(x,            u, Vc);
    k2 = ode_descent(x + 0.5*h*k1, u, Vc);
    k3 = ode_descent(x + 0.5*h*k2, u, Vc);
    k4 = ode_descent(x +     h*k3, u, Vc);
    x  = x + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
end
x_next = x;
```

- **Riga 1**: firma. `x` (5x1), `u` (2x1, la spinta tenuta costante), `dt`
  (lunghezza dell'intervallo ZOH), `Vc`, `n_sub` (numero di sotto-passi RK4
  dentro l'intervallo). Il valore usato in produzione e' `n_sub = 2`
  (`main_task2.m` riga 37).

- **Righe 2-13**: docstring. Riga 12-13: *"No arguments validation by design:
  hot-loop propagator inside fmincon; validate at the call site."*

- **Riga 15**: `h = dt / n_sub`. Il passo di integrazione **non** e' l'intervallo
  ZOH: l'intervallo viene suddiviso in `n_sub` sotto-passi. Il controllo pero'
  resta lo stesso in tutti i sotto-passi -- il ZOH e' sull'intervallo `dt`, il
  raffinamento e' solo numerico. Con tf_nd = 2.173 e N = 50 si ha dt = 0.0443,
  quindi h = 0.0222 con n_sub = 2.

- **Righe 16-22**: il loop di sotto-passi. Notare che `x` viene **riassegnato in
  place** (riga 21): la variabile d'ingresso viene sovrascritta a ogni giro. In
  MATLAB e' innocuo (semantica per valore, copy-on-write), ed evita di allocare
  un secondo buffer nel hot loop.

- **Righe 17-20**: i quattro stadi. **Il punto centrale del file e' che `u` viene
  passato identico a tutti e quattro.** In un RK4 non autonomo standard gli stadi
  andrebbero valutati a tempi diversi (t_k, t_k+h/2, t_k+h/2, t_k+h) e il
  controllo andrebbe interpolato a quei tempi. Qui il controllo e' costante
  sull'intervallo, quindi non c'e' niente da interpolare e la dinamica e'
  **autonoma**: e' esattamente per questo che `ode_descent` non ha l'argomento
  `t`. Lo ZOH non e' un dettaglio dell'integratore, e' cio' che rende l'integrando
  autonomo.

- **Riga 21**: la combinazione con i pesi di Butcher classici (1, 2, 2, 1)/6.

- **Riga 23**: `x_next = x`, lo stato a fine intervallo.

### Derivazione dei quattro stadi

Per un problema autonomo x_dot = f(x), il RK4 classico (tableau di Kutta 1901)
e':

    k1 = f(x_k)
    k2 = f(x_k + (h/2)*k1)
    k3 = f(x_k + (h/2)*k2)
    k4 = f(x_k + h*k3)
    x_{k+1} = x_k + (h/6)*(k1 + 2*k2 + 2*k3 + k4)

L'idea e' quella di una **quadratura di Simpson sulla derivata**: se f dipendesse
solo dal tempo, l'aggiornamento si ridurrebbe a
x + (h/6)*(f(t) + 4*f(t+h/2) + f(t+h)), che e' esattamente Simpson -- k2 e k3 sono
due stime distinte della derivata a meta' intervallo, mediate con peso 2
ciascuna. Il fatto che k2 usi k1 e k3 usi k2 (invece di ripartire da x_k) e' cio'
che rende il metodo capace di catturare anche la dipendenza da x, portando la
consistenza al quarto ordine.

**Ordine di accuratezza**: RK4 soddisfa le condizioni di ordine fino al quarto,
quindi l'errore locale di troncamento e' O(h^5) e l'errore **globale**
sull'intervallo e' O(h^4). Con n_sub sotto-passi l'errore sull'intero intervallo
ZOH scala come (dt/n_sub)^4: raddoppiare n_sub divide l'errore per 16. Il test
`testFourthOrderConvergence` (`rk4ZohTest.m` riga 34) misura l'ordine empirico
confrontando n_sub = 1, 2, 4, 8 contro un riferimento `ode45` a RelTol 1e-12 e
verifica che l'ordine osservato sia > 3.5.

### La riga della massa e' integrata esattamente

Osservazione non banale, ed e' un test a se' (`testMassRowIsExact`, riga 26). Con
`u` costante, la quinta componente di `ode_descent` vale

    m_dot = -Vc * ||u||   =   costante sull'intervallo

Applicare RK4 a un'equazione x_dot = const da' k1 = k2 = k3 = k4 = const, e
l'aggiornamento diventa x + (h/6)*(1+2+2+1)*const = x + h*const: **esatto**, per
qualunque numero di sotto-passi, anche n_sub = 1. Quindi tutto l'errore di RK4 in
questo problema vive nelle righe di posizione e velocita', mai nella massa.
Il test lo verifica a 1e-14: `m_end = m0 - Vc*||u||*dt`.

Questo e' un altro regalo del ZOH: **se il controllo fosse piecewise-lineare**
(FOH), la spinta T(t) sarebbe affine in t ma il suo **modulo** ||T(t)|| non lo
sarebbe (la norma di una funzione affine e' una iperbole, non una retta), e la
riga della massa introdurrebbe errore di integrazione. Con lo ZOH questo non
succede mai.

> **Possibile domanda d'esame** -- Perche' usare un RK4 a passo fisso dentro il
> vincolo, invece di chiamare direttamente `ode45` che e' piu' accurato?
> *Risposta:* Perche' `fmincon` calcola i gradienti dei vincoli alle differenze
> finite, e ha bisogno che `nonlcon` sia una funzione **liscia e deterministica**
> delle variabili di decisione. Un integratore adattivo cambia il numero e la
> lunghezza dei passi in modo discontinuo quando gli input variano di poco: il
> vincolo risultante ha un "rumore numerico" dell'ordine della tolleranza, che
> distrugge la stima alle differenze finite (il rapporto incrementale con
> delta ~ 1e-7 misura il rumore, non la derivata). RK4 a passo fisso e' invece una
> composizione finita di operazioni lisce: la mappa (x_k, u_k) -> x_{k+1} e'
> C-infinito (a parte il kink di ||u|| in u = 0) e i gradienti FD sono affidabili.
> In piu' e' molto piu' veloce: `ode45` a RelTol 1e-10 su un intervallo costa
> decine di valutazioni del RHS, RK4 con n_sub = 2 ne costa 8.

---

## Perche' ZOH e non un'interpolazione di ordine superiore

Cinque ragioni, in ordine di forza:

1. **Realismo di implementazione.** Un computer di bordo aggiorna il comando di
   spinta a istanti discreti e lo *tiene* fino all'aggiornamento successivo. Lo
   ZOH e' letteralmente cio' che l'attuatore fara'. Una soluzione trascritta con
   controllo piecewise-lineare non e' direttamente eseguibile: bisognerebbe
   ricampionarla, reintroducendo un errore. Il "replay gap" del trapezoidale
   (4.3 m dal pad, README) e' esattamente questo fenomeno.

2. **La dinamica diventa autonoma.** Con u costante non serve interpolare il
   controllo agli stadi intermedi di RK4, e `ode_descent` non ha nemmeno bisogno
   dell'argomento tempo.

3. **La discretizzazione lineare ha forma chiusa.** Per la variante LTV
   (Appendice A) le matrici discrete si ottengono da
   B_k = integrale di Phi(t_{k+1}, s)*B(s) ds su un solo hold. Con un
   **first-order hold** servirebbero **due** matrici di ingresso per intervallo
   (una pesata sul valore a inizio intervallo, una su quello a fine intervallo),
   raddoppiando il sistema aumentato di `ltv_aug_rhs` e la dimensione delle
   matrici da propagare. Il costo del guadagno di accuratezza non e' gratis.

4. **La riga della massa resta esatta** (vedi sopra): con FOH la norma della
   spinta non sarebbe piu' costante e il canale di massa introdurrebbe errore.

5. **La soluzione ottima e' bang-off-bang.** Il profilo di spinta ottimo e'
   max-coast-max (README: switch a t ~ 14.0 s e t ~ 33.1 s): a meno degli istanti
   di switch e' **gia' piecewise-constant in modulo**. Un'interpolazione di ordine
   superiore non aggiungerebbe niente di utile, mentre soffrirebbe di ringing
   attorno alle discontinuita'.

**Il prezzo da pagare (onesta').** Lo ZOH e' una rappresentazione di ordine 0 del
controllo: rispetto all'ottimo *continuo* la soluzione ZOH e' sub-ottima di O(dt),
e i tempi di switch sono risolti solo con la granularita' della griglia. La
trascrizione riproduce con altissima fedelta' (1.4e-8) **il proprio** plant
discreto -- cioe' e' *consistente*, non necessariamente piu' *ottima*. Infatti le
masse finali delle quattro varianti in HM2 sono molto simili: nessuna e' "piu'
ottima", cambiano il costo computazionale e la fedelta' del replay.

---

## Costo computazionale e assenza del blocco `arguments`

Righe 12-13 dichiarano l'assenza della validazione. Il conto e' brutale:

- Una valutazione di `zoh_nonlcon` = (N-1) = 49 chiamate a `rk4_zoh`, ciascuna
  con n_sub = 2 sotto-passi da 4 stadi = **392 valutazioni di `ode_descent`**.
- `fmincon` in SQP con gradienti FD su 7N = 350 variabili ricostruisce lo
  Jacobiano dei vincoli con ~351 valutazioni di `nonlcon` per iterazione maggiore
  = **~1.4e5 chiamate a `ode_descent` per iterazione**.
- Con `MaxIterations = 1000` (`fmincon_opts`, riga 733) si arriva nell'ordine di
  1e8 valutazioni del RHS per una singola risoluzione.

Un blocco `arguments` con `mustBeFinite` / controlli di dimensione aggiungerebbe
un overhead confrontabile con il costo dell'aritmetica stessa (il benchmark stima
`ode_descent` intorno ai 60 ns). Da qui la scelta, dichiarata anche in
`CLAUDE.md`. La validazione sta nelle funzioni di frontiera: `solve_zoh`,
`nondim`, `fwd_integrate`.

> **Possibile domanda d'esame** -- Come si sceglie `n_sub`?
> *Risposta:* E' un compromesso fra fedelta' e costo, entrambi lineari-ish in
> n_sub (il costo cresce come n_sub, l'errore cala come n_sub^-4). Con n_sub = 2 e
> h = 0.022 nondim l'errore di nodo osservato e' 1.4e-8, gia' molto sotto la
> `ConstraintTolerance` di 1e-6 di `fmincon`: raffinare oltre sarebbe sprecato,
> perche' l'errore dominante diventerebbe quello del solver, non quello
> dell'integratore. La regola pratica e': si sceglie n_sub minimo tale che
> l'errore di integrazione sia almeno un ordine di grandezza sotto la tolleranza
> sui vincoli.

---

## Possibili domande d'esame

**D: Che differenza c'e' fra il difetto trapezoidale e il difetto ZOH+RK4?**
R: Il trapezoidale impone x_{k+1} - x_k = (dt/2)*(f_k + f_{k+1}): e'
un'approssimazione dell'integrale, accurata a O(dt^2), e presuppone un controllo
piecewise-lineare (i due f coinvolgono u_k e u_{k+1}). Il ZOH impone
x_{k+1} = RK4(x_k, u_k, dt): non approssima l'integrale, lo **simula** sotto un
controllo che e' esattamente quello che sara' comandato, con errore O(h^4)
riducibile a piacere. Il primo e' *collocation* (implicito, tutti i nodi
accoppiati due a due), il secondo e' *multiple shooting* (esplicito, ogni
intervallo propagato in avanti). Entrambi producono un NLP con le stesse variabili
di decisione, ma il secondo ha un errore di replay molto piu' piccolo.

**D: Perche' il controllo `u` viene passato identico ai quattro stadi di RK4?**
R: Perche' e' tenuto costante sull'intervallo (zero-order hold). In un RK4 non
autonomo gli stadi campionerebbero il controllo a t_k, t_k+h/2, t_k+h/2 e t_k+h;
essendo u costante quei quattro valori coincidono. Corollario: la dinamica e'
autonoma sull'intervallo, ed e' per questo che `ode_descent` non prende l'
argomento tempo. Se il controllo fosse piecewise-lineare bisognerebbe valutarlo a
tempi diversi e il RHS dovrebbe accettare `t`.

**D: Qual e' l'ordine di accuratezza e come si verifica?**
R: RK4 classico ha errore locale O(h^5) e globale O(h^4). Il test
`testFourthOrderConvergence` lo misura empiricamente: calcola l'errore rispetto a
un riferimento `ode45` (RelTol 1e-12) per n_sub = 1, 2, 4, 8, prende
log2(err_i / err_{i+1}) e verifica che sia > 3.5 (cioe' che dimezzando il passo
l'errore cali di circa un fattore 16). Attenzione: la riga della massa e'
integrata *esattamente* per costruzione, quindi l'ordine misurato riguarda solo le
componenti di posizione e velocita'.

**D: Perche' RK4 e non un metodo implicito o simplettico?**
R: Il problema non e' rigido (non ci sono scale temporali molto separate: massa,
posizione e velocita' evolvono tutte sulla scala di t_ref = 17.5 s) quindi non
serve un metodo implicito, che costerebbe una risoluzione non lineare dentro ogni
sotto-passo, dentro ogni valutazione del vincolo. Non e' nemmeno un problema
hamiltoniano conservativo (c'e' dissipazione di massa e un controllo esterno),
quindi un simplettico non porta vantaggi. RK4 esplicito e' il punto ottimo: alto
ordine, zero risoluzioni implicite, mappa liscia per i gradienti FD.

**D: Il file non ha il blocco `arguments`. Non e' un rischio?**
R: E' una scelta di design dichiarata (righe 12-13) e allineata a `CLAUDE.md`.
Il propagatore vive dentro il loop di `fmincon`, dove viene chiamato ~1e5 volte
per iterazione: la validazione costerebbe piu' del calcolo. Il contratto (x 5x1,
u 2x1, dt > 0, n_sub intero positivo) e' garantito dal chiamante `zoh_nonlcon`,
che a sua volta e' costruito da `solve_zoh`, che **ha** la validazione. Il costo
di questa scelta e' che un bug di dimensioni si manifesterebbe come un errore
criptico di MATLAB dentro `ode_descent`, non come un messaggio chiaro: e' il
prezzo che si paga.
