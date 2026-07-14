# HM2_powered_descent/main_task1.m

## Ruolo del file nel progetto

Questo e' lo script "Task 1" di HM2: risolve il problema di **powered descent
and landing** a combustibile minimo con un **metodo diretto**, cioe'
trascrivendo il problema di controllo ottimo continuo in un **NLP a dimensione
finita** e passandolo a `fmincon`. E' il file di riferimento dell'homework: le
quattro varianti del Task 2 (ZOH+RK4 nonlineare, SCvx LTV con fmincon, SCvx LTV
con YALMIP/ECOS, GFOLD log-massa) sono tutte confronti contro questa
trascrizione, e due di esse partono da una sua soluzione come warm start.

La trascrizione scelta e' la **collocazione diretta trapezoidale** su una griglia
uniforme di `N = 50` nodi, con tempo finale **fissato** (`tf = 38 s`, dato dalla
traccia). Il vettore delle incognite impila stato e controllo a ogni nodo; i
vincoli di collocazione ("defect") impongono che la traiettoria discreta soddisfi
le equazioni del moto; il costo massimizza la massa finale, che e' esattamente
"minimizza il carburante" perche' la massa iniziale e' fissata.

Il file e' **autocontenuto tranne per un'unica dipendenza esterna**:
`ode_descent.m`, il right-hand side della dinamica adimensionale (riga 314,
tramite il wrapper `dyn_rhs`). Tutto il resto -- non dimensionalizzazione, guess
iniziale, bounds, condizioni al contorno, vincoli di percorso, diagnostiche KKT,
studio di convergenza sulla griglia, plot, export figure -- vive dentro
`main_task1.m` come local function. Nessuna dipendenza da toolbox esterni
(YALMIP/ECOS servono solo alle varianti (c)/(d) del Task 2).

Il flusso e': dati SI -> non dimensionalizzazione -> sweep di sensitivita' su
`tf` (tre solve indipendenti a 0.95/1.00/1.05 volte il nominale) -> plot ed
export -> studio di convergenza di griglia (`N = 25/50/100`) che misura la
fedelta' rigiocando il controllo ottimo attraverso `ode45`.

---

## Header e setup (righe 1-12)

```matlab
%  Solved non-dim: L_ref = y0, a_ref = g,
%  t_ref = sqrt(L_ref/g), V_ref = sqrt(g*L_ref),
%  m_ref = m0, T_ref = m0*g. Only residual
%  parameter is V_c = V_ref/c.
```

- Righe 1-10: il blocco di commento dichiara tutte le scelte di progetto --
  trapezoidale, durata fissa, minimo carburante, punto materiale 2D senza
  aerodinamica su Terra piatta, solver `fmincon` con algoritmo `sqp`.
- Riga 12: `clear; close all; clc;`. Il `close all` e' funzionale: il blocco di
  export figure (righe 62-78) fa `findobj(groot, 'Type', 'figure')` e salva
  *tutte* le figure aperte, quindi lo script deve partire da un desktop pulito.

**Il problema di controllo ottimo continuo.** Il modello e' un punto materiale
2D in un campo di gravita' uniforme, senza aerodinamica. Stato
`x = [x, y, vx, vy, m]`, controllo `u = [Tx, Ty]` (vettore spinta cartesiano).
In forma dimensionale:

    x_dot  = vx
    y_dot  = vy
    vx_dot = Tx / m
    vy_dot = Ty / m - g
    m_dot  = -|T| / c,     con |T| = sqrt(Tx^2 + Ty^2),  c = Isp * g0

Il problema e':

    min  int_0^tf |T(t)| dt        (equivalente: max m(tf))
    s.t. dinamica di sopra
         x(0) = (1000, 3000) m,  v(0) = (300, -200) m/s,  m(0) = 2000 kg
         x(tf) = 0, y(tf) = 0, vx(tf) = 0, vy(tf) = 0     (pinpoint + soft)
         Tmin <= |T(t)| <= Tmax
         |x(t)| <= tan(theta_max) * y(t)                  (glide slope)
         tf = 38 s fissato

Non c'e' **nessun angolo di pitch** nel modello: essendo un punto materiale,
l'assetto non e' uno stato e non esiste un vincolo di puntamento della spinta
(del tipo `T_y >= |T| cos(theta_pointing)`). L'unico vincolo angolare e' il
glide slope, che e' un vincolo **geometrico sulla posizione**, non sull'assetto.
Questo e' fedele alla traccia dell'homework, non una semplificazione arbitraria.

> **Possibile domanda d'esame** -- Perche' `min int |T| dt` e' equivalente a
> `max m(tf)`?
> *Risposta:* perche' `m_dot = -|T|/c` con `c` costante, quindi integrando
> `m(tf) = m0 - (1/c) * int_0^tf |T| dt`. Con `m0` fissata, massimizzare `m(tf)`
> e' identico a minimizzare l'integrale del modulo della spinta. Il vantaggio
> pratico e' che la massa e' gia' uno stato del sistema: il costo diventa di tipo
> **Mayer** (una funzione del solo stato finale), quindi nel NLP e' semplicemente
> `-z(i_massa_ultimo_nodo)` -- lineare nelle variabili di decisione, con gradiente
> costante. Non serve nessuna quadratura del costo.

---

## `Problem data (Table 1, dimensional)` (righe 14-29)

- Righe 15-26: i dati della Tabella 1 della traccia, in SI. Stato iniziale
  `(1000, 3000) m`, `(300, -200) m/s`, `m0 = 2000 kg`; `Isp = 225 s`;
  `Tmin = 0`, `Tmax = 70 kN`; `theta_mx = 60 deg` (semiapertura del cono di
  glide slope).
- Riga 20 vs riga 22: **due gravita' diverse, e non e' un errore**. `data.g =
  9.81` e' la gravita' locale che entra nelle equazioni del moto; `data.g0 =
  9.80665` e' la gravita' standard, usata **solo** alla riga 23 per convertire
  l'impulso specifico in velocita' di efflusso, `c = Isp * g0 = 2206.5 m/s`. La
  definizione di `Isp` e' per convenzione legata a `g0` standard, non alla
  gravita' del sito.
- Riga 24: `Tmin = 0`. E' il valore della traccia, e ha una conseguenza pesante
  sulla struttura del problema (vedi `trap_nonlcon`): con `Tmin = 0` l'insieme
  ammissibile della spinta e' il **disco** `|T| <= Tmax`, che e' convesso. Se
  fosse `Tmin > 0` sarebbe una **corona circolare**, non convessa -- il caso
  classico che motiva la lossless convexification (ticket T006, ancora aperto).
- Riga 28: `tf_nom = 38` s, tempo di volo **fissato dalla traccia**, non ottimizzato.
- Riga 29: `N = 50` nodi di collocazione.

---

## `Non-dimensionalisation` (righe 31-35)

- Riga 32: chiama `nondim(data)` (definita a riga 109) e stampa le scale.
- Righe 33-35: `fprintf` diagnostici. Il numero che conta e' `V_c = V_ref/c`,
  perche' e' l'**unico parametro residuo** del problema adimensionale.

**Perche' e' l'unico parametro residuo.** Con `L_ref = y0`, `a_ref = g`,
`t_ref = sqrt(L_ref/g)`, `V_ref = sqrt(g*L_ref) = L_ref/t_ref`, `m_ref = m0`,
`T_ref = m0*g`, le equazioni si riscrivono cosi' (tilde = adimensionale):

    dx~/dt~  = vx~                          (coefficiente 1 per costruzione)
    dvx~/dt~ = (Tx/m) / a_ref
             = (Tx~ * T_ref) / (m~ * m_ref) / g
             = Tx~ / m~                     perche' T_ref = m_ref * g
    dvy~/dt~ = Ty~/m~ - 1                   perche' g/a_ref = 1
    dm~/dt~  = -(|T|/c) * t_ref/m_ref
             = -|T~| * (m_ref*g*t_ref) / (c*m_ref)
             = -|T~| * (g*t_ref)/c
             = -|T~| * V_ref/c              perche' g*t_ref = sqrt(g*L_ref) = V_ref

Cioe': tutte le equazioni hanno coefficienti unitari **tranne** quella della
massa, che porta il gruppo `V_c = V_ref/c` (un numero di Tsiolkovsky). Questa e'
esattamente la riga 14 di `ode_descent.m`:

```matlab
dx = [ x(3); x(4); u(1)/x(5); u(2)/x(5) - 1; -Vc*Tmag ];
```

Coi numeri della traccia: `t_ref = 17.49 s`, `V_ref = 171.6 m/s`,
`T_ref = 19.62 kN`, `V_c = 0.0778`. Le condizioni iniziali adimensionali
diventano `x0~ = 1/3`, `y0~ = 1`, `vx0~ = 1.75`, `vy0~ = -1.17`, `m0~ = 1`, e
`Tmax~ = 70000/19620 = 3.57` -- che si legge subito come **rapporto
spinta/peso massimo 3.57**. Il tempo di volo nominale diventa `tf~ = 2.17`.

Il beneficio per il solver e' il **condizionamento**: in SI le variabili di
decisione andrebbero da `1e-1` (velocita' finali) a `7e4` (spinta), cinque
ordini di grandezza. `fmincon/sqp` usa differenze finite con passo relativo
fisso e un test di ottimalita' su una norma mista di gradiente e moltiplicatori:
con variabili tutte O(1) il passo di differenza finita e' significativo per
tutte, e la tolleranza `OptimalityTolerance = 1e-5` ha lo stesso peso su ogni
componente.

> **Possibile domanda d'esame** -- La `ConstraintTolerance = 1e-6` e'
> adimensionale: quanto vale in metri?
> *Risposta:* i defect di posizione sono normalizzati su `L_ref = 3000 m`, quindi
> `1e-6` non-dim = `3 mm` di errore di posizione per intervallo; i defect di
> velocita' sono su `V_ref = 171.6 m/s`, quindi `1e-6` non-dim = `0.17 mm/s`. La
> tolleranza e' fisicamente molto stretta. Attenzione a non confondere questo con
> l'errore di **trascrizione**, che e' molto piu' grande (vedi la sezione sullo
> studio di griglia): il defect misura quanto bene si soddisfa la dinamica
> *discreta*, non quanto la dinamica discreta approssima quella *vera*.

---

## `Sensitivity sweep on flight time` (righe 37-57)

- Riga 38: `tf_list = tf_nom * [0.95, 1.00, 1.05]` -> `[36.10, 38.00, 39.90] s`.
- Righe 41-57: tre solve **completamente indipendenti**. Non c'e' continuazione /
  warm start fra un `tf` e il successivo: ogni run riparte dalla guess lineare
  costruita dentro `solve_trapcol`. E' una scelta discutibile (la README lo
  ammette: la continuazione e' elencata fra i possibili rimedi allo stallo di
  convergenza), ma rende i tre risultati indipendenti fra loro.
- Riga 43: `tf_nd = tf_list(k) / ref.t` -- il solve avviene tutto in
  adimensionale.
- Riga 45: `dim_sol` riporta la soluzione in SI **solo per stampa e plot**.
- Righe 48-54: diagnostiche post-solve (tempi di switch, coast, margine di glide
  slope, attivita' KKT).
- Righe 55-56: stampa `iters`, `firstorderopt`, `exitflag` di `fmincon`.

**Attenzione: questo NON e' un'ottimizzazione a `tf` libero.** `tf` non e' una
variabile di decisione del NLP: e' un **parametro** che entra solo tramite
`dt = tf/(N-1)`. Lo sweep e' un'analisi di sensitivita' a tre punti, non una
ricerca dell'ottimo. Dai risultati riportati nella README (593.67 / 596.80 /
601.18 kg di carburante per `tf` crescente) il consumo **cresce
monotonicamente** con `tf` nella finestra esplorata: significa che, dentro quella
finestra, `tf = 38 s` non e' il tempo di volo a minimo carburante, e che
l'ottimo (se esiste) sta a sinistra del bordo inferiore, sotto i 36.10 s. Lo
script non lo cerca, e non lo puo' dire.

L'interpretazione fisica del trend e' la gravity loss: piu' a lungo si sta in
volo, piu' impulso serve semplicemente per sostenere il peso durante il coast.
Ma il trade-off ha un altro lato che questa finestra non vede: accorciando
troppo `tf` la decelerazione richiesta supererebbe `Tmax` e il problema
diventerebbe infeasible. Lo sweep esplora solo il ramo monotono.

> **Possibile domanda d'esame** -- Se volessi rendere `tf` libero, cosa
> cambieresti?
> *Risposta:* due strade. (1) Aggiungere `tf` come 7N+1-esima variabile di
> decisione, con `dt = tf/(N-1)` che compare dentro i defect: i defect diventano
> **nonlineari anche in `tf`** e la Jacobiana acquista una colonna densa. Il
> costo `-m_N` resta lineare. (2) Tenere `tf` come parametro e mettere una
> ricerca 1D esterna (bisezione / golden section) sul valore ottimo, sfruttando il
> fatto che ogni solve interno resta piu' semplice. Nel piano del ticket T006 la
> seconda strada e' quella scelta per la variante SOCP, proprio perche' con la
> lossless convexification il problema interno e' convesso e la ricerca su `tf`
> conserva un certificato di ottimalita' globale per ogni `tf` provato.

---

## `Plots` / `Export figures` / `Summary table` (righe 59-85)

- Riga 60: chiama `plot_results` (definita a riga 389).
- Righe 63-78: export. `slugify` (riga 65) converte il `Name` della figura in un
  nome file; `theme(fig, 'light')` (riga 71) forza il tema chiaro anche se il
  desktop MATLAB e' in dark mode, col `catch` alla riga 73 che ricade su
  `Color = 'w'` per MATLAB pre-R2025a. Export in PNG a 200 dpi in
  `HM2_powered_descent/figures/`.
- Righe 81-85: la tabella di sintesi `tf | m_f | fuel`.

Nota di onesta': l'export salva **tutte** le figure attualmente aperte, non solo
quelle create da `plot_results`. Funziona perche' la riga 12 ha fatto
`close all`, ma e' fragile se lo script venisse eseguito per sezioni.

---

## `Grid-convergence study` (righe 87-101)

- Riga 90: `N_list = [25, 50, 100]`, tutti al `tf` nominale.
- Riga 95: risolve di nuovo il NLP (nota: il caso `N = 50` viene ricalcolato da
  zero, era gia' stato risolto nello sweep -- spreco innocuo).
- Riga 97: `fwd_integrate_pwl` rigioca il controllo ottimo attraverso `ode45`.
- Riga 98: `err = max(node_err(...))` -- massimo errore nodale su posizione e
  velocita', in adimensionale.
- Riga 99-100: stampa `N`, massa finale in kg, errore, wall time.

**Che cosa misura davvero questa metrica.** Il NLP e' (quasi) esattamente
ammissibile per la dinamica **discreta**: i defect sono azzerati a `1e-6`. Ma la
dinamica discreta e' un'approssimazione `O(dt^2)` di quella vera. Il replay con
`ode45` (tolleranze `1e-10 / 1e-12`, quindi praticamente esatto) misura la
distanza fra le due: e' l'**errore di trascrizione**, non una violazione di
vincolo. La README riporta che l'errore cala di un fattore ~4 a ogni raddoppio di
`N` -- cioe' il fattore atteso `(dt/2)^2 / dt^2 = 1/4`. E' la conferma numerica
dell'ordine 2 della regola trapezoidale.

E' anche il motivo per cui, pur avendo `y(tf) = 0` imposto come vincolo di
uguaglianza **esatto**, il replay open-loop atterra a 4.3 m dal pad a 0.11 m/s
(numeri dalla README). Il vincolo e' soddisfatto -- ma sulla traiettoria
*discreta*.

> **Possibile domanda d'esame** -- Il tuo NLP dice `y_N = 0` esattamente, ma
> rigiocando il controllo con `ode45` atterri 4 metri fuori. La tua soluzione e'
> sbagliata?
> *Risposta:* no, e' *ammissibile per il problema che ho scritto*. Il NLP impone
> la dinamica **trapezoidale**, non quella continua; il gap di 4.3 m e' l'errore
> di quadratura `O(dt^2)` con `dt = 0.78 s`. Si chiude in tre modi: infittendo la
> griglia (`N = 100` lo riduce di ~4x), passando a uno schema di ordine
> superiore (Hermite-Simpson, `O(dt^4)`), oppure -- che e' quello che si fa in
> volo -- non usando mai il piano open-loop: si ri-risolve il problema a ogni ciclo
> di guida (MPC) o si mette un tracker in retroazione sopra la traiettoria
> nominale.

---

## `nondim` (righe 109-134)

- Riga 109: firma `[ref, dnd] = nondim(d)`. Chiamata una volta sola (riga 32),
  quindi ha un blocco `arguments` (righe 116-118): la validazione sta sui
  **boundary helper**, non sulle funzioni hot-loop -- convenzione della repo.
- Righe 119-124: costruisce le scale. `ref.L = y0` (quota iniziale),
  `ref.g = g`, `ref.t = sqrt(L/g)`, `ref.V = sqrt(g*L)`, `ref.m = m0`,
  `ref.T = m0*g`. La scelta `T_ref = m_ref * g` e' quella che fa collassare i
  coefficienti nelle equazioni di velocita' (vedi derivazione sopra).
- Righe 125-133: applica le scale. `dnd.m0 == 1` per costruzione (riga 129, il
  commento lo dice). `dnd.Vc = ref.V / d.c` (riga 132) e' il gruppo residuo.
- Riga 133: `theta_mx` **non viene scalato** -- e' gia' un angolo, adimensionale
  per natura.

---

## `dim_sol` (righe 136-164)

- Righe 147-158: rimoltiplica ogni campo per la sua scala. `sol.Tmag` (riga 155)
  viene ricalcolata in SI da `Tx, Ty` invece di essere scalata: equivalente,
  perche' la norma e' omogenea di grado 1.
- Riga 158: `fuel = (m0~ - m_f~) * m_ref` -- il carburante consumato.
- Righe 159-163: **i moltiplicatori di Lagrange restano adimensionali** e il
  commento lo dichiara esplicitamente. E' corretto non scalarli con una scala di
  stato (avrebbero dimensioni miste: `d(costo)/d(vincolo)`), ma significa che i
  numeri stampati da `diagnostics` (riga 348, `max_gs_mult`) vanno letti solo in
  senso relativo (zero / non zero), mai come grandezze fisiche.

---

## `solve_trapcol` (righe 166-267) -- il cuore della trascrizione

E' la funzione che costruisce e risolve il NLP. Vale la pena smontarla pezzo per
pezzo.

### Il vettore delle variabili di decisione (righe 181-183)

```matlab
dt  = tf / (N - 1);
nz  = 7 * N;
idx = @(i) (i-1)*7 + (1:7);
```

- Riga 181: griglia **uniforme**, `dt = tf/(N-1)`. Con `tf = 38 s` e `N = 50`:
  `dt = 0.776 s` (`dt~ = 0.0443`).
- Riga 182: `nz = 7*N = 350` incognite.
- Riga 183: l'ordinamento. `idx(i)` restituisce le 7 posizioni del nodo `i`:

      z = [ x_1 y_1 vx_1 vy_1 m_1 Tx_1 Ty_1 |
            x_2 y_2 vx_2 vy_2 m_2 Tx_2 Ty_2 | ... | ...nodo N ]

  Cioe' **stacking node-major**: prima tutto il nodo 1, poi tutto il nodo 2, ecc.
  Il controllo e' discretizzato **ai nodi**, non agli intervalli -- coerente con la
  trapezoidale, che implica un controllo **piecewise-linear**.

**Perche' proprio quest'ordine?** Perche' rende la Jacobiana dei vincoli
**a banda**. Il defect dell'intervallo `k` coinvolge solo i nodi `k` e `k+1`,
cioe' 14 componenti **contigue** di `z` (le posizioni `7(k-1)+1 ... 7(k+1)`). Con
l'ordinamento alternativo -- tutte le `x`, poi tutte le `y`, poi tutte le `vx`,
ecc. -- le stesse 14 componenti sarebbero sparse per tutto il vettore, e la
Jacobiana perderebbe la struttura a banda. Nota pero' l'onesta' del caso: **qui il
codice non sfrutta questa proprieta'**, perche' non dichiara nessuno sparsity
pattern e non fornisce Jacobiane analitiche (vedi sotto). L'ordinamento e'
comunque quello giusto, e sarebbe l'unico da usare se si passasse a `interior-point`
con `ConstraintGradient` sparso.

Un altro vantaggio pratico: `reshape(z, 7, N)` (riga 280, dentro `trap_nonlcon`)
recupera la matrice nodi-per-colonna **senza copie ne' permutazioni**, perche'
MATLAB e' column-major. Con l'ordinamento alternativo servirebbe una
trasposizione. (Attenzione a non confondere con l'unpack finale a riga 251,
`Z = reshape(z_opt, 7, N).'`: li' la trasposizione c'e', ed e' voluta, perche'
serve la matrice nodi-per-**riga** da cui estrarre le colonne `sol.x`, `sol.y`,
ecc. Ma e' fatta una volta sola a fine solve, non nel hot loop.)

### La guess iniziale (righe 185-197)

```matlab
a = (i-1)/(N-1);
z0(s(1)) = (1-a) * d.x0;      % ... idem y, vx, vy
z0(s(5)) = d.m0 * (1 - 0.3*a);
z0(s(6)) = 0;
z0(s(7)) = d.m0;              % hover
```

- Righe 188-193: **interpolazione lineare** dallo stato iniziale a **zero**
  (perche' `(1-a)` va da 1 a 0). Non e' una propagazione: e' la corda che unisce
  i due estremi. Siccome i target finali di `x, y, vx, vy` sono tutti nulli, la
  corda arriva esattamente sulle condizioni finali.
- Riga 194: massa da `m0` a `0.7*m0`: si indovina il 30% di frazione di
  carburante. (La soluzione vera consuma ~597 kg su 2000, cioe' 29.8% -- la guess
  e' sospettosamente centrata, verosimilmente tarata a posteriori.)
- Righe 195-196: controllo costante `T = (0, m0)` adimensionale, cioe' spinta
  verticale pari al **peso iniziale** (perche' `T_ref = m0*g`): la guess di
  hover. Il commento "gravity = 1 nondim" e' corretto.

**Proprieta' notevole della guess** (che vale la pena saper difendere): questa
guess soddisfa **tutti i vincoli tranne i defect**. Le BC lineari sono esatte per
costruzione. Il vincolo di spinta e' rispettato (`|T| = 1 < Tmax~ = 3.57`). E il
glide slope pure: lungo la corda `x/y = x0/y0` e' **costante** (entrambi sono
scalati dallo stesso `(1-a)`), quindi l'angolo resta `atan(1000/3000) = 18.4 deg
< 60 deg` su tutti i nodi. Quindi l'unica infeasibility che `fmincon` deve
sanare e' quella dinamica. E' esattamente la situazione in cui SQP lavora bene.

### Bounds (righe 199-211)

- Riga 204: `lb(y) = 0` -- niente volo sottoterra.
- Righe 205-206: `1e-3 <= m <= m0`. Il lower bound e' **essenziale**, non
  cosmetico: il RHS ha `Tx/m` e `Ty/m`, quindi `m -> 0` fa esplodere la dinamica.
  Il bound superiore `m <= m0` e' ridondante fisicamente (la massa non puo'
  crescere: dal defect di massa, `m_{k+1} - m_k = -0.5*dt*Vc*(|T_k| + |T_{k+1}|)
  <= 0` sempre) ma da' a SQP una scatola finita su cui lavorare.
- Righe 207-210: `-Tmax <= Tx, Ty <= Tmax`. E' una **scatola**, non un disco: da
  sola permetterebbe `|T|` fino a `sqrt(2)*Tmax`. E' il vincolo nonlineare in
  `trap_nonlcon` a imporre il disco. La scatola serve solo a dare bounds finiti
  a tutte le variabili.

Nota di onesta': `lb(Ty) = -Tmax` **permette spinta verso il basso**. Fisicamente
un lander col motore sotto non puo' spingere all'ingiu'. La traccia non impone un
vincolo di puntamento, e all'ottimo la spinta verso il basso non viene mai usata
(sprecherebbe carburante e violerebbe il glide slope), ma il modello lo
consentirebbe. E' una liberta' non fisica lasciata aperta.

### Condizioni al contorno (righe 213-225)

```matlab
Aeq = sparse(9, nz);
Aeq(1, s1(1)) = 1;  beq(1) = d.x0;   % ... 5 iniziali
Aeq(6, sN(1)) = 1;  beq(6) = 0;      % ... 4 finali
```

- Righe 217-221: **5** BC iniziali (`x, y, vx, vy, m` al nodo 1).
- Righe 222-225: **4** BC finali (`x, y, vx, vy` al nodo N a zero: pinpoint +
  soft landing). **La massa finale NON e' vincolata** -- e' l'obiettivo.
- Sono messe in `Aeq / beq`, cioe' fra le **equalita' lineari**, non dentro
  `nonlcon`. E' la scelta giusta: `fmincon` le tratta esattamente (le usa nel QP
  senza linearizzarle) e soprattutto **non le differenzia numericamente**.

Piccola stonatura: `Aeq` viene costruita con `sparse` (riga 215) e poi passata a
`fmincon` con `full(Aeq)` (riga 243), il che annulla il beneficio. Con una 9x350
non cambia nulla in pratica, ma la `sparse` e' inutile.

### Il costo (righe 227-229)

```matlab
iN_m  = (N-1)*7 + 5;
f_obj = @(z) -z(iN_m);
```

- Riga 228: `iN_m` e' l'indice della **massa all'ultimo nodo**. Con `N = 50`:
  `49*7 + 5 = 348`.
- Riga 229: minimizza `-m_N`, cioe' massimizza la massa finale.

Il costo e' **lineare** nelle variabili di decisione. Il suo gradiente e' il
vettore `-e_348`, esatto e costante. Ne segue una conseguenza importante: nel
NLP **tutta la nonlinearita' vive nei vincoli** (i defect e il modulo della
spinta), non nel costo.

### Le opzioni di `fmincon` (righe 234-244)

```matlab
opts = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...
    'Display',   'final', ...
    'MaxIterations',          1000, ...
    'MaxFunctionEvaluations', 1e6, ...
    'OptimalityTolerance',    1e-5, ...
    'ConstraintTolerance',    1e-6, ...
    'StepTolerance',          1e-10);
```

- `Algorithm = 'sqp'`: Sequential Quadratic Programming. A ogni iterazione
  costruisce un QP locale (Hessiano del Lagrangiano approssimato BFGS, vincoli
  linearizzati) e ne prende la soluzione come direzione di ricerca. Gestisce bene
  i punti di partenza **infeasible** -- che e' esattamente il caso della guess
  (vedi sopra). Lo script non motiva la scelta.
- `ConstraintTolerance = 1e-6`, `OptimalityTolerance = 1e-5`,
  `StepTolerance = 1e-10`, `MaxIterations = 1000`. Anche `Display = 'final'`
  (riga 236) e' esposta esplicitamente: convenzione della repo, tutte le opzioni
  che determinano il comportamento del solver stanno scritte nello script, cosi'
  il run e' riproducibile.
- Righe 246-248: se `exitflag <= 0` emette una `warning`. **Questo warning
  scatta davvero**: la README documenta che i run a `tf = 36.10 s` e
  `tf = 38.00 s` finiscono col cap di `MaxIterations` (`exitflag = 0`) e
  first-order optimality bloccata a `1e-3`-`1e-4` non-dim. La feasibility resta
  invece sotto `1e-6`.

**Il buco nero: nessun gradiente analitico.** Le opzioni **non** contengono
`SpecifyObjectiveGradient` ne' `SpecifyConstraintGradient`, e nessuno sparsity
pattern viene dichiarato. Quindi `fmincon` calcola **tutte** le derivate per
differenze finite in avanti. Il conto: 350 variabili -> 351 valutazioni di
`trap_nonlcon` per ogni Jacobiana, ognuna delle quali fa un loop di 50 chiamate a
`ode_descent`. E' il motivo per cui il singolo solve costa ~1-1.5 minuti. Sono
tutte derivate che si scriverebbero a mano in poche righe:

- gradiente del costo: `-e_{iN_m}`, esatto e costante;
- Jacobiana dei defect: `d(zeta_k)/d(z_k) = -I - (dt/2)*df/dx|_k` e
  `d(zeta_k)/d(z_{k+1}) = +I - (dt/2)*df/dx|_{k+1}`, con `df/dx` un blocco 5x7
  analitico (le uniche derivate non banali sono `d(T/m)/dm = -T/m^2` e
  `d|T|/dT = T/|T|`);
- Jacobiana dei vincoli di percorso: quelli di glide slope sono **lineari**,
  quelli di spinta hanno gradiente `+/- T/|T|`.

Sarebbe una vittoria di uno-due ordini di grandezza in tempo di calcolo. Non e'
stato fatto: e' un limite reale dello script, non una scelta.

---

## `trap_nonlcon` (righe 269-304) -- defect e vincoli di percorso

```matlab
Z = reshape(z, 7, N);
for i = 1:N
    f(:,i) = dyn_rhs(Z(:,i), d.Vc);
end
for k = 1:N-1
    defs(:,k) = Z(1:5,k+1) - Z(1:5,k) ...
                - 0.5*dt*(f(:,k) + f(:,k+1));
end
```
*(indentazione e spaziatura ricompattate rispetto al sorgente, righe 280-292.)*

- Riga 279: il commento dichiara che **non c'e'** blocco `arguments`, per scelta:
  e' una funzione hot-loop dentro `fmincon`. Convenzione della repo, non una
  dimenticanza.
- Riga 280: `reshape(z, 7, N)` -- nodi in colonna. Gratis, per l'ordinamento
  node-major.
- Righe 283-286: valuta il RHS a **tutti** i nodi. Nota che `f` e' 5xN mentre `Z`
  e' 7xN: la dinamica ha 5 stati, le altre 2 righe sono il controllo.
- Righe 289-292: i **defect**.

### Derivazione della regola trapezoidale

Su un intervallo `[t_k, t_{k+1}]`, integrando esattamente l'equazione di stato:

    x_{k+1} - x_k = int_{t_k}^{t_{k+1}} f(x(t), u(t)) dt

L'integrale non e' calcolabile (non conosciamo `x(t)`). Lo si approssima con una
quadratura. Se si approssima l'**integrando** `f` con la sua interpolante
**lineare** fra i due estremi noti `f_k = f(x_k, u_k)` e `f_{k+1}`:

    f(t) ~ f_k + (f_{k+1} - f_k) * (t - t_k)/dt

e si integra questa retta su `[t_k, t_{k+1}]`, si ottiene l'area del trapezio:

    int f dt ~ (dt/2) * (f_k + f_{k+1})

Da cui il **defect** (il "difetto", cioe' il residuo dell'equazione discretizzata):

    zeta_k = x_{k+1} - x_k - (dt/2)*(f_k + f_{k+1})  =  0,   k = 1 ... N-1

Sono **vincoli di uguaglianza** che il solver deve azzerare: quando `zeta_k = 0`
per ogni `k`, la sequenza di nodi `{x_k}` e' una traiettoria valida del sistema
**discreto**. Se il solver li lasciasse diversi da zero, i nodi sarebbero
semplicemente numeri che non stanno su nessuna traiettoria fisica. Da qui il nome:
il "difetto" e' quanto manca alla traiettoria per essere ammissibile.

Osservazioni chiave:

- **Il metodo e' implicito.** `f_{k+1}` dipende da `x_{k+1}`, che e' incognito. In
  uno schema di integrazione classico questo richiederebbe di risolvere
  un'equazione nonlineare a ogni passo. Nella collocazione diretta **non serve**:
  tutte le `x_k` sono variabili di decisione, e l'implicitezza e' semplicemente
  parte del sistema di vincoli che il NLP risolve *simultaneamente*. E' il
  vantaggio strutturale della trascrizione (a volte detta "simultaneous
  approach", in opposizione al "sequential"/shooting).
- **Ordine di accuratezza.** L'errore della quadratura trapezoidale su un
  intervallo e' `-(dt^3/12) * f''(xi)`: locale `O(dt^3)`, quindi globale
  `O(dt^2)`. Il metodo coincide col trapezio implicito (Crank-Nicolson):
  A-stabile, simmetrico, secondo ordine. Lo studio di griglia (righe 87-101) lo
  conferma empiricamente col fattore 4 per raddoppio di `N`.
- **Modello implicito del controllo.** Poiche' l'integrando e' approssimato
  linearmente, anche il controllo e' implicitamente **piecewise-linear** fra i
  nodi. Per questo `fwd_integrate_pwl` (riga 351) rigioca il controllo con
  interpolazione **lineare** e non ZOH: rigiocarlo con uno ZOH misurerebbe la
  fedelta' della trascrizione sbagliata.

**Confronto con la trascrizione ZOH/RK4 del Task 2.** Le due sono
concettualmente diverse:

| | trapezoidale (Task 1) | ZOH + RK4 (Task 2a) |
|---|---|---|
| controllo | piecewise-linear ai nodi | costante a tratti sull'intervallo |
| vincolo di uguaglianza | `x_{k+1} - x_k - (dt/2)(f_k+f_{k+1}) = 0` | `x_{k+1} - RK4(x_k, u_k, dt) = 0` |
| `x_{k+1}` compare | implicitamente (dentro `f_{k+1}`) | solo a sinistra |
| ordine | 2 | 4 (con substep, `rk4_zoh.m`) |
| fedelta' al replay | ~4.3 m (README) | ~`1.4e-8` non-dim (README) |

La colonna di destra e' quello che fa `rk4_zoh.m`: propaga `ode_descent` con
Runge-Kutta 4 tenendo `u` **costante** sull'intervallo. E' anche un vincolo di
uguaglianza, ma di tipo **multiple shooting** (il difetto e' fra lo stato
propagato e il nodo successivo), non di collocazione. Attenzione al terzo file
della famiglia: `ode_descent_uacc.m` tiene costante l'**accelerazione**
`u = T/m` invece del vettore spinta -- e' la convenzione ZOH nativa della
formulazione GFOLD log-massa (Task 2d), e non e' intercambiabile con
`ode_descent.m`.

### I vincoli di percorso (righe 295-303)

```matlab
Tmag     = sqrt(Z(6,:).^2 + Z(7,:).^2).';
g_thr_lo = d.Tmin - Tmag;          % <= 0
g_thr_hi = Tmag   - d.Tmax;        % <= 0
tt       = tan(d.theta_mx);
g_gs_pos = ( Z(1,:).' - tt*Z(2,:).');
g_gs_neg = (-Z(1,:).' - tt*Z(2,:).');
c_ineq   = [g_thr_lo; g_thr_hi; g_gs_pos; g_gs_neg];
```

**Vincolo di spinta** (righe 296-298). Scritto come due disuguaglianze scalari
per nodo:

    Tmin - |T_k| <= 0        e        |T_k| - Tmax <= 0

E' la trascrizione letterale di `Tmin <= |T| <= Tmax`. La domanda cruciale e':
**e' convesso?**

- L'insieme `{T : |T| <= Tmax}` e' un **disco**: convesso. Il vincolo
  `|T| - Tmax <= 0` e' una funzione convessa `<= 0`, quindi definisce un insieme
  convesso. OK.
- L'insieme `{T : |T| >= Tmin}` con `Tmin > 0` e' il **complemento** di un disco:
  **non convesso**. Il vincolo `Tmin - |T| <= 0` e' `-|T| + Tmin <= 0`, cioe' una
  funzione **concava** `<= 0`. La coppia definisce una **corona circolare**, che
  non e' convessa.
- **Ma qui `Tmin = 0`** (riga 24). Il vincolo inferiore diventa `-|T| <= 0`, che
  e' sempre vero: **e' vacuo**. L'insieme ammissibile della spinta e' quindi il
  solo disco, ed e' convesso.

Quindi: il vincolo di spinta *come scritto* e' convesso, per il valore di `Tmin`
della traccia. **Ma il NLP nel suo complesso resta non convesso**, perche' i
defect contengono `Tx/m` e `Ty/m` (bilineari nelle incognite) e `|T|` dentro la
dinamica di massa. `fmincon` restituisce un **minimo locale**, senza alcun
certificato di ottimalita' globale. Questo e' il punto che motiva l'intera linea
del Task 2 e il ticket T006 (lossless convexification): col cambio di variabili
log-massa `z = ln m`, `u = T/m`, piu' lo slack `Gamma >= |u|`, l'intero problema
diventa **un solo SOCP** convesso, risolto in meno di un secondo con certificato
globale -- e la convexification e' *lossless* nel senso che l'ottimo del problema
rilassato e' dimostrabilmente anche ottimo del problema originale.

**Il vero problema numerico: `|T|` non e' differenziabile in `T = 0`.** La norma
euclidea ha un cono (una punta) nell'origine; il suo gradiente `T/|T|` li' non
esiste. E la soluzione ottima **ci passa in mezzo**: durante il coast arc la
spinta e' esattamente zero. `fmincon`, che sta usando differenze finite in avanti
con passo `~1.5e-8`, in quel punto calcola una pendenza unilatera e sistematicamente
di modulo ~1, indipendentemente dalla direzione -- un gradiente **sbagliato**. `|T|`
compare sia nel vincolo `g_thr_lo/hi` sia (via `ode_descent`) nel defect di massa,
quindi la contaminazione e' doppia. **Questa e' la causa radice dello stallo di
convergenza** documentato nella README (first-order optimality bloccata a
`1e-3`-`1e-4`, run che finiscono a `MaxIterations`). Non e' un bug: e' una
patologia strutturale della formulazione min-fuel con norma nella dinamica.

**Vincolo di glide slope** (righe 299-301). Il vincolo fisico e'

    |x| <= tan(theta_max) * y,    theta_max = 60 deg

cioe': il veicolo deve stare dentro un cono di semiapertura 60 gradi **misurata
dalla verticale**, con vertice sul pad. Il valore assoluto lo rende non
differenziabile in `x = 0`, ma si spacca esattamente nella **coppia di
disuguaglianze lineari**:

    +x - tan(theta_max)*y <= 0
    -x - tan(theta_max)*y <= 0

Le due insieme sono equivalenti a `|x| <= tan(theta)*y` -- e sono **lineari**,
quindi convesse e differenziabili ovunque. E' lo stesso trucco che si usa per
linearizzare una norma-1. (Il fatto che `y >= 0` sia garantito dal bound alla
riga 204 e' quello che rende l'equivalenza corretta: con `y < 0` le due
disuguaglianze sarebbero infeasible, il che e' comunque il comportamento voluto.)

Nota di onesta': essendo **lineari**, questi vincoli dovrebbero stare negli slot
`A`, `b` di `fmincon` (che alla riga 243 sono passati vuoti, `[]`), non dentro
`nonlcon`. Cosi' come sono, `fmincon` li differenzia numericamente insieme a
tutto il resto -- `2N = 100` righe di Jacobiana calcolate a differenze finite
quando sarebbero note esattamente e costanti. Innocuo in accuratezza, sprecato in
tempo.

**Un limite strutturale della collocazione, valido per entrambi i vincoli.** I
vincoli di percorso sono imposti **solo ai nodi**. Fra un nodo e l'altro lo stato
e' una spline (quadratica, per la trapezoidale) e nulla vieta che sconfini dal
corridoio o che `|T|` superi `Tmax`. La `gs_margin` calcolata in `diagnostics` e'
un margine **nodale**, non un certificato sull'intera traiettoria. Con `N = 50` e
un margine di 1-2 gradi (README), il rischio e' basso ma non nullo.

**Ordine di impilamento di `c_ineq`** (riga 303): `[thr_lo; thr_hi; gs_pos;
gs_neg]`, `N` elementi ciascuno, `4N = 200` in totale. **Quest'ordine e' un
contratto**: `diagnostics` (righe 347-348) slicia `lambda.ineqnonlin` esattamente
su questi blocchi (`N+1:2N` per il bound superiore di spinta, `2N+1:4N` per il
glide slope). Riordinare qui e non li' rompe silenziosamente le diagnostiche KKT.

> **Possibile domanda d'esame** -- Quanti vincoli e quante incognite ha il tuo
> NLP, e quanti gradi di liberta' restano?
> *Risposta:* con `N = 50`: **350** variabili (`7 x 50`); **245** uguaglianze
> nonlineari (5 defect x 49 intervalli) piu' **9** uguaglianze lineari (5 BC
> iniziali + 4 finali) = 254 uguaglianze; **200** disuguaglianze (`4 x 50`), piu'
> i bound di scatola. Restano `350 - 254 = 96` gradi di liberta' sulla varieta'
> ammissibile, sui quali si minimizza `-m_N` rispettando le 200 disuguaglianze.

---

## `dyn_rhs` (righe 306-315)

```matlab
dx = ode_descent(s(1:5), s(6:7), Vc);
```

- Riga 314: e' solo un **adattatore**. `ode_descent` vuole stato e controllo
  separati; il NLP li tiene impilati in un unico vettore di 7 elementi per nodo.
  `dyn_rhs` spacca `s(1:5)` (stato) da `s(6:7)` (controllo).
- Riga 313: nessuna validazione `arguments`, di nuovo per scelta dichiarata (hot
  loop dentro `fmincon` **e** dentro `ode45`, riga 372).

Il fatto che l'unica dipendenza esterna dello script passi da qui e' importante:
`ode_descent.m` e' l'**unica** definizione della dinamica in tutto HM2 (la usano
anche `rk4_zoh.m` e i replay di Task 2), quindi non c'e' rischio che Task 1 e
Task 2 stiano confrontando modelli diversi.

---

## `diagnostics` (righe 317-349)

### Tempi di switch (righe 333-340)

```matlab
thr  = 0.5 * d.Tmax;
i_dn = find(Tm(1:end-1) >= thr & Tm(2:end) <  thr, 1, 'first');
i_up = find(Tm(1:end-1) <  thr & Tm(2:end) >= thr, 1, 'last');
```

- Riga 333: soglia **euristica** a meta' della spinta massima.
- Riga 335: `i_dn` = **primo** attraversamento in discesa -> fine della prima
  bruciata.
- Riga 336: `i_up` = **ultimo** attraversamento in salita -> inizio della seconda.
- Riga 337: `cross(i)` interpola **linearmente** `|T|` fra i due nodi per stimare
  l'istante esatto. E' coerente col modello: la trapezoidale *dice* che il
  controllo e' piecewise-linear.
- Riga 340: `coast = t_sw2 - t_sw1`.

**Perche' esiste un coast arc -- la teoria dietro il numero.** Applicando il
principio del minimo di Pontryagin al problema continuo, l'Hamiltoniano e'

    H = lambda_x*vx + lambda_y*vy + lambda_vx*(Tx/m)
        + lambda_vy*(Ty/m - 1) - lambda_m*Vc*|T|

Scrivendo `T = |T| * u_hat` con `u_hat` versore, la parte che dipende dal
controllo diventa

    H_u = |T| * [ (lambda_v . u_hat)/m - Vc*lambda_m ]

Due conseguenze:

1. La **direzione** ottima allinea `u_hat` con `-lambda_v` (il *primer vector* di
   Lawden): la spinta punta lungo il costato di velocita', cambiato di segno.
2. `H` e' **affine in `|T|`**: il coefficiente e' la **switching function**
   `S = -|lambda_v|/m - Vc*lambda_m`. Una funzione affine su un intervallo
   `[Tmin, Tmax]` ha il minimo **su un estremo**. Quindi:

       |T| = Tmax  se S < 0
       |T| = Tmin  se S > 0
       arco singolare  se S == 0 su un intervallo

Da qui la struttura **bang-off-bang** (`max - coast - max`), esattamente quella
osservata numericamente: `t_sw1 ~ 14.0 s`, `t_sw2 ~ 33.1 s` per il caso nominale
(README). **Lo script non calcola mai i costati** -- sono le condizioni KKT del
NLP a fare da controparte discreta al PMP, e i moltiplicatori di `fmincon` sono,
a meno della scalatura per `dt`, i costati discretizzati. E' il ponte fra HM2
(metodi diretti) e HM1 (metodi indiretti), e vale la pena saperlo dire.

Nota di onesta' sui limiti del codice: `find(..., 1, 'first')` e
`find(..., 1, 'last')` **assumono** proprio la struttura max-coast-max. Se il
profilo avesse piu' archi, `i_dn` e `i_up` catturerebbero il coast piu' esterno,
mascherando la struttura interna; se `|T|` non attraversasse mai `0.5*Tmax` con
quel pattern, sarebbero **vuoti** e le diagnostiche degenererebbero
silenziosamente (nessun guard, nessun controllo `isempty`).

### Margine di glide slope (righe 342-344)

- Riga 342: `ok = sol.y > 1` -- maschera i nodi sotto 1 metro di quota. Il motivo e'
  scritto nel commento alle righe 325-326: al pad `atan2(|x|, y)` e' `0/0`, e i
  residui numerici entro tolleranza (millimetri) producono angoli **arbitrari**
  fra 0 e 90 gradi. Senza la maschera il grafico e il margine sarebbero rumore.
- Riga 343-344: `gs_margin = theta_max_deg - max(theta_deg)`. Positivo = corridoio
  mai violato **sui nodi mascherati**.

### Attivita' KKT (righe 346-348)

```matlab
lam = sol.lambda.ineqnonlin;
dg.n_thr_active = sum(lam(N+1:2*N) > 1e-6);
dg.max_gs_mult  = max(lam(2*N+1:4*N));
```

- Riga 347: conta i nodi in cui il **moltiplicatore** del bound superiore di
  spinta e' strettamente positivo. La lettura e' la **complementarita'**: nel KKT
  vale `lambda_i * g_i(z) = 0` con `lambda_i >= 0`, quindi `lambda_i > 0`
  implica `g_i = 0`, cioe' vincolo **attivo**. Il codice deduce l'attivita' dal
  moltiplicatore, non dal residuo -- che e' la lettura corretta (piu' robusta
  numericamente, e assume complementarita' stretta).
- Riga 348: il massimo moltiplicatore fra tutti quelli di glide slope. La README
  riporta che e' **numericamente zero**: il corridoio non e' mai attivo, la
  soluzione ci sta dentro **spontaneamente**. Cioe': la traiettoria min-fuel non
  e' vincolata dal glide slope, in questa geometria. Il solo vincolo di percorso
  attivo e' `|T| <= Tmax`.
- Il valore `1e-6` come soglia di positivita' e' arbitrario e non commentato.
- Ricordare (dalla riga 159) che questi `lambda` sono **adimensionali** e non
  vanno letti come grandezze fisiche.

---

## `fwd_integrate_pwl` (righe 351-377) e `node_err` (righe 379-387)

```matlab
u_fcn = @(tt) [
  sol.Tx(k) + (sol.Tx(k+1)-sol.Tx(k))*(tt-t_k)/(t_kp-t_k);
  sol.Ty(k) + (sol.Ty(k+1)-sol.Ty(k))*(tt-t_k)/(t_kp-t_k)];
[~, Y] = ode45(rhs_t, [t_k, t_kp], X(k,:).', opts);
X(k+1,:) = Y(end,:);
```

- Riga 365: parte dalle **condizioni iniziali esatte**, non dal nodo 1 del NLP
  (identiche comunque, essendo la BC un'uguaglianza).
- Riga 366: tolleranze `RelTol = 1e-10`, `AbsTol = 1e-12` -- `ode45` e' qui il
  "verita' di riferimento", quindi deve essere praticamente esatto (convenzione
  della repo per il lavoro shooting/indiretto).
- Righe 369-371: ricostruisce il controllo **piecewise-linear** fra i nodi `k` e
  `k+1`. Questa e' la scelta consistente con la trapezoidale (vedi sopra).
- Righe 373-374: **la propagazione riparte da `X(k,:)`, non dal nodo NLP
  `sol.x(k)`** (la condizione iniziale passata a `ode45` alla riga 373 e'
  `X(k,:).'`, cioe' il risultato del passo precedente, memorizzato alla riga 374).
  Cioe' e' un replay **open-loop** genuino: gli errori si **accumulano**
  intervallo dopo intervallo. Non e' un errore locale per intervallo -- e' la
  domanda "se accendessi i motori esattamente cosi', dove finirei davvero?". La
  risposta (README): a 4.3 m dal pad, a 0.11 m/s.
- Riga 386 (`node_err`): norma euclidea della differenza su **posizione e
  velocita' soltanto** (`X(:,1:4)`), la massa e' esclusa. Il codice non dice
  perche'; ragionevolmente perche' la massa e' l'obiettivo e non un errore di
  tracking, e perche' la sua scala e' diversa. Ma va detto: la metrica **non
  cattura** un eventuale errore sulla massa consumata.

---

## `plot_results` (righe 389-460)

Quattro figure, con le opzioni utili da saper spiegare:

- Righe 404-418, **traiettoria**: `axis equal` (obbligatorio, altrimenti il cono
  di glide slope apparirebbe deformato), le due rette tratteggiate del corridoio
  costruite come `xx = tan(theta_max)*yy` (riga 407), il marker triangolare nero
  sul pad.
- Righe 421-430, **modulo della spinta**: `yline(Tmax)` per rendere visibile che
  la soluzione **satura** il bound. E' il grafico che mostra il bang-off-bang.
- Righe 433-441, **massa**: monotona decrescente con un plateau durante il coast
  (dove `m_dot = -Vc*|T| = 0`).
- Righe 446-459, **angolo di glide slope**: stessa maschera `y > 1` di
  `diagnostics` (righe 444-445 lo commentano), e alla riga 456 un `ylim` forzato
  per tenere la linea di `theta_max` dentro l'inquadratura anche quando la curva
  le sta molto sotto -- altrimenti la scala automatica la butterebbe fuori e non si
  vedrebbe il margine.

---

## Possibili domande d'esame

**D: Che cos'e' esattamente un "defect constraint" e perche' e' un vincolo di
uguaglianza?**
R: E' il residuo dell'equazione di stato discretizzata:
`zeta_k = x_{k+1} - x_k - (dt/2)*(f_k + f_{k+1})`. Nella collocazione diretta i
valori dello stato ai nodi sono **variabili di decisione indipendenti**: nulla, di
per se', le lega alla dinamica. Il defect e' proprio il vincolo che le lega. Se
`zeta_k = 0` per ogni `k`, i nodi giacciono su una traiettoria del sistema
discreto; se non lo fosse, sarebbero numeri senza significato fisico. Deve essere
un'**uguaglianza** (non una disuguaglianza) perche' la dinamica non e' negoziabile:
non esiste un "quasi soddisfatta in un verso". Sono `5 x (N-1) = 245` uguaglianze
nonlineari con `N = 50`.

**D: Perche' collocazione diretta e non shooting (come in HM1)?**
R: Nel single shooting le uniche incognite sono i parametri iniziali, e lo stato
si ottiene propagando: il problema e' piccolo ma **estremamente mal condizionato**
(la sensitivita' dello stato finale rispetto alle condizioni iniziali cresce
esponenzialmente col tempo di volo, ed e' la ragione per cui HM1 ha bisogno di
continuazione). Nella collocazione diretta lo stato e' **tutto** incognito, quindi
il problema e' grande (350 variabili) ma **sparso e ben condizionato**: ogni
vincolo tocca solo due nodi adiacenti, quindi un errore locale non si propaga
lungo tutta la traiettoria. In piu' i vincoli di percorso (glide slope, bound di
spinta) si impongono **direttamente** ai nodi, mentre in un metodo indiretto
richiederebbero moltiplicatori di percorso e una struttura di archi da indovinare
a priori. Il prezzo e' l'accuratezza: `O(dt^2)` invece della precisione
dell'integratore.

**D: Il tuo problema di ottimizzazione e' convesso?**
R: **No**, anche se molti dei suoi pezzi lo sono. Il costo e' lineare
(`-m_N`); il glide slope e' scritto come coppia di disuguaglianze **lineari**; il
bound di spinta e' `|T| <= Tmax`, che con `Tmin = 0` e' un **disco** convesso (con
`Tmin > 0` sarebbe una corona, non convessa -- ed e' il caso interessante).
La non convessita' viene tutta dai **defect**, che contengono `Tx/m` e `Ty/m`
(bilineari) e `|T|` nell'equazione di massa. Quindi `fmincon/sqp` restituisce un
**minimo locale** senza certificato globale, e la risposta dipende (in linea di
principio) dalla guess. E' proprio questo che motiva la variante GFOLD del Task 2:
col cambio di variabili log-massa `z = ln m`, `u = T/m`, la dinamica diventa
**esattamente LTI** e con lo slack `Gamma >= |u|` l'intero problema collassa in
**un solo SOCP** convesso, con ottimo globale certificato.

**D: Il vincolo `Tmin <= |T| <= Tmax` non e' non convesso? Come lo hai gestito?**
R: Nel caso generale si', e' una **corona circolare** (l'insieme
`{|T| >= Tmin}` e' il complemento di un disco). Nel mio problema pero' la traccia
da' `Tmin = 0`, quindi il vincolo inferiore `-|T| <= 0` e' **vacuo** e l'insieme
si riduce al disco convesso. Il codice lo scrive comunque in forma generica
(riga 297), quindi con `Tmin > 0` girerebbe ancora -- ma `fmincon` lo tratterebbe
come un vincolo nonlineare qualsiasi, cioe' cercherebbe un minimo locale su un
insieme non convesso, senza garanzie. La cura corretta e' la **lossless
convexification**: si rilassa a `Tmin <= Gamma <= Tmax` con `|T| <= Gamma`, e si
dimostra che all'ottimo il rilassamento e' **attivo** (`Gamma = |T|`), quindi la
soluzione del problema convesso e' anche soluzione dell'originale.

**D: `fmincon` ti dice `exitflag = 0` su due run su tre. La tua soluzione vale
qualcosa?**
R: `exitflag = 0` significa **MaxIterations raggiunto**, non "fallito". La
**feasibility** e' comunque sotto `1e-6` (i vincoli sono soddisfatti); a stallare
e' la **first-order optimality**, ferma a `1e-3`-`1e-4` non-dim. Quindi ho una
traiettoria valida, ma non ho la prova che sia stazionaria. La causa radice e'
identificata: la norma `|T| = sqrt(Tx^2 + Ty^2)` **non e' differenziabile in
`T = 0`**, e la soluzione ottima ci passa attraverso per tutto il coast arc. Con
gradienti a differenze finite (che e' quello che `fmincon` sta usando, perche' lo
script non fornisce Jacobiane analitiche) il gradiente li' e' spazzatura. Rimedi,
in ordine di eleganza: (1) lo slack `Gamma >= |T|` della lossless convexification,
che **toglie la norma dalla dinamica**; (2) gradienti analitici con norma
regolarizzata `sqrt(|T|^2 + eps)`; (3) continuazione, warm-startando ogni run
dello sweep dalla soluzione nominale.

**D: Perche' il vettore delle incognite e' impilato nodo per nodo e non
variabile per variabile?**
R: Perche' con lo stacking node-major (`[stato_1; u_1; stato_2; u_2; ...]`) il
defect dell'intervallo `k` tocca **14 componenti contigue** di `z` -- quelle dei
nodi `k` e `k+1`. La Jacobiana dei vincoli risulta **a banda** (block-bidiagonale),
ed e' esattamente la struttura che un solver sparso (`interior-point` con
`ConstraintGradient`, o IPOPT) sfrutta. Con l'ordinamento alternativo (tutte le
`x`, poi tutte le `y`, ...) le stesse 14 componenti sarebbero disperse per tutto
il vettore. In piu', essendo MATLAB column-major, `reshape(z, 7, N)` recupera i
nodi in colonna **senza copie**. Onesta': lo script **non** sfrutta questa
sparsita' -- non dichiara sparsity pattern ne' Jacobiane analitiche, e `fmincon`
calcola tutto per differenze finite (351 valutazioni di `nonlcon` per Jacobiana).
E' il principale limite di performance del Task 1.

**D: Come costruisci la guess iniziale, e perche' funziona?**
R: **Interpolazione lineare** (una corda, non una propagazione) fra lo stato
iniziale e lo stato finale per `x, y, vx, vy`; massa lineare da `m0` a `0.7*m0`
(30% di frazione di carburante indovinata); controllo costante di **hover**,
`T = (0, m0*g)`, cioe' `(0, 1)` in adimensionale. La proprieta' che la rende buona
e' che **soddisfa gia' tutti i vincoli tranne i defect**: le BC per costruzione,
il bound di spinta (`|T| = 1 < 3.57`), e anche il glide slope, perche' lungo la
corda il rapporto `x/y` resta **costante** a `atan(1000/3000) = 18.4 deg`, ben
dentro il cono di 60 gradi. Quindi l'unica infeasibility che SQP deve sanare e'
quella dinamica -- che e' proprio il caso in cui SQP e' forte (gestisce bene punti
di partenza infeasible).

**D: Cosa ti aspetteresti di guadagnare passando a Hermite-Simpson?**
R: Hermite-Simpson e' `O(dt^4)` invece di `O(dt^2)`, a costo di aggiungere un
punto di collocazione a meta' intervallo (lo stato al midpoint e' ricostruito con
un'interpolante cubica di Hermite, il controllo al midpoint diventa una nuova
variabile). A parita' di `N` l'errore di replay crollerebbe di ordini di
grandezza; equivalentemente, a parita' di accuratezza servirebbero molti meno
nodi, quindi un NLP piu' piccolo. Nel mio caso specifico pero' il collo di
bottiglia **non e' l'ordine dello schema**: e' la non differenziabilita' di `|T|`
in `T = 0` combinata con le differenze finite, che nessun aumento di ordine
risolve. Prima si sistema il gradiente, poi ha senso alzare l'ordine.
