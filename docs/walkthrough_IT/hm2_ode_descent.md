# HM2_powered_descent/ode_descent.m

## Ruolo del file nel progetto

Questo file contiene la **dinamica continua** del problema di powered descent di
HM2, scritta in forma adimensionale. E' il "cuore fisico" di tutta l'homework:
ogni trascrizione diretta (trapezoidale in Task 1, ZOH+RK4, LTV+SCvx in Task 2)
alla fine chiama questa funzione per sapere come si muove il razzo. Il modello e'
un **punto materiale 2D su Terra piatta, senza aerodinamica**: niente drag,
niente portanza, gravita' costante. La massa e' uno stato, non un parametro,
perche' il razzo brucia propellente e la sua accelerazione dipende da m(t).

La parametrizzazione del controllo e' **il vettore di spinta T = [Tx; Ty]**: le
componenti di T sono le variabili di decisione dell'NLP, e sull'intervallo ZOH
e' T (non l'accelerazione) a essere tenuto costante. Questo e' il punto che
distingue questo file dal gemello `ode_descent_uacc.m`, che invece tiene
costante l'accelerazione u = T/m.

Chi lo chiama:
- `main_task1.m` riga 314, attraverso il wrapper `dyn_rhs`, per costruire i
  difetti trapezoidali;
- `main_task2.m` riga 773 (`trap_nonlcon`, difetti trapezoidali del baseline);
- `main_task2.m` riga 295 (`ltv_aug_rhs`), dove serve a costruire l'offset
  affine `c_off = f - A*x_ref - B*u_ref` della linearizzazione di Appendice A;
- `main_task2.m` riga 852 (`fwd_integrate`), come RHS del replay `ode45` che
  misura la fedelta' della trascrizione;
- `rk4_zoh.m` righe 17-20, come RHS dei quattro stadi RK4;
- i test `tests/odeDescentTest.m`, `tests/rk4ZohTest.m`,
  `tests/descentDynamicsPerformanceTest.m`.

Il file e' lungo 16 righe di cui 2 di codice vero: la densita' sta tutta nella
derivazione delle due righe.

---

## `ode_descent` (righe 1-16)

```matlab
function dx = ode_descent(x, u, Vc)
...
Tmag = sqrt(u(1)^2 + u(2)^2);
dx = [ x(3); x(4); u(1)/x(5); u(2)/x(5) - 1; -Vc * Tmag ];
```

- **Riga 1**: firma. Tre soli argomenti: lo stato `x` (5x1), il controllo `u`
  (2x1) e lo scalare `Vc`. Notare **l'assenza dell'argomento tempo `t`**: la
  dinamica e' *autonoma*. Non lo e' per caso: con il controllo tenuto costante
  sull'intervallo (ZOH) il tempo non compare mai esplicitamente nel RHS, e questo
  e' esattamente cio' che permette a `rk4_zoh.m` di passare lo stesso `u` a tutti
  e quattro gli stadi senza interpolare nulla. Dove serve un RHS con il tempo
  (il replay `ode45` in `fwd_integrate`) e' il chiamante a costruire il wrapper
  `@(tt, x) ode_descent(x, u_fcn(tt), d.Vc)`.

- **Righe 2-11**: la docstring. Dichiara il layout dello stato
  `x = [x; y; vx; vy; m]`, il controllo `u = [Tx; Ty]`, e la convenzione
  chiave: **"Gravity = -1"**. Adimensionalizzato, il modulo della gravita' vale
  esattamente 1. Righe 10-11 dichiarano esplicitamente la scelta di design di
  non mettere il blocco `arguments`.

- **Riga 13**: `Tmag = sqrt(u(1)^2 + u(2)^2)` -- il modulo del vettore spinta,
  ||T||. Scritto a mano invece che con `norm(u)` per velocita' (siamo in un hot
  loop; `norm` fa dispatch e controlli). E' l'unica non-linearita' *nel
  controllo* dell'intero RHS, ed e' anche il punto dolente numerico: `sqrt` di
  una somma di quadrati e' **non differenziabile in T = 0**, cioe' esattamente
  dove si trova l'arco di coast della soluzione bang-off-bang. Il README di HM2
  attribuisce a questo kink lo stallo del criterio di ottimalita' di `fmincon`
  a 1e-3..1e-4 sui run di sensitivita'. Da notare che **qui la norma non e'
  regolarizzata**: e' la funzione `jacobians` in `main_task2.m` (riga 266) ad
  aggiungere `+1e-6` sotto radice per tenere finita df/du in T = 0. Le due
  funzioni quindi non sono perfettamente coerenti fra loro -- il RHS e' il vero
  modello, lo Jacobiano e' una versione ammorbidita.

- **Riga 14**: le cinque righe della dinamica, in ordine.

### Derivazione della dinamica

Il modello dimensionale (punto materiale, Terra piatta, no drag) e':

    x_dot  = vx
    y_dot  = vy
    vx_dot = Tx/m
    vy_dot = Ty/m - g
    m_dot  = -||T|| / c        con c = Isp * g0   (velocita' efficace di scarico)

La quinta riga e' la legge di consumo del razzo: la spinta e'
||T|| = m_dot_propellente * c, quindi il flusso di massa vale ||T||/c. Nel codice
non compaiono ne' Isp ne' g0: sono gia' collassati dentro `Vc`. In
`main_task2.m` riga 31 si ha `data.c = data.Isp * data.g0` e riga 226
`dnd.Vc = ref.V / d.c`.

**Adimensionalizzazione** (funzione `nondim`, `main_task2.m` righe 213-227):

    L_ref = y0 = 3000 m
    g_ref = g  = 9.81 m/s^2
    t_ref = sqrt(L_ref/g)        -> 17.49 s
    V_ref = sqrt(g*L_ref)        -> 171.55 m/s
    m_ref = m0 = 2000 kg
    T_ref = m_ref*g = 19620 N

Sostituendo x = L_ref*x', t = t_ref*t', v = V_ref*v', m = m_ref*m', T = T_ref*T'
e usando le identita' t_ref*V_ref = L_ref e t_ref/V_ref = 1/g:

- posizione: dx'/dt' = (t_ref*V_ref/L_ref) * vx' = **vx'** (coefficiente 1);
- velocita': dvx'/dt' = (t_ref/V_ref)*(T_ref/m_ref) * Tx'/m' = (1/g)*(g) * Tx'/m'
  = **Tx'/m'**;
- gravita': (t_ref/V_ref)*g = 1, da cui il **-1** costante nella quarta riga;
- massa: dm'/dt' = -(t_ref*T_ref)/(m_ref*c) * ||T'|| = -(t_ref*g/c)*||T'||
  = -(V_ref/c)*||T'|| = **-Vc * ||T'||**.

Quindi `Vc = V_ref/c` (il "numero di Tsiolkovsky") e' l'**unico** parametro
fisico rimasto nel RHS: tutta la fisica del propulsore e' compressa in uno
scalare. Con i dati di HM2 vale Vc = 171.55/2206.5 = 0.0777 (il valore usato nei
test come costante). Il messaggio implicito e' potente: due razzi con Isp e
scale diverse ma stesso Vc hanno **la stessa traiettoria ottima adimensionale**.

### Cosa NON c'e' nella dinamica

- **Nessun drag.** Il RHS non ha termini in rho, S, Cd, ne' V^2. Il modello e'
  valido solo perche' HM2 lo dichiara esplicitamente (2D point-mass, no
  aerodynamics). In una discesa vera da 3 km il drag non e' trascurabile; e'
  un'ipotesi dell'assegnazione, non una verita' fisica.
- **Gravita' costante**, non inverso del quadrato: su 3 km di quota la
  variazione di g e' ~1e-3 relativo, quindi l'approssimazione e' innocua.
- **Nessun vincolo di puntamento della spinta** e nessuna dinamica di assetto:
  T puo' ruotare istantaneamente. La direzione della spinta e' un controllo
  libero.
- **Nessuna guardia su m > 0.** Se un'iterata di `fmincon` porta `x(5)` a zero
  o sotto, le righe 3 e 4 esplodono o cambiano segno silenziosamente. La
  protezione e' delegata al chiamante: `box_bounds` in `main_task2.m` riga 691
  impone `lb = 1e-3` sulla massa.

> **Possibile domanda d'esame** -- La dinamica e' affine nel controllo?
> *Risposta:* No, non del tutto. Le righe di posizione non dipendono da u; le
> righe di velocita' sono *lineari* in u a stato fissato (u1/m, u2/m), ma il
> coefficiente 1/m dipende dallo stato, quindi il prodotto e' **bilineare** in
> (u, 1/m) e non convesso. La riga della massa contiene ||u||, che e' convessa ma
> non affine. Percio' il vincolo di difetto x_{k+1} = F(x_k, u_k) resta
> un'uguaglianza non convessa. E' proprio questo che costringe le varianti (b) e
> (c) a linearizzare e iterare (SCvx), e che la variante (d) rimuove con il
> cambio di variabili z = ln m, u = T/m.

---

## Perche' due RHS diversi (`ode_descent` vs `ode_descent_uacc`)

Il confronto vale la pena farlo qui perche' e' la domanda piu' probabile
all'orale. Le due funzioni descrivono **la stessa fisica**, ma differiscono per
*quale grandezza viene tenuta costante sull'intervallo ZOH*:

| | `ode_descent` | `ode_descent_uacc` |
|---|---|---|
| Controllo tenuto costante | spinta T = [Tx; Ty] | accelerazione u = T/m |
| vx_dot | Tx/m (varia: m cala) | ux (costante) |
| m_dot | -Vc*norm(T) (costante sull'intervallo) | -Vc*m*norm(u) |
| d(ln m)/dt | -Vc*norm(T)/m (varia) | -Vc*norm(u) (costante) |

Con T tenuta costante, l'**accelerazione cresce** durante l'intervallo perche' la
massa cala: T/m(t) e' crescente. Con u tenuta costante, e' la **spinta a
crescere/calare** seguendo la massa, T(t) = m(t)*u. Sono due ZOH diversi, quindi
due plant discreti diversi: non e' lecito confrontare i due set di controlli
nodo per nodo pretendendo che descrivano la stessa traiettoria. Infatti
`main_task2.m` usa due funzioni di replay distinte (`fwd_integrate` riga 852 con
`ode_descent`, `fwd_integrate_uacc` riga 1298 con `ode_descent_uacc`) proprio
per non barare nella validazione.

Osservazione (mia, non nel codice): con u costante la traiettoria dell'intervallo
ha **soluzione analitica banale** (doppio integratore + log-massa lineare),
mentre con T costante la velocita' segue una legge logaritmica di tipo
Tsiolkovsky, vx(t) = vx_k - (Tx/(Vc*||T||)) * ln(m(t)/m_k), che degenera in 0/0
quando ||T|| -> 0 (arco di coast). Questo e' un buon motivo pratico per
integrare `ode_descent` con RK4 invece che con la sua forma chiusa: RK4 e'
uniformemente valido anche a spinta nulla.

---

## Nessun blocco `arguments`: e' voluto

Righe 10-11 lo dichiarano: *"No arguments validation by design: hot-loop RHS
inside ode45/fmincon; validate at the call site."*

Il conto: nella variante (a) `zoh_nonlcon` chiama `rk4_zoh` per ciascuno dei
N-1 = 49 intervalli, ognuno con `n_sub = 2` sotto-passi da 4 stadi -> 49*2*4 =
392 valutazioni di `ode_descent` per **una sola** valutazione del vincolo. Con
`fmincon` in modalita' SQP e gradienti alle differenze finite su 7N = 350
variabili, ogni iterazione maggiore ricostruisce lo Jacobiano dei vincoli con
~351 valutazioni di `nonlcon`, cioe' circa **1.4e5 chiamate a `ode_descent` per
iterazione**. Un blocco `arguments` (che fa validazione di classe e dimensione a
ogni chiamata) costerebbe piu' del RHS stesso -- il benchmark
`descentDynamicsPerformanceTest` (riga 39) stima una singola valutazione intorno
ai **60 ns**, tanto bassa da dover essere misurata a batch di 1000.

La validazione vive invece nelle funzioni di frontiera chiamate una volta per
run: `nondim`, `dim_sol`, `compute_ltv_zoh`, `fwd_integrate`, `lti_zoh`. E' la
convenzione dichiarata anche in `CLAUDE.md` ("Hot-loop functions carry no
arguments validation by design").

> **Possibile domanda d'esame** -- Se il RHS non valida gli input, come si evita
> di passargli uno stato con la massa negativa?
> *Risposta:* Non lo si evita nel RHS, lo si evita a monte. I bound di scatola
> dell'NLP (`box_bounds`, `main_task2.m` riga 691) impongono m in [1e-3, m0] su
> ogni nodo, quindi `fmincon` non valuta mai il vincolo fuori da quella scatola
> (con l'algoritmo SQP i bound sono rispettati a ogni iterata). Il rischio
> residuo e' negli *stadi interni* di RK4, dove uno stato intermedio potrebbe in
> teoria uscire dai bound: e' un'ipotesi implicita che il codice non verifica.

---

## Possibili domande d'esame

**D: Perche' la gravita' compare come un "-1" secco e non come 9.81?**
R: Perche' il problema e' adimensionalizzato scegliendo a_ref = g. Con le scale
L_ref = y0, t_ref = sqrt(L_ref/g), V_ref = sqrt(g*L_ref) si ottiene
(t_ref/V_ref)*g = 1, quindi il termine gravitazionale nella quarta riga vale
esattamente -1. Il vantaggio e' che tutti gli stati sono O(1): posizioni e
velocita' iniziali diventano (0.33, 1) e (1.75, -1.17), la massa iniziale 1. Con
uno scaling cosi' l'NLP e' molto meglio condizionato e le tolleranze di `fmincon`
(1e-5 su ottimalita', 1e-6 sui vincoli) hanno lo stesso significato su tutte le
componenti.

**D: Cosa rappresenta fisicamente `Vc` e perche' e' l'unico parametro rimasto?**
R: Vc = V_ref/c = sqrt(g*L_ref)/(Isp*g0). E' il rapporto fra la velocita'
caratteristica del problema e la velocita' efficace di scarico del motore -- un
"numero di Tsiolkovsky". Regola quanto costa in massa una data manovra: dalla
riga 5 nondimensionale, la massa consumata su un intervallo e' Vc volte
l'impulso adimensionale. Con i dati di HM2, Vc = 0.0777. Se Vc -> 0 (motore
infinitamente efficiente) il razzo non consuma e il problema di minimo consumo
degenera.

**D: Perche' `m_dot` dipende solo dal modulo della spinta e non dalla sua
direzione?**
R: Perche' la portata di propellente e' fissata dal regime del motore, e la
spinta prodotta e' ||T|| = m_dot_prop * c indipendentemente da dove la si punta.
Ruotare il vettore di spinta (nel modello, gimbal ideale istantaneo) non costa
massa. Il test `testMassFlowDependsOnlyOnThrustMagnitude` (`odeDescentTest.m`
riga 39) verifica proprio questa invarianza confrontando T = [1;0] e T = [0;-1].

**D: Il modulo ||T|| e' non differenziabile in T = 0. Che problema crea e come
lo si aggira?**
R: La soluzione ottima e' bang-off-bang (max-coast-max): l'arco di coast si
piazza *esattamente* sul kink T = 0. `fmincon` con gradienti alle differenze
finite vede una derivata mal definita li' e il criterio di ottimalita' di
prim'ordine si blocca a 1e-3..1e-4 invece di scendere a 1e-6 (documentato nel
README di HM2). Il codice mitiga solo lo Jacobiano (`jacobians`, riga 266,
regolarizza con `sqrt(Tx^2 + Ty^2 + 1e-6)`), non il RHS. Le cure vere sono: (i)
il rilassamento lossless con slack sigma >= ||T||, che toglie la norma dalla
dinamica e la trasforma in un vincolo conico (la strada della variante GFOLD);
(ii) gradienti analitici con norma smussata; (iii) continuation.

**D: In `ltv_aug_rhs` (riga 295-296) `ode_descent` e' usata per calcolare
`c_off = f - A*x_ref - B*u_ref`. Perche' serve questo termine?**
R: Perche' la linearizzazione di una dinamica non lineare attorno a una
traiettoria di riferimento non e' lineare ma **affine**: f(x,u) ~ f(x_r,u_r) +
A*(x-x_r) + B*(u-u_r) = A*x + B*u + [f(x_r,u_r) - A*x_r - B*u_r]. Il termine in
parentesi e' l'offset c_off, che non e' zero perche' il riferimento non e' un
punto di equilibrio (il razzo sta accelerando). Se lo si dimenticasse, il
sistema linearizzato "perderebbe" gravita' e la parte di dinamica non catturata
dai Jacobiani, e SCvx non convergerebbe alla soluzione non lineare.

**D: Questa dinamica e' convessa? Se no, cosa la rende non convessa?**
R: No. Due sorgenti di non convessita': (i) il termine u/m, bilineare fra
controllo e reciproco della massa (il vincolo di difetto e' un'uguaglianza non
lineare, quindi non convessa a prescindere dal segno); (ii) il vincolo di
**lower bound** sulla spinta ||T|| >= Tmin, che in HM2 e' inattivo perche'
Tmin = 0 (`main_task2.m` riga 32) ma che in un caso realistico e' un vincolo di
tipo "esterno a una palla", palesemente non convesso. La riga della massa
-Vc*||T|| e' invece convessa in u, ma comparendo dentro un'uguaglianza non
aiuta. Il cambio di variabili di Acikmese-Blackmore (z = ln m, u = T/m,
sigma >= ||u||) elimina la prima e convessifica la seconda.
