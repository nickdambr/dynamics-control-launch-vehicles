# HM2_powered_descent/proto_gfold_logmass.m

## Ruolo del file nel progetto

Questo file e' il **prototipo standalone della trascrizione GFOLD log-mass**, cioe'
la quarta (e ultima) trascrizione del Task 2 di HM2. E' uno script: si lancia da
solo, risolve l'intero problema di powered descent in pochi secondi e stampa a
schermo un blocco di validazione. L'header (righe 1-9) dichiara esplicitamente che
il codice e' stato poi promosso dentro `main_task2.m` come "variant (d)"
(`solve_gfold_scvx` / `solve_gfold_socp`), e che i due kernel numerici sono stati
estratti in file autonomi -- `lti_zoh.m` e `ode_descent_uacc.m` -- condivisi fra il
prototipo e `main_task2.m`. Il prototipo e' quindi **codice vivo ma ridondante**:
serve da demo minimale, e l'header stesso dice "Safe to delete if the integrated
version is enough". Onesta': nessun altro file lo chiama, e i test
(`tests/gfoldLogMassTest.m`) testano i due kernel condivisi, *non* questo script.

Il problema risolto e' sempre quello della traccia HM2: punto materiale 2D, massa
variabile, atterraggio pinpoint e soft da (1000, 3000) m con velocita' iniziale
(300, -200) m/s, `tf = 38 s` fisso, `0 <= |T| <= 70 kN`, corridoio glide-slope a
60 gradi, `Isp = 225 s`, costo = massimizzare la massa finale. Le altre tre
trascrizioni (trapezoidale, ZOH+RK4 nonlineare, LTV+SCvx) attaccano il problema
come un **NLP non convesso**. Questa lo attacca cambiando variabili in modo che il
problema diventi (quasi) un **programma conico del second'ordine**.

L'idea, dovuta ad Acikmese/Ploen/Blackmore, e' in due mosse. Prima si sostituisce
il controllo: invece del vettore spinta `T` si usa l'**accelerazione comandata**
`u = T/m`, e invece della massa si usa il suo logaritmo `z = ln(m)`. Con questa
sostituzione la dinamica diventa **esattamente lineare tempo-invariante** in
`(r, v, z)` con ingresso `(u, sigma)` -- non linearizzata: *esattamente* lineare.
Seconda mossa: si introduce uno **slack scalare** `sigma` con `|u| <= sigma`
(vincolo conico, convesso) che rimpiazza la norma `|u|` dentro l'equazione della
massa. Il rilassamento e' **lossless**: all'ottimo `sigma = |u|` e quindi la
soluzione del problema convesso rilassato risolve anche l'originale.

Il file usa questi due mattoni condivisi:

- `lti_zoh.m` -- costruisce le matrici ZOH discrete `(Abar, Bbar, cbar)` con un
  singolo `expm` (trucco di van Loan). Poiche' il sistema e' LTI, le matrici sono
  **costanti lungo tutta la griglia** e si calcolano **una volta sola** (riga 181).
- `ode_descent_uacc.m` -- il right-hand side nonlineare "vero" con l'accelerazione
  `u` tenuta costante sull'intervallo; serve solo come *ground truth* per il replay
  ode45 (righe 241-254), che misura la fedelta' della trascrizione.

Cio' che il codice **non** fa, e va detto subito: non risolve un unico SOCP. Il
paper GFOLD linearizza il bound di spinta superiore attorno a un profilo di massa
a-priori *fissato* e risolve **un solo** problema convesso, con certificato di
ottimalita' globale. Qui invece la linearizzazione viene ri-centrata a ogni
iterazione dentro un loop SCvx (righe 179-239). Il README di HM2 lo riconosce: il
ticket "Lossless convexification -> single SOCP" e' ancora aperto.

---

## Header e cambio di variabile (righe 1-27)

```matlab
% Change of variables (GFOLD, Acikmese/Blackmore):
%   state    xi = [x; y; vx; vy; z],  z = ln(m)
%   control  w  = [ux; uy; sigma],    u = T/m, sigma >= ||u||
```

- Righe 11-23: il commento dichiara il cuore teorico del file. Vale la pena
  derivarlo per intero, perche' e' *il* punto d'orale.

**Punto di partenza.** La dinamica non-dim in variabili originali (quella delle
altre tre trascrizioni, `ode_descent.m`) e':

    x'  = vx
    y'  = vy
    vx' = Tx/m
    vy' = Ty/m - 1
    m'  = -Vc * |T|

dove `|T| = sqrt(Tx^2 + Ty^2)`, la gravita' vale 1 in unita' non-dim (vedi la
scelta delle scale, righe 103-113) e `Vc = V_rif/c` con `c = Isp*g0` la velocita'
efficace di scarico. Due nonlinearita': la **divisione per m** nelle righe di
velocita' (bilineare in `(T, 1/m)`) e la **norma** nella riga di massa. Sono
entrambe *algebriche*, non dinamiche.

**Mossa 1 -- cambio di controllo.** Si definisce l'accelerazione comandata

    u = T/m      ->      T = m*u

La sostituzione e' invertibile finche' `m > 0`, quindi non si perde nulla. Le righe
di velocita' diventano

    vx' = ux
    vy' = uy - 1

cioe' **non contengono piu' la massa**. Questa e' esattamente la proprieta' che il
test `testAccelerationIsDirect` verifica.

**Mossa 2 -- cambio di stato.** Si definisce `z = ln(m)` (diffeomorfismo fra
`m > 0` e `z` reale). Allora

    z' = m'/m = -Vc*|T|/m = -Vc*|m*u|/m = -Vc*|u|

La massa **sparisce anche dalla riga di massa**. In variabili dimensionali lo
stesso conto da'

    d(ln m)/dt = -|u| / (Isp*g0) = -alpha*|u|,     alpha = 1/(Isp*g0)

e `Vc` non e' altro che `alpha` reso non-dim (`Vc = alpha * V_rif`, vedi riga 111:
`dnd.Vc = ref.V/d.c`).

**Mossa 3 -- slack.** Resta la norma `|u|`, che e' convessa ma **nonlineare**: una
*uguaglianza* che contiene una funzione convessa nonlineare e' un vincolo non
convesso. Si introduce lo scalare `sigma` e si scrive

    z' = -Vc*sigma          (lineare)
    |u| <= sigma            (cono del second'ordine, convesso)

Il sistema risultante e' affine (LTI + termine costante di gravita'):

    xi' = A*xi + B*w + c
    A: x<-vx, y<-vy   B: vx<-ux, vy<-uy, z<--Vc*sigma   c = [0;0;0;-1;0]

**Perche' e' LTI e non "linearizzato".** Non c'e' nessuno sviluppo di Taylor,
nessun termine trascurato, nessun punto di lavoro. `A`, `B`, `c` sono matrici
costanti che valgono per **qualunque** stato e qualunque controllo, non solo vicino
a una traiettoria di riferimento. Le variabili `(u, sigma)` sono controlli
*genuini* (l'inversa `T = m*u` e' esplicita), quindi la riscrittura e' un cambio di
coordinate esatto sullo spazio stato-controllo. Contrasto con la variante (c),
dove la linearizzazione LTV e' un'approssimazione al prim'ordine attorno alla
traiettoria di riferimento, valida solo dentro la trust region: li' l'errore di
linearizzazione esiste e il ratio test SCvx deve misurarlo.

- Righe 20-23: il commento e' preciso nel dire che **l'unica non convessita'
  residua** e' il bound superiore di spinta `|u| <= Tmax*exp(-z)`. Su questo torno
  alla sezione `solve_gfold_socp`.

> **Possibile domanda d'esame** -- Se il cambio di variabile e' esatto, dove e'
> finita la non convessita' del problema originale?
> *Risposta:* Si e' spostata dalla **dinamica** ai **vincoli**. Nel problema
> originale la non convessita' sta nelle equazioni di stato (`T/m` bilineare, `|T|`
> dentro un'uguaglianza). Dopo la sostituzione la dinamica e' convessa (affine), il
> costo e' lineare, glide-slope e quota sono lineari: resta solo
> `|u| <= Tmax*exp(-z)`, cioe' la regione **sotto** il grafico di una funzione
> convessa di `z`, che e' non convessa. Il guadagno e' netto: una singola
> disuguaglianza scalare non convessa al posto di 5N vincoli di uguaglianza non
> convessi.

---

## Dati del problema e non-dimensionalizzazione (righe 31-48)

- Righe 32-40: dati di Tabella 1 della traccia in SI. `Tmin = 0` (riga 38): la
  traccia non impone una spinta minima. Questo ha una conseguenza importante che il
  codice non commenta -- vedi sotto.
- Riga 42: `[ref, d] = nondim(d_si)` -- scale di riferimento e struct non-dim.
- Riga 43-44: `tf = 38/ref.t`, `dt = tf/(N-1)` con `N = 50`. Griglia uniforme,
  `N-1 = 49` intervalli ZOH.
- Righe 46-48: stampa diagnostica delle scale.

Con i numeri della riga 32-40 le scale valgono `L = 3000 m`, `V = sqrt(g*L) ~ 171.6
m/s`, `t = sqrt(L/g) ~ 17.49 s`, quindi `Vc ~ 0.0777` e `dt ~ 0.0444` -- sono
esattamente i due valori hard-coded come costanti nel test (`gfoldLogMassTest.m`,
righe 7-8).

> **Possibile domanda d'esame** -- `Tmin = 0`: il rilassamento lossless serve ancora
> a qualcosa?
> *Risposta:* Si', ma per un motivo diverso da quello del paper. Nel paper lo slack
> serve soprattutto a convessificare il bound **inferiore** `|T| >= Tmin`, che e' la
> non convessita' piu' cattiva (un lower bound su una norma non e' mai convesso).
> Qui `Tmin = 0` e quel vincolo non c'e'. Lo slack serve comunque perche'
> l'*uguaglianza* della massa `z' = -Vc*|u|` contiene una norma: rilassarla in
> `z' = -Vc*sigma` con `|u| <= sigma` e' cio' che rende il vincolo dinamico affine.
> Onesta': questo codice **non esercita la meta' piu' spettacolare** della lossless
> convexification.

---

## Controllo disponibilita' solver (righe 50-53)

- Righe 51-53: se YALMIP (`sdpvar`) o ECOS non sono sul path, lo script muore con
  un `error`. Nessun fallback: a differenza di `main_task2.m` (che secondo il README
  salta con grazia le varianti coniche se mancano i pacchetti), qui la dipendenza e'
  dura. E' accettabile per un prototipo.

---

## Chiamata al solver SCvx (righe 55-57)

- Riga 56: `max_iter = 20`, `tol = 1e-3`.
- Riga 57: `sol = solve_gfold_scvx(...)`. Il secondo output `hist` della funzione
  (righe 197-198) esiste ma **viene scartato**: nessuno lo usa.

---

## Validazione in variabili originali (righe 59-86)

```matlab
[~, X_rep] = fwd_integrate_uacc(sol, d);   % replay nonlineare u-ZOH
pos_err = norm(X_rep(end,1:2)) * ref.L;
cone_gap = max(abs(sol.sig(1:N-1) - ...
                vecnorm([sol.ux(1:N-1) sol.uy(1:N-1)],2,2)));
```

- Riga 60: `dim_sol` riporta la soluzione in SI.
- Riga 61: **replay**. La soluzione ottima (le accelerazioni `u_k`) viene rigiocata
  con `ode45` attraverso la dinamica **nonlineare** in `(x,y,vx,vy,m)`. E' il test
  di fedelta' della trascrizione: se il modello discreto mente, il replay atterra
  altrove.
- Righe 62-65: errori di touchdown in posizione [m] e velocita' [m/s], drift di
  massa finale modello-vs-replay [kg], errore nodo per nodo. Attenzione alla riga
  65: `node_e` confronta **solo le colonne 1:4** (posizione e velocita'), la massa
  e' esclusa -- esattamente come `node_err` in `main_task2.m` (righe 859-868). Il
  canale di massa e' coperto da `dmf` (riga 64), non dall'errore di nodo.
- Riga 69: margine glide-slope `tan(theta_mx)*y - |x| >= 0`.
- Riga 70: monotonia della massa (`diff(m) <= 1e-9`).
- **Riga 71**: `cone_gap` -- la riga piu' importante del file. Misura
  `max_k |sigma_k - |u_k||`, cioe' **quanto e' attivo il rilassamento lossless**.
  Se questo numero e' ~1e-9 o meno, allora `sigma = |u|` a tutti i nodi e la
  soluzione del SOCP rilassato **e' anche soluzione del problema originale non
  convesso**. Nota che l'indice si ferma a `N-1`: il nodo `N` non ha controllo (i
  controlli sono `N-1`, uno per intervallo) ed e' stato riempito con zeri alla riga
  173, quindi includerlo falserebbe il massimo.
- Righe 73-86: stampa. Riga 76: valori di riferimento delle altre trascrizioni,
  **hard-coded a mano** nel `fprintf` (non ricalcolati) -- vanno letti come commento,
  non come output.
- Righe 79-80: check spinta con tolleranza `+1 N`. Riga 81-84: glide-slope e quota
  con tolleranza `-1e-6` non-dim.

**Perche' il replay deve tornare a macchina.** Nelle varianti (a)/(b)/(c) il replay
ode45 ha un errore residuo (la trascrizione approssima). Qui no: la mappa discreta
e' la soluzione esatta della dinamica LTI con `w` costante, quindi *se* la dinamica
vera vista dal replay e' quella con `u` costante (ed e' proprio cio' che fa
`ode_descent_uacc`), l'unica discrepanza possibile fra modello e replay e'
`sigma - |u|`. In altre parole: `dmf` (riga 64) e `cone_gap` (riga 71) misurano
**la stessa cosa** da due angoli diversi. Il README riporta per questa variante un
errore di nodo di `7.3e-12` non-dim, cioe' il floor dell'integratore.

> **Possibile domanda d'esame** -- Come fai a sapere che il rilassamento e' lossless
> e non ti sei semplicemente risolto un problema *diverso* (piu' facile)?
> *Risposta:* Perche' `cone_gap` e' zero a precisione numerica: `sigma_k = |u_k|`
> per ogni k. Se la disuguaglianza `|u| <= sigma` e' attiva ovunque, il vincolo
> `z' = -Vc*sigma` coincide con `z' = -Vc*|u|`, cioe' la vera equazione della massa,
> e il bound `sigma <= Tmax*exp(-z)` coincide con `|u| <= Tmax*exp(-z)`, cioe' il
> vero bound di spinta. Ogni vincolo del problema rilassato collassa sul vincolo
> originale, quindi la soluzione e' ammissibile per l'originale -- ed essendo
> l'ottimo di un rilassamento (feasible set piu' grande), il suo costo e' un lower
> bound: coincidendo con un punto ammissibile dell'originale, e' l'ottimo globale
> dell'originale.

---

## `nondim` (righe 103-113)

- Riga 105-107: scale -- `L = y0 = 3000 m`, `g = 9.81`, `t = sqrt(L/g)`,
  `V = sqrt(g*L)`, `m = m0`, `T = m0*g`. La scala di accelerazione implicita e'
  `L/t^2 = g`, quindi `u` non-dim e' **letteralmente in "g"**.
- Riga 111: `dnd.Tmax = Tmax/(m0*g)` = 70000/19620 ~ 3.57, cioe' 3.57 g di
  accelerazione massima **alla massa iniziale**. E `dnd.Vc = ref.V/c`.
- Riga 110: `dnd.Tmin = 0/ref.T = 0`. **Variabile morta**: `d.Tmin` non compare mai
  nel SOCP (righe 125-177). Con `Tmin = 0` il vincolo sarebbe comunque vacuo
  (`sigma >= 0` e' gia' implicato da `sigma >= |u| >= 0`), ma il lettore va avvisato:
  se domani si volesse `Tmin > 0`, non basta cambiare `d_si.Tmin` -- bisogna
  *aggiungere* il vincolo `sigma >= Tmin*exp(-z)` al SOCP.
- Con `m0` non-dim = 1, si ha `z0 = ln(1) = 0` (usato alla riga 132).

---

## `dim_sol` (righe 115-123)

- Righe 117-122: ri-dimensionalizzazione. Riga 119-120: `Tx = m*ux`, `Ty = m*uy`
  ricostruiti a monte (riga 174) e qui riscalati; `Tmag = sqrt(Tx^2 + Ty^2)`.
- Attenzione al significato di `Tmag`: `T = m_k * u_k` usa la massa al nodo
  **sinistro** dell'intervallo. Dentro l'intervallo `u` e' costante ma `m(t)`
  cala, quindi la spinta vera `T(t) = m(t)*u` **decresce**. Il valore riportato e'
  dunque il **massimo** della spinta sull'intervallo. Conseguenza elegante: se il
  bound e' rispettato ai nodi, e' rispettato **con continuita'** su tutto
  l'intervallo -- non serve nessun vincolo inter-nodo per il bound superiore. (Il
  contrario varrebbe per un eventuale bound *inferiore* `Tmin > 0`: li' il punto
  critico sarebbe la fine dell'intervallo, e il check ai nodi non basterebbe.)

---

## `solve_gfold_socp` (righe 125-177)

E' il sotto-problema convesso: **un SOCP** costruito con YALMIP e risolto da ECOS.

```matlab
XI = sdpvar(5, N,   'full');   % [x; y; vx; vy; z]
W  = sdpvar(3, N-1, 'full');   % [ux; uy; sigma] per intervallo ZOH
```

- Righe 127-128: variabili di decisione. `5N + 3(N-1)` scalari: 250 + 147 = 397.
  Nota il **disallineamento voluto** stato/controllo: N nodi di stato, N-1 controlli
  (uno per intervallo), com'e' giusto per un ZOH.
- Riga 130: `z0 = log(d.m0)` = 0.
- Righe 132-133: condizione iniziale su tutte e 5 le componenti; condizione finale
  **solo su `XI(1:4,N) = 0`** (pinpoint + soft landing). `z_N` resta libera: e' la
  quantita' che si massimizza.

### Dinamica (riga 136)

    XI(:,k+1) == Abar*XI(:,k) + Bbar*W(:,k) + cbar

Vincolo di **uguaglianza lineare**: e' cio' che rende il problema un SOCP e non un
NLP. `Abar, Bbar, cbar` arrivano da `lti_zoh` e sono **le stesse per ogni k** (LTI).
In un problema LTV avresti `Abar_k`, e dovresti ricalcolarle a ogni iterazione
integrando l'ODE aumentata dell'Appendice A.

### Cono lossless (riga 137)

    norm(W(1:2,k)) <= W(3,k)          ->      |u_k| <= sigma_k

E' un vincolo conico del second'ordine (SOC) in forma standard. Questo e' il
rilassamento: al posto dell'uguaglianza non convessa `z' = -Vc*|u|` abbiamo
l'uguaglianza lineare `z' = -Vc*sigma` (dentro `Bbar`) piu' questa disuguaglianza
convessa.

**Perche' e' lossless (l'argomento serio).** Il costo (riga 164) e' `-z_N`, e
sommando la riga di massa lungo la griglia

    z_N = z_0 - Vc*dt*sum_k sigma_k

quindi massimizzare `z_N` significa **minimizzare `sum sigma_k`**. Lo slack ha
dunque un costo strettamente positivo: l'ottimizzatore lo vuole il piu' piccolo
possibile, e il suo unico lower bound e' `|u_k|`. L'intuizione "quindi
`sigma = |u|` all'ottimo" e' giusta, ma **attenzione**: non e' un argomento di
scambio banale, perche' abbassare `sigma_k` alza `z_j` per `j > k`, il che
*restringe* il bound superiore ai nodi successivi (piu' massa = meno accelerazione
disponibile). L'argomento rigoroso e' via PMP sul problema rilassato continuo:
l'Hamiltoniana e'

    H = lambda_r . v + lambda_v . (u + g) - alpha*lambda_z*sigma

e il controllo `u` compare **solo** nel termine `lambda_v . u`. Minimizzando H
rispetto a `u` a `sigma` fissato, con `|u| <= sigma`:

    u* = -sigma * lambda_v / |lambda_v|      ->      |u*| = sigma

ogni volta che il primer vector `lambda_v` e' diverso da zero. Il cono e' quindi
attivo quasi ovunque, **indipendentemente** da come viene scelto `sigma`. Il caso
patologico `lambda_v = 0` su un insieme di misura positiva viene escluso nel paper
con un argomento di controllabilita'. Nota che questo argomento **non e' nel
codice**: il codice si limita a *verificare a posteriori* la losslessness stampando
`cone_gap` (riga 86).

Interpretazione fisica bonus: `sum sigma_k * dt ~ integrale di |u| dt = integrale
di |T|/m dt`, che per Tsiolkovsky vale `ln(m0/mf)/alpha`, cioe' **il delta-V
ideale**. Minimizzare il carburante e' *esattamente* minimizzare il delta-V -- non
un'analogia, un'identita'.

### Bound superiore di spinta linearizzato (righe 138-140)

```matlab
ezr = exp(-z_ref(k));
cstr = [cstr, W(3,k) <= d.Tmax*ezr*(1 - (XI(5,k) - z_ref(k)))];
```

Il vincolo **vero** e' `|T| <= Tmax`, cioe', diviso per `m = exp(z)`:

    sigma <= Tmax * exp(-z)

L'insieme `{(z, sigma) : sigma <= Tmax*exp(-z)}` e' la regione **sotto** il grafico
di una funzione convessa: **non convesso** (prova: presi due punti sul grafico con
`z1 != z2`, il loro punto medio ha `sigma = Tmax*(e^{-z1}+e^{-z2})/2 >
Tmax*e^{-(z1+z2)/2}` per stretta convessita', quindi cade fuori).

Il codice lo sostituisce con la **retta tangente** a `exp(-z)` in `z_ref`:

    sigma <= Tmax * e^{-z_ref} * (1 - (z - z_ref))

che e' il Taylor **al prim'ordine**. Due osservazioni cruciali:

1. **E' un'approssimazione conservativa (inner), non una violazione.** Una funzione
   convessa sta sempre **sopra** la sua tangente: `exp(-z) >= e^{-z_ref}*(1 - (z -
   z_ref))` per ogni `z`. Quindi ogni punto che soddisfa il vincolo linearizzato
   soddisfa **anche** quello vero: la spinta `|T| = m*|u| <= m*sigma <=
   m*Tmax*e^{-z_ref}(1-(z-z_ref)) <= m*Tmax*e^{-z} = Tmax`. **Il bound di 70 kN non
   puo' essere violato**, per costruzione, a qualunque iterazione. Cio' che si puo'
   perdere e' ottimalita' (si sta risolvendo un problema piu' stretto del vero).
2. **La conservativita' e' nulla in `z = z_ref` ed e' tanto peggiore quanto piu'
   `z` si allontana da `z_ref`.** Da qui la necessita' di ri-centrare la tangente,
   che e' esattamente il compito del loop SCvx.

**Differenza dal paper.** Acikmese/Blackmore usano `z_ref(t) = ln(m0 - alpha*Tmax*t)`
(profilo di svuotamento a spinta massima) **fissato**, usano lo stesso Taylor al
prim'ordine per il bound superiore, e un Taylor al **second'ordine**
(`mu1*[1 - dz + dz^2/2]`, un vincolo quadratico convesso, quindi SOC-rappresentabile)
per il bound **inferiore** `sigma >= Tmin*exp(-z)`; poi risolvono **un solo** SOCP.
Questo codice: (a) non ha il bound inferiore (`Tmin = 0`, quindi il Taylor al
second'ordine **non compare da nessuna parte**); (b) parte da quel `z_ref` del paper
(riga 186) ma poi lo **aggiorna** iterativamente. Se qualcuno all'orale cita "Taylor
al second'ordine come nel paper", la risposta onesta e': nel paper c'e', **in questo
codice no**, perche' serve solo per `Tmin > 0`.

### Vincoli di percorso (righe 143-147)

- Righe 143-145: glide-slope. `|x| <= tan(theta_mx)*y` scritto come **coppia di
  disuguaglianze lineari** `x <= tan*y` e `-x <= tan*y`. Con `theta_mx = 60 deg`,
  `tan = 1.732`. Convesso e lineare -- nessun costo aggiuntivo per il SOCP.
- Riga 146: `y >= 0` (quota).
- Riga 147: `log(1e-3) <= z <= 0`, cioe' `1e-3*m0 <= m <= m0`. Il bound superiore
  e' fisico (non si guadagna massa). Il bound inferiore e' un **floor numerico
  arbitrario** (2 kg): impedisce a `z` di divergere verso `-inf` ma **non e' una
  massa a secco**. Il problema, cosi' com'e' scritto, **non ha un vincolo di massa
  minima strutturale**: nulla garantisce che il propellente consumato sia
  fisicamente disponibile. Il paper usa invece un lower bound *tempo-variante*,
  `z >= ln(m0 - alpha*Tmax*t)`, che serve anche a limitare l'errore
  dell'approssimazione di Taylor.
- **Limite comune a tutte le trascrizioni per collocazione**: glide-slope e `y >= 0`
  sono imposti **solo ai nodi**. Fra due nodi la traiettoria e' una parabola e in
  linea di principio puo' uscire dal corridoio. Il codice non lo controlla.

### Obiettivo e risoluzione (righe 164-167)

- Riga 164: `optimize(cstr, -XI(5,N), ...)` -- minimizza `-z_N`, cioe' massimizza
  `z_N`. **Perche' `z_N` e non `m_N`?** Perche' `m_N = exp(z_N)` sarebbe una
  funzione nonlineare della variabile di decisione, mentre il SOCP vuole un
  obiettivo lineare (o al piu' convesso). `exp` e' monotona crescente, quindi
  `argmax z_N = argmax m_N`: massimizzare il log-massa e' *equivalente* a
  massimizzare la massa. Il cambio di variabile regala anche questo.
- Righe 165-167: se ECOS restituisce un flag != 0 il codice **si limita a un
  warning e prosegue**. Su un problema infeasible `value(XI)` restituisce `NaN`, che
  poi si propaga silenziosamente nel replay e nei check. Fragilita' nota: dovrebbe
  essere un `error` o almeno un rifiuto dell'iterazione.

### Estrazione della soluzione (righe 169-176)

- Riga 173: i controlli vengono **paddati con uno zero** al nodo N (`sol.ux =
  [Wv(1,:).'; 0]`), perche' ci sono N-1 controlli e N nodi. Conseguenze: (i) il
  `cone_gap` alla riga 71 deve escludere il nodo N (e lo fa); (ii) nel grafico
  `stairs` della spinta (riga 92) l'ultimo gradino **crolla a zero** -- artefatto
  cosmetico, non un vero coast finale.
- Riga 174: `T = m .* u` con `m = exp(z)` nodo per nodo.

---

## `solve_gfold_scvx` (righe 179-239)

Loop esterno di successive convexification attorno al SOCP.

- **Riga 181**: `[Abar, Bbar, cbar] = lti_zoh(tf/(N-1), d.Vc)` -- **una volta sola,
  fuori dal loop**. E' la firma della trascrizione LTI: nelle varianti LTV le
  matrici discrete vanno ricostruite a ogni iterazione, integrando l'ODE aumentata
  dell'Appendice A per ognuno dei 49 intervalli.
- Righe 185-192: **warm start**. `m_apri = max(m0 - Vc*Tmax*t, 1e-2)` e
  `ref.z = log(m_apri)`: e' proprio il profilo `z0(t) = ln(m0 - alpha*rho2*t)` del
  paper, cioe' il consumo che si avrebbe **spingendo al massimo per tutto il volo**.
  E' quindi una *sottostima* della massa vera (la soluzione ottima e' un
  bang-off-bang con un arco di coast -- vedi README: max-coast-max), il che rende la
  prima tangente parecchio conservativa. Stato di riferimento: interpolazione
  lineare fra le condizioni al contorno; controllo di riferimento: `u = [0; 1]`
  (hover-ish, un g verso l'alto).
- Riga 194: dimensioni base della trust region (non-dim): `pos 0.5` (~1500 m),
  `vel 1.0` (~172 m/s), `lz 0.4`, `u 4.0` (4 g), `sig 4.0`. Sono **box molto larghi**.
- Riga 195: `rho = 1`, `rho_max = 1`, `rho_min = 1e-3`, `eta_l = 0.25`, `eta_h = 0.7`.
  Nota: `rho` **parte gia' al massimo**, quindi la riga 228 (`rho = min(rho_max,
  2*rho)`) non puo' mai farla crescere. La trust region puo' solo **restringersi**.
- Righe 202-212: alla **prima iterazione** il SOCP viene risolto **senza trust
  region** (`ref_sol = []`, `trust = []`). Il commento (righe 203-206) spiega il
  perche' ed e' corretto: la dinamica e' esatta, quindi non serve restare vicini al
  riferimento per fidarsi del modello; e il warm start e' cinematicamente incoerente
  (interpolazione lineare), quindi un box stretto attorno ad esso renderebbe il
  problema infeasible.
- Righe 214-217: **ratio test**. `J_pred = m_f(modello) - m_f(riferimento)`,
  `J_act = m_f(replay ode45) - m_f(riferimento)`, `eta = J_act/J_pred`.
- Riga 219: `delta` = norma della variazione di stato rispetto al riferimento.
- Righe 221-235: accetta se `eta >= 0.25`; se accettato aggiorna il riferimento,
  eventualmente allarga `rho` (di fatto mai, vedi sopra) e dichiara convergenza se
  `delta < tol = 1e-3`. Se rifiutato dimezza `rho` e si ferma se `rho < 1e-3`.

**Il ratio test qui e' quasi vuoto, e va detto.** L'unica sorgente di discrepanza
fra modello discreto e replay nonlineare e' il gap del cono (`sigma - |u|`): le
righe di posizione e velocita' sono riprodotte **esattamente** dal replay (perche'
`vx' = ux` e' vera anche nel modello nonlineare), e la riga di massa coincide non
appena `sigma = |u|`. Ma `sigma = |u|` e' garantito **a ogni ottimo del SOCP**, per
la struttura del costo. Quindi `eta ~ 1` gia' dalla prima iterazione, ogni passo
viene accettato, `rho` non si restringe mai e la trust region non morde. In pratica
**quello che fa convergere il loop non e' l'SCvx ma la ri-linearizzazione del bound
di spinta**, e il criterio che conta e' `delta < tol` (punto fisso della tangente),
non `eta`. Peggio: `eta` e' **strutturalmente cieco** all'errore che il loop sta
davvero correggendo, perche' l'errore di linearizzazione del bound `sigma <=
Tmax*exp(-z)` non influenza la *fedelta' dinamica* (il replay non conosce i
vincoli), solo la *feasibility/ottimalita'*. In sostanza questo e' un **successive
linear programming con safeguard**, non un SCvx nel senso pieno. Il README riporta
convergenza in 3 iterazioni (~5 s), coerente con questa lettura.

> **Possibile domanda d'esame** -- Se la dinamica e' esatta, perche' c'e' una trust
> region?
> *Risposta:* Non serve per la dinamica (nessun errore di linearizzazione li'), ma
> in linea di principio per l'unico vincolo linearizzato, `sigma <= Tmax*exp(-z)`:
> la tangente e' affidabile solo vicino a `z_ref`, quindi limitare `|z - z_ref|`
> (`trust.lz`) limita la conservativita'/il salto fra iterazioni. In questo codice
> pero' il box e' largo e il ratio test non misura quell'errore, quindi la trust
> region e' di fatto inerte. Nel GFOLD "puro" non c'e' ne' trust region ne' loop: si
> accetta la conservativita' della tangente attorno al profilo a-priori e si risolve
> un SOCP solo.

---

## `fwd_integrate_uacc` (righe 241-254)

- Righe 244-245: parte dallo stato iniziale **esatto** (non da quello del SOCP).
- Riga 246: `opts = odeset('RelTol',1e-10,'AbsTol',1e-12)` -- la tolleranza
  stretta prevista dalle convenzioni della repo.
- Righe 247-252: per ogni intervallo tiene **l'accelerazione** `u_k` costante e
  integra `ode_descent_uacc` con `ode45` usando quelle tolleranze.
- La convenzione ZOH e' `u` costante, **non** `T` costante: dentro l'intervallo
  `T(t) = m(t)*u` "galleggia" al calare della massa. E' la convenzione nativa del
  GFOLD ed e' esattamente quella che rende la mappa discreta esatta. Se si rigiocasse
  la soluzione tenendo `T` costante (`ode_descent.m`), il replay **non** tornerebbe:
  sarebbe un errore di consistenza, non un errore numerico.

---

## `ternary` (righe 256-258)

- Helper cosmetico per gli `fprintf`. Valuta entrambi i rami (qui sono stringhe o
  scalari gia' calcolati, quindi innocuo). Usato anche alla riga 217 per proteggere
  la divisione `J_act/J_pred` quando `|J_pred| < 1e-10`.

---

## Limiti noti, in chiaro

1. **Non e' il GFOLD del paper.** Il paper risolve **un SOCP** con certificato di
   ottimalita' globale; qui c'e' un loop SCvx di ri-linearizzazione. Il ticket
   "single SOCP" e' aperto nel README.
2. **`Tmin = 0`**: il bound inferiore di spinta (la non convessita' piu' interessante,
   quella per cui la lossless convexification e' famosa) **non e' implementato**.
   `d.Tmin` e' calcolato e mai usato.
3. **Nessuna massa a secco**: il floor `z >= log(1e-3)` e' numerico, non fisico.
4. **Vincoli solo ai nodi**: glide-slope e `y >= 0` non sono garantiti fra i nodi.
5. **Errore di solver ignorato**: flag ECOS != 0 produce solo un warning; i `NaN`
   si propagano.
6. **`rho` non puo' crescere** (parte a `rho_max`); `hist` e' calcolato e scartato;
   il grafico `stairs` della spinta mostra un finto azzeramento all'ultimo nodo per
   via del padding.
7. **Il ratio test non misura l'errore che il loop corregge** (vedi sopra).

---

## Possibili domande d'esame

**D: Deriva il cambio di variabile GFOLD e spiega perche' la dinamica risultante e'
LTI e non semplicemente linearizzata.**
R: Si pone `u = T/m` (accelerazione comandata) e `z = ln(m)`. Le righe di velocita'
`v' = T/m + g` diventano `v' = u + g`: la massa scompare. La riga di massa
`m' = -alpha*|T|` diventa `z' = m'/m = -alpha*|T|/m = -alpha*|u|`: la massa scompare
di nuovo. Aggiungendo lo slack `sigma >= |u|` si ha `z' = -alpha*sigma`, e il sistema
in `(r, v, z)` con ingresso `(u, sigma)` e' **affine a coefficienti costanti**:
`xi' = A*xi + B*w + c`, con `A`, `B`, `c` indipendenti da stato e controllo. Non c'e'
Taylor, non c'e' punto di lavoro, non c'e' resto trascurato: e' un cambio di
coordinate esatto (invertibile finche' `m > 0`), valido globalmente. Per questo la
discretizzazione ZOH e' esatta e le matrici discrete sono **le stesse per tutti gli
intervalli**, calcolate una volta con un `expm`.

**D: Il rilassamento con lo slack sigma: perche' e' *lossless*?**
R: Perche' all'ottimo il cono e' attivo, `|u| = sigma`, e quindi tutti i vincoli
rilassati collassano su quelli originali. La dimostrazione pulita e' via PMP:
nell'Hamiltoniana il controllo `u` compare solo nel termine `lambda_v . u`, quindi a
`sigma` fissato il minimo su `{|u| <= sigma}` si ottiene sul **bordo**,
`u* = -sigma*lambda_v/|lambda_v|`, cioe' `|u*| = sigma`, ogni volta che il primer
vector `lambda_v` non e' nullo (escluso su insiemi di misura positiva da un argomento
di controllabilita'). Intuizione complementare: `sigma` costa carburante (il costo e'
`-z_N = -z_0 + alpha*dt*sum sigma_k`), quindi l'ottimizzatore lo schiaccia verso il
basso finche' non sbatte contro `|u|`. Conseguenza pratica: la soluzione del problema
**convesso** rilassato e' ammissibile per il problema **non convesso** originale e, in
quanto ottimo di un rilassamento, ne e' l'ottimo **globale**. Nel codice questo si
verifica alla riga 71/86: `cone_gap = max|sigma_k - |u_k||` deve essere zero a
precisione numerica.

**D: Che fine ha fatto il vincolo non convesso `Tmin <= |T| <= Tmax`? Come lo tratta
il codice?**
R: Diviso per la massa diventa `Tmin*exp(-z) <= |u| <= Tmax*exp(-z)`. Il lower bound
su una norma non e' mai convesso: e' li' che lo slack fa il vero lavoro (si scrive
`sigma >= Tmin*exp(-z)`, convesso in `(z, sigma)` perche' epigrafo di funzione
convessa, e la losslessness restituisce `|u| = sigma`). **Ma questo codice ha
`Tmin = 0`, quindi quel pezzo non e' implementato**: `d.Tmin` e' una variabile morta.
Resta l'upper bound `sigma <= Tmax*exp(-z)`, non convesso perche' e' la regione sotto
il grafico di una funzione convessa. Il codice (riga 140) lo sostituisce con la
**tangente** in `z_ref`: `sigma <= Tmax*e^{-z_ref}*(1 - (z - z_ref))`. Essendo
`exp(-z)` convessa, la tangente le sta sempre **sotto**: il vincolo linearizzato e'
piu' stretto del vero, quindi **conservativo ma mai violato** (i 70 kN sono garantiti
per costruzione). Il loop SCvx (righe 179-239) ri-centra `z_ref` sulla soluzione
corrente per azzerare la conservativita' al punto fisso.

**D: Il paper usa un Taylor al second'ordine. Dov'e' nel codice?**
R: Non c'e'. Il Taylor al second'ordine del paper serve per il bound **inferiore**
`sigma >= Tmin*exp(-z)`: si approssima `exp(-z)` dall'alto con
`e^{-z0}*(1 - dz + dz^2/2)`, che e' un vincolo quadratico convesso (SOC-rappresentabile).
Con `Tmin = 0` quel vincolo non esiste e il codice usa **solo** la tangente
(prim'ordine) sul bound superiore. Dire il contrario all'orale sarebbe raccontare il
paper, non il codice.

**D: Perche' `lti_zoh` (discretizzazione esatta con `expm`) e non un RK4 come nella
variante (b)?**
R: Perche' qui il sistema e' davvero LTI e il controllo e' ZOH: la soluzione
dell'intervallo e' nota in forma chiusa,
`xi_{k+1} = e^{A*dt}*xi_k + (integrale_0^dt e^{A*s} ds)*(B*w_k + c)`, ottenibile con
un solo `expm` del sistema aumentato (van Loan). Un RK4 sarebbe un'**approssimazione
di qualcosa che si sa risolvere esattamente**: errore locale O(dt^5) invece di zero,
piu' costo, piu' Jacobiani. In piu' `A` e' nilpotente (`A^2 = 0`), quindi la serie si
tronca e le matrici discrete sono le formule elementari del moto uniformemente
accelerato: `Abar = I + A*dt`, `Bbar = dt*B + (dt^2/2)*A*B`, `cbar = dt*c +
(dt^2/2)*A*c`. Nella variante (b), invece, con **la spinta** (non l'accelerazione)
tenuta costante, `v' = T/m(t)` varia dentro l'intervallo: la mappa ZOH non e' piu'
affine nelle variabili di decisione, serve integrazione numerica o una linearizzazione
LTV, e le Jacobiane si portano dietro la massa. Guardando `main_task2.m:jacobians`
(righe 255-277): le **righe di velocita'** di `df/dx` vanno come `1/m^2`
(`A_jac(3,5) = -Tx/m^2`, `A_jac(4,5) = -Ty/m^2`, righe 270-271) e come `1/m` in
`df/du`; la **riga di massa** `m' = -Vc*|T|` non contiene affatto `m` (Jacobiano
rispetto a `m` identicamente nullo), ma e' **non differenziabile in `T = 0`**, ed e'
quella la singolarita' che il codice tampona con
`Tmag_reg = sqrt(Tx^2 + Ty^2 + 1e-6)` (riga 266) prima di usarla in
`B_jac(5,1:2)`. E' l'epsilon della "epsilon-regularised singular mass row" citata
nell'header (righe 17-19), e la trascrizione GFOLD la elimina alla radice: con
`z' = -Vc*sigma` la riga di massa e' lineare e non c'e' nulla da regolarizzare.

**D: Perche' l'obiettivo e' `-z_N` e non `-m_N`?**
R: Perche' `m_N = exp(z_N)` e' una funzione nonlineare della variabile di decisione, e
il SOCP richiede un obiettivo lineare/convesso. Essendo `exp` strettamente monotona,
`argmax z_N = argmax m_N`, quindi massimizzare il log-massa e' *equivalente* a
massimizzare la massa finale. In piu' `z_N = z_0 - Vc*dt*sum_k sigma_k` rende
esplicito che minimizzare il carburante equivale a minimizzare `sum sigma_k`, cioe'
(per Tsiolkovsky) il **delta-V ideale**.

**D: Come verifichi che la trascrizione sia fedele, e cosa ti aspetti di vedere?**
R: Con un replay open-loop: si prendono le accelerazioni ottime `u_k`, si integrano
con `ode45` (tolleranze 1e-10/1e-12) attraverso la dinamica **nonlineare** in massa
(`ode_descent_uacc`, riga 249) e si confrontano gli stati nodo per nodo. Poiche' la
mappa discreta e' esatta, ci si aspetta un errore al **floor dell'integratore** (il
README riporta 7.3e-12 non-dim), contro l'O(dt^2) della trascrizione trapezoidale.
Attenzione pero' a *cosa* misura ciascuna metrica. L'errore di nodo (riga 65) guarda
solo posizione e velocita', e quelle righe non contengono `sigma`: sono riprodotte
esattamente dal replay **comunque**, quindi sono cieche a un cono non attivo. La
patologia `sigma > |u|` (rilassamento non lossless) si manifesta nell'unico canale che
la vede, la massa: `dmf` (riga 64) e `cone_gap` (riga 71) misurano quella, da due
angoli diversi.
