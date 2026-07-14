# HM2_powered_descent/main_task2.m

## Ruolo del file nel progetto

Questo e' lo script piu' lungo e piu' denso della repo (1420 righe). Implementa il
**Task 2 (opzionale)** dell'Homework 2: la traccia chiede una trascrizione
alternativa del problema di powered descent basata su **Zero-Order Hold (ZOH)**,
seguendo l'Appendice A del PDF, e di validarla ri-propagando i controlli ottimizzati
con un integratore ad alta fedelta'.

Lo script non si limita a una trascrizione: ne implementa **quattro**, tutte sullo
stesso problema fisico (stessi dati, stessa griglia, stesso costo), piu' la baseline
trapezoidale del Task 1 come termine di paragone. Le quattro varianti sono, in ordine
di "quanta struttura convessa espongono al solver":

- **(a) ZOH nonlineare + RK4** (`solve_zoh`, righe 626-650): la dinamica resta
  nonlineare, il difetto di collocazione e' `x_{k+1} - RK4(x_k, u_k, dt) = 0`. E'
  un multiple shooting. Risolto con `fmincon`/SQP.
- **(b) ZOH LTV + SCvx con fmincon** (`solve_scvx`, righe 429-537): si linearizza la
  dinamica attorno a una traiettoria di riferimento, si costruiscono le matrici
  discrete `Abar_k, Bbar_k, cbar_k` con l'ODE aumentata dell'Appendice A, si risolve
  il sottoproblema LTV, si aggiorna il riferimento, si itera.
- **(c) ZOH LTV + SCvx con YALMIP/ECOS** (`solve_scvx_yalmip`, righe 1042-1121):
  stesso loop esterno di (b), ma il sottoproblema interno e' scritto come **SOCP**
  (Second-Order Cone Program) e dato a un solver conico invece che a un NLP generico.
- **(d) GFOLD log-mass** (`solve_gfold_scvx`, righe 1198-1283): cambio di variabili
  `z = ln(m)`, `u = T/m`, slack `sigma >= |u|`. La dinamica diventa **esattamente
  LTI**, quindi la ZOH dell'Appendice A si riduce a un singolo esponenziale di
  matrice; l'unica cosa che resta da linearizzare e' il bound superiore di spinta.

Il file dipende da quattro mattoni esterni condivisi con `main_task1.m` e con la test
suite: `ode_descent.m` (RHS nonlineare, controllo = spinta), `ode_descent_uacc.m`
(RHS con controllo = accelerazione, convenzione GFOLD), `rk4_zoh.m` (propagatore RK4
a passo fisso) e `lti_zoh.m` (ZOH esatta via `expm` per il sistema log-mass).
Tutto il resto -- Jacobiane, ODE aumentata, assemblaggio degli NLP, loop SCvx, replay
di validazione, plotting -- vive in **local functions** dentro questo file (dalla riga
203 in poi).

Il flusso e' lineare: dati SI -> non-dimensionalizzazione -> risolvi le 5 trascrizioni
-> ri-propaga i controlli con `ode45` -> misura l'errore ai nodi e la dispersione al
touchdown -> tabelle e figure.

---

## Il problema di partenza e perche' NON e' convesso

La dinamica continua (non-dim, gravita' = 1, stato `x = [x; y; vx; vy; m]`,
controllo `T = [Tx; Ty]`) e' quella di `ode_descent.m`:

    x_dot  = vx
    y_dot  = vy
    vx_dot = Tx / m
    vy_dot = Ty / m - 1
    m_dot  = -Vc * |T|,        Vc = V_ref / c

con `c = Isp * g0` velocita' efficace di scarico. Costo: massimizzare `m(tf)`
(equivalente a minimizzare il carburante, perche' `m0` e' fissa).

Ci sono **due** sorgenti di non convessita', ed e' importante non confonderle:

1. **La dinamica e' bilineare / nonlineare.** I termini `Tx/m` e `Ty/m` sono un
   prodotto fra controllo e (inverso della) massa: bilineari. Il termine `-Vc*|T|`
   e' una norma dentro un vincolo di **uguaglianza** (l'equazione di stato). Anche se
   `|T|` e' una funzione convessa, un vincolo `m_dot = -Vc*|T|` con la funzione
   convessa **al secondo membro di un'uguaglianza** definisce un insieme non convesso.
   Quindi il set di traiettorie ammissibili non e' convesso.

2. **Il vincolo di spinta minima.** In generale il motore ha `T_min <= |T| <= T_max`:
   e' un **guscio anulare** (in 2D una corona circolare), il classico insieme non
   convesso. In *questo* homework `T_min = 0` (riga 32), quindi l'insieme e' un disco
   pieno, gia' convesso -- ma la macchineria per gestire `T_min > 0` (lossless
   convexification) va comunque capita, perche' e' quello che la traccia e la
   letteratura GFOLD hanno in mente e perche' lo slack `sigma` serve comunque, come
   vedremo, per un'altra ragione.

Nota bene: i **vincoli di percorso** dell'homework sono invece gia' convessi.
Il glide-slope `|x| <= tan(theta_max) * y` e' una coppia di semipiani (righe 817-818),
il bound `|T| <= T_max` e' un cono di Lorentz, i box su massa e quota sono lineari.
Quindi tutta la non convessita' e' concentrata **nella dinamica** (piu' `T_min` se
fosse > 0). Questo e' esattamente il motivo per cui la ricetta SCvx funziona bene qui:
basta convessificare la dinamica.

---

## Teoria: la trascrizione ZOH per sistemi LTV (Appendice A)

Sia il sistema **lineare tempo-variante**

    x_dot(t) = A(t) x(t) + B(t) u(t) + c(t)

su una griglia uniforme `t_k`, con passo `dt`, e controllo tenuto costante
sull'intervallo (**zero-order hold**): `u(t) = u_k` per `t in [t_k, t_{k+1})`.

La formula di variazione delle costanti da'

    x(t) = Phi(t, t_k) x_k + int_{t_k}^{t} Phi(t, s) [B(s) u_k + c(s)] ds

dove `Phi(t, s)` e' la matrice di transizione di stato del sistema omogeneo. Valutando
in `t = t_{k+1}` si ottiene la ricorsione discreta

    x_{k+1} = Abar_k x_k + Bbar_k u_k + cbar_k

    Abar_k = Phi(t_{k+1}, t_k)
    Bbar_k = int_{t_k}^{t_{k+1}} Phi(t_{k+1}, s) B(s) ds
    cbar_k = int_{t_k}^{t_{k+1}} Phi(t_{k+1}, s) c(s) ds

Il punto e' che `u_k` esce dall'integrale **solo perche' e' costante**: e' tutta qui la
ragion d'essere della ZOH. Con un first-order hold (controllo lineare a tratti)
servirebbe un secondo integrale ausiliario.

**Come si calcolano gli integrali.** L'Appendice A usa l'identita' di semigruppo

    Phi(t_{k+1}, s) = Phi(t_{k+1}, t_k) * Phi(t_k, s)

che permette di portare fuori dall'integrale il fattore fisso `Phi(t_{k+1}, t_k)`.
Definendo la transizione **all'indietro** `Psi(s) := Phi(t_k, s) = Phi(s, t_k)^{-1}` e
i due integrali (riferiti a `t_k`)

    beta(t)  = int_{t_k}^{t} Psi(s) B(s) ds
    gamma(t) = int_{t_k}^{t} Psi(s) c(s) ds

si ha

    Bbar_k = Phi(t_{k+1}, t_k) * beta(t_{k+1})
    cbar_k = Phi(t_{k+1}, t_k) * gamma(t_{k+1})

e le quantita' `(Phi, Psi, beta, gamma)` obbediscono a un sistema di **ODE lineari**
che si integrano tutte insieme su un intervallo:

    Phi_dot   =  A Phi,          Phi(t_k)   = I
    Psi_dot   = -Psi A,          Psi(t_k)   = I
    beta_dot  =  Psi B,          beta(t_k)  = 0
    gamma_dot =  Psi c,          gamma(t_k) = 0

L'equazione per `Psi` si ricava derivando l'identita' `Phi^{-1} Phi = I`:

    d(Phi^{-1})/dt = -Phi^{-1} * Phi_dot * Phi^{-1} = -Phi^{-1} A Phi Phi^{-1}
                   = -Psi A

Le equazioni per `beta` e `gamma` sono **quadrature pure**: il loro integrando dipende
solo da `s`, non da `beta`/`gamma` stessi. Nessun feedback -- ecco perche' nel codice
(riga 292-293) si dice esplicitamente che le loro derivate non dipendono da loro.

Questo e' letteralmente cio' che fanno `ltv_aug_rhs` (righe 279-303) e
`compute_ltv_zoh` (righe 305-339).

> **Possibile domanda d'esame** -- Perche' portarsi dietro `Psi = Phi^{-1}` come stato
> aumentato, invece di calcolare `Phi(t_{k+1}, s)` direttamente?
> *Risposta:* perche' `Phi(t_{k+1}, s)` ha l'**estremo finale fisso e quello iniziale
> variabile**: integrandolo in avanti da `t_k` non lo si ottiene mai in un colpo solo.
> Riferendo tutto a `t_k` con l'identita' di semigruppo, tutte e quattro le quantita'
> partono da condizioni iniziali note in `t_k` e si integrano in **una sola passata in
> avanti**; si paga un prodotto matriciale finale (righe 336-337) per ricostruire
> `Bbar_k`, `cbar_k`. E' il prezzo della convenienza numerica.

---

## Teoria: Successive Convexification (SCvx)

### Linearizzazione

Il problema di partenza e' `x_dot = f(x, u)` con `f` nonlineare. SCvx lo affronta cosi':
si prende una **traiettoria di riferimento** `(x_ref(t), u_ref(t))` e si espande `f` al
primo ordine attorno ad essa:

    f(x, u) ~= f(x_ref, u_ref) + A (x - x_ref) + B (u - u_ref)

con

    A(t) = df/dx |_ref     (5x5)
    B(t) = df/du |_ref     (5x2)

Raccogliendo, si ottiene la forma LTV **in variabili assolute** (non in deviazioni):

    x_dot ~= A(t) x + B(t) u + c(t)

    c(t) = f(x_ref, u_ref) - A x_ref - B u_ref        <- termine affine / residuo

Il termine `c(t)` e' il pezzo che spesso si dimentica: e' quello che fa passare la
retta tangente **per il punto di riferimento**. Se lo togli, la linearizzazione passa
per l'origine e il modello e' completamente sbagliato. Nel codice e' la riga 296
(`c_off = f_val - A_jac*x_ref - B_jac*u_k`).

Per il nostro `f`, le Jacobiane si calcolano a mano (funzione `jacobians`, righe
255-277):

    A = [ 0 0 1 0      0
          0 0 0 1      0
          0 0 0 0  -Tx/m^2
          0 0 0 0  -Ty/m^2
          0 0 0 0      0    ]

    B = [   0        0
            0        0
          1/m        0
            0      1/m
      -Vc*Tx/|T|  -Vc*Ty/|T| ]

Le righe 3 e 4 di `A` vengono da `d/dm (Tx/m) = -Tx/m^2`. La riga 5 di `A` e' nulla
perche' `m_dot = -Vc*|T|` non dipende dallo stato. La riga 5 di `B` e'
`d/dT (-Vc*|T|) = -Vc * T/|T|`, cioe' `-Vc` per il **versore** di spinta: e' singolare
in `T = 0` (il gradiente della norma non esiste nell'origine). Il codice regolarizza con
`sqrt(Tx^2 + Ty^2 + 1e-6)` (riga 266). Torneremo su questo hack.

Con la dinamica linearizzata, tutti i vincoli del problema diventano convessi
(uguaglianze lineari + semipiani + un cono per nodo), quindi il **sottoproblema e' un
SOCP**: ha un ottimo globale unico raggiungibile in tempo polinomiale.

### Il ciclo

    1. parti da un riferimento (x_ref, u_ref)
    2. calcola A, B, c lungo il riferimento; discretizza (ZOH) -> Abar, Bbar, cbar
    3. risolvi il sottoproblema convesso dentro una trust region
    4. la soluzione diventa il nuovo riferimento
    5. torna al punto 2 finche' ||x^{j+1} - x^{j}|| < tol

Il fatto che il sottoproblema sia risolto **globalmente** non significa che SCvx trovi
l'ottimo globale del problema originale: converge a un punto stazionario (KKT) del
problema nonlineare. La convessita' serve a rendere ogni passo affidabile e veloce,
non a certificare l'ottimalita' globale.

### La trust region: perche' serve e com'e' implementata qui

La linearizzazione e' accurata **solo vicino al riferimento**. Se il sottoproblema
convesso e' lasciato libero, succede il fenomeno dell'**artificial unboundedness**: il
modello lineare "promette" un guadagno che la dinamica vera non consegna. Nel nostro
caso il colpevole tipico e' proprio la riga della massa: con `|T|` piccolo, `df5/dT`
e' quasi singolare, il modello lineare crede di poter guadagnare massa e il solver ci
si butta dentro (nel report questo e' descritto come "mass creation", `m_dot > 0`).

Nel codice la trust region e':

- **Hard**, non una penalita' soft. E' un box in norma infinito, per-variabile e
  per-nodo, **intersecato con i box bounds** gia' esistenti (`apply_trust`, righe
  571-600): `lb = max(lb, ref - raggio)`, `ub = min(ub, ref + raggio)`.
- **Anisotropa**: raggi diversi per posizione, velocita', massa, spinta
  (riga 49: `pos 0.17, vel 0.6, mass 0.1, thrust 1.0` in unita' non-dim -- circa
  510 m, 103 m/s, 200 kg, 19620 N).
- **Adattiva**, scalata da un fattore `rho in [1e-3, 1]` (righe 464-468).

L'aggiornamento e' il **ratio test** classico (righe 487-530). Definendo:

    J_pred = m_f(candidato secondo il modello LTV) - m_f(riferimento)
    J_act  = m_f(candidato ri-propagato con ode45) - m_f(riferimento)
    eta    = J_act / J_pred

cioe' **guadagno reale / guadagno predetto**, si ha:

    eta <  eta_l = 0.25    -> il modello lineare mente: RIFIUTA, rho <- rho/2
    0.25 <= eta < eta_h=0.7 -> accetta, rho invariato
    eta >= eta_h = 0.7     -> il modello e' ottimo: accetta, rho <- min(1, 2*rho)

Se `rho` scende sotto `rho_min = 1e-3`, il loop si ferma ("trust region collapsed",
riga 527). Notare che `rho_max = 1.0`: la trust region non puo' mai crescere oltre i
raggi base.

**Un dettaglio importante e specifico di questo codice:** il "guadagno reale" e'
misurato ri-propagando i controlli candidati attraverso la **dinamica nonlineare vera**
con `ode45` (riga 489). Non e' un semplice ricalcolo del costo: e' un test di
fedelta' dinamica. Questo e' un punto forte, ed e' anche il motivo per cui il loop
sopravvive senza virtual control (vedi sotto).

### Virtual control / slack: perche' servono e perche' QUI NON CI SONO

Il secondo modo di fallire di SCvx e' l'**artificial infeasibility**: il sottoproblema
linearizzato puo' risultare **infeasible anche se il problema nonlineare e' feasible**.
Succede quando la dinamica linearizzata, combinata con le condizioni al contorno
(qui: stato iniziale completamente assegnato **e** posizione/velocita' finali nulle),
non ammette nessuna soluzione: il modello lineare semplicemente non riesce a "chiudere"
il problema ai due estremi.

La cura canonica e' rilassare i vincoli che possono diventare vuoti:

- **virtual control** `nu_k`: si aggiunge un termine libero alla dinamica discreta

      x_{k+1} = Abar_k x_k + Bbar_k u_k + cbar_k + nu_k

  Con `nu` libero, la dinamica e' **sempre** soddisfacibile: il sottoproblema non puo'
  piu' essere infeasible per colpa della dinamica.

- **buffer / slack** `w >= 0` sui vincoli di percorso linearizzati:

      g(x_ref) + dg/dx (x - x_ref) + w <= 0

- entrambi vengono **penalizzati in costo** con una penalita' esatta L1:

      J_nu = lambda_nu * ||nu||_1 + lambda_w * ||w||_1

  I pesi vengono alzati finche' `nu -> 0` e `w -> 0` alla convergenza. La penalita' L1
  e' *esatta*: per `lambda` abbastanza grande, se esiste una soluzione con `nu = 0` il
  minimo la trova, senza dover mandare `lambda -> infinito`.

**Nel codice di `main_task2.m` NON ci sono ne' virtual control ne' buffer.** E' la
deviazione piu' importante rispetto alla ricetta canonica. Il loop se la cava per due
motivi: (i) il warm start dalla soluzione trapezoidale mette il riferimento gia' vicino
all'ottimo; (ii) quando il sottoproblema *diventa* infeasible, il solver interno
restituisce comunque il suo iterato meno infeasible, e il ratio test con `ode45` lo
vaglia prima di accettarlo. Il codice si limita a stampare un warning
(riga 410 per fmincon, riga 1028 per YALMIP). E' un comportamento fragile: un run
cold-started, o una trust region piu' larga, potrebbero rompersi.

> **Possibile domanda d'esame** -- Trust region e virtual control curano due patologie
> diverse. Quali?
> *Risposta:* la **trust region** cura l'*artificial unboundedness*: impedisce al
> sottoproblema linearizzato di allontanarsi dal riferimento fin dove il modello
> lineare promette guadagni che la dinamica vera non conferma. I **virtual control**
> curano l'*artificial infeasibility*: rendono il sottoproblema sempre risolvibile,
> aggiungendo alla dinamica un termine libero penalizzato in costo, che si annulla alla
> convergenza. Sono complementari: la prima restringe, i secondi rilassano.

---

## Teoria: lossless convexification e il cambio di variabili GFOLD

Questo e' il cuore teorico dell'homework, ed e' la variante (d).

### Il cambio di variabili

Seguendo Acikmese & Ploen (2007) e Acikmese & Blackmore (GFOLD -- *Guidance for
Fuel-Optimal Large Diverts*), si introducono:

    z     = ln(m)                 <- log-massa (nuovo stato)
    u     = T / m                 <- accelerazione di spinta (nuovo controllo)
    sigma >= |u|                  <- slack scalare (nuovo controllo ausiliario)

Sostituiamo nella dinamica:

- `vx_dot = Tx/m = ux`  -> **lineare** (il prodotto bilineare e' sparito per costruzione)
- `vy_dot = Ty/m - 1 = uy - 1` -> **lineare** (piu' un termine affine costante)
- massa: da `m_dot = -Vc*|T|` e `z = ln(m)`,

      z_dot = m_dot / m = -Vc * |T| / m = -Vc * |u|

  e **sostituendo `|u|` con lo slack `sigma`**:

      z_dot = -Vc * sigma        -> **lineare**

Risultato: il sistema in `xi = [x; y; vx; vy; z]` con controllo `w = [ux; uy; sigma]`
e' **esattamente LTI**:

    xi_dot = A xi + B w + c

    A = [0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0]
    B = [0 0 0; 0 0 0; 1 0 0; 0 1 0; 0 0 -Vc]
    c = [0; 0; 0; -1; 0]                      <- gravita'

Questo e' esattamente `lti_zoh.m` (righe 23-25 di quel file).

**Attenzione al ruolo di `sigma`.** In questo homework `T_min = 0`, quindi il vincolo
di spinta e' gia' convesso e -- a rigore -- lo slack "non servirebbe" per convessificare
il bound. Ma `sigma` serve **comunque**, per un motivo diverso: senza di lui la riga
della massa sarebbe `z_dot = -Vc*|u|`, cioe' una **uguaglianza con una norma dentro**:
non lineare, non convessa. Introdurre `sigma` con la **disuguaglianza** `|u| <= sigma`
(che e' un cono di Lorentz, convesso) e mettere `sigma` -- non `|u|` -- nella dinamica
e' quello che rende il sistema LTI. E' questo il punto che spesso si perde.

### Perche' il rilassamento e' LOSSLESS

Il rilassamento sostituisce l'insieme `{ (u, sigma) : sigma = |u| }` (non convesso) con
`{ (u, sigma) : sigma >= |u| }` (convesso). E' un rilassamento: si e' allargato
l'insieme ammissibile. **Lossless** significa che all'ottimo il vincolo e' **attivo**,
cioe' `sigma* = |u*|`, quindi la soluzione del problema rilassato e' anche soluzione
del problema originale -- niente e' stato perso.

Nel caso di questo codice la dimostrazione e' quasi banale ed e' istruttivo farla, perche'
la struttura ZOH la rende esplicita. Poiche' la riga 5 di `A` e' nulla, l'esponenziale
di matrice ha `exp(A*s)` con quinta riga `[0 0 0 0 1]`, e quindi la discretizzazione
esatta della riga della massa e':

    z_{k+1} = z_k - Vc * dt * sigma_k

cioe' `Bbar(5,3) = -Vc*dt` e `cbar(5) = 0`. Iterando:

    z_N = z_0 - Vc * dt * sum_{k=1}^{N-1} sigma_k

Il **costo** e' `max z_N` (righe 1183: `optimize(cstr, -XI(5,N))`). Ma `z_0` e' fissato
dalla condizione iniziale, quindi

    max z_N   <=>   min sum_k sigma_k

Ogni `sigma_k` compare quindi nel costo con **coefficiente positivo** (lo si vuole il
piu' piccolo possibile), e negli altri vincoli solo come:
- `sigma_k >= |u_k|`   (limite inferiore),
- `sigma_k <= T_max * e^{-z}` (limite superiore, che abbassare non viola).

Quindi all'ottimo ogni `sigma_k` viene **schiacciato sul suo limite inferiore**:
`sigma_k = |u_k|`. Il cono e' attivo. Il rilassamento e' lossless. Nessun `sigma`
"sprecato" darebbe vantaggio: brucerebbe carburante senza produrre accelerazione.

Nota che questo argomento **non e' quello generale**. Il risultato forte di
Acikmese & Ploen e' che il rilassamento resta lossless anche con `T_min > 0`, dove il
vincolo `T_min <= sigma <= T_max` **non** ha piu' `|u|` come lower bound naturale --
li' la dimostrazione passa per il Principio del Massimo di Pontryagin (un cono
inattivo forzerebbe `lambda_m(tf) = 0`, in contraddizione con la condizione di
trasversalita'). Nel codice `T_min = 0`, quindi la parte "difficile" del teorema non
viene esercitata: siamo nel caso facile.

### L'unica non convessita' superstite: il bound superiore di spinta

Il vincolo fisico e' `|T| <= T_max`. Con `T = m*u = e^z * u`:

    |u| <= T_max * e^{-z}      cioe'      sigma <= T_max * e^{-z}

`e^{-z}` e' **convessa**, quindi `sigma <= (funzione convessa di z)` definisce un
insieme **non convesso** (il sopragrafo di una funzione convessa non e' convesso).
Ecco perche' serve ancora un loop SCvx, sia pure per un solo vincolo scalare.

Si linearizza attorno a `z_ref`:

    e^{-z} ~= e^{-z_ref} * (1 - (z - z_ref))

da cui il vincolo del codice (righe 1158-1159):

    sigma_k <= T_max * e^{-z_ref,k} * (1 - (z_k - z_ref,k))

Poiche' `e^{-z}` e' convessa, la sua **tangente sta sempre sotto** la funzione: la
retta linearizzata **sottostima** il bound vero. Quindi il vincolo linearizzato e'
**conservativo** (approssimazione interna): ogni iterato accettato e' **feasible per il
problema originale**. E' una proprieta' molto comoda -- il conservatorismo si riduce
man mano che `z_ref` insegue la soluzione, ma non si viola mai il bound fisico.

> **Possibile domanda d'esame** -- Il rilassamento `sigma >= |u|` allarga l'insieme
> ammissibile: perche' non "bara"?
> *Risposta:* perche' `sigma` entra nella dinamica della massa con segno negativo
> (`z_dot = -Vc*sigma`) e il costo e' `max z_N`, quindi il costo e' una funzione
> **monotona decrescente** di ogni `sigma_k`. L'ottimizzatore ha percio' un incentivo
> stretto a comprimere ogni `sigma_k` sul suo limite inferiore `|u_k|`: il cono e'
> attivo all'ottimo e la soluzione rilassata coincide con quella originale. Con
> `T_min > 0` l'argomento e' meno immediato e richiede il PMP (Acikmese & Ploen 2007).

---

## `%% Problem data (Table 1, dimensional)` (righe 25-39)

- Righe 26-33: i dati della Tabella 1 della traccia, in SI. Stato iniziale
  `(1000, 3000) m`, velocita' `(300, -200) m/s`, `m0 = 2000 kg`, `Isp = 225 s`,
  `Tmax = 70 kN`, glide-slope `theta_max = 60 deg`.
- Riga 31: `c = Isp * g0` con `g0 = 9.80665` (gravita' standard, non `g = 9.81` del
  pianeta). E' la distinzione classica: `g0` e' una **costante di conversione** per
  passare da Isp a velocita' efficace di scarico, `g` e' l'accelerazione di gravita'
  del corpo su cui si atterra. Confonderle e' un errore comune.
- Riga 32: **`data.Tmin = 0`**. Come discusso, questo elimina la non convessita'
  anulare del vincolo di spinta. Tutta la macchineria GFOLD implementata qui e'
  quindi esercitata nel suo caso "facile".
- Riga 35: `tf = 38 s` **fissato** (non e' una variabile di decisione). Questo e' cio'
  che tiene le trascrizioni lineari: con `tf` libero il problema diventa nonlineare
  anche nel caso GFOLD (bisognerebbe fare una ricerca 1-D su `tf`).
- Riga 36: `N = 50` nodi -> `dt = 38/49 = 0.7755 s`.
- Riga 37: `n_sub = 2` sotto-passi RK4 per intervallo ZOH (variante a).
- Righe 38-39: `scvx_max = 15` iterazioni esterne, `scvx_tol = 1e-3` (non-dim). Il
  commento dice "(variant b)" ma i valori sono passati a **tutti e tre** i loop SCvx
  (righe 79, 90, 104): il commento e' fuorviante.

## `%% Non-dimensionalisation` (righe 41-49)

- Riga 42: chiama `nondim` (righe 203-228). Lo schema di riferimento e':

      L_ref = y0 = 3000 m          (quota iniziale)
      g_ref = g  = 9.81 m/s^2
      t_ref = sqrt(L/g) = 17.487 s
      V_ref = sqrt(g*L) = 171.55 m/s
      m_ref = m0 = 2000 kg
      T_ref = m0*g = 19620 N       (peso iniziale)
      Vc    = V_ref / c = 0.0777   (numero di Tsiolkovsky)

  Con questa scelta la gravita' non-dim vale **esattamente 1** (ecco perche' in
  `ode_descent.m` la riga di `vy_dot` e' `u(2)/x(5) - 1`), e `m0_nd = 1`,
  `Tmax_nd = 70000/19620 = 3.568` (cioe' il motore puo' dare ~3.57 g a massa piena).
  `tf_nd = 38/17.487 = 2.173`.
- Riga 46: `tf_nd = tf / ref.t`.
- Riga 49: i **raggi base della trust region**, in unita' non-dim:

  ```matlab
  trust = struct('pos', 0.17, 'vel', 0.6, ...
                 'mass', 0.1, 'thrust', 1.0);
  ```

  In SI: `0.17*3000 = 510 m`, `0.6*171.55 = 103 m/s`, `0.1*2000 = 200 kg`,
  `1.0*19620 = 19620 N`. Sono raggi **larghi**: il warm start e' gia' vicino
  all'ottimo, il ratio test si occupa di stringerli. Questi raggi valgono solo per le
  varianti (b) e (c): la (d) ha i suoi, hard-coded alla riga 1231.

## `%% YALMIP + ECOS availability check` (righe 51-55)

- Riga 52: `yalmip_ok = exist('yalmip','file') && exist('ecos','file')`. Se manca uno
  dei due, le varianti (c) e (d) vengono **saltate con grazia** e lo script cade su
  `plot_compare3`. E' la ragione per cui la repo dichiara "no external dependency" per
  le varianti (a)/(b): sono pure `fmincon`.

## `%% Solve all three transcriptions` (righe 57-110)

Il banner dice "three" ma le trascrizioni risolte sono cinque (trapezoidale + quattro
varianti ZOH): il commento e' rimasto indietro rispetto al codice.

- Righe 58-64: baseline **trapezoidale** (Task 1), `solve_trap`. Serve a due cose: (i)
  termine di paragone, (ii) **warm start** per le varianti (b) e (c).
- Righe 66-72: variante **(a)**, `solve_zoh` con `n_sub = 2`.
- Righe 74-85: variante **(b)**, `solve_scvx`, **warm-startata da `sol_trap_nd`**
  (riga 80). Il commento alle righe 76-77 e' esplicito: la linearizzazione LTV e'
  accurata solo localmente, quindi serve un riferimento gia' quasi ottimo piu' una
  trust region hard.
- Righe 87-96: variante **(c)**, `solve_scvx_yalmip`, stesso warm start, stessi raggi.
- Righe 98-109: variante **(d)**, `solve_gfold_scvx`. Notare la firma: **non riceve
  ne' `sol_trap_nd` ne' `trust`** (riga 104). E' **self-starting** -- si costruisce da
  sola il riferimento (riga 1223). E' la conseguenza diretta del fatto che la sua
  dinamica e' esatta: non c'e' errore di modello da cui difendersi con un buon warm
  start.

Ogni solve e' cronometrato con `tic`/`toc` e la soluzione riportata in SI con
`dim_sol`. `sol.iter` (righe 83, 94, 107) raccoglie il numero di iterazioni SCvx.

## `%% Forward-integration validation` (righe 112-137)

E' il cuore della **validazione** richiesta dalla traccia: si prendono i controlli
ottimizzati e li si "rivola" attraverso la dinamica continua vera con `ode45`
(`RelTol = 1e-10`, `AbsTol = 1e-12`, riga 839), campionando ai nodi della griglia.

- Righe 114-116: `fwd_integrate` con la convenzione di interpolazione **coerente con la
  trascrizione**: `'pwl'` (piecewise-linear) per la trapezoidale, `'zoh'`
  (piecewise-constant) per le ZOH. Questo e' fondamentale: replayare i controlli
  trapezoidali come ZOH darebbe un errore artificialmente enorme, e viceversa.
- Righe 124-127: la variante GFOLD richiede un replay **diverso**
  (`fwd_integrate_uacc`, righe 1285-1303), perche' tiene costante **l'accelerazione**
  `u = T/m`, non la spinta `T`. Con `u` costante la spinta vera
  `T(t) = m(t) * u` **fluttua** con la massa che si consuma. E' la convenzione ZOH
  nativa di GFOLD, e sbagliarla vanificherebbe tutta la fedelta' guadagnata.
- Righe 118-120, 123, 127: `node_err` calcola la norma dell'errore ai nodi.

> **Possibile domanda d'esame** -- Perche' GFOLD ha bisogno di un replay diverso dalle
> altre varianti ZOH?
> *Risposta:* perche' il suo controllo non e' la spinta `T` ma l'accelerazione
> `u = T/m`. Tenere `u` costante su un intervallo significa che `T = m(t)*u` varia con
> la massa che si esaurisce. Se lo si replayasse tenendo `T` costante si simulerebbe
> una legge di controllo diversa da quella ottimizzata, e l'errore di fedelta' sarebbe
> dominato da questo mismatch di convenzione, non dalla qualita' della trascrizione.

## `%% Replay landing accuracy + wall-time summary` (righe 139-165)

- Righe 142-144: la funzione anonima `land` estrae tre numeri dal replay:

  ```matlab
  land = @(s_nd, X) [norm(X(end,1:2)) * ref.L, ...
                     norm(X(end,3:4)) * ref.V, ...
                     (s_nd.m_f - X(end,5)) * ref.m];
  ```

  cioe' **dispersione di posizione** al touchdown (norma della posizione finale
  replayata rispetto al bersaglio `(0,0)`), **dispersione di velocita'** (norma della
  velocita' finale, bersaglio `(0,0)`) e **drift della massa finale** (differenza fra
  la `m_f` che la trascrizione crede di avere e quella che il replay consegna).
- Righe 149-165: tabella formattata. E' la figura di merito piu' operativa
  dell'homework: *se volassi questa legge di controllo open-loop, dove atterrerei?*

## `%% Plots` e `%% Export figures` (righe 167-193)

- Righe 168-175: `plot_compare5` se YALMIP c'e', `plot_compare3` altrimenti.
- Righe 178-193: esporta **tutte** le figure aperte in `figures/`, forzando il tema
  chiaro (riga 186, `theme(fig,'light')`, con fallback `Color = 'w'` per MATLAB
  pre-R2025a). I nomi file sono slugificati dal `'Name'` della figura (riga 180).

---

## `nondim` (righe 203-228) e `dim_sol` (righe 230-253)

Boundary helper: qui l'`arguments` validation **c'e'** (righe 210-212, 237-240), a
differenza delle funzioni hot-loop. E' una scelta esplicita documentata nel banner alle
righe 195-201: validare una volta per run al confine, mai dentro i loop di `fmincon` /
`ode45` (che chiamano quelle funzioni ~1e6 volte).

`dim_sol` (righe 241-252) riporta la soluzione in SI e calcola `sol.fuel = (m0 - m_f) *
m_ref`. Nota che **non** copia `sol.iter`: per questo lo script deve ricopiarlo a mano
alle righe 83, 94, 107. Piccola scomodita' del design.

## `jacobians` (righe 255-277)

```matlab
Tmag_reg = sqrt(Tx^2 + Ty^2 + 1e-6);
...
A_jac(3,5) = -Tx / m^2;
A_jac(4,5) = -Ty / m^2;
B_jac(3,1) = 1/m;
B_jac(4,2) = 1/m;
B_jac(5,1) = -Vc * Tx / Tmag_reg;
B_jac(5,2) = -Vc * Ty / Tmag_reg;
```

- Riga 265: legge massa e spinta dallo stato/controllo.
- Riga 266: **la regolarizzazione**. `d/dT (-Vc*|T|) = -Vc * T/|T|` non esiste in
  `T = 0`. Sostituendo `|T| -> sqrt(|T|^2 + eps^2)` con `eps^2 = 1e-6` (quindi
  `eps = 1e-3` non-dim, circa 20 N) la derivata resta finita e tende a zero
  linearmente vicino all'origine. **Il costo:** una distorsione `O(eps)` proprio dove
  `|T| ~ eps` -- cioe' **sull'arco di coast**, che e' esattamente dove il modello
  lineare fatica di piu'. E' l'hack piu' rilevante del file, ed e' proprio il motivo
  per cui la variante (d) (che elimina quella riga singolare per costruzione) e'
  cosi' piu' pulita.
- Righe 267-271: `A = df/dx`. La quinta riga e' **identicamente nulla**: la portata di
  massa non dipende dallo stato. Le righe 3-4 catturano l'unico accoppiamento
  stato-stato non banale: la sensitivita' dell'accelerazione alla massa.
- Righe 272-276: `B = df/du`. Le righe 3-4 sono `1/m`; la riga 5 e' `-Vc` per il
  versore di spinta.

Nessuna `arguments` validation: sta dentro il loop di `ode45` chiamato da
`compute_ltv_zoh`, quindi viene invocata migliaia di volte per intervallo.

## `ltv_aug_rhs` (righe 279-303)

E' l'ODE aumentata dell'Appendice A, nella forma **beta-gamma** (riferita a `t_k`).
Lo stato aumentato `z` e' un vettore 70x1:

    z(1:5)    = x_ref     (5)     traiettoria di riferimento
    z(6:30)   = vec(Phi)  (25)    matrice di transizione
    z(31:55)  = vec(Psi)  (25)    transizione inversa, Psi = Phi^{-1}
    z(56:65)  = vec(Beta) (10)    integrale di Psi*B     (5x2)
    z(66:70)  = Gamma     (5)     integrale di Psi*c

```matlab
[A_jac, B_jac] = jacobians(x_ref, u_k, Vc);
f_val = ode_descent(x_ref, u_k, Vc);
c_off = f_val - A_jac * x_ref - B_jac * u_k;
dx_ref =  f_val;
dPhi   =  A_jac * Phi;
dPsi   = -Psi * A_jac;
dBeta  =  Psi * B_jac;
dGamma =  Psi * c_off;
```

- Righe 289-291: spacchetta lo stato aumentato.
- Riga 294: **le Jacobiane sono rivalutate a ogni passo interno di `ode45`**, sulla
  `x_ref` che sta evolvendo dentro l'intervallo (riga 297: `dx_ref = f_val`). Quindi
  `A(t)` e `B(t)` sono genuinamente **tempo-varianti dentro l'intervallo**, non
  congelate al valore nodale. E' la lettura "letterale" dell'Appendice A -- alcune
  implementazioni si accontentano di congelare `A`, `B` a `t_k` (e allora la ZOH
  diventa un `expm` per intervallo), qui no.
- Riga 296: il **termine affine** `c = f(x_ref, u_ref) - A x_ref - B u_ref` derivato
  sopra. E' cio' che rende la linearizzazione esatta nel punto di riferimento.
- Riga 299: `dPsi = -Psi * A`, con la derivazione fatta sopra.
- Righe 300-301: quadrature pure, nessun feedback.

**Sottigliezza da capire (e da dire all'orale):** `x_ref` dentro l'intervallo e' la
**propagazione nonlineare** dello stato nodale di riferimento sotto `u_k`. Se il
riferimento non e' dinamicamente consistente con la ZOH (e la soluzione trapezoidale
del warm start **non lo e'**), allora `x_ref(t_{k+1})` calcolato qui **non coincide**
con `ref.x(k+1)`. La linearizzazione e' quindi costruita attorno a una traiettoria
leggermente diversa da quella nodale. E' una delle ragioni per cui i primi passi SCvx
vengono rifiutati.

## `compute_ltv_zoh` (righe 305-339)

- Riga 322: `dt = tf/(N-1)`.
- Riga 326: **tolleranze `RelTol 1e-8`, `AbsTol 1e-10`** -- piu' lasche di quelle del
  replay (1e-10 / 1e-12, riga 839). Ragionevole: qui si sta costruendo un modello che
  verra' comunque rifatto alla prossima iterazione SCvx, li' si sta misurando la
  verita'.
- Riga 331: condizioni iniziali dell'ODE aumentata: `Phi(t_k) = Psi(t_k) = I`,
  `Beta = 0`, `Gamma = 0`.
- Riga 332: **una `ode45` per intervallo** -- quindi 49 integrazioni di un sistema 70-D
  per **ogni** iterazione SCvx. E' il collo di bottiglia computazionale delle varianti
  (b) e (c) (e la ragione per cui la (d), che fa un solo `expm`, e' ~10x piu' veloce).
- Righe 334-337: lettura dei valori terminali e ricostruzione

      Abar_k = Phi(t_{k+1})
      Bbar_k = Phi(t_{k+1}) * Beta(t_{k+1})
      cbar_k = Phi(t_{k+1}) * Gamma(t_{k+1})

  E' il prodotto finale imposto dall'identita' di semigruppo.

## `solve_ltv_nlp` (righe 341-412) e `ltv_nonlcon` (righe 414-427)

E' il **sottoproblema interno** della variante (b), risolto con `fmincon`.

- Righe 362-363: `nz = 7N`, con il layout per nodo `[x; y; vx; vy; m; Tx; Ty]`.
  Stesso layout di tutte le altre varianti fmincon (e di `main_task1.m`).
- Righe 364-368: warm start = il riferimento corrente (`ref_to_z`), oppure la guess
  lineare se non c'e' riferimento.
- Riga 369: box bounds con `zero_uN = true` (il controllo del nodo N e' inchiodato a 0:
  sotto ZOH non ha nessun intervallo su cui agire).
- Righe 370-372: intersezione con la trust region.
- Righe 375-405: **assemblaggio della matrice di uguaglianza sparsa via triplette**
  `(rows, cols, vals)`. Le righe di `Aeq` sono `9 + 5*(N-1)`:
  - righe 1-5: stato iniziale completo (`x0, y0, vx0, vy0, m0`);
  - righe 6-9: posizione e velocita' finali nulle (**la massa finale e' libera**: e'
    proprio la variabile che si massimizza);
  - righe 10 in poi, a blocchi di 5: la **dinamica LTV**
    `x_{k+1} - Abar_k x_k - Bbar_k u_k = cbar_k`. Ogni riga ha 8 entrate non nulle:
    1 (per `x_{k+1}`) + 5 (per `A`) + 2 (per `B`). Da qui `nnz_dyn = 8*5*(N-1)`
    (riga 376).
- Riga 405: `Aeq = sparse(...)`, poi riga 408 la converte in `full(Aeq)` per `fmincon`.
  Costruire in sparso e passare in denso e' un po' contorto (il commento a riga 374 dice
  che serviva a evitare un warning di indicizzazione sparsa); con `N = 50` la matrice e'
  `254 x 350`, quindi non e' un problema.
- Righe 407-409: obiettivo `-z(iN_m)` con `iN_m = (N-1)*7 + 5`, cioe' **massimizza la
  massa al nodo N**. E' un obiettivo **lineare**. Vincoli: uguaglianze lineari
  (dinamica + BC), box (bounds + trust), e `ltv_nonlcon` per i vincoli di percorso.
- `ltv_nonlcon` (righe 414-427): **solo `path_ineq`**, `c_eq = []`. La dinamica non e'
  piu' qui, e' gia' dentro `Aeq`. Questo e' il salto qualitativo rispetto a `solve_zoh`:
  la dinamica e' passata da vincolo **nonlineare** a vincolo **lineare**.

> **Possibile domanda d'esame** -- Il sottoproblema di `solve_ltv_nlp` e' convesso.
> Perche' allora e' dato a `fmincon`, che e' un solver per NLP non convessi?
> *Risposta:* funziona (fmincon trova lo stesso ottimo globale, il problema e'
> convesso), ma e' uno spreco: SQP approssima il cono di spinta `|T| <= Tmax` con
> successive linearizzazioni quadratiche invece di trattarlo come primitiva, e paga una
> Jacobiana alle differenze finite. La variante (c) e' proprio la stessa cosa fatta
> bene: stesso loop esterno, sottoproblema riscritto come SOCP e dato a ECOS -- circa
> 3x piu' veloce end-to-end e con fedelta' per-passo un ordine di grandezza migliore.

## `solve_scvx` (righe 429-537)

Il **loop esterno SCvx** con fmincon. Il docstring alle righe 442-447 riassume esattamente
la logica del ratio test.

- Righe 458-462: riferimento iniziale = `init_ref` (il trapezoidale) o, se assente, la
  guess lineare. Nella pratica e' sempre il trapezoidale (riga 80).
- Righe 464-468: parametri del ratio test.

  ```matlab
  rho     = 1.0;     rho_min = 1e-3;   rho_max = 1.0;
  eta_l   = 0.25;    eta_h   = 0.7;
  ```

  Nota `rho_max = 1.0`: la trust region puo' solo **restringersi** rispetto ai raggi
  base, mai allargarsi oltre. E' una scelta conservativa coerente col fatto che i raggi
  base sono gia' generosi.
- Righe 470-474: `conv_hist` traccia per-iterazione `m_f`, `delta_x`, `rho`, `eta`, e
  se il passo e' stato accettato.
- Righe 484-485: ricalcola le matrici LTV **attorno al riferimento corrente** e risolve
  il sottoproblema. Questo e' il "successive" di successive convexification.
- Righe 487-496: il **ratio test**.

  ```matlab
  J_pred = sol_cand.m_f - ref.m_f;      % guadagno predetto (LTV)
  [~, X_act] = fwd_integrate(sol_cand, d, 'zoh');
  m_f_actual = X_act(end, 5);
  J_act      = m_f_actual - ref.m_f;    % guadagno reale (ode45)
  eta = J_act / J_pred;                 % (con guardia se J_pred ~ 0)
  ```

  Il costo e' `-m_f`, quindi "riduzione del costo" = "aumento di `m_f`": il commento a
  riga 487 lo dice.

  **Onesta':** numeratore e denominatore usano **la stessa baseline** `ref.m_f`, che e'
  la massa finale *secondo l'NLP* del riferimento -- non la sua massa finale *replayata*.
  Poiche' il riferimento e' a sua volta un candidato accettato la cui `m_f` NLP differisce
  (di poco) dal suo replay, `eta` porta dentro l'errore di trascrizione del riferimento.
  Un ratio test rigoroso userebbe `m_f_replay(ref)` come baseline in entrambi. L'effetto
  e' piccolo ma non nullo, e spiega parte del comportamento a scacchiera di `eta`
  (accettati/rifiutati alternati) riportato nel report.
- Righe 498-501: `delta_x` = norma della variazione dell'iterato di **stato** (posizione,
  velocita', massa). **I controlli non entrano** nel criterio di convergenza: due
  soluzioni con lo stesso stato ma spinte diverse (possibile solo in coast) sarebbero
  considerate convergenti. In pratica non succede.
- Riga 503: `accepted = (eta >= eta_l)`. Notare: `eta` puo' essere **negativo** (il
  modello predice un guadagno e la realta' consegna una perdita) -- in quel caso e'
  giustamente rifiutato.
- Righe 513-530: l'aggiornamento.
  - accettato: `sol_best = sol_cand; ref = sol_cand;` e, se `eta > eta_h`, `rho`
    raddoppia (cappato a `rho_max`). Convergenza se `delta_x < tol`.
  - rifiutato: **il riferimento NON viene aggiornato** (questo e' essenziale: si
    ri-risolve lo stesso sottoproblema con una trust region piu' stretta), `rho`
    dimezza. Se `rho < rho_min` si esce.
- Riga 535: `sol = sol_best`. **Attenzione al nome:** `sol_best` e' l'**ultimo passo
  accettato**, non quello con la `m_f` migliore. Se un passo accettato peggiorasse
  `m_f` (possibile: `eta >= 0.25` non implica `J_act > 0` se `J_pred < 0`), verrebbe
  comunque tenuto. Il nome e' fuorviante.

## `ternary` (righe 539-548)

```matlab
function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
```

Innocuo di per se', ma **attenzione a come e' usato**. Alla riga 1091 (e 1253):

    eta = ternary(abs(J_pred) < 1e-10, 1, J_act / J_pred);

MATLAB valuta **entrambi** gli argomenti prima di chiamare la funzione: `J_act/J_pred`
viene calcolato **anche** quando `J_pred ~ 0`, producendo `Inf` o `NaN` (che poi viene
scartato). Non e' un bug -- MATLAB non solleva errore sulla divisione per zero -- ma la
guardia **non guardia' un bel niente**. `solve_scvx` (righe 492-496) usa invece un
`if/else` vero, che e' corretto. Incoerenza fra le tre implementazioni dello stesso loop.

## `ref_to_z` (righe 550-569) e `apply_trust` (righe 571-600)

- `ref_to_z`: impacchetta la struct di riferimento nel vettore di decisione `7N x 1`.
  **Piccola verruca:** copia anche `ref.Tx(N)`, `ref.Ty(N)` (righe 566-567). Ma il
  riferimento e' la soluzione **trapezoidale**, che NON inchioda il controllo del nodo N
  a zero. Poiche' `box_bounds(..., zero_uN=true)` impone `lb = ub = 0` su quelle due
  variabili, il punto iniziale passato a `fmincon` e' **fuori dai bound**, e `fmincon` lo
  proietta silenziosamente dentro. Innocuo ma sporco.
- `apply_trust` (righe 581-599): intersezione **in norma infinito** dei bound esistenti
  con `[ref - raggio, ref + raggio]`. La trust region e' quindi un **box**, non una palla
  euclidea, e non e' una penalita' soft. E' la variante "adaptive hard trust region" della
  tassonomia SCvx. I controlli sono vincolati solo per `i < N` (riga 593), coerentemente
  con la convenzione ZOH.

## `solve_trap` (righe 602-624) e `solve_zoh` (righe 626-650)

Stessa struttura, differiscono solo nel `nonlcon` e nella convenzione sul controllo
terminale.

- `solve_trap` (baseline Task 1): difetti trapezoidali, controllo piecewise-linear
  (tutti gli `N` controlli sono significativi -> `zero_uN = false`).
- `solve_zoh` (variante a): difetti ZOH via RK4, controllo piecewise-constant ->
  `zero_uN = true`, `u_N` inchiodato a zero. Il vettore di decisione ha la stessa
  dimensione `7N`: si spreca un controllo, ma il layout resta identico e riusabile.
- Entrambi: obiettivo `-z(iN_m)` (massimizza `m_N`), `Aeq/beq` dalle `bcs`, opzioni da
  `fmincon_opts()` con `Display = 'final'`. `solve_ltv_nlp` invece usa
  `fmincon_opts('off')` (riga 409), per non stampare 15 volte dentro il loop SCvx.

## `init_guess` (righe 652-676) e `box_bounds` (righe 678-700)

- `init_guess`: interpolazione **lineare** fra stato iniziale e origine per posizione e
  velocita'; massa che scende linearmente del 30% (riga 668); controllo `Tx = 0`,
  `Ty = m0` -- cioe' **spinta di hover** (riga 673), che in non-dim vale esattamente
  `m*g = 1*1 = 1` (ecco perche' il codice scrive `d.m0`).
- `box_bounds`: `y >= 0` (non si scava), massa in `[1e-3, m0]` (**strettamente
  positiva**: `m` compare a denominatore in `ode_descent`, `m = 0` farebbe esplodere
  tutto), componenti di spinta in `[-Tmax, Tmax]`. Nota che questo box e' un
  **quadrato** che contiene il disco `|T| <= Tmax`: il vincolo circolare vero e'
  imposto separatamente in `path_ineq`. Il box e' ridondante ma aiuta il solver a
  restare in una regione limitata.

## `bcs` (righe 702-720)

9 uguaglianze lineari: 5 per lo stato iniziale completo, 4 per posizione e velocita'
finali nulle. **La massa finale non e' vincolata** -- e' l'obiettivo.

## `fmincon_opts` (righe 722-736)

```matlab
opts = optimoptions('fmincon', ...
    'Algorithm', 'sqp', 'Display', display_mode, ...
    'MaxIterations', 1000, 'MaxFunctionEvaluations', 1e6, ...
    'OptimalityTolerance', 1e-5, 'ConstraintTolerance', 1e-6, ...
    'StepTolerance', 1e-10);
```

- `sqp`: buono per NLP con vincoli densi e problemi di dimensione media come questo.
- `MaxFunctionEvaluations = 1e6`: serve, perche' **le Jacobiane sono alle differenze
  finite**. Con 350 variabili, ogni gradiente costa ~350 valutazioni del `nonlcon`, e
  ogni valutazione del `nonlcon` della variante (a) costa 49 propagazioni RK4. E' il
  costo dominante e la ragione per cui la variante (b) e' lenta pur avendo dinamica
  lineare (la nonlinearita' residua in `path_ineq` obbliga comunque fmincon a fare
  differenze finite sull'intero vettore).
- Il README ammette che le run del Task 1 a `tf` corto arrivano al cap di 1000
  iterazioni con ottimalita' del primo ordine ferma a `1e-3`-`1e-4`: la causa e' proprio
  il punto **non differenziabile** `T = 0` della norma `|T|` in `path_ineq`, che cade
  esattamente sull'arco di coast.

## `unpack` (righe 738-757)

Da vettore di decisione a struct. Calcola `sol.Tmag`, `sol.m_f = sol.m(end)`.

## `trap_nonlcon` (righe 759-781) e `zoh_nonlcon` (righe 783-803)

Le due trascrizioni "a difetti":

- **trapezoidale** (riga 777):

      def_k = x_{k+1} - x_k - (dt/2) * (f_k + f_{k+1})

  Regola del trapezio: ordine `O(dt^2)` globale. Il controllo e' implicitamente
  piecewise-linear (entra sia in `f_k` che in `f_{k+1}`).

- **ZOH/RK4** (riga 798):

      def_k = x_{k+1} - RK4(x_k, u_k, dt, Vc, n_sub)

  Qui il "difetto" e' la differenza rispetto al **flusso esatto** (approssimato a
  ordine `O(h^4)` con `h = dt/n_sub`). E' un **multiple shooting**, non collocazione.
  Con `n_sub = 2`, `h = dt/2 = 0.0222` non-dim, l'errore globale e' ~`1e-8` -- e infatti
  il replay `ode45` della variante (a) misura `1.4e-8` (README, riga 110). I due numeri
  **coincidono per costruzione**: il difetto ZOH *e'* la mappa di flusso RK4, quindi
  l'errore di trascrizione e' esattamente l'errore RK4-vs-ode45.

## `path_ineq` (righe 805-820)

```matlab
Tmag     = sqrt(Z(6,:).^2 + Z(7,:).^2).';
g_thr_lo = d.Tmin - Tmag;          % Tmin - |T| <= 0
g_thr_hi = Tmag - d.Tmax;          % |T| - Tmax <= 0
tt       = tan(d.theta_mx);
g_gs_pos = ( Z(1,:).' - tt*Z(2,:).');   %  x - tan(th)*y <= 0
g_gs_neg = (-Z(1,:).' - tt*Z(2,:).');   % -x - tan(th)*y <= 0
```

- Il glide-slope: la coppia di semipiani equivale a `|x| <= tan(theta_max) * y`, cioe' il
  veicolo deve restare **dentro un cono** con vertice sulla piazzola e semiapertura
  `theta_max = 60 deg` misurata **dalla verticale**. E' un vincolo **convesso** (un cono
  di Lorentz in 2D degenerato in due semipiani). Il punto iniziale ha
  `x0/y0 = 1000/3000 = 0.33`, ben dentro `tan(60) = 1.73`.
- `g_thr_lo` con `Tmin = 0` diventa `-|T| <= 0`: **sempre soddisfatto**. E' un vincolo
  ridondante -- ma non gratuito: e' **non differenziabile in `T = 0`**, e siccome
  l'arco di coast sta esattamente li', fmincon ci sbatte contro con le sue differenze
  finite. Il README lo identifica come la causa radice della convergenza lenta.

## `fwd_integrate` (righe 822-857)

Il replay di validazione. Per ogni intervallo, integra la dinamica **nonlineare vera**
con `ode45` a tolleranze strette (riga 839: `1e-10 / 1e-12`), partendo **sempre dallo
stato replayato precedente** `X(k,:)` -- non dallo stato NLP. E' quindi un replay
**open-loop puro**: gli errori si accumulano, ed e' esattamente cio' che si vuole
misurare (dove atterrerei davvero?).

Due modalita' (righe 842-851): `'zoh'` (controllo costante) e `'pwl'` (interpolazione
lineare fra `u_k` e `u_{k+1}`). Scegliere la modalita' **coerente con la trascrizione**
e' obbligatorio, altrimenti si misura il mismatch di convenzione, non la fedelta'.

## `node_err` (righe 859-868)

```matlab
e = vecnorm([sol.x sol.y sol.vx sol.vy] - X(:,1:4), 2, 2);
```

Norma euclidea dell'errore su **posizione e velocita'**, per nodo, in unita' non-dim.
**La massa e' esclusa** (il docstring alla riga 860 lo dice). Ma il messaggio stampato
alla riga 130 dice "max grid-node nondim **state** error": e' una piccola imprecisione,
lo stato completo avrebbe 5 componenti. Il canale della massa e' monitorato
separatamente, tramite il drift di `m_f` (riga 144) e tramite il ratio test.

## `plot_compare3` (righe 870-963)

Percorso di fallback quando YALMIP/ECOS non ci sono. Quattro figure (traiettoria con
corridoio glide-slope, `|T|` con `stairs` per le ZOH e `plot` per la trapezoidale --
scelta corretta, la ZOH e' davvero una scala --, massa, fedelta' in scala semilog) piu'
il pannello di convergenza SCvx (righe 928-962: `delta_x`, `rho` e `m_f` sopra, barre di
`eta` con le soglie `0.25` / `0.7` sotto, e i passi rifiutati colorati in arancione).

---

## `solve_ltv_nlp_yalmip` (righe 965-1040)

Il sottoproblema interno della variante **(c)**: lo stesso problema di `solve_ltv_nlp`,
ma scritto come **SOCP** e dato a ECOS.

- Righe 988-989: le variabili di decisione **modellate correttamente**:

  ```matlab
  X = sdpvar(5, N,   'full');   % stato per nodo
  U = sdpvar(2, N-1, 'full');   % controllo per INTERVALLO
  ```

  Notare `U` ha `N-1` colonne, non `N`: sotto ZOH esistono solo `N-1` controlli. E' piu'
  pulito della versione fmincon, che porta un `u_N` fantasma inchiodato a zero. Il
  padding a `N` avviene solo in uscita (riga 1032).
- Righe 993-994: BC iniziali (5) e finali (4, massa libera).
- Righe 996-998: **dinamica LTV come uguaglianze lineari**.
- Righe 1000-1002: **il vincolo conico**:

  ```matlab
  cstr = [cstr, norm(U(:,k)) <= d.Tmax];
  ```

  Questo e' un **cono di Lorentz** (SOC). YALMIP lo riconosce come tale e lo passa a
  ECOS in forma nativa. Il commento (riga 1000) nota giustamente che con `Tmin = 0` il
  bound inferiore e' banalmente soddisfatto e non serve modellarlo.
- Righe 1004-1007: glide-slope, lineare.
- Righe 1009-1010: quota `>= 0`, massa in `[1e-3, m0]`.
- Righe 1012-1024: trust region, **lineare** (box in norma infinito). Anche la trust
  region resta dentro la classe SOCP.
- Riga 1026:

  ```matlab
  res = optimize(cstr, -X(5,N), ...
                 sdpsettings('solver','ecos','verbose',0));
  ```

  Obiettivo **lineare** (`-m_N`), vincoli lineari + coni: e' un SOCP puro. ECOS e' un
  solver interior-point per coni (lineari, di secondo ordine, esponenziali): risolve in
  tempo polinomiale con **certificato di ottimalita' globale del sottoproblema**.
- Righe 1027-1029: se `res.problem ~= 0` stampa un warning ma **prosegue comunque**, e
  `value(X)` restituisce l'ultimo iterato. E' quello che salva il loop nelle iterazioni
  in cui il sottoproblema linearizzato risulta infeasible (perche' mancano i virtual
  control). Fragile, ma funzionante.
- Righe 1031-1039: spacchetta e ricostruisce la struct nel formato standard.

**Nota di performance:** i vincoli sono accumulati con `cstr = [cstr, ...]` dentro dei
loop (crescita quadratica) e non c'e' nessun `yalmip('clear')` fra una chiamata e
l'altra, quindi il contatore interno di `sdpvar` cresce per tutte le 15 iterazioni. Con
`N = 50` non e' un problema (29 s end-to-end contro 82 s della variante fmincon), ma non
scalerebbe.

## `solve_scvx_yalmip` (righe 1042-1121)

**Copia carbone** di `solve_scvx`, con la sola sostituzione della riga 485
(`solve_ltv_nlp`) con la riga 1086 (`solve_ltv_nlp_yalmip`). Stesso ratio test, stessi
`eta_l = 0.25` / `eta_h = 0.7`, stesso warm start trapezoidale, stessi raggi di trust
region. E' proprio questo che rende il confronto (b) vs (c) **pulito**: l'unica variabile
che cambia e' il solver interno.

Duplicazione di codice notevole (~80 righe quasi identiche). Un refactor con un handle
al solver interno passato come argomento sarebbe stato piu' elegante; il codice ha
scelto la chiarezza espositiva.

Unica differenza sostanziale: riga 1091 usa `ternary` (con il problema di valutazione
eager discusso sopra) invece dell'`if/else` di `solve_scvx`.

## `solve_gfold_socp` (righe 1123-1196)

Il sottoproblema interno della variante **(d)**. E' il pezzo teoricamente piu' denso.

- Righe 1147-1148: le nuove variabili.

  ```matlab
  XI = sdpvar(5, N,   'full');   % [x; y; vx; vy; z],  z = ln(m)
  W  = sdpvar(3, N-1, 'full');   % [ux; uy; sigma]
  ```

  Il controllo ha **tre** componenti: le due dell'accelerazione piu' lo **slack**
  `sigma`. Lo slack e' una variabile di decisione a tutti gli effetti.
- Riga 1150: `z0 = log(d.m0)` che vale **esattamente 0**, perche' `m0_nd = 1`.
- Righe 1152-1153: BC. Notare il commento: **`z_N` e' libero** -- e' l'obiettivo.
- Riga 1156: **dinamica LTI**:

  ```matlab
  cstr = [cstr, XI(:,k+1) == Abar*XI(:,k) + Bbar*W(:,k) + cbar];
  ```

  `Abar`, `Bbar`, `cbar` sono **le stesse per tutti i `k`** (senza indice!): il sistema
  e' tempo-invariante. Calcolate una volta sola alla riga 1218 con `lti_zoh`.

  Vale la pena esplicitare cosa sono davvero. In `lti_zoh.m`, `A` e' **nilpotente**
  (`A^2 = 0`), quindi `expm(A*s) = I + A*s` **esattamente** e gli integrali si fanno a
  mano:

      Abar = I + A*dt              -> doppio integratore ZOH
      Bbar = (I*dt + A*dt^2/2) * B
      cbar = (I*dt + A*dt^2/2) * c

  che, componente per componente, danno

      x_{k+1}  = x_k  + dt*vx_k + (dt^2/2)*ux_k
      y_{k+1}  = y_k  + dt*vy_k + (dt^2/2)*(uy_k - 1)
      vx_{k+1} = vx_k + dt*ux_k
      vy_{k+1} = vy_k + dt*(uy_k - 1)
      z_{k+1}  = z_k  - Vc*dt*sigma_k

  Cioe': **le formule del moto uniformemente accelerato**, esatte, piu' una massa che
  decade esponenzialmente in modo esatto. Non c'e' **nessun errore di discretizzazione**.
  Ecco perche' la fedelta' del replay e' `7e-12` (il floor dell'integratore) invece di
  `1e-4`-`1e-8`.

- Riga 1157: **il cono lossless**:

  ```matlab
  cstr = [cstr, norm(W(1:2,k)) <= W(3,k)];   % |u| <= sigma
  ```

- Righe 1158-1159: **il bound superiore linearizzato** -- l'unica non convessita'
  residua:

  ```matlab
  ezr = exp(-z_ref(k));
  cstr = [cstr, W(3,k) <= d.Tmax*ezr*(1 - (XI(5,k) - z_ref(k)))];
  ```

  cioe' `sigma_k <= Tmax * e^{-z_ref,k} * (1 - (z_k - z_ref,k))`, la tangente di
  `Tmax*e^{-z}` in `z_ref`. Sotto-stima il bound vero -> **conservativa** -> ogni
  iterato accettato e' feasible per il problema originale. E' l'**unico** motivo per cui
  serve ancora un loop SCvx.
- Righe 1162-1166: glide-slope (lineare) e bound sulla massa **in log**:
  `log(1e-3) <= z <= 0`. **Fragile:** l'upper bound `0` e' hard-coded assumendo
  `m0_nd = 1`; sarebbe piu' robusto scrivere `<= z0`.
- Righe 1168-1181: trust region, con raggi **su variabili diverse**: `pos`, `vel`, `lz`
  (in log-massa!), `u` (accelerazione), `sig`. Non sono i raggi della riga 49 -- sono
  quelli della riga 1231.
- Riga 1183: `optimize(cstr, -XI(5,N))` -> **massimizza `z_N`**, che per monotonia di
  `exp` equivale a massimizzare `m_N`.
- Righe 1184-1186: **`res.problem ~= 0 && res.problem ~= 1`**. Il flag `1` di YALMIP e'
  "infeasible problem": il codice lo **sopprime deliberatamente** e prosegue con
  l'iterato che ECOS restituisce. E' l'ammissione implicita che manca il virtual
  control. Da dichiarare all'orale: e' l'hack piu' discutibile del file.
- Righe 1188-1195: **back-transform** alle variabili originali:

      m = exp(z),      T = m * u

  con l'ultimo nodo paddato a zero.

> **Possibile domanda d'esame** -- Se la dinamica GFOLD e' esatta, perche' serve ancora
> un loop SCvx?
> *Risposta:* perche' resta una non convessita': il bound superiore di spinta
> `|u| <= Tmax * e^{-z}`. La funzione `e^{-z}` e' convessa, quindi vincolare `sigma`
> **sotto** di essa definisce un insieme non convesso. Lo si linearizza attorno a
> `z_ref` e si itera. Ma e' **un solo vincolo scalare per nodo**, non l'intera dinamica:
> per questo il loop converge in 3 iterazioni invece delle 15 delle varianti (b)/(c), e
> ogni passo e' accettato con `eta ~= 1`.

## `solve_gfold_scvx` (righe 1198-1283)

Il loop esterno della variante (d). Struttura identica alle altre due, con **tre
differenze sostanziali**:

1. **Riga 1218:** `[Abar, Bbar, cbar] = lti_zoh(tf/(N-1), d.Vc)` -- **una volta sola**,
   fuori dal loop. Un solo `expm` 9x9 contro le `15 iterazioni x 49 intervalli = 735`
   integrazioni `ode45` 70-D delle varianti (b)/(c). E' da qui che viene il fattore ~10
   sul tempo di calcolo (5 s contro 82 s / 29 s).

2. **Righe 1221-1229: e' self-starting.** Nessun warm start trapezoidale. Il riferimento
   iniziale e':

   ```matlab
   m_apri = max(d.m0 - d.Vc*d.Tmax*t_grid, 1e-2);
   ref.z  = log(m_apri);
   ```

   cioe' il **profilo di massa analitico a spinta massima costante**:
   `m(t) = m0 - Vc*Tmax*t` (integrale banale di `m_dot = -Vc*Tmax`). E' il profilo di
   consumo **piu' pessimistico possibile** -- serve solo come punto attorno a cui
   linearizzare il bound di spinta, e siccome quel bound e' l'unica cosa linearizzata, un
   riferimento grossolano basta. Posizione e velocita' sono interpolate linearmente e
   `ref.ux/uy/sig` sono valori segnaposto (`0, 1, 1`) che di fatto non vengono mai usati
   (all'iterazione 1 non c'e' trust region; dalla 2 in poi `ref` e' gia' il candidato
   accettato).

3. **Righe 1240-1248: la prima iterazione e' risolta SENZA trust region.**

   ```matlab
   if iter == 1
       cand = solve_gfold_socp(tf,N,d,Abar,Bbar,cbar,ref.z,[],[]);
   else
       ...
   end
   ```

   Il commento (righe 1241-1242) lo motiva: la dinamica e' esatta, quindi lasciare il
   SOCP libero di trovare una traiettoria dinamicamente feasible e' sicuro. E' l'opposto
   di quello che si farebbe con dinamica linearizzata (dove il primo passo senza trust
   region e' esattamente quello che esplode).

- Riga 1231: i raggi base **specifici di GFOLD**:
  `pos 0.5, vel 1.0, lz 0.4, u 4.0, sig 4.0` -- molto piu' larghi di quelli delle
  varianti (b)/(c), coerentemente col fatto che qui non c'e' errore di modello dinamico
  da cui difendersi. `lz = 0.4` in log-massa corrisponde a un fattore `e^0.4 = 1.49`.
- Riga 1251: il ratio test usa **`fwd_integrate_uacc`** (replay ad accelerazione
  costante), non `fwd_integrate`.
- **Onesta' sull'iterazione 1:** `ref.m_f = exp(ref.z(N))` con il profilo a spinta
  massima da' circa `0.40` non-dim, mentre l'ottimo vero e' `~0.70`. Quindi al primo
  passo sia `J_pred` sia `J_act` sono grandi e positivi e `eta ~= 1` **quasi per
  costruzione**: il ratio test della prima iterazione non sta davvero verificando nulla.
  Diventa informativo dalla seconda in poi.

## `fwd_integrate_uacc` (righe 1285-1303)

Gemello di `fwd_integrate`, ma il RHS e' `ode_descent_uacc` e il controllo tenuto
costante e' `u = T/m`, non `T`. La dinamica replayata e':

    vx_dot = ux
    vy_dot = uy - 1
    m_dot  = -Vc * m * |u|          <- NOTA: dipende da m!

Con `u` costante, le righe di velocita' e posizione sono **integrabili esattamente** in
forma polinomiale (moto uniformemente accelerato), e la massa decade
**esponenzialmente**: `m(t) = m_k * exp(-Vc*|u|*(t - t_k))`, cioe'
`z(t) = z_k - Vc*|u|*(t - t_k)` -- **lineare in `z`**. Che e' esattamente cio' che la
riga discreta `z_{k+1} = z_k - Vc*dt*sigma_k` prevede, **a patto che `sigma_k = |u_k|`**.

Questo chiude il cerchio: **la fedelta' `7e-12` del replay GFOLD e' anche una verifica
numerica indipendente che il cono lossless e' attivo.** Se `sigma_k > |u_k|`, il SOCP
farebbe decadere `z` piu' in fretta di quanto il replay consegni, e l'errore sarebbe
visibile. Non lo e'. Bel risultato da citare all'orale.

## `plot_compare5` (righe 1305-1420)

Come `plot_compare3` piu' le varianti YALMIP (viola) e GFOLD (verde acqua). Le figure di
convergenza SCvx sono **tre**, generate in un ciclo (righe 1377-1419), una per ciascun
loop SCvx, cosi' da poterle confrontare a colpo d'occhio.

---

## Cosa misura davvero il confronto, e cosa emerge

Lo script produce **quattro** metriche per ogni trascrizione:

1. **Massa finale `m_f`** (quindi carburante) -- e' l'ottimalita'.
2. **Errore massimo ai nodi** (`node_err`) -- e' la **fedelta' di trascrizione**: quanto
   la traiettoria discretizzata si discosta da un'integrazione `ode45` degli stessi
   controlli.
3. **Dispersione al touchdown** (posizione, velocita', drift di massa) dal replay
   open-loop -- e' la metrica **operativa**.
4. **Wall time e numero di iterazioni SCvx**.

I numeri (dal report, Tabelle `tab:hm2_compare4` e `tab:hm2_replay`, `tf = 38 s`,
`N = 50`):

| Trascrizione            | m_f [kg] | Err. max nodi | Pos. touchdown | Wall [s] | Iter |
|-------------------------|----------|---------------|----------------|----------|------|
| Trapezoidale (PWL)      | 1403.20  | 1.55e-3       | 4.3 m          | 90       | --   |
| ZOH nonlineare + RK4    | 1403.37  | 1.40e-8       | 6.1e-6 m       | 100      | --   |
| LTV + SCvx (fmincon)    | 1399.75  | 4.01e-4       | 0.18 m         | 82       | 15   |
| LTV + SCvx (YALMIP/ECOS)| 1400.84  | 4.09e-5       | 0.019 m        | 29       | 15   |
| GFOLD log-mass          | 1402.77  | 7.30e-12      | 2.1e-8 m       | 5        | 3    |

Cosa emerge, in ordine di importanza:

- **Tutte e cinque trovano la stessa soluzione fisica.** `m_f` varia di ~4 kg su 1400,
  cioe' meno dello 0.3%. Le differenze sono sotto la precisione ingegneristica del
  modello (niente aerodinamica, niente assetto). Tutte esibiscono la struttura
  **max-coast-max** (bang-off-bang) prevista dal PMP.
- **Fedelta' e ordine di convergenza sono cose diverse dall'ottimalita'.** La
  trapezoidale ha `m_f` fra le migliori ma la fedelta' peggiore (`1e-3`): il suo
  "ottimo" e' calcolato su una dinamica che vale solo a `O(dt^2)`, quindi e'
  leggermente **ottimistico**. La ZOH-RK4 e' `O(h^4)` e va a `1e-8`.
- **Esporre la struttura conica paga.** Varianti (b) e (c) sono **matematicamente lo
  stesso sottoproblema**: (c) e' 3x piu' veloce e ha fedelta' un ordine di grandezza
  migliore, solo perche' il cono `|T| <= Tmax` e' passato a ECOS come primitiva invece
  di essere approssimato da SQP con Jacobiane alle differenze finite.
- **Il cambio di variabili batte l'iterazione numerica.** GFOLD e' l'unica variante che
  attacca la **causa** della non convessita' invece di gestirne gli **effetti**:
  dinamica esatta -> nessun errore di modello -> ogni passo accettato con `eta ~= 1` ->
  3 iterazioni contro 15 -> 5 s contro 82 s -> fedelta' al floor dell'integratore.
- **Entrambi i loop LTV+SCvx (b)/(c) NON convergono**: arrivano al cap di 15 iterazioni
  con `||delta x|| ~ 1e-2`, sopra la tolleranza `1e-3`. Il report lo dice apertamente.
  Non e' un disastro (ogni passo accettato ha superato il ratio test contro `ode45`,
  quindi l'iterato restituito e' realizzabile), ma non e' convergenza.

---

## Limiti, hack e deviazioni dalla teoria -- da dichiarare

Elenco onesto, in ordine di gravita':

1. **Nessun virtual control, nessun buffer di slack sui vincoli di percorso.** E' la
   deviazione piu' seria dalla ricetta SCvx canonica. Il codice riconosce
   l'infeasibility del sottoproblema (righe 410, 1028) ma si limita a un warning; la
   variante GFOLD la **sopprime esplicitamente** (riga 1184: `res.problem ~= 1`). Il
   loop sopravvive solo grazie al warm start e al fatto che i solver restituiscono
   comunque il loro iterato meno infeasible.
2. **La regolarizzazione `sqrt(|T|^2 + 1e-6)` nella riga della massa della Jacobiana**
   (riga 266). Introduce un bias `O(eps)` proprio dove `|T| ~ eps`, cioe' sull'arco di
   coast, che e' anche dove il modello lineare fatica di piu'. La variante (d) elimina
   il problema alla radice.
3. **Il vincolo ridondante e non differenziabile `Tmin - |T| <= 0`** con `Tmin = 0`
   (riga 814). Non vincola nulla, ma introduce un punto di non differenziabilita' in
   `T = 0` che `fmincon` deve attraversare con differenze finite. E' la causa radice
   della convergenza lenta documentata nel README.
4. **La baseline del ratio test e' `ref.m_f` (valore NLP) sia al numeratore che al
   denominatore.** Rigorosamente, il numeratore dovrebbe usare la `m_f` **replayata** del
   riferimento. L'errore di trascrizione del riferimento si infiltra in `eta`.
5. **`ternary` valuta eagerly entrambi i rami** (righe 1091, 1253): la divisione
   `J_act/J_pred` viene calcolata anche quando `J_pred ~ 0`. Innocuo in MATLAB, ma la
   guardia non guardia niente. `solve_scvx` usa un `if/else` corretto: le tre
   implementazioni dello stesso loop non sono coerenti.
6. **`sol_best` non e' il migliore**, e' l'ultimo accettato (riga 535).
7. **Duplicazione massiccia**: `solve_scvx`, `solve_scvx_yalmip` e `solve_gfold_scvx`
   sono ~80 righe quasi identiche ciascuna. Idem `plot_compare3` / `plot_compare5`.
8. **Warm start fuori dai bound**: `ref_to_z` copia `ref.Tx(N)`, `ref.Ty(N)` dal
   riferimento trapezoidale, ma `box_bounds` li inchioda a `[0, 0]`. `fmincon` proietta
   silenziosamente.
9. **`z <= 0` hard-coded** in `solve_gfold_socp` (riga 1166): assume `m0_nd = 1`.
   Fragile se si cambiasse lo schema di non-dimensionalizzazione.
10. **Commenti disallineati**: il banner alla riga 57 dice "three transcriptions" (sono
    cinque); il commento alla riga 38 attribuisce `scvx_max` alla sola variante (b) (vale
    per tutte e tre); il messaggio alla riga 130 dice "state error" ma `node_err` esclude
    la massa.
11. **`tf` e' fissato.** Con `tf` libero il problema tornerebbe nonlineare anche in
    GFOLD; la strada canonica e' una ricerca 1-D su `tf` risolvendo un SOCP per ogni
    valore. Non implementata.
12. **Il ticket T006 (single-SOCP lossless) resta aperto.** L'obiettivo finale --
    eliminare del tutto il loop SCvx e risolvere **un solo SOCP** con certificato di
    ottimalita' globale -- non e' stato raggiunto: la variante (d) resta un loop SCvx
    perche' il bound superiore di spinta e' ancora linearizzato. Per chiuderlo servirebbe
    la sostituzione canonica (bound superiore approssimato con un'espansione al
    secondo ordine, o riformulazione che rende il vincolo convesso senza iterare).

---

## Possibili domande d'esame

**D: Il problema di powered descent e' non convesso. Da dove viene esattamente la non convessita', e come la affronti?**
R: Da due sorgenti distinte. (1) **La dinamica**: i termini `Tx/m`, `Ty/m` sono bilineari
e il termine di massa `m_dot = -Vc*|T|` mette una norma dentro un'**uguaglianza**, il che
rende non convesso l'insieme delle traiettorie ammissibili. (2) **Il vincolo di spinta
minima** `T_min <= |T| <= T_max`, che e' un guscio anulare (nel mio caso `T_min = 0`,
quindi questa seconda sorgente non c'e'). I vincoli di percorso -- glide-slope, bound
superiore di spinta, box -- sono invece **gia' convessi**. Quindi il bersaglio e' la
dinamica: o la si linearizza e si itera (SCvx, varianti b/c), o la si rende lineare per
**cambio di variabili** (GFOLD, variante d).

**D: Deriva la linearizzazione LTV e spiega perche' il termine affine `c(t)` e' indispensabile.**
R: Espando `f(x,u)` al primo ordine attorno a `(x_ref, u_ref)`:
`f ~= f(x_ref,u_ref) + A*(x - x_ref) + B*(u - u_ref)` con `A = df/dx`, `B = df/du`
valutate sul riferimento. Raccogliendo in variabili **assolute** ottengo
`x_dot ~= A*x + B*u + c` con `c = f(x_ref,u_ref) - A*x_ref - B*u_ref`. Il termine `c`
e' la costante che fa passare l'iperpiano tangente **per il punto di riferimento**: se
lo ometto, la retta tangente passa per l'origine e il modello e' completamente sbagliato
(ad esempio la gravita', che e' un termine costante di `f`, sparirebbe). Nel codice e' la
riga 296, e viene propagato nell'integrale `gamma` dell'ODE aumentata.

**D: A cosa serve la trust region, come e' implementata nel tuo codice, e come viene aggiornata?**
R: Serve contro l'**artificial unboundedness**: la linearizzazione e' accurata solo
localmente, e se lascio il sottoproblema convesso libero, questo "promette" un guadagno
di massa che la dinamica vera non consegna (tipicamente sfruttando la riga singolare
della massa, dove `df5/dT` e' quasi indefinita a spinta nulla). Nel mio codice e' una
trust region **hard** (non una penalita' soft): un box in norma infinito, per-variabile e
per-nodo, intersecato con i bound gia' esistenti (`apply_trust`, righe 571-600), con
raggi anisotropi `(0.17, 0.6, 0.1, 1.0)` non-dim e un fattore di scala adattivo
`rho`. L'aggiornamento e' il **ratio test**: calcolo
`eta = (guadagno reale in m_f) / (guadagno predetto dal modello LTV)`, dove il guadagno
reale viene da un replay `ode45` della dinamica **nonlineare** con i controlli candidati.
Se `eta < 0.25` rifiuto e dimezzo `rho`; se `0.25 <= eta < 0.7` accetto e tengo `rho`; se
`eta >= 0.7` accetto e raddoppio `rho` (cappato a 1).

**D: Cosa sono i virtual control, perche' servono, e li hai implementati?**
R: Servono contro l'**artificial infeasibility**: il sottoproblema linearizzato puo'
risultare infeasible anche se il problema nonlineare e' feasible, perche' la dinamica
lineare non riesce a "chiudere" fra stato iniziale assegnato e stato finale assegnato. La
cura canonica e' aggiungere alla dinamica discreta un termine libero `nu_k`
(`x_{k+1} = Abar*x_k + Bbar*u_k + cbar + nu_k`) e dei buffer `w >= 0` sui vincoli di
percorso, penalizzandoli in costo con una penalita' esatta L1
(`lambda_nu*||nu||_1 + lambda_w*||w||_1`), i cui pesi si alzano finche' `nu` e `w` si
annullano alla convergenza. **Nel mio codice non ci sono.** Il loop sopravvive per due
motivi: il warm start dalla soluzione trapezoidale mette il riferimento gia' vicino
all'ottimo, e quando il sottoproblema diventa infeasible il solver restituisce comunque
il suo iterato meno infeasible, che poi il ratio test contro `ode45` vaglia prima di
accettarlo. E' un comportamento fragile e una versione di produzione dovrebbe portarli.

**D: Deriva il cambio di variabili GFOLD e spiega perche' rende la dinamica esattamente LTI.**
R: Pongo `z = ln(m)`, `u = T/m` (accelerazione di spinta) e introduco lo slack
`sigma >= |u|`. Allora: `vx_dot = Tx/m = ux` e `vy_dot = Ty/m - 1 = uy - 1` diventano
lineari **per costruzione** (il prodotto bilineare `T/m` e' la nuova variabile). Per la
massa: `z_dot = m_dot/m = -Vc*|T|/m = -Vc*|u|`, e **sostituendo `|u|` con `sigma`**
ottengo `z_dot = -Vc*sigma`, lineare. Il sistema in `[x, y, vx, vy, z]` con controllo
`[ux, uy, sigma]` ha `A`, `B`, `c` **costanti**, quindi la ZOH dell'Appendice A si riduce
a un singolo esponenziale di matrice (van Loan), calcolato una volta per tutta la
griglia (`lti_zoh.m`). E siccome `A` e' nilpotente, `expm(A*dt) = I + A*dt` esattamente:
le equazioni discrete sono le formule del moto uniformemente accelerato, **senza nessun
errore di discretizzazione**. E' per questo che il replay dei controlli GFOLD ha errore
`7e-12`, cioe' il floor di `ode45`, invece dei `1e-3`-`1e-8` delle altre trascrizioni.

**D: Perche' il rilassamento `sigma >= |u|` e' "lossless"?**
R: Perche' all'ottimo il vincolo e' **attivo**, `sigma = |u|`, quindi la soluzione del
problema rilassato (convesso) e' anche soluzione del problema originale. Nel mio setup lo
si vede in modo diretto: la riga discreta della massa e' `z_{k+1} = z_k - Vc*dt*sigma_k`
(esatta, perche' la quinta riga di `A` e' nulla), quindi
`z_N = z_0 - Vc*dt*sum(sigma_k)`. Il costo e' `max z_N`, cioe' `min sum(sigma_k)`: ogni
`sigma_k` entra nel costo con coefficiente **positivo**, e negli altri vincoli compare
solo come `sigma_k >= |u_k|` (limite inferiore) e `sigma_k <= Tmax*e^{-z}` (limite
superiore, che abbassare non viola). L'ottimizzatore lo schiaccia quindi sul suo limite
inferiore `|u_k|`. Uno `sigma` "in eccesso" brucerebbe carburante senza produrre
accelerazione: e' strettamente sub-ottimo. **Attenzione:** questo argomento e' il caso
facile (`T_min = 0`). Il risultato forte di Acikmese & Ploen (2007) e' che il
rilassamento resta lossless anche con `T_min > 0`, e li' la dimostrazione passa per il
Principio del Massimo (un cono inattivo forzerebbe `lambda_m(tf) = 0`, in contraddizione
con la condizione di trasversalita').

**D: Se la dinamica GFOLD e' esatta e il cono e' lossless, perche' il tuo codice fa ancora un loop SCvx? Cosa servirebbe per arrivare a un singolo SOCP?**
R: Resta **una** non convessita': il bound superiore di spinta. `|T| <= Tmax` con
`T = e^z * u` diventa `sigma <= Tmax * e^{-z}`, e `e^{-z}` e' **convessa**: vincolare
`sigma` sotto una funzione convessa da' un insieme non convesso. Lo linearizzo attorno a
`z_ref` (riga 1159), il che e' **conservativo** (la tangente di una funzione convessa sta
sempre sotto la funzione, quindi ogni iterato accettato e' feasible per il problema
vero), e itero finche' `z_ref` insegue la soluzione. Ma e' **un solo vincolo scalare per
nodo**, non l'intera dinamica: infatti il loop converge in 3 iterazioni con `eta ~= 1` a
ogni passo, contro le 15 iterazioni non convergenti delle varianti (b)/(c). Per arrivare
a un **singolo SOCP** (ottimo globale certificato, nessun warm start, nessuna iterazione)
servirebbe la formulazione canonica GFOLD completa -- e' il ticket T006, ancora aperto.

**D: Le varianti (b) e (c) risolvono lo stesso identico sottoproblema. Perche' allora danno risultati diversi?**
R: Perche' lo risolvono in modo diverso. Il sottoproblema e' un SOCP: obiettivo lineare
(`max m_N`), uguaglianze lineari (dinamica LTV discreta + condizioni al contorno),
disuguaglianze lineari (glide-slope, box, trust region) e **un cono di Lorentz per nodo**
(`|T_k| <= Tmax`). La variante (c) lo passa a ECOS, che tratta il cono come **primitiva**
e lo risolve con un metodo interior-point in tempo polinomiale. La variante (b) lo passa a
`fmincon`/SQP, che **butta via** la struttura convessa: approssima il cono con successive
linearizzazioni quadratiche e calcola le Jacobiane alle **differenze finite** (350
variabili -> ~350 valutazioni del `nonlcon` per gradiente). Il risultato: (c) e' ~3x piu'
veloce end-to-end (29 s vs 82 s) e la fedelta' per-passo e' un ordine di grandezza
migliore (`4e-5` vs `4e-4`). Non e' un difetto di `fmincon` -- e' che lo si sta usando
fuori dalla sua classe.

**D: Le tue due varianti LTV+SCvx non convergono. E' un problema?**
R: E' un limite dichiarato. Entrambe arrivano al cap di 15 iterazioni con
`||delta x|| ~ 1e-2`, sopra la tolleranza di `1e-3`. Le cause sono tre: (i) manca il
virtual control, quindi quando il sottoproblema diventa infeasible il loop procede su un
iterato non del tutto pulito; (ii) la regolarizzazione `sqrt(|T|^2 + 1e-6)` della riga di
massa introduce un bias `O(eps)` proprio sull'arco di coast; (iii) il ratio test usa come
baseline il valore NLP di `m_f` del riferimento invece del suo replay, quindi `eta` porta
dentro l'errore di trascrizione del riferimento. Detto questo, **ogni passo accettato ha
superato il ratio test contro un'integrazione `ode45` della dinamica nonlineare vera**,
quindi l'iterato restituito e' fisicamente realizzabile: replayando i controlli si atterra
a 0.18 m (fmincon) e 0.019 m (ECOS) dalla piazzola. La `m_f` finale differisce di 2-4 kg
su 1400 dalle altre trascrizioni, cioe' lo 0.2%: e' un ottimo locale vicino, non un
risultato sbagliato. La variante GFOLD, che rimuove le cause (ii) e in gran parte (i),
converge in 3 iterazioni e atterra a `2e-8` m.
