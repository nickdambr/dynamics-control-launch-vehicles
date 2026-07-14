# HM2_powered_descent/lti_zoh.m

## Ruolo del file nel progetto

Questo file discretizza **esattamente** il sistema LTI che si ottiene dopo il
cambio di variabili GFOLD (Acikmese & Blackmore): z = ln m, u = T/m, con lo slack
sigma >= ||u||. E' il cuore numerico della **variante (d) di Task 2**, ed e' il
file dove si vede materializzato tutto il guadagno del cambio di variabili: la
dinamica di powered descent, che in coordinate (m, T) e' non lineare e non
convessa, diventa qui un sistema **lineare, tempo-invariante e con offset
costante**, discretizzabile con un unico esponenziale di matrice.

Il contrasto con le altre varianti e' netto:

| Variante | Come si ottengono le matrici discrete | Costo |
|---|---|---|
| (a) ZOH+RK4 | non si linearizza affatto: `rk4_zoh` propaga la dinamica non lineare | 392 valutazioni RHS per vincolo |
| (b)/(c) LTV+SCvx | `compute_ltv_zoh` integra con `ode45` un sistema aumentato a 70 stati (Phi, Psi, Beta, Gamma) su ciascuno dei 49 intervalli, **a ogni iterazione SCvx** | 49 chiamate `ode45` per iterazione |
| (d) GFOLD | `lti_zoh`: un `expm` di una matrice 9x9, **una volta sola** | 1 `expm`, punto |

`main_task2.m` riga 1218: `[Abar, Bbar, cbar] = lti_zoh(tf/(N-1), d.Vc);` con il
commento *"exact LTI ZOH, once"* -- fuori dal ciclo SCvx, perche' le matrici non
dipendono ne' dal nodo (tempo-invarianza) ne' dalla traiettoria di riferimento
(non c'e' nulla da linearizzare nella dinamica). Le tre matrici finiscono poi nel
vincolo di dinamica del SOCP (`solve_gfold_socp`, riga 1156):

    XI(:,k+1) == Abar*XI(:,k) + Bbar*W(:,k) + cbar

che e' un'**uguaglianza lineare** nelle variabili di decisione -- ed e' esattamente
questo che permette al sotto-problema di essere un SOCP genuino, risolvibile da
ECOS in frazioni di secondo con certificato di ottimalita' globale.

Testato da `tests/gfoldLogMassTest.m` righe 55-83 (forma chiusa analitica +
confronto con `ode45`).

---

## `lti_zoh` (righe 1-28)

```matlab
function [Abar, Bbar, cbar] = lti_zoh(dt, Vc)
...
A = zeros(5);    A(1,3) = 1;  A(2,4) = 1;
B = zeros(5,3);  B(3,1) = 1;  B(4,2) = 1;  B(5,3) = -Vc;
c = [0; 0; 0; -1; 0];
E = expm([A, B, c; zeros(4,9)] * dt);
Abar = E(1:5, 1:5);   Bbar = E(1:5, 6:8);   cbar = E(1:5, 9);
```

- **Riga 1**: firma. Restituisce **tre** oggetti: la matrice di transizione
  `Abar` (5x5), la matrice di ingresso `Bbar` (5x3) e l'**offset affine** `cbar`
  (5x1). Il terzo non e' un dettaglio: la gravita' e' un termine costante, non un
  controllo, quindi il sistema non e' xi_dot = A*xi + B*w ma
  xi_dot = A*xi + B*w + c, e la sua discretizzazione richiede un terzo pezzo.

- **Righe 2-18**: la docstring, che dichiara esplicitamente lo stato
  `xi = [x; y; vx; vy; z]` con z = ln m, il controllo `w = [ux; uy; sigma]` e la
  dinamica continua. Righe 11-13 sottolineano che le matrici sono **CONSTANT
  across the grid** perche' il sistema e' tempo-invariante. Righe 15-18 anticipano
  il trucco della matrice a blocchi.

- **Righe 19-22: il blocco `arguments`.** *Attenzione: contrariamente a
  `ode_descent.m`, `ode_descent_uacc.m` e `rk4_zoh.m`, questo file **ha** la
  validazione degli input* (`dt` positivo e finito, `Vc` finito). Non e' una
  dimenticanza ne' un'incoerenza: `lti_zoh` **non e' una funzione di hot loop**.
  Viene chiamata una volta sola per risoluzione (riga 1218 di `main_task2.m`,
  fuori dal ciclo SCvx), quindi ricade nella categoria "boundary helper" della
  convenzione di `CLAUDE.md`, dove la validazione e' dovuta. Questa e' proprio la
  dimostrazione che la regola "niente `arguments` nel hot loop" e' applicata con
  criterio e non a tappeto.
  Nota minore: `Vc` e' validato solo come finito, non come positivo. Vc = 0
  (Isp infinito, nessun consumo) e' un caso limite fisicamente sensato e resta
  ammesso; Vc < 0 sarebbe non fisico ma non viene intercettato.

- **Riga 23**: `A(1,3) = 1; A(2,4) = 1` -- le due sole entrate non nulle: x_dot =
  vx e y_dot = vy. E' un **doppio integratore** in 2D. Il resto di A e' zero:
  nessun accoppiamento con la massa, nessun drag, nessun gradiente di gravita'.

- **Riga 24**: `B(3,1) = 1; B(4,2) = 1; B(5,3) = -Vc` -- l'accelerazione comandata
  entra direttamente nelle righe di velocita' (coefficiente 1, non 1/m: e' il
  guadagno del cambio di variabili u = T/m), e lo slack sigma entra nella riga di
  log-massa con coefficiente -Vc. Cioe' **z_dot = -Vc*sigma**: la riga della
  massa non dipende piu' dallo stato.

- **Riga 25**: `c = [0; 0; 0; -1; 0]` -- la gravita' adimensionale, che agisce
  solo su vy_dot. E' un **termine affine**, non un controllo: e' noto e costante,
  e va discretizzato come tale.

- **Riga 26**: l'esponenziale della matrice aumentata (vedi sotto).

- **Riga 27**: estrazione dei tre blocchi dal risultato 9x9.

### Derivazione della discretizzazione ZOH esatta

Il sistema continuo e'

    xi_dot(t) = A*xi(t) + B*w + c        con w costante su [t_k, t_k + dt]  (ZOH)

La soluzione esatta di un sistema lineare con ingresso costante e' la formula di
variazione delle costanti:

    xi(t_k + dt) = expm(A*dt)*xi(t_k) + integrale_0^dt expm(A*s) ds * (B*w + c)

da cui, per confronto con xi_{k+1} = Abar*xi_k + Bbar*w_k + cbar:

    Abar = expm(A*dt)
    Bbar = [ integrale_0^dt expm(A*s) ds ] * B
    cbar = [ integrale_0^dt expm(A*s) ds ] * c

Il fatto che l'ingresso sia **costante** e' cio' che permette di portarlo fuori
dall'integrale. Se w(t) variasse (per esempio FOH), l'integrale sarebbe
integrale expm(A*(dt - s))*B*w(s) ds e servirebbero due matrici di ingresso
distinte.

### Il trucco della matrice a blocchi (van Loan)

Calcolare l'integrale di un esponenziale di matrice a mano e' scomodo (e diventa
singolare se A non e' invertibile -- qui A **non lo e'**, e' addirittura
nilpotente, quindi la scorciatoia `A^-1*(expm(A*dt) - I)*B` **non e' utilizzabile**).
Il trucco di van Loan lo evita: si costruisce la matrice aumentata (9x9)

    M = [ A   B   c ]  * dt
        [ 0   0   0 ]        <- 4 righe di zeri (3 per w, 1 per la costante)

e si dimostra che

    expm(M) = [ Abar   Bbar   cbar ]
              [  0      I3     0   ]
              [  0      0      1   ]

**Dimostrazione** (serve saperla): sia M = [[F, G], [0, 0]] con F = A*dt e
G = [B c]*dt. Per induzione, M^n = [[F^n, F^(n-1)*G], [0, 0]] per ogni n >= 1
(il blocco in basso e' nullo, quindi le potenze non lo "risvegliano" mai).
Sommando la serie:

    expm(M) = I + somma_{n>=1} M^n/n!
            = [ somma F^n/n!        (somma_{n>=1} F^(n-1)/n!) * G ]
              [ 0                    I                            ]

Il blocco in alto a sinistra e' expm(F) = expm(A*dt) = **Abar**. Il blocco in alto
a destra: somma_{n>=1} F^(n-1)/n! = integrale_0^1 expm(F*s) ds (basta integrare la
serie termine a termine). Quindi il blocco vale
integrale_0^1 expm(A*dt*s) * [B c] * dt ds, e con la sostituzione tau = dt*s
(dtau = dt*ds) diventa

    integrale_0^dt expm(A*tau) * [B c] dtau  =  [ Bbar   cbar ]

che e' esattamente cio' che cercavamo. Il codice lo scrive nella docstring alla
riga 17 e lo implementa alla riga 26. Un solo `expm` di una 9x9 restituisce
**tutte e tre** le matrici in un colpo, senza integrare nulla e senza invertire A.

### Verifica in forma chiusa

Qui la struttura di A e' cosi' semplice che la forma chiusa si scrive a mano -- ed
e' quello che fa il test `testZohClosedForm` (`gfoldLogMassTest.m` riga 55).

**A e' nilpotente di indice 2**: (A^2)(1,j) = riga 3 di A = 0, (A^2)(2,j) = riga 4
di A = 0, quindi **A^2 = 0**. La serie esponenziale si tronca dopo due termini:

    expm(A*s) = I + A*s
    integrale_0^dt (I + A*s) ds = dt*I + (dt^2/2)*A

Da cui:

    Abar = I + A*dt         ->  Abar(1,3) = dt,  Abar(2,4) = dt, resto identita'
    Bbar = dt*B + (dt^2/2)*A*B
    cbar = dt*c + (dt^2/2)*A*c

Svolgendo (A*B ha riga 1 = riga 3 di B e riga 2 = riga 4 di B):

    Bbar = [ dt^2/2   0        0      ]        cbar = [ 0       ]
           [ 0        dt^2/2   0      ]               [ -dt^2/2 ]
           [ dt       0        0      ]               [ 0       ]
           [ 0        dt       0      ]               [ -dt     ]
           [ 0        0       -Vc*dt  ]               [ 0       ]

Che sono **letteralmente** `B_exp` e `c_exp` del test (righe 60-65). La lettura
fisica e' immediata:

- righe 1-2 di Bbar: lo spostamento dovuto a un'accelerazione costante,
  (1/2)*a*dt^2;
- righe 3-4: la variazione di velocita', a*dt;
- riga 2 e 4 di cbar: la **caduta libera** su un intervallo, -(1/2)*g*dt^2 di
  quota e -g*dt di velocita' verticale (con g = 1 adimensionale);
- riga 5 di Bbar: **z_{k+1} = z_k - Vc*sigma_k*dt**, la depletion di log-massa,
  perfettamente lineare.

Quest'ultima riga chiude anche il cerchio sulla funzione obiettivo: poiche'
z_N = z_0 - Vc*dt*somma(sigma_k), massimizzare z_N (il costo del SOCP e'
`-XI(5,N)`, riga 1183) equivale a **minimizzare la somma degli sigma**, cioe' a
minimizzare il consumo. L'obiettivo e' lineare nelle variabili di decisione.

> **Possibile domanda d'esame** -- Perche' non calcolare Bbar con la formula
> A^-1 * (expm(A*dt) - I) * B, che si trova su tutti i libri?
> *Risposta:* Perche' qui **A e' singolare** (e' nilpotente: A^2 = 0, quindi tutti
> gli autovalori sono nulli e det(A) = 0), e quella formula richiede A
> invertibile. E' un caso tutt'altro che raro: qualunque doppio integratore ha A
> nilpotente. Il trucco della matrice a blocchi funziona **sempre**, invertibile o
> no, ed e' anche numericamente piu' stabile (non moltiplica per l'inversa di una
> matrice mal condizionata).

---

## "Esatta": in che senso, e quando conviene davvero rispetto a RK4

**In che senso e' esatta.** La mappa xi_{k+1} = Abar*xi_k + Bbar*w_k + cbar
riproduce la soluzione della ODE continua sotto ZOH **senza alcun errore di
troncamento**, per qualunque dt: non c'e' un O(h^p) da controllare, non c'e' un
`n_sub` da scegliere. L'unico errore e' quello di macchina di `expm` (~1e-16). Il
test `testZohMatchesOde45` (riga 71) lo conferma confrontando un passo discreto
con un'integrazione `ode45` a RelTol 1e-12 (tolleranza 1e-9), e il replay non
lineare completo raggiunge un errore di nodo di 7.3e-12 nondim (README).

**Onesta' -- in questo caso specifico anche RK4 sarebbe esatto.** Vale la pena
dirlo perche' e' una domanda che un professore puo' fare. Poiche' A^2 = 0, la
soluzione esatta su un intervallo e' un **polinomio di grado 2 in dt**, e RK4 e'
esatto sui polinomi fino al grado 4. Facendo il conto (con b = B*w + c e A^2 = 0):

    k1 = A*xi + b
    k2 = k3 = k1 + (h/2)*A*b
    k4 = k1 + h*A*b
    xi + (h/6)*(k1 + 2*k2 + 2*k3 + k4) = xi + h*(A*xi + b) + (h^2/2)*A*b

che coincide con (I + A*h)*xi + (h*I + (h^2/2)*A)*b, cioe' con la formula esatta.
Quindi **su questo particolare sistema** RK4 anche con n_sub = 1 darebbe lo stesso
risultato. Il vantaggio dell'`expm` **non e' l'accuratezza in questo caso**: e'
strutturale.

**I veri motivi per cui `lti_zoh` e' preferibile:**

1. **Costo.** Le matrici sono costanti: un `expm` 9x9, una volta sola per l'intera
   risoluzione. La variante LTV (b/c) deve integrare con `ode45` un sistema
   aumentato a 70 stati su 49 intervalli **a ogni iterazione SCvx**; la variante
   (a) ricalcola 392 valutazioni del RHS a ogni valutazione di vincolo. Il GFOLD
   converge in 3 iterazioni SCvx in ~5 s di wall time (README).

2. **Convessita'.** Con Abar, Bbar, cbar **costanti e note**, il vincolo di
   dinamica e' un'uguaglianza **lineare** nelle variabili di decisione: il
   sotto-problema e' un SOCP genuino, che ECOS risolve a ottimalita' globale. Con
   RK4 sulla dinamica non lineare in (m, T) il difetto e' un'uguaglianza non
   convessa e serve `fmincon`, con tutte le sue patologie (minimi locali,
   gradienti FD, stallo sul kink di ||T||).

3. **Niente riga di massa singolare** (dichiarato alla riga 18 della docstring).
   Nelle coordinate originali la linearizzazione ha A(3,5) = -Tx/m^2 e
   A(4,5) = -Ty/m^2 (`jacobians`, `main_task2.m` righe 270-271): termini che
   esplodono quando m -> 0 e che rendono la linearizzazione fragile a fine
   traiettoria, dove il razzo e' leggero. In log-massa questi termini **non
   esistono**: le righe di velocita' non dipendono affatto dallo stato.

4. **Generalita'.** Se il modello venisse arricchito (drag lineare, frame
   rotante con Coriolis, gradiente di gravita'), A **non sarebbe piu' nilpotente**
   e la soluzione esatta diventerebbe genuinamente trascendente: li' RK4 avrebbe
   un errore O(h^4) reale, mentre `expm` resterebbe esatto. La formulazione a
   blocchi si adatta senza modifiche.

5. **Nessun parametro da sintonizzare.** Nessun `n_sub`, nessuna analisi di
   convergenza di griglia sull'integratore.

---

## Cosa resta non convesso (e perche' c'e' ancora una SCvx)

E' il limite piu' importante da esplicitare. La **dinamica** e' esatta e lineare,
ma **un vincolo non lo e'**: il bound superiore di spinta ||T|| <= Tmax, riscritto
in coordinate GFOLD, diventa

    sigma <= Tmax * exp(-z)

Il membro destro e' **convesso in z**, quindi "variabile <= funzione convessa" e'
un vincolo **non convesso**. Il codice lo linearizza con Taylor al prim'ordine
attorno al riferimento z_ref (`main_task2.m` righe 1158-1159):

    sigma_k <= Tmax * exp(-z_ref_k) * (1 - (z_k - z_ref_k))

e itera su z_ref nel ciclo SCvx, con trust region adattiva e ratio test. **Quindi
la variante (d) non e' un singolo SOCP**: e' una successione di SOCP. Il vero
single-shot (lossless convexification completa) e' ancora nel backlog come ticket
T006, come dichiara il README.

Nota: l'altro vincolo tipicamente non convesso -- il **lower bound** ||T|| >= Tmin
(vincolo "fuori da una palla") -- in HM2 non si pone, perche' `data.Tmin = 0`
(`main_task2.m` riga 32). Con Tmin > 0 servirebbe la maggiorazione quadratica di
Acikmese-Blackmore, ed e' proprio li' che il teorema di *lossless convexification*
diventa non banale.

---

## Possibili domande d'esame

**D: Deriva Abar e Bbar per un sistema LTI con ingresso ZOH.**
R: Dalla variazione delle costanti,
xi(t_k + dt) = expm(A*dt)*xi_k + integrale_0^dt expm(A*(dt - tau))*B*w(tau) dtau.
Con w costante (ZOH) l'ingresso esce dall'integrale e, con il cambio di variabile
s = dt - tau, si ottiene
Bbar = [integrale_0^dt expm(A*s) ds] * B, e Abar = expm(A*dt). L'offset costante c
(la gravita') si tratta identicamente: cbar = [integrale_0^dt expm(A*s) ds] * c.
Sono **esatti**, non approssimati: nessun errore di troncamento, per qualunque dt.

**D: Spiega il trucco della matrice a blocchi. Perche' funziona?**
R: Si costruisce M = [[A, B, c], [0, 0, 0]] * dt (9x9, con 4 righe di zeri).
Poiche' M e' triangolare a blocchi con blocco diagonale inferiore nullo, si ha
M^n = [[F^n, F^(n-1)*G], [0, 0]] con F = A*dt, G = [B c]*dt. Sommando la serie
esponenziale, il blocco in alto a sinistra da' expm(A*dt) = Abar e quello in alto
a destra da' (somma_{n>=1} F^(n-1)/n!)*G = integrale_0^1 expm(F*s) ds * G, che con
la sostituzione tau = dt*s e' esattamente [Bbar cbar]. Un solo `expm` restituisce
tutte e tre le matrici. Il vantaggio pratico: funziona **anche se A e' singolare**
(qui A e' nilpotente, quindi la formula A^-1*(expm(A*dt) - I)*B **non** sarebbe
applicabile).

**D: Perche' la dinamica GFOLD e' esattamente LTI, mentre quella originale non lo
e'?**
R: Per il doppio cambio di variabili di Acikmese & Blackmore. (i) Prendendo come
controllo l'accelerazione u = T/m invece della spinta, le righe di velocita'
perdono il fattore 1/m e diventano un doppio integratore: vx_dot = ux,
vy_dot = uy - 1. (ii) Prendendo come stato z = ln m invece di m, la riga della
massa m_dot = -Vc*m*||u|| diventa z_dot = -Vc*||u||, indipendente dallo stato.
(iii) Introducendo lo slack sigma >= ||u|| (vincolo conico convesso) si ha
z_dot = -Vc*sigma, lineare. Tutte le non linearita' della dinamica sono state
spostate in un vincolo **convesso**. Nel codice, il punto (iii) e' `B(5,3) = -Vc`
(riga 24) e il cono e' `norm(W(1:2,k)) <= W(3,k)` (riga 1157 di `main_task2.m`).

**D: Il costo del SOCP e' `-XI(5,N)`. Perche' e' equivalente a minimizzare il
consumo, e come si sa che sigma = ||u|| all'ottimo?**
R: XI(5,N) = z_N = ln(m_f), quindi massimizzare z_N e' massimizzare m_f, cioe'
minimizzare il propellente (m0 e' fissata). Dalla riga 5 di Bbar,
z_N = z_0 - Vc*dt*somma(sigma_k): massimizzare z_N significa minimizzare la somma
degli sigma, e l'ottimizzatore li spinge quindi contro il loro **unico lower
bound**, che e' il cono sigma_k >= ||u_k||. Il cono si chiude da solo, senza
bisogno di imporlo: e' il nucleo dell'argomento di *lossless convexification*.
Il test `testMassRowConsistency` (`gfoldLogMassTest.m` riga 85) verifica
numericamente che, con sigma = ||u||, l'aggiornamento LTI di z coincide con la
depletion di massa non lineare integrata da `ode45`.

**D: Se la discretizzazione e' esatta e la dinamica lineare, perche' serve ancora
un ciclo SCvx?**
R: Perche' e' la **dinamica** a essere esatta, non tutti i vincoli. Il bound
superiore ||T|| <= Tmax diventa sigma <= Tmax*exp(-z), che e' non convesso
(maggiorante convesso). Il codice lo linearizza attorno a z_ref e itera. Con la
lossless convexification completa (ticket T006) si potrebbe fissare il punto di
espansione su un profilo di massa a priori e risolvere **un solo** SOCP; qui si e'
scelta la strada SCvx, che comunque converge in 3 iterazioni perche' la dinamica
esatta non introduce errore di linearizzazione.

**D: Questo file ha il blocco `arguments` mentre gli altri tre no. Incoerenza?**
R: No, e' esattamente la regola applicata correttamente. `ode_descent`,
`ode_descent_uacc` e `rk4_zoh` sono chiamati ~1e5 volte per iterazione dentro
`fmincon`/`ode45`, e per loro la validazione costerebbe piu' del calcolo (una
valutazione di `ode_descent` sta sui 60 ns). `lti_zoh` invece e' chiamata **una
volta per risoluzione** (riga 1218, fuori dal ciclo SCvx), quindi e' un "boundary
helper" e la validazione e' dovuta e gratuita. La convenzione di `CLAUDE.md` dice
proprio questo: niente `arguments` nel hot loop, validazione al confine.

**D: Cosa succederebbe se dimenticassi `cbar`?**
R: Il razzo cadrebbe in assenza di gravita'. `cbar = [0; -dt^2/2; 0; -dt; 0]` e'
la discretizzazione esatta del termine affine c = [0;0;0;-1;0], cioe' la caduta
libera: -(1/2)*g*dt^2 di quota e -g*dt di velocita' verticale per intervallo (con
g = 1 adimensionale). Non e' assorbibile in Bbar perche' la gravita' **non e' un
controllo**: e' un ingresso noto e costante. E' la ragione per cui il sistema si
scrive xi_dot = A*xi + B*w + c e non semplicemente xi_dot = A*xi + B*w, e per cui
la matrice aumentata ha una colonna in piu' (la nona).
