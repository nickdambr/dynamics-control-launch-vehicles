# HM1/ode_burn.m

## Ruolo del file nel progetto

`ode_burn.m` e' il **cuore numerico di tutta HM1**. E' il right-hand side (RHS)
dell'arco propulso: la funzione che `ode45` chiama a ogni step per propagare in
avanti lo stato del lanciatore **con la legge di controllo ottima gia' sostituita
dentro**. Non e' un integratore "della dinamica" e basta: e' la dinamica *chiusa
sulle condizioni necessarie di Pontryagin*, cioe' il sistema Hamiltoniano ridotto
che nasce dal PMP dopo aver eliminato il controllo `phi` tramite la condizione di
massimizzazione.

Il file esiste perche' il metodo indiretto trasforma il problema di controllo
ottimo in un **Boundary Value Problem (BVP)**: dati i costati iniziali incogniti,
si integra in avanti e si guarda quanto si sbaglia sulle condizioni al contorno
finali. `ode_burn` e' esattamente il pezzo "si integra in avanti". E' condiviso da
**tutti e quattro i task**: lo chiamano `main_task1.m` (dentro `shooting1` e per
la ri-integrazione delle soluzioni convergenti), `main_task2.m` (arco propulso
dopo la salita verticale), `main_task3.m` (arco propulso prima del coast),
`main_task4.m` e `validate_staging_corner.m` (i due archi separati dallo staging),
piu' le suite `HM1/tests/odeBurnTest.m` e `odeBurnPerformanceTest.m`.

Implementa le equazioni del moto planare (Introduction.tex, eq. xdot-mdot) piu' le
equazioni di Eulero-Lagrange dei costati (eq. lx-lm), con la legge bilinear-tangent
(eq. bilinear) inserita al posto del controllo. Il modello e' 2D, Terra piatta,
gravita' costante, **niente atmosfera, niente drag, spinta sempre accesa e di
modulo costante** `T = c*Q`.

Dipende solo dalla struct `p` che gli passa il chiamante: `p.T`, `p.Q`, `p.c` (in
realta' `p.c` **non viene mai letto** da questa funzione -- vedi sotto) e i tre
parametri dei costati `p.lam_vx0`, `p.lam_vy0`, `p.lam_y`, che sono esattamente
tre delle quattro incognite dello shooting.

---

## `ode_burn` -- intestazione e contratto (righe 1-14)

```matlab
function dz = ode_burn(t, z, p)
% Powered-flight RHS, linear-tangent steering.
%   z - state [x; y; vx; vy; m; lam_m]
%   p - struct: T, Q, c, lam_vx0, lam_vy0, lam_y
```

- Riga 1: firma. Tre argomenti nell'ordine imposto da `ode45`: tempo scalare `t`,
  vettore di stato `z` (6x1), struct di parametri `p` passata via closure
  (`@(t,z) ode_burn(t,z,pp)`). Restituisce `dz`, la derivata 6x1.
- Riga 5: il vettore integrato e' `[x; y; vx; vy; m; lam_m]`. Attenzione: **non e'**
  lo stato fisico e basta (che sarebbe 5 componenti, `x, y, vx, vy, m`), e **non e'**
  nemmeno stato + tutti i costati (che sarebbero 5+5 = 10). E' un ibrido: 5 stati
  fisici piu' **un solo costato**, `lam_m`. Il motivo e' spiegato sotto (righe 18-21).
- Righe 10-11: il commento dichiara la struttura analitica dei costati:
  `lam_x = 0`, `lam_y = p.lam_y` (costante), `lam_vx = p.lam_vx0` (costante),
  `lam_vy = p.lam_vy0 - p.lam_y*t` (rampa lineare). Sono le soluzioni in forma
  chiusa delle equazioni di Eulero-Lagrange -- quindi **non serve integrarle
  numericamente**.
- Riga 11: `phi = atan2(lam_vy, lam_vx)` -- la legge di controllo, dichiarata gia'
  nell'header.
- Righe 13-14: **niente blocco `arguments` per scelta di design**. Questa non e'
  pigrizia: la funzione sta nel loop piu' interno di `ode45`, che a sua volta sta
  nel loop di `fsolve`, che a sua volta sta nel loop di continuazione su `Q`
  (80 valori) e su `yf` (3 valori). Il numero di chiamate e' dell'ordine di 1e6+.
  Un blocco `arguments` con `mustBeNumeric`/`mustBeScalar` costa qualche
  microsecondo a chiamata e diventerebbe una frazione dominante del tempo totale.
  La validazione degli input vive **al boundary** (una volta per run), non nel hot
  loop. E' una convenzione dichiarata nel `CLAUDE.md` della repo e va difesa, non
  "aggiustata".

> **Possibile domanda d'esame** -- Perche' il vettore integrato ha 6 componenti e
> non 10 (5 stati + 5 costati)?
> *Risposta:* Perche' quattro delle cinque equazioni dei costati si risolvono in
> forma chiusa. Da `lambda_dot = -dH/dx` si ottiene `lam_x_dot = 0` (e la
> trasversalita' su `x` libero da' `lam_x = 0` identicamente), `lam_y_dot = 0`
> (costante), `lam_vx_dot = -lam_x = 0` (costante), `lam_vy_dot = -lam_y` (rampa
> lineare). Solo `lam_m_dot` dipende dallo stato (`m`) e va integrato. Quindi si
> integrano 5 stati fisici + `lam_m` = 6 equazioni, e i tre parametri
> `lam_vx0, lam_vy0, lam_y` bastano a ricostruire analiticamente tutti gli altri
> costati a qualunque `t`. Si risparmiano 4 ODE e, cosa piu' importante, si elimina
> l'errore di integrazione su di esse.

---

## Unpack dello stato (riga 16)

```matlab
vx = z(3); vy = z(4); m = z(5);
```

- Riga 16: si estraggono solo le tre componenti che servono al RHS. `x = z(1)` e
  `y = z(2)` non compaiono a destra di nessuna equazione: la dinamica **non dipende
  dalla posizione** (Terra piatta, `g` costante, niente atmosfera). E' proprio questo
  che rende `lam_x` e `lam_y` costanti (`dH/dx = dH/dy = 0`). Se ci fosse il drag,
  `dH/dy ~= 0` (densita' esponenziale in `y`), `lam_y` non sarebbe piu' costante,
  `lam_vy` non sarebbe piu' lineare in `t`, e **tutta la struttura linear-tangent
  cadrebbe**: bisognerebbe integrare anche i costati numericamente.
- Anche `z(6) = lam_m` non viene letto: la sua derivata (riga 33) dipende da `m` e
  dai costati di velocita', non da se stesso. L'equazione di `lam_m` e' quindi una
  **pura quadratura accodata**, accoppiata in una sola direzione.

---

## Costati: ricostruzione analitica (righe 18-21)

```matlab
lam_vx = p.lam_vx0;
lam_vy = p.lam_vy0 - p.lam_y * t;
lam_v_norm = sqrt(lam_vx^2 + lam_vy^2);
```

**Da dove viene questa forma.** L'Hamiltoniana del problema (Mayer, massimizzare
`m(tf)`) e', con `g = 1` in nondimensionale:

    H = lam_x*vx + lam_y*vy + lam_vx*(T/m)*cos(phi)
        + lam_vy*((T/m)*sin(phi) - 1) - lam_m*Q

Le equazioni di Eulero-Lagrange `lambda_dot = -dH/dx` danno, componente per
componente:

    lam_x_dot  = -dH/dx  = 0            -> lam_x  = cost = 0
    lam_y_dot  = -dH/dy  = 0            -> lam_y  = cost
    lam_vx_dot = -dH/dvx = -lam_x = 0   -> lam_vx = cost = lam_vx0
    lam_vy_dot = -dH/dvy = -lam_y       -> lam_vy = lam_vy0 - lam_y*t
    lam_m_dot  = -dH/dm  = +(T/m^2)*|lam_v|

- Riga 19: `lam_vx` costante. Deriva da `lam_vx_dot = -lam_x` e da `lam_x = 0`.
  Il fatto che `lam_x = 0` non e' un'ipotesi: e' la **condizione di trasversalita'**
  su `x(tf)` libero e non presente nel costo (`dphi/dx = 0`), combinata con
  `lam_x_dot = 0` che la propaga all'indietro su tutto l'arco.
- Riga 20: `lam_vy` rampa lineare. Il segno meno viene da `lam_vy_dot = -lam_y`.
  **Interpretazione fisica del segno:** con `lam_y > 0` (guess iniziale del Task 1:
  `lam_y = 14`), `lam_vy` *decresce* nel tempo, quindi `phi` scende da quasi
  verticale a orizzontale (e persino sotto l'orizzonte). E' il **pitch-over**: la
  legge lo produce da sola, non e' programmato a mano. Se si mettesse `lam_y < 0`
  il razzo alzerebbe il muso durante la salita -- non ottimo.
- Riga 21: `lam_v_norm = |lam_v|` e' il modulo del **primer vector**. Compare in due
  posti: nella derivata di `lam_m` (riga 33) e -- implicitamente -- nella funzione di
  switching `S = |lam_v|/m - lam_m/c` usata dal Task 3 per il cutoff.
- **Nota di onesta':** `t` qui e' il tempo **globale** dell'arco, cioe' lo stesso
  tempo rispetto a cui `lam_vy0` e' definito. In `main_task4.m` il secondo arco
  viene integrato su `[ts tf]` (non su `[0 tf-ts]`) proprio perche' la rampa deve
  restare continua. Se si passasse un `t` locale ripartito da 0 su un arco che non
  inizia a `t = 0`, la legge di controllo verrebbe silenziosamente sbagliata.
  Il codice non ha alcun assert su questo: e' un contratto implicito.

> **Possibile domanda d'esame** -- Cosa significa fisicamente che `lam_vy` e' lineare
> in `t` mentre `lam_vx` e' costante?
> *Risposta:* Significa che `tan(phi) = lam_vy(t)/lam_vx` e' una funzione lineare del
> tempo: e' la **legge linear-tangent** (bilinear-tangent). E' la firma
> caratteristica dei problemi di ascesa senza atmosfera e con gravita' costante, e
> vale perche' la dinamica non dipende dalla posizione (`lam_y = cost`) e
> l'accelerazione di gravita' e' costante. Non appena si aggiunge il drag, o una
> gravita' che varia con la quota, la linearita' si perde.

---

## Legge di controllo ottima (righe 23-24)

```matlab
% Optimal thrust angle
phi = atan2(lam_vy, lam_vx);
```

**Derivazione.** Raccogliendo i termini in `phi` nell'Hamiltoniana:

    H = [lam_x*vx + lam_y*vy - lam_vy - lam_m*Q]
        + (T/m) * [ lam_vx*cos(phi) + lam_vy*sin(phi) ]

Il PMP dice: `phi*` **massimizza** `H` puntualmente. Poiche' `T/m > 0` (spinta
positiva, massa positiva), massimizzare `H` equivale a massimizzare

    g(phi) = lam_vx*cos(phi) + lam_vy*sin(phi)

che e' il prodotto scalare tra il versore di spinta `[cos(phi), sin(phi)]` e il
vettore `lam_v = [lam_vx, lam_vy]`. Un prodotto scalare tra un versore e un vettore
fisso e' massimo quando il versore e' **allineato** al vettore:

    [cos(phi*), sin(phi*)] = lam_v / |lam_v|      ->   g_max = +|lam_v|

La condizione di stazionarieta' `dH/dphi = 0` da' `tan(phi*) = lam_vy/lam_vx`, che ha
**due radici per periodo** (differiscono di `pi`). La radice sbagliata da'
`g = -|lam_v|`, cioe' il **minimo** di `H`: spinta anti-allineata al primer vector,
razzo che frena. `atan2(lam_vy, lam_vx)` -- a differenza di `atan(lam_vy/lam_vx)` --
usa i **segni separati** dei due argomenti per selezionare il quadrante giusto,
quindi sceglie automaticamente il ramo di massimo. Usare `atan` qui sarebbe un bug
silenzioso ogni volta che `lam_vx < 0`.

- **Attenzione al segno (punto delicato, spesso confuso).** In molti testi la
  direzione di spinta ottima si scrive `-lambda_v/|lambda_v|`. Qui **no**: il codice
  usa `+lam_v/|lam_v|`. Non e' un errore, e' la convenzione. In HM1 il problema e'
  posto come **massimizzazione** di `m(tf)` con il funzionale aumentato in forma
  `(f - x_dot)`; questa scelta da' `lam_m(tf) = +1` e `lam_m_dot = +(T/m^2)*|lam_v| >= 0`,
  e la spinta va **con** `lam_v`. Con la convenzione opposta (minimizzazione del
  consumo, funzionale in forma `(x_dot - f)`) tutti i costati cambiano segno e la
  spinta risulta lungo `-lambda_v`. **Le due formulazioni danno la stessa
  traiettoria fisica**; quello che conta e' che segno dell'Hamiltoniana, segno di
  `lam_m_dot` e segno del controllo siano coerenti tra loro -- e nel codice lo sono
  (righe 24 e 33 usano entrambe il `+`).
- **Vincolo di controllo assente.** `phi` e' libero su tutto il cerchio: nessuna
  saturazione, nessun limite di rateo di beccheggio, nessun vincolo di angolo
  d'attacco. E' questo che rende la massimizzazione dell'Hamiltoniana un problema
  algebrico banale. Con un vincolo del tipo `|phi - psi| <= alpha_max` la soluzione
  diventerebbe un problema con archi vincolati (le cose diventano *molto* piu'
  difficili).

> **Possibile domanda d'esame** -- Perche' `atan2` e non `atan`?
> *Risposta:* Perche' la condizione di stazionarieta' `tan(phi) = lam_vy/lam_vx`
> determina `phi` solo modulo `pi`, e delle due radici una massimizza e l'altra
> minimizza l'Hamiltoniana. `atan` restituisce sempre un angolo in `(-pi/2, pi/2)`
> perche' vede solo il rapporto e perde i segni individuali; `atan2` vede i segni di
> numeratore e denominatore separatamente e restituisce il quadrante corretto,
> selezionando il ramo che allinea la spinta con il primer vector (il massimo).
> In Task 1 con `lam_vx0 > 0` i due coinciderebbero, ma la funzione e' condivisa da
> tutti i task e non c'e' garanzia a priori sul segno di `lam_vx0`.

---

## Derivate dello stato (righe 26-33)

```matlab
dz = zeros(6,1);
dz(1) = vx;                          % dx/dt
dz(2) = vy;                          % dy/dt
dz(3) = (p.T / m) * cos(phi);        % dvx/dt
dz(4) = (p.T / m) * sin(phi) - 1;    % dvy/dt  (g = 1 nondim)
dz(5) = -p.Q;                         % dm/dt
dz(6) = (p.T / m^2) * lam_v_norm;    % dlam_m/dt
```

- Riga 27: preallocazione `zeros(6,1)`. `ode45` pretende un vettore **colonna**;
  restituire una riga fa fallire l'integrazione.
- Righe 28-29: cinematica. Banali, ma sono loro a rendere `x` e `y` non presenti nel
  RHS (vedi sopra).
- Riga 30: `dvx/dt = (T/m)*cos(phi)`. Nessuna forza orizzontale oltre alla spinta:
  niente drag, niente Coriolis (Terra non rotante).
- Riga 31: `dvy/dt = (T/m)*sin(phi) - 1`. **Il `-1` e' la gravita'.** Nello schema di
  nondimensionalizzazione `a_rif = g`, quindi `g/a_rif = 1` esatto. Se si togliesse
  quel `-1` sparirebbe la gravity loss, `lam_vy` resterebbe comunque lineare (perche'
  `-dH/dvy` non contiene `g`: la gravita' e' un termine costante nell'Hamiltoniana,
  non dipende dallo stato) ma la traiettoria non salirebbe piu' contro nulla e il
  problema perderebbe il suo trade-off centrale. **Sottigliezza:** il termine `-1`
  entra in `H` come `-lam_vy` (vedi il termine `lam_vy*(-1)`), ed e' esattamente
  quel `-lam_vy` che compare nella condizione `H(0) = -lam_vy0 + T*(|lam_v0| - 1/c)`
  usata come quarto residuo dello shooting.
- Riga 32: `dm/dt = -Q`, con `Q` **costante**. Conseguenza notevole: la massa e'
  **lineare** nel tempo, `m(t) = m0 - Q*t = 1 - Q*t`. Quindi `mf = 1 - Q*tf` in forma
  chiusa. Il codice invece la legge da `Z(end,5)` dopo l'integrazione -- funziona, ed
  e' anche un buon test di consistenza dell'integratore, ma e' informazione ridondante.
  Nota: **non c'e' nessun guard su `m > 0`**. Con `Q = 7` (estremo alto dello sweep)
  basta `tf > 1/7 ~ 0.143` per far diventare la massa negativa e `T/m` singolare.
  In pratica non succede perche' la soluzione ottima ha `tf` ben sotto quella soglia,
  ma la protezione e' assente per costruzione.
- Riga 33: `dlam_m/dt = (T/m^2)*|lam_v|`. **Derivazione del segno:**

      dH/dm = lam_vx*T*(-1/m^2)*cos(phi) + lam_vy*T*(-1/m^2)*sin(phi)
            = -(T/m^2) * [lam_vx*cos(phi) + lam_vy*sin(phi)]
            = -(T/m^2) * |lam_v|        (sostituendo phi ottimo)

  e quindi `lam_m_dot = -dH/dm = +(T/m^2)*|lam_v| >= 0`. **`lam_m` cresce sempre**
  lungo un arco propulso (e' costante solo dove `T = 0`, cioe' in coast).
  Interpretazione: `lam_m` e' la sensibilita' del costo (massa finale) alla massa
  istantanea; piu' si e' leggeri, piu' ogni chilo in piu' "vale", perche' `T/m`
  e' maggiore.
- **Nota di onesta' importante:** in `main_task1.m` (e in Task 2 e 4) `z(6)` viene
  integrato ma il suo valore **non entra mai nel residuo dello shooting**. Il motivo
  e' la normalizzazione `lam_m0 = 1`: fissando la scala dei costati all'istante
  iniziale si *rinuncia* alla condizione di trasversalita' `lam_m(tf) = 1`. Quindi
  `Z(end,6)` a fine integrazione **non vale 1** (vale piu' di 1, dato che `lam_m`
  cresce): e' semplicemente un multiplo positivo dei costati "veri". Questo e'
  innocuo perche' la traiettoria dipende dai costati solo tramite la **direzione**
  `phi` (un rapporto, invariante di scala) e perche' `H` e' omogenea di grado 1 in
  `lambda`, quindi `H = 0` e' anch'essa invariante di scala. `z(6)` e' portato in
  giro perche' **serve al Task 3**, dove la normalizzazione `lam_m0 = 1` viene
  abbandonata: `lam_m0` diventa una **quinta incognita vera** dello shooting
  (`z0(4)` in `main_task3.m`) e il residuo di trasversalita' `lam_m(tc) = 1` legge
  direttamente `zf(6)`, cioe' il valore integrato di `z(6)` al cutoff
  (`main_task3.m` riga 302). E' **questo residuo** -- non la funzione di switching --
  a rendere indispensabile l'integrazione di `z(6)`: la switching function teorica e'
  effettivamente `S = |lam_v|/m - lam_m/c`, ma nel codice, una volta che il residuo ha
  fissato `lam_mc = 1`, viene implementata nella forma equivalente
  `S = |lam_v|/mc - 1/c` (riga 298), in cui `lam_m` **non compare**. La scelta e'
  deliberata e commentata nel sorgente: la forma letterale `lam_mc/c` fa scivolare
  `fsolve` sulla radice spuria con `vyc < 0`, mentre le due forme coincidono al
  cutoff.
- **Nota:** `p.c` e' documentato nell'header (riga 6) come campo della struct ma
  `ode_burn` **non lo usa mai**. La relazione `T = c*Q` viene applicata dal
  chiamante (`p.T = c * p.Q` in `main_task1.m`). Il commento e' quindi leggermente
  fuorviante: `p.c` e' richiesto dal *contratto complessivo* della pipeline, non da
  questa funzione.

> **Possibile domanda d'esame** -- Il sistema e' autonomo (l'Hamiltoniana e'
> costante), eppure `ode_burn` dipende esplicitamente da `t`. Contraddizione?
> *Risposta:* No. Il sistema **stato + costati** e' autonomo, ed e' per questo che
> `H = cost = 0`. La dipendenza esplicita da `t` in `ode_burn` e' un **artefatto
> dell'aver risolto in forma chiusa l'equazione di `lam_vy`**: invece di integrare
> `lam_vy_dot = -lam_y` come settima equazione, si sostituisce direttamente
> `lam_vy(t) = lam_vy0 - lam_y*t`. E' una riduzione esatta, non un'approssimazione:
> l'informazione temporale che vediamo e' quella del costato eliminato.

---

## Possibili domande d'esame

**D: Da dove viene esattamente `dz(6) = (T/m^2)*|lam_v|` e perche' e' sempre positiva?**
R: E' l'equazione di Eulero-Lagrange per il costato di massa, `lam_m_dot = -dH/dm`.
Nell'Hamiltoniana la massa compare solo nei termini di spinta `lam_vx*(T/m)*cos(phi)`
e `lam_vy*(T/m)*sin(phi)`; derivando rispetto a `m` esce `-(T/m^2)*(lam_vx*cos + lam_vy*sin)`,
e sostituendo il controllo ottimo la parentesi vale `|lam_v|`. Il doppio segno meno
(`-dH/dm` di un termine gia' negativo) lascia il segno `+`. E' sempre >= 0 perche'
`T >= 0`, `m > 0` e `|lam_v| >= 0`: lungo un arco propulso `lam_m` cresce
monotonicamente, cioe' il valore marginale di un chilo di massa aumenta man mano che
il veicolo si alleggerisce.

**D: Come si passa da "massimizzare H" a `phi = atan2(lam_vy, lam_vx)` in due righe?**
R: Solo due termini di `H` contengono `phi`, e si raccolgono come
`(T/m)*(lam_vx*cos(phi) + lam_vy*sin(phi))`, cioe' `(T/m)` per il prodotto scalare tra
il versore di spinta e il vettore `lam_v`. Con `T/m > 0` massimizzare `H` significa
massimizzare quel prodotto scalare, che e' massimo quando il versore di spinta e'
allineato con `lam_v`. Quindi `cos(phi*) = lam_vx/|lam_v|`, `sin(phi*) = lam_vy/|lam_v|`,
cioe' `phi* = atan2(lam_vy, lam_vx)`, e il valore massimo del termine di controllo e'
`(T/m)*|lam_v|`. Il vettore `lam_v` e' il **primer vector** di Lawden.

**D: Perche' questa funzione non ha un blocco `arguments` di validazione, quando il
resto della repo lo usa?**
R: Perche' e' una hot-loop function. Sta dentro `ode45`, dentro `fsolve`, dentro il
doppio loop di continuazione su `Q` (80 punti) e `yf` (3 valori): l'ordine di
grandezza e' 1e6+ chiamate per run. Con `RelTol = 1e-10` l'integratore fa migliaia di
step per traiettoria, e ogni step richiede piu' valutazioni del RHS (Runge-Kutta 4/5
ne fa 6). Un blocco `arguments` aggiungerebbe overhead a ogni chiamata senza dare
nessuna informazione nuova (gli input arrivano sempre dallo stesso chiamante interno).
La validazione e' collocata al boundary del run, una volta sola. E' una scelta
documentata nel `CLAUDE.md` della repo.

**D: Cosa succede se aggiungi il drag aerodinamico a questo RHS?**
R: Si rompe tutta la struttura analitica. Il drag dipende dalla densita' `rho(y)` e
dalla velocita', quindi `dH/dy ~= 0` e `dH/dvx, dH/dvy` acquistano termini nuovi.
Conseguenze a cascata: `lam_y` non e' piu' costante, quindi `lam_vy` non e' piu'
lineare in `t`, quindi **la legge linear-tangent non vale piu'**. Bisognerebbe
integrare tutti e cinque i costati numericamente (stato aumentato 10x1 invece di 6x1),
il controllo resterebbe `phi = atan2(lam_vy, lam_vx)` solo se la spinta e' l'unica
forza controllata e il drag non dipende da `phi` -- cosa che in generale e' falsa
(il drag dipende dall'angolo d'attacco). Il costo per chiamata almeno raddoppierebbe
e il basin di convergenza dello shooting si restringerebbe ulteriormente.

**D: `z(6)` (cioe' `lam_m`) viene integrato ma in Task 1 non serve. Perche' non toglierlo?**
R: Perche' `ode_burn.m` e' condiviso da tutti e quattro i task, e il **Task 3 ne ha
bisogno**: li' la normalizzazione `lam_m0 = 1` viene abbandonata, `lam_m0` diventa
una quinta incognita vera dello shooting e il residuo di trasversalita'
`lam_m(tc) = 1` legge direttamente `zf(6)` (`main_task3.m` riga 302) -- senza
integrare `z(6)` quel residuo non esisterebbe. Attenzione a non raccontarla male
all'orale: la funzione di switching teorica e' `S = |lam_v|/m - lam_m/c`, ma nel
codice, dopo che il residuo ha fissato `lam_mc = 1`, e' implementata nella forma
equivalente `S = |lam_v|/mc - 1/c` (riga 298), dove `lam_m` non compare -- la forma
letterale `lam_mc/c` porterebbe `fsolve` sulla radice spuria `vyc < 0`.
Portare 6 equazioni invece di 5 in Task 1 costa pochissimo (una
quadratura accodata, accoppiata in una sola direzione: nessun altro `dz(i)` dipende da
`z(6)`) ed evita di duplicare il RHS. Va pero' detto chiaramente che in Task 1 il
valore di `z(6)` a `tf` **non soddisfa** la trasversalita' `lam_m(tf) = 1`, perche'
quella condizione e' stata sostituita dalla normalizzazione `lam_m0 = 1`.

**D: Perche' la spinta e' lungo `+lambda_v` e non `-lambda_v` come in molti libri?**
R: E' una questione di convenzione, non di fisica. HM1 pone il problema come
**massimizzazione** di `m(tf)` e costruisce il funzionale aumentato con
`lambda^T*(f - x_dot)`. Con questa scelta la trasversalita' da' `lam_m(tf) = +1`, la
PMP richiede di **massimizzare** `H`, e la spinta risulta allineata con `+lam_v`. Con
la convenzione duale (minimizzare il consumo, funzionale con `lambda^T*(x_dot - f)`)
tutti i costati cambiano segno e la spinta va lungo `-lambda_v`. Le due traiettorie
ottime coincidono. L'unica cosa che conta e' la **coerenza interna**: nel codice il `+`
di riga 24 (controllo) e il `+` di riga 33 (`lam_m_dot`) appartengono alla stessa
convenzione, e la stessa convenzione e' usata in `H(0) = -lam_vy0 + T*(|lam_v0| - 1/c)`
nel residuo dello shooting.
