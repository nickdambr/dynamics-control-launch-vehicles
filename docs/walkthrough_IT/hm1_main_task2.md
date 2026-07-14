# HM1/main_task2.m

## Ruolo del file nel progetto

`main_task2.m` risolve la **seconda variante** della missione di HM1: prima di
accendere il gravity turn ottimo, il lanciatore deve compiere una **salita
verticale** fino a una quota assegnata `y1 = 1e-4` (adimensionale; circa 620 m
dimensionali secondo il README di HM1, con `L_rif = V_rif^2/g ~ 6200 km`). La
missione diventa quindi a due fasi: `salita verticale -> arco propulso ottimo`.
Lo scopo dichiarato nell'intestazione (righe 1-3) e' confrontare il payload con
quello del Task 1 (nessuna salita verticale).

La differenza sostanziale rispetto a `main_task1.m` **non e' nella matematica
del problema di controllo ottimo**, ma nel **punto di partenza dell'arco
ottimizzato**. Nel Task 1 l'arco propulso parte da `(x,y,vx,vy,m) = (0,0,0,0,1)`;
qui parte da `(0, y1, 0, vy1, m1)`, dove `vy1` e `m1` sono il risultato
dell'integrazione della salita verticale. La fase 1 **non ha gradi di liberta'**
(l'angolo di spinta e' imposto a `phi = 90 deg` e la durata e' fissata
dall'evento `y = y1`), quindi non c'e' nulla da ottimizzare in fase 1 e i suoi
costati sono irrilevanti: non esiste nessuna condizione di continuita' dei
costati alla giunzione. Solo lo **stato** viene trasferito da fase 1 a fase 2.

Di conseguenza il vettore delle incognite dello shooting resta **identico** a
quello del Task 1 -- `[lam_vx0, lam_vy0, lam_y, t_burn]`, 4 incognite e 4
residui -- e cambia solo (a) la condizione iniziale passata a `ode45`, e (b) la
forma algebrica della condizione `H(0) = 0`, che ora contiene i termini
`lam_y*vy1` e `T/m1` (nel Task 1 valevano zero e `T`, perche' `vy0 = 0` e
`m0 = 1`).

Il file dipende da `ode_burn.m` (RHS condiviso stato+costato dell'arco propulso)
e definisce tre funzioni locali: `ode_vertical`, `event_altitude`, `shooting2`.
Non e' chiamato da nessun altro script: e' un entry point standalone che, come
tutti i `main_task*.m` di HM1, esporta le figure in `HM1/figures/` con prefisso
`task2_`.

---

## Richiamo: le condizioni necessarie (PMP) usate dall'arco propulso

Sono le stesse per tutti e quattro i task e sono cablate dentro `ode_burn.m`.
Dinamica adimensionale (piana, Terra piatta, senza drag, `g = 1`, `T = c*Q`
costante):

    x_dot  = vx
    y_dot  = vy
    vx_dot = (T/m)*cos(phi)
    vy_dot = (T/m)*sin(phi) - 1
    m_dot  = -Q

Costo: massimizzare `m(tf)`. Hamiltoniana (convenzione del **massimo**, quella
usata dal codice):

    H = lam_x*vx + lam_y*vy + lam_vx*(T/m)*cos(phi)
        + lam_vy*((T/m)*sin(phi) - 1) - lam_m*Q

Equazioni dei costati (`lam_dot = -dH/dstato`):

    lam_x_dot  = 0
    lam_y_dot  = 0
    lam_vx_dot = -lam_x
    lam_vy_dot = -lam_y
    lam_m_dot  = +(T/m^2)*|lam_v|          con |lam_v| = sqrt(lam_vx^2+lam_vy^2)

Conseguenze usate dal codice:

- `x(tf)` e' **libera** (nessuna condizione terminale su `x`) -> trasversalita'
  `lam_x(tf) = 0` -> essendo `lam_x` costante, `lam_x == 0` ovunque -> `lam_vx`
  e' **costante** `= lam_vx0`.
- `lam_y` e' costante (incognita), quindi `lam_vy(t) = lam_vy0 - lam_y*t` e'
  **lineare** nel tempo. Questo e' esattamente cio' che fa `ode_burn.m` (riga 20).
- Massimizzando `H` rispetto a `phi`: `cos(phi) = lam_vx/|lam_v|`,
  `sin(phi) = lam_vy/|lam_v|`, cioe'

      tan(phi) = (lam_vy0 - lam_y*t) / lam_vx0

  cioe' la **legge della tangente lineare** (caso degenere della bilinear-tangent
  perche' `lam_x = 0`).
- `m(tf)` e' libera e coincide col costo -> `lam_m(tf) = 1`.
- `tf` libero + sistema autonomo -> `H(t) = 0` **per ogni t** (non solo in `tf`).

---

## Intestazione e parametri (righe 1-17)

```matlab
c   = 0.6;
eta = 0.1;
y1  = 0.0001;   % vertical climb altitude
yf  = 0.04;     % target orbit altitude
Q   = 2;
```

- Righe 8-12: dati fissi. `c = 0.6` velocita' efficace di scarico adimensionale,
  `eta = 0.1` coefficiente strutturale `ms/mp`, `yf = 0.04` quota di iniezione.
  `Q = 2` e' il valore nominale della traccia -- **non** il `Q*` ottimo trovato
  nel Task 1 (che per `yf = 0.04` vale circa 2.52 secondo il README): il Task 2
  non riesegue lo sweep su `Q`.
- Riga 10: `y1 = 1e-4` e' un **dato di progetto imposto**, non una variabile di
  ottimizzazione. Se lo si volesse ottimizzare servirebbe una condizione di
  interior-point sui costati alla giunzione -- qui non serve proprio perche' `y1`
  e' fisso.
- Righe 14-17: tolleranze. `ode45` a `RelTol=1e-10, AbsTol=1e-12` e `fsolve` a
  `1e-10` su funzione e passo: e' la ricetta standard del metodo indiretto, dove
  i residui dello shooting sono estremamente sensibili ai costati iniziali (il
  sistema aggiunto e' instabile in avanti) e un'integrazione lasca produce
  gradienti numerici privi di senso.

---

## Fase 1 -- salita verticale (righe 19-37)

```matlab
T = c * Q;
ic_vert = [0; 0; 1];
opts_vert = odeset('RelTol', 1e-12, 'AbsTol', 1e-14, ...
    'Events', @(t,z) event_altitude(t, z, y1));
[T_vert, Z_vert] = ode45(@(t,z) ode_vertical(t,z,T,Q), ...
                         [0 1], ic_vert, opts_vert);
```

- Riga 20: `T = c*Q = 1.2`. In unita' adimensionali il peso iniziale e'
  `m0*g = 1`, quindi `T` **coincide numericamente col rapporto T/W al liftoff**:
  per questo la riga 22 stampa `T` sotto l'etichetta `T/W`. E' corretto solo a
  `t = 0`; man mano che `m` cala, `T/W` cresce. Condizione di decollo:
  `T > 1` cioe' `Q > 1/c = 1.667` (con `Q = 2` siamo dentro).
- Righe 24-25: stato ridotto `z = [y; vy; m]` -- durante la salita verticale `x` e
  `vx` restano identicamente nulli, quindi si integra solo il sottospazio utile.
- Righe 26-27: tolleranze ancora **piu' strette** (`1e-12/1e-14`) che nell'arco
  propulso. Motivo: `vy1` e `m1` sono le condizioni iniziali del BVP a valle;
  un errore su di esse si propaga in un errore sui residui, ed e' un errore che
  `fsolve` non puo' compensare (non e' una delle sue incognite).
- Riga 28: `tspan = [0 1]` con evento terminale. **Punto fragile**: il codice non
  verifica che l'evento sia effettivamente scattato. Se `T <= 1` (`Q < 1/c`) il
  veicolo non decolla, l'evento non si attiva mai e lo script prosegue
  silenziosamente con `t1 = 1` e uno stato assurdo. Non c'e' guardia.
- Righe 30-33: si estraggono `t1`, `y_1 = Z_vert(end,1)` (la quota **davvero
  raggiunta**, non `y1` nominale), `vy1`, `m1`. Nota: `main_task3.m` in questo
  punto usa invece `y1` nominale; e' un'incoerenza tra i due script, di impatto
  numerico trascurabile (dell'ordine della tolleranza dell'evento).

> **Possibile domanda d'esame** -- perche' la salita verticale non produce
> nessuna condizione sui costati alla giunzione `t1`?
> *Risposta:* perche' in fase 1 non ci sono gradi di liberta': `phi` e' imposto a
> `90 deg` e la durata `t1` e' determinata dal vincolo `y(t1) = y1`. I costati
> servono a esprimere la stazionarieta' rispetto a un controllo libero; se il
> controllo non e' libero, i moltiplicatori della fase 1 non entrano in nessuna
> condizione necessaria. La giunzione trasferisce quindi solo lo **stato**
> `(0, y1, 0, vy1, m1)`, e i costati dell'arco propulso ripartono come incognite
> nuove di zecca. Sarebbero servite condizioni di interior-point solo se `y1` (o
> `t1`) fosse stata una variabile da ottimizzare.

---

## Fase 2 -- impostazione e soluzione del BVP (righe 39-64)

```matlab
p.x0  = 0;   p.y0  = y_1;
p.vx0 = 0;   p.vy0 = vy1;   p.m0  = m1;

z_guess = [0.6; 3.8; 14; 0.30];
[z_sol, ~, ef] = fsolve(@(z) shooting2(z,p,opts_ode), ...
                        z_guess, opts_fs);
```

- Righe 44-52: la struct `p` porta i parametri fisici **e** lo stato iniziale
  dell'arco propulso. E' l'unico canale attraverso cui la fase 1 comunica con la
  fase 2.
- Righe 41-42 (commento): incognite `[lam_vx0, lam_vy0, lam_y, t_burn]`, residui
  `[y(tf)-yf, vx(tf)-1, vy(tf), H(0)]`. **Identico al Task 1 come cardinalita'**:
  la normalizzazione `lam_m0 = 1` rimuove una incognita e la condizione
  `lam_m(tf) = 1` viene sostituita da essa (vedi sotto, `shooting2`).
- Riga 55: `z_guess = [0.6; 3.8; 14; 0.30]` e' **letteralmente lo stesso guess
  cablato in `main_task1.m` (riga 40)**. Questo e' l'unico "warm start" presente:
  un travaso manuale del guess da un task all'altro. **Non c'e' continuazione**
  in questo script -- nessun loop che riparte dalla soluzione precedente. Funziona
  perche' `y1 = 1e-4` e' minuscolo (1/400 della quota finale), quindi la
  soluzione del Task 2 e' una perturbazione piccola di quella del Task 1. Nota
  pero' che il guess era tarato su `Q ~ 3` (nel Task 1 lo sweep parte dal `Q` piu'
  vicino a 3) mentre qui `Q = 2`.
- Righe 60-64: se `fsolve` fallisce (`ef <= 0`) si riprova **una sola volta** con
  un secondo guess `[0.4; 3.0; 10; 0.35]`. E' un multi-start rudimentale, non una
  continuazione: se anche il secondo fallisce, lo script stampa un errore
  (riga 155) e non produce risultati.
- Lettura fisica del guess: con `lam_vy0 = 3.8` e `lam_y = 14`, dopo
  `t = 0.30` si ha `lam_vy = 3.8 - 14*0.30 = -0.4`, quindi
  `phi = atan2(lam_vy, 0.6)` scende da circa `+81 deg` a circa `-34 deg`: e'
  proprio la pitch-over di un gravity turn. Un guess con `lam_y < 0` darebbe una
  spinta che si alza nel tempo -- non convergerebbe mai.

> **Possibile domanda d'esame** -- se il metodo indiretto e' cosi' sensibile,
> perche' qui basta un guess cablato e nel Task 1 serviva la continuazione?
> *Risposta:* nel Task 1 bisogna risolvere il BVP per 80 valori di `Q` e 3 valori
> di `yf`: quando `Q` si allontana dal valore "amichevole" il bacino di
> convergenza si restringe e un guess fisso fallisce, quindi si usa la
> continuazione (warm start dalla soluzione al `Q` vicino). Nel Task 2 c'e' **una
> sola** soluzione da trovare, a `Q = 2` fissato, e la salita verticale e' cosi'
> corta che la soluzione e' vicinissima a quella del Task 1: il guess del Task 1
> cade gia' dentro il bacino di attrazione.

---

## Estrazione della soluzione e risultati (righe 66-87)

- Righe 68-71: i tre costati vengono ricopiati dentro `pp`, che e' la struct che
  `ode_burn` si aspetta. La riga 72 estrae invece `t_burn` come variabile locale
  ordinaria: **non** e' un campo di `pp` (`ode_burn` legge solo `T`, `Q`, `c`,
  `lam_vx0`, `lam_vy0`, `lam_y`), serve come estremo del `tspan` alla riga 75.
- Riga 74: `ic2 = [x0; y0; vx0; vy0; m0; 1]` -- il sesto stato e' `lam_m`, inizializzato
  a **1**: e' la normalizzazione, non una condizione fisica.
- Riga 75: reintegrazione dell'arco su una griglia densa `linspace(0,t_burn,500)`
  solo per i grafici (la soluzione era gia' nota da `fsolve`).
- Riga 78: la durata totale e' `tf_total = t1 + t_burn` (le righe 80-84 si
  limitano a stamparla insieme agli altri risultati). `t1` **non** e' un'incognita
  dello shooting: e' un dato prodotto dalla fase 1.
- Riga 87: payload `mu = mf*(1+eta) - eta`. Derivazione: a fine missione la massa
  residua e' `mf = ms + mu`, con struttura `ms = eta*mp` e propellente
  `mp = m0 - mf = 1 - mf` (adimensionale, `m0 = 1`). Sostituendo:
  `mu = mf - eta*(1 - mf) = mf*(1+eta) - eta`. E' la stessa formula del Task 1,
  quindi il confronto di payload richiesto dall'intestazione e' un confronto
  diretto di `mf`.

---

## Plots (righe 89-153)

- Righe 99-106: traiettoria `y(x)`: rosso il tratto verticale (`x == 0`), blu
  l'arco propulso.
- Righe 108-126: **inset di zoom**. Serve perche' `y1 = 1e-4` e' 1/400 di
  `yf = 0.04`: a scala piena la salita verticale e' un punto. L'inset mostra il
  ginocchio di pitch-over. Il rettangolo tratteggiato (righe 117-118) marca la
  regione ingrandita.
- Righe 129-140: ricostruzione degli angoli lungo l'arco propulso.
  `lam_vy_k = z_sol(2) - z_sol(3)*T2(kk)` ricalcola il costato **analiticamente**
  (e' lineare: non serve integrarlo), poi `phi = atan2(lam_vy_k, lam_vx0)` e
  `psi = atan2(vy, vx)` (angolo del vettore velocita'). La guardia `V_k > 1e-10`
  evita l'`atan2(0,0)` al primo istante, dove `vx = 0` e `vy = vy1` (in realta'
  `vy1 > 0`, quindi `V_k` non e' nullo: la guardia e' un residuo del Task 1, dove
  la velocita' iniziale e' davvero zero).
- Righe 146-147: nella fase verticale si tracciano `phi = 90 deg` e `psi = 90 deg`
  sovrapposte. Il salto visibile in `t1` tra `90 deg` e `phi(0) = atan2(lam_vy0,
  lam_vx0)` e' **fisico e atteso**: il controllo e' discontinuo alla giunzione
  perche' la fase 1 e' imposta, non ottima. Solo se `y1 -> 0` il salto si
  annulla.

---

## Export figure (righe 158-174)

- Righe 159-160: cartella `HM1/figures/`, creata se manca.
- Righe 166-171: `theme(fig,'light')` forza il tema chiaro (altrimenti in dark
  mode le PNG escono con sfondo nero e finiscono cosi' nel report LaTeX);
  fallback su `Color='w'` per MATLAB pre-R2025a.
- Riga 172-173: `exportgraphics` a 200 dpi, nome `task2_<nome-figura-slugificato>.png`.

---

## `ode_vertical` (righe 178-187)

```matlab
function dz = ode_vertical(~, z, T, Q)
    vy = z(2); m = z(3);
    dz = [vy; T/m - 1; -Q];
end
```

- Riga 178: firma. Chiamata solo da `ode45` alla riga 28.
- Riga 186: e' la dinamica con `phi = pi/2` **sostituito a mano**:
  `vx_dot = (T/m)*cos(90 deg) = 0` (quindi `x` e `vx` non compaiono),
  `vy_dot = (T/m)*sin(90 deg) - g = T/m - 1` (con `g = 1`), `m_dot = -Q`.
- Non ci sono costati: la fase non e' ottimizzata.
- Questa funzione e' **duplicata identica** in `main_task3.m` (righe 221-230).
  Nessuna delle due e' fattorizzata in un file condiviso come `ode_burn.m`.

---

## `event_altitude` (righe 189-199)

```matlab
value = z(1) - y_target;
isterminal = 1;
direction = 1;
```

- Riga 196: la funzione evento si annulla quando `y = y1`.
- Riga 197: `isterminal = 1` -> l'integrazione si ferma.
- Riga 198: `direction = 1` -> si intercetta solo l'attraversamento **in salita**.
  Se il veicolo non decollasse (`T < 1`), `y` non salirebbe mai e l'evento non
  scatterebbe: da qui la fragilita' segnalata sopra.
- Anche questa e' duplicata in `main_task3.m` (righe 232-242).

---

## `shooting2` (righe 201-247)

E' il cuore del file: la funzione residuo che `fsolve` azzera.

- Riga 201: firma `res = shooting2(z0, p, opts_ode)`. `z0` sono le 4 incognite,
  `p` porta parametri **e** stato iniziale dell'arco. Nessun blocco `arguments`:
  come dichiarato ai commenti (riga 211) e nel `CLAUDE.md` del repo, la funzione
  gira dentro il loop di `fsolve` (migliaia di chiamate) e la validazione starebbe
  fuori posto.
- Righe 213-216: unpack delle incognite `[lam_vx0; lam_vy0; lam_y; t_burn]`.
- Righe 218-221: **box constraint via penalita'**. Se `t_burn <= 0` o `t_burn > 2`
  si restituisce `1e6*ones(4,1)`. Serve a impedire a `fsolve` (che e'
  unconstrained) di esplorare durate assurde o negative. Effetto collaterale
  onesto da dichiarare: il residuo diventa **discontinuo** su quella frontiera, e
  la Jacobiana per differenze finite che `fsolve` costruisce li' e' priva di
  senso. Funziona in pratica perche' il guess parte lontano dalla frontiera.
- Riga 228: `ic = [x0; y0; vx0; vy0; m0; 1]` -- di nuovo `lam_m0 = 1`.
- Righe 230-236: `try/catch` attorno a `ode45`; un fallimento dell'integratore
  produce residui `1e6` invece di far crashare lo sweep.
- Righe 240-241: **la condizione di Hamiltoniana nulla**, il pezzo di matematica
  specifico di questo task.

```matlab
Lam0 = sqrt(lam_vx0^2 + lam_vy0^2);
H0 = lam_y * p.vy0 + (p.T / p.m0) * Lam0 ...
     - lam_vy0 - p.T / p.c;
```

  Derivazione. Valutiamo `H` all'inizio dell'arco propulso, dove lo stato e'
  **noto** (`x=0, y=y1, vx=0, vy=vy1, m=m1`) e i costati valgono
  `lam_x = 0`, `lam_vx = lam_vx0`, `lam_vy = lam_vy0`, `lam_m = 1`:

      H(0) = lam_x*vx0 + lam_y*vy0 + (T/m0)*|lam_v0|
             - lam_vy0*g - lam_m0*Q

  con `lam_x = 0`, `g = 1`, `lam_m0 = 1` e `Q = T/c` (perche' `T = c*Q`):

      H(0) = lam_y*vy1 + (T/m1)*sqrt(lam_vx0^2 + lam_vy0^2)
             - lam_vy0 - T/c

  che e' esattamente la riga 241. Il termine `(T/m0)*|lam_v0|` viene dal fatto
  che il massimo su `phi` di `lam_vx*cos(phi) + lam_vy*sin(phi)` e' `|lam_v|`.

  **Perche' in `t = 0` e non in `tf`?** Perche' il sistema e' autonomo (`H` non
  dipende esplicitamente dal tempo), quindi `H` e' **costante** lungo l'arco: la
  condizione di trasversalita' `H(tf) = 0` (tempo finale libero) equivale a
  `H(0) = 0`. Ma in `t = 0` lo stato e' esatto, mentre in `tf` sarebbe il
  risultato dell'integrazione numerica -- imporla in `0` la rende **algebrica**,
  senza errore di integrazione, e migliora sensibilmente il condizionamento
  dello shooting.

  **Confronto col Task 1** (`main_task1.m` riga 268):

      H(0)_task1 = -lam_vy0 + T*( |lam_v0| - 1/c )

  E' il caso particolare di questa formula con `vy0 = 0` (il primo termine
  sparisce) e `m0 = 1` (`T/m0 -> T`). Task 2 non introduce una condizione nuova:
  **generalizza** la stessa condizione a uno stato iniziale non banale.

- Righe 243-246: i 4 residui:

      res(1) = y(t_burn) - yf     % quota di iniezione
      res(2) = vx(t_burn) - 1     % velocita' orbitale (Vrif = vxf)
      res(3) = vy(t_burn)         % iniezione orizzontale
      res(4) = H(0)               % tempo finale libero

  Notare che `x(tf)` non compare: e' libera, ed e' precisamente la ragione per cui
  `lam_x = 0` e quindi `lam_vx` e' costante. Notare anche che `lam_m(tf)` non
  compare: grazie alla normalizzazione `lam_m0 = 1`, `lam_m` e' un puro quadratura
  (in `ode_burn.m`, `dz(6)` dipende solo da `m` e `|lam_v|`, e nessuna altra
  derivata dipende da `lam_m`), quindi **non entra mai nei residui**.

> **Possibile domanda d'esame** -- la trasversalita' corretta e' `lam_m(tf) = 1`;
> perche' il codice impone invece `lam_m(0) = 1`?
> *Risposta:* le condizioni necessarie (equazioni dei costati, legge del controllo,
> `H = 0`) sono **omogenee di grado 1** nel vettore dei costati, e la legge del
> controllo `phi = atan2(lam_vy, lam_vx)` e' addirittura **invariante di scala**.
> Quindi si puo' riscalare l'intero vettore dei costati per una costante positiva
> senza cambiare la traiettoria: fissare `lam_m(0) = 1` invece di `lam_m(tf) = 1`
> e' solo una scelta di normalizzazione, che pero' **elimina un'incognita**
> (`lam_m0`) e una condizione (`lam_m(tf) = 1`), portando lo shooting da 5x5 a
> 4x4. Attenzione al fine print: il riscalamento e' lecito solo con fattore
> **positivo**, cioe' la normalizzazione assume implicitamente `lam_m(0) > 0` alla
> soluzione. Non e' una proprieta' garantita in generale -- infatti nel Task 3
> `lam_m0` torna a essere un'incognita, e i guess usati li' sono **negativi**.

---

## Limiti noti / punti onesti da dichiarare

- **Nessuna continuazione.** Contrariamente al Task 1, qui c'e' solo un guess
  cablato (identico a quello del Task 1) + un fallback. Se il guess fallisse, non
  c'e' una strategia di recupero sistematica. Una continuazione naturale sarebbe
  su `y1`: partire da `y1 = 0` (soluzione del Task 1) e aumentarlo gradualmente
  fino a `1e-4`, warm-startando ogni volta.
- **Nessuna verifica che l'evento sia scattato** (righe 26-30): con `Q < 1/c`
  lo script produrrebbe numeri senza segnalare nulla.
- **`Q = 2` non e' il `Q*` ottimo** del Task 1: il confronto di payload
  Task 1 vs Task 2 e' quindi corretto solo se fatto a parita' di `Q` (e nel Task 1
  a `Q = 2` la soluzione esiste, essendo `Q_vec` che parte da 1.8).
- **Codice duplicato**: `ode_vertical` e `event_altitude` sono copiate identiche in
  `main_task3.m` invece di stare in un file condiviso.
- **La guardia `V_k > 1e-10`** (riga 135) e' inutile in questo task (all'inizio
  dell'arco `vy = vy1 > 0`): e' un residuo del Task 1.

---

## Possibili domande d'esame

**D: Qual e' esattamente la differenza tra Task 1 e Task 2 nel problema di
controllo ottimo?**
R: Nessuna nella struttura delle condizioni necessarie. Il problema ottimizzato e'
sempre e solo l'arco propulso, con le stesse 4 incognite
`[lam_vx0, lam_vy0, lam_y, t_burn]` e gli stessi 4 residui
`[y(tf)-yf, vx(tf)-1, vy(tf), H(0)]`. Cambia solo la **condizione iniziale**
dell'arco, che nel Task 1 e' `(0,0,0,0,1)` e nel Task 2 e' `(0, y1, 0, vy1, m1)`
prodotta dalla salita verticale. Come conseguenza cambia la forma algebrica di
`H(0)`, che acquisisce il termine `lam_y*vy1` e sostituisce `T` con `T/m1`.

**D: Perche' non serve imporre la continuita' dei costati alla giunzione
salita-verticale / arco propulso?**
R: Perche' la salita verticale non ha gradi di liberta': il controllo e' imposto
(`phi = 90 deg`) e la durata e' fissata dal vincolo `y(t1) = y1`. Non essendoci
un controllo libero, non ci sono condizioni di stazionarieta' da scrivere in fase
1 e i suoi moltiplicatori non compaiono in nessuna condizione necessaria. Se
invece `y1` (o `t1`) fosse una variabile da ottimizzare, servirebbe una condizione
di interior-point: continuita' di `lam_x, lam_vx, lam_vy, lam_m`, salto di `lam_y`
pari al moltiplicatore del vincolo interno `y(t1) - y1 = 0`, e continuita' di `H`.

**D: Perche' `H = 0` viene imposta all'inizio dell'arco e non alla fine?**
R: Il sistema e' autonomo, quindi `H` e' un integrale primo lungo l'arco: `H(0) =
H(tf)`. Imporla in `tf` significherebbe valutarla su uno stato ottenuto per
integrazione numerica (con il suo errore); imporla in `0` la rende un'equazione
**algebrica esatta** nelle incognite, perche' lo stato iniziale e' noto. Meglio
condizionata e piu' economica.

**D: Il codice usa la continuazione (warm start)?**
R: No. Il Task 2 usa un guess cablato `[0.6; 3.8; 14; 0.30]` -- che e' letteralmente
il guess iniziale del Task 1 -- piu' un unico fallback `[0.4; 3.0; 10; 0.35]`. La
continuazione vera e propria (warm start a catena su un parametro che varia) e'
usata solo nel Task 1, per lo sweep su `Q`. Qui non serve perche' c'e' una sola
soluzione da trovare e la salita verticale e' cosi' corta (`y1 = 1e-4` contro
`yf = 0.04`) che la soluzione e' una perturbazione piccola di quella del Task 1.

**D: Perche' `lam_vx` e' costante e `lam_vy` lineare?**
R: `lam_vx_dot = -lam_x` e `lam_vy_dot = -lam_y`, e sia `lam_x` che `lam_y` sono
costanti (`H` non dipende esplicitamente da `x` e `y`). Poiche' `x(tf)` e' libera,
la trasversalita' da' `lam_x(tf) = 0`, e quindi `lam_x == 0` ovunque:
`lam_vx = lam_vx0` costante. `lam_y` invece resta un'incognita (`y(tf) = yf` e'
vincolata), quindi `lam_vy(t) = lam_vy0 - lam_y*t` e' lineare. Da qui la legge
della **tangente lineare** `tan(phi) = (lam_vy0 - lam_y*t)/lam_vx0`.

**D: Cosa succede fisicamente al payload aggiungendo la salita verticale, e
perche'?**
R: Peggiora leggermente. La salita verticale e' un tratto in cui la spinta e'
allineata contro gravita' e non produce nessuna velocita' orizzontale: e' tutta
perdita gravitazionale pura, e nessuna traiettoria non ottima puo' battere
l'ottimo del Task 1. La si accetta per realismo (sicurezza della base di lancio,
clearance dalla torre, carichi aerodinamici); essendo `y1` piccolissima, la
penalita' e' contenuta.

**D: Il payload si calcola come `mf*(1+eta) - eta`. Da dove viene?**
R: A fine missione `mf = ms + mu`. Il modello strutturale del corso da'
`ms = eta*mp` con `mp = m0 - mf`. In adimensionale `m0 = 1`, quindi `mp = 1 - mf`
e `ms = eta*(1-mf)`. Sostituendo: `mu = mf - eta*(1-mf) = mf*(1+eta) - eta`.
Massimizzare `mf` equivale quindi a massimizzare `mu` (la mappa e' affine e
crescente), il che giustifica l'uso di `m(tf)` come funzionale di costo.
