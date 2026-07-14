# HM1/main_task1.m

## Ruolo del file nel progetto

`main_task1.m` e' lo **script di ingresso del Task 1** di HM1: l'ascesa planare a
**singolo arco propulso**, dal decollo da fermo fino allo stato di iniezione
orbitale, risolta con il metodo indiretto (PMP + single shooting). E' il piu'
semplice dei quattro programmi di missione ed e' quello che stabilisce
l'infrastruttura riusata da tutti gli altri: la struct di parametri `p`, il
residuo di shooting, la strategia di continuazione, le tolleranze.

Lo script risponde a tre domande della traccia, che sono anche le tre figure che
produce:
- **(a)** come varia la massa finale `mf` al variare della portata `Q`, per tre
  quote target `yf = 0.04, 0.05, 0.06`;
- **(b)** come si decompone la perdita di velocita' (misalignment `Wd` + gravita'
  `Wg`) in funzione di `Q`, a `yf = 0.04`;
- **(c)** che aspetto hanno traiettoria e legge di controllo alla `Q` ottima.

Il punto centrale e' che **`Q` non e' un dato ma un parametro di progetto**: c'e'
un massimo interno di `mf(Q)`. A `Q` bassa il razzo e' lento, la spinta e' appena
sopra il peso e le gravity loss divorano il delta-V; a `Q` alta il razzo brucia in
fretta ma con una spinta enorme e mal allineata alla velocita', e domina la
misalignment loss. Il compromesso e' `Q*`.

Lo script contiene, oltre al corpo principale, **tre local function** in fondo:
`shooting1` (il residuo del BVP), `set_costates` (utility), `ode_burn_losses` (una
copia estesa di `ode_burn` con due integratori di perdita accodati). Dipende da
`HM1/ode_burn.m`, che vive in un file a parte perche' e' condiviso con gli altri
task.

---

## Header e `%% Parameters` (righe 1-18)

```matlab
c   = 0.6;          % effective exhaust velocity (nondim)
eta = 0.1;          % structural coefficient ms/mp
yf_vec = [0.04, 0.05, 0.06];
Q_vec  = linspace(1.8, 7, 80);   % Q > 1/c ~ 1.667 for liftoff
```

- Riga 6: `clear; close all; clc` -- questo e' uno **script**, non una funzione. Le
  variabili finiscono nel workspace base. E' il motivo per cui il `close all` conta:
  l'export delle figure (righe 207-223) rastrella con `findobj(groot,'Type','figure')`
  **tutte** le figure aperte, quindi partire con la lavagna pulita e' necessario per
  non esportare roba di sessioni precedenti.
- Riga 9: `c = 0.6`. E' la velocita' di efflusso efficace **nondimensionale**. Lo
  schema di nondimensionalizzazione (documentato in `HM1/README.md` e nel report) e':

      a_rif = g = 9.81 m/s^2
      V_rif = vxf = 7800 m/s       (velocita' orbitale al target)
      m_rif = m0 = 1e6 kg
      L_rif = V_rif^2 / g  ~ 6200 km
      t_rif = V_rif / g    ~ 795 s

  Quindi `c = 0.6` significa `c_dim = 0.6 * 7800 = 4680 m/s`, cioe' `Isp ~ 477 s`.
  In queste unita' `g = 1`, `m0 = 1`, e la condizione terminale sulla velocita'
  orizzontale diventa semplicemente `vx(tf) = 1`.
- Riga 10: `eta = 0.1` -- coefficiente strutturale `ms/mp`. **Non entra nella
  dinamica**: e' usato solo *a posteriori* (riga 104) per convertire la massa finale
  in payload: `mu = mf*(1+eta) - eta`. Deriva da `mf = mu + ms` con `ms = eta*mp` e
  `mp = 1 - mf`, da cui `mu = mf - eta*(1-mf) = mf*(1+eta) - eta`. Conseguenza: il
  payload puo' venire **negativo** (e in effetti a `yf = 0.06` viene negativo, secondo
  i risultati riportati nel README) -- l'ottimizzatore non lo sa e non se ne cura,
  perche' massimizza `mf`, non `mu`. Massimizzare `mf` e massimizzare `mu` sono pero'
  equivalenti, perche' `mu` e' una funzione **affine crescente** di `mf`.
- Riga 13: `Q_vec = linspace(1.8, 7, 80)`. Il limite inferiore non e' arbitrario:
  a `t = 0` si ha `T/W = c*Q/m0 = c*Q`, quindi il decollo richiede `c*Q > 1`, cioe'
  `Q > 1/c = 1.667`. `Q = 1.8` da' `T/W = 1.08`, appena sopra soglia. Il decollo e' il
  caso vincolante perche' `T/W = c*Q/m` **cresce** man mano che la massa scende.
- Riga 15: `opts_ode = odeset('RelTol', 1e-10, 'AbsTol', 1e-12)`. Tolleranze
  estremamente strette. Il perche' e' il punto piu' importante dell'intero script --
  vedi il box qui sotto.
- Righe 16-18: opzioni `fsolve` esposte **esplicitamente** (`Display`,
  `MaxIterations`, `MaxFunctionEvaluations`, `FunctionTolerance`, `StepTolerance`)
  per riproducibilita'. `FunctionTolerance = StepTolerance = 1e-10`.
  Nota: **l'algoritmo non e' specificato**, quindi `fsolve` usa il default. Con 4
  residui e 4 incognite il sistema e' quadrato, e il default e' `trust-region-dogleg`,
  con **Jacobiano alle differenze finite** (nessun Jacobiano analitico e' fornito).

> **Possibile domanda d'esame** -- Perche' `RelTol = 1e-10` e `AbsTol = 1e-12`? Non
> sono esagerate per una simulazione di ascesa?
> *Risposta:* Non sono le tolleranze di una *simulazione*, sono le tolleranze di un
> *residuo dentro un Newton*. `fsolve` deve annullare il residuo a `1e-10`
> (`FunctionTolerance`) e costruisce il Jacobiano per **differenze finite**,
> perturbando le incognite di circa `sqrt(eps) ~ 1.5e-8` in relativo. Se il residuo
> fosse affetto da rumore di integrazione dell'ordine di `1e-6` (tolleranza ode45 di
> default), la derivata alle differenze finite sarebbe **rumore puro**: numeratore
> `1e-6`, denominatore `1e-8`, rapporto senza significato. La regola pratica e': la
> tolleranza dell'integratore deve essere **piu' stretta** della tolleranza sul
> residuo, non solo confrontabile. In piu' lo shooting sui costati e' intrinsecamente
> mal condizionato (piccole variazioni di `lambda_0` vengono amplificate
> esponenzialmente dall'integrazione in avanti nel residuo terminale), quindi il
> numero di condizionamento del Jacobiano e' grande e amplifica ulteriormente il
> rumore. Le tolleranze strette sono il prezzo del metodo indiretto.

---

## `%% Storage` (righe 20-24)

- Righe 23-24: `mf_results = nan(nQ, nyf)` e `sol_results = cell(nQ, nyf)`.
  L'inizializzazione a **NaN** (non a zero) e' deliberata: i punti in cui `fsolve`
  non converge restano NaN e vengono poi filtrati nei plot con la maschera
  `valid = ~isnan(...)` (righe 112, 147). Con zeri, i fallimenti apparirebbero come
  punti a massa finale nulla -- un artefatto grafico indistinguibile da un risultato.
- `sol_results` e' un cell array perche' ogni elemento e' il vettore delle 4 incognite
  convergenti `[lam_vx0; lam_vy0; lam_y; tf]`, che serve poi (i) come warm start per
  il `Q` successivo e (ii) per la ri-integrazione nei plot 1b e 1c.

---

## `%% Solve BVP for each yf, sweeping Q with continuation` (righe 26-94)

E' il blocco che fa il lavoro. La struttura e' un **doppio loop con continuazione a
due livelli**.

### Punto di partenza (righe 27-63)

```matlab
[~, idx0] = min(abs(Q_vec - 3));      % start near Q = 3
if jj == 1
    z_guess = [0.6; 3.8; 14; 0.30];
else
    z_guess = sol_results{idx0, jj-1};   % warm start dal yf precedente
end
```

- Righe 32-33: la struct `p` viene (ri)costruita a ogni `jj` con `c` e `yf`. `p` e'
  una variabile **di script**, non locale a una funzione: viene mutata piu' volte e
  ricostruita anche nei blocchi dei plot (riga 131, riga 162). Funziona perche' le
  closure `@(z) shooting1(z, p, opts_ode)` catturano `p` **per valore** al momento
  della creazione dell'handle, ma e' un pattern fragile: basta dimenticare di
  reimpostare un campo e si risolve silenziosamente il problema sbagliato.
- Riga 36: si parte da `Q` piu' vicina a 3. Sulla griglia `linspace(1.8, 7, 80)` il
  passo e' `5.2/79 = 0.0658`, quindi il nodo piu' vicino a 3 e' `Q ~ 2.985`
  (indice 19). Perche' proprio 3? Perche' `T/W = c*Q = 0.6*3 = 1.8` e' un rapporto
  spinta-peso "amichevole": abbastanza sopra 1 da non essere marginale, abbastanza
  basso da non produrre traiettorie estreme. E' il punto in cui il BVP e' meglio
  condizionato, quindi e' da li' che si parte a fare continuazione.
- Riga 39-40: **il guess iniziale a freddo** `[0.6; 3.8; 14; 0.30]`. E' l'unico numero
  "magico" dello script -- e va difeso, perche' i costati **non hanno significato
  fisico diretto** e non si possono indovinare da considerazioni ingegneristiche
  ovvie. Cosa si puo' dire di questo guess:
  - `phi(0) = atan2(3.8, 0.6) ~ 81 deg`: spinta quasi verticale al decollo. Corretto.
  - `lam_vy(tf) = 3.8 - 14*0.30 = -0.4`, quindi `phi(tf) = atan2(-0.4, 0.6) ~ -34 deg`:
    la spinta finisce **sotto** l'orizzonte. Non e' assurdo: verso fine burn il razzo
    ha ancora `vy > 0` e deve azzerarla per soddisfare `vy(tf) = 0`, quindi deve
    spingere leggermente verso il basso.
  - `mf = 1 - Q*tf = 1 - 2.985*0.30 ~ 0.105`, dello stesso ordine dei valori finali
    riportati (`mf* ~ 0.118` a `yf = 0.04`).
  - Il quarto residuo con questi numeri vale
    `H0 = -3.8 + T*(|lam_v0| - 1/c)` con `T = c*Q ~ 1.79` e `|lam_v0| = sqrt(0.6^2 + 3.8^2) = 3.847`,
    cioe' `H0 ~ -3.8 + 1.79*(3.847 - 1.667) ~ +0.10` -- piccolo rispetto ai termini
    individuali (~3.9). Il guess e' gia' quasi ammissibile sulla condizione di
    Hamiltoniana.

  (Questi valori li ho ricavati a mano dai numeri nel codice; lo script non li stampa.)
- Riga 43: per `jj >= 2` (`yf = 0.05, 0.06`) il guess non e' piu' a freddo: si prende
  la soluzione convergente **allo stesso `Q`** ma alla `yf` precedente. E' un secondo
  livello di continuazione, questa volta nella quota target.
- Righe 50-52: si costruisce `p.T = c * p.Q` (la spinta e' `T = c*Q`, coerente con
  `dm/dt = -Q` e `dm/dt = -T/c`) e si lancia il primo `fsolve`.
- Riga 54: `if ef > 0` -- `ef` e' l'**exit flag**. `fsolve` restituisce `ef > 0` solo se
  ha convergito. Se fallisce sul punto di partenza (riga 61) l'intero `yf` viene
  saltato con `continue`: senza il punto di ancoraggio non c'e' continuazione possibile.
- Righe 55-57: **ri-integrazione**. Una volta convergito, si integra di nuovo con
  `ode_burn` per estrarre la massa finale `Z(end,5)`. Nota: dato che `dm/dt = -Q` e'
  costante, si avrebbe `mf = 1 - Q*tf` in forma chiusa. La ri-integrazione e' quindi
  ridondante per la sola massa -- ma serve comunque il `Z` completo, ed e' un check
  implicito di consistenza dell'integratore.
- Riga 56: le condizioni iniziali `[0;0;0;0;1;1]` -- razzo fermo nell'origine, `m0 = 1`,
  e **`lam_m0 = 1`** (la normalizzazione dei costati, sesta componente).

### Sweep forward e backward (righe 65-93)

```matlab
z_prev = z_sol;
for ii = idx0+1:nQ            % Q crescente
    ...
    [z_sol, ~, ef] = fsolve(@(z) shooting1(z,p,opts_ode), z_prev, opts_fs);
    if ef > 0
        ...
        z_prev = z_sol;       % <- aggiornato SOLO se converge
    end
end
```

- Righe 67-78 e 82-93: la **continuazione (omotopia)**, il vero segreto per far
  convergere questo shooting. Non si risolve ogni `Q` da zero: si risolve `Q_{k+1}`
  partendo dalla soluzione convergente di `Q_k`. Poiche' i due problemi differiscono
  di poco (`Delta Q ~ 0.066`), il guess cade dentro il basin di attrazione di Newton.
  Un cold start su `Q = 7` quasi certamente divergerebbe.
- Si sweepa in **entrambe le direzioni** dal punto di ancoraggio: in avanti verso
  `Q = 7`, all'indietro verso `Q = 1.8`. La riga 81 resetta `z_prev` alla soluzione
  di ancoraggio prima di partire all'indietro -- se non lo facesse, ripartirebbe dalla
  soluzione a `Q = 7`, lontanissima.
- Righe 76 e 91: `z_prev = z_sol` sta **dentro** `if ef > 0`. Questo e' un dettaglio
  di robustezza importante: se un `Q` fallisce, la catena **non viene avvelenata** --
  il `Q` successivo riparte comunque dall'ultimo guess *buono*, non da un vettore
  divergente. Costa un buco (NaN) nella curva, non un crollo dell'intero sweep.
- **Nota di onesta':** a differenza del punto di partenza (riga 61), i fallimenti
  dentro gli sweep sono **completamente silenziosi**: nessun `fprintf`, nessun
  warning. Se una porzione della curva `mf(Q)` fosse mancante, lo si scoprirebbe solo
  guardando il grafico. E' un limite reale dello script.

> **Possibile domanda d'esame** -- Cos'e' la continuazione e perche' e' *necessaria*,
> non solo comoda, nel metodo indiretto?
> *Risposta:* La continuazione (o omotopia) consiste nell'immergere il problema in una
> famiglia parametrizzata (qui da `Q`) e risolvere una catena di problemi vicini,
> usando ogni soluzione come guess iniziale del successivo. E' necessaria perche' nel
> metodo indiretto le incognite sono i costati, che non hanno interpretazione fisica:
> non esiste un modo ragionevole di indovinarli. Inoltre l'integrazione in avanti
> amplifica esponenzialmente le perturbazioni sui costati iniziali, quindi il basin di
> convergenza di Newton e' piccolissimo. La continuazione garantisce che ogni solve
> parta gia' **dentro** il basin. E' la ragione per cui lo script risolve prima a
> `Q ~ 3` (dove il problema e' ben condizionato) e poi si muove a passi di 0.066.

---

## `%% Optimal Q* per target altitude` (righe 96-105)

- Righe 101-102: `[mf_star(jj), idx] = max(mf_results(:, jj))` -- l'ottimo su `Q` non
  e' trovato con un ottimizzatore, ma per **enumerazione sulla griglia**. La
  risoluzione su `Q*` e' quindi limitata dal passo della griglia, `~0.066`. E'
  accettabile perche' la curva `mf(Q)` e' piatta vicino al massimo (il report lo dice
  esplicitamente), ma va detto: `Q*` riportato ha 4 cifre decimali (riga 103) che
  **non sono giustificate** dalla risoluzione della griglia.
- Riga 104: il payload `mf*(1+eta) - eta`, calcolato a posteriori (vedi sopra).
  Massimizzare `mf` e' equivalente a massimizzare il payload perche' la mappa e'
  affine crescente in `mf` -- quindi lo stesso `Q*` ottimizza entrambi.

---

## `%% PLOT 1a: Final mass vs Q` (righe 107-118)

- Riga 112: `valid = ~isnan(mf_results(:,jj))` -- la maschera che esclude i `Q` non
  convergenti. Buca la curva invece di falsificarla.
- Il risultato atteso (e riportato nel README) e' un **massimo interno** in `Q` per
  ciascuna `yf`, che si sposta verso `Q` piu' basse al crescere della quota target.

---

## `%% PLOT 1b: Velocity losses vs Q for yf = 0.04` (righe 120-153)

E' il blocco piu' "fisico" dello script: decompone il delta-V perso.

```matlab
ic8 = [0; 0; 0; 0; 1; 1; 0; 0];
[~, Z8] = ode45(@(t,z) ode_burn_losses(t,z,pp), [0 tf], ic8, opts_ode);
Wd_vec(ii) = Z8(end,7);
Wg_vec(ii) = Z8(end,8);
Wt_vec(ii) = c * log(1/mf_ii) - 1;   % total = DV_ideal - V_final
```

**Derivazione dell'identita' di decomposizione.** Il delta-V ideale (Tsiolkovsky) e'

    DV_ideal = c * ln(m0/mf) = c * ln(1/mf) = integrale di (T/m) dt

(l'ultima uguaglianza segue da `dm/dt = -T/c`). La velocita' effettivamente ottenuta
e' pero' `v_f = 1` (nondim). La differenza e' la perdita totale. Proiettando le
equazioni del moto lungo il versore di velocita', con `psi = atan2(vy, vx)` angolo di
traiettoria e `phi` angolo di spinta:

    d|v|/dt = (T/m)*cos(phi - psi) - sin(psi)

Integrando da fermo (`|v|(0) = 0`) fino a `|v|(tf) = v_f` e sottraendo da `DV_ideal`:

    DV_ideal - v_f  =  integrale di (T/m)*[1 - cos(phi - psi)] dt   (= W_d)
                     + integrale di sin(psi) dt                      (= W_g)

- `W_d` (**misalignment loss**, riga 321 in `ode_burn_losses`): la frazione di spinta
  che non spinge lungo la velocita'. Si annulla **solo se** `phi = psi` identicamente,
  cioe' se la spinta e' sempre allineata alla velocita' -- cosa che l'ottimo *non* fa,
  perche' deve anche curvare la traiettoria.
- `W_g` (**gravity loss**, riga 322): `integrale di sin(psi) dt`. E' la componente di
  gravita' opposta al moto. Massima quando si vola verticale (`psi = 90 deg`,
  `sin(psi) = 1`, si perde 1 g per ogni secondo), nulla quando si vola orizzontale.
  E' per questo che una `Q` bassa e' penalizzante: allunga il tempo di volo mentre si
  e' ancora ripidi.
- Riga 142: `Wt_vec(ii) = c*log(1/mf_ii) - 1`. Il `-1` e' `v_f = 1`, hard-coded.
  E' legittimo perche' `vx(tf) = 1, vy(tf) = 0` sono **condizioni al contorno imposte**
  dallo shooting; ma se una soluzione fosse mal convergita, `Wt` verrebbe calcolata
  come se `v_f` fosse esattamente 1 e la discrepanza si nasconderebbe.
  **Diagnostico gratis:** nel grafico, `Wd + Wg` deve coincidere con `Wt`. Se le curve
  non si sovrappongono, o l'integrazione o la convergenza sono da rifare. Le tre curve
  sono calcolate per **strade diverse** (le prime due per quadratura lungo la
  traiettoria, la terza da una formula chiusa sulla massa), quindi il loro accordo e'
  una verifica non banale.

> **Possibile domanda d'esame** -- Perche' `mf(Q)` ha un massimo interno? Spiegalo con
> `Wd` e `Wg`.
> *Risposta:* La massa finale e' `mf = exp(-DV_ideal/c)` e `DV_ideal = v_f + Wd + Wg`,
> quindi massimizzare `mf` equivale a **minimizzare la perdita totale** `Wd + Wg`. A
> `Q` bassa la spinta e' appena sopra il peso, il burn dura a lungo e il razzo resta
> ripido a lungo: `Wg = integrale di sin(psi) dt` esplode. A `Q` alta il burn e'
> brevissimo, `Wg` crolla, ma la spinta enorme deve essere fortemente disallineata
> dalla velocita' per curvare la traiettoria nel poco tempo disponibile, e `Wd`
> cresce. I due contributi vanno in direzioni opposte, quindi la somma ha un minimo
> interno: `Q*`.

---

## `%% PLOT 1c: Trajectory and angles for optimal Q` (righe 155-205)

- Righe 157-158: si ricalcola `[mf_opt, idx_opt] = max(mf_results(:, jj_ref))`,
  duplicando quello che le righe 101-102 avevano gia' fatto. Ridondanza innocua ma
  reale.
- Riga 169: integrazione **densa** su `linspace(0, tf, 500)`. Nota: passare un vettore
  di tempi a `ode45` **non cambia i passi di integrazione** (che restano adattivi e
  guidati da `RelTol`/`AbsTol`): cambia solo dove l'output viene interpolato
  (interpolazione densa di ordine 4). Quindi l'accuratezza non degrada, si ottengono
  solo curve lisce da plottare.
- Righe 180-189: ricostruzione degli angoli **a posteriori**, senza reintegrare i
  costati:

      lam_vy_k = z_opt(2) - z_opt(3) * T_sol(kk);
      phi_traj(kk) = atan2(lam_vy_k, z_opt(1));

  cioe' `phi(t) = atan2(lam_vy0 - lam_y*t, lam_vx0)`. E' esattamente la stessa formula
  che `ode_burn` usa internamente (righe 20 e 24 di `ode_burn.m`) -- e' la legge
  **linear-tangent**. Il fatto di poterla ricostruire fuori dall'integratore, con tre
  soli numeri, e' la manifestazione pratica del fatto che i costati di velocita' hanno
  soluzione in forma chiusa.
- Righe 184-188: guard su `V_k > 1e-10`. A `t = 0` la velocita' e' **esattamente zero**
  e `psi = atan2(0,0)` non e' definito (in MATLAB darebbe 0, che sarebbe fuorviante:
  suggerirebbe volo orizzontale a un razzo fermo sulla rampa). Il fallback e'
  `psi = phi`, che e' la scelta fisicamente sensata: nell'istante in cui il veicolo si
  stacca, la velocita' nasce nella direzione dell'accelerazione netta, che e'
  dominata dalla spinta. **Un guard analogo esiste anche in `ode_burn_losses`**
  (righe 307-312, ma con soglia `1e-12` invece di `1e-10`), ed e' li' che davvero
  conta: senza, `Wd` avrebbe un integrando indefinito all'istante iniziale.
- Righe 191-203: le due figure (traiettoria `x-y` con `axis equal`, e i due angoli
  `phi`/`psi` in gradi vs tempo). Il comportamento atteso di `phi` e' la discesa
  monotona da quasi verticale a quasi orizzontale (ed eventualmente sotto lo zero):
  e' la legge linear-tangent, `tan(phi)` lineare in `t`.

---

## `%% EXPORT FIGURES` (righe 207-223)

- Riga 208: `fig_dir` costruito con `fileparts(mfilename('fullpath'))` -- quindi le
  figure finiscono in `HM1/figures/` **indipendentemente dalla working directory**.
- Riga 211: `findobj(groot, 'Type', 'figure')` prende **tutte** le figure aperte. Vedi
  la nota sul `close all` di riga 6: e' quello che garantisce che siano solo le proprie.
- Righe 215-220: `theme(fig, 'light')` dentro un `try/catch`, con fallback
  `fig.Color = 'w'` per MATLAB pre-R2025a. Serve a forzare figure su fondo chiaro anche
  se il desktop MATLAB e' in dark mode -- altrimenti le PNG esportate finirebbero nel
  report con sfondo scuro.
- Riga 210 + 222: `slugify` produce nomi tipo `task1_task_1a_final_mass_vs_q.png`
  (il "task1_" del prefisso piu' il "Task 1a" slugificato dal `Name` della figura --
  da cui la ripetizione). Sono i nomi che il report LaTeX si aspetta.
- **Deviazione dalle convenzioni della repo:** il `CLAUDE.md` chiede SVG (vettoriale) +
  PNG (preview). Qui si esporta **solo PNG a 200 dpi**. Nel report va bene lo stesso,
  ma vale la pena saperlo.

---

## `shooting1` -- il residuo del BVP (righe 227-274)

E' la funzione **piu' importante dello script**: e' la definizione operativa del
Boundary Value Problem.

```matlab
function res = shooting1(z0, p, opts_ode)
    lam_vx0 = z0(1);  lam_vy0 = z0(2);
    lam_y   = z0(3);  tf      = z0(4);
    ...
    ic = [0; 0; 0; 0; 1; 1];   % state + lam_m0 = 1 (normalization)
    [~, Z] = ode45(@(t,z) ode_burn(t, z, pp), [0 tf], ic, opts_ode);
    zf = Z(end,:);

    Lam0 = sqrt(lam_vx0^2 + lam_vy0^2);
    H0   = -lam_vy0 + p.T * (Lam0 - 1/p.c);

    res = [zf(2) - p.yf;    % y(tf)  = yf
           zf(3) - 1;       % vx(tf) = 1
           zf(4);           % vy(tf) = 0
           H0];             % H(0)   = 0  (free final time)
end
```

- Riga 227: firma. `z0` sono le **4 incognite dello shooting**
  `[lam_vx0; lam_vy0; lam_y; tf]`; `p` porta `c, T, Q, yf`. Chiamata solo da `fsolve`.
- Riga 238: di nuovo **nessun blocco `arguments`, per scelta** -- `shooting1` gira
  dentro il loop di `fsolve` (decine di valutazioni per iterazione, per costruire il
  Jacobiano alle differenze finite, moltiplicato per 240 solve).
- Righe 245-248: **guard sul tempo finale.** Se `tf <= 0` o `tf > 2`, si restituisce
  un residuo enorme `1e6*ones(4,1)`.
  **Nota di onesta': e' un hack, non una barriera liscia.** Il residuo diventa
  discontinuo sul bordo della regione ammissibile; se `fsolve` ci finisce dentro con
  una perturbazione di differenze finite, la colonna corrispondente del Jacobiano
  diventa spazzatura (`~1e6/1e-8`). In pratica non da' problemi perche' la
  continuazione tiene `tf` ben all'interno (dell'ordine di 0.1-0.35 nondim), ma e'
  fragile. Nota anche che il limite superiore `tf > 2` (cioe' `~1590 s` dimensionali)
  e' **arbitrario**: una soluzione con `tf > 2` verrebbe rifiutata a priori.
- Righe 250-253: si copiano le tre incognite dei costati dentro `pp`, che diventa la
  struct passata a `ode_burn`.
- Riga 255: le condizioni iniziali. Le prime cinque sono la fisica (`x = y = vx = vy = 0`,
  `m0 = 1`). La sesta e' `lam_m0 = 1`: **la normalizzazione dei costati**.
- Righe 257-263: `try/catch` attorno a `ode45` con lo stesso residuo-penalita' `1e6`.
  Protegge dai casi in cui l'integratore esplode (per esempio `m -> 0` se `Q*tf > 1`).
  Stessa critica di prima: e' una discontinuita', non una barriera.
- Righe 265-268: **la condizione di Hamiltoniana nulla, imposta a `t0`.**
  Derivazione: partendo dall'Hamiltoniana ottima

      H* = lam_y*vy + (T/m)*|lam_v| - lam_vy - lam_m*Q

  e valutandola a `t = 0` con `vx0 = vy0 = 0`, `m0 = 1`, `lam_m0 = 1` e `Q = T/c`:

      il termine lam_y*vy  -> 0        (perche' vy0 = 0)
      il termine (T/m)*|lam_v| -> T*|lam_v0|   (perche' m0 = 1)
      il termine -lam_vy   -> -lam_vy0
      il termine -lam_m*Q  -> -T/c     (perche' lam_m0 = 1 e Q = T/c)

  da cui

      H(0) = -lam_vy0 + T*|lam_v0| - T/c
           = -lam_vy0 + T*( |lam_v0| - 1/c ) = 0

  che e' **esattamente** la riga 268. Notare che la parentesi `(|lam_v0| - 1/c)` e' il
  valore iniziale della **funzione di switching** `S = |lam_v|/m - lam_m/c`, valutata
  con `m0 = 1` e `lam_m0 = 1`. `S(0) > 0` significa "motore acceso" -- coerente con un
  arco propulso che parte da terra. (`S` diventa una quantita' di primo piano nel
  Task 3, dove il coast inizia quando `S = 0`.)
- **Perche' `H = 0` a `t0` e non a `tf`?** La condizione di trasversalita' per tempo
  finale libero e' `H(tf) = 0`. Ma il sistema e' **autonomo**, quindi `H` e' costante
  lungo la traiettoria ottima, e imporre `H(0) = 0` e' equivalente. Il vantaggio e'
  enorme: a `t0` la formula e' **puramente algebrica** e non contiene errore di
  integrazione, mentre `H(tf)` andrebbe valutata sullo stato integrato (e su `lam_m(tf)`,
  che ha accumulato errore per tutto l'arco).
- Righe 270-273: i **quattro residui**:
  - `zf(2) - p.yf` -> vincolo di quota `y(tf) = yf`
  - `zf(3) - 1`    -> vincolo di velocita' orizzontale `vx(tf) = 1` (velocita' orbitale)
  - `zf(4)`        -> vincolo `vy(tf) = 0` (iniezione orizzontale, apogeo)
  - `H0`           -> condizione di tempo finale libero

  **Sistema quadrato: 4 incognite, 4 residui.** Notare cosa *non* c'e': non c'e' nessun
  vincolo su `x(tf)` (la posizione orizzontale e' libera: `lam_x = 0` per trasversalita'),
  e non c'e' nessun residuo su `lam_m(tf) = 1` (assorbito dalla normalizzazione
  `lam_m0 = 1`). Le due "migliorie" numeriche -- normalizzazione dei costati e `H`
  imposta a `t0` -- tolgono ciascuna **una incognita e una condizione**, portando il
  BVP da 6x6 a 4x4. E `lam_m` **non compare mai nel residuo**: `zf(6)` non e' letto.

> **Possibile domanda d'esame** -- Dove finiscono le condizioni di trasversalita'
> `lam_x(tf) = 0` e `lam_m(tf) = 1`? Non compaiono nel residuo.
> *Risposta:* Entrambe sono state usate *analiticamente* prima di scrivere il codice,
> non numericamente. `lam_x(tf) = 0` (perche' `x(tf)` e' libero e non entra nel costo),
> combinata con `lam_x_dot = 0`, da' `lam_x = 0` **su tutto l'arco**: sparisce
> dall'Hamiltoniana e da `ode_burn`. `lam_m(tf) = 1` e' invece **sostituita** dalla
> normalizzazione `lam_m0 = 1`: poiche' l'Hamiltoniana e' omogenea di grado 1 in
> `lambda`, i costati sono determinati solo a meno di un fattore di scala positivo, e
> la trasversalita' su `lam_m` serviva solo a fissarlo. Fissarlo all'istante iniziale
> invece che finale e' equivalente e piu' comodo (rimuove una incognita **e** una
> condizione). Conseguenza da sapere: `lam_m(tf)` calcolato dal codice **non vale 1**.

---

## `set_costates` (righe 276-287)

- Utility di tre righe: copia `z_sol(1:3)` nei campi `lam_vx0`, `lam_vy0`, `lam_y`
  della struct. Esiste solo per evitare di ripetere le tre assegnazioni nei cinque
  punti dello script che ri-integrano una soluzione convergente (righe 55, 72, 87, 133,
  164). Nota che **ignora `z_sol(4)` (`tf`)**, che va letto separatamente dal chiamante
  (righe 134, 165).

---

## `ode_burn_losses` (righe 289-323)

- E' `ode_burn` con **due integratori di perdita accodati**. Le prime sei componenti
  (righe 314-320) sono una **copia letterale** di `ode_burn.m` (righe 27-33).
- Riga 321: `dz(7) = (T/m)*(1 - cos(phi - psi))` -> `Wd`, la misalignment loss.
  Nota che `1 - cos(.) >= 0` sempre, quindi `Wd` e' **monotona crescente**: e' una
  perdita, non puo' recuperare. E' nulla solo dove `phi = psi`.
- Riga 322: `dz(8) = sin(psi)` -> `Wg`, la gravity loss. Non e' monotona: se il razzo
  scendesse (`psi < 0`), `sin(psi) < 0` e la "perdita" diventerebbe un guadagno.
  Nell'ascesa `psi > 0` quasi ovunque, quindi in pratica cresce.
  Notare che **`Wg` non dipende dalla spinta ne' dalla massa**: e' pura geometria della
  traiettoria -- quanto tempo si passa a volare ripidi.
- Righe 307-312: il guard `V > 1e-12` con fallback `psi = phi` (stessa struttura del
  guard di riga 184, ma con soglia `1e-12` invece di `1e-10`: le due soglie sono
  scelte indipendenti, non una costante condivisa). Qui e' **necessario**, non
  cosmetico: senza, `psi` sarebbe indefinito a `t = 0` e l'integrando di `Wd`
  inizierebbe da un valore arbitrario.
- **Nota di onesta' -- duplicazione del codice.** `ode_burn_losses` replica il RHS di
  `ode_burn` invece di chiamarlo. Se un giorno si cambiasse `ode_burn.m` (per esempio
  aggiungendo il drag), questa copia **non** si aggiornerebbe e i due divergerebbero
  silenziosamente. La ragione della scelta e' presumibilmente la performance (evitare
  una chiamata di funzione e la ricostruzione del vettore nel hot loop), ma il costo e'
  un rischio di manutenzione reale, e va dichiarato.

---

## Possibili domande d'esame

**D: Elenca le quattro incognite dello shooting e i quattro residui, e spiega perche'
il sistema e' quadrato.**
R: Incognite: `lam_vx0`, `lam_vy0`, `lam_y`, `tf`. Residui: `y(tf) - yf`, `vx(tf) - 1`,
`vy(tf)`, `H(0)`. Il sistema di partenza avrebbe 5 costati iniziali + `tf` = 6
incognite, e 6 condizioni (3 vincoli terminali di stato + `lam_x(tf) = 0` +
`lam_m(tf) = 1` + `H(tf) = 0`). Due riduzioni analitiche lo portano a 4x4:
(i) `lam_x = 0` identicamente (per trasversalita' + `lam_x_dot = 0`), che elimina
un'incognita e la sua condizione; (ii) la normalizzazione `lam_m0 = 1`, lecita perche'
`H` e' omogenea di grado 1 in `lambda`, che elimina `lam_m0` come incognita e rende
superflua la trasversalita' `lam_m(tf) = 1`. Restano 4 e 4, e `lam_m` non entra mai nel
residuo.

**D: Perche' `H = 0` viene imposta a `t = 0` e non a `t = tf`, dove la condizione di
trasversalita' la richiede?**
R: Perche' il sistema e' autonomo (l'Hamiltoniana non dipende esplicitamente dal tempo),
quindi `H` e' un integrale primo: `H(t) = cost`. Se `H(tf) = 0` allora `H(0) = 0`, e
viceversa: le due condizioni sono **equivalenti**. Ma a `t = 0` lo stato e' noto
esattamente (`vx = vy = 0`, `m = 1`, `lam_m = 1`), quindi `H(0)` si scrive in forma
**algebrica chiusa**, `H(0) = -lam_vy0 + T*(|lam_v0| - 1/c)`, senza errore di
integrazione. Valutarla a `tf` significherebbe usare lo stato integrato -- inquinato
dall'errore numerico accumulato -- e anche `lam_m(tf)`, che e' il costato con l'errore
maggiore. Imporre la condizione dove e' esatta migliora il condizionamento del residuo.

**D: Perche' `RelTol = 1e-10 / AbsTol = 1e-12` e non le tolleranze di default di ode45?**
R: Perche' l'output di `ode45` **e' il residuo** che `fsolve` deve annullare a `1e-10`,
e `fsolve` costruisce il Jacobiano per differenze finite perturbando le incognite di
`~sqrt(eps) ~ 1.5e-8`. Se il residuo avesse rumore di integrazione `~1e-6` (default
`RelTol = 1e-3`), la derivata numerica sarebbe `1e-6/1e-8 = 100` volte il rumore: il
Jacobiano sarebbe puro rumore e Newton non convergerebbe (o convergerebbe a una
soluzione fasulla). La regola: la tolleranza dell'integratore deve stare **sotto** la
tolleranza sul residuo, con margine. Si aggiunge il fatto che lo shooting sui costati e'
mal condizionato per natura, quindi il Jacobiano amplifica ulteriormente ogni rumore.

**D: Che cos'e' la strategia di continuazione qui, e cosa succederebbe senza?**
R: Si risolve prima il BVP a `Q ~ 3` (`T/W = 1.8`, punto ben condizionato) partendo dal
guess a freddo `[0.6; 3.8; 14; 0.30]`, poi si sweepa `Q` **in avanti e all'indietro**
usando ogni soluzione convergente come guess della successiva (passo `~0.066`). Senza
continuazione bisognerebbe indovinare i costati a freddo per ognuno degli 80 `Q` e 3
`yf` -- e i costati non hanno significato fisico, quindi non si indovinano. Il basin di
convergenza di Newton per uno shooting indiretto e' molto piccolo (l'integrazione in
avanti amplifica esponenzialmente le perturbazioni su `lambda_0`), quindi la maggior
parte dei cold start divergerebbe. Il codice fa continuazione anche su `yf`: la prima
soluzione per `yf = 0.05` e' presa da quella convergente a `yf = 0.04` alla stessa `Q`.

**D: Come si decompone la perdita di velocita' e perche' `Wd` e `Wg` vanno in direzioni
opposte con `Q`?**
R: `DV_ideal = c*ln(1/mf) = v_f + Wd + Wg`, con
`Wd = int (T/m)*[1 - cos(phi - psi)] dt` (misalignment) e `Wg = int sin(psi) dt`
(gravita'), ottenute proiettando le equazioni del moto sul versore di velocita'.
Al crescere di `Q` il burn si accorcia (perche' `T/W` cresce), quindi `Wg`, che e'
un integrale nel tempo di una quantita' limitata da 1, **diminuisce**. Ma per curvare
la traiettoria dalla verticale all'orizzontale nel poco tempo disponibile, la spinta
deve essere fortemente disallineata dalla velocita', quindi `1 - cos(phi - psi)` cresce
e `Wd` **aumenta**. Il minimo della somma e' `Q*`, che e' anche il massimo di `mf`
(perche' `mf = exp(-DV_ideal/c)` e' decrescente nella perdita totale).

**D: Quali sono i punti deboli / gli hack di questo script?**
R: (i) I guard `tf <= 0 || tf > 2` e il `try/catch` attorno a `ode45` restituiscono un
residuo costante `1e6`: sono **barriere discontinue**, non lisce, e possono corrompere
il Jacobiano alle differenze finite se `fsolve` ci finisce dentro. In pratica la
continuazione le evita. (ii) Il limite `tf > 2` e' arbitrario. (iii) I fallimenti di
`fsolve` dentro gli sweep sono **silenziosi** (lasciano solo un NaN). (iv) `Q*` e'
trovato per enumerazione su griglia (passo `0.066`), quindi le 4 cifre decimali stampate
non sono giustificate. (v) `ode_burn_losses` **duplica** il RHS di `ode_burn` invece di
chiamarlo: rischio di divergenza in manutenzione. (vi) La massa finale si potrebbe
ottenere in forma chiusa (`mf = 1 - Q*tf`, perche' `Q` e' costante) invece di
ri-integrare. (vii) Non c'e' nessun guard su `m > 0`. (viii) Le figure sono esportate
solo in PNG, mentre il `CLAUDE.md` chiede anche SVG.
