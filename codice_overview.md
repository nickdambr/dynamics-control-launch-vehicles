# Documentazione del codice — HM1, HM2 & HM3

Schema delle funzioni MATLAB scritte per HM1 (ottimizzazione indiretta
dell'ascesa), HM2 (powered descent via collocazione diretta + SCvx) e HM3
(controllo d'assetto LV a max-q in frequenza). Ogni
sezione segue lo stesso schema:

1. struttura del file (cosa fa lo script al top-level),
2. tabella delle funzioni locali con firma, ruolo e — quando rilevante —
   formato dell'input/output.

Tutti i collegamenti puntano al sorgente (file:linea).

---

## HM1 — Ottimizzazione indiretta dell'ascesa

Pattern comune a tutti i task:

```
script principale (main_taskN.m)
   ├── definisce parametri (c, Q, η, y_f, …) e tolleranze
   ├── eventuale fase verticale ad eventi (Task 2/3)
   ├── fsolve(shooting_residual, z0)   ← single-shooting BVP
   ├── re-integrazione fine della soluzione ottima
   └── plot + export figure
```

I 5 stati nondimensionali sono `[x, y, vx, vy, m]`; il sesto stato aggiunto
in ODE è il costato `λ_m` (gli altri costati hanno forma chiusa lineare).
Le condizioni terminali sono fisse: `vx(t_f)=1, vy(t_f)=0, y(t_f)=y_f`.
**Migliorie numeriche** (note del corso): i costati sono normalizzati fissando
`λ_m0=1`, e la condizione di tempo libero `H=0` è imposta a `t₀` in forma
algebrica. Lo shooting ha così **4 unknown** `[λ_vx0, λ_vy0, λ_y, t_f]` e `λ_m`
non entra nel residuo (Task 1/2/4); Task 3 mantiene `λ_m0` per lo switch del coast.

### Funzione condivisa

| Funzione | Posizione | Ruolo |
|---|---|---|
| `ode_burn` | [HM1/ode_burn.m:1](HM1/ode_burn.m#L1) | RHS della fase di volo propulso: integra lo **stato esteso** `[x,y,vx,vy,m,λ_m]` con la **bilinear-tangent law** `φ = atan2(λ_vy, λ_vx)`. `λ_vx` costante, `λ_vy` lineare in `t`. Equazione di `λ_m`: `dλ_m/dt = (T/m²)·‖λ_v‖`. |

### `main_task1.m` — arco singolo, sweep su `Q`

Sweep di `Q ∈ [1.8, 7]` (sopra la soglia di liftoff `1/c≈1.67`) per tre quote
target `y_f ∈ {0.04, 0.05, 0.06}`, con **continuazione** (warm-start della
successiva risoluzione dalla precedente, seme a `Q≈3`) per garantire
convergenza dello shooting su tutto l'intervallo. Gli `Q*` ottimi risultano
`2.52 / 2.33 / 2.19`.
Dalla soluzione ottima vengono inoltre calcolate le perdite di gravità e
sterzo per decomposizione.

| Funzione | Riferimento | Ruolo |
|---|---|---|
| `shooting1` | [main_task1.m:227](HM1/main_task1.m#L227) | Residuo di shooting (4 componenti) sulle unknown `z = [λ_vx0, λ_vy0, λ_y, t_f]` (con `λ_m0≡1`). Integra `ode_burn` su `[0, t_f]` e impone: `y(t_f)=y_f`, `vx(t_f)=1`, `vy(t_f)=0`, `H(0)=−λ_vy0+T(‖λ_v0‖−1/c)=0` (tempo libero, imposto a `t₀` in forma algebrica). |
| `set_costates` | [main_task1.m:275](HM1/main_task1.m#L275) | Pack helper: copia `λ_vx0, λ_vy0, λ_y` dentro la struct `p` per passarli a `ode_burn`. |
| `ode_burn_losses` | [main_task1.m:283](HM1/main_task1.m#L283) | Versione estesa di `ode_burn` (8 stati) che integra in parallelo `dW_d/dt = (T/m)(1−cos(φ−ψ))` (perdita di sterzo) e `dW_g/dt = sin(ψ)` (perdita gravitazionale). Usata solo nel post-processing. |

### `main_task2.m` — salita verticale + burn

Sequenza vertical-climb (fino a `y₁ = 10⁻⁴`) seguita da arco di burn
ottimo. La fase verticale è integrata con event detection per agganciare
la quota target; lo shooting BVP riparte poi dallo stato finale del climb.

| Funzione | Riferimento | Ruolo |
|---|---|---|
| `ode_vertical` | [main_task2.m:179](HM1/main_task2.m#L179) | RHS della salita verticale: `φ=π/2`, stato ridotto `[y, vy, m]`. |
| `event_altitude` | [main_task2.m:185](HM1/main_task2.m#L185) | Evento `ode45` che termina l'integrazione quando `y = y_target`. |
| `shooting2` | [main_task2.m:191](HM1/main_task2.m#L191) | Stessa logica (migliorata) di `shooting1` ma con I.C. **non nulle** (uscita del climb) via `p.x0,…,p.m0`. Unknown: `[λ_vx0, λ_vy0, λ_y, t_burn]` (`λ_m0≡1`); `H=0` valutato sullo stato post-climb (`H0=λ_y·vy1+(T/m1)‖λ_v0‖−λ_vy0−T/c`). |

### `main_task3.m` — salita verticale + burn + coast

Aggiunge l'arco balistico finale dopo il cutoff: durante coast `λ_m` è
costante, `λ_vy` rampa lineare, quindi si **chiude analiticamente** la
condizione di iniezione con `y_c + ½ vy_c² = y_f` invece di integrare.

| Funzione | Riferimento | Ruolo |
|---|---|---|
| `ode_vertical` | [main_task3.m:200](HM1/main_task3.m#L200) | Identica a Task 2. |
| `event_altitude` | [main_task3.m:205](HM1/main_task3.m#L205) | Identica a Task 2. |
| `shooting3` | [main_task3.m:211](HM1/main_task3.m#L211) | Residuo 5-d alla fine del **burn** (Task 3 mantiene `λ_m0` per lo switch). Impone: `vx_c=1`, `y_c + ½ vy_c² = y_f` (match balistico), `λ_m(t_c)=1`, `λ_vy(t_c) = λ_y · vy_c` (ottimalità coast), `S(t_c) = ‖λ_v‖/m_c − λ_m/c = 0` (switching; `λ_m=1` al cutoff). Una guess a burn corto seleziona la radice fisica `vy_c>0`. |

### `main_task4.m` — staging

Due burn arcs separati da un evento di jettison della struttura del primo
stadio (`m_s¹ = η·Q·t_s`). Sweep esterno su `t_s`, BVP interno con un
warm-start sulla soluzione single-stage di Task 1.

| Funzione | Riferimento | Ruolo |
|---|---|---|
| `shooting_single` | [main_task4.m:215](HM1/main_task4.m#L215) | Identico a `shooting1` (4 unknown, `H0=0`); riferimento single-stage per warm-start. |
| `shooting_twostage` | [main_task4.m:241](HM1/main_task4.m#L241) | Integra `ode_burn` su `[0, t_s]`, applica `m⁺ = m⁻ − η·Q·t_s` (`λ_m` continua), integra di nuovo su `[t_s, t_f]` e impone gli stessi 4 residui di `shooting1` (`H0=0`). `t_s` è parametro fissato per ciascuna iterazione dello sweep. |

---

## HM2 — Powered descent

Il problema OCP è formulato in forma **non-dimensionale** con scale
`L_ref=y₀, t_ref=√(L_ref/g), V_ref=√(g·L_ref), m_ref=m₀, T_ref=m₀·g`.
L'unico parametro residuo è `V_c = V_ref/c`. Tutto il lavoro interno è in
non-dim e viene riportato in SI solo per stampa/plot.

### `main_task1.m` — collocazione trapezoidale (baseline)

Trascrizione trapezoidale su `N=50` nodi, vettore decisionale
`z = [x,y,vx,vy,m,Tx,Ty]` per nodo (lunghezza `7·N`). Solver: `fmincon`
con `'sqp'`. Sweep di sensibilità su `t_f ∈ {0.95, 1.00, 1.05}·38 s`.

| Funzione | Riferimento | Ruolo |
|---|---|---|
| `nondim` | [main_task1.m:79](HM2_powered_descent/main_task1.m#L79) | Calcola le scale di riferimento e ritorna `dnd` con tutti i dati problema convertiti in non-dim. |
| `dim_sol` | [main_task1.m:99](HM2_powered_descent/main_task1.m#L99) | Riporta una `sol` non-dim in SI (`t·t_ref, x·L_ref, …`) e aggiunge `Tmag` e fuel consumato. |
| `solve_trapcol` | [main_task1.m:115](HM2_powered_descent/main_task1.m#L115) | Assembla e risolve l'NLP: initial guess (interpolazione lineare + hover), box bounds, B.C. lineari, obiettivo `−z(m_N)`, vincoli non-lineari via `trap_nonlcon`. Ritorna la soluzione spacchettata. |
| `trap_nonlcon` | [main_task1.m:203](HM2_powered_descent/main_task1.m#L203) | Costruisce: (a) `5·(N−1)` defects trapezoidali `x_{k+1}−x_k − ½dt(f_k+f_{k+1})` come `c_eq`, (b) vincoli di path (bounds spinta + glide-slope `±x − tan(θ)·y ≤ 0`) come `c_ineq`. |
| `dyn_rhs` | [main_task1.m:230](HM2_powered_descent/main_task1.m#L230) | Dinamica continua non-dim. Stato `s=[x,y,vx,vy,m,Tx,Ty]`, output `ṡ=[vx, vy, Tx/m, Ty/m−1, −V_c·‖T‖]`. |
| `plot_results` | [main_task1.m:242](HM2_powered_descent/main_task1.m#L242) | Figure trajectory + thrust + mass + glide-slope, una curva per ciascun `t_f`. |

### `main_task2.m` — ZOH multiple-shooting + SCvx + variante YALMIP/ECOS

Script più articolato: risolve in **parallelo** quattro trascrizioni con
gli stessi dati, le stessa initial conditions e (dove rilevante) lo stesso
warm-start. Le quattro varianti:

1. **Trapezoidale** (riproduce il baseline di Task 1, usata come riferimento e initial-ref).
2. **ZOH non-lineare con RK4**: multiple shooting `x_{k+1} = RK4(x_k, u_k, dt)` con `n_sub` sub-step.
3. **LTV-ZOH + SCvx con fmincon**: linearizza attorno alla reference corrente via l'ODE aumentata di Appendice A, costruisce le matrici `(Ā, B̄, c̄)` per ogni intervallo e risolve l'NLP linearizzato con `fmincon`.
4. **LTV-ZOH + SCvx con YALMIP/ECOS**: stesso outer-loop ma sottoproblema convesso scritto come SOCP (`norm(U(:,k)) ≤ Tmax`) e risolto da ECOS.

Il file termina con validazione tramite **forward-integration** (i.e., si
re-propaga il controllo ottimo discreto con `ode45` e si misura
l'errore di nodo). I tre/quattro errori sono poi confrontati nei plot.

#### Mapping dei nomi di funzione

Le funzioni si dividono in cinque famiglie. Tutte sono definite localmente
in [main_task2.m](HM2_powered_descent/main_task2.m).

**Dinamica & linearizzazione**

| Funzione | Riferimento | Ruolo |
|---|---|---|
| `rhs` | [main_task2.m:171](HM2_powered_descent/main_task2.m#L171) | RHS non-dim 5-d. Identica a `dyn_rhs` di Task 1 ma scritta su `(x, u, V_c)` separati per usarla nella linearizzazione. |
| `jacobians` | [main_task2.m:177](HM2_powered_descent/main_task2.m#L177) | Jacobiani `A = ∂f/∂x`, `B = ∂f/∂u` al punto di riferimento. Termine `−V_c · T/‖T‖` regolarizzato con `‖T‖_reg = √(T² + 10⁻⁶)` per evitare singolarità in `T=0`. |
| `rk4_zoh` | [main_task2.m:193](HM2_powered_descent/main_task2.m#L193) | Avanza lo stato di un intervallo `dt` con RK4 a `n_sub` sub-step a controllo costante (ZOH non-lineare). |
| `ltv_aug_rhs` | [main_task2.m:205](HM2_powered_descent/main_task2.m#L205) | RHS dell'ODE aumentata di Appendice A: stato `z = [x_ref(t), vec(Φ), vec(B̂), ĉ]` (45 componenti). Integrando su un singolo intervallo `[0,dt]` si recuperano `Ā_k = Φ(dt)`, `B̄_k = B̂(dt)`, `c̄_k = ĉ(dt)` della discretizzazione ZOH lineare. |
| `compute_ltv_zoh` | [main_task2.m:225](HM2_powered_descent/main_task2.m#L225) | Loop su `k=1…N−1`: imposta `z₀` corretto, chiama `ode45` su `ltv_aug_rhs` e stipa in `(Ā[:,:,k], B̄[:,:,k], c̄[:,k])` le matrici discrete. |

**Setup NLP / supporto**

| Funzione | Riferimento | Ruolo |
|---|---|---|
| `init_guess` | [main_task2.m:466](HM2_powered_descent/main_task2.m#L466) | Initial guess: stato lineare in `α=(i−1)/(N−1)`; massa `m₀(1−0.3α)`; controllo costante `(0, m₀·g)` (hover). Se `zero_uN=true` annulla il controllo al nodo finale (ZOH ha un solo controllo in meno). |
| `box_bounds` | [main_task2.m:484](HM2_powered_descent/main_task2.m#L484) | Bounds: `y≥0`, `m∈[10⁻³, m₀]`, `|T_{x,y}| ≤ T_max`. Annulla il controllo del nodo finale se ZOH. |
| `bcs` | [main_task2.m:500](HM2_powered_descent/main_task2.m#L500) | Costruisce le 9 equazioni lineari di B.C. (5 iniziali + 4 finali per posizione/velocità). |
| `fmincon_opts` | [main_task2.m:512](HM2_powered_descent/main_task2.m#L512) | Opzioni standard `fmincon` (`'sqp'`, `MaxIterations = 1000`, tolleranze `1e-5` / `1e-6` / `1e-10`). |
| `unpack` | [main_task2.m:521](HM2_powered_descent/main_task2.m#L521) | Da vettore `z` a struct `sol` con tutti i campi (`x, y, vx, vy, m, Tx, Ty, Tmag, tf, m_f`). |
| `ref_to_z` | [main_task2.m:404](HM2_powered_descent/main_task2.m#L404) | Inversa di `unpack` (utile per warm-start). |
| `apply_trust` | [main_task2.m:418](HM2_powered_descent/main_task2.m#L418) | Stringe le box-bounds attorno alla `ref` corrente: `lb_i = max(lb_i, ref_i − Δ_i)`, idem `ub`. Trust region "hard" (intersezione con i box esistenti). |
| `ternary` | [main_task2.m:400](HM2_powered_descent/main_task2.m#L400) | `if cond, a else, b`. |
| `nondim`, `dim_sol` | [main_task2.m:138](HM2_powered_descent/main_task2.m#L138), [:156](HM2_powered_descent/main_task2.m#L156) | Identici a Task 1. |

**Vincoli (nonlcon callbacks)**

| Funzione | Riferimento | Ruolo |
|---|---|---|
| `path_ineq` | [main_task2.m:559](HM2_powered_descent/main_task2.m#L559) | Vincoli di path **condivisi** tra tutte le trascrizioni: bounds spinta + cono di glide-slope `±x ≤ tan(θ)·y`. |
| `trap_nonlcon` | [main_task2.m:534](HM2_powered_descent/main_task2.m#L534) | Defects trapezoidali + `path_ineq`. Riprodotta qui per indipendenza da Task 1. |
| `zoh_nonlcon` | [main_task2.m:548](HM2_powered_descent/main_task2.m#L548) | Defects ZOH non-lineari `x_{k+1} − RK4(x_k, u_k, dt)` + `path_ineq`. |
| `ltv_nonlcon` | [main_task2.m:296](HM2_powered_descent/main_task2.m#L296) | Solo `path_ineq` (la dinamica LTV è già lineare ed entra come `Aeq, beq`). |

**Solver principali**

| Funzione | Riferimento | Ruolo |
|---|---|---|
| `solve_trap` | [main_task2.m:442](HM2_powered_descent/main_task2.m#L442) | Risolve la trascrizione trapezoidale (gemello di `solve_trapcol` di Task 1 ma con utility-funzioni condivise). |
| `solve_zoh` | [main_task2.m:454](HM2_powered_descent/main_task2.m#L454) | Risolve la trascrizione ZOH-RK4. Stessa struttura, `nonlcon = zoh_nonlcon`. |
| `solve_ltv_nlp` | [main_task2.m:243](HM2_powered_descent/main_task2.m#L243) | Risolve il sotto-problema LTV per **una** iterazione SCvx: la dinamica `x_{k+1} = Ā_k·x_k + B̄_k·u_k + c̄_k` entra come blocco di equazioni lineari sparse (triplet-format per evitare i warning di `sparse`); resta solo `path_ineq` come `c_ineq`. |
| `solve_ltv_nlp_yalmip` | [main_task2.m:689](HM2_powered_descent/main_task2.m#L689) | Stessa cosa di `solve_ltv_nlp` ma in YALMIP: bound spinta come `norm(U(:,k)) ≤ T_max` (SOCP), risolto da ECOS. |
| `solve_scvx` | [main_task2.m:302](HM2_powered_descent/main_task2.m#L302) | **Outer-loop SCvx** con trust-region adattiva (vedi sotto). Inner solver: `solve_ltv_nlp` (fmincon). |
| `solve_scvx_yalmip` | [main_task2.m:747](HM2_powered_descent/main_task2.m#L747) | Stessa logica di `solve_scvx` ma inner solver = `solve_ltv_nlp_yalmip`. |

**SCvx outer-loop in dettaglio** ([solve_scvx](HM2_powered_descent/main_task2.m#L302)):

```
init ref = warm-start dalla trapezoidale
ρ = 1                            ; ρ ∈ [10⁻³, 1]
per ogni iterazione:
   1. trust-region scaling: Δ_i ← ρ·Δ_i^base
   2. linearizza dinamica attorno a ref  (compute_ltv_zoh)
   3. risolvi sotto-problema convesso     (solve_ltv_nlp[_yalmip])
   4. predetto: J_pred = m_f^cand − m_f^ref
   5. attuale : J_act  = m_f(fwd_integrate(cand)) − m_f^ref
   6. η = J_act / J_pred
   7. η < η_l = 0.25  →  rifiuta, ρ ← ρ/2
      η_l ≤ η < η_h  →  accetta
      η ≥ η_h = 0.7   →  accetta, ρ ← min(ρ_max, 2ρ)
   8. accept → ref ← cand; convergenza se ‖Δx‖ < tol
```

**Validazione & plot**

| Funzione | Riferimento | Ruolo |
|---|---|---|
| `fwd_integrate` | [main_task2.m:569](HM2_powered_descent/main_task2.m#L569) | Ri-propaga la dinamica continua (`ode45`, `RelTol=1e-10`) con il controllo ottimo discreto, modalità `'zoh'` (controllo costante) o `'pwl'` (lineare a tratti, per la trapezoidale). |
| `node_err` | [main_task2.m:596](HM2_powered_descent/main_task2.m#L596) | Norma 2 del mismatch sui 4 stati di posizione/velocità per ogni nodo (la massa è esclusa perché in scala diversa). |
| `plot_compare3`, `plot_compare4` | [main_task2.m:601](HM2_powered_descent/main_task2.m#L601), [:810](HM2_powered_descent/main_task2.m#L810) | Confronto multi-traccia: traiettoria, spinta, massa, glide-slope, errore di nodo, convergenza SCvx (`m_f`, `‖Δx‖`, `ρ`, `η` per iterazione). La versione `4` aggiunge la traccia YALMIP/ECOS. |

---

## HM3 — Controllo d'assetto LV a max-q (frequenza)

A differenza di HM1/HM2 (script monolitici con funzioni locali), HM3 è
**modulare**: builder riusabili in file separati, condivisi dai tre
`main_task*.m`. Il flusso è sempre *plant → controller → loop → margini +
simulazione vento*. Il loop è **condizionatamente stabile** (airframe instabile,
polo a `+√A6`), quindi i target `|GM|≈6 dB, |PM|≈30°` si leggono in modulo.

```
main_taskN.m
   ├── load_hw3_params           ← parametri t=72 s (LPV data / Tabella 1)
   ├── build_plant_rigid/full    ← ss [delta, alpha_w] → misure + plot
   ├── design_controller         ← tuning PD su Nichols (fminsearch)
   ├── assemble_loop             ← connect → L (Nichols) + T (simulazione)
   ├── load_wind_profile + simulate_gust_response
   └── plot + export figure (figures/taskN_*.png, theme light)
```

| Funzione | Sorgente | Ruolo |
|----------|----------|-------|
| `load_hw3_params` | [load_hw3_params.m:1](HM3/load_hw3_params.m#L1) | Struct `p` con i parametri a t=72 s; legge le grandezze time-varying da `GreensiteLPV_DATA.mat` (interp) con fallback alla Tabella 1. Opzioni `mu_alpha_scale`/`mu_c_scale` per i corner di Task 3. |
| `build_plant_rigid` | [build_plant_rigid.m:1](HM3/build_plant_rigid.m#L1) | Plant rigido 4 stati `[z,ż,θ,θ̇]` (Eq. 1 senza bending), ingressi `[δ,α_w]`, uscite misure + plot. |
| `build_plant_full` | [build_plant_full.m:1](HM3/build_plant_full.m#L1) | Plant completo 6 stati (+ `η,η̇`); misura INS (Eq. 2) che accoppia il bending nei sensori (`'ins'` default) o feedback vero (`'true'`). |
| `build_tvc` | [build_tvc.m:1](HM3/build_tvc.m#L1) | `W_TVC(s)` 2° ordine (Eq. 3) × ritardo 20 ms via `pade` (ordine 3). |
| `build_notch_filter` | [build_notch_filter.m:1](HM3/build_notch_filter.m#L1) | Filtro Eq. 4: `sgn=-1` lead-lag (fase, non-min-fase), `sgn=+1` notch min-fase (gain-stab). |
| `assemble_loop` | [assemble_loop.m:1](HM3/assemble_loop.m#L1) | Chiude il loop PD (`connect`): ritorna `L` (SISO al break `delta`, convenzione `1+L`) e `T` (`{α_w,θ_ref}→{θ,z,ż,δ}`). |
| `design_controller` | [design_controller.m:1](HM3/design_controller.m#L1) | Tuning `Kp_θ,Kd_θ` con `fminsearch` su `(|GM|-6)²+(|PM|-30)²`, vincolo CL stabile; gain laterali fissi negativi ~1e-3. |
| `load_wind_profile` | [load_wind_profile.m:1](HM3/load_wind_profile.m#L1) | Raffica deterministica 1-cosine; ampiezza da `drywind.mat` (severe) alla quota corrente; `α_w=v_w/V`. |
| `simulate_gust_response` | [simulate_gust_response.m:1](HM3/simulate_gust_response.m#L1) | `lsim` del closed-loop `T` con `α_w(t)`; ritorna `θ,z,ż,δ` + picchi. |
| `main_task1/2/3` | [main_task1.m:1](HM3/main_task1.m#L1), [main_task2.m:1](HM3/main_task2.m#L1), [main_task3.m:1](HM3/main_task3.m#L1) | Entry point: rigido; completo (Step B senza notch instabile → Step C col notch); robustezza 4 corner ±30%. |
| `init_simulink_hm3` | [init_simulink_hm3.m:1](HM3/init_simulink_hm3.m#L1) | Calcola e spinge nel base workspace tutte le variabili (matrici, guadagni, tf TVC/notch, `wind_ts`) per il modello a blocchi. |
| `run_simulink_closed_loop` | [run_simulink_closed_loop.m:1](HM3/run_simulink_closed_loop.m#L1) | Simula `models/hm3_closed_loop.slx` (costruito a mano, vedi `SIMULINK_GUIDE.md`) e sovrappone allo script. |

---

## Convenzioni globali (HM1+HM2)

- **Stato esteso vs ridotto.** Lo stato MATLAB ha sempre i costati (HM1)
  o i controlli (HM2) impacchettati nello stesso vettore di lavoro per
  semplificare lo storage; le funzioni RHS sanno spacchettarlo.
- **Hot-start sistematico.** In HM1 lo sweep di `Q` è risolto per
  continuazione; in HM2 SCvx parte dalla soluzione trapezoidale, e la
  variante YALMIP parte (come prima ref) dalla stessa trapezoidale.
- **Tolleranze.** ODE in HM1: `RelTol=1e-10, AbsTol=1e-12`; in HM2
  `compute_ltv_zoh` usa `1e-8/1e-10` (linearizzazione locale, non serve
  oltre); `fwd_integrate` rimette `1e-10/1e-12` perché è la verifica.
- **Output.** Tutti gli script salvano i PNG in `HM<N>/figures/` con
  prefisso `task<n>_…`. I sorgenti LaTeX in `HM<N>/report/` puntano a
  questi figures via `\graphicspath{{figures/}{../figures/}}`.
