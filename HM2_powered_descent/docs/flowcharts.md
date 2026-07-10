# HM2 — Control-Flow Diagrams (ISO 5807)

Detailed, color-coded step-by-step flowcharts of the two entry-point scripts
(Task 1 trapezoidal collocation, Task 2 ZOH transcriptions) and every shared
routine in this folder, obtained by static reading of the source. Symbols
follow **ISO 5807** (terminator, process, predefined process, decision,
preparation, data I/O) and are mapped onto Mermaid node shapes — and colored by
category — so the diagrams render natively on GitHub.

Task 2 now assembles **five transcriptions** of the same fixed-time
minimum-fuel problem (one trapezoidal baseline plus ZOH variants **a–d**),
which share a non-dimensionalisation, boundary/path structure, and an
`ode45` fidelity replay, but differ in how they discretise the dynamics and in
which convex solver drives the inner step. The diagrams below trace each path
down to the constraint block, the discretisation kernel, and the trust-region
logic.

## Symbol & color legend

| ISO 5807 symbol | Meaning | Mermaid node | Color |
|---|---|---|---|
| Terminator (stadium) | Start / End / return | `([ ... ])` | grey |
| Preparation (hexagon) | Setup / loop init | `{{ ... }}` | amber |
| Process / predefined | Operation / solver call | `[ ... ]` · `[[ ... ]]` | blue |
| Decision (rhombus) | Conditional branch | `{ ... }` | yellow |
| Data I/O (parallelogram) | print / figure / export | `[/ ... /]` | green |
| — | Error / abort / warning path | — | red |

---

## Transcription roster (Task 2)

The five transcriptions and their defining choices — the decision vector,
how the dynamics enter the optimisation, what is held constant over a ZOH
interval, the solver, the warm start, and the `ode45` replay convention used
for the fidelity check.

| # | Name · local fn | Dynamics in the NLP/convex problem | Control hold | Inner solver | Warm start | Replay |
|---|---|---|---|---|---|---|
| — | **Trapezoidal** · `solve_trap` | trapezoidal defects (nonlinear equalities) | PWL thrust `T` | `fmincon` SQP | linear interp + hover | `pwl` |
| a | **Nonlinear ZOH + RK4** · `solve_zoh` | `x_{k+1} = RK4(x_k,u_k)` multiple-shooting defects (nonlinear eq) | ZOH thrust `T` | `fmincon` SQP | linear interp + hover, `u_N=0` | `zoh` |
| b | **LTV + SCvx** · `solve_scvx` | `x_{k+1}=Ā x + B̄ u + c̄` (linear eq, re-linearised each iter) | ZOH thrust `T` | `fmincon` SQP (inner) | trapezoidal baseline | `zoh` |
| c | **LTV + SCvx (YALMIP)** · `solve_scvx_yalmip` | same LTV linear eq; thrust bound as SOC | ZOH thrust `T` | ECOS SOCP (inner) | trapezoidal baseline | `zoh` |
| d | **GFOLD log-mass + SCvx** · `solve_gfold_scvx` | exact **LTI** (one matrix exp); only thrust upper bound linearised | ZOH accel `u=T/m` | ECOS SOCP (inner) | analytic max-thrust profile (self-start) | `u-zoh` |

Variants **c** and **d** are built only when both YALMIP and ECOS are on the
path; otherwise Task 2 runs the trapezoidal baseline plus variants **a** and
**b** and plots three curves instead of five.

---

## 0 · Code architecture (call graph)

How the two entry points sit on the numerical engines (`fmincon` for the NLP
paths, `YALMIP`+`ECOS` for the SOCP inner steps, `expm` for the LTI
discretisation, `ode45` for the augmented-ODE discretisation and every replay),
the per-transcription defect/path locals, and the two continuous right-hand
sides (`ode_descent`, `ode_descent_uacc`).

```mermaid
flowchart TD
  subgraph entry["Entry points — one script per task"]
    m1["main_task1.m<br/>trapezoidal collocation · tf sweep + grid study"]:::entryC
    m2["main_task2.m<br/>5 transcriptions: trap baseline + ZOH a/b/c/d"]:::entryC
  end

  fmincon(["fmincon() · SQP<br/>NLP solver"]):::engine
  yalmip(["YALMIP + ECOS<br/>SOCP solver"]):::engine
  ode45(["ode45()<br/>integrator"]):::engine
  expmE(["expm()<br/>matrix exponential"]):::engine

  subgraph drv["Solve drivers — outer"]
    st["solve_trapcol / solve_trap<br/>baseline"]:::proc
    sz["solve_zoh<br/>variant a"]:::proc
    sc["solve_scvx<br/>variant b"]:::proc
    scy["solve_scvx_yalmip<br/>variant c"]:::proc
    sg["solve_gfold_scvx<br/>variant d"]:::proc
  end

  subgraph inr["Convex inner subproblems"]
    sl["solve_ltv_nlp<br/>LTV eq · fmincon"]:::proc
    sly["solve_ltv_nlp_yalmip<br/>LTV eq · SOCP"]:::proc
    sgs["solve_gfold_socp<br/>LTI eq · SOCP"]:::proc
  end

  subgraph nlc["nonlcon / defect + path"]
    tn["trap_nonlcon"]:::proc
    zn["zoh_nonlcon"]:::proc
    ln["ltv_nonlcon"]:::proc
    pi["path_ineq<br/>thrust + glide-slope"]:::proc
  end

  subgraph rhs["Dynamics RHS / discretization"]
    od[["ode_descent.m<br/>T-hold point-mass RHS"]]:::io
    odu[["ode_descent_uacc.m<br/>u-hold (accel) RHS"]]:::io
    rk[["rk4_zoh.m<br/>RK4 step, T held"]]:::io
    lt["compute_ltv_zoh<br/>ltv_aug_rhs + jacobians"]:::io
    li["lti_zoh.m<br/>van Loan matrix exp"]:::io
  end

  m1 --> st
  m2 --> st & sz & sc & scy & sg
  st & sz --> fmincon
  sc -->|"per iter"| sl --> fmincon
  scy -->|"per iter"| sly --> yalmip
  sg -->|"per iter"| sgs --> yalmip
  st --> tn
  sz --> zn
  sl --> ln
  tn & zn & ln --> pi
  sc & scy -->|"per iter"| lt --> ode45
  sg -->|"once"| li --> expmE
  zn -->|"@rk4_zoh"| rk --> od
  tn --> od
  m1 & m2 -->|"fwd_integrate replay"| ode45
  sg -->|"fwd_integrate_uacc"| ode45
  ode45 --> od & odu

  classDef entryC fill:#e7d4f7,stroke:#6f42c1,color:#111
  classDef engine fill:#d7d7d7,stroke:#333,color:#111
  classDef proc fill:#cfe2ff,stroke:#1c5d99,color:#111
  classDef io fill:#d1e7dd,stroke:#2e7d4f,color:#111
```

---

## 1 · Task 1 — Trapezoidal direct collocation, tf sensitivity sweep

Fixed-duration minimum-fuel NLP solved non-dim with `fmincon` (SQP); a
three-point sweep on flight time (`tf` nominal ±5%), per-solve diagnostics, then
a grid-convergence study replaying the PWL control through `ode45`. The NLP
assembly is factored into §1.1 and the post-solve diagnostics into §1.2.

```mermaid
flowchart TD
  start(["Start Task 1"]):::term --> setup{{"Setup: data (Table 1), tf_nom = 38 s, N = 50"}}:::setup
  setup --> nd["[ref, dnd] = nondim(data)<br/>L=y0, g, t=sqrt(L/g), V=sqrt(gL), m=m0, T=m0 g; Vc=V/c"]:::proc
  nd --> ndp[/"print non-dim reference scales (L,V,t,m,T,Vc)"/]:::io
  ndp --> sweep{{"tf_list = tf_nom * [0.95, 1.00, 1.05]"}}:::setup
  sweep --> loopK{{"for k = 1 : numel(tf_list)"}}:::setup
  loopK -->|"body"| solve[["solve_trapcol(tf_nd, N, dnd) — see §1.1"]]:::proc
  solve --> dimsol["dim_sol: scale non-dim solution back to SI"]:::proc
  dimsol --> diag[["diagnostics(sol, data, N) — see §1.2"]]:::proc
  diag --> rep[/"print m_f, fuel, switch/coast times, glide-slope margin, KKT activity, fmincon stats"/]:::io
  rep --> loopK
  loopK -->|"all tf done"| plots[/"plot_results: trajectory + corridor, |T|, mass, glide-slope angle"/]:::io
  plots --> export[/"force light theme · exportgraphics -> task1_*.png"/]:::io
  export --> tbl[/"print sensitivity summary table (tf, m_f, fuel)"/]:::io
  tbl --> gridL{{"Grid convergence: N_list = [25, 50, 100], nominal tf"}}:::setup
  gridL -->|"body"| gsolve["solve_trapcol(tf_nom/ref.t, N_k) · tic/toc"]:::proc
  gsolve --> replay["fwd_integrate_pwl: replay PWL control via ode45 (RelTol 1e-10, AbsTol 1e-12)"]:::proc
  replay --> nerr["node_err: max pos+vel node error (nondim)"]:::proc
  nerr --> gridL
  gridL -->|"done"| gprint[/"print N, m_f, max err, wall time"/]:::io
  gprint --> done(["End Task 1"]):::term

  classDef term fill:#d7d7d7,stroke:#333,color:#111
  classDef setup fill:#fff3cd,stroke:#b8860b,color:#111
  classDef proc fill:#cfe2ff,stroke:#1c5d99,color:#111
  classDef io fill:#d1e7dd,stroke:#2e7d4f,color:#111
  classDef decision fill:#ffe08a,stroke:#b8860b,color:#111
  classDef err fill:#f8d7da,stroke:#b02a37,color:#111
```

### 1.1 · NLP assembly (`solve_trapcol` / `solve_trap` / `solve_zoh`)

The three fmincon-driven drivers share the same skeleton: decision vector
`z = [x;y;vx;vy;m;Tx;Ty]` stacked node-by-node (length `7N`, `idx(i)=(i-1)*7+(1:7)`),
a linear-interp state + hover-thrust initial guess, box bounds, a 9-row linear
equality block for the boundary conditions, the maximise-final-mass objective,
and a transcription-specific nonlinear-constraint handle. They differ only in
that constraint handle (and in whether the last node's control is pinned to
zero for the ZOH variants).

```mermaid
flowchart TD
  a(["Enter solve_trapcol / solve_trap / solve_zoh"]):::term --> dt["dt = tf/(N-1); nz = 7N; idx(i) = (i-1)*7 + (1:7)"]:::proc
  dt --> ig["init_guess: linear state interp, m = m0(1 - 0.3a),<br/>Tx=0, Ty=m0 (hover, gravity=1 nondim); ZOH -> u(N)=0"]:::proc
  ig --> bnd["box_bounds: y &gt;= 0, m in [1e-3, m0],<br/>|Tx|,|Ty| &lt;= Tmax; ZOH -> u(N) pinned to 0"]:::proc
  bnd --> aeq["bcs: Aeq (9 rows) = full state at node 1 (5) + pos/vel = 0 at node N (4)"]:::proc
  aeq --> obj["objective f(z) = -z(iN_m), iN_m = (N-1)*7 + 5   (maximize m_N)"]:::proc
  obj --> nlc{"which transcription?"}:::decision
  nlc -->|"trap"| tnl["nonlcon = trap_nonlcon (trapezoidal defects + path) — see §3.2"]:::proc
  nlc -->|"ZOH variant a"| znl["nonlcon = zoh_nonlcon (RK4 shooting defects + path) — see §3.2"]:::proc
  tnl --> run
  znl --> run
  run[["fmincon (SQP): MaxIter 1e3, MaxFunEval 1e6,<br/>OptTol 1e-5, ConTol 1e-6, StepTol 1e-10"]]:::proc
  run --> ef{"exitflag &lt;= 0 ?"}:::decision
  ef -->|"yes"| warn[/"warning: fmincon did not converge cleanly"/]:::err --> up
  ef -->|"no"| up["unpack z_opt -> sol (t, x..m, Tx, Ty, Tmag, m_f)<br/>+ exitflag, iters, first-order optimality, lambda"]:::proc
  up --> r(["return sol"]):::term

  classDef term fill:#d7d7d7,stroke:#333,color:#111
  classDef proc fill:#cfe2ff,stroke:#1c5d99,color:#111
  classDef decision fill:#ffe08a,stroke:#b8860b,color:#111
  classDef err fill:#f8d7da,stroke:#b02a37,color:#111
```

### 1.2 · Post-solve diagnostics (`diagnostics`)

Runs on the SI solution: burn/coast structure from `|T|` crossings, the
glide-slope margin, and KKT activity read from the `fmincon` inequality
multipliers.

```mermaid
flowchart TD
  d(["Enter diagnostics(sol, data, N)"]):::term --> thr["thr = 0.5 * Tmax"]:::proc
  thr --> sw["find first |T| down-crossing (burn 1 end) and last up-crossing (burn 2 start)"]:::proc
  sw --> cross["linear-interpolate crossing times -> t_sw1, t_sw2;  coast = t_sw2 - t_sw1"]:::proc
  cross --> gs["mask nodes y &gt; 1 m; th = atan2(|x|, y);<br/>gs_margin = theta_mx - max(th)   [deg]"]:::proc
  gs --> kkt["KKT from lambda.ineqnonlin, stacked [thr_lo; thr_hi; gs_pos; gs_neg] (N each):<br/>n_thr_active = #(thr_hi mult &gt; 1e-6);  max_gs_mult = max(gs mult)"]:::proc
  kkt --> r(["return dg (t_sw1, t_sw2, coast, gs_margin, n_thr_active, max_gs_mult)"]):::term

  classDef term fill:#d7d7d7,stroke:#333,color:#111
  classDef proc fill:#cfe2ff,stroke:#1c5d99,color:#111
```

> **NLP structure recap.** Objective `-z(iN_m)` maximises final mass; the linear
> equality `Aeq` fixes the full state at node 1 and pos+vel = 0 at node N; box
> bounds enforce `y>=0`, `m` in `[1e-3, m0]`, `|Tx|,|Ty| <= Tmax`; the nonlinear
> constraints come from `trap_nonlcon` (collocation defects `= 0`, thrust +
> glide-slope path `<= 0`) or, for variant a, `zoh_nonlcon` (RK4 shooting
> defects `= 0`, same path block).

---

## 2 · Task 2 — Five ZOH transcriptions

Five transcriptions solved on the same grid and compared: the trapezoidal
baseline, then (a) nonlinear ZOH with RK4 multiple-shooting defects; (b)
LTV-linearised ZOH inside a successive-convexification (SCvx) outer loop with an
adaptive trust region (`fmincon` inner NLP); (c) the same SCvx loop with a
YALMIP/ECOS SOCP inner subproblem; (d) the GFOLD log-mass change of variables,
whose dynamics become exactly LTI so the discretisation is a single matrix
exponential and only the thrust upper bound is linearised (SOCP via
YALMIP/ECOS). Variants a–c warm-start from the trapezoidal baseline; variant d
self-starts. Variants c and d run only if YALMIP + ECOS are on the path.

```mermaid
flowchart TD
  s2(["Start Task 2"]):::term --> p2{{"Setup: data, tf=38 s, N=50, n_sub=2, scvx_max=15, scvx_tol=1e-3"}}:::setup
  p2 --> nd2["[ref, dnd] = nondim; tf_nd = tf/ref.t;<br/>trust = (pos .17, vel .6, mass .1, thrust 1.0)"]:::proc
  nd2 --> base[["baseline: solve_trap -> sol_trap_nd   (§1.1)"]]:::proc
  base --> va[["variant a: solve_zoh(tf,N,dnd,n_sub) — RK4 shooting NLP   (§1.1)"]]:::proc
  va --> vb[["variant b: solve_scvx(... sol_trap_nd, trust) — LTV+SCvx, fmincon   (§2.1)"]]:::proc
  vb --> ycheck{"YALMIP + ECOS on path?"}:::decision
  ycheck -->|"no"| ywarn[/"warning: variants c & d skipped"/]:::err --> valid
  ycheck -->|"yes"| vc[["variant c: solve_scvx_yalmip(... sol_trap_nd, trust) — LTV+SCvx, SOCP   (§2.1)"]]:::proc
  vc --> vd[["variant d: solve_gfold_scvx(tf,N,dnd,scvx_max,scvx_tol) — GFOLD log-mass, SOCP   (§2.2)"]]:::proc
  vd --> valid
  valid["fwd_integrate replay per transcription: 'pwl' trap, 'zoh' a/b/c,<br/>u-ZOH (fwd_integrate_uacc) for d"]:::proc
  valid --> ferr["node_err: max grid-node nondim state error per transcription"]:::proc
  ferr --> fprint[/"print transcription fidelity table"/]:::io
  fprint --> land["landing accuracy: pos/vel error norms + m_f drift at t_f (SI)"]:::proc
  land --> lprint[/"print replay landing accuracy + wall time"/]:::io
  lprint --> pbranch{"yalmip_ok?"}:::decision
  pbranch -->|"yes"| p5[/"plot_compare5: 5 transcriptions + 3 SCvx traces (fmincon, YALMIP, GFOLD)"/]:::io
  pbranch -->|"no"| p3[/"plot_compare3: 3 transcriptions + 1 SCvx trace"/]:::io
  p5 --> export2
  p3 --> export2
  export2[/"force light theme · exportgraphics -> task2_*.png"/]:::io
  export2 --> e2(["End Task 2"]):::term

  classDef term fill:#d7d7d7,stroke:#333,color:#111
  classDef setup fill:#fff3cd,stroke:#b8860b,color:#111
  classDef proc fill:#cfe2ff,stroke:#1c5d99,color:#111
  classDef io fill:#d1e7dd,stroke:#2e7d4f,color:#111
  classDef decision fill:#ffe08a,stroke:#b8860b,color:#111
  classDef err fill:#f8d7da,stroke:#b02a37,color:#111
```

### 2.1 · SCvx outer loop (variants b and c)

Successive convexification with an adaptive trust region. `solve_scvx` uses an
`fmincon` inner NLP (LTV dynamics as linear equalities); `solve_scvx_yalmip` is
identical but solves an SOCP inner subproblem with YALMIP/ECOS. Each iteration
re-linearises about the current reference (`compute_ltv_zoh`, §3.3), solves the
convex subproblem inside the trust region, then validates the step against the
**nonlinear** dynamics by forward integration. The trust ratio
`eta = (actual m_f gain)/(predicted m_f gain)` drives accept/reject and
trust-region resizing.

```mermaid
flowchart TD
  sc([Enter solve_scvx / solve_scvx_yalmip]):::term --> init{{"Init: ref = init_ref (trap baseline) or init_guess;<br/>rho=1, rho_min=1e-3, rho_max=1; eta_l=0.25, eta_h=0.7"}}:::setup
  init --> loop{{"for iter = 1 : max_iter"}}:::setup
  loop -->|"body"| scaled["scaled trust = rho * base_trust (pos,vel,mass,thrust)"]:::proc
  scaled --> ltv["compute_ltv_zoh(ref): integrate Appendix A ODE -> Abar, Bbar, cbar   (§3.3)"]:::proc
  ltv --> inner[["Inner convex solve: solve_ltv_nlp (fmincon) OR solve_ltv_nlp_yalmip (SOCP)"]]:::proc
  inner --> jpred["J_pred = sol_cand.m_f - ref.m_f   (LTV-predicted gain)"]:::proc
  jpred --> fint["fwd_integrate(sol_cand,'zoh'): nonlinear ode45 replay -> m_f_actual"]:::proc
  fint --> jact["J_act = m_f_actual - ref.m_f; eta = J_act / J_pred  (1 if |J_pred|&lt;1e-10)"]:::proc
  jact --> dx["delta_x = ||[x;y;vx;vy;m]_cand - ref||"]:::proc
  dx --> trace[/"print iter, rho, eta, delta_x, m_f, ACCEPTED/rejected"/]:::io
  trace --> acc{"eta &gt;= eta_l ?  (accept step)"}:::decision
  acc -->|"no"| shrink["reject: rho = 0.5 * rho"]:::proc
  shrink --> coll{"rho &lt; rho_min ?  (trust collapsed)"}:::decision
  coll -->|"yes"| stop1[/"print: trust region collapsed · break"/]:::err --> ret
  coll -->|"no"| loop
  acc -->|"yes"| keep["accept: sol_best = ref = sol_cand"]:::proc
  keep --> grow{"eta &gt; eta_h ?"}:::decision
  grow -->|"yes"| up["rho = min(rho_max, 2*rho)"]:::proc --> conv
  grow -->|"no"| conv{"delta_x &lt; tol ?"}:::decision
  conv -->|"yes"| done2[/"print: SCvx converged · break"/]:::io --> ret
  conv -->|"no"| loop
  loop -->|"iter == max_iter"| cap[/"print: hit iteration cap"/]:::io
  cap --> ret([return sol_best, iter, conv_hist]):::term

  classDef term fill:#d7d7d7,stroke:#333,color:#111
  classDef setup fill:#fff3cd,stroke:#b8860b,color:#111
  classDef proc fill:#cfe2ff,stroke:#1c5d99,color:#111
  classDef io fill:#d1e7dd,stroke:#2e7d4f,color:#111
  classDef decision fill:#ffe08a,stroke:#b8860b,color:#111
  classDef err fill:#f8d7da,stroke:#b02a37,color:#111
```

> **Inner LTV subproblem.** `solve_ltv_nlp` / `solve_ltv_nlp_yalmip` treat the
> discrete LTV dynamics `x_{k+1} = Abar_k x_k + Bbar_k u_k + cbar_k` as **linear
> equality constraints**, so `ltv_nonlcon` carries only the path constraints
> (`path_ineq`). The fmincon variant keeps the thrust bound as the nonlinear
> path inequality `|T| <= Tmax`; the YALMIP variant writes it as a second-order
> cone `norm(U_k) <= Tmax` and the glide-slope as linear inequalities. Box
> bounds are intersected with the per-variable trust region via `apply_trust`
> (fmincon) or inline trust inequalities (YALMIP).

### 2.2 · GFOLD log-mass SCvx loop (variant d)

The change of variables `z = ln(m)`, `u = T/m`, with slack `sigma >= ||u||`,
makes the translational dynamics **exactly LTI** (`ode_descent_uacc`), so the
Appendix A ZOH collapses to one matrix exponential (`lti_zoh`, §3.4) computed
**once** — no per-iteration re-linearisation of the dynamics and no singular
mass row. Only the thrust *upper* bound `sigma <= Tmax·e^{-z}` is nonconvex; it
is linearised about the current `z_ref` each iteration. The loop self-starts
from an analytic max-thrust mass profile (no trapezoidal warm start), leaves the
first solve free of the trust region, and validates each step with a
`u = T/m`-hold replay (`fwd_integrate_uacc`).

```mermaid
flowchart TD
  g(["Enter solve_gfold_scvx(tf, N, d, max_iter, tol)"]):::term --> lti["[Abar,Bbar,cbar] = lti_zoh(dt, Vc) — exact LTI ZOH, computed ONCE   (§3.4)"]:::proc
  lti --> ref0{{"self-start reference: m_apri = max(m0 - Vc·Tmax·t, 1e-2), z = ln(m_apri),<br/>linear pos/vel, ux=0, uy=1, sig=1 — NO trapezoidal warm start"}}:::setup
  ref0 --> init{{"base trust = (pos .5, vel 1.0, lz .4, u 4.0, sig 4.0);<br/>rho=1, rho_min=1e-3, rho_max=1; eta_l=0.25, eta_h=0.7"}}:::setup
  init --> loop{{"for iter = 1 : max_iter"}}:::setup
  loop -->|"iter == 1"| free["solve_gfold_socp(... ref.z, [], []) — FREE of trust region"]:::proc
  loop -->|"iter &gt; 1"| trr["solve_gfold_socp(... ref.z, ref, rho*base) — inside trust region"]:::proc
  free --> jp
  trr --> jp
  jp["J_pred = cand.m_f - ref.m_f"]:::proc
  jp --> fint["fwd_integrate_uacc(cand): nonlinear u-ZOH ode45 replay -> m_f_actual   (§3.1)"]:::proc
  fint --> eta["J_act = m_f_actual - ref.m_f; eta = J_act/J_pred  (1 if |J_pred|&lt;1e-10)"]:::proc
  eta --> dx["delta_x = ||[x;y;vx;vy;z]_cand - ref||"]:::proc
  dx --> trace[/"print iter, rho, eta, delta_x, m_f, ACCEPTED/rejected"/]:::io
  trace --> acc{"eta &gt;= eta_l ?"}:::decision
  acc -->|"no"| shrink["reject: rho = 0.5 * rho"]:::proc --> coll{"rho &lt; rho_min ?"}:::decision
  coll -->|"yes"| stop[/"print: trust region collapsed · break"/]:::err --> ret
  coll -->|"no"| loop
  acc -->|"yes"| keep["accept: sol_best = ref = cand"]:::proc --> grow{"eta &gt; eta_h ?"}:::decision
  grow -->|"yes"| up2["rho = min(rho_max, 2*rho)"]:::proc --> conv
  grow -->|"no"| conv{"delta_x &lt; tol ?"}:::decision
  conv -->|"yes"| done[/"print: SCvx-GFOLD converged · break"/]:::io --> ret
  conv -->|"no"| loop
  loop -->|"iter == max_iter"| cap[/"print: hit iteration cap"/]:::io --> ret
  ret(["return sol_best, iter, conv_hist"]):::term

  classDef term fill:#d7d7d7,stroke:#333,color:#111
  classDef setup fill:#fff3cd,stroke:#b8860b,color:#111
  classDef proc fill:#cfe2ff,stroke:#1c5d99,color:#111
  classDef io fill:#d1e7dd,stroke:#2e7d4f,color:#111
  classDef decision fill:#ffe08a,stroke:#b8860b,color:#111
  classDef err fill:#f8d7da,stroke:#b02a37,color:#111
```

> **GFOLD inner SOCP (`solve_gfold_socp`).** State `XI = [x;y;vx;vy;z]`
> (`z = ln m`), control `W = [ux;uy;sigma]` (`u = T/m`). Constraints: initial
> condition `XI(:,1) = [x0;y0;vx0;vy0; z0=ln m0=0]` and terminal `XI(1:4,N)=0`
> (`z_N` free); LTI dynamics `XI_{k+1}=Abar·XI_k+Bbar·W_k+cbar` (equalities);
> the **lossless** cone `||u_k|| <= sigma_k` (exact SOC); the **linearised**
> upper thrust bound `sigma_k <= Tmax·e^{-z_ref}(1-(z_k - z_ref))` — the only
> nonconvexity; glide-slope cone (linear), `y>=0`, `z in [ln 1e-3, 0]`, and an
> optional trust box `(pos,vel,lz,u,sig)`. Objective `-XI(5,N)` maximises the
> terminal log-mass `z_N = ln(m_f)`. Recovered as `m = e^z`, `T = m·u`, with the
> last node's control padded to zero.

---

## 3 · Shared subroutines

The continuous right-hand sides, the RK4 propagator, the common defect/path
pattern, and the two discretisation kernels — all shared across Task 1, Task 2,
and the test suite.

### 3.1 · Continuous RHS (two ZOH conventions) + RK4

Two point-mass RHS files differ only in what is held constant over a ZOH
interval: `ode_descent` holds the **thrust vector** `T` (so acceleration
`T/m` grows as mass depletes), while `ode_descent_uacc` holds the
**acceleration** `u = T/m` (so `T = m(t)·u` floats). Both use gravity `= -1`
nondim. `rk4_zoh` is the fixed-step RK4 propagator used by the variant-a
shooting defects, holding `T` constant.

```mermaid
flowchart LR
  subgraph A["ode_descent(x,u,Vc) — thrust T held constant"]
    od(["entry"]):::term --> odm["Tmag = sqrt(Tx^2 + Ty^2)"]:::proc
    odm --> odd["dx = [vx; vy; Tx/m; Ty/m - 1; -Vc*Tmag]"]:::proc
    odd --> odr(["return dx (5x1)"]):::term
  end
  subgraph B["ode_descent_uacc(x,u,Vc) — accel u=T/m held constant"]
    ou(["entry"]):::term --> oum["umag = sqrt(ux^2 + uy^2)"]:::proc
    oum --> oud["dx = [vx; vy; ux; uy - 1; -Vc*m*umag]"]:::proc
    oud --> our(["return dx (5x1)"]):::term
  end
  subgraph C["rk4_zoh(x,u,dt,Vc,n_sub) — T held"]
    rk(["entry"]):::term --> sub{{"for ii = 1 : n_sub  (h = dt/n_sub)"}}:::setup
    sub -->|"body"| k1["k1..k4 via ode_descent, u constant"]:::proc
    k1 --> step["x = x + (h/6)(k1 + 2k2 + 2k3 + k4)"]:::proc
    step --> sub
    sub -->|"done"| rkr(["return x_next (5x1)"]):::term
  end

  classDef term fill:#d7d7d7,stroke:#333,color:#111
  classDef setup fill:#fff3cd,stroke:#b8860b,color:#111
  classDef proc fill:#cfe2ff,stroke:#1c5d99,color:#111
```

### 3.2 · nonlcon defect/path pattern

Every fmincon transcription reshapes `z` into the `7 x N` state/control grid,
forms the dynamics defects (equality `c_eq`), and appends the shared
`path_ineq` block (inequality `c_ineq`). The LTV inner NLP contributes no
defect (its dynamics are linear equalities in `Aeq`), so `c_eq = []`.

```mermaid
flowchart TD
  nc([Entry: trap_nonlcon / zoh_nonlcon / ltv_nonlcon]):::term --> resh["Z = reshape(z, 7, N)"]:::proc
  resh --> which{"which defect ?"}:::decision
  which -->|"trap"| dft["c_eq = Z(:,k+1)-Z(:,k) - 0.5 dt (f_k + f_(k+1))   via ode_descent"]:::proc
  which -->|"zoh (variant a)"| dfz["c_eq = Z(:,k+1) - rk4_zoh(Z_k, u_k, dt, Vc, n_sub)"]:::proc
  which -->|"ltv (variant b)"| dfl["c_eq = []   (dynamics are linear equalities in Aeq)"]:::proc
  dft --> pth
  dfz --> pth
  dfl --> pth["path_ineq(Z, d)"]:::proc
  pth --> pp["c_ineq = [Tmin-|T|; |T|-Tmax; x - tan(theta)*y; -x - tan(theta)*y]"]:::proc
  pp --> ncr([return c_ineq, c_eq]):::term

  classDef term fill:#d7d7d7,stroke:#333,color:#111
  classDef proc fill:#cfe2ff,stroke:#1c5d99,color:#111
  classDef decision fill:#ffe08a,stroke:#b8860b,color:#111
```

### 3.3 · LTV discretization (`compute_ltv_zoh`, Appendix A)

Used by the SCvx variants b and c. It does **not** use `rk4_zoh`: it integrates
the Appendix A augmented ODE (`ltv_aug_rhs`) over each interval with `ode45` in
the beta-gamma form, carrying the reference state, the transition matrix `Phi`,
the inverse transition `Psi = Phi^-1`, and the start-referenced integrals
`Beta = ∫ Psi·B` and `Gamma = ∫ Psi·c`. The Jacobians `df/dx`, `df/du`
(`jacobians`) feed the augmented RHS; `df/du` regularises `|T|` with a small
`1e-6` term to stay finite at `T = 0`.

```mermaid
flowchart TD
  c(["Enter compute_ltv_zoh(ref, tf, N, d)"]):::term --> loop{{"for k = 1 : N-1   (interval [t_k, t_{k+1}], dt = tf/(N-1))"}}:::setup
  loop -->|"body"| z0["augment z0 = [x_k; vec I(5); vec I(5); zeros(10); zeros(5)]<br/>= (state, Phi=I, Psi=I, Beta=0, Gamma=0)"]:::proc
  z0 --> ode["ode45 ltv_aug_rhs over [0, dt] (RelTol 1e-8, AbsTol 1e-10):<br/>dx=f, dPhi=A·Phi, dPsi=-Psi·A, dBeta=Psi·B, dGamma=Psi·c_off"]:::proc
  ode --> jac["ltv_aug_rhs calls jacobians(x_ref, u_k, Vc):<br/>A=df/dx, B=df/du (|T| reg. by 1e-6); c_off = f - A·x - B·u"]:::proc
  jac --> read["read interval end: Phi_f = reshape(zf);<br/>Abar_k = Phi_f, Bbar_k = Phi_f·Beta, cbar_k = Phi_f·Gamma"]:::proc
  read --> loop
  loop -->|"done"| r(["return Abar (5x5x(N-1)), Bbar (5x2x(N-1)), cbar (5x(N-1))"]):::term

  classDef term fill:#d7d7d7,stroke:#333,color:#111
  classDef setup fill:#fff3cd,stroke:#b8860b,color:#111
  classDef proc fill:#cfe2ff,stroke:#1c5d99,color:#111
```

### 3.4 · LTI discretization (`lti_zoh`, GFOLD log-mass)

Used by variant d. Because `z = ln(m)`, `u = T/m` linearises the dynamics
*exactly*, the discretisation is a single van Loan block matrix exponential
computed **once** and reused for every interval (the system is time-invariant).
No per-interval integration, no singular mass row.

```mermaid
flowchart TD
  l(["Enter lti_zoh(dt, Vc)"]):::term --> mats["A: x&lt;-vx, y&lt;-vy;  B: vx&lt;-ux, vy&lt;-uy, z&lt;--Vc·sigma;  c = [0;0;0;-1;0]"]:::proc
  mats --> exp["E = expm([A B c; zeros(4,9)] * dt)   (van Loan augmented block)"]:::proc
  exp --> split["Abar = E(1:5,1:5);  Bbar = E(1:5,6:8);  cbar = E(1:5,9)"]:::proc
  split --> r(["return Abar (5x5), Bbar (5x3), cbar (5x1) — CONSTANT across the grid (LTI)"]):::term

  classDef term fill:#d7d7d7,stroke:#333,color:#111
  classDef proc fill:#cfe2ff,stroke:#1c5d99,color:#111
```

> **Why two discretisers.** The LTV path (§3.3) linearises the *nonlinear*
> thrust dynamics about a moving reference, so its matrices change every SCvx
> iteration and every interval, and are obtained by integrating the augmented
> STM ODE. The GFOLD path (§3.4) works in coordinates where the dynamics are
> already linear and time-invariant, so a single `expm` gives the exact ZOH
> matrices and the SCvx loop only has to chase the one linearised inequality
> (`sigma <= Tmax·e^{-z}`).
