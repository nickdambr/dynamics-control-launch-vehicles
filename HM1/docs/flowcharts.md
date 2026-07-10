# HM1 — Control-Flow Diagrams (ISO 5807)

Detailed, color-coded step-by-step flowcharts of the four entry-point scripts,
the two Appendix validation scripts, and the shared routines in this folder,
obtained by static reading of the source. Symbols follow **ISO 5807** (terminator, process, predefined process,
decision, preparation, data I/O) and are mapped onto Mermaid node shapes —
and colored by category — so the diagrams render natively on GitHub.

## Symbol & color legend

| ISO 5807 symbol | Meaning | Mermaid node | Color |
|---|---|---|---|
| Terminator (stadium) | Start / End | `([ ... ])` | grey |
| Preparation (hexagon) | Setup / loop init | `{{ ... }}` | amber |
| Process / predefined | Operation / solver call | `[ ... ]` · `[[ ... ]]` | blue |
| Decision (rhombus) | Conditional branch | `{ ... }` | yellow |
| Data I/O (parallelogram) | print / figure / export | `[/ ... /]` | green |
| — | Error / abort path | — | red |

---

## 0 · Code architecture (call graph)

How the four entry points sit on the two numerical engines (`fsolve`, `ode45`)
and the dynamic right-hand sides.

```mermaid
flowchart TD
  subgraph entry["Entry points (one script per task)"]
    m1["main_task1.m<br/>single burn · Q sweep"]:::entryC
    m2["main_task2.m<br/>vertical + burn"]:::entryC
    m3["main_task3.m<br/>vertical + burn + coast"]:::entryC
    m4["main_task4.m<br/>two-stage staging"]:::entryC
  end

  fsolve(["fsolve()<br/>BVP root-finder"]):::engine
  ode45(["ode45()<br/>integrator"]):::engine

  subgraph resid["Shooting residuals (local functions)"]
    sh1["shooting1"]:::proc
    sh2["shooting2"]:::proc
    sh3["shooting3"]:::proc
    shs["shooting_single"]:::proc
    sht["shooting_twostage"]:::proc
  end

  subgraph rhs["Dynamic right-hand sides"]
    odeb[["ode_burn.m<br/>powered flight · linear tangent"]]:::io
    odev["ode_vertical"]:::io
    odel["ode_burn_losses"]:::io
  end

  subgraph valid["Appendix validation scripts"]
    vr["validate_rocket_sled.m<br/>min-energy sled · closed-form check"]:::entryC
    vs["validate_staging_corner.m<br/>staging corner condition · 5-unknown solve"]:::entryC
    sled[["sled_ode<br/>sled RHS · u = lam_v/2"]]:::io
  end

  m1 -->|"@shooting1"| fsolve
  m2 -->|"@shooting2"| fsolve
  m3 -->|"@shooting3"| fsolve
  m4 -->|"@shooting_single / @shooting_twostage"| fsolve
  fsolve --> sh1 & sh2 & sh3 & shs & sht
  sh1 & sh2 & sh3 & shs & sht -->|"@ode_burn"| ode45
  ode45 --> odeb
  m2 & m3 -->|"@ode_vertical"| ode45
  m1 -->|"@ode_burn_losses"| ode45
  vr -->|"@sled_residual"| fsolve
  vs -->|"@shooting_inner / @shooting_corner"| fsolve
  vs -.->|"@ode_burn"| odeb
  vr -.-> sled

  classDef entryC fill:#e7d4f7,stroke:#6f42c1,color:#111
  classDef engine fill:#d7d7d7,stroke:#333,color:#111
  classDef proc fill:#cfe2ff,stroke:#1c5d99,color:#111
  classDef io fill:#d1e7dd,stroke:#2e7d4f,color:#111
```

---

## 1 · Task 1 — Single burn arc, Q sweep with continuation

Bidirectional sweep (forward/backward) starting from `Q ≈ 3` with warm-start of
the previous solution; post-processing into three plots.

```mermaid
flowchart TD
  start([Start Task 1]):::term --> setup{{"Setup: c, eta, yf_vec=[0.04,0.05,0.06], Q_vec=80 pts, ODE/fsolve opts"}}:::setup
  setup --> alloc["Preallocate mf_results, sol_results (nQ x nyf)"]:::proc
  alloc --> loopYf{{"for jj = 1 : nyf  (target altitude yf)"}}:::setup
  loopYf -->|body| setp["Set p.c, p.yf; find idx0 (Q closest to 3)"]:::proc
  setp --> guessDec{"jj == 1 ?"}:::decision
  guessDec -->|yes| g1["fixed z_guess = [0.6; 3.8; 14; 0.30]"]:::proc
  guessDec -->|no| g2["z_guess = previous-yf solution at idx0 (fallback to fixed if empty)"]:::proc
  g1 --> solve0
  g2 --> solve0
  solve0[["fsolve(shooting1) at Q(idx0)"]]:::proc --> conv0{"converged? (ef &gt; 0)"}:::decision
  conv0 -->|no| failed[/"print FAILED · skip this yf"/]:::err --> nextYf
  conv0 -->|yes| int0["ode45(ode_burn): store mf, sol at idx0"]:::proc
  int0 --> fwd{{"Sweep forward: for ii = idx0+1 : nQ"}}:::setup
  fwd -->|body| fsolveF[["fsolve(shooting1) warm-start from z_prev"]]:::proc
  fsolveF --> convF{"converged?"}:::decision
  convF -->|yes| storeF["ode45: store mf, sol; update z_prev"]:::proc --> fwd
  convF -->|no| fwd
  fwd -->|done| bwd{{"Sweep backward: for ii = idx0-1 : -1 : 1"}}:::setup
  bwd -->|body| fsolveB[["fsolve(shooting1) warm-start"]]:::proc
  fsolveB --> convB{"converged?"}:::decision
  convB -->|yes| storeB["ode45: store mf, sol; update z_prev"]:::proc --> bwd
  convB -->|no| bwd
  bwd -->|done| nextYf["next jj"]:::setup
  nextYf --> loopYf
  loopYf -->|all yf done| qstar["Q* = argmax mf for each yf"]:::proc
  qstar --> plot1a[/"Plot 1a: mf vs Q"/]:::io
  plot1a --> losses["For each Q (yf=0.04): ode45(ode_burn_losses) -> Wd, Wg, Wt"]:::proc
  losses --> plot1b[/"Plot 1b: velocity losses vs Q"/]:::io
  plot1b --> opt["Q_opt = argmax mf (yf=0.04)"]:::proc
  opt --> dense["Dense ode45; compute phi(t), psi(t)"]:::proc
  dense --> plot1c[/"Plot 1c: trajectory + thrust/velocity angles"/]:::io
  plot1c --> export[/"Export figures to PNG (light theme)"/]:::io
  export --> done([End Task 1]):::term

  classDef term fill:#d7d7d7,stroke:#333,color:#111
  classDef setup fill:#fff3cd,stroke:#b8860b,color:#111
  classDef proc fill:#cfe2ff,stroke:#1c5d99,color:#111
  classDef io fill:#d1e7dd,stroke:#2e7d4f,color:#111
  classDef decision fill:#ffe08a,stroke:#b8860b,color:#111
  classDef err fill:#f8d7da,stroke:#b02a37,color:#111
```

---

## 2 · Task 2 — Vertical climb + optimal burn arc

Phase 1 with an altitude event; Phase 2 BVP with a fallback guess if the first
solve fails.

```mermaid
flowchart TD
  s2([Start Task 2]):::term --> p2{{"Setup: c, eta, y1=1e-4, yf=0.04, Q=2"}}:::setup
  p2 --> ph1["Phase 1 · vertical climb: ode45(ode_vertical) + altitude event @ y1"]:::proc
  ph1 --> st1["Extract t1, vy1, m1 at y = y1"]:::proc
  st1 --> setp2["Phase 2 · set p with initial state (0, y1, 0, vy1, m1)"]:::proc
  setp2 --> guess2["z_guess = [0.6; 3.8; 14; 0.30]"]:::proc
  guess2 --> solve2[["fsolve(shooting2)"]]:::proc
  solve2 --> conv2{"converged? (ef &gt; 0)"}:::decision
  conv2 -->|no| retry["Retry with alternative guess [0.4; 3.0; 10; 0.35]"]:::proc
  retry --> solve2b[["fsolve(shooting2)"]]:::proc
  solve2b --> conv2b{"converged?"}:::decision
  conv2 -->|yes| good2
  conv2b -->|yes| good2["ode45(ode_burn): compute mf, tf_total, payload"]:::proc
  conv2b -->|no| err2[/"print ERROR: BVP not converged"/]:::err
  good2 --> plots2[/"Plots: trajectory (with vertical-climb inset), thrust/velocity angles"/]:::io
  plots2 --> exp2[/"Export figures to PNG"/]:::io
  err2 --> exp2
  exp2 --> e2([End Task 2]):::term

  classDef term fill:#d7d7d7,stroke:#333,color:#111
  classDef setup fill:#fff3cd,stroke:#b8860b,color:#111
  classDef proc fill:#cfe2ff,stroke:#1c5d99,color:#111
  classDef io fill:#d1e7dd,stroke:#2e7d4f,color:#111
  classDef decision fill:#ffe08a,stroke:#b8860b,color:#111
  classDef err fill:#f8d7da,stroke:#b02a37,color:#111
```

---

## 3 · Task 3 — Vertical + burn + coast (physical-root search)

The four guesses are scanned, accepting the first root with `vyc > 0` (physical
ballistic coast); the coast phase is analytic.

```mermaid
flowchart TD
  s3([Start Task 3]):::term --> p3{{"Setup: c, eta, y1, yf, Q, T=c*Q"}}:::setup
  p3 --> ph1["Phase 1 · vertical climb: ode45 + event @ y1 -> t1, vy1, m1"]:::proc
  ph1 --> setp3["Phases 2+3 · set p with state (0, y1, 0, vy1, m1)"]:::proc
  setp3 --> glist["Build guess_list (4 guesses · short burn -> physical root)"]:::proc
  glist --> loopG{{"for each guess in guess_list"}}:::setup
  loopG -->|body| fs3[["fsolve(shooting3)"]]:::proc
  fs3 --> c3{"converged?"}:::decision
  c3 -->|no| loopG
  c3 -->|yes| chk["ode45(ode_burn): evaluate vyc at cutoff"]:::proc
  chk --> phys{"physical root? (vyc &gt; 0)"}:::decision
  phys -->|no| loopG
  phys -->|yes| accept["Accept z_sol · break"]:::proc
  loopG -->|exhausted| c3f{"root found?"}:::decision
  accept --> c3f
  c3f -->|no| err3[/"print ERROR: not converged"/]:::err
  c3f -->|yes| reint["Re-integrate burn (ode_burn) -> cutoff state xc,yc,vxc,vyc,mc"]:::proc
  reint --> coast["Analytic coast: t_coast=vyc; ballistic x,y,vx,vy"]:::proc
  coast --> res3["Compute tf_total, mf=mc, payload; print summary"]:::proc
  res3 --> plots3[/"Plots: trajectory (vertical+burn+coast), angles, mass profile"/]:::io
  plots3 --> exp3[/"Export figures"/]:::io
  err3 --> exp3
  exp3 --> e3([End Task 3]):::term

  classDef term fill:#d7d7d7,stroke:#333,color:#111
  classDef setup fill:#fff3cd,stroke:#b8860b,color:#111
  classDef proc fill:#cfe2ff,stroke:#1c5d99,color:#111
  classDef io fill:#d1e7dd,stroke:#2e7d4f,color:#111
  classDef decision fill:#ffe08a,stroke:#b8860b,color:#111
  classDef err fill:#f8d7da,stroke:#b02a37,color:#111
```

---

## 4 · Task 4 — Optimal staging (two stages)

Single-stage reference, then a bidirectional sweep of the staging time `ts`
outward from the middle of the range.

```mermaid
flowchart TD
  s4([Start Task 4]):::term --> p4{{"Setup: c, eta, yf, Q, T"}}:::setup
  p4 --> ref[["Single-stage reference: fsolve(shooting_single)"]]:::proc
  ref --> refc{"converged?"}:::decision
  refc -->|no| stop[/"error(): abort"/]:::err
  refc -->|yes| refint["ode45(ode_burn) -> mf_single, payload_single"]:::proc
  refint --> sweep{{"Setup sweep: ts_vec=50 pts, ts_max; find idx_mid"}}:::setup
  sweep --> passes{{"for pass = 1:2 (forward, then backward from idx_mid)"}}:::setup
  passes -->|body| loopTs{{"for ii in idx_range"}}:::setup
  loopTs -->|body| setp4["Set p4 with ts; warm-start guess (z_prev or z_ref)"]:::proc
  setp4 --> fs4[["fsolve(shooting_twostage)"]]:::proc
  fs4 --> c4{"converged?"}:::decision
  c4 -->|no| loopTs
  c4 -->|yes| stage["ode45 stage 1 -> staging (drop ms1) -> ode45 stage 2"]:::proc
  stage --> store4["Compute mf, payload; store; update z_prev"]:::proc
  store4 --> loopTs
  loopTs -->|done| passes
  passes -->|both done| findopt["ts_opt = argmax payload"]:::proc
  findopt --> anyv{"any valid solution?"}:::decision
  anyv -->|no| err4[/"print ERROR"/]:::err
  anyv -->|yes| res4["Print optimal staging; compare vs single-stage"]:::proc
  res4 --> plots4[/"Plots: payload vs ts, mf vs ts, optimal trajectory, mass profile"/]:::io
  plots4 --> exp4[/"Export figures"/]:::io
  err4 --> exp4
  exp4 --> e4([End Task 4]):::term

  classDef term fill:#d7d7d7,stroke:#333,color:#111
  classDef setup fill:#fff3cd,stroke:#b8860b,color:#111
  classDef proc fill:#cfe2ff,stroke:#1c5d99,color:#111
  classDef io fill:#d1e7dd,stroke:#2e7d4f,color:#111
  classDef decision fill:#ffe08a,stroke:#b8860b,color:#111
  classDef err fill:#f8d7da,stroke:#b02a37,color:#111
```

---

## 5 · Shared subroutines

Left: the dynamic RHS (`ode_burn`, linear-tangent steering). Right: the common
pattern of every shooting residual — a guard on `tf`, a `try/catch` around the
integration, and the free-time condition `H = 0` imposed algebraically at the
initial instant.

```mermaid
flowchart TD
  subgraph A["ode_burn(t, z, p)"]
    ob([Entry]):::term --> unpack["Unpack vx, vy, m from z"]:::proc
    unpack --> cost["Costates: lam_vx=const; lam_vy=lam_vy0 - lam_y*t; |lam_v|"]:::proc
    cost --> phi["Optimal angle phi = atan2(lam_vy, lam_vx)"]:::proc
    phi --> deriv["dz = [vx; vy; (T/m)cos phi; (T/m)sin phi - 1; -Q; (T/m^2)|lam_v|]"]:::proc
    deriv --> ret([return dz]):::term
  end

  subgraph B["shooting*(z0, p) — common pattern"]
    sh([Entry: unknowns lam_v0, lam_y, (lam_m0), tf]):::term --> guard{"tf in valid range?"}:::decision
    guard -->|no| pen[/"res = 1e6 (penalty)"/]:::err --> rret
    guard -->|yes| pack["Pack costates into pp; ic with lam_m0=1"]:::proc
    pack --> tryint[["try: ode45(ode_burn) 0 -> tf -> zf"]]:::proc
    tryint --> caught{"integration ok?"}:::decision
    caught -->|no| pen
    caught -->|yes| h0["Compute H0 at t0 (free-time condition)"]:::proc
    h0 --> resid["res = [y(tf)-yf; vx(tf)-1; vy(tf); H0]"]:::proc
    resid --> rret([return res]):::term
  end

  classDef term fill:#d7d7d7,stroke:#333,color:#111
  classDef proc fill:#cfe2ff,stroke:#1c5d99,color:#111
  classDef decision fill:#ffe08a,stroke:#b8860b,color:#111
  classDef err fill:#f8d7da,stroke:#b02a37,color:#111
```

> **Note.** `shooting3` departs from the common pattern: 5 unknowns (`lam_m0`
> kept free) plus switching/coast conditions at cutoff instead of `H = 0` at
> `t0`. It is the only variant not reducible to the 4-unknown scheme above.

---

## 6 · Appendix validation scripts

Two standalone checks that back the appendix claims: `validate_rocket_sled.m`
(Appendix A) recovers a closed-form optimum, and `validate_staging_corner.m`
(Appendix C) reproduces the swept staging optimum of Task 4 and shows that
dropping the burnout reference misplaces it.

```mermaid
flowchart TD
  sv([Start validate_rocket_sled]):::term --> pv{{"Setup: tf=2, rf=1/2, vf=0, ODE/fsolve opts"}}:::setup
  pv --> fsv[["fsolve(sled_residual) from lam = [0; 0]"]]:::proc
  fsv --> rep["Recover lam_r0, lam_v0; terminal residual"]:::proc
  rep --> cc["Cross-check: ode45(sled_ode) -> u_num vs u* = 3/4(1 - t)"]:::proc
  cc --> passv{"ef &gt; 0 and ||lam - 1.5|| &lt; 1e-6 ?"}:::decision
  passv -->|yes| okv[/"print PASS (recovers 3/2, 3/2)"/]:::io
  passv -->|no| failv[/"print FAIL"/]:::err
  okv --> ev([End]):::term
  failv --> ev

  subgraph SR["sled_residual(L) · sled_ode"]
    sr([Entry]):::term --> sri["ode45(sled_ode): integrate [0; 0; Lr; Lv] to tf"]:::proc
    sri --> sro(["res = [r(tf)-rf; v(tf)-vf]"]):::term
  end

  classDef term fill:#d7d7d7,stroke:#333,color:#111
  classDef setup fill:#fff3cd,stroke:#b8860b,color:#111
  classDef proc fill:#cfe2ff,stroke:#1c5d99,color:#111
  classDef io fill:#d1e7dd,stroke:#2e7d4f,color:#111
  classDef decision fill:#ffe08a,stroke:#b8860b,color:#111
  classDef err fill:#f8d7da,stroke:#b02a37,color:#111
```

```mermaid
flowchart TD
  ss([Start validate_staging_corner]):::term --> ps{{"Setup: c, eta, yf, Q, T (same as main_task4)"}}:::setup
  ps --> warm[["Step 1 · warm start: fsolve(shooting_inner) at ts0=0.33 (4 unknowns)"]]:::proc
  warm --> wconv{"warm start converged?"}:::decision
  wconv -->|no| werr[/"error(): abort"/]:::err
  wconv -->|yes| aug[["Step 2 · augmented: fsolve(shooting_corner) 5 unknowns + corner residual"]]:::proc
  aug --> ext["propagate -> mf; payload = mf - eta*Q*(tf - ts)"]:::proc
  ext --> passs{"ef &gt; 0 and |ts - 0.336| &lt; 5e-3 ?"}:::decision
  passs -->|yes| oks[/"print PASS (ts=0.336, tf=0.424, mu=0.068)"/]:::io
  passs -->|no| fails[/"print FAIL"/]:::err
  oks --> caut[["Cautionary: fsolve(shooting_corner_wrong) · un-referenced term"]]:::proc
  fails --> caut
  caut --> cautp[/"print spurious ts ≈ 0.225"/]:::io
  cautp --> es([End]):::term

  subgraph PR["propagate(w) + shooting_corner residual"]
    pr([Entry]):::term --> pr1["ode45(ode_burn): stage 1 over [0, ts]"]:::proc
    pr1 --> pr2["jettison: m+ = m- - eta*Q*ts (lam_m continuous)"]:::proc
    pr2 --> pr3["ode45(ode_burn): stage 2 over [ts, tf]"]:::proc
    pr3 --> pr4["res = [y(tf)-yf; vx(tf)-1; vy(tf); H0; corner]"]:::proc
    pr4 --> pro(["corner = eta(lam_m(tf)-lam_m(ts)) - c|lam_v|(1/m+ - 1/m-)"]):::term
  end

  classDef term fill:#d7d7d7,stroke:#333,color:#111
  classDef setup fill:#fff3cd,stroke:#b8860b,color:#111
  classDef proc fill:#cfe2ff,stroke:#1c5d99,color:#111
  classDef io fill:#d1e7dd,stroke:#2e7d4f,color:#111
  classDef decision fill:#ffe08a,stroke:#b8860b,color:#111
  classDef err fill:#f8d7da,stroke:#b02a37,color:#111
```
