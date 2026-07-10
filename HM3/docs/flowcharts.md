# HM3 — Control-Flow Diagrams (ISO 5807)

Detailed, color-coded step-by-step flowcharts of the three task entry-point
scripts, the Monte-Carlo robustness study, the shared plant/filter builders,
and the LPV full-ascent extension in this folder, obtained by static reading of
the source. Symbols follow **ISO 5807** (terminator, process, predefined
process, decision, preparation, data I/O) and are mapped onto Mermaid node
shapes — and colored by category — so the diagrams render natively on GitHub.

Unlike HM1 (indirect ODE shooting), HM3 is a **classical frequency-domain
control** study: the numerical engines are Control System Toolbox functions
(`tf`/`ss`/`connect`, `margin`/`allmargin`, `nichols`, `getLoopTransfer`,
`lsim`) plus a base-MATLAB `fminsearch` margin-matching tuner and a Monte-Carlo
loop. Every task script follows the same spine: **load params → build plant
(rigid/full) → design PD + notch → assemble open/closed loop → frequency
margins / Nichols → gust time response → plots/export**.

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

How the three task scripts and the Monte-Carlo study sit on a shared set of
builders/helpers and on the Control System Toolbox engines. The three task
scripts reuse the *same* PD design (`design_controller`) and the *same* loop
assembly (`assemble_loop`); they differ only in which plant builder they call
and how much of the filter chain they bolt on.

```mermaid
flowchart TD
  subgraph entry["Entry points"]
    m1["main_task1.m<br/>rigid · PD · Nichols"]:::entryC
    m2["main_task2.m<br/>full · TVC+delay+bending · filter trade"]:::entryC
    m3["main_task3.m<br/>parametric corners · fixed controller"]:::entryC
    mc["main_montecarlo.m<br/>probabilistic margins"]:::entryC
  end

  subgraph build["Builders / helpers (shared)"]
    lp[["load_hw3_params.m<br/>Table 1 / LPV @ t=72 s"]]:::proc
    lw[["load_wind_profile.m<br/>1-cosine gust · alpha_w"]]:::proc
    gr[["build_plant_rigid.m<br/>4-state ss"]]:::proc
    gf[["build_plant_full.m<br/>6-state ss + bending/INS"]]:::proc
    bt[["build_tvc.m<br/>2nd-order + Pade delay"]]:::proc
    bn[["build_notch_filter.m<br/>Eq.4 notch / lead-lag"]]:::proc
    dc[["design_controller.m<br/>fminsearch margin match"]]:::proc
    al[["assemble_loop.m<br/>connect -&gt; L, T"]]:::proc
    sg[["simulate_gust_response.m<br/>lsim closed loop"]]:::proc
  end

  subgraph eng["Control System Toolbox engines"]
    cst(["tf / ss / connect / getLoopTransfer"]):::engine
    mar(["margin / allmargin / isstable / pole"]):::engine
    nic(["nichols / nicholsplot / bode / freqresp"]):::engine
    lsi(["lsim"]):::engine
    fms(["fminsearch (base MATLAB)"]):::engine
  end

  m1 & m2 & m3 & mc --> lp & lw
  m1 --> gr
  m2 & m3 & mc --> gf
  m2 & m3 & mc --> bt & bn
  m1 & m2 & m3 & mc --> dc
  m1 & m2 & m3 & mc --> al
  m1 & m2 & m3 & mc --> sg

  gr & gf -->|"ss(A,B,C,D)"| cst
  bt & bn -->|"tf / pade"| cst
  dc -->|"cost = GM/PM margin miss"| fms
  dc --> al
  al -->|"connect + getLoopTransfer"| cst
  al -->|"margins read by caller"| mar
  m1 & m2 & m3 & mc -->|"nichols overlay"| nic
  sg -->|"y = lsim(T,u,t)"| lsi

  classDef entryC fill:#e7d4f7,stroke:#6f42c1,color:#111
  classDef engine fill:#d7d7d7,stroke:#333,color:#111
  classDef proc fill:#cfe2ff,stroke:#1c5d99,color:#111
  classDef io fill:#d1e7dd,stroke:#2e7d4f,color:#111
```

> **Simulink mirror (brief).** `init_simulink_hm3.m` precomputes all block
> parameters (split plant matrices, gains, `tvc_num/den`, `notch_num/den`, wind
> timeseries) and `run_simulink_closed_loop.m` simulates the hand-built
> `models/hm3_closed_loop.slx` and overlays it on the script response. The
> scripts are the source of truth; the `.slx` only replays the same controller
> and plant. Both return cleanly with a message if the model file is absent.

---

## 1 · Task 1 — Rigid LV attitude control at max-q

Single design point (`t = 72 s`, max-q). Rigid 4-state plant, ideal actuator,
PD pitch + weak negative drift feedback tuned to `|GM| ≈ 6 dB` / `|PM| ≈ 30°`,
checked on Nichols and against a wind-gust time response.

```mermaid
flowchart TD
  s1([Start Task 1]):::term --> p1[["load_hw3_params() -> p (A6, K1, V, wBM, ...)"]]:::proc
  p1 --> echo[/"print param source · unstable airframe pole +sqrt(A6)"/]:::io
  echo --> gr1[["build_plant_rigid(p) -> G (4-state ss)"]]:::proc
  gr1 --> dc1[["design_controller(G, []) -> K, m  (ideal actuator)"]]:::proc
  dc1 --> pp["Pole-placement cross-check: wc_eq, zeta_eq vs course range"]:::proc
  pp --> al1[["assemble_loop(G, K) -> L (open), T (closed)"]]:::proc
  al1 --> mrg[/"print stability margins: GM_dB, PM_deg, wc, stable flag"/]:::io
  mrg --> w1[["load_wind_profile(p) -> w (severe 1-cosine gust)"]]:::proc
  w1 --> sg1[["simulate_gust_response(T, w) -> r (theta,z,zdot,delta,alpha)"]]:::proc
  sg1 --> gsum[/"print peak theta, z, delta, alpha -> peak qbar*alpha"/]:::io
  gsum --> f1[/"Fig 1: Nichols chart of L + GM/PM markers (bode at wc)"/]:::io
  f1 --> f2[/"Fig 2: gust response tiles (theta, z, zdot, delta)"/]:::io
  f2 --> f3[/"Fig 3: angle-of-attack budget + aerodynamic load qbar*alpha"/]:::io
  f3 --> exp1[/"Export PNG (light theme, 200 dpi) to figures/"/]:::io
  exp1 --> e1([End Task 1]):::term

  classDef term fill:#d7d7d7,stroke:#333,color:#111
  classDef setup fill:#fff3cd,stroke:#b8860b,color:#111
  classDef proc fill:#cfe2ff,stroke:#1c5d99,color:#111
  classDef io fill:#d1e7dd,stroke:#2e7d4f,color:#111
  classDef decision fill:#ffe08a,stroke:#b8860b,color:#111
  classDef err fill:#f8d7da,stroke:#b02a37,color:#111
```

---

## 2 · Task 2 — Full model: TVC, transport delay, bending notch trade

Extends Task 1 to the 6-state plant (lightly damped first bending mode +
2nd-order TVC + 20 ms Pade delay + INS coupling). Four steps: reuse the rigid
PD (A), expose the resonance instability (B), sweep four bending-filter
candidates (C), then test off-nominal `wBM` (D).

```mermaid
flowchart TD
  s2([Start Task 2]):::term --> mute["warning off MarginUnstable"]:::proc
  mute --> p2[["load_hw3_params() -> p"]]:::proc
  p2 --> stepA["Step A: build_plant_rigid -> design_controller -> assemble_loop -> Trigid (reuse Task-1 PD)"]:::proc
  stepA --> stepB[["Step B: build_plant_full('ins') + build_tvc(p,3); assemble_loop(G,K,Wtvc)"]]:::proc
  stepB --> bchk[/"print |L(wBM)| dB, isstable(Tb) -> resonance destabilises"/]:::io
  bchk --> cprep{{"Step C setup: grids wx=wBM±4, zN=0.1:0.3, zD=0.4:0.6"}}:::setup
  cprep --> cLL{{"C-LL sweep: for each (wx,zN,zD) lead-lag (sgn=-1)"}}:::setup
  cLL -->|body| cLLa[["assemble_loop(G,K,Wtvc*Hc); track isstable, max Re(pole)"]]:::proc
  cLLa --> cLL
  cLL -->|done| cLLr[/"print C-LL: how many stabilise · least-unstable combo"/]:::io
  cLLr --> cN["C-N: deep min-phase notch at wBM (zN=0.002, zD=0.7, sgn=+1) = retained"]:::proc
  cN --> cT["C-T: triplet notch at 0.9/1.0/1.1 wBM"]:::proc
  cT --> cNLL{{"C-NLL sweep: notch Hn * lead-lag partner"}}:::setup
  cNLL -->|body| cNLLa[["assemble_loop(G,K,Wtvc*Hn*Hc); if stable, allmargin -> DelayMargin"]]:::proc
  cNLLa --> keepdm{"stable? higher DM than best?"}:::decision
  keepdm -->|yes| cNLLu["update bestC (max DM)"]:::proc --> cNLL
  keepdm -->|no| cNLL
  cNLL -->|done| ctab["Comparison table: for each candidate assemble_loop + allmargin/margin"]:::proc
  ctab --> ctabp[/"print rigidGM, minGM, PM, DM, |L(wBM)|, stable per candidate"/]:::io
  ctabp --> decide["Decision: deep notch retained (Wfull = Wtvc*Hn) -> Lc, Tfull"]:::proc
  decide --> dsum[/"print retained design margins (GM/PM preserved, extra delay margin)"/]:::io
  dsum --> stepD{{"Step D: scales = 0.90:0.05:1.10 (true wBM off-nominal, filters fixed)"}}:::setup
  stepD -->|body| dscan[["per candidate C-N/C-T/C-NLL: rescale wBM, build_plant_full, assemble_loop, isstable"]]:::proc
  dscan --> stepD
  stepD -->|done| dtab[/"print stability table vs wBM scale"/]:::io
  dtab --> g2[["load_wind_profile + simulate_gust_response(Tfull) and (Trigid)"]]:::proc
  g2 --> g2p[/"print full-model gust peaks"/]:::io
  g2p --> f2a[/"Fig 1: Nichols trade (no filter vs lead-lag vs deep notch)"/]:::io
  f2a --> f2b[/"Fig 2: full-model gust tiles"/]:::io
  f2b --> f2c[/"Fig 3: rigid (Task 1) vs full response comparison"/]:::io
  f2c --> exp2[/"Export PNG to figures/"/]:::io
  exp2 --> e2([End Task 2]):::term

  classDef term fill:#d7d7d7,stroke:#333,color:#111
  classDef setup fill:#fff3cd,stroke:#b8860b,color:#111
  classDef proc fill:#cfe2ff,stroke:#1c5d99,color:#111
  classDef io fill:#d1e7dd,stroke:#2e7d4f,color:#111
  classDef decision fill:#ffe08a,stroke:#b8860b,color:#111
  classDef err fill:#f8d7da,stroke:#b02a37,color:#111
```

---

## 3 · Task 3 — Parametric robustness (corner cases)

Controller **fixed** at the Task-2 design (no re-tuning). Sweep nine cases:
nominal, four `±30 %` box vertices `V1..V4` (the assignment corners), and four
one-at-a-time sensitivities `S1..S4`. Per case: rebuild loop, read GM/PM/DM,
run the gust sim.

```mermaid
flowchart TD
  s3([Start Task 3]):::term --> mute3["warning off MarginUnstable"]:::proc
  mute3 --> fix["Fixed controller: build_plant_rigid(p0) -> design_controller -> K; notch fixed @ p0.wBM"]:::proc
  fix --> fxp[/"print fixed Kp_th, Kd_th, notch params"/]:::io
  fxp --> cases{{"Build 9-case list: Nominal, V1..V4 (corners), S1..S4 (sensitivities)"}}:::setup
  cases --> w3[["load_wind_profile(p0) -> w (same gust for all)"]]:::proc
  w3 --> loop3{{"for i = 1 : nC"}}:::setup
  loop3 -->|body| pscale[["load_hw3_params(mu_alpha_scale, mu_c_scale) -> perturbed p"]]:::proc
  pscale --> gbuild[["build_plant_full('ins') + build_tvc*build_notch_filter -> Wf"]]:::proc
  gbuild --> al3[["assemble_loop(Gf, K, Wf) -> L{i}, T"]]:::proc
  al3 --> mar3["allmargin/margin: rigidGM (low-freq xover), minGM, PM, DM"]:::proc
  mar3 --> sg3[["simulate_gust_response(T, w) -> res{i}"]]:::proc
  sg3 --> row3[/"print case row: mu_a, mu_c, margins, peak theta/z, stable"/]:::io
  row3 --> loop3
  loop3 -->|done| f3a[/"Fig 1: Nichols overlay over Nominal + V1..V4"/]:::io
  f3a --> f3b[/"Fig 2: gust response overlay (pitch + lateral drift)"/]:::io
  f3b --> exp3[/"Export PNG to figures/"/]:::io
  exp3 --> e3([End Task 3]):::term

  classDef term fill:#d7d7d7,stroke:#333,color:#111
  classDef setup fill:#fff3cd,stroke:#b8860b,color:#111
  classDef proc fill:#cfe2ff,stroke:#1c5d99,color:#111
  classDef io fill:#d1e7dd,stroke:#2e7d4f,color:#111
  classDef decision fill:#ffe08a,stroke:#b8860b,color:#111
  classDef err fill:#f8d7da,stroke:#b02a37,color:#111
```

---

## 4 · Monte-Carlo — probabilistic margin robustness

Probabilistic counterpart of Task 3. Five uncertain factors are pre-sampled
(seeded, outside the `parfor`), each draw reassembles `L` with the **fixed**
Task-2 controller and notch, and margins/metrics are collected into
distributions, quantiles, a Nichols cloud, and sensitivity scatters.

```mermaid
flowchart TD
  smc([Start Monte-Carlo]):::term --> cfg{{"Config: N=1500, Nsub=150, wgrid, seed=2026; unc spec (5 factors)"}}:::setup
  cfg --> fixc["Fixed controller: build_plant_rigid(p0) -> design_controller -> K; notch fixed @ nominal wBM"]:::proc
  fixc --> nom["Nominal loop Ln = assemble_loop(Gn,K,Wtvc*Wnotch); nichols_branch -> mn (alignment ref)"]:::proc
  nom --> seed["rng(seed); pre-sample fa,fc,fw,fz,ft = sample_factor(unc, N)"]:::proc
  seed --> alloc["Preallocate rigidGM, minGM, PM, DM, peakTh, peakZ, stab (N x 1)"]:::proc
  alloc --> pf{{"parfor i = 1 : N  (serial without PCT)"}}:::setup
  pf -->|body| pert["Scale p: A6,K1,wBM,zBM,tau by fa,fc,fw,fz,ft(i)"]:::proc
  pert --> reasm[["build_plant_full('ins') + build_tvc*Wnotch; assemble_loop -> L, T"]]:::proc
  reasm --> amm["allmargin(L): minGM, rigidGM (0.2-1 rad/s xover), DM"]:::proc
  amm --> stb["isstable(T); margin(L) -> PM (cap 180 if no phase xover)"]:::proc
  stb --> sgm[["simulate_gust_response(T,w) -> peakTh, peakZ"]]:::proc
  sgm --> pf
  pf -->|done| toc[/"print Monte-Carlo wall time"/]:::io
  toc --> stat["Statistics: Pstab, P(minGM&gt;=3dB), P(PM&gt;=30deg); local_quantile p05/p50/p95; worst draw = argmin minGM"]:::proc
  stat --> statp[/"print robustness summary + percentile table + worst-draw factors"/]:::io
  statp --> f4a[/"Fig 1: 6 histograms (minGM, rigidGM, PM, DM, peakTh, peakZ) w/ target lines"/]:::io
  f4a --> cloud{{"Fig 2 setup: nominal branch shift; subset sub = linspace(1,N,Nsub)"}}:::setup
  cloud -->|body| crebuild[["rebuild Lk per draw; nichols_branch aligned to mn; plot gray=stable / red=unstable"]]:::proc
  crebuild --> cloud
  cloud -->|done| f4b[/"Fig 2: Nichols cloud + nominal + critical points"/]:::io
  f4b --> f4c[/"Fig 3: sensitivity scatter (a) mu plane + Task-3 box, (b) wBM detuning vs minGM"/]:::io
  f4c --> save4[/"Export PNG + save mc_results.mat (factors + per-draw metrics)"/]:::io
  save4 --> emc([End Monte-Carlo]):::term

  classDef term fill:#d7d7d7,stroke:#333,color:#111
  classDef setup fill:#fff3cd,stroke:#b8860b,color:#111
  classDef proc fill:#cfe2ff,stroke:#1c5d99,color:#111
  classDef io fill:#d1e7dd,stroke:#2e7d4f,color:#111
  classDef decision fill:#ffe08a,stroke:#b8860b,color:#111
  classDef err fill:#f8d7da,stroke:#b02a37,color:#111
```

> **Local helpers (`main_montecarlo`).** `sample_factor` draws `gauss`
> (truncated by resampling), `uniform`, or `lognorm` multiplicative factors;
> `nichols_branch` returns aligned `(gain dB, phase deg)` and shifts the phase
> by `360k` to keep the cloud on one branch; `print_pct` / `local_quantile`
> give Statistics-Toolbox-free percentiles; `hist_metric` draws one histogram
> tile with a median line and an optional reference target.

---

## 5 · LTV full-ascent extension (`LTV_FULL_ASCENT/`)

Portfolio showcase beyond the assignment (the deliverable is the max-q point
only). The frozen HM3 design is lifted to a **time-varying (LPV)** plant from
`GreensiteLPV_DATA.mat` so the wind generator (`strong_wind.slx`, 0–140 s) and
the dynamics share one clock. `ode45` is the source of truth; two
script-authored `.slx` mirrors reproduce it.

### 5.0 Call graph

```mermaid
flowchart TD
  subgraph lentry["LPV entry points"]
    mfa["main_full_ascent.m<br/>frozen vs gain-scheduled PD"]:::entryC
    mqs["main_q_scheduling.m<br/>schedule on q(t) vs t"]:::entryC
    mfx["main_flex.m<br/>fixed vs varying bending notch"]:::entryC
  end

  isl[["init_simulink_lpv.m<br/>LPV grids + interpolants (fc1..fc7, fKp/fKd, flex set)"]]:::proc
  oa[["ode_lpv_ascent.m<br/>rigid LTV RHS (ideal actuator PD)"]]:::io
  of[["ode_lpv_flex.m<br/>flex LTV RHS (bending + INS + notch + TVC)"]]:::io

  subgraph reuse["Reused HM3 builders"]
    rr[["build_plant_rigid / build_plant_full"]]:::proc
    rp[["load_hw3_params(t_ref)"]]:::proc
    ra[["assemble_loop / design_controller / build_tvc / build_notch_filter"]]:::proc
  end

  subgraph slx["Simulink mirrors (script-authored)"]
    ba[["build_hm3_full_ascent.m -> hm3_full_ascent.slx"]]:::proc
    bf[["build_hm3_full_ascent_flex.m -> ..._flex.slx"]]:::proc
  end

  mfa & mqs & mfx --> isl
  mfa & mqs --> oa
  mfx --> of
  mfa & mqs & mfx -->|"frozen-time margin sweep"| rr & rp & ra
  oa & of -->|"ode45 inner loop"| ode(["ode45"]):::engine
  isl -.->|"base vars"| ba & bf

  classDef entryC fill:#e7d4f7,stroke:#6f42c1,color:#111
  classDef engine fill:#d7d7d7,stroke:#333,color:#111
  classDef proc fill:#cfe2ff,stroke:#1c5d99,color:#111
  classDef io fill:#d1e7dd,stroke:#2e7d4f,color:#111
```

### 5.1 LPV mains (higher-level flow)

All three share the spine: `init_simulink_lpv` → `ode45` over the ascent (one
run per controller variant) → frozen-time margin sweep (freeze the plant at
each instant, read margins) → summary → figures/export.

```mermaid
flowchart TD
  sl([Start LPV main]):::term --> ini[["init_simulink_lpv() -> S (grids, interpolants, K0, schedule, flex set)"]]:::proc
  ini --> branch{"which script?"}:::decision

  branch -->|main_full_ascent| fa1["ode45(ode_lpv_ascent): frozen-gain run + scheduled-gain run"]:::proc
  fa1 --> fa2["consistency check @ t=72 s: LPV loop reduces to HM3 frozen response"]:::proc

  branch -->|main_q_scheduling| qs1["Build q-keyed gain lookup from ASCENDING branch; quantify hysteresis on descending branch"]:::proc
  qs1 --> qs2["ode45(ode_lpv_ascent) x3: frozen / t-scheduled / q-scheduled"]:::proc

  branch -->|main_flex| fx1["Frozen-time detuning sweep: |L(omega(t))| + isstable, fixed notch @omega(72) vs varying @omega(t)"]:::proc
  fx1 --> fx2["ode45(ode_lpv_flex) x2: fixed-notch run + varying-notch run (find first unstable instant)"]:::proc

  fa2 --> swp["Frozen-time margin sweep: per t, build_plant_*(load_hw3_params(t)) + assemble_loop -> GM/PM"]:::proc
  qs2 --> swp
  fx2 --> swp
  swp --> sum[/"print summary: peak theta/z/delta, qbar*alpha, margins, consistency/hysteresis"/]:::io
  sum --> figs[/"Figures: response, load indicator / detuning / hysteresis, margin sweep -> PNG"/]:::io
  figs --> el([End LPV main]):::term

  classDef term fill:#d7d7d7,stroke:#333,color:#111
  classDef setup fill:#fff3cd,stroke:#b8860b,color:#111
  classDef proc fill:#cfe2ff,stroke:#1c5d99,color:#111
  classDef io fill:#d1e7dd,stroke:#2e7d4f,color:#111
  classDef decision fill:#ffe08a,stroke:#b8860b,color:#111
  classDef err fill:#f8d7da,stroke:#b02a37,color:#111
```

> **LPV RHS files.** `ode_lpv_ascent` (4-state rigid) and `ode_lpv_flex`
> (13-state: 6 plant + 2 notch + TVC) carry **no `arguments` validation by
> design** — they sit in the `ode45` inner loop. The PD is closed inside the
> RHS (`theta_ref = 0`); coefficients are `griddedInterpolant` handles
> evaluated at `t`. `build_hm3_full_ascent[_flex].m` author the `.slx` mirrors
> block-by-block (1-D lookups on flight time, a read-only copy of the
> professor's `strong_wind/Subsystem`, integrator chains); they never modify the
> source generator.

---

## 6 · Shared subroutines / builders

The builders are pure constructors (params → LTI object); `assemble_loop` is
the one place the loop topology lives; `simulate_gust_response` is the only
time-domain call. Left: how `assemble_loop` closes the loop and breaks it for
analysis. Right: the `design_controller` margin-matching search.

```mermaid
flowchart TD
  subgraph A["assemble_loop(G, K, Wact)"]
    al([Entry: plant G, gains K, actuator chain Wact]):::term --> aw{"Wact empty?"}:::decision
    aw -->|yes| aw1["Wact = tf(1)  (ideal actuator)"]:::proc --> kc
    aw -->|no| kc["Build static controller gain block Kc (named IO): u_pd from theta_ref, theta_m, thetadot_m, z_m, zdot_m"]:::proc
    kc --> chain["Actuator/filter chain Wa: u_pd -> delta"]:::proc
    chain --> conn[["T = connect(G, Kc, Wa) keeping 'delta' as analysis point"]]:::proc
    conn --> lt[["L = getLoopTransfer(T,'delta',-1); L = minreal(tf(L))"]]:::proc
    lt --> aret(["return L (open), T (closed), info"]):::term
  end

  subgraph B["design_controller(G, Wact)"]
    dc([Entry: targets GM=6 dB, PM=30 deg; K0 guess; fixed lateral gains]):::term --> mw["mute MarginUnstable warning (onCleanup restore)"]:::proc
    mw --> fs[["fminsearch(cost, log(K0))"]]:::proc
    fs --> cost["cost(x): Kp_th,Kd_th = exp(x); assemble_loop; margin -> (|GM|-6)^2+(|PM|-30)^2; +1e4 if unstable; +1e6 if non-finite"]:::proc
    cost --> fs
    fs --> kbuild["K = exp(xopt); re-assemble; margin -> m struct"]:::proc
    kbuild --> dret([return K, m]):::term
  end

  classDef term fill:#d7d7d7,stroke:#333,color:#111
  classDef proc fill:#cfe2ff,stroke:#1c5d99,color:#111
  classDef decision fill:#ffe08a,stroke:#b8860b,color:#111
```

```mermaid
flowchart TD
  subgraph P["Plant builders"]
    pr(["build_plant_rigid(p)"]):::term --> pra["A,Bd,Bw from a1,a3,a4,A6,K1,V; C = 4 meas + 3 plot rows"]:::proc
    pra --> prr(["ss(A,[Bd Bw],C,D): 4 states [z zdot theta thetadot]"]):::term
    pf(["build_plant_full(p, meas)"]):::term --> pfa["Augment with bending [eta etadot]: -wBM^2, -2 zBM wBM; phi_tvc*Tc forcing"]:::proc
    pfa --> pfm{"meas = 'ins'?"}:::decision
    pfm -->|ins| pfi["Cm leaks bending: sigma_ins, phi_ins into theta_m/z_m (Eq. 2)"]:::proc
    pfm -->|true| pft["Cm = clean states"]:::proc
    pfi --> pfr
    pft --> pfr(["ss: 6 states"]):::term
  end

  subgraph F["Filter / actuator builders"]
    bt2(["build_tvc(p, padeOrder)"]):::term --> bta["Wact = wTVC^2/(s^2+2 zTVC wTVC s+wTVC^2)"]:::proc
    bta --> btp["[nd,dd] = pade(tau, order); Wtvc = Wact * tf(nd,dd)"]:::proc
    btp --> btr(["return Wtvc"]):::term
    bn2(["build_notch_filter(wx, zN, zD, numSign)"]):::term --> bnn["num=[1 sgn*2 zN wx wx^2]; den=[1 2 zD wx wx^2] (Eq. 4)"]:::proc
    bnn --> bnr(["return Hx (sgn=-1 NMP lead-lag / +1 min-phase notch)"]):::term
  end

  subgraph G["Inputs + time-domain"]
    lp2(["load_hw3_params(opt)"]):::term --> lpd{"GreensiteLPV_DATA.mat present?"}:::decision
    lpd -->|yes| lpi["interp LPV coeffs @ t_ref"]:::proc --> lps
    lpd -->|no| lpl["Table 1 literals"]:::proc --> lps
    lps["derive rho, qbar; apply mu_alpha/mu_c scaling"]:::proc --> lpr(["return p"]):::term
    lw2(["load_wind_profile(p, o)"]):::term --> lwg["gust/step/doublet alpha_w = vw/V (or strong_wind.slx windowed)"]:::proc
    lwg --> lwr(["return w (t, vw, alphaw)"]):::term
    sg2(["simulate_gust_response(T, w)"]):::term --> sgl[["y = lsim(T, [alpha_w, theta_ref=0], t)"]]:::proc
    sgl --> sga["alpha = theta + zdot/V + alpha_w; peak_* metrics"]:::proc
    sga --> sgr(["return r"]):::term
  end

  classDef term fill:#d7d7d7,stroke:#333,color:#111
  classDef proc fill:#cfe2ff,stroke:#1c5d99,color:#111
  classDef decision fill:#ffe08a,stroke:#b8860b,color:#111
```

> **Note.** The airframe is open-loop unstable (`+sqrt(A6)`), so the loop is
> *conditionally* stable: `|GM|`/`|PM|` are matched in **magnitude** (Nichols),
> and the binding stability indicator throughout is the `isstable()` flag, not
> the sign of the gain margin. The lateral drift gains `Kp_z, Kd_z` are held
> fixed (small, negative) — only the pitch PD pair is tuned.
