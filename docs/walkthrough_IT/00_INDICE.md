# Walkthrough del codice -- materiale di studio per l'orale (DCLV)

Documento **locale, per la preparazione dell'orale** -- non fa parte della consegna.
Un file `.md` per ogni file di codice dei quattro homework: ruolo del file nel
progetto, poi spiegazione funzione-per-funzione / riga-per-riga, con la matematica
derivata (non solo nominata) e riquadri **"Possibili domande d'esame"** sui punti
delicati.

**Come usarlo.** Se il professore indica un file o una riga, apri il `.md`
corrispondente qui sotto e cerca il sotto-titolo della funzione: ogni sezione
riporta i **numeri di riga reali** del sorgente.

> **Copertura:** 55 file (~8.700 righe di MATLAB) su HM0, HM1, HM2, HM3 e
> l'estensione LPV. Scritto contro il commit **`b673640`**. I numeri di riga
> invecchiano quando il codice cambia: se una sezione non torna, ricontrolla il
> sorgente prima di fidarti della pagina.

> **Onesta'.** Le pagine documentano il codice **com'e'**, non come dovrebbe
> essere. Dove il codice ha un bug, un pezzo stale, un commento che contraddice
> l'implementazione o una deviazione dalla teoria, la pagina lo dice. Sono
> esattamente i punti su cui un esaminatore mette il dito: meglio arrivarci
> sapendolo.

---

## 1 - HM0 -- Falcon 9, ascesa 3-DoF (`HM0_falcon9_ascent/`)

Propagazione in avanti, niente ottimizzazione: e' il banco su cui si impara il modello.

| file | cosa spiega |
|---|---|
| [main.m](hm0_main.md) | EOM in coordinate sferiche derivate da zero, termini di trasporto vs Coriolis, spinta corretta per contropressione, gravity turn, max-Q |
| [main2.m](hm0_main2.md) | adimensionalizzazione, riparametrizzazione in `tau` a tre archi con chain rule, e l'event function inerte smascherata |

## 2 - HM1 -- ascesa planare, metodo indiretto (`HM1/`)

PMP, costati, funzione di switching, BVP risolto con shooting. La teoria pesante.

| file | cosa spiega |
|---|---|
| [ode_burn.m](hm1_ode_burn.md) | **il cuore**: Hamiltoniana, Eulero-Lagrange, primer vector, legge linear-tangent; perche' niente `arguments` block |
| [main_task1.m](hm1_main_task1.md) | shooting a 4 incognite, `H(0)=0` algebrica, continuazione su `Q` e `yf`, tolleranze strette |
| [main_task2.m](hm1_main_task2.md) | salita verticale imposta + arco propulso: stesso shooting 4x4, `H(0)=0` generalizzata, e perche' NON serve una condizione di giunzione |
| [main_task3.m](hm1_main_task3.md) | il coast balistico: shooting 5x5, switching function `S=0`, corner conditions di Weierstrass-Erdmann, la radice spuria |
| [main_task4.m](hm1_main_task4.md) | staging: salto di massa, costati continui, Hamiltoniana **discontinua**, sweep sull'istante di separazione |
| [validate_rocket_sled.m](hm1_validate_rocket_sled.md) | validazione dello shooting contro il rocket sled a energia minima (soluzione analitica nota) |
| [validate_staging_corner.m](hm1_validate_staging_corner.md) | corner condition come quinto residuo, invarianza di gauge di `lam_m`, confronto con lo sweep |

## 3 - HM2 -- powered descent, metodi diretti (`HM2_powered_descent/`)

Dalla collocazione trapezoidale alla convessificazione lossless. La parte piu' moderna.

| file | cosa spiega |
|---|---|
| [ode_descent.m](hm2_ode_descent.md) | RHS punto-materiale 2D con spinta ZOH: derivazione adimensionale, `Vc = V_ref/c`, il kink di `abs(T)` in `T=0` |
| [ode_descent_uacc.m](hm2_ode_descent_uacc.md) | RHS con **accelerazione** `u = T/m` tenuta costante: omogeneita' in massa, il ponte verso il log-mass |
| [rk4_zoh.m](hm2_rk4_zoh.md) | RK4 a 4 stadi su intervallo ZOH: ordine quattro, riga di massa esatta, perche' passo fisso |
| [lti_zoh.m](hm2_lti_zoh.md) | ZOH **esatto** via `expm`: trucco di van Loan a blocchi, `A` nilpotente, forma chiusa |
| [main_task1.m](hm2_main_task1.md) | collocazione trapezoidale: decision vector, **defect constraints**, costo min-fuel, glide slope, `fmincon`/SQP, diagnostiche KKT |
| [main_task2.m](hm2_main_task2.md) | **il file piu' grosso della repo** (1420 righe): quattro trascrizioni ZOH -- RK4 nonlineare, LTV+SCvx (`fmincon` e YALMIP/ECOS), GFOLD log-mass |
| [proto_gfold_logmass.m](hm2_proto_gfold_logmass.md) | **convessificazione lossless**: cambio di variabile log-mass, dinamica esattamente LTI, slack `sigma`, tangente sul bound di spinta |

## 4 - HM3 -- il modello del lanciatore a max-q (`HM3/`)

Pitch-plane, corpo rigido instabile, bending, attuatore. Le fondamenta di tutto il resto.

| file | cosa spiega |
|---|---|
| [load_hw3_params.m](hm3_load_hw3_params.md) | i parametri Greensite a `t = 72 s`: Table 1 vs dataset LPV, derivazione di `mu_alpha`/`mu_c`, `qbar`, scaling del Task 3 |
| [build_plant_rigid.m](hm3_build_plant_rigid.md) | EOM pitch-plane a 4 stati, il polo instabile `+sqrt(A6)`, il polo lento di deriva, **la convenzione di segno di `alpha_w`** |
| [build_plant_full.m](hm3_build_plant_full.md) | il modo di bending forzato dal TVC e **la contaminazione INS** (Eq. 2): da dove escono i +29 dB di risonanza nel loop |
| [build_tvc.m](hm3_build_tvc.md) | attuatore TVC del 2o ordine + ritardo di 20 ms via **Pade** (all-pass a fase non minima) |
| [build_notch_filter.m](hm3_build_notch_filter.md) | notch / lead-lag (Eq. 4): **gain vs phase stabilization**, il prezzo del notch profondo, il detuning |
| [assemble_loop.m](hm3_assemble_loop.md) | chiusura dell'anello PD: `L` per Nichols e `T` per la simulazione; **l'integratore di deriva** che avvelena i margini |
| [load_wind_profile.m](hm3_load_wind_profile.md) | la raffica 1-coseno, la conversione `alpha_w = v_w/V`, perche' max-`qbar` e' la condizione di carico |
| [simulate_gust_response.m](hm3_simulate_gust_response.md) | risposta alla raffica, budget di `alpha`, indicatore di carico `qbar*alpha` |

## 5 - HM3 -- margini e progetto: **il cuore metodologico**

Il lanciatore e' open-loop instabile, quindi l'anello e' **condizionalmente stabile** e
un `margin()` secco mente. E' la storia che devi saper raccontare.

| file | cosa spiega |
|---|---|
| [classify_margins.m](hm3_classify_margins.md) | **perche' `margin()` e' privo di senso** qui, e i margini classificati per banda (aero / rigid / flex / delay) via `allmargin` |
| [design_controller.m](hm3_design_controller.md) | i guadagni canonici `Kp0 = 2*A6/K1`, `Kd0 = sqrt(A6)/K1` **derivati** dai requisiti 6 dB / 30 gradi; auto-tuner `fminsearch` sul loop pieno |
| [plot_nichols_lv.m](hm3_plot_nichols_lv.md) | Nichols col punto critico a `(-180, 0 dB)` da `1+L=0` (convenzione del corso), la curva che "viene dall'alto", marker per banda |
| [make_pz_figures.m](hm3_make_pz_figures.md) | mappe poli-zeri: polo aero instabile, bending quasi sull'asse, zero RHP del Pade |

## 6 - HM3 -- gli entry point (quelli che il professore fa girare)

| file | cosa spiega |
|---|---|
| [main_task1.m](hm3_main_task1.md) | **Task 1**: corpo rigido, attuatore ideale, PD sul loop pieno verso 6 dB / 30 gradi, Nichols, risposta alla raffica |
| [main_task2.m](hm3_main_task2.md) | **Task 2**: TVC + ritardo + bending; il trade fra quattro filtri, il notch profondo, e la **ri-sintonizzazione** del PD |
| [main_task3.m](hm3_main_task3.md) | **Task 3**: robustezza `+/-30%` su `mu_alpha` e `mu_c`, i quattro vertici, il margine aerodinamico come vincolo che stringe |
| [main_montecarlo.m](hm3_main_montecarlo.md) | *(extra)* Monte Carlo N=1500 su cinque incertezze -- **e la sua staleness rispetto ai Task 1-3, verificata** |

## 7 - HM3 -- track Simulink

| file | cosa spiega |
|---|---|
| [init_simulink_hm3.m](hm3_init_simulink.md) | il ponte script -> Simulink: ogni parametro di blocco precalcolato dagli stessi moduli dei `main_task*` |
| [run_simulink_closed_loop.m](hm3_run_simulink_closed_loop.md) | driver di validazione: modello vs baseline analitica degli script, e cosa significherebbe uno scostamento |

## 8 - LTV_FULL_ASCENT -- l'estensione LPV (oltre la traccia)

Il design congelato a max-q sollevato all'intera ascesa (0-140 s), con gain scheduling.

| file | cosa spiega |
|---|---|
| [ode_lpv_ascent.m](lpv_ode_lpv_ascent.md) | RHS LTV rigida: coefficienti interpolati **dentro** il loop di `ode45`, PD chiuso in linea, **la fallacia del frozen-time** |
| [ode_lpv_flex.m](lpv_ode_lpv_flex.md) | RHS LTV flessibile: bending con `omega(t)`, misure INS contaminate, notch variabile |
| [build_hm3_full_ascent.m](lpv_build_full_ascent.md) | authoring del `.slx` rigido da script: **il trucco dei coefficienti effettivi** (l'interpolazione non commuta con la moltiplicazione) |
| [build_hm3_full_ascent_flex.m](lpv_build_full_ascent_flex.md) | modello flessibile da script: notch che insegue `omega(t)`, e la discrepanza tabelle-vs-RHS |
| [main_full_ascent.m](lpv_main_full_ascent.md) | propagazione dell'ascesa completa e il limite del design congelato lontano da max-q |
| [main_flex.m](lpv_main_flex.md) | il gemello flessibile: `omega_BM` **non e' una costante**, e cosa succede al notch fisso |
| [main_q_scheduling.m](lpv_main_q_scheduling.md) | **schedulare in `qbar` invece che nel tempo**: perche' il tempo e' barare |
| [init_simulink_lpv.m](lpv_init_simulink.md) | tabelle dei coefficienti effettivi, gain schedule, generatore di vento, contratto col base workspace |
| [run_full_ascent_simulink.m](lpv_run_full_ascent_simulink.md) | overlay Simulink vs `ode45`: cosa deve coincidere **esattamente** per arrivare a 1e-7 rad |
| [run_flex_simulink.m](lpv_run_flex_simulink.md) | overlay flessibile a 13 stati e il mismatch fra tabelle e RHS dietro il residuo |

## 9 - Test e benchmark

| file | cosa verifica |
|---|---|
| [HM0 falcon9AscentTest.m](hm0_test_falcon9Ascent.md) | il pattern run+harvest contro il `clear` iniziale degli script, e cosa resta **non testato** |
| [HM1 odeBurnTest.m](hm1_test_odeBurn.md) | `ode_burn`: linear tangent law, costato di massa, limiti balistico e di Tsiolkovsky |
| [HM1 odeBurnPerformanceTest.m](hm1_test_odeBurnPerformance.md) | benchmark della RHS hot-loop; **giustifica l'assenza dell'`arguments` block** |
| [HM2 odeDescentTest.m](hm2_test_odeDescent.md) | la RHS del descent: derivata analitica vs differenze finite, segni fisici |
| [HM2 rk4ZohTest.m](hm2_test_rk4Zoh.md) | RK4: verifica dell'**ordine di convergenza** (dimezzo `dt`, l'errore cala di 16x) |
| [HM2 gfoldLogMassTest.m](hm2_test_gfoldLogMass.md) | il cambio di variabile log-mass: velocita' indipendente dalla massa, `z` lineare, ZOH esatta |
| [HM2 descentDynamicsPerformanceTest.m](hm2_test_descentDynamicsPerformance.md) | `matlab.perftest` sull'hot loop -- e perche' chiamarlo *test* di regressione e' generoso |
| [HM3 hm3PlantTest.m](hm3_test_plant.md) | congela la **fisica**: il polo a `+sqrt(A6)`, il bending a `omega_BM`, `zeta_BM` |
| [HM3 hm3FilterTest.m](hm3_test_filter.md) | congela i blocchi di compensazione: profondita' del notch, fase del Pade |
| [HM3 hm3LoopTest.m](hm3_test_loop.md) | congela **le conclusioni**: i margini del Task 1, l'instabilita' senza filtro, il notch che salva l'anello (golden values) |
| [HM3 hm3LoopPerformanceTest.m](hm3_test_loopPerformance.md) | il costo dell'assemblaggio del loop: hot path dell'auto-tuner e del Monte Carlo |

---

## Consiglio di studio

Non leggerli in ordine. I file ad alto rendimento, nell'ordine in cui conviene
affrontarli:

1. **[hm3_build_plant_rigid.md](hm3_build_plant_rigid.md)** -- da qui parte tutto. Se
   sai derivare `theta_ddot = A6*theta + (A6/V)*zdot + K1*delta` e spiegare **perche'
   `A6 > 0` significa instabilita' statica** (centro di pressione davanti al
   baricentro) e perche' il polo cade esattamente in `+sqrt(A6)`, hai in mano meta'
   di HM3.

2. **[hm3_classify_margins.md](hm3_classify_margins.md)** -- il punto metodologico su
   cui l'homework si gioca. Il velivolo e' open-loop instabile, quindi l'anello e'
   **condizionalmente stabile**: la curva di Nichols deve passare *fra* i due punti
   critici e un singolo numero di `margin()` non significa niente. Saper dire questo,
   e perche' i margini vanno classificati per banda, e' la differenza fra 28 e 30.

3. **[hm3_build_notch_filter.md](hm3_build_notch_filter.md)** + **[hm3_main_task2.md](hm3_main_task2.md)**
   -- **gain stabilization vs phase stabilization** del modo di bending, e il prezzo
   che si paga per il notch profondo (devi conoscere `omega_BM` con precisione, e il
   ritardo di fase che aggiunge ti mangia il phase margin rigido, tanto da costringere
   a ri-sintonizzare il PD).

4. **[hm2_proto_gfold_logmass.md](hm2_proto_gfold_logmass.md)** -- la **convessificazione
   lossless**. Il cambio di variabile `z = ln(m)`, `u = T/m`, `sigma = abs(T)/m` che
   rende la dinamica esattamente LTI, e il rilassamento `abs(u) <= sigma` che
   all'ottimo e' attivo. E' il pezzo di teoria piu' elegante di tutta la repo, ed e'
   quello su cui e' piu' facile fare bella figura.

5. **[hm1_ode_burn.md](hm1_ode_burn.md)** -- Hamiltoniana, Eulero-Lagrange, primer
   vector. Poche righe di codice, tutta la teoria del controllo ottimo indiretto.
