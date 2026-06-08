# Domande per il docente — HM1 & HM2

Lista di dubbi / punti da chiarire emersi durante l'implementazione.
Ordinati per priorità all'interno di ciascuna sezione.

---

## HM1 — Ottimizzazione indiretta dell'ascesa

### Metodologia / formulazione

1. **Strategia di continuazione (Task 1).** Per far convergere lo shooting
   sullo sweep in `Q` ho dovuto usare una continuazione (parto da `Q ≈ 1.8`
   e tengo come warm-start la soluzione precedente, propagando avanti e
   indietro). È un approccio accettato per la consegna o ci si aspetta una
   convergenza "cold" tramite buone initial guess analitiche?

2. **Quota della salita verticale (Task 2/3).** Ho fissato `y₁ = 10⁻⁴`
   (≈ 620 m dimensionali). Il valore è suggerito altrove o va dimensionato
   in base ad un criterio fisico (es. clearance dalla rampa, fine effetto
   suolo)?

3. **Arco di coast (Task 3).** Durante il coast `λ_m = cost`, `λ_vy` rampa
   lineare e ho usato una **chiusura analitica balistica**
   (`y(tc) + ½ vy(tc)² = yf`) invece di re-integrare numericamente. È la
   modalità preferita o sarebbe meglio mantenere coerenza propagando tutte
   le fasi con `ode45` per uniformità di errore numerico?

4. **Condizione di switching al cutoff (Task 3).** Sto imponendo
   `S(t_c) = |λ_v|/m − 1/c = 0`. Ci sono casi (parametri "estremi" del
   problema) in cui il cutoff può essere a `S < 0` con vincolo di
   non-negatività della spinta attivo? Vogliamo che il codice gestisca
   anche quel branch o è fuori dallo scope?

5. **Staging senza salita verticale (Task 4).** Il task chiede staging
   con sequenza burn-staging-burn (no vertical climb). È una scelta
   didattica voluta o vogliamo anche valutare staging + climb come
   estensione? Nel caso, lo stesso `t_s` ottimo del caso senza climb
   resta significativo?

6. **`Q` costante tra gli stadi (Task 4).** Nella mia formulazione i due
   stadi hanno lo **stesso** `Q` (e quindi stessa spinta). Vogliamo
   esplorare il caso `Q₁ ≠ Q₂` (parametro aggiuntivo dell'outer loop)?
   Il guadagno di payload (+182 % nel caso `Q=2`, `y_f=0.04`) sembra
   sovrastimare la realtà proprio per via di questa semplificazione +
   assenza di drag.

7. **Modello di payload con coefficiente strutturale.** Sto usando
   `m_payload = m_f · (1+η) − η`. Confermo che è la definizione attesa
   (i.e. payload = utile, escludendo struttura jettisonata)?

### Codice / numerica

8. **Tolleranze `fsolve`.** `FunctionTolerance = StepTolerance = 1e-10`
   in tutti i task. Sono coerenti con le tolleranze ODE (`1e-10/1e-12`)
   o stiamo "tirando" oltre la precisione macchina dell'integrazione?

9. **`λ_m` come incognita anche con `Q` fissato.** Nel single-shooting
   Task 1 le incognite sono `[λ_vx0, λ_vy0, λ_y, λ_m0, t_f]`. Includere
   `λ_m0` quando `Q` è fissato è ridondante o serve per la transversality
   (`H(t_f)=0` + free t_f)?

---

## HM2 — Powered descent

### Formulazione del problema

10. **Tempo di volo `t_f = 38 s`.** Il valore è dato (Tabella 1 della
    consegna). Per la sensitivity sweep ho usato `t_f ∈ {0.95, 1.00, 1.05}·38`.
    Ci aspetta una formulazione free-time (`t_f` come decision variable
    aggiuntiva) o resta fixed-time per tutti i task della HM2?

11. **Glide-slope.** L'ho scritto come coppia lineare
    `±x − tan(θ_max)·y ≤ 0` (convessa, ancorata all'origine).
    Confermo che l'angolo di 60° si misura dalla **verticale** locale e
    che il vertice del cono è la **posizione di landing target (0,0)**?

12. **Vincolo di magnitudine spinta.** `0 ≤ |T| ≤ T_max` è
    intrinsecamente SOC (palla in `T_x, T_y`). Per il task con SCvx
    convesso ho usato `norm(T) ≤ T_max` (cono di Lorentz); per fmincon
    sto imponendo `Tx² + Ty² − T_max² ≤ 0`. Il vincolo `T_min ≤ |T|` è
    non-convesso: nella consegna `T_min = 0`, quindi non si attiva, ma
    se in HM3/estensione `T_min > 0` lo dobbiamo gestire via "lossless
    convexification" (Açikmeşe-Ploen)?

13. **Dinamica della massa.** `ṁ = −|T|/c` è non-smooth in `T = 0`.
    Per fmincon-SQP è un problema accettabile; per la versione SCvx la
    linearizzazione del termine `|T|` introduce un termine non-differenziabile
    in `B(x_ref, u_ref)`. Nel codice ho regolarizzato con
    `‖T‖_reg = √(T² + ε)` (con `ε = 1e-6` nondim). Va bene o c'è una
    formulazione "canonica" usata in classe?

### Task 2 (ZOH, Appendice A) e SCvx

14. **Variante richiesta.** L'Appendice A propone la discretizzazione
    LTV-ZOH costruita tramite ODE ausiliaria su
    `[Φ, B̂, ĉ]`. La consegna richiede *solo* la versione LTV+SCvx o
    accetta anche la **versione ZOH non-lineare** (multiple-shooting con
    RK4) come Task 2? Nel mio codice ho implementato entrambe per
    confronto: trapezoidale (Task 1), ZOH+RK4 e LTV-ZOH+SCvx — la
    convergenza fmincon dei tre risultati è coerente a `1e-4` non-dim.

15. **Trust region SCvx.** Sto usando una box dura sulle variabili di
    stato/controllo (in unità non-dim) centrata sulla reference
    corrente, con scaling adattivo `ρ ∈ [10⁻³, 1]` aggiornato in base al
    rapporto `η = ΔJ_actual / ΔJ_predicted`. È accettabile o è preferito
    lo schema con **virtual control** (slack `ν_k` sui defects) +
    penalità `‖ν‖₁` come da Mao/Açikmeşe?

16. **Convergenza SCvx.** Criterio attuale: `‖Δx‖ < 10⁻³` non-dim.
    Sono ~12-15 iterate. Esiste un benchmark di riferimento sul numero
    di iterate / `m_f` finale atteso? Volevo verificare di non essere
    bloccato in un minimo locale spurio della linearizzazione.

17. **YALMIP + ECOS.** Ho aggiunto come quarto track una variante
    SCvx con il sottoproblema convesso modellato in YALMIP e risolto da
    ECOS (replicando `cvx_sled_class_2026.m`). I risultati combaciano
    con fmincon ma il tempo è significativamente maggiore. È il path
    "preferito" dal corso oppure è equivalente?

### Numerica / convergenza fmincon

18. **Convergenza Task 1.** Con `MaxIterations = 500` e algoritmo `'sqp'`
    una delle tre run di sensitivity (`t_f = 0.95·38`) termina al cap
    con first-order optimality ≈ 5×10⁻⁴. Soglia accettabile per la
    consegna o richiede gradienti analitici / `IPOPT` via CasADi?

19. **Initial guess.** Per fmincon sto usando interpolazione lineare
    stato + hover costante `(0, m₀·g)` per il controllo. Funziona, ma
    in presenza di vincolo di glide-slope attivo all'inizio sembra
    sub-ottimale. Suggerimenti su una guess migliore (es. propagazione
    a thrust nullo + correzione)?

20. **N = 50 nodi.** Convergenza pulita a `N = 50` (Δm_f < 0.5 kg fra
    `N = 50` e `N = 100`). Lasciamo `N = 50` o ci aspetta uno studio di
    grid convergence?

---

## Domande trasversali

21. **Stile di report.** Per HM1 ho usato `subfiles` LaTeX con un
    capitolo per task. Va bene mantenere lo stesso schema per HM2/HM3
    o preferisce un report unico per omogeneità?

22. **Confronto fra metodi.** Per HM2 ho confronto trapezoidale vs.
    ZOH-RK4 vs. LTV+SCvx (e variante YALMIP). Lo presentiamo come
    estensione o lo limito a una sola transcrizione (quella richiesta)
    nel report principale, mettendo le altre in appendice?
