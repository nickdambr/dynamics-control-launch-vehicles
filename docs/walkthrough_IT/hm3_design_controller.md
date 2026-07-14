# HM3/design_controller.m

## Ruolo del file nel progetto

Questa funzione e' l'**auto-tuner** del PD d'assetto di HM3. Riceve l'impianto
(`build_plant_rigid` per il Task 1, `build_plant_full` per i Task 2/3), la catena
attuatore/filtro `Wact` (vuota per il Task 1, `TVC * notch` per il Task 2), e
restituisce la coppia di guadagni `(Kp_th, Kd_th)` che porta i **margini
classificati** ai target d'assegnazione: |GM| ~ 6 dB e |PM| ~ 30 deg. E' chiamata
da `main_task1.m` riga 22, da `main_task2.m` riga 25 (Step A, rigido) e riga 153
(re-tune sul loop pieno), da `main_task3.m` riga 23 e da `main_montecarlo.m`
riga 42; la usano anche `make_pz_figures.m` (righe 40 e 55), `init_simulink_hm3.m`
(riga 34), `LTV_FULL_ASCENT/init_simulink_lpv.m` (righe 72 e 82) e i test.

Il file esiste per risolvere un problema molto concreto. La formula analitica di
libro (D'Antuono, Trotta) per un lanciatore aerodinamicamente instabile da'
in forma chiusa la coppia

    Kp0 = 2*A_6/K_1     Kd0 = sqrt(A_6)/K_1

ma quelle formule sono ricavate sulla **dinamica rotazionale disaccoppiata**,
cioe' su theta_ddot = A_6*theta + K_1*delta. Il loop vero di HM3 non e' quello:
c'e' anche il feedback di deriva laterale (Kp_z*z_m + Kd_z*zdot_m) e c'e'
l'accoppiamento (A_6/V)*zdot dentro l'impianto. Il risultato e' che i guadagni
canonici, applicati al loop pieno, **non centrano piu' i target**: eseguendo il
codice, il seed canonico (Kp0 = 1.4817, Kd0 = 0.4029) sul loop rigido pieno da'
|GM_aero| = 4.20 dB e PM_rigid = 27.33 deg -- e non 6 dB / 30 deg. La docstring
(righe 7-9) dichiara esattamente questo: *"the lateral-drift feedback erodes the
aerodynamic gain margin (the canonical decoupled 6 dB drops to ~4 dB on the full
loop)"*. Verificato.

La strategia adottata e' quindi in due tempi: **(1)** parti dalla coppia
canonica come punto iniziale (riga 55), **(2)** ri-sintonizza con un
`fminsearch` che minimizza lo scarto quadratico dai target, dove i margini sono
letti sul **loop pieno** e **classificati per banda di frequenza** da
`classify_margins`. La classificazione e' obbligatoria, non un vezzo: il loop
e' condizionalmente stabile e `margin()` da solo restituisce un numero senza
significato (vedi la pagina di `classify_margins`).

Dipendenze: `assemble_loop` (costruisce L e T), `classify_margins` (legge i
margini per banda), `isstable` per il verdetto finale. Nessun toolbox oltre al
Control System Toolbox; il tuner usa `fminsearch` di base MATLAB (nessun
Optimization Toolbox richiesto).

---

## Firma, docstring e contratto (righe 1-23)

```matlab
function [K, m] = design_controller(G, Wact, o)
```

- Riga 1: firma. `G` = impianto LTI (4 stati nel Task 1, 6 nel Task 2/3),
  `Wact` = blocco in serie controllore->impianto (attuatore TVC + ritardo +
  notch), `o` = name-value. Restituisce `K` (struct dei 4 guadagni) e `m` (i
  margini classificati, arricchiti con `stable`, `L`, `T`).
- Righe 2-11: la docstring dichiara il metodo -- seed canonico D'Antuono
  Eq. 3.6-3.7, poi re-tune sul loop pieno, margini classificati per banda, e il
  verdetto di stabilita' preso da `isstable(T)` e **non** dal segno di un margine.
  Questa distinzione e' il cuore dell'homework.
- Riga 15: `Wact = []` o `tf(1)` significa **attuatore ideale** (Task 1).

## Blocco `arguments` (righe 25-37)

- Righe 28-29: `Kp_z`, `Kd_z` di default **-1e-3** entrambi. Sono i guadagni di
  deriva laterale: piccoli e **negativi**. Non sono ottimizzati dal tuner --
  restano fissi. Il segno negativo e' il classico *load relief* / *drift-minimum*
  (il lanciatore lascia un po' di deriva pur di scaricare l'incidenza aerodinamica
  a max-qbar). Il canale z/delta ha uno zero a destra (eseguendo il codice:
  `z_m/delta = 20.609*(s+3.553)*(s-3.553) / [s*(s+1.861)*(s-1.816)*(s-0.02909)]`),
  quindi e' **non-minimo di fase**: un feedback di posizione laterale forte
  destabilizzerebbe. Da qui l'ordine di grandezza 1e-3.
  - **Il segno e' l'intenzione, non il risultato.** A 1e-3 l'effetto di load relief
    e' trascurabile e il bilancio di incidenza resta dominato dal PD d'assetto:
    misurando la risposta alla raffica nominale, il picco di
    `|alpha| = |theta + zdot/V - alpha_w|` vale **0.577 deg** contro i **0.390 deg**
    che il vento produrrebbe da solo. Cioe' l'anello **aggrava** il carico invece di
    alleggerirlo -- vedi la domanda d'esame in fondo alla pagina. Va detto cosi',
    senza spacciare l'intento per un risultato.
- Righe 30-31: i target `GM = 6` dB e `PM = 30` deg -- i numeri della traccia.
- Riga 32: **`K0` e' accettato ma non usato** (il commento lo dice: *"accepted,
  unused"*). E' un residuo di una versione precedente in cui il tuner poteva
  essere warm-startato. Conseguenza reale e onesta: **il re-tune del Task 2 non
  riparte dai guadagni del Task 1**, ma sempre dalla coppia canonica. Va detto
  all'orale se qualcuno chiede "e il warm start?".
- Righe 33-35: le tre frequenze di banda passate al classificatore. Default
  `Inf/Inf/NaN` = "nessun modo flessibile", che e' il caso Task 1.
- Riga 38: `[]` -> `tf(1)`, l'alias per attuatore ideale.

## Estrazione di A_6, K_1 e definizione delle bande (righe 40-48)

```matlab
iTh = strcmp(G.StateName, 'theta');
iTd = strcmp(G.StateName, 'thetadot');
iDe = strcmp(G.InputName, 'delta');
A6  = G.A(iTd, iTh);
K1  = G.B(iTd, iDe);
w_drift = 0.3*sqrt(A6);
```

- Righe 41-45: i coefficienti **non** vengono passati come argomento, ma **letti
  per nome** dalle matrici dell'impianto. `A6 = A(thetadot, theta)` e
  `K1 = B(thetadot, delta)` sono per costruzione (vedi `build_plant_rigid`
  righe 11-16) esattamente i coefficienti della riga theta_ddot. Il vantaggio:
  la stessa funzione lavora sull'impianto a 4 stati e su quello a 6 stati senza
  modifiche, e nei corner del Task 3 (dove `load_hw3_params` scala A6 e K1) legge
  automaticamente i valori **scalati**. Nessuna possibilita' di disallineamento
  fra impianto e guadagni.
- Riga 46: `w_drift = 0.3*sqrt(A6)`. E' il confine fra la **banda di deriva** e la
  **banda rigida** usato dal classificatore. La scala naturale del corpo rigido e'
  sqrt(A_6) = 1.839 rad/s (il polo instabile aerodinamico); il lobo di deriva sta
  ben sotto. 0.3*sqrt(A6) = **0.5517 rad/s** al nominale. Nota che il default
  interno di `classify_margins` e' 0.5 fisso: qui viene **sovrascritto** con un
  valore che **scala con l'instabilita' dell'aeroshell**, il che e' essenziale nei
  corner del Task 3 (con mu_alpha = 1.3, w_drift diventa 0.629 rad/s).
- Righe 47-48: le bande in formato cell array, poi espanse con `bands{:}` nelle
  chiamate a `classify_margins`. Task 1: default (nessun bending). Task 2/3: il
  chiamante passa `w_flex = 0.6*wBM = 11.34`, `w_flex_hi = 1.5*wBM = 28.35`,
  `w_bending = wBM = 18.9`.

> **Possibile domanda d'esame** -- perche' A_6 e K_1 sono letti dalle matrici
> dell'impianto invece di essere passati come parametri?
> *Risposta:* perche' cosi' i guadagni canonici di partenza sono garantiti
> coerenti con l'impianto che si sta effettivamente chiudendo. Nei corner del
> Task 3 `load_hw3_params` scala mu_alpha e mu_c del +/-30 %: se A_6 e K_1
> venissero passati a mano si rischierebbe di seminare il tuner con i valori
> nominali su un impianto perturbato. Leggendoli per nome di stato non c'e' modo
> di sbagliare.

## Silenziamento del warning `MarginUnstable` (righe 50-52)

```matlab
warnState   = warning('off', 'Control:analysis:MarginUnstable');
restoreWarn = onCleanup(@() warning(warnState));
```

- Riga 51: MATLAB emette `Control:analysis:MarginUnstable` **a ogni chiamata** di
  `allmargin`/`margin` su un loop open-loop instabile -- ed e' il nostro caso
  (l'impianto ha il polo a +1.82 rad/s). Con `fminsearch` che valuta il costo
  centinaia di volte, il command window verrebbe sommerso.
- Riga 52: `onCleanup` **ripristina** lo stato del warning all'uscita, anche in
  caso di errore. E' testato esplicitamente: `tests/hm3LoopTest.m` righe 88-96
  (`testDesignControllerRestoresWarningState`). Non e' pedanteria: silenziare un
  warning globalmente e non ripristinarlo e' un effetto collaterale che
  inquinerebbe tutte le sessioni successive dell'utente.

## Il re-tune: seed canonico + `fminsearch` (righe 54-57)

```matlab
x0 = log([2*A6/K1, sqrt(A6)/K1]);         % D'Antuono Eq. 3.6-3.7
xo = fminsearch(@cost, x0, ...
                optimset('Display','off','TolX',1e-4, ...
                         'TolFun',1e-3,'MaxFunEvals',400));
```

Qui c'e' tutta la teoria dell'homework compressa in due righe. **Derivo le due
formule canoniche**, perche' e' la domanda d'orale piu' probabile su questo file.

**Dinamica rotazionale disaccoppiata.** Si trascura l'accoppiamento con la deriva
laterale (il termine (A_6/V)*zdot, che pesa 1/V ~ 1e-3) e resta

    theta_ddot = A_6*theta + K_1*delta

Con la legge PD (theta_ref = 0): delta = -Kp*theta - Kd*theta_dot, il loop aperto
e'

    L_theta(s) = K_1*(Kp + Kd*s) / (s^2 - A_6)

e la caratteristica a ciclo chiuso e'

    s^2 + K_1*Kd*s + (K_1*Kp - A_6) = 0

da cui, per confronto con s^2 + 2*zeta*omega_n*s + omega_n^2,

    omega_n = sqrt(K_1*Kp - A_6)      zeta = K_1*Kd / (2*omega_n)

**(a) Kp0 = 2*A_6/K_1 viene dal requisito di 6 dB di gain margin.** Si moltiplica
il loop per un guadagno k e si riapplica Routh alla caratteristica
s^2 + k*K_1*Kd*s + (k*K_1*Kp - A_6): la stabilita' richiede il termine noto
positivo, cioe'

    k > k_min = A_6 / (K_1*Kp)

Attenzione: e' un **limite inferiore** su k. Il margine e' una **riduzione** di
guadagno, non un aumento -- ed e' la ragione per cui in dB il margine aerodinamico
e' **negativo** (cfr. Barrows & Orr, cap. 9: *gain-reduction margin*). In dB:

    GM_aero = 20*log10(A_6/(K_1*Kp))

Imporre |GM_aero| = 6 dB significa k_min = 1/2, cioe' K_1*Kp = 2*A_6:

    Kp0 = 2*A_6/K_1

Lettura equivalente e piu' geometrica: |L_theta(0)| = K_1*Kp/A_6 = 2 = +6.02 dB,
e la fase di L_theta a DC e' **esattamente la fase critica** (il denominatore
s^2 - A_6 vale -A_6 < 0 in continua, e un segno meno e' -180 deg di fase --
ossia +180: stesso punto, fase mod 360). Quindi la curva di Nichols **parte
esattamente sulla fase critica, 6 dB sopra il punto critico** (-180 deg, 0 dB
nella convenzione del corso, da 1 + L = 0; D'Antuono etichetta lo stesso punto
+180 sulla carta rietichettata di +360 -- e' la convenzione che il codice usava
fino a poco fa): abbassa il guadagno di 6 dB e ci finisci dentro. Verifica
numerica eseguendo il codice: `allmargin(L_theta)` restituisce GM = **-6.02 dB
alla frequenza 0**.

**(b) Kd0 = sqrt(A_6)/K_1 viene dai 30 deg di phase margin.** Con Kp = Kp0 si ha
K_1*Kp = 2*A_6 e K_1*Kd = sqrt(A_6). La frequenza di attraversamento a 0 dB
risolve |L_theta(j*w)| = 1, cioe'
K_1*sqrt(Kp^2 + Kd^2*w^2) = w^2 + A_6. Elevando al quadrato:

    w_c^4 + A_6*w_c^2 - 3*A_6^2 = 0   ->   w_c^2 = A_6*(sqrt(13)-1)/2 = 1.3028*A_6

La fase di L_theta e' -180 deg + atan(Kd*w/Kp) (il numeratore ha fase positiva,
il denominatore reale negativo contribuisce la fase critica -180 deg -- +180 se
si preferisce l'altro ramo, e' lo stesso mod 360), quindi il phase margin
rispetto al punto critico a -180 deg e'

    PM = atan(Kd*w_c/Kp) = atan( w_c / (2*sqrt(A_6)) ) = atan( sqrt(1.3028)/2 )
       = 29.71 deg

**che non dipende ne' da A_6 ne' da K_1.** Ecco perche' la coppia canonica e' *la*
coppia canonica: da' 6 dB e ~30 deg su qualunque lanciatore instabile, quali che
siano i suoi coefficienti. Verifica numerica: sul loop disaccoppiato con
Kp0 = 1.4817, Kd0 = 0.4029 si ottiene GM = -6.02 dB @ w = 0 e **PM = 29.71 deg
@ 2.099 rad/s**, con la coppia a ciclo chiuso omega_n = sqrt(A_6) = 1.839 rad/s e
**zeta = 0.5 esatto**.

**Perche' allora il re-tune?** Perche' il loop vero non e' quello disaccoppiato.
Eseguendo il codice sul loop rigido **pieno** (4 stati + feedback di deriva) con
gli stessi guadagni canonici si ottiene |GM_aero| = **4.20 dB** e PM_rigid =
**27.33 deg**. Spegnendo i soli guadagni di deriva (Kp_z = Kd_z = 0) si risale a
5.56 dB / 28.74 deg: quindi l'erosione e' **in parte** dovuta all'accoppiamento
gia' presente nell'impianto (i termini a1, a4, A_6/V) e **in parte maggiore** al
feedback di deriva. Il tuner recupera la differenza spostando i guadagni a
**Kp_th = 1.7845, Kd_th = 0.4433**, che sul loop pieno danno esattamente
|GM_aero| = 6.000 dB @ 0.593 rad/s e PM_rigid = 30.00 deg @ 2.455 rad/s.

- Riga 55: `x0 = log([...])`. La **parametrizzazione logaritmica** e' il modo di
  imporre Kp > 0 e Kd > 0 con un ottimizzatore **non vincolato** come
  `fminsearch`: si ottimizza su x = log(K), e K = exp(x) e' positivo per
  costruzione, qualunque cosa faccia il simplesso. Guadagni negativi
  invertirebbero il segno del feedback e renderebbero il problema privo di senso
  fisico.
- Righe 56-57: `fminsearch` = Nelder-Mead, **derivative-free**. E' la scelta
  giusta qui, non un ripiego: la funzione di costo contiene una
  **classificazione** (quale attraversamento appartiene a quale banda) e due
  penalita' a gradino (righe 86 e 89). Il costo e' quindi discontinuo a tratti e
  il suo gradiente non esiste dove la classificazione cambia ramo: `fsolve` o
  `fmincon` con gradienti alle differenze finite sarebbero fragili. Tolleranze:
  `TolX = 1e-4` sui log-guadagni, `TolFun = 1e-3` sul costo (che ha unita' dB^2 +
  deg^2), `MaxFunEvals = 400`.

> **Possibile domanda d'esame** -- il problema ha 2 incognite e 2 target: perche'
> minimizzare uno scarto quadratico invece di risolvere il sistema con `fsolve`?
> *Risposta:* il problema e' effettivamente quadrato, e infatti il costo viene
> annullato (i margini finali sono 6.000 dB e 30.00 deg, non "vicini a"). Ma la
> mappa (Kp, Kd) -> (GM_aero, PM_rigid) passa attraverso `allmargin` +
> classificazione per banda: e' liscia solo finche' gli attraversamenti restano
> nelle rispettive bande. Formulare come minimizzazione permette di **penalizzare
> in modo continuo** il fallimento (NaN -> 1e6, instabile -> +1e4) invece di far
> divergere un root-finder; e Nelder-Mead non ha bisogno di derivate. E' la
> scelta pragmatica corretta.

## Applicazione dei guadagni e verdetto di stabilita' (righe 59-69)

```matlab
K.Kp_th = exp(xo(1));
K.Kd_th = exp(xo(2));
...
[L, T]   = assemble_loop(G, K, Wact);
L        = minreal(L, 1e-6);
m        = classify_margins(L, bands{:});
m.stable = isstable(T);
```

- Righe 59-62: si esce dai log e si riattaccano i guadagni di deriva **fissi**
  (non ottimizzati).
- Riga 64: `assemble_loop` ricostruisce il loop una volta finale, con i guadagni
  vincitori.
- Riga 65: `minreal(L, 1e-6)` **cancella i modi non minimi** introdotti da
  `connect`/`getLoopTransfer`. Serve: senza cancellazione, `allmargin` lavora su
  una realizzazione con poli/zeri quasi coincidenti e i margini si sporcano.
- Riga 67: **il verdetto di stabilita' e' `isstable(T)`, non un segno di margine.**
  Questa e' la riga concettualmente piu' importante del file. Su un sistema
  condizionalmente stabile il segno del gain margin non dice se il loop chiuso e'
  stabile (il criterio di Nyquist richiede il **conteggio degli avvolgimenti**
  rispetto ai poli a destra dell'impianto, e qui i poli instabili sono **due**:
  +1.8165 e +0.0291 rad/s, oltre all'integratore in 0). L'unico test onesto e'
  guardare i poli del ciclo chiuso.

## Stampa diagnostica (righe 71-76)

- Riga 75: stampa `abs(m.aeroGM_dB)`. **La `abs()` non e' cosmetica**: il margine
  aerodinamico e' *negativo* in dB per costruzione (e' una gain-reduction margin),
  e stamparlo come |GM| = 6.0 dB e' cio' che rende confrontabile il numero con il
  requisito "6 dB" della traccia. Il segno resta pero' nella struct: chi legge
  `m.aeroGM_dB` trova -6.000.

## Funzione annidata `cost` (righe 78-90)

```matlab
function c = cost(x)
    Kt.Kp_th = exp(x(1));  Kt.Kd_th = exp(x(2));
    Kt.Kp_z  = o.Kp_z;     Kt.Kd_z  = o.Kd_z;
    [Lt, Tt] = assemble_loop(G, Kt, Wact);
    Lt = minreal(Lt, 1e-6);
    mt = classify_margins(Lt, bands{:});
    if isnan(mt.aeroGM_dB) || isnan(mt.rigidPM_deg)
        c = 1e6;  return;
    end
    c = (abs(mt.aeroGM_dB) - o.GM)^2 + (mt.rigidPM_deg - o.PM)^2;
    if ~isstable(Tt), c = c + 1e4; end
end
```

- Riga 78: funzione **annidata** (non sotto-funzione): vede `G`, `Wact`, `o` e
  `bands` per chiusura lessicale. E' il motivo per cui `fminsearch(@cost, x0)`
  puo' essere chiamata con un solo argomento.
- Righe 80-84: a ogni valutazione si **riassembla l'intero loop** e si
  **riclassificano** i margini. E' costoso (un `connect` + un `minreal` + un
  `allmargin` per iterazione, fino a 400 volte) ma e' l'unico modo per
  ottimizzare sui margini *veri* del loop pieno e non su una loro approssimazione.
- Righe 85-87: **penalita' dura 1e6** se una delle due bande **non ha piu' un
  attraversamento** (`classify_margins` restituisce NaN). Traduzione fisica: se il
  guadagno e' tale che la curva non attraversa piu' la fase critica in banda
  aerodinamica, il "margine aerodinamico" non esiste come numero -- non e' che sia
  infinito, e' che il concetto e' perso. Il `return` anticipato evita di calcolare
  un costo su NaN (che si propagherebbe e farebbe uscire `fminsearch` con un
  simplesso degenere).
- Riga 88: **la funzione di costo vera e propria.**

      c = (|GM_aero| - 6)^2 + (PM_rigid - 30)^2

  Tre osservazioni oneste:
  1. **Le unita' sono mescolate**: dB^2 sommati a deg^2, con peso implicito 1:1.
     Un errore di 1 dB "costa" quanto un errore di 1 deg. E' arbitrario, ma
     essendo il problema quadrato (2 incognite, 2 target) il minimo e' comunque
     lo zero e la pesatura non cambia la soluzione -- influenza solo il cammino
     del simplesso.
  2. Su `aeroGM_dB` c'e' `abs()` (perche' e' negativo), su `rigidPM_deg` **no**.
     Voluto: un PM negativo (curva dalla parte sbagliata) verrebbe cosi'
     penalizzato pesantemente, ((-8) - 30)^2, e spinto verso +30.
  3. Il tuner insegue il **target**, non lo massimizza: cerca |GM| = 6 dB, non
     |GM| >= 6 dB. Se un guadagno desse 9 dB di margine aerodinamico, il costo
     salirebbe. E' coerente con la traccia (che chiede quei valori) ma va detto:
     non e' un'ottimizzazione di robustezza, e' un *matching* di specifica.
- Riga 89: **penalita' additiva 1e4** se il ciclo chiuso e' instabile. Additiva e
  non sostitutiva: il costo dei margini resta dentro, quindi anche nella regione
  instabile il simplesso ha ancora un gradiente numerico che lo guida verso i
  target. Se fosse `c = 1e4` secco, la regione instabile sarebbe un plateau
  piatto e Nelder-Mead ci si perderebbe.

> **Possibile domanda d'esame** -- perche' la penalita' di instabilita' e'
> *additiva* (`c = c + 1e4`) mentre quella dei NaN e' *sostitutiva* (`c = 1e6`)?
> *Risposta:* nel caso instabile i margini classificati esistono ancora e portano
> informazione utile ("quanto sei lontano dai target"), quindi si somma la
> penalita' e si conserva il gradiente. Nel caso NaN i margini **non esistono**:
> non c'e' nulla da sommare, e l'unica cosa sensata e' un valore piatto e molto
> alto che dice al simplesso "torna indietro". Le due situazioni sono
> qualitativamente diverse e sono trattate diversamente.

---

## Possibili domande d'esame

**D: Da dove escono Kp = 2*A_6/K_1 e Kd = sqrt(A_6)/K_1? Ricavale.**
R: Dalla dinamica rotazionale disaccoppiata theta_ddot = A_6*theta + K_1*delta con
PD. La caratteristica a ciclo chiuso e' s^2 + K_1*Kd*s + (K_1*Kp - A_6). Scalando
il loop di k, Routh richiede k > A_6/(K_1*Kp): e' un limite **inferiore**, cioe'
una gain-*reduction* margin. Imporre 6 dB (k_min = 1/2) da' K_1*Kp = 2*A_6.
Con quel Kp, la frequenza di crossover risolve w_c^4 + A_6*w_c^2 - 3*A_6^2 = 0
cioe' w_c = sqrt(1.3028*A_6), e la fase di L e' -180 deg + atan(Kd*w/Kp) (la
fase critica viene dal denominatore reale negativo), da cui il margine rispetto
al punto critico a -180 deg e'
PM = atan(w_c/(2*sqrt(A_6))) = 29.71 deg **indipendentemente da A_6 e K_1** se
Kd = sqrt(A_6)/K_1. Quindi 6 dB e 30 deg escono automaticamente dalla coppia
canonica. La coppia da' anche omega_n = sqrt(A_6) e zeta = 0.5 esatti a ciclo
chiuso.

**D: Se le formule canoniche danno gia' 6 dB e 30 deg, perche' il codice fa
comunque un `fminsearch`?**
R: Perche' quelle formule valgono sul loop **disaccoppiato**. Il loop vero
include il feedback di deriva laterale (Kp_z, Kd_z = -1e-3) e l'accoppiamento
(A_6/V)*zdot dentro l'impianto. Eseguendo il codice, la coppia canonica sul loop
pieno da' 4.20 dB / 27.33 deg invece di 6 / 30: il margine aerodinamico e' eroso
di quasi 2 dB. Il `fminsearch` ricompra quei 2 dB muovendo i guadagni a
Kp = 1.7845, Kd = 0.4433 (Task 1). Nel Task 2, con TVC + ritardo + notch in serie,
il re-tune sposta soprattutto **Kd** (0.44 -> 0.69) per rimpiazzare la fase
mangiata dall'attuatore e dal ritardo.

**D: Perche' `fminsearch` e non `fmincon` o `fsolve`?**
R: (a) La funzione di costo contiene una classificazione per banda e due penalita'
a gradino, quindi e' discontinua a tratti e senza gradiente definito dove un
attraversamento passa da una banda all'altra; Nelder-Mead e' derivative-free e non
si spezza. (b) La positivita' dei guadagni e' gia' imposta dalla parametrizzazione
`x = log(K)`, quindi non serve un ottimizzatore vincolato. (c) Non serve alcun
toolbox aggiuntivo: `fminsearch` e' base MATLAB, e HM3 dichiara di dipendere dal
solo Control System Toolbox.

**D: Perche' il verdetto finale di stabilita' e' `isstable(T)` e non il segno del
gain margin?**
R: Perche' il lanciatore e' open-loop instabile: eseguendo il codice, il loop
rigido ha **due** poli a parte reale positiva (+1.8165 e +0.0291 rad/s) piu'
l'integratore di deriva in 0. Il criterio di Nyquist su un sistema con P > 0 non
si legge dal segno di un margine ma dal conteggio degli avvolgimenti attorno al
punto critico (-180 deg, 0 dB sulla carta di Nichols). Il loop e'
**condizionalmente stabile**: nel Task 2 (loop pieno con
TVC + ritardo + notch) la banda di guadagni stabili verificata numericamente e'
k in [0.50, 2.39] -- sia ridurre sia aumentare il guadagno destabilizza, e i due
estremi sono esattamente l'Aero GM (-6.0 dB) e il Rigid GM (+7.56 dB): il punto
critico resta incastrato fra i due attraversamenti della fase critica. L'unico
test che non mente e' guardare i poli di T.

**D: Nel Task 1 la banda di guadagni stabili e' davvero a due estremi?**
R: **No, e questo va detto onestamente.** Con l'attuatore *ideale* (Task 1) il
loop e' stabile per k in [0.50, +inf): verificato fino a k = 1e6. La fase del loop
tende a -90 deg alle alte frequenze e **non riattraversa** mai la fase critica,
quindi non c'e' un estremo superiore. Coerentemente, `classify_margins` restituisce
`rigidGM_dB = NaN` nel Task 1 (e il commento nel codice, righe 10-11 di
`classify_margins`, lo dice: *"absent for an ideal actuator"*). La conditional
stability a **due** lati compare solo nel Task 2, quando il ritardo di 20 ms e la
dinamica TVC del 2o ordine aggiungono il ritardo di fase che crea
l'attraversamento superiore a 11.1 rad/s.

**D: I guadagni di deriva Kp_z, Kd_z: perche' negativi, perche' piccoli, e perche'
il tuner non li tocca?**
R: Negativi perche' *intendono* realizzare il *load relief* (drift-minimum): il
veicolo accetta un po' di deriva laterale pur di ridurre l'incidenza aerodinamica

    alpha = theta + zdot/V - alpha_w

e quindi il carico qbar*alpha a max-qbar. (Il **meno** su alpha_w lo impone il
plant: l'Eq. (1) ha la colonna di vento `Bw = [0; -a1*V; 0; -A6]`, quindi alpha e'
l'incidenza rispetto alla velocita' **relativa all'aria**.) Piccoli (1e-3) perche'
il canale z/delta e' non-minimo di fase (zero a +3.553 rad/s): un feedback di
posizione laterale forte destabilizzerebbe. Il tuner non li tocca perche' sono una
**scelta di missione** fissata dalle guidelines della traccia, non un grado di
liberta' di stabilizzazione: si sintonizza il PD d'assetto *dato* il canale di
load relief. Il prezzo lo si paga in margine aerodinamico (i ~2 dB persi) e lo si
recupera con il re-tune.

**Ma va detta la verita' sui numeri**: a 1e-3 quel load relief **non funziona**.
Misurando la risposta alla raffica nominale, il picco di `|alpha|` vale **0.577
deg** mentre il vento da solo produrrebbe solo **0.390 deg**: l'assetto non elide
il vento, **ci si somma**, e `qbar*alpha` sale a 46.8 kPa*deg. La ragione e' che
per tenere `theta_ref = 0` il loop deve opporsi al momento aerodinamico, e cosi'
facendo **becca il muso dentro il vento relativo**. **Un puro attitude-hold e'
load-aggravating, non load-relieving.** E non c'e' effetto banderuola a salvarlo:
con `A_6 > 0` il centro di pressione sta davanti al baricentro, quindi il momento
aerodinamico e' **divergente** -- il veicolo non si allinea al vento, ci diverge
contro, ed e' esattamente per questo che serve il PD. Un load relief vero
richiederebbe guadagni di deriva molto piu' aggressivi (o un accelerometro), pagati
in deriva laterale e in margine.

> *Onesta' storica, se l'esaminatore chiede della versione precedente:* fino a poco
> fa il **post-processing** di `simulate_gust_response` scriveva
> `alpha = theta + zdot/V + alpha_w` (con il `+`), e dava 0.255 deg / 20.7 kPa*deg
> -- facendo sembrare il progetto load-relieving. Era un bug **di
> post-processing**, non di modello: il plant integrato (`build_plant_rigid`,
> `ode_lpv_ascent`) ha sempre avuto il meno, quindi `theta`, `z`, `delta` e **tutti
> i margini** non ne erano toccati; a sbagliare era solo il canale diagnostico
> `alpha`. Corretto in tutti i siti.
