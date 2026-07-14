# HM3/tests/hm3LoopTest.m

## Ruolo del file nel progetto

E' la suite di **regressione sulle conclusioni dell'homework**. Mentre
`hm3PlantTest` congela la fisica e `hm3FilterTest` congela i blocchi di
compensazione, questo file congela i **risultati**: che il PD del Task 1
raggiunga i target di margine assegnati, che il modello completo senza filtro sia
**instabile** per colpa del bending, e che il notch profondo lo **stabilizzi in
guadagno** (Task 2). Sono esattamente le tre affermazioni che l'homework va a
difendere.

Copre quattro file: `load_wind_profile.m`, `assemble_loop.m`,
`design_controller.m` e `simulate_gust_response.m`. La struttura logica dei test
segue l'ordine narrativo del report:

1. **Il disturbo** (test 1-4): la raffica 1-coseno ha la forma giusta e
   l'ampiezza giusta.
2. **Il progetto del Task 1** (test 5-7): i guadagni riproducono |GM_aero| = 6 dB
   e PM_rigido = 30 deg, e l'auto-tuner li ritrova da solo.
3. **Il problema del Task 2** (test 8): il modello completo nudo e' instabile.
4. **La soluzione del Task 2** (test 9): il notch lo stabilizza.
5. **Il post-processing** (test 10): il bilancio di incidenza aerodinamica.

**Attenzione, questo e' il file dove vivono i golden value.** Diversi test
confrontano contro numeri hard-coded (`Kref`, 6.0 dB, 30.0 deg, 6.38 m/s) che
**non hanno derivazione analitica**: sono l'output di un `fminsearch` su un
progetto specifico. Vanno capiti per quello che sono -- vedi la sezione dedicata
piu' sotto.

---

## `Kref` -- i golden value del progetto (righe 7-13)

```matlab
properties (Constant)
    % Task-1 PD design: gains re-tuned on the FULL loop so the classified
    % Aero GM = 6 dB / Rigid PM = 30 deg (canonical D'Antuono start, then the
    % lateral-drift aero-GM erosion is compensated). Pinned for regression.
    Kref = struct('Kp_th', 1.7845, 'Kd_th', 0.4433, ...
                  'Kp_z', -1e-3,   'Kd_z', -1e-3)
end
```

- Righe 11-12: quattro numeri, usati da **cinque** dei dieci test. Sono i guadagni
  del **progetto del Task 1**, e coincidono con i valori del README (Kp = 1.78,
  Kd = 0.44). Verificato in MATLAB: `design_controller` restituisce oggi
  Kp_th = 1.784455, Kd_th = 0.443323 (entrambi arrotondano ai valori pinnati).
- `Kp_z = Kd_z = -1e-3` sono i guadagni di **deriva laterale**, fissati (non
  ottimizzati) come prescrive la traccia. Sono **negativi e piccoli**: e' il
  cosiddetto *load relief*. Un guadagno positivo su `z` cercherebbe di riportare
  il lanciatore sulla traiettoria nominale contro il vento, aumentando l'incidenza
  e quindi il carico qbar*alpha; un guadagno **negativo** lascia che il veicolo
  derivi *con* il vento, riducendo l'incidenza. E' il classico compromesso
  "seguire la traiettoria" contro "sopravvivere al carico", e a max-qbar vince il
  secondo.
  - **Ma attenzione a non sovravendere il load relief.** Quella e' l'*intenzione*
    del segno; l'*effetto*, a 1e-3, e' trascurabile. Misurando il bilancio di
    incidenza sull'anello effettivo (vedi
    `testGustResponseAngleOfAttackBudget` piu' sotto) il picco di `|alpha|` vale
    **0.577 deg** contro i **0.390 deg** della raffica da sola: l'anello chiuso
    **aggrava** il carico invece di alleggerirlo. Il PD d'assetto e' dominante e
    becca il muso dentro il vento relativo; i -1e-3 sulla deriva non bastano a
    invertire il bilancio. Va detto cosi', senza abbellirlo.

**Da dove vengono Kp = 1.78 e Kd = 0.44 (la derivazione che serve all'orale).**
Non sono numeri magici: `design_controller.m` (riga 55) parte da una coppia in
**forma chiusa** e poi la raffina. Congelando la deriva laterale, la dinamica di
beccheggio e' `theta_ddot = A_6*theta + K_1*delta`; con la legge PD
`delta = -Kp*theta - Kd*theta_dot` si ottiene

    theta_ddot + K_1*Kd*theta_dot + (K_1*Kp - A_6)*theta = 0

da cui, per confronto con `s^2 + 2*zeta*omega_n*s + omega_n^2`:

    omega_n = sqrt(K_1*Kp - A_6)          zeta = K_1*Kd / (2*omega_n)

**Il termine `-A_6` e' tutto il problema**: l'instabilita' aerodinamica
**sottrae** rigidezza all'anello chiuso. Serve `K_1*Kp > A_6` **solo per non
divergere**. Le gains canoniche di D'Antuono (Eq. 3.6-3.7, citate alla riga 55 del
sorgente) sono

    Kp0 = 2*A_6/K_1        Kd0 = sqrt(A_6)/K_1

e sostituendo si vede subito cosa fanno:

    omega_n = sqrt(2*A_6 - A_6) = sqrt(A_6) = 1.839 rad/s
    zeta    = K_1*(sqrt(A_6)/K_1) / (2*sqrt(A_6)) = 1/2 = 0.5

cioe' **piazzano la coppia di beccheggio in anello chiuso esattamente alla
frequenza del polo instabile, con smorzamento 0.5**. E' una scelta elegante: la
banda di controllo e' fissata dalla velocita' dell'instabilita' che deve
combattere. Numericamente: Kp0 = 1.4816, Kd0 = 0.4029.

Da li' `fminsearch` raffina **sull'anello completo** (che include la deriva
laterale, che erode il margine aerodinamico) e arriva a Kp = 1.7845, Kd = 0.4433.
La coppia equivalente di beccheggio risultante e'

    omega_n = sqrt(4.5647*1.7845 - 3.3818) = 2.18 rad/s      zeta = 0.46

che e' il "equivalent pitch pair omega_c = 2.18, zeta = 0.46" del README
(verificato). Restare in banda 1-4 rad/s e' il criterio tipico del corso.

> **Possibile domanda d'esame** -- Perche' i guadagni canonici non bastano e devi
> raffinarli?
> *Risposta:* Perche' sono derivati sulla dinamica di beccheggio **disaccoppiata**,
> mentre l'anello reale include la deriva laterale (i termini `a_1`, `a_4`, `A_6/V`
> della matrice A) e i guadagni `Kp_z`, `Kd_z`. Il feedback di deriva **erode il
> margine di guadagno aerodinamico**: come dice il commento di `design_controller.m`
> (righe 8-10), i 6 dB canonici della coppia disaccoppiata scendono a ~4 dB
> sull'anello completo. Il `fminsearch` li ricompensa alzando Kp da 1.48 a 1.78 e
> Kd da 0.40 a 0.44, e riporta la coppia classificata a 6 dB / 30 deg.

---

## Setup (righe 15-35)

```matlab
methods (TestMethodSetup)
    function muteConditionallyStableMarginWarning(testCase)
        ws = warning('off', 'Control:analysis:MarginUnstable');
        testCase.addTeardown(@() warning(ws));
    end
end
```

- Righe 15-28 (il resto del setup, non riportato sopra): la proprieta' `p`
  (righe 15-17) e il blocco `TestClassSetup` (righe 19-28), eseguito **una volta
  per classe**. `addHm3ToPath` aggiunge `HM3/` al path con una `PathFixture`
  (che il framework rimuove da sola a fine suite, anche in caso di errore), e
  `loadNominalParams` riempie `testCase.p` con `load_hw3_params()`. E' da li' che
  arriva il `pp = testCase.p` usato da tutti i test.
- Righe 30-35: **questo blocchetto e' un indizio fisico**, non un dettaglio
  igienico. `Control:analysis:MarginUnstable` e' l'avviso che la Control System
  Toolbox emette quando le si chiede un margine su un anello aperto **instabile**.
  Il nostro anello aperto e' instabile per costruzione (due poli a parte reale
  positiva, vedi `hm3PlantTest`), quindi la warning arriverebbe a ogni singola
  chiamata di `margin`/`allmargin`, riempiendo l'output dei test. Va mutata, ma
  con `addTeardown`, che **ripristina lo stato precedente anche se il test
  fallisce**.
- La warning e' il sintomo della **stabilita' condizionata**: il sistema e' stabile
  in anello chiuso solo se il guadagno sta **dentro una banda**, non semplicemente
  "sotto un massimo". Abbassando troppo il guadagno l'anello torna instabile,
  perche' smette di controbilanciare il polo aerodinamico. E' il motivo per cui
  esiste un **gain margin di riduzione** (aeroGM, negativo in dB) oltre al
  classico gain margin di aumento.

---

## `testGustProfileShape` (righe 38-47)

```matlab
w  = load_wind_profile(pp, 'Vg', 8.0, 'Tg', 3.0, 't0', 1.0);
testCase.verifyEqual(max(abs(w.vw(w.t < 1.0))), 0, 'AbsTol', 1e-15);
[vwPeak, iPeak] = max(w.vw);
testCase.verifyEqual(vwPeak, 8.0, 'AbsTol', 1e-6);
testCase.verifyEqual(w.t(iPeak), 1.0 + 1.5, 'AbsTol', w.t(2)-w.t(1));
testCase.verifyEqual(w.alphaw, w.vw/pp.V, 'AbsTol', 1e-15);
```

- Righe 42-46: quattro invarianti della raffica **1-coseno** (`load_wind_profile.m`
  riga 56):

      v_w(t) = 0.5*V_g * (1 - cos(2*pi*(t - t0)/T_g)),   t0 <= t <= t0 + T_g
      v_w(t) = 0                                          altrove

  - **Causalita'** (riga 42): `v_w = 0` prima dell'istante di innesco `t0`. Se il
    vento fosse non nullo a t = 0, il sistema partirebbe da una condizione forzata
    e la risposta al transitorio si mescolerebbe con quella alla raffica.
  - **Ampiezza di picco** (riga 44): il massimo vale esattamente `V_g`. Deriva dal
    fatto che `1 - cos(x)` ha massimo 2, quindi `0.5*V_g*2 = V_g`. Il fattore 0.5
    e' quello che rende `V_g` interpretabile come *velocita' di picco* e non come
    ampiezza da raddoppiare -- un classico errore di fattore 2.
  - **Istante di picco** (riga 45): a meta' raffica, `t0 + T_g/2 = 2.5 s`. La
    tolleranza e' `w.t(2)-w.t(1)`, cioe' **un passo di campionamento** (dt = 5 ms):
    e' la scelta corretta, perche' il picco viene individuato con `max` su una
    griglia discreta e la sua posizione e' nota solo a meno di un campione.
  - **Definizione di alpha_w** (riga 46): `alpha_w = v_w / V`. E' la
    linearizzazione a piccoli angoli dell'angolo di incidenza indotto da una
    componente di vento **laterale** `v_w` su un veicolo che avanza a `V`:
    `alpha_w = atan(v_w/V) ~ v_w/V` per `v_w << V`. Con V_g = 6.38 m/s e
    V = 937.7 m/s, `alpha_w = 6.8e-3 rad = 0.39 deg`: siamo profondamente in
    regime lineare, l'approssimazione e' irreprensibile.

---

## `testStepProfileShape` (righe 49-53) e `testWindProfileRejectsUnknownProfile` (righe 61-65)

- Righe 51-52: il profilo `'step'` e' zero prima di `t0` e vale `V_g` a fine
  orizzonte. E' il disturbo "peggiore" in senso di persistenza (non rientra mai),
  utile come controllo diagnostico.
- Righe 62-64: un profilo sconosciuto (`'sinusoid'`) solleva
  `load_wind_profile:profile` (riga 64 del sorgente). Come in
  `build_plant_full`, si testa l'**identificatore**, non il messaggio.

---

## `testDefaultGustAmplitudeFromDrywind` (righe 55-59)

```matlab
% Severe dry-wind dispersion at 15.1 km is the documented default
w = load_wind_profile(testCase.p);
testCase.verifyEqual(w.Vg, 6.38, 'AbsTol', 0.05);
```

- Riga 58: **golden value numero uno**. `6.38 m/s` non e' una scelta di progetto:
  e' il risultato dell'interpolazione della dispersione di vento `severe` del file
  `General/hw3-v3/drywind.mat` alla quota `Alt/1000 = 15.143 km`
  (`load_wind_profile.m` righe 39-48). Verificato: `Vg = 6.3785`, e il picco di
  `alpha_w` che ne deriva vale 0.3897 deg (il "0.39 deg" del README).
- **Questo test e' anche, implicitamente, un test di presenza del data file.** Se
  `drywind.mat` mancasse, `load_wind_profile` cadrebbe sul fallback `Vg = 8.0`
  (riga 47) **senza errori**, e il test fallirebbe (8.0 contro 6.38 +/- 0.05). E' un
  effetto collaterale desiderabile: il fallback silenzioso e' pericoloso, perche'
  produrrebbe un'analisi di raffica del 25 % piu' severa senza avvisare nessuno.
  Il test lo intercetta. Vale la pena saperlo dire all'orale.
- La `AbsTol` di 0.05 e' larga rispetto al valore (0.8 %): assorbe eventuali
  differenze nella griglia di interpolazione ma non un cambio di severita'
  (`light`/`moderate` darebbero valori ben diversi).

---

## `testRigidLoopMeetsMarginTargets` (righe 67-76)

```matlab
G = build_plant_rigid(testCase.p);
[L, T] = assemble_loop(G, testCase.Kref);
mm = classify_margins(minreal(L, 1e-6), ...
        'w_drift', 0.3*sqrt(testCase.p.A6));
testCase.verifyTrue(isstable(T));
testCase.verifyEqual(abs(mm.aeroGM_dB), 6.0,  'AbsTol', 0.3);
testCase.verifyEqual(mm.rigidPM_deg,    30.0, 'AbsTol', 0.7);
```

E' il test centrale del Task 1, e il piu' denso di contenuto teorico.

**Riga 72 -- perche' `w_drift = 0.3*sqrt(A_6)`.** E' il confine fra la banda
"deriva laterale" e la banda "corpo rigido" passato a `classify_margins`. Non e'
una costante arbitraria: e' **scalato sulla frequenza dell'instabilita'**. La
logica e': il crossover di controllo deve stare all'ordine di `sqrt(A_6)`
(altrimenti non batte il polo instabile), mentre i fenomeni di deriva laterale
sono molto piu' lenti (l'integratore di posizione `z`). Prendere il 30 % di
`sqrt(A_6)` mette la soglia a `0.3*1.839 = 0.552 rad/s`, comodamente sopra le
crossings di deriva e sotto il crossover rigido. Ed e' esattamente lo stesso
valore che `design_controller.m` (riga 46) calcola internamente: i due devono
coincidere, altrimenti il test misurerebbe una cosa e il tuner ne ottimizzerebbe
un'altra.

**Perche' `classify_margins` e non `margin()`.** Questo e' il punto piu'
importante da saper difendere. L'anello aperto ha **due poli instabili** e un
integratore, quindi la curva di Nichols passa **fra** i due punti critici invece
che sotto (stabilita' condizionata) e ci sono **piu' attraversamenti** dello 0 dB
e dei -180 deg. Verificato con `allmargin` sull'anello progettato:

    attraversamenti di fase (0 dB crossings):
        -133.12 deg  @  0.161 rad/s      <- artefatto di deriva
         -40.53 deg  @  0.222 rad/s      <- artefatto di deriva
         +30.00 deg  @  2.455 rad/s      <- IL margine di fase rigido

    attraversamento di guadagno (-180 deg crossing):
          -6.00 dB   @  0.593 rad/s      <- margine di guadagno AERODINAMICO

Le due crossings a 0.161 e 0.222 rad/s hanno margini di fase **negativi**, ma il
sistema **e' stabile** (`isstable(T) = 1`, `allmargin.Stable = 1`): sono artefatti
del lobo di deriva, non indicatori di instabilita'. `classify_margins` le esclude
con la maschera `pf > w_drift` (riga 49 del sorgente) e poi seleziona con la
strategia `'maxv'`. `mm.drift_w` le conserva comunque, etichettate come tali.

Il margine di guadagno e' **negativo** (-6 dB) ed e' una **riduzione**: dice che
l'anello si destabilizza se il guadagno **scende** di 6 dB. E' la firma della
stabilita' condizionata, e fisicamente e' il *margine aerodinamico*: se il
controllo perde meta' della sua autorita', non riesce piu' a contrastare il polo
a +1.84 rad/s. `classify_margins` lo isola con la maschera `gmdb < 0` (riga 44).

Nota di onesta': su **questo specifico** progetto, `margin()` grezza restituisce
per caso gli stessi numeri (Gm = 0.5012 = -6.00 dB @ 0.593, Pm = 30.00 deg @
2.455) -- verificato. Il commento di `classify_margins.m` (righe 13-15) sostiene
che `margin()` "prenderebbe una delle crossings di deriva", e su questo punto e'
troppo assertivo. Ma la robustezza del ragionamento non cambia: la scelta di
`margin()` **dipende dal progetto**. Con un'altra coppia di guadagni -- per esempio
`Kp = 1.98, Kd = 1.40`, che sono i guadagni **pre-ritaratura** rimasti a lungo
pinnati per sbaglio nel perf-test (vedi la pagina su `hm3LoopPerformanceTest`) --
`margin()` restituisce `Pm = -30.00 deg` su un anello che **e' stabile**:
verificato. Cioe' `margin()` e' inaffidabile su questo tipo di anello e a volte da'
la risposta giusta per fortuna. **`classify_margins` esiste per non dipendere dalla
fortuna**, e questo va detto cosi'.

**Righe 74-75 -- i target 6 dB / 30 deg.** Sono i valori richiesti dalla traccia,
non un'ottimizzazione. Le tolleranze (0.3 dB, 0.7 deg) sono **strette**: sono
tolleranze di **regressione**, non di ingegneria. Servono a far fallire il test se
cambia una qualunque cosa a monte -- il seed del tuner, `TolX`/`TolFun` di
`fminsearch`, la logica di banda del classificatore, i coefficienti letti
dall'LPV. Non dicono "il progetto e' buono entro 0.3 dB": dicono "il progetto e'
**questo**".

**Riga 73 -- `isstable(T)`.** E' il verdetto vero, ed e' logicamente **prioritario**
rispetto ai margini: su un anello condizionatamente stabile, l'unica affermazione
non ambigua e' "tutti i poli dell'anello chiuso hanno parte reale negativa". I
margini sono **metriche di robustezza** rispetto a quel verdetto, non sostituti.

---

## `testDesignControllerMeetsTargets` (righe 78-86)

```matlab
[K, m] = design_controller(G, [], 'verbose', false);
testCase.verifyEqual(abs(m.aeroGM_dB), 6.0,  'AbsTol', 0.3);
testCase.verifyEqual(m.rigidPM_deg,    30.0, 'AbsTol', 0.7);
testCase.verifyTrue(m.stable);
testCase.verifyEqual(K.Kp_th, testCase.Kref.Kp_th, 'AbsTol', 1e-2);
testCase.verifyEqual(K.Kd_th, testCase.Kref.Kd_th, 'AbsTol', 1e-2);
```

- E' il test **piu' forte** della suite, e il piu' fragile. Le righe 81-83
  verificano che il tuner **raggiunga i target**; le righe 84-85 verificano che ci
  arrivi **passando esattamente per Kref**.
- La differenza e' sostanziale. Le prime tre asserzioni testano una **proprieta'**
  (l'ottimizzatore converge dove deve). Le ultime due pinnano il **percorso**: dato
  che il costo (righe 78-90 di `design_controller.m`) e'

      c = (|aeroGM| - 6)^2 + (rigidPM - 30)^2   [+ 1e4 se instabile]

  e' un problema a 2 incognite e 2 obiettivi: la soluzione e' (localmente) unica e
  `Kref` e' quella soluzione. Pinnarla a 1e-2 significa pinnare **il seed
  D'Antuono, la tolleranza di `fminsearch`, e la definizione delle bande** tutti
  insieme. Se domani si cambiasse il seed (per esempio partendo da `[1 1]`),
  `fminsearch` potrebbe convergere a un minimo locale diverso che soddisfa
  comunque i target: le prime tre asserzioni passerebbero, le ultime due no. Il
  test e' progettato per accorgersene, ed e' la scelta giusta in un contesto di
  regressione.
- Nota sulla parametrizzazione: il tuner ottimizza in `log` (riga 55:
  `x0 = log([2*A6/K1, sqrt(A6)/K1])`). E' un dettaglio importante -- garantisce
  che i guadagni restino **positivi** durante la ricerca, senza vincoli espliciti,
  e rende la ricerca naturalmente **moltiplicativa** (un passo di `fminsearch` e'
  una variazione percentuale, non assoluta), che e' la scala giusta per dei
  guadagni. `fminsearch` e' Nelder-Mead, senza vincoli: senza il `log` potrebbe
  proporre guadagni negativi e produrre anelli assurdi.

---

## `testDesignControllerRestoresWarningState` (righe 88-96)

```matlab
warning('on', 'Control:analysis:MarginUnstable');
G = build_plant_rigid(testCase.p);
design_controller(G, [], 'verbose', false);
st = warning('query', 'Control:analysis:MarginUnstable');
testCase.verifyEqual(st.state, 'on');
```

- E' un **test di regressione su un bug di igiene**, non sulla fisica. Il commento
  (righe 89-90) e' esplicito: `design_controller` muta la warning al proprio
  interno (righe 51-52 del sorgente) e **deve ripristinare lo stato del chiamante**
  all'uscita.
- Perche' e' importante e non pedanteria: `warning('off', ...)` e' uno **stato
  globale della sessione MATLAB**. Se `design_controller` lo lasciasse spento, ogni
  successiva chiamata a `margin()` **in tutta la sessione dell'utente** smetterebbe
  di avvisare che l'anello e' instabile. Su un progetto dove la stabilita'
  condizionata e' il fenomeno centrale, sopprimere silenziosamente quell'avviso e'
  esattamente il tipo di errore che porta a conclusioni sbagliate.
- Il sorgente lo fa correttamente con `onCleanup` (riga 52), che garantisce il
  ripristino anche se `fminsearch` solleva un'eccezione. Questo test verifica che
  il pattern funzioni davvero.
- Sottigliezza: il test **si affida** al fatto che `TestMethodSetup` (righe 30-35)
  spenga la warning e la riaccenda in teardown. La riga 91 la riaccende
  esplicitamente *dentro* il test, perche' altrimenti il setup l'avrebbe appena
  spenta e il test verificherebbe `'off' == 'on'`.

---

## `testBareFullModelIsBendingUnstable` (righe 98-104)

```matlab
% Task 2, Step B: TVC + delay with no bending filter -> unstable
Gf = build_plant_full(pp, 'ins');
[~, T] = assemble_loop(Gf, testCase.Kref, build_tvc(pp));
testCase.verifyFalse(isstable(T));
```

- Riga 103: **un `verifyFalse` che vale come mezzo homework**. E' la
  *dimostrazione eseguibile* che il Task 2 e' un problema vero: si prende il
  progetto del Task 1 (stesso `Kref`), lo si mette sul modello a 6 stati con INS
  contaminato e con il TVC reale, **senza filtro**, e l'anello chiuso e'
  **instabile**. Verificato: `isstable(T) = 0`.
- La causa e' il picco di risonanza di **+29 dB** a omega_BM = 18.9 rad/s (README).
  Il modo di bending, di per se', e' stabile (zeta_BM = 0.005 > 0): e' la
  **chiusura dell'anello** attraverso il percorso `delta -> eta` (forzamento TVC)
  e `eta -> theta_m` (contaminazione INS) a renderlo instabile, come spiegato
  nella pagina su `hm3PlantTest`.
- **Un test negativo di questo tipo e' prezioso** e spesso trascurato: senza di
  esso, si potrebbe aggiungere il notch "per sicurezza" e non avere alcuna prova
  che serva. Qui la prova c'e', ed e' automatica.

---

## `testDeepNotchStabilisesFullModel` (righe 106-114)

```matlab
Gf = build_plant_full(pp, 'ins');
Hn = build_notch_filter(pp.wBM, 0.002, 0.7, +1);
[L, T] = assemble_loop(Gf, testCase.Kref, build_tvc(pp)*Hn);
testCase.verifyTrue(isstable(T));
testCase.verifyLessThan(20*log10(abs(freqresp(L, pp.wBM))), -10);
```

- Riga 110: il **notch ritenuto** del Task 2, con la parametrizzazione esatta del
  progetto: centrato su `omega_BM`, `zeta_N = 0.002`, `zeta_D = 0.7`, **a fase
  minima** (`+1`). Profondita' -50.9 dB (vedi la pagina su `hm3FilterTest`).
- Riga 112: l'anello ora **e' stabile**. Verificato.
- Riga 113: **questa e' l'asserzione con contenuto ingegneristico**. Non basta che
  l'anello sia stabile: si verifica che lo sia per **gain stabilisation**, cioe'
  che il guadagno d'anello alla frequenza di bending sia sceso **sotto -10 dB**.
  Verificato: `|L(j*omega_BM)| = -21.87 dB` con questi guadagni.
  - La soglia -10 dB e' un requisito ingegneristico, non un numero preciso: la
    letteratura chiede tipicamente **>= 12 dB** di attenuazione per un modo
    gain-stabilised (il README lo dice). Il test usa 10 dB come soglia
    conservativa, cioe' verifica una condizione **piu' debole** di quella
    dichiarata nel report. Va bene per un test di regressione (non deve essere
    fragile), ma bisogna sapere che il requisito vero e' piu' stringente.
  - Perche' e' l'invariante giusto: `isstable` da solo non distingue una
    stabilizzazione **in guadagno** (robusta: il lobo e' cosi' in basso che la sua
    fase e' irrilevante) da una **in fase** (fragile: dipende dal passare dalla
    parte giusta del punto critico). Aggiungere la soglia sul modulo pinna **quale
    strategia** e' stata usata, non solo che ha funzionato.

**Nota di onesta' importante -- questo test NON pinna il progetto finale.** Usa
`Kref`, cioe' i guadagni del **Task 1** (1.7845 / 0.4433), non i guadagni
**ri-tarati** del Task 2 che il README dichiara come progetto consegnato
(`Kp = 1.73`, `Kd = 0.69`). La configurazione testata e' quella intermedia della
tabella del report ("deep notch (retained)" con i guadagni del Task 1: aero GM
6.08 dB, **rigid PM 14.6 deg**, |L(omega_BM)| = -21.9 dB), cioe' proprio la
configurazione con il **margine di fase collassato** che ha costretto al
re-tuning. **Non esiste alcun test che pinni i guadagni finali del Task 2**, ne'
i suoi margini (30 deg @ 3.2 rad/s, DM 165 ms, |L(omega_BM)| = -18 dB). E' la
lacuna di copertura piu' significativa della suite.

---

## `testGustResponseAngleOfAttackBudget` (righe 116-128)

```matlab
% alpha = theta + zdot/V - alpha_w (the plant's own convention:
% Eq. (1) has Bw = [0; -a1*V; 0; -A6]), and peaks match the histories
r  = simulate_gust_response(T, w);
testCase.verifyEqual(r.alpha, r.theta + r.zdot/pp.V - r.alphaw, ...
    'AbsTol', 1e-15);
testCase.verifyEqual(r.peak_theta, max(abs(r.theta)), 'AbsTol', 1e-15);
testCase.verifyEqual(r.theta(1), 0, 'AbsTol', 1e-15);
```

**Le tre asserzioni sono, nell'ordine, debolissime.** Vanno esaminate con onesta'.

- Righe 124-125: verifica che `r.alpha == r.theta + r.zdot/V - r.alphaw`. Ma quella
  **e' letteralmente la riga 29 di `simulate_gust_response.m`**. Il test
  **ricopia la formula** e verifica che la formula sia se stessa: e' una
  tautologia. Il suo unico valore e' da *change-detector* (se qualcuno riscrive la
  riga 29, il test fallisce e lo costringe a giustificarsi).
- Riga 126: `peak_theta == max(abs(theta))` -- di nuovo, e' la riga 31 del
  sorgente riscritta. Tautologia.
- Riga 127: `theta(1) == 0`. Sembra un test di condizione iniziale nulla, ma passa
  **automaticamente**: `lsim` parte da stato nullo per default, `D = 0` (nessun
  feedthrough) e in ogni caso `alphaw(1) = 0` perche' la raffica comincia a
  t0 = 1 s. Non puo' fallire se non con un errore grossolano.

### Perche' il segno e' il MENO (il punto da difendere all'orale)

Il segno di `alpha_w` non e' una convenzione libera: **lo detta il plant**. Guardiamo
la colonna di disturbo di `build_plant_rigid.m` (riga 17), che ricopia l'Eq. (1)
della traccia:

      Bw = [0; -a1*V; 0; -A6]

La riga di beccheggio del plant e' quindi

      theta_ddot = A_6*theta + (A_6/V)*z_dot + K_1*delta - A_6*alpha_w
                 = A_6 * ( theta + z_dot/V - alpha_w ) + K_1*delta

e la riga laterale, sviluppata, da'

      z_ddot = a_1*V*( theta + z_dot/V - alpha_w ) + a_4*theta + a_3*delta

**Entrambe le righe forzano sulla stessa incidenza aerodinamica**

      alpha = theta + z_dot/V  -  alpha_w        (segno MENO)

Il significato fisico e' preciso: `alpha` e' l'incidenza rispetto alla velocita'
**relativa all'aria**, non rispetto al suolo. `theta + z_dot/V` e' l'orientamento
del veicolo rispetto alla traiettoria inerziale; `alpha_w = v_w/V` e' la rotazione
del vettore velocita' relativa prodotta dal vento laterale. Un vento che spinge
nella stessa direzione in cui il muso e' gia' inclinato **riduce** l'incidenza
vista dal veicolo: da qui la sottrazione. Un `+` descriverebbe un plant diverso da
quello che si sta integrando -- sarebbe in contraddizione con l'Eq. (1) stessa.

Il codice usa oggi il meno **in tutti i siti**: `simulate_gust_response.m` riga 29,
questo test alla riga 124, e il post-processing LPV di
`LTV_FULL_ASCENT/main_full_ascent.m` riga 153. Il plant (sia `build_plant_rigid`,
sia `ode_lpv_ascent`) lo ha sempre avuto.

> **Nota onesta, da raccontare com'e' andata.** Fino a poco fa il **post-processing**
> aveva un `+`: `r.alpha = theta + zdot/V + alpha_w`. Era un bug di
> post-processing, non di modello -- il plant integrato era ed e' sempre stato
> quello giusto, quindi `theta`, `z` e `delta` non ne erano toccati; a sbagliare
> era **solo** il canale diagnostico `alpha` (e quindi `qbar*alpha`). Ora e'
> corretto ovunque. L'effetto sui numeri e' grosso:
>
>       picco |alpha| col vecchio `+`  = 0.255 deg  ->  qbar*alpha = 20.7 kPa*deg
>       picco |alpha| col corretto `-` = 0.577 deg  ->  qbar*alpha = 46.8 kPa*deg
>
> **Un fattore 2.3 sull'indicatore di carico**, ed e' il numero con cui si
> argomentava che il progetto fosse load-relieving.

### La conseguenza fisica: l'attitude hold e' *load-aggravating*

Questo e' il punto d'orale vero, e vale la pena saperlo enunciare da solo.

Il picco di incidenza vale **0.577 deg**, mentre la raffica **da sola** produrrebbe
solo `alpha_w` di picco = **0.390 deg**. Cioe' l'anello chiuso **peggiora** il
bilancio di incidenza invece di alleggerirlo: il contributo dell'assetto non elide
il vento, **ci si somma**.

La ragione e' immediata una volta scritto il segno giusto. Il vento fa ruotare il
vettore velocita' relativa; per **tenere l'assetto** (`theta_ref = 0`) il loop deve
opporsi al momento aerodinamico, e cosi' facendo becca il muso **dentro** il vento
relativo. Il termine `theta + z_dot/V` cresce nello stesso verso in cui `-alpha_w`
gia' spinge, e i due si sommano in modulo. **Una legge di puro attitude-hold e'
load-aggravating, non load-relieving** -- il contrario di quello che il segno
sbagliato lasciava credere.

Non c'e' effetto banderuola che salvi la situazione: con `A_6 > 0` il centro di
pressione sta **davanti** al baricentro, quindi il momento aerodinamico e'
**divergente**. Il veicolo non si allinea da solo al vento -- ci diverge contro, e
il controllo deve spendere autorita' proprio per impedirglielo.

> **Possibile domanda d'esame** -- Perche' nel bilancio di incidenza `alpha_w` entra
> con il segno meno, e cosa ti dice il numero che ne esce?
> *Risposta:* Perche' l'Eq. (1) ha la colonna di vento `Bw = [0; -a1*V; 0; -A6]`,
> quindi entrambi i termini aerodinamici del plant forzano su
> `alpha = theta + z_dot/V - alpha_w`: `alpha` e' l'incidenza rispetto alla
> velocita' **relativa all'aria**, e il vento ruota quel vettore. Scrivere un `+`
> descriverebbe un plant diverso da quello che sto integrando. Il numero che ne
> esce e' istruttivo: il picco di `|alpha|` e' **0.577 deg** contro i **0.390 deg**
> del vento da solo. L'assetto **non** compensa il vento, **si somma**: per tenere
> `theta ~ 0` il loop becca il muso dentro il vento relativo. Quindi il mio Task 1
> e' un attitude-hold **load-aggravating**, e i guadagni di deriva a -1e-3 sono
> troppo piccoli per invertire il bilancio. Con `A_6 > 0` il centro di pressione e'
> davanti al baricentro: nessuna banderuola, il momento e' divergente. Un vero
> load relief richiederebbe guadagni di deriva molto piu' aggressivi (o un
> accelerometro), pagati in deriva laterale.

**Nota metodologica: il test, per come e' scritto, non poteva accorgersene.**
Congelando la formula invece di confrontarla con la dinamica, la tautologia ha
pinnato il segno sbagliato per tutto il tempo in cui c'e' stato. Un test non
tautologico avrebbe confrontato `r.alpha` con il segnale che **effettivamente
guida** la dinamica -- per esempio verificando che
`theta_ddot - K_1*delta = A_6*alpha` sulla risposta simulata. Questa e' la lezione
generale da portare all'orale: **un change-detector non e' un test di
correttezza**, e va detto quando se ne scrive uno.

---

## I golden value di questo file, in chiaro

| Valore | Riga | Da dove viene | Cosa succede se il progetto cambia |
|---|---|---|---|
| `Kp_th = 1.7845` | 11 | output di `fminsearch` in `design_controller` sul plant rigido, seed D'Antuono | falliscono 5 test (tutti quelli che usano `Kref`) |
| `Kd_th = 0.4433` | 11 | idem | idem |
| `Kp_z = Kd_z = -1e-3` | 12 | **fissati dalla traccia** (load relief), non ottimizzati | non cambiano |
| `6.0 dB` | 74, 81 | **target della traccia** (|GM| aerodinamico) | e' un requisito, non un risultato |
| `30.0 deg` | 75, 82 | **target della traccia** (PM rigido) | idem |
| `6.38 m/s` | 58 | interpolazione di `drywind.mat` (severe) a 15.143 km | cambia se cambia quota o severita' |
| `-10 dB` | 113 | soglia conservativa di gain stabilisation (il requisito vero e' >= 12 dB) | soglia scelta a mano |
| `0.002 / 0.7` | 110 | parametri del notch ritenuto | scelti dal trade a 4 filtri del Task 2 |

**Come leggere i golden value.** I 6 dB e i 30 deg sono **requisiti**: se il
progetto cambia, restano. `Kref`, `6.38`, `0.002/0.7` sono **fotografie di un
progetto specifico**: se il progetto cambia (altro seed, altra tolleranza del
tuner, altro data set, altro notch), **devono** cambiare, e la giusta reazione a un
fallimento non e' allargare la tolleranza ma **capire perche' il progetto si e'
spostato** e, se lo spostamento e' voluto, aggiornare il numero pinnato **insieme
al report**. Un golden value non e' una verita': e' un contratto con la versione
del progetto documentata nel report.

---

## Possibili domande d'esame

**D: Perche' non usi semplicemente `margin()`? Cosa aggiunge `classify_margins`?**
R: Perche' l'anello aperto ha due poli instabili e un integratore, quindi e'
**condizionatamente stabile**: la curva di Nichols passa *fra* i due punti critici
e ci sono piu' attraversamenti dello 0 dB. Su questo anello `allmargin` restituisce
tre margini di fase: -133 deg a 0.161 rad/s, -40 deg a 0.222 rad/s e +30 deg a 2.455
rad/s -- eppure l'anello chiuso e' stabile. I primi due sono artefatti del lobo di
deriva laterale. `margin()` restituisce **un** numero e quale numero scelga dipende
dal progetto: sui guadagni del Task 1 azzecca (-6 dB / 30 deg), su altri guadagni
che ho provato restituisce `Pm = -30 deg` su un anello stabile. `classify_margins`
separa i margini **per banda fisica** -- GM aerodinamico alla crossing di bassa
frequenza (dove il gain margin e' di **riduzione**, cioe' negativo in dB), PM rigido
al crossover di controllo, GM/PM flessionale attorno a omega_BM -- e non dipende
dalla fortuna.

**D: Perche' il margine di guadagno e' negativo (-6 dB)? Non dovrebbe dire di
quanto posso alzare il guadagno?**
R: Su un sistema condizionatamente stabile ci sono **due** limiti al guadagno: uno
superiore e uno **inferiore**. Il -6 dB dice che l'anello si destabilizza se il
guadagno **scende** di 6 dB, cioe' se il controllo perde meta' della sua autorita'.
Fisicamente e' il **margine aerodinamico**: sotto quella soglia il controllo non
riesce piu' a contrastare il polo instabile a +1.84 rad/s e il lanciatore diverge.
E' l'unico margine di guadagno che esiste sul modello rigido del Task 1 (con
attuatore ideale non c'e' crossing di aumento: `rigidGM_dB = NaN`, verificato).
Nel Task 2, con TVC e ritardo, ricompare anche il margine di aumento (7.56 dB
secondo il README), perche' l'attuatore aggiunge fase negativa in alta frequenza.

**D: Che cosa sono i tuoi golden value e cosa fai se un test fallisce?**
R: Sono numeri pinnati che fotografano **un progetto specifico**: `Kref =
(1.7845, 0.4433)` e' l'output di `fminsearch` a partire dal seed D'Antuono sul
plant rigido; `6.38 m/s` e' l'ampiezza di raffica interpolata dal file
`drywind.mat` alla mia quota. Non hanno derivazione analitica. I 6 dB e i 30 deg
invece **non** sono golden value: sono i **requisiti** della traccia. Se un test su
`Kref` fallisce, la reazione giusta **non** e' allargare la tolleranza: e' capire
cosa e' cambiato a monte (il seed, la tolleranza del tuner, il data set LPV, le
bande del classificatore) e, se il cambiamento e' voluto, aggiornare il numero
pinnato **insieme al report**, perche' i due devono raccontare la stessa storia.

**D: Come dimostri che il notch serve, e non l'hai messo per sicurezza?**
R: Con un test negativo. `testBareFullModelIsBendingUnstable` prende i guadagni del
Task 1, li mette sul modello a 6 stati con INS contaminato e TVC reale **senza
filtro**, e verifica che l'anello chiuso sia **instabile** (`verifyFalse(isstable(T))`).
Poi `testDeepNotchStabilisesFullModel` aggiunge il notch a fase minima e verifica
due cose: che l'anello torni stabile, **e** che il guadagno d'anello a omega_BM sia
sceso sotto -10 dB (misurato: -21.9 dB). La seconda asserzione e' quella che
qualifica la strategia: non e' stabile "per fase", e' stabile "per guadagno", cioe'
il lobo di bending e' cosi' basso che la sua fase e' irrilevante.

**D: Perche' il tuner ottimizza il logaritmo dei guadagni?**
R: Per due ragioni. Primo, garantisce che i guadagni restino **positivi** durante
tutta la ricerca senza dover imporre vincoli: `fminsearch` e' Nelder-Mead, un
metodo **non vincolato**, e senza il `log` potrebbe proporre Kp o Kd negativi,
producendo anelli con feedback positivo. Secondo, rende la ricerca
**moltiplicativa**: un passo del simplesso e' una variazione percentuale, non
assoluta, che e' la scala naturale per dei guadagni (passare da 0.4 a 0.5 e da 4 a
5 sono la stessa cosa in log). Il seed e' `log([2*A6/K1, sqrt(A6)/K1])`, cioe' la
coppia canonica di D'Antuono.

**D: I guadagni di deriva sono negativi. Non e' un errore di segno?**
R: No, e' *load relief*, ed e' voluto (la traccia li fissa a -1e-3). Un guadagno
**positivo** su `z` cercherebbe di riportare il lanciatore sulla traiettoria
nominale contro il vento, e per farlo dovrebbe **aumentare l'incidenza**
aerodinamica, cioe' aumentare il carico qbar*alpha -- proprio nel punto di volo in
cui il carico e' massimo. Un guadagno **negativo** lascia che il veicolo derivi
*con* il vento, riducendo l'incidenza al prezzo di qualche metro di deriva
laterale. A max-qbar il carico strutturale e' il vincolo dominante e la deriva no,
quindi si accetta il compromesso. **Ma il segno e' l'intenzione, non il risultato**:
a 1e-3 i guadagni sono cosi' piccoli che il bilancio di incidenza resta dominato
dal PD d'assetto, e misurando si trova un picco di `|alpha|` di 0.577 deg contro i
0.390 deg del vento da solo -- cioe' l'anello **aggrava** il carico. E' una
correzione debole, non un anello di guida, e non va spacciata per un load relief
che funziona.

**D: Il picco di incidenza e' piu' grande del vento da solo. Non e' un errore?**
R: No, ed e' il punto fisico piu' interessante dell'homework. Il bilancio corretto
e' `alpha = theta + z_dot/V - alpha_w` (lo impone l'Eq. (1) con
`Bw = [0; -a1*V; 0; -A6]`): `alpha` e' l'incidenza rispetto alla velocita'
**relativa all'aria**. Per tenere l'assetto (`theta_ref = 0`) il loop deve opporsi
al momento aerodinamico, e cosi' facendo **becca il muso dentro il vento
relativo**: il termine d'assetto si **somma** a quello di vento invece di eliderlo.
Risultato misurato: picco `|alpha|` = 0.577 deg contro `alpha_w` di picco = 0.390
deg, e `qbar*alpha` = 46.8 kPa*deg. **Un puro attitude-hold e' load-aggravating**,
non load-relieving. E non c'e' effetto banderuola a salvarlo: con `A_6 > 0` il
centro di pressione sta davanti al baricentro, quindi il momento aerodinamico e'
**divergente** -- il veicolo non si allinea al vento, ci diverge contro. (Onesta'
storica: fino a poco fa il post-processing aveva un `+` invece del `-` e dava
0.255 deg / 20.7 kPa*deg, facendo sembrare il progetto load-relieving. Era un bug
di post-processing -- il plant integrato aveva sempre il segno giusto, quindi
`theta`, `z` e `delta` non erano toccati -- ed e' stato corretto.)
