# HM3/tests/hm3PlantTest.m

## Ruolo del file nel progetto

E' la suite `matlab.unittest` che **congela la fisica** di HM3. Copre tre file:
`load_hw3_params.m` (i coefficienti di Tabella 1 letti a t = 72 s dal data set
LPV), `build_plant_rigid.m` (il modello a 4 stati del Task 1) e
`build_plant_full.m` (il modello a 6 stati con il primo modo di bending, Task
2/3). Nessuna di queste funzioni fa controllo: costruiscono soltanto la matrice
`A`, `B`, `C` dell'Eq. (1) e dell'Eq. (2) della traccia. Se sono sbagliate,
tutto il resto dell'homework (Nichols, notch, margini, Monte Carlo) e' costruito
su sabbia -- da qui la scelta di testarle per prime e in modo molto letterale.

La domanda che questa suite si pone e': *il plant che ho in mano e' davvero il
lanciatore Greensite a max-qbar, o e' un sistema qualunque a 4/6 stati?* La
risposta viene data pinnando quattro classi di invarianti:

1. **Invarianti sui dati** -- i valori interpolati dal file
   `General/hw3-v3/GreensiteLPV_DATA.mat` a t = 72 s coincidono con i letterali
   di Tabella 1 (test 1), e il carico dinamico qbar e' coerente con l'atmosfera
   esponenziale (test 2).
2. **Invarianti di interfaccia** -- dimensioni e, soprattutto, i **nomi** dei
   segnali (test 4). Non e' cosmetica: `assemble_loop.m` chiude l'anello con
   `connect(...)`, che collega i blocchi **per nome**. Un `InputName` rinominato
   non rompe `build_plant_rigid`, rompe silenziosamente il loop.
3. **L'invariante fisico centrale** -- il polo instabile del corpo rigido cade a
   `+sqrt(A_6)` (test 5). E' la firma analitica dell'instabilita' aerodinamica
   ed e' il motivo per cui esiste l'homework.
4. **Invarianti sulla struttura di misura** -- la contaminazione INS del bending
   (Eq. 2) entra con i segni giusti (test 8), e sparisce nella variante `'true'`
   (test 9). E' la contaminazione che destabilizza l'anello a omega_BM e che
   motiva il notch del Task 2.

La suite gira in due minuti scarsi e non richiede nulla oltre alla Control
System Toolbox.

---

## Setup di classe (righe 7-20)

```matlab
properties
    p    % nominal parameter struct, loaded once per class
end

methods (TestClassSetup)
    function addHm3ToPath(testCase)
        hm3 = fileparts(fileparts(mfilename('fullpath')));
        testCase.applyFixture( ...
            matlab.unittest.fixtures.PathFixture(hm3));
    end
```

- Righe 7-9: la property `p` tiene lo struct dei parametri nominali. Essendo
  popolata in `TestClassSetup` (non in `TestMethodSetup`), `load_hw3_params()`
  viene chiamata **una volta sola per classe**, non una volta per test. E'
  corretto perche' lo struct e' immutabile e i test non lo modificano: quelli
  che vogliono parametri diversi (test 3) ricostruiscono uno struct locale.
- Righe 12-15: `mfilename('fullpath')` e' il path di questo file di test
  (`.../HM3/tests/hm3PlantTest.m`); il doppio `fileparts` risale di due livelli
  e restituisce `.../HM3`. La `PathFixture` aggiunge quella cartella al path
  MATLAB **e la rimuove alla fine della classe**. E' il modo pulito di rendere
  la suite eseguibile da qualunque cwd (`runtests('HM3/tests')` funziona anche
  dalla root della repo) senza sporcare il path dell'utente.
- Righe 17-19: caricamento dei parametri nominali (scala di incertezza = 1).

> **Possibile domanda d'esame** -- Perche' `PathFixture` e non un `addpath`
> secco?
> *Risposta:* `addpath` e' un effetto collaterale globale: se la suite fallisce
> a meta', il path resta sporco. La `PathFixture` e' gestita dal framework, che
> garantisce il ripristino anche in caso di eccezione (e' lo stesso principio
> dell'`onCleanup` usato in `design_controller.m` per le warning). Nei test
> serve isolamento: uno stato residuo puo' far passare o fallire il test
> successivo per ragioni che non c'entrano nulla col codice sotto test.

---

## `testParamsMatchTable1` (righe 23-31)

```matlab
pp = testCase.p;
testCase.verifyEqual(pp.A6,  3.3818, 'AbsTol', 5e-3);
testCase.verifyEqual(pp.K1,  4.5647, 'AbsTol', 5e-3);
testCase.verifyEqual(pp.V,   937.70, 'AbsTol', 0.5);
testCase.verifyEqual(pp.wBM, 18.9,   'AbsTol', 0.05);
testCase.verifyEqual(pp.a4,  -27.2710, 'AbsTol', 5e-2);
```

- Righe 26-30: cinque valori attesi hard-coded. Sono i **letterali di Tabella 1**
  della traccia:

      A_6 = mu_alpha = 3.3818 s^-2      (momento aerodinamico destabilizzante)
      K_1 = mu_c     = 4.5647 s^-2      (efficacia del TVC)
      V              = 937.70 m/s
      omega_BM       = 18.9 rad/s
      a_4            = -27.2710 s^-2

  Il punto sottile e' *perche' questo test non e' una tautologia*. In
  `load_hw3_params.m` (righe 40-48) quegli stessi numeri sono scritti come
  letterali, ma alle righe 59-80 **vengono sovrascritti** dal data set LPV (il ramo
  LPV riscrive i soli coefficienti tempo-varianti -- A6, K1, a1, a3, a4, V, Tc,
  sigma_ins, phi_ins, phi_tvc, wBM; `zBM`, `wTVC`, `zTVC` e `tau` restano i
  letterali di Tabella 1):

      p.A6 = interp1(L.A6.Time, squeeze(L.A6.Data), 72)

  Quindi il test confronta il valore **interpolato a t = 72 s** con il valore
  **stampato in tabella**. E' un cross-check fra due fonti indipendenti: se il
  file `.mat` non c'e', se `t_ref` fosse sbagliato, se il mapping dei campi LPV
  fosse invertito (`L.a3` scambiato con `L.a4`), o se l'interpolazione cadesse
  su un istante diverso, i numeri divergerebbero e il test fallirebbe. Verificato
  in MATLAB: la sorgente effettiva e' `GreensiteLPV_DATA.mat @ t=72 s` e i valori
  interpolati sono A6 = 3.3818, K1 = 4.5647, V = 937.71, wBM = 18.900,
  a4 = -27.2710. Combaciano.
- Le tolleranze sono asimmetriche e questo e' voluto: `5e-3` sui coefficienti
  (~0.1 %) perche' la tabella li stampa a 4 decimali, `0.5` su `V` (~0.05 %)
  perche' la tabella arrotonda a 937.70 mentre il dato LPV da' 937.71, `0.05`
  su `omega_BM` che in tabella e' a una sola cifra decimale. Le tolleranze
  seguono la **precisione di stampa della fonte**, non un criterio numerico.
- Nota di onesta': `a_4` merita attenzione. Il commento in `load_hw3_params.m`
  (righe 36-38) segnala che dai valori di massa e spinta della tabella si
  ricaverebbe `-(T_t - D)/m = -23.17 s^-2`, mentre la tabella dichiara
  `a_4 = -27.2710`. Il data set LPV concorda con -27.2710, e **quello e' il
  valore che entra nella dinamica**. Il test pinna -27.2710, cioe' pinna
  l'incoerenza cosi' com'e'. E' la scelta giusta (il test deve riflettere il
  codice), ma all'orale bisogna saperla difendere: e' una discrepanza nota della
  traccia, non un errore di implementazione.

> **Possibile domanda d'esame** -- Il tuo test verifica `a4 = -27.2710`, ma dai
> dati di Tabella 1 verrebbe -23.17. Quale usi e perche'?
> *Risposta:* Uso -27.2710, perche' e' il valore che compare sia nella riga
> `a_4` di Tabella 1 sia nel data set LPV di riferimento del professore, che sono
> due fonti indipendenti in accordo fra loro. Il -23.17 nasce dal ricalcolare
> `-(T_t - D)/m` con i valori di massa e spinta arrotondati della stessa tabella,
> e la discrepanza (~15 %) e' plausibilmente dovuta a quegli arrotondamenti o a
> un istante di riferimento leggermente diverso per T e D. `a_4` entra solo nella
> riga della dinamica laterale (accelerazione `zddot`), non nel momento di
> beccheggio, quindi non tocca il polo instabile ne' i margini: sposta solo
> leggermente la deriva laterale.

---

## `testDynamicPressureSelfConsistent` (righe 33-38)

```matlab
qbarExpected = 0.5 * 1.225 * exp(-pp.Alt/8000) * pp.V^2;
testCase.verifyEqual(pp.qbar, qbarExpected, 'AbsTol', 1e-6);
```

- Righe 36-37: ricalcola la pressione dinamica con l'**atmosfera esponenziale**

      rho(h) = rho_0 * exp(-h / H_scale),   rho_0 = 1.225 kg/m^3, H_scale = 8000 m
      qbar   = 0.5 * rho(h) * V^2

  A h = 15143 m e V = 937.7 m/s viene qbar = 81.1 kPa (verificato: 81129 Pa),
  che e' il valore di max-qbar citato nel README.
- **Onesta' totale**: questo test e' semi-tautologico. Riscrive la **stessa
  formula** delle righe 83-84 di `load_hw3_params.m`, con le stesse costanti
  hard-coded (1.225 e 8000). Non puo' quindi accorgersi di un errore *fisico*
  (per esempio: se `H_scale` fosse 7500 m, il test verrebbe aggiornato allo
  stesso modo). Cattura solo due cose: (a) che `p.qbar` sia effettivamente
  derivato e non un letterale scollegato da `V` e `Alt`, (b) che una futura
  modifica alla formula in `load_hw3_params.m` sia deliberata e non accidentale.
  Va classificato come *change-detector*, non come test di correttezza. Da
  notare inoltre che `p.Alt` resta un letterale (15143 m): non e' letto
  dall'LPV, quindi la quota non e' cross-checkata da nessun test.
- La `AbsTol` di `1e-6` su un numero dell'ordine di 8e4 e' di fatto una
  richiesta di identita' bit-a-bit modulo l'ordine delle operazioni in
  floating point.

---

## `testUncertaintyScalingAppliesToA6K1` (righe 40-47)

```matlab
ps = load_hw3_params('mu_alpha_scale', 1.3, 'mu_c_scale', 0.7);
testCase.verifyEqual(ps.A6, 1.3*pp.A6, 'AbsTol', 1e-12);
testCase.verifyEqual(ps.K1, 0.7*pp.K1, 'AbsTol', 1e-12);
testCase.verifyEqual(ps.a3, pp.a3,     'AbsTol', 1e-12);
```

- Riga 43: usa esattamente il vertice **V3** del box di incertezza del Task 3
  (mu_alpha +30 %, mu_c -30 %), che il README identifica come il caso peggiore
  (0.91 dB / 18.0 deg). Non e' una scelta casuale di numeri.
- Righe 44-46: l'invariante e' **quali coefficienti l'incertezza tocca e quali
  no**. `A_6` (= mu_alpha) e `K_1` (= mu_c) scalano; `a_3` (efficacia della
  forza laterale del TVC) **resta identico**. La terza asserzione e' quella
  interessante: e' un test *negativo*, garantisce che lo scaling non sia stato
  applicato con un ciclo su tutti i campi dello struct.
- Perche' e' l'invariante giusto: la traccia chiede di perturbare +/- 30 % *i due
  parametri aerodinamico e di controllo*, non l'intero modello. Se lo scaling
  si propagasse anche ad `a_1`, `a_3`, `a_4`, lo studio di robustezza del Task 3
  risponderebbe a una domanda diversa da quella posta, e i risultati non
  sarebbero confrontabili con quelli attesi. Nota fisica: mu_c scala `K_1` (il
  momento di controllo) ma **non** `a_3` (la forza laterale prodotta dalla
  stessa deflessione), il che e' fisicamente discutibile -- entrambe derivano
  dalla stessa spinta deviata. Il test pinna la scelta della traccia, non una
  verita' fisica.
- `AbsTol 1e-12` perche' e' una moltiplicazione esatta: `1.3*p.A6` calcolato nel
  test e nel codice deve dare lo stesso double.

> **Possibile domanda d'esame** -- Perche' lo scaling di robustezza non tocca
> `a_1`, `a_3`, `a_4`?
> *Risposta:* Perche' la traccia definisce la scatola di incertezza sui due
> parametri che dominano l'anello di beccheggio: mu_alpha (l'instabilita'
> aerodinamica, che fissa il polo a +sqrt(A_6) e quindi il margine di guadagno
> aerodinamico) e mu_c (l'autorita' di controllo, che fissa il guadagno d'anello
> in alta frequenza). Sono anche i due piu' incerti in volo reale, perche'
> dipendono dal centro di pressione e dalla spinta effettiva. I coefficienti
> laterali influenzano la deriva, non la stabilita' dell'anello di assetto,
> quindi restano nominali.

---

## `testRigidPlantDimensionsAndNames` (righe 49-55)

```matlab
G = build_plant_rigid(testCase.p);
testCase.verifySize(G.A, [4 4]);
testCase.verifySize(G.B, [4 2]);
testCase.verifyEqual(G.InputName,  {'delta'; 'alpha_w'});
testCase.verifyEqual(G.OutputName(1), {'theta_m'});
```

- Righe 51-52: 4 stati `[z, zdot, theta, thetadot]` e 2 ingressi
  `[delta, alpha_w]`. Il secondo ingresso non e' un comando: e' il **disturbo di
  vento**, trattato come un ingresso esogeno del plant. E' la scelta che permette
  poi a `assemble_loop` di restituire una `T` con ingressi `{alpha_w, theta_ref}`
  e di simulare la raffica con un semplice `lsim`.
- Righe 53-54: **questo e' il test importante della funzione**, ed e' facile
  scambiarlo per cosmesi. `assemble_loop.m` (righe 26-27, 31-32, 35) costruisce
  l'anello con

      T = connect(G, Kc, Wa, {'alpha_w','theta_ref'}, ...
                  {'theta','z','zdot','delta'}, {'delta'});

  e `connect` collega i blocchi **per nome di segnale**. Se `build_plant_rigid`
  rinominasse `theta_m` in `thetam`, il plant continuerebbe a essere corretto,
  ma `connect` non troverebbe piu' l'ingresso corrispondente del controllore e
  l'anello si chiuderebbe male (o `connect` lascerebbe un ingresso scollegato,
  producendo un `L` privo del ramo di feedback). Il test protegge il **contratto
  di interfaccia** fra il plant e tutto il resto della catena.
- Nota: `verifyEqual(G.InputName, {'delta'; 'alpha_w'})` usa un cell **colonna**
  (punto e virgola), perche' e' cosi' che `ss` normalizza `InputName`. Un cell
  riga farebbe fallire il confronto.

---

## `testRigidAirframeUnstablePole` (righe 57-63)

```matlab
% Aerodynamically unstable airframe: dominant pole at ~ +sqrt(A6)
% (the a1/a4 drift coupling shifts it ~1% off the pitch-only value)
G = build_plant_rigid(testCase.p);
testCase.verifyEqual(max(real(pole(G))), sqrt(testCase.p.A6), ...
    'RelTol', 0.02);
```

Questo e' **il** test della suite: cattura in una riga la ragione d'essere di
tutto HM3.

**Da dove viene `+sqrt(A_6)`.** La riga di beccheggio dell'Eq. (1) e'

    theta_ddot = A_6*theta + (A_6/V)*z_dot + K_1*delta - A_6*alpha_w

Se si isola la sola rotazione (si congela la deriva laterale: z_dot = 0,
alpha_w = 0, delta = 0), resta

    theta_ddot = A_6 * theta   ->   s^2 - A_6 = 0   ->   s = +/- sqrt(A_6)

cioe' **una coppia di poli reali simmetrici rispetto all'asse immaginario**, uno
dei quali nel semipiano destro. Non e' un'oscillazione smorzata male: e' una
divergenza esponenziale pura.

**Perche' A_6 e' positivo.** A_6 = mu_alpha = N_alpha * l_alpha / I_yy. E'
positivo perche' sul lanciatore **il centro di pressione sta davanti al centro
di massa** (l_alpha > 0 misurato in avanti). Un incremento di incidenza alpha
genera una forza normale che, applicata davanti al baricentro, produce un momento
che **aumenta ulteriormente alpha**. E' l'esatto opposto di un aereo (che ha il
CP dietro il CG ed e' staticamente stabile): il lanciatore e' un pendolo
rovesciato aerodinamico. Con A_6 = 3.3818 s^-2:

    sqrt(A_6) = 1.839 rad/s   ->   tempo di raddoppio = ln(2)/1.839 = 0.377 s

In meno di mezzo secondo l'errore d'assetto raddoppia: **il feedback non e' un
miglioramento delle prestazioni, e' una condizione di sopravvivenza**. E questo
fissa anche l'ordine di grandezza della banda: il crossover di controllo deve
stare comodamente sopra 1.84 rad/s (il progetto lo mette a 2.45 rad/s).

**Perche' la tolleranza e' relativa e vale 2 %.** Il plant completo non e' la
sola rotazione: la matrice `A` 4x4 accoppia beccheggio e deriva (il termine
`A_6/V` sulla riga 4, i termini `a_1` e `a_1*V + a_4` sulla riga 2). L'autovalore
esatto quindi **non** e' esattamente sqrt(A_6). Verificato in MATLAB, i poli del
plant rigido sono

    0,   -1.8610,   +0.0291,   +1.8165

Il polo dominante e' +1.8165 contro sqrt(A_6) = 1.8390: uno scarto dell'**1.22 %**,
compatibile con la `RelTol` di 0.02 e con quanto dichiara il commento alle righe
58-59. La `RelTol` (non `AbsTol`) e' la scelta giusta perche' lo scarto e'
proporzionale, non additivo: se nel Task 3 si scala A_6 di +30 %, il polo si
sposta a ~sqrt(1.3*A_6) e la stessa tolleranza relativa continua a valere.

**Il dettaglio che vale all'orale.** Nella lista dei poli ci sono **due poli a
parte reale positiva** (+1.8165 e +0.0291) e **uno nell'origine** (lo stato di
posizione `z`, che e' un integratore puro: nessuna forza richiama il lanciatore
verso la traiettoria nominale). L'accoppiamento con la deriva ha spezzato la
coppia simmetrica `+/- sqrt(A_6)` in una quaterna. Che il loop abbia due poli
instabili in anello aperto e' esattamente il motivo per cui il sistema e'
**condizionatamente stabile** e per cui un singolo numero di `margin()` non
significa nulla -- vedi `classify_margins.m` e la pagina su `hm3LoopTest`.

> **Possibile domanda d'esame** -- Perche' verifichi `max(real(pole(G)))` e non
> confronti tutti i poli con valori attesi?
> *Risposta:* Perche' l'invariante fisico che voglio proteggere e' *uno solo*: il
> polo dominante instabile, che e' la firma dell'instabilita' aerodinamica ed e'
> l'unica quantita' che ha una forma chiusa (+sqrt(A_6)) derivabile dalla teoria.
> Gli altri tre poli (l'integratore in zero, il polo stabile a -1.86, il polo
> lento a +0.029) sono artefatti dell'accoppiamento con la deriva laterale, non
> hanno un'espressione analitica pulita, e pinnarli significherebbe scrivere dei
> golden value senza contenuto teorico. Il test resta cosi' leggibile: dice
> "l'aereo diverge a sqrt(A_6)", che e' la frase che voglio poter difendere.

---

## `testRigidMeasurementsEqualTrueStates` (righe 65-70)

```matlab
G = build_plant_rigid(testCase.p);
testCase.verifyEqual(G.C(1,:), G.C(5,:), 'AbsTol', 1e-15);  % theta
testCase.verifyEqual(G.C(3,:), G.C(6,:), 'AbsTol', 1e-15);  % z
```

- Righe 68-69: la `C` di `build_plant_rigid` ha 7 righe, nell'ordine
  `[theta_m, thetadot_m, z_m, zdot_m, theta, z, zdot]`. Le prime quattro sono le
  **misure** che vanno al controllore; le ultime tre sono segnali di **plot**
  (gli stati veri). Il test verifica che, in assenza di bending, riga 1
  (`theta_m`) e riga 5 (`theta`) siano identiche, e cosi' riga 3 (`z_m`) e riga 6
  (`z`).
- Perche' e' un invariante e non una banalita': stabilisce il **punto di
  riferimento** rispetto a cui il Task 2 e' un cambiamento. Nel modello rigido
  misura e stato coincidono; nel modello completo con INS **non coincidono piu'**
  (test 8), ed e' proprio quello scarto che destabilizza l'anello a omega_BM. Il
  test fissa il "prima" del confronto prima/dopo.
- `AbsTol 1e-15`: sono zeri e uni esatti, non numeri calcolati. Si poteva usare
  `isequal`; la tolleranza e' innocua.

---

## `testFullPlantBendingMode` (righe 72-80)

```matlab
G  = build_plant_full(pp);
[wn, zeta] = damp(G);
[~, k] = min(abs(wn - pp.wBM));
testCase.verifyEqual(wn(k),   pp.wBM, 'AbsTol', 1e-9);
testCase.verifyEqual(zeta(k), pp.zBM, 'AbsTol', 1e-9);
```

- Righe 76-77: `damp(G)` restituisce frequenze naturali e smorzamenti di **tutti**
  i modi; la riga 77 seleziona quello con `wn` piu' vicina a `omega_BM = 18.9`
  rad/s. E' un modo robusto di dire "prendi il modo di bending" senza dipendere
  dall'ordinamento con cui `damp` restituisce gli autovalori (che non e'
  garantito). Funziona perche' gli altri modi stanno a 0 / 1.8 / 1.9 rad/s: sono
  un ordine di grandezza sotto, non c'e' ambiguita'.
- Righe 78-79: **l'invariante e' l'esattezza**, con `AbsTol 1e-9`. Non "il modo
  sta circa a 18.9", ma "sta a 18.9 con nove cifre". Perche' si puo' pretendere
  tanto? Perche' la matrice `A` di `build_plant_full` (righe 20-25) e' **diagonale
  a blocchi**:

      A = [ A_rigid(4x4)   0(4x2)     ]
          [ 0(2x4)         A_bend(2x2) ]

  Le righe 1-4 hanno colonne 5-6 nulle (il corpo rigido non "sente" eta); le
  righe 5-6 hanno colonne 1-4 nulle (il bending non e' forzato dagli stati
  rigidi). Quindi **lo spettro del plant e' l'unione esatta dei due spettri**, e
  il blocco di bending e' l'oscillatore canonico

      A_bend = [ 0        1       ]     ->  s^2 + 2*zeta_BM*omega_BM*s + omega_BM^2 = 0
               [ -w^2   -2*z*w    ]

  i cui autovalori sono per costruzione a frequenza naturale omega_BM e
  smorzamento zeta_BM. Il test verifica quindi che `build_plant_full` scriva
  correttamente quel blocco 2x2, e la tolleranza 1e-9 lascia margine solo
  all'errore numerico di `eig`. Verificato: `damp` restituisce esattamente
  18.9000 / 0.0050.
- **Il punto concettuale, che e' il piu' importante di tutta la pagina**: il
  bending **non** e' accoppiato al corpo rigido nella matrice `A`. L'accoppiamento
  passa da altre due strade:
  - via `B`: `Bd(6) = -phi_tvc * T_c` (riga 27) -- il TVC, deviando la spinta,
    eccita il modo flessionale;
  - via `C`: eta contamina le misure INS (righe 33-36) -- il sensore, montato sul
    corpo flessibile, legge assetto vero **piu'** deformazione.

  In **anello aperto** questi due percorsi non si chiudono, e infatti il plant
  completo ha esattamente gli stessi autovalori del plant rigido piu' la coppia
  di bending: il modo flessionale, con zeta_BM = 0.005, e' *poco smorzato ma
  stabile*. E' **la chiusura dell'anello** (delta -> eta -> misura -> delta) che
  crea il percorso di retroazione a omega_BM e lo rende instabile. Questo spiega
  perche' il test sul plant non puo' rilevare l'instabilita' di bending: quella
  vive in `hm3LoopTest` (`testBareFullModelIsBendingUnstable`).

> **Possibile domanda d'esame** -- Se il modo di bending e' stabile nel plant,
> perche' l'anello chiuso e' instabile?
> *Risposta:* Perche' l'instabilita' nasce dall'anello, non dal modo. In anello
> aperto la matrice A e' diagonale a blocchi e il bending e' un oscillatore
> stabile a zeta = 0.005. Chiudendo l'anello si crea il percorso
> delta -> eta (forzamento TVC, `B(6,1) = -phi_tvc*T_c`) -> theta_m
> (contaminazione INS, `C(1,5) = sigma_ins`) -> u_pd -> delta. A omega_BM lo
> smorzamento e' talmente basso che il guadagno d'anello ha un picco di +29 dB:
> se la fase in quel punto e' sfavorevole, il criterio di Nyquist da' instabilita'.
> La cura e' abbassare quel picco sotto 0 dB con un notch (gain stabilisation),
> che e' quello che fa il Task 2.

---

## `testInsMeasurementsContaminatedByBending` (righe 82-88)

```matlab
G  = build_plant_full(pp, 'ins');
testCase.verifyEqual(G.C(1,5),  pp.sigma_ins, 'AbsTol', 1e-15);
testCase.verifyEqual(G.C(3,5), -pp.phi_ins,   'AbsTol', 1e-15);
```

- Righe 86-87: pinna l'**Eq. (2)** della traccia, cioe' il modello di sensore:

      theta_m = theta + sigma_ins * eta        (sigma_ins = 0.178 rad/m)
      z_m     = z     - phi_ins   * eta        (phi_ins   = 0.8)

  `sigma_ins` e' la **pendenza** della forma modale nel punto in cui e' montata
  la piattaforma inerziale (rad di rotazione locale per metro di deformazione
  generalizzata): il giroscopio non misura l'assetto del corpo rigido, misura
  l'assetto **della sezione in cui e' avvitato**, che ruota anche per flessione.
  `phi_ins` e' invece l'**ampiezza** della forma modale nello stesso punto
  (spostamento laterale per unita' di eta).
- **L'invariante e' il segno**, e non e' un dettaglio contabile. Il segno di
  `C(1,5)` decide la **fase** con cui il lobo di bending compare nel guadagno
  d'anello `L`, e quindi se il modo si presenta sulla carta di Nichols dalla
  parte "buona" o "cattiva" del punto critico. Un errore di segno qui
  produrrebbe un anello che sembra stabilizzarsi senza notch (o che si
  destabilizza in modo diverso), e tutto il trade dei quattro filtri del Task 2
  risponderebbe a un problema che non esiste. I due segni opposti (+ su theta,
  - su z) sono l'unica cosa che questi due `verifyEqual` proteggono, ed e'
  esattamente cio' che serve.
- Nota: il test verifica solo le colonne di `eta` (colonna 5). Le colonne di
  `etadot` (colonna 6, righe 2 e 4 della `C`: `thetadot_m = thetadot +
  sigma_ins*etadot`, `zdot_m = zdot - phi_ins*etadot`) **non sono testate**. E'
  una lacuna minore ma reale: un errore di segno su `C(2,6)` passerebbe. Dato che
  il feedback derivativo `Kd_th * thetadot_m` e' proprio quello che tocca il
  bending in alta frequenza, sarebbe il termine da coprire per primo se si
  volesse rinforzare la suite.
- Nota 2: nemmeno il termine di forzamento `Bd(6) = -phi_tvc * T_c` (riga 27 di
  `build_plant_full.m`) e' coperto da un test. E' l'altra meta' dell'accoppiamento
  di bending.

---

## `testTrueMeasurementsBypassBending` (righe 90-94)

```matlab
G = build_plant_full(testCase.p, 'true');
testCase.verifyEqual(G.C(1:4, 5:6), zeros(4, 2), 'AbsTol', 1e-15);
```

- Riga 93: nella variante `'true'` il blocco di misura (righe 1-4 della `C`) ha
  **colonne di bending identicamente nulle**: il controllore riceve gli stati
  veri, come se il sensore fosse montato su un corpo rigido ideale.
- A cosa serve questa variante: e' il **controfattuale diagnostico**. Permette di
  dimostrare che l'instabilita' del Task 2 e' colpa della *contaminazione del
  sensore* e non della *dinamica flessibile in se'*. Con `'true'` il bending c'e'
  ancora nel plant (viene ancora eccitato dal TVC via `B(6,1)`), ma non torna
  indietro nel controllore, l'anello non si chiude su omega_BM e resta stabile.
  E' il modo pulito di isolare la causa. Un solo `verifyEqual` su una sottomatrice
  4x2 basta a pinnare l'intera semantica dell'opzione.

---

## `testFullPlantRejectsUnknownMeas` (righe 96-99)

```matlab
testCase.verifyError(@() build_plant_full(testCase.p, 'bogus'), ...
    'build_plant_full:meas');
```

- Righe 97-98: verifica che un valore non riconosciuto di `meas` sollevi l'errore
  con **identificatore** `build_plant_full:meas` (riga 43 del sorgente), non un
  errore generico. L'`arguments` block (righe 13-16) valida solo che `meas` sia
  testo scalare; e' lo `switch ... otherwise` a rifiutare i valori fuori dominio.
- Perche' conta: senza il ramo `otherwise`, uno `switch` MATLAB su una stringa
  sconosciuta **non fa nulla** e `Cm` resterebbe indefinita -> errore criptico piu'
  a valle, oppure (peggio, se `Cm` fosse preinizializzata) un plant silenziosamente
  sbagliato. Testare l'**identificatore** e non il messaggio e' la scelta corretta:
  i messaggi si riscrivono, gli ID sono il contratto d'errore.

---

## Cosa NON e' coperto (limiti noti)

Per onesta' e perche' all'orale puo' venire chiesto "quanto ti fidi di questa
suite":

- **Coefficienti `a_1`, `a_3` nella matrice `A`.** Nessun test verifica che le
  righe 2 e 4 della `A` (l'accoppiamento deriva-beccheggio) siano scritte
  correttamente. Un `a_1*V + a_4` diventato `a_1*V - a_4` passerebbe tutti i test
  di questa suite (sposterebbe il polo dominante, ma probabilmente ancora entro
  la `RelTol` del 2 %). L'unico presidio indiretto e' che i margini del Task 1
  cambierebbero e farebbe fallire `hm3LoopTest`.
- **Colonna del disturbo `Bw`.** Nessun test tocca `Bw = [0; -a1*V; 0; -A6]`. Il
  segno negativo davanti ad `A6` implica che l'incidenza aerodinamica interna al
  plant sia `alpha = theta + zdot/V - alpha_w`; vedi la nota di onesta' nella
  pagina su `hm3LoopTest`, dove questa convenzione **non** coincide con quella
  usata da `simulate_gust_response.m`.
- **`Bd(6) = -phi_tvc * T_c`** e le colonne di `etadot` nella `C`, come detto
  sopra.
- **`p.Alt`** e' un letterale mai cross-checkato con l'LPV.

---

## Possibili domande d'esame

**D: Il lanciatore ha un polo a +1.84 rad/s. Che cosa significa in pratica e come
lo hai verificato?**
R: Significa che, senza controllo, un errore di assetto cresce come e^(1.84*t):
raddoppia in ln(2)/1.84 = 0.38 s. E' l'instabilita' aerodinamica dovuta al centro
di pressione davanti al centro di massa (mu_alpha = A_6 > 0). Analiticamente,
congelando la deriva laterale la dinamica di beccheggio si riduce a
`theta_ddot = A_6*theta`, con poli +/- sqrt(A_6). Il test
`testRigidAirframeUnstablePole` confronta `max(real(pole(G)))` con `sqrt(A_6)` a
meno del 2 %: il plant reale, che include l'accoppiamento con la deriva, da'
+1.8165 contro sqrt(A_6) = 1.8390, cioe' l'1.2 % di scarto. Il polo instabile
esiste, sta dove la teoria dice, e fissa la banda minima di controllo.

**D: Perche' il tuo test sul modo di bending pretende una tolleranza di 1e-9,
mentre quello sul polo rigido si accontenta del 2 %?**
R: Perche' i due invarianti hanno natura diversa. Il modo di bending e' **esatto
per costruzione**: la matrice A e' diagonale a blocchi, il blocco flessionale e'
l'oscillatore canonico [0 1; -omega_BM^2, -2*zeta_BM*omega_BM], e i suoi
autovalori sono per definizione a frequenza omega_BM e smorzamento zeta_BM. Non
c'e' approssimazione, solo errore di `eig`, quindi posso pretendere 1e-9. Il polo
rigido invece e' `+sqrt(A_6)` solo nel **limite disaccoppiato**: nel plant a 4
stati l'accoppiamento con la deriva laterale (i termini a_1, a_4, A_6/V) lo sposta
dell'1.2 %. La tolleranza del 2 % e' la larghezza di quella approssimazione, non
sciatteria.

**D: Che cosa rende instabile l'anello nel Task 2, se il modo di bending nel
plant e' stabile?**
R: La retroazione. In anello aperto il bending e' disaccoppiato dal corpo rigido
nella matrice A e ha smorzamento 0.005 > 0, quindi e' stabile. Chiudendo l'anello
si crea il ciclo delta -> eta (il TVC eccita il modo, `B(6,1) = -phi_tvc*T_c`) ->
theta_m (l'INS misura la deformazione, `C(1,5) = sigma_ins = 0.178`) -> u_pd ->
delta. A omega_BM = 18.9 rad/s lo smorzamento bassissimo produce un picco di
risonanza di +29 dB nel guadagno d'anello, ben sopra 0 dB, e con fase sfavorevole:
Nyquist da' instabilita'. `hm3PlantTest` pinna gli ingredienti (il forzamento e la
contaminazione), `hm3LoopTest` pinna la conseguenza (l'anello nudo e' instabile).

**D: A cosa serve l'opzione `meas = 'true'` di `build_plant_full`, visto che il
lanciatore reale ha sempre un INS contaminato?**
R: E' uno strumento diagnostico, non un modello di volo. Serve a isolare la causa
dell'instabilita': con `'true'` il bending e' ancora nel plant e viene ancora
eccitato dal TVC, ma non torna nel controllore, e l'anello resta stabile. Questo
dimostra sperimentalmente che la colpa e' della **contaminazione della misura**
(Eq. 2) e non della dinamica flessibile in se'. Il test
`testTrueMeasurementsBypassBending` verifica che nella variante `'true'` le
colonne di eta ed etadot nel blocco di misura siano identicamente nulle.

**D: I test sui nomi dei segnali sembrano cosmetici. Perche' ci sono?**
R: Non sono cosmetici: sono un test di **contratto**. `assemble_loop.m` chiude
l'anello con `connect()`, che collega i blocchi **per nome** dei segnali
(`'delta'`, `'theta_m'`, `'alpha_w'`, ...). Se rinominassi un ingresso del plant,
`build_plant_rigid` continuerebbe a restituire un `ss` corretto e tutti i test
numerici passerebbero, ma `connect` non troverebbe piu' il segnale, e l'anello si
chiuderebbe in modo sbagliato o incompleto. Testare i nomi e' l'unico modo di
proteggere quell'accoppiamento implicito.

**D: `testDynamicPressureSelfConsistent` non e' una tautologia? Riscrive la stessa
formula del codice.**
R: In buona parte si', ed e' giusto ammetterlo: e' un *change-detector*, non un
test di correttezza fisica. Se `H_scale` fosse sbagliato, il test non se ne
accorgerebbe, perche' la costante 8000 e' hard-coded in entrambi i posti. Quello
che il test garantisce e' che `p.qbar` sia effettivamente **derivato** da `V` e
`Alt` e non sia un letterale scollegato: se domani cambio la quota di riferimento
o riscalo `V`, `qbar` si muove di conseguenza. Ha valore perche' `qbar` entra
nell'indicatore di carico qbar*alpha del Task 1, che e' il numero che uso per
giustificare il progetto a max-qbar.
