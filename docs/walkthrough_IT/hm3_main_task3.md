# HM3/main_task3.m

## Ruolo del file nel progetto

E' il **terzo entry point** di HM3 e copre il Task 3 della traccia ("Robust
Control Design", dichiarato *optional*). La domanda e' precisa: il coefficiente di
momento aerodinamico **mu_alpha = A_6** e l'efficacia di controllo
**mu_c = K_1** sono incerti; si valuti la robustezza del **controllore gia'
progettato** (quindi **senza ri-sintonizzarlo**) sui **quattro casi d'angolo**
ottenuti variando indipendentemente i due parametri di **+/-30%**.

Lo script e' concettualmente il piu' semplice dei tre -- e' un ciclo su nove casi
-- ma e' quello che produce il messaggio ingegneristico piu' forte
dell'homework: **il vincolo che stringe non e' il margine di fase, e' il margine
aerodinamico**, e il vertice peggiore V3 lo porta a **0.91 dB**, cioe' a un capello
dalla perdita di stabilita'.

Struttura:

    1. Ricostruisce il controllore FISSO del Task 2 (PD ri-sintonizzato + notch)
    2. Cicla su 9 casi:  Nominale + 4 vertici V1-V4 + 4 sensibilita' S1-S4
    3. Per ognuno: margini classificati per banda + risposta alla stessa raffica
    4. Due figure: overlay Nichols sui vertici, overlay risposta alla raffica

La distinzione **V** vs **S** e' il punto metodologico:

- I **V** sono i **vertici della scatola di incertezza** -- le combinazioni
  (mu_alpha, mu_c) = (0.7, 0.7), (0.7, 1.3), (1.3, 0.7), (1.3, 1.3). Sono i
  quattro corner cases che la traccia chiede letteralmente.
- Gli **S** sono **sensibilita' one-at-a-time**: si muove un parametro alla volta
  tenendo l'altro al nominale. Sono i **punti medi degli spigoli**, non i vertici.
  Il commento di header (righe 5-7) dice che "miss the worst combo (V3)" -- vero
  nel senso che nessuna singola riga S e' cattiva quanto V3. Ma come mostro sotto,
  le perdite in dB degli S **si sommano quasi esattamente** a quella di V3, quindi
  gli S sono in realta' un ottimo predittore lineare dei vertici.

---

## `%% Fixed controller (Task 2 retained design)` (righe 16-26)

```matlab
p0     = load_hw3_params();
notch  = struct('wx',p0.wBM,'zN',0.002,'zD',0.7,'sgn',+1);
Gfull0 = build_plant_full(p0,'ins');
Wact0  = build_tvc(p0,3) * build_notch_filter(notch.wx,notch.zN,notch.zD,notch.sgn);
K = design_controller(Gfull0, Wact0, 'w_flex',0.6*p0.wBM, 'w_flex_hi',1.5*p0.wBM, ...
                      'w_bending',p0.wBM, 'verbose',false);
```

- Righe 19-24: lo script **non hard-codifica** Kp = 1.732 / Kd = 0.687: li
  **ri-deriva** rilanciando `design_controller` sull'anello completo nominale,
  esattamente come fa `main_task2.m` alla riga 153. Essendo `fminsearch`
  deterministico da un seed fisso, riproduce gli stessi valori (verificato:
  Kp_th = 1.732, Kd_th = 0.687). Costa un run del tuner ma rende lo script
  **autosufficiente**: puo' essere lanciato da solo, senza aver prima eseguito
  `main_task2`. E' la scelta giusta per uno script che il professore fara' girare
  isolatamente.
- Riga 20: il notch e' **ricopiato a mano** come struct letterale
  (`zN = 0.002, zD = 0.7, sgn = +1`) invece di essere importato da un file di
  configurazione condiviso. **Limite di manutenibilita' da segnalare**: se si
  cambiasse il notch in `main_task2.m`, il Task 3 continuerebbe a usare quello
  vecchio senza avvisare. Non c'e' una singola sorgente di verita' per il
  progetto del filtro.
- Riga 25: la stampa conferma che il controllore congelato e' quello del **Task 2
  ri-sintonizzato**, non quello del Task 1. E' importante: usare i guadagni del
  Task 1 (Kd = 0.443) darebbe un margine di fase nominale di 14.6 gradi, e la
  campagna di robustezza partirebbe gia' in deficit.

---

## `%% Cases: box vertices + one-at-a-time sensitivities` (righe 28-67)

```matlab
cases = {
    'Nominal', 1.00, 1.00
    'V1',      0.70, 0.70      % vertici della scatola
    'V2',      0.70, 1.30
    'V3',      1.30, 0.70      % peggiore: max instabilita', min autorita'
    'V4',      1.30, 1.30
    'S1',      0.70, 1.00      % sensibilita' one-at-a-time (extra)
    'S2',      1.30, 1.00
    'S3',      1.00, 0.70
    'S4',      1.00, 1.30 };
```

- Riga 41: `nPlot = 5` -> **solo Nominale + V1..V4 vengono disegnati**. Gli S sono
  calcolati e tabulati ma non compaiono in nessuna figura.
- Riga 43: **la stessa raffica per tutti i casi** (`load_wind_profile(p0, ...)`,
  costruita su `p0` nominale). Corretto: il disturbo e' un dato d'ambiente, non
  deve cambiare con l'incertezza del veicolo. `Vg` dipende solo da quota e
  severita'.

### Il ciclo (righe 48-65)

```matlab
p  = load_hw3_params('mu_alpha_scale',cases{i,2}, 'mu_c_scale',cases{i,3});
Gf = build_plant_full(p,'ins');
Wf = build_tvc(p,3) * build_notch_filter(notch.wx,notch.zN,notch.zD,notch.sgn);
[L{i}, T] = assemble_loop(Gf, K, Wf);
mm{i} = classify_margins(L{i}, 'w_drift',0.3*sqrt(p.A6), ...);
r = simulate_gust_response(T, w);
```

- Riga 49: lo scaling e' applicato **dentro** `load_hw3_params` (righe 89-90 di
  quel file): `p.A6 = p.A6 * mu_alpha_scale`, `p.K1 = p.K1 * mu_c_scale`. **Nient'
  altro viene scalato**: a1, a3, a4, V, omega_BM, sigma_ins, phi_tvc restano
  nominali.

  **Osservazione fisica onesta**: A_6 = N_alpha*l_alpha/I_yy e
  a1 = -N_alpha/(m*V) **condividono N_alpha**. Una dispersione fisica reale del
  +30% sulla derivata di portanza muoverebbe **entrambi**. Il codice li disaccoppia.
  Non e' un errore: e' letteralmente cio' che la traccia chiede ("mu_alpha = A6 and
  mu_c = K1 are uncertain"), quindi l'incertezza e' definita **sui coefficienti**,
  non sulle grandezze fisiche a monte. Ma va detto: e' un'astrazione parametrica,
  non una dispersione fisica.

- Riga 51: `Wf` viene **ricostruito ad ogni iterazione** ma e' **identico in tutti
  i casi**: `build_tvc` dipende solo da wTVC/zTVC/tau (non scalati) e il notch e'
  centrato su `notch.wx = p0.wBM` (nominale). E' una **ridondanza innocua** ma
  merita di essere capita: significa che **il notch resta centrato sul nominale**,
  che e' proprio il senso di "controllore fisso".
- Righe 56-57: `classify_margins` riceve `w_drift = 0.3*sqrt(p.A6)` **calcolato sul
  corner corrente**. Sensato: il confine fra artefatti di deriva e corpo rigido si
  sposta con l'instabilita'. Le altre bande (`w_flex`, `w_flex_hi`, `w_bending`)
  usano `p.wBM`, che **non e' scalato**, quindi sono le stesse in tutti i casi.

### Risultati reali (run completo)

    Case       mu_a   mu_c |  AeroGM RigidPM  RigidGM  DM[ms] | peakTh  peakZ  stab
    Nominal    1.00   1.00 |   6.00    30.0     7.56     165 |  0.231   2.27    1
    V1         0.70   0.70 |   5.70    26.6     8.70     196 |  0.241   2.28    1
    V2         0.70   1.30 |  11.06    30.8     6.50     121 |  0.097   1.57    1
    V3         1.30   0.70 |   0.91    18.0     8.75     209 |  0.878   4.92    1
    V4         1.30   1.30 |   6.17    30.9     6.57     135 |  0.222   2.26    1
    S1         0.70   1.00 |   8.79    30.6     7.53     156 |  0.142   1.80    1
    S2         1.30   1.00 |   3.93    28.9     7.59     176 |  0.346   2.83    1
    S3         1.00   0.70 |   2.95    23.6     8.72     206 |  0.417   3.19    1
    S4         1.00   1.30 |   8.26    30.9     6.54     128 |  0.154   1.90    1

**Tutti e nove i casi restano stabili.** Il controllore fisso del Task 2 sopravvive
alla scatola +/-30%. Ma i margini si muovono moltissimo.

---

## La matematica dietro i numeri (derivazione, non nel codice)

### Perche' mu_c muove il margine di guadagno di esattamente 20*log10(mu_c)

K_1 compare **solo nella matrice B** (colonna di delta, riga di thetadot_dot).
Non tocca la matrice A, quindi **non sposta i poli del plant**. Scalare K_1 per
mu_c scala quasi esattamente il guadagno dell'intero anello aperto L. E un
margine di guadagno, misurato in **dB**, trasla di conseguenza:

    Delta(AeroGM_dB) = 20*log10(mu_c)

Verifica sui dati:

    S4 (mu_c = 1.30): previsto 20*log10(1.30) = +2.28 dB   |  misurato +2.26 dB
    S3 (mu_c = 0.70): previsto 20*log10(0.70) = -3.10 dB   |  misurato -3.05 dB

Accordo quasi perfetto. (Non e' esatto al 100% perche' la colonna di delta
(`build_plant_full.m`, riga 27: `Bd = [0; a3; 0; K1; 0; -phi_tvc*Tc]`) contiene
anche `a3` sulla riga di zdot e `-phi_tvc*Tc` sulla riga di etadot -- la sesta,
cioe' l'equazione del modo flessibile -- e **nessuno dei due viene scalato**:
solo il canale di assetto scala. Ma quel canale domina.)

### Perche' mu_alpha e' diverso

A_6 compare nella **matrice A** (riga di thetadot_dot) **e** nella colonna del
disturbo. Scalarlo **sposta i poli**: il polo instabile passa da
+sqrt(A_6) = 1.84 a +sqrt(1.3*A_6) = 2.10 rad/s. Non e' quindi una pura traslazione
di guadagno, e infatti l'accordo con -20*log10(mu_alpha) e' solo approssimativo:

    S2 (mu_alpha = 1.30): -20*log10(1.30) = -2.28 dB previsti  |  -2.07 dB misurati
    S1 (mu_alpha = 0.70): -20*log10(0.70) = +3.10 dB previsti  |  +2.79 dB misurati

### La scoperta: le perdite in dB si sommano

Se i due effetti fossero indipendenti e additivi in dB, il margine di un vertice
sarebbe prevedibile dagli spigoli:

    AeroGM(V) ~ AeroGM(nom) + Delta_S(mu_alpha) + Delta_S(mu_c)

Verifica su tutti e quattro i vertici (Delta_S presi dalle righe S misurate):

    V1 (0.7, 0.7):  6.00 + 2.79 - 3.05 =  5.74   |  misurato  5.70
    V2 (0.7, 1.3):  6.00 + 2.79 + 2.26 = 11.05   |  misurato 11.06
    V3 (1.3, 0.7):  6.00 - 2.07 - 3.05 =  0.88   |  misurato  0.91
    V4 (1.3, 1.3):  6.00 - 2.07 + 2.26 =  6.19   |  misurato  6.17

**Accordo entro 0.05 dB su tutti e quattro.** Il margine aerodinamico e'
essenzialmente **separabile e additivo in dB** nei due scaling logaritmici. Sono
tre conseguenze importanti:

1. **Il vertice peggiore e' prevedibile**: V3 e' il corner in cui entrambe le
   perdite sono negative. Non c'e' bisogno di cercarlo, si sa in anticipo che sara'
   (mu_alpha alto, mu_c basso).
2. **Gli S non "mancano" V3** come dice il commento di header: lo **predicono**
   con precisione, semplicemente vanno **sommati**. Il commento e' vero solo in
   senso letterale (nessuna singola riga S vale 0.91 dB), ma fuorviante sul
   metodo.
3. **La ricerca sui vertici e' giustificata a posteriori**: se la dipendenza e'
   monotona in entrambi i parametri, il peggior caso su una scatola sta
   necessariamente su un vertice. **Il codice pero' non lo dimostra** -- campiona 4
   vertici e 4 spigoli, e da questo non si puo' escludere formalmente un caso
   interno peggiore. La quasi-additivita' osservata lo rende molto improbabile, ma
   e' un argomento empirico, non una prova. Una vera prova richiederebbe
   un'analisi mu (structured singular value), che non e' nel codice.

### Perche' mu_c alto peggiora invece Rigid GM e DM

Guardando le colonne RigidGM e DM: mu_c = 1.30 (V2, V4, S4) le **abbassa**
(6.50-6.57 dB e 121-135 ms contro i 7.56 dB / 165 ms nominali), mentre mu_c = 0.70
(V1, V3, S3) le **alza** (8.70-8.75 dB, 196-209 ms). E' il rovescio esatto della
medaglia: il Rigid GM e' un margine di **aumento** di guadagno, quindi piu'
guadagno d'anello (mu_c alto) lo **erode**; l'Aero GM e' un margine di
**riduzione**, quindi piu' guadagno lo **migliora**.

**Ecco il vero trade-off di HM3**: i due margini si muovono in **direzioni
opposte** rispetto a mu_c. Non esiste un valore di guadagno che li massimizzi
entrambi -- e' esattamente la firma della **stabilita' condizionale**. Il progetto
sta in mezzo alla banda ammessa, e l'incertezza lo spinge verso un bordo o
verso l'altro.

### Perche' V3 amplifica anche la risposta alla raffica

    peak theta:  0.231 deg (nom)  ->  0.878 deg (V3)     = 3.8x
    peak z:      2.27 m           ->  4.92 m             = 2.2x

L'ampiezza dell'eccitazione da vento e' proporzionale ad **A_6** (il termine
-A_6*alpha_w nella colonna del disturbo), mentre la capacita' di reagire e'
proporzionale a **K_1**. V3 alza il primo del 30% e abbassa il secondo del 30%:
il rapporto A_6/K_1 passa da 0.74 a **1.37**, cioe' **+85%**. La risposta cresce
di conseguenza. Il caso opposto V2 (mu_alpha basso, mu_c alto) porta il picco a
0.097 deg, cioe' a meno della meta' del nominale.

---

## `%% Figures` (righe 69-105)

### f1 -- Nichols overlay sui corner (righe 77-92)

- Riga 79: `wv = logspace(-2, log10(30), 3000)` -- la griglia arriva **solo a
  30 rad/s**.
- Riga 81: `sh0 = 360*round((-180 - interp1(wv, ph1, mm{1}.rigidPM_w))/360)` -- lo
  **stesso** shift di fase (derivato dal caso nominale) e' applicato a **tutte** le
  curve, e ancora il crossover nominale sul ramo che contiene il punto critico a
  **-180 deg** (convenzione del corso, da 1 + L = 0; il commento riscritto alle
  righe 71-76 la dichiara). Essendo un multiplo esatto di 360 deg e' una
  rietichettatura del ramo, non una rotazione: la dispersione fra i corner resta
  confrontabile. Stesso trucco di `plot_nichols_lv` (riga 49) e di `main_task2.m`
  (riga 226). (Nota storica: fino a poco fa lo shift era riferito a +180, la
  rietichettatura D'Antuono della stessa carta spostata di +360; riallineato al
  corso, margini invariati.)
- Riga 88: `plot(ax, -180, 0, 'r+')` -- il punto critico e' marcato
  esplicitamente a (-180, 0 dB).
- Riga 89: `xlim [-270 -90]`, `ylim [-15 20]` -- **zoom sulla sola regione
  rigida**, centrato sul punto critico a -180 e tutto dentro il foglio nativo
  [-360, 0] della griglia di `ngrid` (riga 78).

  **Conseguenza da conoscere**: con `ylim` a [-15 20], il **lobo di bending non e'
  visibile** (sta a -18 dB). La figura non dice nulla sulla robustezza del modo
  flessibile. E' pero' una scelta **difendibile**, e vale la pena saper spiegare
  perche': lo scaling di mu_alpha e mu_c **quasi non tocca il bending**. Il modo e'
  forzato da -phi_tvc*Tc (non scalato) e rientra via sigma_ins (non scalato);
  omega_BM e zeta_BM non sono scalati. A_6 e K_1 muovono solo il **canale rigido**,
  che a 18.9 rad/s e' gia' rotolato via e conta poco contro la risonanza. Quindi
  |L(omega_BM)| e' praticamente lo stesso in tutti i corner, e non c'e' niente da
  vedere. Il posto dove il bending **viene** perturbato e' `main_montecarlo.m`, che
  disperde anche omega_BM, zeta_BM e il ritardo.

- La curva di **V3** e' quella che riga piu' vicino al punto critico a -180 deg:
  0.91 dB di margine sono, sul grafico, un capello.

### f2 -- Gust response overlay (righe 95-105)

Due pannelli (theta e z) con Nominale + V1..V4 sovrapposti. Mostra visivamente
l'escursione di V3 (0.88 deg / 4.9 m) contro il resto del gruppo.

### Export (righe 107-118)

Identico agli altri due script: `theme(f,'light')` in `try/catch`, PNG a 200 dpi in
`HM3/figures/`, prefisso `task3_`.

---

## Cosa il Task 3 NON copre (limiti dichiarati)

- **Nessuna incertezza su omega_BM, zeta_BM o sul ritardo.** Eppure lo Step D del
  Task 2 mostra che il notch profondo va **instabile a +5% su omega_BM**. Il Task 3
  potrebbe quindi dare un falso senso di sicurezza: dice che il controllore e'
  robusto a +/-30% sui due coefficienti aerodinamici, **ma il vero punto fragile del
  progetto e' altrove**. La combinazione (incertezza parametrica) x (detuning del
  bending) e' esplorata solo in `main_montecarlo.m`, che non e' un deliverable
  della traccia.
- **Nessuna prova formale che il caso peggiore stia su un vertice.** Si campionano
  4 vertici e 4 spigoli. La quasi-additivita' in dB (verificata sopra a 0.05 dB) e'
  un forte indizio di monotonia, ma non e' una dimostrazione.
- **Nessun gain scheduling.** Il messaggio implicito del task (V3 a 0.91 dB e'
  troppo poco per un progetto di volo) e' che servirebbe schedulare i guadagni;
  lo script non lo fa. Lo fa l'estensione `LTV_FULL_ASCENT/`, fuori traccia.

---

## Possibili domande d'esame

**D: Perche' il vertice V3 e' il peggiore, e si poteva prevederlo senza simulare?**
R: Si', si poteva. V3 e' (mu_alpha = 1.3, mu_c = 0.7): **velivolo piu' instabile**
(A_6 alto -> polo aerodinamico da 1.84 a 2.10 rad/s) **e meno autorita' di
controllo** (K_1 basso -> meno guadagno d'anello). Il margine aerodinamico e' un
margine di **riduzione** di guadagno: misura quanto guadagno puoi togliere prima di
non dominare piu' l'instabilita'. V3 aumenta cio' che va dominato e diminuisce cio'
con cui lo domini: i due effetti si **sommano**, e in dB lo fanno quasi
esattamente. Numericamente: la riga S2 (solo mu_alpha alto) costa 2.07 dB, la riga
S3 (solo mu_c basso) costa 3.05 dB; la somma e' 5.12 dB, e V3 perde 5.09 dB
(6.00 -> 0.91). Coincidenza entro 0.03 dB.

**D: Il commento nel codice dice che le sensibilita' S1-S4 "mancano" il caso
peggiore V3. E' vero?**
R: Letteralmente si', ma metodologicamente e' fuorviante. E' vero che **nessuna
singola riga S** e' cattiva quanto V3 (la peggiore, S3, si ferma a 2.95 dB contro
0.91). Ma le **perdite in dB delle righe S si sommano** a quella di V3 con un
errore di 0.03 dB, e lo stesso vale per tutti e quattro i vertici (accordo entro
0.05 dB). Quindi gli S non mancano affatto il vertice: lo **predicono** in modo
quasi esatto, purche' si sommino invece di leggerli isolati. La ragione e' che
l'Aero GM in dB e' **separabile** nei due scaling logaritmici -- perche' K_1 e'
un puro guadagno d'anello (Delta_dB = 20*log10(mu_c), verificato: -3.05 dB
misurati contro -3.10 previsti per mu_c = 0.7).

**D: Perche' mu_c alto migliora il margine aerodinamico ma peggiora quello rigido?**
R: Perche' sono margini di **segno opposto**. L'anello e' **condizionalmente
stabile**: c'e' una banda di guadagni ammessi, con un bordo sotto e uno sopra.
L'Aero GM misura la distanza dal **bordo inferiore** (quanto guadagno puoi
togliere), il Rigid GM dal **bordo superiore** (quanto puoi aggiungerne). Alzare
mu_c alza il guadagno d'anello: ti allontani dal bordo inferiore (Aero GM: 6.00 ->
8.26 dB) e ti avvicini a quello superiore (Rigid GM: 7.56 -> 6.54 dB; DM: 165 ->
128 ms). E' il trade-off strutturale di tutto HM3: **non esiste un guadagno che
massimizzi entrambi**, si sta in mezzo.

**D: Perche' la figura di Nichols del Task 3 non mostra il modo di bending?**
R: Perche' e' zoomata su ylim = [-15 20] dB e il lobo di bending sta a -18 dB. Ma
la scelta e' difendibile: **lo scaling di mu_alpha e mu_c non tocca praticamente il
bending**. Il modo e' forzato dal termine -phi_tvc*Tc*delta e rientra nella
retroazione via sigma_ins -- **nessuno dei due e' scalato**, e nemmeno omega_BM o
zeta_BM lo sono. A_6 e K_1 muovono solo il canale rigido, che a 18.9 rad/s e' gia'
rotolato via e non compete con la risonanza. Quindi |L(omega_BM)| resta
essenzialmente lo stesso in tutti i corner: non c'e' dispersione da mostrare. La
robustezza del bending si studia altrove (Step D del Task 2 per omega_BM,
`main_montecarlo.m` per tutto il resto).

**D: Il controllore sopravvive a tutti e nove i casi. Il progetto e' allora
robusto?**
R: **No, e questa e' la risposta corretta.** Sopravvive, ma V3 lascia **0.91 dB** di
margine aerodinamico e **18.0 gradi** di margine di fase: nessuno dei due e'
accettabile per un progetto di volo (le soglie tipiche sono 6 dB e 30 gradi, che
sono i target su cui il controllore e' stato tunato al nominale). Inoltre V3
triplica l'escursione di beccheggio (0.88 gradi contro 0.23) e raddoppia la deriva
(4.9 m contro 2.3). Il messaggio ingegneristico del Task 3 e' proprio questo: il
progetto a punto fisso **non ha abbastanza margine** per coprire una dispersione
+/-30%, ed e' l'argomento quantitativo a favore del **gain scheduling**. E c'e' un
secondo motivo, ancora piu' serio: il Task 3 **non disperde omega_BM**, e lo Step D
del Task 2 ha gia' mostrato che bastano **+5%** su omega_BM per far cadere l'anello.
La fragilita' vera del progetto non e' nemmeno in questa tabella.

**D: Perche' lo script ri-esegue il tuner invece di leggere i guadagni del Task 2?**
R: Per **autosufficienza**: cosi' `main_task3` puo' essere lanciato da solo, senza
aver prima eseguito `main_task2`. `fminsearch` parte da un seed analitico fisso
(Kp0 = 2*A_6/K_1, Kd0 = sqrt(A_6)/K_1) ed e' deterministico, quindi riproduce
esattamente Kp = 1.732 / Kd = 0.687. Il prezzo e' un run del tuner in piu' e, come
limite di manutenibilita', il fatto che il **notch e' ricopiato a mano** come
struct letterale alla riga 20: non c'e' una sorgente unica di verita' per il
progetto del filtro, quindi una modifica in `main_task2.m` non si propagherebbe
qui.
