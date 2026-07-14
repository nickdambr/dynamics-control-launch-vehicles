# HM3/main_montecarlo.m

## Ruolo del file nel progetto

`main_montecarlo.m` e' un **extra** rispetto alla traccia di HM3: la traccia si
ferma al Task 3, cioe' alla robustezza valutata sui **quattro vertici** della
scatola di incertezza +/-30 % su mu_alpha (= A6) e mu_c (= K1). Questo script
sostituisce i vertici con **variabili aleatorie**: estrae N = 1500 campioni su
**cinque** fattori incerti (mu_alpha, mu_c, omega_BM, zeta_BM, ritardo TVC),
riassembla il loop aperto L e il loop chiuso T ad ogni estrazione, e produce
distribuzioni di margini, una probabilita' di stabilita', una nuvola di Nichols
e due scatter di sensibilita'.

La catena di dipendenze e' la stessa dei Task 1-3: `load_hw3_params` per i
parametri Greensite a t = 72 s, `build_plant_full(p,'ins')` per il modello a 6
stati con la contaminazione INS del modo flessionale, `build_tvc(p,3)` per
l'attuatore TVC piu' il ritardo approssimato con Pade di ordine 3,
`build_notch_filter` per il notch profondo, `assemble_loop` per chiudere il PD e
restituire (L, T), `simulate_gust_response` per la risposta alla raffica.
Le uscite sono `HM3/figures/mc_margins_hist.png`, `HM3/figures/mc_nichols_cloud.png`,
`HM3/figures/mc_sensitivity.png` e `HM3/mc_results.mat`.

Il file **non e' citato dal report LaTeX** (`HM3/report/`): compare solo nel
README di HM3 e in `HM3/docs/flowcharts.md`. E' materiale di contorno, non
consegnato.

---

## RIQUADRO DI ONESTA': lo script E' STALE rispetto ai Task 1-3

Ho verificato entrambi i punti sul codice. **Entrambe le incoerenze ci sono**,
piu' una terza.

### (a) I guadagni NON sono quelli del Task 2

L'intestazione (riga 2) dichiara: *"Task-2 controller (rigid PD + bending notch
fixed at nominal wBM) held fixed"*, e il banner di riga 39 dice
*"Fixed controller (Task-2 design)"*. Ma il codice reale e':

```matlab
Grigid = build_plant_rigid(p0);
K      = design_controller(Grigid, [], 'verbose', false);
```
(righe 41-42)

Questa e' **esattamente** la chiamata del Task 1 (`main_task1.m` riga 22:
`[K, m] = design_controller(G, []);`). Il secondo argomento `[]` viene
convertito in `tf(1)` da `design_controller` (riga 38), cioe' **attuatore
ideale**, e il plant e' quello **rigido a 4 stati**. Nessuna banda `w_flex` /
`w_bending` viene passata, quindi il classificatore interno gira con i default
"nessun modo flessionale".

Il Task 2 e il Task 3 fanno tutt'altro. `main_task3.m` righe 23-24:

```matlab
K = design_controller(Gfull0, Wact0, 'w_flex',0.6*p0.wBM, ...
      'w_flex_hi',1.5*p0.wBM, 'w_bending',p0.wBM, ...
      'verbose',false);
```

cioe' ri-sintonizzano il PD **sul loop pieno** (plant a 6 stati + TVC + ritardo
+ notch). Secondo il README di HM3 i due set di guadagni sono:

| origine | Kp_theta | Kd_theta |
|---|---|---|
| Task 1 (plant rigido, attuatore ideale) -- **quello che usa il Monte Carlo** | 1.78 | 0.44 |
| Task 2 ri-sintonizzato sul loop pieno -- quello dei Task 2 e 3 | 1.73 | 0.69 |

Il **notch** invece e' giusto: riga 43 costruisce `zN = 0.002`, `zD = 0.7`,
`sgn = +1` a `wx = p0.wBM`, identico al notch profondo ritenuto nel Task 2
(`main_task2.m` righe 74-77). Quindi il Monte Carlo monta il **notch del Task 2
sopra i guadagni del Task 1**: e' precisamente la configurazione che
`main_task2.m` chiama "BEFORE re-tuning" (righe 142-149) e che **scarta**,
perche' il ritardo piu' l'attuatore piu' il lag del notch fanno collassare il
margine di fase rigido da 30 gradi a 14.6 gradi, con delay margin a 98 ms
(numeri dal README, tabella del trade e testo del Task 2).

**Conseguenza operativa**: il "nominale" attorno a cui il Monte Carlo disperde
NON e' il progetto finale di HM3. E' il progetto intermedio bocciato. Le
distribuzioni di margine, la P(stabile) e le percentuali stampate a video
**non sono confrontabili** con la tabella del Task 3 del README (Nominale:
Aero |GM| = 6.00 dB, Rigid PM = 30.0 gradi, DM = 165 ms).

### (b) I margini NON passano da `classify_margins`

Nel corpo del `parfor` (righe 90-111) i margini sono letti cosi':

```matlab
am  = allmargin(L);
gmv = 20*log10(am.GainMargin);
gf  = am.GMFrequency;
...
minGM(i) = min(abs(gmv));
idx = find(gf>0.2 & gf<1, 1);
if ~isempty(idx), rigidGM(i) = abs(gmv(idx)); end
...
[~,Pm] = margin(L);
```

Nessuna chiamata a `classify_margins`. I problemi, in ordine di gravita':

1. **`margin(L)` restituisce il margine di fase PIU' PICCOLO** fra tutti gli
   attraversamenti a 0 dB. Ma il loop ha un lobo di drift laterale a bassa
   frequenza che genera attraversamenti spuri: e' scritto nero su bianco in
   `classify_margins.m` (righe 13-15) --
   *"Crossings below w_drift are lateral-drift artifacts ... Taking margin()'s
   default instead would pick one of these."* Il Monte Carlo fa esattamente
   quello che il resto di HM3 ha smesso di fare. Il vettore `PM` quindi non e'
   il "Rigid PM" delle tabelle dei Task 1-3.
2. **`minGM = min(abs(gmv))` mescola bande e butta via il segno.** Nel loop di
   un lanciatore instabile il segno del gain margin e' fisica, non convenzione:
   `gmdb < 0` e' il margine di **riduzione** di guadagno (bordo inferiore della
   banda condizionalmente stabile, il margine *aerodinamico*), `gmdb > 0` e' il
   margine di **aumento** di guadagno (bordo superiore, il margine *rigido*).
   `classify_margins` li separa proprio con i test `gmdb < 0` e `gmdb > 0`
   (righe 44-45). Prendendo il valore assoluto e il minimo su tutto, `minGM`
   e' un numero senza banda di appartenenza.
3. **`rigidGM` e' mal nominato.** La finestra hard-coded `gf > 0.2 & gf < 1`
   e il commento di riga 97 (*"low-freq aerodynamic crossover"*) dicono che il
   numero e' il margine **aerodinamico**, non quello rigido. Ma la variabile si
   chiama `rigidGM` e l'istogramma la etichetta `'rigid |GM| [dB]'` (riga 153)
   con la linea di riferimento a 6 dB, che e' il target *aerodinamico*. Il
   numero e' plausibilmente giusto, il nome e l'etichetta sono sbagliati.
   In piu' la finestra [0.2, 1] rad/s e' fissa mentre A6 disperde di +/-30 %:
   `find(..., 1)` prende la prima voce in finestra senza controllare il segno
   di `gmv`, cosa che `classify_margins` invece fa.
4. `PM(i) = 180*double(isStab)` (riga 108) quando non c'e' attraversamento di
   fase: e' un tappo arbitrario a 180 gradi. Va bene come segnaposto, ma
   inquina l'istogramma di `PM` e la statistica `P(|PM| >= 30)`.

### (c) Terzo scostamento: orizzonte della raffica

Riga 45: `w = load_wind_profile(p0);` -- **senza** `Tend`. Il default di
`load_wind_profile` e' `Tend = 12` s (riga 25 di quel file). Task 1, Task 2 e
Task 3 chiamano tutti `load_wind_profile(p, Tend=80)`, con un commento esplicito
(main_task1 riga 44) che 80 s servono perche' il modo lento di drift ha
tau ~ 18-20 s. Quindi `peakTh` e `peakZ` del Monte Carlo sono presi su
**12 s di simulazione**: il picco di theta (che avviene presto) e' probabilmente
catturato, ma il picco di z (drift, che cresce lentamente) puo' essere
**sottostimato**. Anche questi due numeri non sono confrontabili con le colonne
`peak th` / `peak z` della tabella del Task 3.

### Cosa e' comunque valido

La **struttura** dell'esperimento (campionamento, propagazione, statistiche,
nuvola di Nichols, scatter di sensibilita') e' corretta e riusabile. Le
**tendenze qualitative** -- omega_BM come driver dominante, instabilita' che
compare all'interno della scatola quando si disperdono anche flessione e
ritardo -- restano indicative. Sono i **valori numerici** a non essere
allineati. Per riallinearlo servirebbero tre modifiche: (1) costruire K con
`design_controller(build_plant_full(p0,'ins'), build_tvc(p0,3)*Wnotch, ...)`
passando le bande, (2) sostituire il blocco `allmargin`/`margin` con
`classify_margins(L, 'w_drift',0.3*sqrt(p.A6), 'w_flex',0.6*p.wBM, ...)`,
(3) mettere `Tend=80` nella raffica.

---

## Intestazione e configurazione (righe 1-37)

- Righe 1-21: docstring. Elenca i cinque fattori incerti e le loro leggi, e
  dichiara la dipendenza dal solo Control System Toolbox. La riga 15 contiene
  l'osservazione progettuale piu' importante: *"Notch NOT retuned when wBM
  disperses (deep notch needs near-exact wBM), so wBM is the dominant driver."*
  Cioe': il notch e' **congelato** a omega_BM nominale, quindi ogni disaccordo
  fra il notch e la vera frequenza di flessione e' pura perdita.
- Riga 24: `warning('off','Control:analysis:MarginUnstable')`. Serve perche'
  `margin`/`allmargin` avvisano ad ogni valutazione di un loop
  condizionalmente stabile. Nota: e' spento sul **client**; i worker di `parfor`
  non ereditano automaticamente lo stato dei warning del client, quindi in
  esecuzione parallela il messaggio puo' ricomparire (le funzioni di libreria
  `classify_margins` e `design_controller` se lo spengono localmente con
  `onCleanup`, qui invece no).
- Riga 27: `N = 1500` estrazioni. Riga 28: `Nsub = 150` sono le sole curve
  disegnate nella nuvola di Nichols (le altre 1350 non vengono ridisegnate).
  Riga 29: `wgrid = logspace(-2,2,600)` rad/s, la griglia su cui si valuta
  Nichols. Riga 30: `seed = 2026`, seme fisso.

### Specifica delle incertezze (righe 32-37)

```matlab
unc.mu_alpha = struct('dist','gauss','sigma',0.10,'trunc',0.30);
unc.mu_c     = struct('dist','gauss','sigma',0.10,'trunc',0.30);
unc.wBM      = struct('dist','gauss','sigma',0.02,'trunc',0.06);
unc.zBM      = struct('dist','lognorm','sigma',0.40);
unc.tau      = struct('dist','uniform','half', 0.25);
```

Tutti i fattori sono **moltiplicativi** sul valore nominale: p.A6 = p0.A6*fa, e
cosi' via. Le scelte:

- **mu_alpha, mu_c**: gaussiane con sigma = 10 % **troncate** a +/-30 %. La
  troncatura e' esattamente la scatola del Task 3, quindi +/-30 % = **3 sigma**.
  Questa e' una scelta di modellazione, non un dato: la traccia da' un intervallo
  (+/-30 %), non una distribuzione. Interpretare l'intervallo come 3 sigma di una
  gaussiana e' molto piu' ottimista che interpretarlo come supporto di una
  uniforme, e cambia le probabilita' che escono. Va detto all'orale.
- **omega_BM**: gaussiana sigma = 2 % troncata a +/-6 % (di nuovo 3 sigma). E'
  molto piu' stretta perche' e' il parametro a cui il notch profondo e'
  ipersensibile: il README documenta che con i filtri fissi il notch tollera
  -10 % ma va instabile a +5 %.
- **zeta_BM**: lognormale, `f = exp(0.40*randn)`. Mediana 1, **media**
  exp(0.5*0.40^2) = exp(0.08) = 1.083 (una lognormale non e' centrata sulla
  mediana). A 2 sigma il fattore sta in [exp(-0.8), exp(+0.8)] = [0.45, 2.23],
  che e' quanto dice il commento di riga 12: il commento e' corretto. La scelta
  lognormale e' fisicamente sensata perche' lo smorzamento e' **positivo per
  costruzione** (zeta_BM = 0.005 nominale) e non puo' diventare negativo: una
  gaussiana lo permetterebbe.
- **tau**: uniforme su +/-25 % (15-25 ms). Uniforme = "non so nulla dentro
  l'intervallo", la scelta agnostica per un ritardo di implementazione.

**Non c'e' Latin Hypercube, non c'e' Sobol, non c'e' correlazione fra i
fattori.** E' un Monte Carlo i.i.d. puro (`randn`/`rand`), con i cinque fattori
**indipendenti**. E' la scelta piu' semplice; costa un fattore ~1/sqrt(N) in
convergenza rispetto a un LHS a parita' di N.

> **Possibile domanda d'esame** -- perche' troncare la gaussiana invece di usare
> direttamente una uniforme sulla scatola +/-30 %?
> *Risposta:* Sono due modelli di ignoranza diversi. La gaussiana troncata a 3
> sigma dice "il valore vero e' quasi certamente vicino al nominale, i +/-30 %
> sono un caso estremo"; la uniforme dice "qualunque valore nella scatola e'
> ugualmente probabile, vertici inclusi". Con la gaussiana i vertici del Task 3
> hanno probabilita' praticamente nulla di essere campionati, quindi il Monte
> Carlo **non riproduce** il caso peggiore del Task 3 e le due analisi restano
> complementari. Con la uniforme le due si avvicinerebbero, al prezzo di una
> ipotesi molto piu' pessimistica.

---

## Controller fisso e modello nominale (righe 39-54)

- Righe 40-45: costruisce `p0` (nominale), il PD `K`, il notch `Wnotch` e il
  profilo di vento `w`. Vedi il riquadro di onesta' sopra per i due difetti
  (guadagni Task 1 anziche' Task 2; `Tend` di default a 12 s).
- Riga 48-50: loop nominale di riferimento.

```matlab
Gn      = build_plant_full(p0,'ins');
Ln      = assemble_loop(Gn, K, build_tvc(p0,3)*Wnotch);
[~, mn] = nichols_branch(Ln, wgrid, []);
```

  `mn` e' il **vettore di fase nominale**: serve solo come riferimento per
  allineare tutte le curve della nuvola sullo stesso ramo di 360 gradi (vedi
  `nichols_branch` piu' avanti). `Gn` usa la misura `'ins'`, cioe' Eq. (2):
  theta_m = theta + sigma_ins*eta, z_m = z - phi_ins*eta. E' proprio questa
  contaminazione del modo flessionale nelle misure a rendere necessario il notch.

---

## Pre-campionamento (righe 56-63)

```matlab
rng(seed);
fa = sample_factor(unc.mu_alpha, N);
fc = sample_factor(unc.mu_c,     N);
fw = sample_factor(unc.wBM,      N);
fz = sample_factor(unc.zBM,      N);
ft = sample_factor(unc.tau,      N);
```

Il commento di riga 57 spiega il perche': *"Sample OUTSIDE the parfor:
deterministic given the seed, loop body RNG-free."* E' il punto tecnico
importante del parallelismo. Se si chiamasse `randn` **dentro** il `parfor`,
ogni worker avrebbe il proprio stream e l'ordine di assegnazione delle
iterazioni ai worker non e' deterministico: il risultato non sarebbe
riproducibile fra un run seriale e uno parallelo, ne' fra pool di dimensioni
diverse. Estraendo tutto **prima**, con un `rng(seed)` singolo sul client, i
cinque vettori di fattori sono fissati una volta per tutte e il corpo del loop
diventa una funzione **deterministica** dell'indice `i`. Stesso seme, stessi
numeri, con o senza Parallel Computing Toolbox.

---

## Propagazione Monte Carlo (righe 65-118)

- Righe 66-72: pre-allocazione di sette vettori N x 1 (`nan` per gli scalari,
  `false(N,1)` per il flag di stabilita').
- Riga 75: `parfor i = 1:N`. **Senza Parallel Computing Toolbox MATLAB esegue
  `parfor` come un normale `for`**: nessun errore, solo esecuzione seriale. E'
  il "degrado seriale" citato nel README. Le variabili `p0, K, Wnotch, w, fa,
  fc, fw, fz, ft` sono *broadcast* (lette da tutti i worker), i sette vettori di
  uscita sono *sliced* (ogni iterazione tocca solo l'elemento `i`) -- e' la
  condizione che rende il loop parallelizzabile senza dipendenze fra iterazioni.
- Righe 77-82: perturbazione. `p = p0;` poi si scalano **cinque** campi:
  `A6`, `K1`, `wBM`, `zBM`, `tau`. Nota che si scala il **plant vero**, non il
  modello usato per il progetto: il controller `K` e il notch `Wnotch` restano
  quelli nominali. E' il senso stesso di un'analisi di robustezza: filtro
  progettato sul nominale, veicolo vero diverso.
- Righe 85-87: riassemblaggio.

```matlab
Gf      = build_plant_full(p,'ins');
Wf      = build_tvc(p,3) * Wnotch;
[L, T]  = assemble_loop(Gf, K, Wf);
```

  `build_tvc(p,3)` viene **ricostruito** con il `tau` perturbato (Pade di ordine
  3 sul nuovo ritardo), mentre `Wnotch` e' il blocco nominale, congelato. Il
  disaccordo notch-flessione nasce qui: `Gf` ha il modo a `p0.wBM*fw(i)`, il
  notch ha lo zero a `p0.wBM`.
- Righe 90-111: lettura dei margini. **E' il blocco problematico**, vedi il
  riquadro di onesta'. In sintesi: `allmargin` + `margin` grezzi, niente
  `classify_margins`.
- Riga 94: se `gmv` e' vuoto, `minGM = Inf` -- nessun attraversamento di guadagno
  significa margine di guadagno infinito, che e' formalmente giusto ma va
  gestito nelle statistiche (infatti `local_quantile` filtra con `isfinite`).
- Righe 103-104: `isStab = isstable(T)`. **Questo e' il vero verdetto di
  stabilita'**, ed e' l'unico che si puo' credere ciecamente. Il codice lo dice
  a voce alta nelle righe 132-134 dello stampato: con un veicolo instabile ad
  anello aperto il loop e' condizionalmente stabile, quindi il **segno** del gain
  margin non e' un indicatore di stabilita' e la sola cosa binaria affidabile e'
  la posizione dei poli del loop chiuso.
- Righe 114-116: risposta alla raffica, sui 12 s di default (vedi punto (c) del
  riquadro). `peak_theta` viene convertito in gradi, `peak_z` resta in metri.

> **Possibile domanda d'esame** -- perche' `isstable(T)` e non il segno del gain
> margin?
> *Risposta:* Il velivolo e' aerodinamicamente instabile: ha un polo di corpo
> rigido in +sqrt(A6) = +1.84 rad/s. Un sistema a ciclo aperto instabile e'
> stabile in ciclo chiuso solo dentro una **banda** di guadagni: esistono un
> margine di riduzione (sotto il quale il polo instabile non viene piu' vinto) e
> un margine di aumento (sopra il quale si eccita l'attraversamento ad alta
> frequenza). La curva di Nichols passa **fra** due punti critici. In questa
> situazione il criterio "GM positivo = stabile" non vale, perche' dipende da
> quanti giri fa la curva attorno al punto critico (Nyquist generalizzato). Il
> conteggio dei poli a parte reale positiva del loop chiuso, cioe' `isstable`,
> e' l'unico verdetto non ambiguo.

---

## Statistiche (righe 120-145)

- Riga 121: `Pstab = mean(stab)`. E' la stima Monte Carlo della probabilita' di
  stabilita': la media di una variabile di Bernoulli su N estrazioni
  indipendenti.
- Righe 122-123: `Pgm3 = mean(minGM >= 3)` e `Ppm30 = mean(PM >= 30)`. Sono
  probabilita' di **soddisfare un requisito**, non margini. Attenzione: entrambe
  ereditano i difetti di lettura del blocco 90-111, quindi vanno lette come
  indicative.
- Riga 126: `[~,iWorst] = min(minGM)` -- l'estrazione con il margine piu' stretto
  viene poi stampata con i suoi cinque fattori (righe 142-144). E' il modo
  corretto di riportare un Monte Carlo: non solo statistiche aggregate, ma anche
  il **campione peggiore**, con i valori dei parametri che l'hanno prodotto, in
  modo da poterlo riprodurre e analizzare a mano.
- Righe 135-141: percentili p05 / p50 / p95 per sei metriche via `print_pct`.
  L'uso di p05 e p95 anziche' media +/- sigma e' la scelta giusta: le
  distribuzioni dei margini non sono gaussiane (sono troncate a sinistra dalla
  perdita di stabilita', e `minGM` ha un atomo a +Inf).

Il commento stampato alle righe 132-134 e' importante e corretto:
`|GM|` e `|PM|` sono **magnitudini**, e l'indicatore vincolante e' il flag di
`isstable()`, non il segno del gain margin.

---

## Figura 1 -- istogrammi (righe 147-157)

Sei riquadri, uno per metrica, prodotti da `hist_metric`: `minGM`, `rigidGM`
(con linea di riferimento a 6 dB), `PM` (linea a 30 gradi), `DM`, `peakTh`,
`peakZ`. Ogni riquadro riporta la **mediana** come linea nera etichettata
(riga 310).

Come si legge un istogramma di margine:

- la **massa a sinistra della linea di riferimento** e' la frazione di
  estrazioni che violano il requisito: e' esattamente la probabilita' di
  violazione stimata dal Monte Carlo;
- la **coda sinistra** conta piu' del picco. Un progetto con mediana 8 dB e coda
  a 0.5 dB e' peggiore di uno con mediana 6 dB e coda a 4 dB;
- l'istogramma di `minGM` e' **troncato** dalla presenza di campioni instabili:
  quelli non hanno un margine sensato, quindi il fatto che l'istogramma sembri
  "sano" non dice nulla sulla P(instabile), che va letta separatamente nel
  titolo (riga 150).

**Nota di coerenza**: l'etichetta `'rigid |GM| [dB]'` con la linea a 6 dB e'
sbagliata, vedi punto (b)(3) del riquadro di onesta': quel numero e' il margine
aerodinamico.

---

## Figura 2 -- nuvola di Nichols (righe 159-186)

- Righe 164-167: allineamento globale.

```matlab
[gn, phn] = nichols_branch(Ln, wgrid, mn);
sc = find(gn(1:end-1).*gn(2:end) <= 0, 1, 'last');
if isempty(sc), shift = 360*round(median(phn)/360);
else,           shift = 360*round((phn(sc)+180)/360); end
```

  `find(gn(1:end-1).*gn(2:end) <= 0, 1, 'last')` cerca il **cambio di segno del
  guadagno in dB**, cioe' l'attraversamento a 0 dB, e prende **l'ultimo** (il
  piu' alto in frequenza). Poi calcola uno spostamento di fase multiplo di 360
  gradi che porta quell'attraversamento vicino a -180 gradi -- lo stesso ruolo
  del `PhaseMatchingValue = -180` di `plot_nichols_lv` (riga 40 di quel file).
  E' cosmetica: sposta
  la curva sul ramo di fase "giusto" del piano di Nichols, senza cambiare la
  fisica (aggiungere 360 gradi alla fase non cambia il sistema).
- **Convenzione di fase**: il punto critico e' a **(-180 gradi, 0 dB)**
  (righe 182-184 marcano le copie -180 + 360k), cioe' la convenzione degli
  appunti del corso, da 1 + L = 0 <=> L = -1. Oggi e' la convenzione di **tutto**
  HM3: `plot_nichols_lv` e i tre main sono stati riallineati a -180, quindi la
  nuvola Monte Carlo e' direttamente confrontabile a occhio con le altre
  Nichols dell'homework. (Nota storica: fino a poco fa i Task 1-3 usavano la
  rietichettatura D'Antuono con il punto critico a +180 -- stessa curva, fase
  mod 360, spostata di un giro -- e questo script era l'**eccezione**; dopo il
  riallineamento non lo e' piu', e infatti **non e' stato toccato**. Il flip e'
  una pura rietichettatura dell'asse: nessun margine cambia.)
- Righe 169-178: si ridisegnano `Nsub = 150` estrazioni prese con
  `round(linspace(1,N,150))`, cioe' un sottoinsieme **regolarmente spaziato negli
  indici**, non un campione casuale del campione -- va bene perche' gli indici
  non sono ordinati per gravita', quindi un passo regolare e' equivalente a
  un'estrazione casuale.
- Riga 176: colore. Grigio trasparente (alpha 0.18) se stabile, rosso (alpha
  0.35) se instabile. Riga 179: la curva nominale in blu spesso sopra tutte.
- **Costo**: il ciclo delle righe 170-178 **riassembla il loop da zero** per
  ognuna delle 150 estrazioni (`build_plant_full`, `build_tvc`,
  `assemble_loop`, `nichols`), duplicando lavoro gia' fatto dentro il `parfor`.
  E' inefficiente ma innocuo (150 assemblaggi contro i 1500 gia' fatti).

**Come si legge la nuvola**: la dispersione dei parametri produce un **fascio**
di curve. Se una parte del fascio circonda diversamente il punto critico rispetto
al nominale, li' c'e' un cambio di stabilita'. Le curve rosse mostrano *dove* la
curva va a sbattere: con il notch congelato e omega_BM disperso, ci si aspetta
che il lobo flessionale -- che nel nominale sta a -18 dB, ben sotto lo 0 dB -- si
alzi e risalga verso il punto critico quando la vera omega_BM esce dal null del
notch.

---

## Figura 3 -- sensibilita' (righe 188-216)

- Righe 196-205, pannello (a): piano (mu_alpha, mu_c), estrazioni stabili in
  verde e instabili con croce rossa, con sovrapposti il **rettangolo +/-30 % del
  Task 3** e i suoi quattro vertici V1..V4.

  Il punto sottile: dato che mu_alpha e mu_c sono gaussiane troncate con
  sigma = 10 % e troncatura a 3 sigma, la nuvola di punti **si concentra al
  centro** del rettangolo e i vertici marcati V1..V4 restano praticamente
  **vuoti**. La probabilita' di finire simultaneamente a 3 sigma su entrambi i
  fattori e' dell'ordine di (1.35e-3)^2 ~ 2e-6: con N = 1500 non capita mai.
  La figura mostra quindi due popolazioni diverse: i vertici (deterministici,
  Task 3) e la nuvola (probabilistica, interno della scatola). Il commento delle
  righe 193-195 lo dice a modo suo: *"Task 3 found all four vertices stable in
  (mu_alpha, mu_c) alone; dispersing bending/delay too makes instability appear
  inside the box."*
- Righe 208-216, pannello (b): `minGM` contro la dispersione percentuale di
  omega_BM. E' la verifica della tesi dichiarata a riga 15: il **detuning del
  notch e' il driver dominante**. Se l'ipotesi e' vera, i punti rossi
  (instabili) si addensano ad un estremo dell'asse x e il margine crolla in modo
  monotono al crescere del disaccordo.

> **Possibile domanda d'esame** -- se il notch e' cosi' sensibile a omega_BM,
> perche' non lo si ri-sintonizza ad ogni estrazione?
> *Risposta:* Perche' sarebbe barare. In volo il controllore e' congelato: la
> vera omega_BM del veicolo non e' misurabile in tempo reale, e' proprio
> l'incertezza che si vuole quantificare. Ri-sintonizzare il notch su ogni
> estrazione risponderebbe alla domanda "quanto e' buono il metodo di progetto",
> non "quanto e' robusto **questo** controllore". La riga 15 del sorgente lo
> dichiara esplicitamente. Il prezzo e' che il notch profondo (zeta_N = 0.002,
> null strettissimo, largo ~2*zeta_N*omega_BM ~ 0.08 rad/s) e' fragile: e' l'argomento a favore
> dell'alternativa notch + lead-lag, che tiene tutta la banda +/-10 % con margini
> nominali piu' magri.

---

## Export (righe 218-237)

- Righe 221-228: forza il tema chiaro (`theme(f,'light')` in `try`, con
  `f.Color = 'w'` come fallback per MATLAB pre-R2025a) ed esporta i PNG a 200
  dpi con il prefisso `mc_`.
- Righe 231-236: salva **tutti i dati grezzi** in `mc_results.mat`: i cinque
  vettori di fattori, le sei metriche per estrazione, il flag di stabilita' e le
  tre probabilita'. E' la cosa giusta da fare: le figure sono ricostruibili dai
  dati, i dati no dalle figure. Con `fa..ft` e `stab` salvati, si puo' rifare
  qualunque post-processing (per esempio ri-leggere i margini con
  `classify_margins`, se non fosse che L non e' salvato).

---

## Funzioni locali (righe 239-315)

### `sample_factor` (righe 240-263)

```matlab
case 'gauss'
    f = 1 + spec.sigma*randn(n,1);
    if isfield(spec,'trunc') && ~isempty(spec.trunc)
        lo = 1-spec.trunc; hi = 1+spec.trunc;
        bad = f<lo | f>hi;
        while any(bad)
            f(bad) = 1 + spec.sigma*randn(nnz(bad),1);
            bad = f<lo | f>hi;
        end
    end
```

- Righe 246-255, `'gauss'`: gaussiana con **troncatura per rifiuto**
  (*rejection sampling*): i campioni fuori dalla scatola vengono ri-estratti
  finche' cadono dentro. Questo produce la **gaussiana troncata** corretta
  (densita' rinormalizzata sul supporto), a differenza del **clipping**
  (`f = min(max(f,lo),hi)`) che accumulerebbe massa di probabilita' sui bordi --
  creerebbe due atomi in lo e hi, cioe' una distribuzione fisicamente assurda.
  La differenza conta: con clipping la scatola avrebbe i bordi "pesanti" e la
  P(violazione) sarebbe sovrastimata.
- Il `while` termina quasi certamente: con sigma = 0.10 e troncatura a 3 sigma la
  probabilita' di rifiuto e' ~0.27 % per campione, quindi il ciclo si chiude
  praticamente al primo giro. Ma **non c'e' guardia sul numero di iterazioni**:
  se qualcuno impostasse `trunc` molto minore di `sigma` il ciclo diventerebbe
  lentissimo (probabilita' di accettazione minuscola). E' un rischio latente,
  non un bug attivo con questi numeri.
- Righe 256-257, `'uniform'`: `f = 1 + half*(2*rand-1)`, supporto [1-half, 1+half].
- Righe 258-259, `'lognorm'`: `f = exp(sigma*randn)`. Mediana 1 esatta, media
  exp(sigma^2/2) > 1, supporto (0, +Inf) -- **strettamente positivo**, che e' il
  motivo per cui e' la legge giusta per uno smorzamento.

### `nichols_branch` (righe 265-278)

```matlab
[mag, phase] = nichols(L, w);
g  = 20*log10(squeeze(mag));
ph = squeeze(phase);
if ~isempty(ref)
    ph = ph - 360*round((ph(1)-ref(1))/360);
end
```

Il problema che risolve: `nichols` restituisce la fase **srotolata** (unwrapped),
ma il ramo su cui finisce dipende dal sistema. Due loop quasi identici possono
uscire su rami distanti 360 gradi, e nella nuvola apparirebbero come due
popolazioni separate anche se fisicamente coincidono. La correzione confronta la
fase **al primo punto della griglia** (w = 1e-2 rad/s, cioe' molto sotto la
dinamica) con quella nominale e sottrae il multiplo di 360 piu' vicino. E'
lecito perche' una fase definita modulo 360 gradi non cambia il sistema.

**Limite**: l'allineamento e' fatto su **un solo punto** (`ph(1)`). Se una
estrazione avesse una fase a bassa frequenza gia' spostata di piu' di 180 gradi
rispetto al nominale, l'arrotondamento potrebbe agganciare il ramo sbagliato. Con
questi livelli di dispersione e' improbabile, ma e' un allineamento fragile.

### `print_pct` e `local_quantile` (righe 280-298)

```matlab
x = sort(x(isfinite(x)));
n  = numel(x);
pos = min(max(p/100*n + 0.5, 1), n);
qv  = interp1(1:n, x, pos, 'linear');
```

Percentili **senza lo Statistics Toolbox**. La regola implementata e' quella del
"midpoint": il k-esimo campione ordinato rappresenta il percentile
(k - 0.5)/n * 100. Invertendo, il percentile p corrisponde alla posizione
`pos = p/100 * n + 0.5`, con interpolazione lineare fra i due campioni adiacenti
e saturazione agli estremi ([1, n]). E' la stessa convenzione di `quantile` di
MATLAB con il metodo di default. Il filtro `isfinite` (riga 293) e' necessario
perche' `minGM` e `DM` possono valere +Inf e `rigidGM` puo' essere NaN quando
l'attraversamento non cade nella finestra hard-coded [0.2, 1] rad/s.

**Effetto collaterale da tenere a mente**: scartando i +Inf e i NaN, i percentili
sono calcolati su un **sottoinsieme** dei campioni, la cui dimensione non viene
stampata. Un p05 di `rigidGM` calcolato su 900 campioni su 1500 e' una statistica
diversa da un p05 su 1500, e il codice non lo segnala.

### `hist_metric` (righe 300-315)

30 bin, mediana come linea nera etichettata (`xline`), linea di riferimento
opzionale (6 dB oppure 30 gradi) tratteggiata. Usa `median(v,'omitnan')` dopo
aver gia' filtrato con `isfinite`, quindi l'`omitnan` e' ridondante -- innocuo.

---

## Senso statistico: vertici (Task 3) contro Monte Carlo

Questo e' il punto teorico da difendere all'orale.

**Cosa fa il Task 3.** Prende la scatola di incertezza
[0.7, 1.3] x [0.7, 1.3] su (mu_alpha, mu_c) e valuta i **quattro vertici**.
La logica implicita e': *se la stabilita' (o il margine) e' peggiore ai vertici,
allora la stabilita' ai vertici implica la stabilita' ovunque nella scatola*.

**Quando questa logica e' valida.** Solo se la dipendenza della metrica dai
parametri e' **monotona** in ciascun parametro (o piu' in generale se il minimo
del margine sulla scatola e' assunto sulla frontiera, e sulla frontiera in un
vertice). Per una funzione monotona in ogni argomento, il minimo su un
ipercubo cade sempre in un vertice. Intuitivamente qui e' plausibile:

    mu_alpha su  -> velivolo piu' instabile -> margine aerodinamico peggiore
    mu_c    giu' -> meno autorita' di controllo -> margine peggiore

quindi il vertice V3 = (1.3, 0.7) e' "il peggiore dei quattro" e infatti il README
lo conferma (0.91 dB / 18.0 gradi). **Ma la monotonia e' una congettura, non un
teorema**: nessuno la dimostra. Con quattro valutazioni non si puo' escludere
che il minimo cada **dentro** la scatola. Basta un'interazione non monotona --
per esempio un attraversamento a 0 dB che cambia banda al variare di un
parametro, o una risonanza che entra o esce da un filtro -- e i vertici non sono
piu' il caso peggiore. Con **cinque** parametri incerti la faccenda peggiora: i
vertici diventano 2^5 = 32, e comunque nessuno garantisce la monotonia in
omega_BM: li' il notch ha un null strettissimo, quindi l'attenuazione ha un
massimo al centro del null e degrada da **entrambi** i lati -- una dipendenza
palesemente **non monotona**, per cui i vertici in omega_BM sarebbero anzi i
punti migliori, non i peggiori.

**Cosa aggiunge il Monte Carlo.**

1. **Campiona l'interno**, non solo i vertici. Se il minimo sta dentro, un
   campionamento denso lo trova (approssimativamente).
2. **Restituisce una probabilita', non un verdetto binario.** Il Task 3 dice
   "tutti e quattro i vertici sono stabili". Il Monte Carlo dice "P(stabile) =
   x %", che e' l'unica cosa con cui si puo' fare un budget di affidabilita' di
   missione. Il caso peggiore assoluto e' spesso inutilizzabile: se il progetto
   e' instabile in un vertice a probabilita' 1e-6, e' un problema o no? Il
   ragionamento sui vertici non sa rispondere, quello probabilistico si'.
3. **Permette di aggiungere fattori** che la traccia non disperde (omega_BM,
   zeta_BM, tau). Ed e' proprio l'aggiunta di questi che fa comparire
   instabilita' *dentro* la scatola (mu_alpha, mu_c) dove il Task 3 non ne
   trovava: il Task 3 non e' sbagliato, sta semplicemente rispondendo a una
   domanda su uno spazio di incertezza **piu' piccolo**.

**Cosa il Monte Carlo NON aggiunge.**

1. **Non e' un certificato.** Non copre il caso peggiore: con code gaussiane
   troncate, i vertici hanno probabilita' ~2e-6 e con N = 1500 non vengono mai
   estratti. Quindi il Monte Carlo **non sostituisce** il Task 3: sono
   complementari. Un certificato vero richiederebbe mu-analisi o LMI su un
   politopo di incertezza.
2. **La risposta dipende dalle distribuzioni assunte**, che sono un'ipotesi di
   modellazione dell'ingegnere, non un dato della traccia. Cambiando le
   gaussiane troncate in uniformi, la P(violazione) cambia di ordini di
   grandezza.
3. **La precisione va come 1/sqrt(N).** Per una probabilita' p stimata su N
   campioni indipendenti, l'errore standard e'

       SE = sqrt( p*(1-p) / N )

   Con N = 1500 e p ~ 0.01, SE ~ 0.26 %, cioe' un errore **relativo** del ~26 %
   sulla stima di una probabilita' dell'1 %. Le probabilita' piccole sono
   proprio quelle mal stimate. Caso limite: se su 1500 estrazioni non se ne
   osservano di instabili, la "regola del tre" da' come limite superiore al 95 %
   di confidenza p <= 3/N = 0.2 %. **Non** p = 0.

> **Possibile domanda d'esame** -- il Monte Carlo trova un'estrazione instabile
> e i quattro vertici del Task 3 sono tutti stabili. Contraddizione?
> *Risposta:* No. Sono spazi di incertezza diversi. Il Task 3 disperde solo
> mu_alpha e mu_c; il Monte Carlo disperde anche omega_BM, zeta_BM e tau, e il
> notch profondo e' congelato alla omega_BM nominale. L'instabilita' che compare
> e' guidata dal detuning del notch, un meccanismo che il Task 3 non esplora
> proprio. Se si fissassero omega_BM, zeta_BM e tau al nominale, la nuvola Monte
> Carlo tornerebbe tutta stabile dentro la scatola.

---

## Possibili domande d'esame

**D: Perche' in questo homework un `margin(L)` secco non va bene, e cosa fa
invece `classify_margins`?**
R: Il loop aperto ha piu' di un attraversamento a 0 dB: uno a bassa frequenza
generato dal canale di drift laterale (che contiene un integratore di posizione),
uno di corpo rigido attorno a sqrt(A6), e potenzialmente uno flessionale a
omega_BM. `margin()` restituisce il **piu' piccolo** dei margini fra tutti gli
attraversamenti, e nel caso del lanciatore pesca l'artefatto di drift, che non ha
alcun significato progettuale. `classify_margins` legge invece `allmargin` e
smista gli attraversamenti per **banda di frequenza**: sotto `w_drift` = artefatti
(scartati), fra `w_drift` e `w_flex` = corpo rigido (Rigid PM, Rigid GM per
gmdb > 0, Aero GM per gmdb < 0), sopra `w_flex` = modo flessionale. E' l'unico
modo di ottenere numeri confrontabili con i target della traccia (6 dB / 30
gradi). **Nota critica: `main_montecarlo.m` NON usa `classify_margins`** ed e'
rimasto alla lettura grezza -- e' il difetto principale dello script.

**D: I guadagni usati dal Monte Carlo sono quelli finali del progetto HM3?**
R: **No.** La riga 42 chiama `design_controller(Grigid, [], ...)`, cioe' progetta
sul plant **rigido** con attuatore **ideale**: sono i guadagni del Task 1
(Kp = 1.78, Kd = 0.44). Il progetto finale di HM3 e' quello ri-sintonizzato sul
loop pieno nel Task 2 (Kp = 1.73, Kd = 0.69), che e' quello che usano
`main_task2.m` e `main_task3.m`. Il notch invece e' quello giusto (Task 2). Il
commento in testa allo script, che parla di "Task-2 controller held fixed", e'
falso. Il Monte Carlo sta quindi valutando la robustezza della configurazione
che il Task 2 chiama esplicitamente "BEFORE re-tuning" e **scarta** (Rigid PM
14.6 gradi contro i 30 richiesti). Le sue statistiche non sono confrontabili con
la tabella del Task 3.

**D: Perche' i campioni casuali vengono estratti PRIMA del `parfor` e non
dentro?**
R: Per la riproducibilita'. Dentro un `parfor` ogni worker ha il proprio stream
di numeri casuali e l'assegnazione delle iterazioni ai worker non e'
deterministica: lo stesso seme darebbe risultati diversi al variare della
dimensione del pool, o fra esecuzione parallela e seriale. Estraendo i cinque
vettori di fattori sul client dopo un singolo `rng(seed)`, il corpo del loop
diventa una funzione deterministica dell'indice e il risultato e' identico con o
senza Parallel Computing Toolbox. Questo e' anche il motivo per cui `parfor`
degrada senza problemi a `for` seriale se il toolbox non c'e'.

**D: Perche' lo smorzamento del modo flessionale e' modellato lognormale e non
gaussiano?**
R: Perche' lo smorzamento e' una quantita' **strettamente positiva** (zeta_BM =
0.005 nominale, gia' vicino a zero). Una gaussiana con dispersione ampia
assegnerebbe probabilita' non nulla a valori negativi, cioe' a un modo
flessionale **auto-eccitato**: fisicamente assurdo e numericamente catastrofico.
La lognormale `exp(sigma*randn)` ha supporto (0, +Inf), mediana esatta 1 e la
giusta asimmetria (piu' facile che lo smorzamento reale sia 2x del nominale che
0.5x). Per mu_alpha e mu_c, che sono coefficienti con segno noto e dispersione
del 10 %, la gaussiana troncata va benissimo.

**D: Perche' la troncatura per rifiuto invece del clipping?**
R: Il clipping (`min(max(f,lo),hi)`) schiaccerebbe tutta la coda esterna sui due
estremi, creando due **atomi** di probabilita' in lo e hi. Il risultato non e'
una gaussiana troncata: e' una gaussiana troncata piu' due masse concentrate sui
bordi della scatola, che nel nostro caso sono proprio i punti piu' critici. La
P(violazione) risulterebbe sovrastimata e l'istogramma avrebbe due picchi
artificiali agli estremi. Il rejection sampling ri-estrae fino a cadere dentro il
supporto e produce la densita' gaussiana correttamente rinormalizzata.

**D: Che cosa aggiunge il Monte Carlo rispetto ai quattro vertici del Task 3?**
R: I vertici sono il caso peggiore **solo se** la metrica dipende in modo
monotono da ciascun parametro: per una funzione monotona su un ipercubo il minimo
cade in un vertice. Qui la monotonia e' plausibile per mu_alpha e mu_c (piu'
instabile, meno autorita' = peggio) ma non e' dimostrata, e per omega_BM e'
palesemente **falsa** (il notch ha un null stretto: l'attenuazione degrada da
entrambi i lati della frequenza nominale). Il Monte Carlo campiona l'**interno**
della scatola, aggiunge tre fattori che il Task 3 non disperde, e soprattutto
restituisce una **probabilita' di violazione** invece di un verdetto binario --
che e' l'unica cosa integrabile in un budget di affidabilita'. Non e' pero' un
certificato: con code troncate, i vertici hanno probabilita' ~2e-6 e con N = 1500
non vengono mai estratti, quindi il Monte Carlo **non sostituisce** l'analisi ai
vertici.

**D: Con N = 1500 quanto e' precisa la P(stabile) stimata?**
R: La stima e' una media di Bernoulli, quindi ha errore standard
SE = sqrt(p*(1-p)/N). Per p ~ 0.5 (caso peggiore) SE ~ 1.3 %; per p ~ 0.01, SE ~
0.26 %, cioe' un errore relativo del 26 % sulla probabilita' piccola -- e le
probabilita' piccole sono esattamente quelle che interessano in affidabilita'. Se
si osservano **zero** fallimenti su 1500, non si puo' concludere p = 0: la regola
del tre da' p <= 3/1500 = 0.2 % al 95 % di confidenza. Per stimare probabilita'
dell'ordine di 1e-4 servono N ~ 1e6, oppure tecniche di riduzione della varianza
(importance sampling, subset simulation) o un metodo deterministico
(mu-analisi/LMI) che copra il caso peggiore per costruzione.

**D: Perche' il notch non viene ri-sintonizzato quando omega_BM disperde?**
R: Perche' il controllore in volo e' congelato: la vera omega_BM non e'
misurabile in tempo reale, e' esattamente l'incertezza che si sta quantificando.
Ri-sintonizzare a ogni estrazione risponderebbe a una domanda diversa ("quanto e'
buono il metodo di progetto") anziche' a quella giusta ("quanto e' robusto questo
controllore"). La riga 15 del sorgente lo dichiara e ne trae la conseguenza:
proprio perche' il notch e' profondo (zeta_N = 0.002) e strettissimo, omega_BM
diventa il driver dominante dell'instabilita' nel Monte Carlo.
