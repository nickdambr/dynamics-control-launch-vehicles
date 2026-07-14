# HM3/main_task2.m

## Ruolo del file nel progetto

E' il **secondo entry point** di HM3 e il piu' ricco dei tre. Implementa il Task 2
della traccia ("Full LV Model"): estendere il progetto del Task 1 al modello
completo, aggiungendo (a) la **dinamica dell'attuatore TVC** di Eq. (3), (b) il
**ritardo puro di 20 ms**, (c) il **primo modo di bending** con la sua
contaminazione della misura inerziale (Eq. 2), e infine inserire il **filtro
notch** e ri-verificare la stabilita'.

Lo script e' organizzato in quattro step espliciti (i banner `%%` sono la
scomposizione logica voluta dall'autore):

    Step A : ricostruisce il PD del Task 1 (attuatore ideale)
    Step B : monta TVC + ritardo + bending SENZA filtro -> l'anello esplode
    Step C : trade study su QUATTRO filtri di bending -> sceglie il notch profondo
             + RI-SINTONIZZA il PD sull'anello completo
    Step D : detuning di omega_BM (+/-10%) a filtri fissi

Il fatto tecnico che genera tutto il resto e' questo: il modo di bending ha
**zeta_BM = 0.005**, cioe' un fattore di merito Q = 1/(2*zeta_BM) = **100**, una
risonanza da +40 dB. Il modo entra nell'anello per due strade contemporaneamente:

- **in avanti**: il TVC forza il bending, riga eta_ddot del plant, con il termine
  -phi_tvc*Tc*delta;
- **all'indietro**: la piattaforma inerziale non misura theta ma
  theta_m = theta + sigma_ins*eta (Eq. 2, sigma_ins = 0.178 rad/m), quindi il
  bending **rientra nella retroazione**.

Il prodotto di questi due percorsi chiude un anello parassita attorno alla
risonanza. Il codice misura che a omega_BM il guadagno d'anello vale **+29 dB**:
ventinove decibel *sopra* la soglia di instabilita', con fase arbitraria.
Senza filtro l'anello chiuso e' instabile (polo a Re = +1.61).

Nota: il bending contamina anche z_m e zdot_m (con -phi_ins = -0.8), ma quei
canali sono chiusi con guadagni da 1e-3, quindi il percorso che conta davvero e'
sigma_ins sul canale di assetto.

**Refuso nel sorgente**: il commento di header alle righe 7-8 parla della
"+39 dB bending resonance". Il valore reale, stampato dal codice stesso alla riga
33, e' **+29.0 dB**. Il README riporta correttamente 29 dB.

---

## `%% Model and parameters` + `%% Step A` (righe 19-26)

```matlab
p = load_hw3_params();
Grigid = build_plant_rigid(p);
[K, mR] = design_controller(Grigid, []);
[~, Trigid] = assemble_loop(Grigid, K);
```

- Riga 17: `warning('off','Control:analysis:MarginUnstable')`. `margin`/`allmargin`
  emettono un warning ogni volta che valutano un anello instabile in anello
  aperto -- che qui e' **sempre**, per costruzione. Senza questo silenziamento la
  console verrebbe sommersa. E' legittimo, non e' nascondere un problema.
- Righe 23-25: **Step A non riusa un file di risultati: ri-esegue il tuner del
  Task 1**. `design_controller(Grigid, [])` e' deterministico (`fminsearch` da un
  seed fisso), quindi riproduce esattamente Kp = 1.7845, Kd = 0.4433. Il commento
  dice "reused from Task 1" ma tecnicamente e' **ricalcolato**. Costa un tuner run
  in piu' ma rende lo script autosufficiente.
- Riga 26: `Trigid` serve solo alla fine, per il confronto rigido-vs-completo
  (figura f3).

---

## `%% Step B - full plant + TVC + delay, NO bending filter` (righe 28-34)

```matlab
Gfull = build_plant_full(p, 'ins');
Wtvc  = build_tvc(p, 3);          % 2nd-order TVC + Pade(20 ms, order 3)
[Lb, Tb] = assemble_loop(Gfull, K, Wtvc);
```

- Riga 29: `'ins'` seleziona la matrice di uscita di **Eq. (2)**, quella
  contaminata dal bending. L'alternativa `'true'` (misure pulite) esiste nel
  modulo ma non e' usata qui: usarla farebbe sparire il problema, che e' proprio
  cio' che si vuole studiare.
- Riga 30: `build_tvc(p, 3)` costruisce

      W_TVC(s) = wTVC^2 / (s^2 + 2*zTVC*wTVC*s + wTVC^2) * Pade(tau, 3)

  con wTVC = 70 rad/s, zTVC = 0.7, tau = 0.020 s. **Perche' Pade ordine 3**: il
  ritardo puro exp(-tau*s) non e' razionale e non si puo' mettere in un
  `ss`/`tf`. L'approssimante di Pade di ordine 3 riproduce fedelmente la fase
  fino a circa 3/tau ~ 150 rad/s, quindi copre bene sia il crossover rigido
  (~3 rad/s) sia la zona di bending (18.9 rad/s). Un ordine 1 sarebbe troppo
  grossolano a omega_BM. **Costo**: Pade introduce zeri a destra (non-minimum
  phase), che e' esattamente il comportamento di un ritardo, quindi e' fedele.
- Righe 32-34: output reale

      |L(omega_BM)| = 29.0 dB  -> closed-loop stable: 0 (max Re pole = 1.61)

  Il modo di bending non e' ne' gain-stabilizzato (sarebbe |L| << 0 dB) ne'
  phase-stabilizzato (la fase a omega_BM e' quella che capita): l'anello e'
  semplicemente instabile.

> **Possibile domanda d'esame** -- I margini di bassa frequenza (Aero GM, Rigid
> PM) sopravvivono all'aggiunta di TVC e ritardo?
> *Risposta:* Si', quasi intatti. Nella tabella di Step C la riga "no filter" da'
> Aero GM = 6.25 dB e Rigid PM = 26.1 deg: partivano da 6.00 e 30.0. Il motivo e'
> che il TVC ha wTVC = 70 rad/s, ben oltre un decennio sopra il crossover rigido
> (~2.5 rad/s), e 20 ms di ritardo a 2.5 rad/s valgono solo 0.020*2.5 = 0.05 rad
> = 3 deg. **Il problema del Task 2 non e' l'attuatore in banda: e' il bending.**
> L'attuatore diventa un problema solo *dopo*, in combinazione con il ritardo di
> fase del notch (vedi Step C).

---

## `%% Step C - bending filter trade study` (righe 36-129)

Quattro candidati, valutati **a guadagni PD fissi** (quelli del Task 1). La
griglia comune (righe 47-49):

    wx_grid = wBM + (-4:2:4) = [14.9  16.9  18.9  20.9  22.9]  rad/s
    zN_grid = 0.10 : 0.05 : 0.30                               (5 valori)
    zD_grid = 0.40 : 0.10 : 0.60                               (3 valori)
    -> 5 * 5 * 3 = 75 combinazioni

Le tre griglie sono **esattamente i range suggeriti dalla traccia**
("zeta_x,N ~ 0.1-0.3, zeta_x,D ~ 0.4-0.6, omega_x around the first bending mode
frequency +/- 4 rad/s").

### La matematica di Eq. (4)

`build_notch_filter(wx, zN, zD, sgn)` costruisce

    Hx(s) = (s^2 + sgn*2*zN*wx*s + wx^2) / (s^2 + 2*zD*wx*s + wx^2)

Due osservazioni importanti:

1. **La traccia ha un refuso**: Eq. (4) e' stampata come
   `(s^2 - 2*zeta_xN*omega_x + omega_x^2)/(...)`, cioe' senza la `s` sul termine
   centrale del numeratore -- dimensionalmente impossibile. Il codice ripristina
   la `s` mancante, che e' l'unica lettura sensata.
2. **Il segno del numeratore e' il vero interruttore di progetto.**
   - `sgn = -1` (Eq. 4 **come stampata**): il numeratore ha zeri a **parte reale
     positiva** -> filtro **non-minimum-phase**. Il modulo e' identico a quello
     del corrispondente minimum-phase, ma la fase e' specchiata: il filtro
     **ruota** la curva senza abbassarla. E' un **phase shaper**: e' esattamente
     il "targeted phase shift at the bending mode frequency in the Nichols plot
     (left-right shift)" di cui parla la traccia.
   - `sgn = +1`: notch **minimum-phase** simmetrico, che **abbassa** il modulo.
     E' un **gain stabiliser**.

   Attenuazione al centro (per entrambi i segni), valutando |Hx(j*wx)|: i termini
   in s^2 e wx^2 si cancellano e resta

       |Hx(j*wx)| = (2*zN*wx*wx)/(2*zD*wx*wx) = zN/zD

   -> profondita' in dB = **20*log10(zN/zD)**. Con i range della traccia il
   massimo che si ottiene e' 20*log10(0.1/0.6) = **-15.6 dB**: molto poco contro
   una risonanza a +29 dB. **Da qui nasce tutta la logica del trade**: i valori
   suggeriti dalla traccia non possono gain-stabilizzare il modo, quindi o si
   phase-stabilizza (usando sgn = -1) o si esce dai range.

### C-LL: lead-lag di Eq. 4 da sola (righe 51-71)

- Righe 54-67: triplo ciclo sulle 75 combinazioni, `sgn = -1`. Per ognuna chiude
  l'anello e registra `max(real(pole(Tc_)))`. Il criterio di selezione (riga 61)
  e' **il minimo del massimo polo reale**, cioe' "il meno instabile" / il piu'
  smorzato. **Non e' una selezione sui margini.**
- Output reale:

      C-LL alone: 45/75 guideline candidates stabilise the loop; least unstable
                  (wx=18.9, zN=0.30, zD=0.60) has max Re(pole) = -0.06

**Contraddizione fra commento e codice, da segnalare**: il commento alla riga 39
dice `"Alone it never stabilises"`. E' **falso**: 45 candidati su 75 stabilizzano
l'anello, e il migliore ha max Re(pole) = -0.06, cioe' e' **stabile** (appena).
Anche la legenda della figura (riga 238) lo etichetta "marginal", che e' piu'
corretto. Il commento va letto come "non stabilizza *bene*", non "non stabilizza".

Riga di tabella corrispondente:

      C-LL  lead-lag only |   5.96    11.4     8.02      23.0      15 | stabile

Il dato chiave e' **|L(omega_BM)| = +23.0 dB**: il modo e' ancora **sopra** 0 dB.
La lead-lag lo ha attenuato solo di 6 dB (= 20*log10(0.30/0.60)), esattamente
come predice la formula. Sopravvive quindi **solo per fase**: il lobo di bending
passa dalla parte giusta del punto critico. E' phase stabilisation pura. Ma il
delay margin e' **15 ms**, cioe' **inferiore ai 20 ms di ritardo gia' modellati**:
il progetto e' sul filo del rasoio e non e' difendibile.

### C-N: notch profondo minimum-phase (righe 73-78) -- il progetto ritenuto

```matlab
notch.wx  = p.wBM;    notch.zN = 0.002;
notch.zD  = 0.7;      notch.sgn = +1;
```

- Profondita': 20*log10(0.002/0.7) = **-50.9 dB** (verificato). Contro i +29 dB
  della risonanza -> |L(omega_BM)| = 29 - 51 = **-21.9 dB**: il modo e'
  **gain-stabilizzato**, cioe' semplicemente non ha piu' abbastanza guadagno per
  destabilizzare, qualunque sia la sua fase.
- **Deviazione dichiarata dalla traccia**: zN = 0.002 e zD = 0.7 sono **fuori dai
  range suggeriti** (zN 0.1-0.3, zD 0.4-0.6). E' una scelta consapevole -- come
  mostrato sopra, dentro quei range la profondita' massima e' -15.6 dB, che non
  basta. Il progetto ritenuto **cambia dispositivo**: non e' piu' la lead-lag di
  Eq. (4), e' un notch stretto e profondo. All'orale va difeso cosi', non
  nascosto.

### C-T: tripletta di notch (righe 80-83)

Tre notch identici a 0.9, 1.0 e 1.1 volte omega_BM. E' la ricetta classica per
coprire l'incertezza sulla frequenza del modo. Risultato:

      C-T   notch triplet |    NaN    -7.3      NaN     -55.8    2469 | INSTABILE

- |L(omega_BM)| = -55.8 dB: il bending e' sepolto.
- Ma **Rigid PM = -7.3 deg**: l'anello e' instabile. **Ogni notch costa fase alle
  frequenze sotto il suo centro**, e tre notch impilati costano circa il triplo:
  al crossover rigido (~2.5 rad/s) mangiano ~30 deg di margine di fase, che e'
  tutto quello che c'era. **E' la lezione piu' istruttiva del trade**: i notch non
  sono gratis, e sovra-filtrare uccide il corpo rigido.
- Il **DM = 2469 ms** in quella riga e' un numero **privo di senso**: `allmargin`
  restituisce margini anche su anelli instabili, dove non significano nulla. Vale
  lo stesso per la riga "no filter" (DM = 41 ms). **L'unica colonna che decide e'
  `stable`.** Il codice stampa i margini comunque, il che e' onesto ma va letto
  con questa avvertenza.

### C-NLL: notch + lead-lag (righe 85-106)

- Secondo sweep sulle stesse 75 combinazioni, questa volta in **serie al notch
  profondo** (`Wtvc*Hn*Hc`). Criterio di selezione (riga 97): **massimo delay
  margin**, valutato solo sui candidati stabili.
- Output: 67/75 stabilizzano; il migliore e' wx = 22.9, zN = 0.10, zD = 0.40.
- Riga di tabella: Aero 5.72 / PM 8.0 / |L(wBM)| -28.1 dB / DM 54 ms. Margini
  nominali **piu' magri** del solo notch, ma vedremo in Step D che e' molto piu'
  robusto al detuning.

### Tabella comparativa (righe 108-129)

- Riga 119: le bande passate a `classify_margins`:

      w_drift   = 0.3*sqrt(A6) = 0.55 rad/s   (sotto: artefatti di deriva)
      w_flex    = 0.6*wBM      = 11.34 rad/s  (confine rigido/flessibile)
      w_flex_hi = 1.5*wBM      = 28.35 rad/s
      w_bending = wBM          = 18.90 rad/s

  **Limite da segnalare**: il Rigid GM del progetto finale cade a **11.11 rad/s**,
  contro un confine `w_flex` a **11.34 rad/s**. Sta *dentro* la banda rigida per
  soli 0.23 rad/s. E' fragile: una piccola variazione di guadagno o di omega_BM lo
  riclassificherebbe come margine flessibile e i numeri stampati cambierebbero
  categoria (non valore). La scelta di w_flex = 0.6*wBM e' euristica e il codice
  non la giustifica.

Tabella completa dalla run reale:

    candidate           |  AeroGM RigidPM  RigidGM  |L(wBM)|  DM[ms] | stable
    B     no filter     |   6.25    26.1      NaN      29.0      41 | 0
    C-LL  lead-lag only |   5.96    11.4     8.02      23.0      15 | 1
    C-N   deep notch    |   6.08    14.6    10.43     -21.9      98 | 1
    C-T   notch triplet |    NaN    -7.3      NaN     -55.8    2469 | 0
    C-NLL notch+leadlag |   5.72     8.0     7.41     -28.1      54 | 1

Si legge cosi': **tutti** i candidati che sopravvivono pagano il bending con il
margine di fase rigido. Nessuno mantiene i 30 deg del Task 1. Ecco perche' serve
lo step successivo.

> **Possibile domanda d'esame** -- Gain stabilisation o phase stabilisation: qual
> e' la differenza e perche' qui si sceglie la prima?
> *Risposta:* **Gain stabilisation** = si abbassa |L| sotto 0 dB alla frequenza
> del modo, cosi' la fase li' non conta piu': il lobo e' troppo piccolo per
> circondare il punto critico. **Phase stabilisation** = si lascia |L| > 0 dB ma
> si ruota la fase in modo che il lobo passi dalla parte sicura del punto critico.
> La seconda e' piu' efficiente (non serve un filtro profondo, quindi si paga meno
> fase in banda) ma e' **intrinsecamente fragile**: dipende dal conoscere la fase
> del modo, che dipende da omega_BM, da zeta_BM, dai coefficienti modali. La prima
> e' robusta in fase ma costosa in fase-in-banda. Qui si sceglie il **gain
> stabilisation** perche' il margine risultante (-18 dB finale) e' ampiamente sopra
> i 12 dB tipicamente richiesti e perche' la lead-lag phase-stabilizzante da sola
> lascia un delay margin di 15 ms, inaccettabile.

---

## `%% Step C decision - deep notch, then RE-TUNE the PD` (righe 131-161)

```matlab
Wfull = Wtvc * Hn;
[Lb1, Tb1] = assemble_loop(Gfull, Ktask1, Wfull);   % PRIMA
mB = classify_margins(Lb1, bands{:});
...
[K, mF] = design_controller(Gfull, Wfull, 'w_flex', 0.6*p.wBM, ...);  % DOPO
```

### Prima del re-tune (righe 142-149)

      Kp=1.784 Kd=0.443 | Aero |GM|=6.08 dB  Rigid PM=14.6 deg @2.59  DM=98 ms

Il margine di fase rigido **crolla da 30.0 a 14.6 gradi** e il delay margin scende
a 98 ms, cioe' **sotto** il requisito tipico di 100 ms citato nel Task 1.

### Da dove viene esattamente la fase persa

Il commento (righe 132-139) attribuisce il crollo a "the actuator + transport
delay + notch phase lag". Ho misurato i tre contributi separatamente alla
frequenza di crossover di quel caso (wc = 2.59 rad/s):

    TVC 2nd-order (wTVC=70, zTVC=0.7) :  -2.97 deg
    Pade(20 ms, ordine 3)             :  -2.97 deg
    notch profondo (zD=0.7 @ 18.9)    : -11.03 deg
    -------------------------------------------------
    totale                            : -16.97 deg

**Il colpevole principale e' il notch, non l'attuatore.** Da solo vale 11 dei
17 gradi persi, cioe' i due terzi. Attuatore e ritardo insieme valgono appena
6 gradi. Il commento del sorgente elenca i tre in ordine, il che suggerisce che
l'attuatore sia il problema; **quantitativamente non lo e'**.

Il motivo e' geometrico: il denominatore del notch e' un secondo ordine a
zD = 0.7 centrato a 18.9 rad/s, e un secondo ordine *inizia a ruotare la fase gia'
una decade sotto* il suo centro. A 2.6 rad/s (una decade sotto) ha gia' speso
11 gradi. Un notch **piu' stretto in denominatore** (zD piu' basso) costerebbe
meno fase in banda -- ma renderebbe il null ancora piu' fragile al detuning.

### Il re-tune (righe 151-161)

      PD design (full loop): Kp_th=1.7318 Kd_th=0.6867
        Aero |GM|=6.00 dB @0.54  Rigid PM=30.0 deg @3.17  DM=165 ms | stable: 1
        |L(omega_BM)| = -18.2 dB
        Rigid GM = 7.56 dB @ 11.11 rad/s

- **Kp scende leggermente** (1.784 -> 1.732), **Kd sale del 55%** (0.443 ->
  0.687). E' la mossa attesa: il termine derivativo Kd*s contribuisce fase
  **positiva** (+90 deg asintotici), quindi alzarlo e' il modo diretto di
  recuperare il ritardo di fase introdotto dal notch. Il crossover si sposta in
  alto (2.45 -> 3.17 rad/s) perche' l'anello ha piu' guadagno ad alta frequenza.
- **Il prezzo, che il codice non commenta:** |L(omega_BM)| peggiora da **-21.9 dB
  a -18.2 dB**. Alzando Kd si alza il guadagno d'anello *anche* a omega_BM (il
  derivativo cresce con la frequenza), quindi il re-tune **baratta 3.7 dB di
  margine di bending per 15 gradi di margine di fase**. Resta comunque un margine
  di bending di 18 dB, sopra i 12 dB usualmente richiesti per un modo
  gain-stabilizzato.

### Deviazione dalla guideline della traccia

La traccia, al punto 2 delle guidelines, mette una nota a pie' di pagina esplicita:
*"In most cases, retuning is not necessary."* **Questo script ri-sintonizza
comunque.** La deviazione e' giustificata dai numeri (PM 14.6 -> 30.0, DM 98 ->
165 ms), ma va dichiarata all'orale come **scelta deliberata contro il
suggerimento**, non come applicazione della traccia. Nota inoltre che la nota
della traccia si riferisce al punto 2 (TVC + ritardo, **senza** notch): e li'
infatti il retune *non* serve (Aero 6.25 / PM 26.1). E' il **notch** ad aver reso
il retune necessario, e il notch e' il punto 3.

---

## `%% Step D - sensitivity to the bending-frequency knowledge` (righe 163-179)

```matlab
scales = 0.90:0.05:1.10;
for i = [3 5]                       % C-N (retained) e C-NLL
    for sc = scales
        ps = load_hw3_params();  ps.wBM = sc*ps.wBM;
        Gs = build_plant_full(ps, 'ins');
        [~, Ts] = assemble_loop(Gs, K, Wtvc*cand{i,2});
        fprintf('  %5d', isstable(Ts));
```

I **filtri restano fissi** (progettati per omega_BM = 18.9), mentre il **modo vero**
si sposta. E' il test di robustezza giusto: nella realta' omega_BM e' noto con
qualche percento di incertezza (dipende da riempimento serbatoi, temperatura,
modellazione FEM).

Risultato reale:

    candidate           |  x0.90  x0.95  x1.00  x1.05  x1.10
    C-N   deep notch    |      1      1      1      0      0
    C-NLL notch+leadlag |      1      1      1      1      1

Il notch profondo **tollera -10% ma va instabile gia' a +5%**. La combinazione
notch+lead-lag tiene tutta la banda.

### La spiegazione del README e' sbagliata -- ecco quella corretta

Il README attribuisce l'asimmetria a "the asymmetric edge of its 0.08 rad/s
null", cioe' alla forma del null di attenuazione. **Non regge**, e si verifica in
due righe. L'attenuazione del notch al modo spostato:

    |Hn| a 0.90*wBM = -16.5 dB
    |Hn| a 1.05*wBM = -23.1 dB      <-- attenua DI PIU' a +5% che a -10%
    |Hn| a 1.10*wBM = -17.4 dB

Il notch attenua **di piu'** a +5% che a -10%, eppure e' -10% a sopravvivere.
L'attenuazione non puo' essere la causa. Misurando il guadagno d'anello al modo
vero:

     sc  | wBM_true | |Hn(wBM_true)| | |L(wBM_true)| | stabile
    0.90 |  17.01   |   -16.5 dB     |   +17.1 dB    |   si
    0.95 |  17.95   |   -22.7 dB     |   +10.4 dB    |   si
    1.00 |  18.90   |   -50.9 dB     |   -18.2 dB    |   si
    1.05 |  19.84   |   -23.1 dB     |    +9.1 dB    |   NO
    1.10 |  20.79   |   -17.4 dB     |   +14.4 dB    |   NO

Il fatto decisivo: **appena il modo si sposta, il lobo di bending torna sopra
0 dB in tutti i casi** (da +9 a +17 dB). La gain stabilisation **evapora
completamente** con qualunque detuning: il null e' largo pochi centesimi di
rad/s, il modo esce e la risonanza da +29 dB riemerge quasi intera. Da quel
momento in poi la sopravvivenza dipende **solo dalla fase** del lobo, cioe' da
quale lato del punto critico passa.

E la fase e' asimmetrica attorno al notch: **sotto** il centro il notch
contribuisce una rotazione favorevole, **sopra** una sfavorevole. Ecco perche' a
+5% l'anello cade pur avendo *meno* guadagno al modo (9.1 dB) di quanto ne abbia a
-5%... anzi, di quanto ne abbia a 0.95 (10.4 dB), che invece e' stabile. E' una
dimostrazione pulita che qui **conta la fase, non il modulo**.

**Conclusione onesta**: il progetto ritenuto e' gain-stabilizzato **solo al punto
nominale**. Fuori dal nominale sopravvive per fortuna (phase stabilisation
accidentale), e da un lato solo. E' il vero punto debole del Task 2, ed e' anche
il motivo per cui la variante C-NLL, pur avendo margini nominali peggiori, e'
l'unica robusta.

### Due incoerenze minori in Step D

1. **I guadagni usati sono quelli RI-SINTONIZZATI.** A riga 153 la variabile `K`
   viene **sovrascritta** con il design del Task 2. Step D (riga 175) usa quindi
   il `K` nuovo. Ma il partner lead-lag di C-NLL (`bestC`, righe 88-102) era stato
   **selezionato con i guadagni del Task 1**. Step D confronta quindi un filtro
   ottimizzato per un set di guadagni con un set di guadagni diverso.
2. **La tabella del trade (righe 122-129) usa i guadagni del Task 1**, mentre il
   progetto finale usa quelli ri-sintonizzati. La riga "C-N deep notch" della
   tabella (PM 14.6) **non e' il progetto finale** (PM 30.0). Lo script lo dice
   ("BEFORE re-tuning"), quindi non e' fuorviante -- ma le righe C-LL, C-T e
   C-NLL non vengono mai ri-valutate con i guadagni finali. Se all'orale si mostra
   la tabella, va detto che e' una **fotografia a guadagni Task-1**.

---

## `%% Time response to a wind gust` (righe 181-189)

Stessa raffica del Task 1 (`Tend=80`, 1-cosine severa, Vg = 6.38 m/s). Risultati:

      peak |theta| = 0.231 deg, |z| = 2.27 m, |delta| = 0.513 deg
      peak |alpha| = 0.565 deg -> peak qbar*alpha = 45.8 kPa deg

Praticamente **identici al Task 1** (theta 0.261 / z 2.27 / delta 0.528, e
alpha 0.577 deg -> 46.8 kPa*deg). Ed e' il risultato giusto: bending, attuatore
e ritardo vivono tutti a frequenze molto sopra la banda della raffica (la raffica
dura 3 s, cioe' contenuto spettrale sotto ~2 rad/s). La risposta temporale a
bassa frequenza **non li vede**. E' esattamente per questo che il bending e' un
problema di **stabilita'** e non di **prestazione**: non lo trovi guardando le
time history, lo trovi solo sul Nichols.

**Il segno di `alpha_w` (righe 188-189).** L'incidenza ricostruita da
`simulate_gust_response.m` (riga 29) e'

    alpha = theta + zdot/V - alpha_w

con il **meno**, che e' l'unico segno coerente con l'Eq. (1): la colonna del
vento del plant e' `Bw = [0; -a1*V; 0; -A6; 0; 0]`
(`build_plant_full.m` riga 28), cioe' i coefficienti di `theta` cambiati di
segno -- il che *e'* la definizione `alpha = theta + zdot/V - alpha_w`. Fisicamente:
l'incidenza si misura rispetto alla velocita' **relativa all'aria**.

La conseguenza e' un punto d'orale: il picco di incidenza (**0.565 deg**)
**supera** quello del vento da solo (0.390 deg). Per tenere l'assetto il loop
becca il muso **dentro** il vento relativo, e quel contributo si **somma** a
quello del vento. Un puro attitude-hold e' quindi **load-aggravating**, non
load-relieving -- e con `A6 > 0` (centro di pressione davanti al baricentro,
momento aerodinamico divergente) non esiste nessuna stabilita' a banderuola che
allevi il carico da sola.

> **Nota storica.** Fino a poco fa `simulate_gust_response.m` usava il segno
> **piu'**, dando un picco di 0.253 deg (20.5 kPa*deg) -- *sotto* il vento da
> solo, e quindi un'apparente azione di load relief. Il plant e' sempre stato
> corretto: il difetto stava solo nella ricostruzione a valle di `alpha`. La
> correzione **non tocca** margini, Nichols, notch, Step C/D o stabilita'
> (`alpha_w` e' un disturbo, non entra in `L(s)`), e nemmeno i picchi di `theta`,
> `z`, `delta`: cambia **solo** la metrica di carico.

---

## `%% Figures` (righe 191-276)

Cinque figure.

- **f1 `nichols`** (righe 196-200): il progetto ritenuto, via `plot_nichols_lv`.
  Convenzione degli appunti del corso: punto critico a **(-180 deg, 0 dB)**, da
  1 + L = 0, finestra `xlim [-720 0]`. (D'Antuono Fig. 3.2 mostra la stessa carta
  rietichettata di +360; il codice usava quella rietichettatura fino a poco fa ed
  e' stato riallineato al corso -- pura rietichettatura dell'asse, margini
  invariati.)
- **f5 `retune`** (righe 205-214): **la figura piu' importante del task**. Sovrappone
  l'anello con i guadagni Task-1 (`Lb1`) e quello ri-sintonizzato (`Lc`), con
  `PhaseMatching` agganciato alla frequenza di crossover del design finale e
  `PhaseMatchingValue = -180` (riga 208) -- stessa convenzione del corso. Mostra
  visivamente il recupero del margine di fase.
- **f4 `nichols_trade`** (righe 222-243): il trade a tre (nessun filtro / lead-lag
  / notch profondo), disegnato **a mano** invece che con `nicholsplot`, perche'
  serve applicare a tutte e tre le curve **lo stesso** shift di fase. La riga 226:

      sh0 = 360*round((-180 - interp1(wv, ph0, mB.rigidPM_w))/360);

  Lo shift e' **un multiplo esatto di 360 deg** -> e' una rietichettatura del ramo
  di fase, non una rotazione. Le tre curve restano confrontabili: le regioni
  rigide si sovrappongono attorno a -180 deg (e' quanto dichiara il commento
  riscritto alle righe 217-221) e i lobi di bending sono direttamente
  confrontabili sullo stesso foglio. Il punto critico e' marcato con una croce
  rossa a -180 (riga 239); il gemello a +180 di riga 238 cade fuori dalla
  finestra `xlim [-720 0]` di riga 240. Bonus della convenzione: la griglia M/N
  di `ngrid` (riga 223) vive nativamente in [-360, 0], quindi ora la regione
  rigida della figura casca **dentro** la griglia (con la vecchia
  rietichettatura a +180 finiva fuori). Stesso trucco in
  `plot_nichols_lv` (riga 49) e in `main_task3.m` (riga 81).
- **f2 `gust_response`** (righe 246-256): theta, z, zdot, delta -- le quattro
  variabili chieste dalla traccia.
- **f3 `comparison_rigid_vs_full`** (righe 259-267): rigido vs completo, che e' il
  "compare the results" richiesto dalla traccia.
- Export (righe 269-280): identico al Task 1, PNG a 200 dpi in `HM3/figures/`,
  prefisso `task2_`.

---

## Possibili domande d'esame

**D: Perche' il modo di bending destabilizza l'anello? Descrivi il percorso.**
R: E' un anello parassita chiuso attorno alla risonanza. In avanti, il TVC forza
il modo: nella riga eta_ddot del plant compare il termine -phi_tvc*Tc*delta.
All'indietro, la piattaforma inerziale non misura l'assetto vero ma
theta_m = theta + sigma_ins*eta (Eq. 2, sigma_ins = 0.178 rad/m): il bending
**rientra nel controllore**, che lo riamplifica e lo rimanda al TVC. Con
zeta_BM = 0.005 il Q vale 100 (+40 dB di risonanza), e il guadagno d'anello a
omega_BM risulta **+29 dB**: molto sopra 0 dB, con fase non controllata ->
encirclement sbagliato -> instabilita' (polo a Re = +1.61). Se si usasse la
misura pulita (`build_plant_full(p,'true')`) il problema sparirebbe: e' la
**contaminazione della misura**, non il modo in se', a chiudere l'anello.

**D: Perche' non basta la lead-lag di Eq. (4), che pure e' quella suggerita dalla
traccia?**
R: Perche' la sua attenuazione al centro vale **zN/zD**, cioe' al massimo
20*log10(0.1/0.6) = -15.6 dB dentro i range suggeriti -- contro una risonanza di
+29 dB. Non puo' gain-stabilizzare. Puo' solo **phase-stabilizzare**: con
sgn = -1 gli zeri sono a destra, il filtro ruota la fase senza abbassare il
modulo, ed e' il "left-right shift" di cui parla la traccia. Nel run 45/75
candidati effettivamente stabilizzano cosi'. Ma il migliore lascia
|L(omega_BM)| = +23 dB e un **delay margin di 15 ms**, cioe' inferiore ai 20 ms
di ritardo gia' presenti nel modello. Non e' un progetto difendibile. (Nota: il
commento nel sorgente dice che la lead-lag "non stabilizza mai" -- e' sbagliato,
45 su 75 stabilizzano.)

**D: La tripletta di notch attenua il bending di -56 dB. Perche' e' il candidato
peggiore?**
R: Perche' e' **instabile**: Rigid PM = -7.3 deg. Ogni notch spende fase alle
frequenze *sotto* il proprio centro (il denominatore e' un secondo ordine che
inizia a ruotare gia' una decade prima), e tre notch impilati spendono circa il
triplo. Al crossover rigido (~2.5 rad/s) si mangiano ~30 gradi, cioe' tutto il
margine di fase disponibile. E' la dimostrazione che **i notch non sono gratis**:
si compra margine di bending pagando in margine di fase rigido, e la tripletta
spende piu' di quanto si possa permettere. Attenzione anche al DM = 2469 ms
stampato su quella riga: e' un numero senza senso, perche' `allmargin` calcola
margini anche su anelli instabili. L'unica colonna che decide e' `stable`.

**D: Il re-tune del PD contraddice la traccia, che dice "in most cases retuning is
not necessary". Come lo difendi?**
R: Distinguendo i due punti della guideline. Al **punto 2** (TVC + ritardo, senza
notch) il retune davvero **non serve**: i margini restano 6.25 dB / 26.1 deg,
praticamente quelli del Task 1, perche' wTVC = 70 rad/s e' oltre una decade sopra
il crossover e 20 ms di ritardo a 2.5 rad/s costano solo 3 gradi. E' il **punto 3**,
cioe' l'inserimento del **notch**, a rendere necessario il retune: il notch
profondo costa da solo **11 gradi** di fase al crossover (misurati), e insieme a
attuatore e ritardo porta il margine di fase da 30 a 14.6 gradi, con delay margin
sotto i 100 ms. Alzando Kd del 55% (0.443 -> 0.687) si recuperano i 30 gradi. La
nota della traccia si riferisce all'attuatore, non al filtro.

**D: Il notch profondo tollera -10% di detuning su omega_BM ma cade a +5%. Perche'
questa asimmetria?**
R: Non e' una questione di profondita' del null, come si potrebbe pensare (e come
sostiene il README): il notch attenua **di piu'** a +5% (-23.1 dB) che a -10%
(-16.5 dB). Il punto e' che **appena il modo si sposta, la gain stabilisation
evapora del tutto**: il lobo di bending torna sopra 0 dB in *tutti* i casi
detunati (da +9 a +17 dB di guadagno d'anello al modo). Il null e' largo pochi
centesimi di rad/s e il modo ne esce subito. Da li' in poi la stabilita' dipende
**solo dalla fase** del lobo, cioe' da che lato del punto critico passa -- e la
fase del notch e' asimmetrica attorno al proprio centro (favorevole sotto,
sfavorevole sopra). Prova decisiva: a 0.95*omega_BM il guadagno al modo e' 10.4 dB
e l'anello e' **stabile**; a 1.05*omega_BM e' 9.1 dB (cioe' **meno**) e l'anello e'
**instabile**. Il modulo non spiega nulla, la fase spiega tutto.

**D: Perche' la risposta alla raffica del Task 2 e' praticamente identica a quella
del Task 1?**
R: Perche' bending (18.9 rad/s), attuatore (70 rad/s) e ritardo (20 ms) vivono
tutti **molto sopra** la banda del disturbo: la raffica 1-cosine dura 3 s, quindi
ha contenuto spettrale sotto ~2 rad/s. La risposta temporale a bassa frequenza
semplicemente non li vede. E' il motivo per cui il bending e' un problema di
**stabilita'** e non di **prestazione**: non lo si scopre guardando le time
history, lo si scopre solo sul Nichols. Chiunque validasse il Task 2 solo con
`lsim` non si accorgerebbe di nulla.

**D: Perche' Pade di ordine 3 e non 1?**
R: Perche' l'approssimante deve essere fedele **fino a omega_BM = 18.9 rad/s**, non
solo al crossover rigido. Pade di ordine n riproduce bene la fase di exp(-tau*s)
fino a circa n/tau: con tau = 20 ms, l'ordine 1 e' affidabile fino a ~50 rad/s in
modulo ma sbaglia gia' la fase molto prima, mentre l'ordine 3 copre bene fino a
~150 rad/s. Dato che tutta la decisione sul filtro si gioca sulla **fase a
omega_BM**, un Pade grossolano falserebbe proprio il numero che conta.
