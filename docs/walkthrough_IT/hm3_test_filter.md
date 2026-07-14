# HM3/tests/hm3FilterTest.m

## Ruolo del file nel progetto

E' la suite `matlab.unittest` che congela i **due blocchi di compensazione**
inseriti in serie fra il controllore PD e il plant: `build_tvc.m` (attuatore TVC
del secondo ordine + ritardo di trasporto approssimato con Pade -- Eq. 3 della
traccia) e `build_notch_filter.m` (la sezione lead-lag / notch -- Eq. 4). Sono
i due elementi che nel Task 2 trasformano un problema risolto (il PD rigido del
Task 1) in un problema difficile: l'attuatore e il ritardo **mangiano fase** alla
frequenza di crossover, il notch **abbassa il guadagno** alla frequenza di
bending. Il compromesso fra queste due cose e' l'intero Task 2.

La difficolta' di testare questi blocchi e' che non hanno un "risultato" ovvio da
confrontare: sono funzioni di trasferimento. La suite risolve il problema
verificando **proprieta' analitiche in forma chiusa**, cioe' quantita' che si
sanno derivare a mano dalla struttura della `tf` e che si possono quindi
pretendere con tolleranze severissime (1e-9), invece che golden value numerici.
Le proprieta' scelte sono quattro:

1. **Guadagno statico unitario** (TVC e notch): un blocco di compensazione non
   deve riscalare il guadagno d'anello in continua, altrimenti sposterebbe
   silenziosamente tutti i margini.
2. **Ordine** della catena: 2 (attuatore) + n (Pade), senza cancellazioni
   accidentali.
3. **Il Pade riproduce davvero un ritardo**: la fase in eccesso rispetto al solo
   attuatore vale esattamente `-omega*tau` nella banda di lavoro.
4. **Profondita' e fase del notch**: `|H(j*omega_x)| = zeta_N / zeta_D` e i
   **zeri nel semipiano giusto** -- il discrimine fra la variante a fase minima
   (gain stabilisation) e quella a fase non minima (phase stabilisation).

Piu' tre test negativi sull'`arguments` block.

---

## Setup di classe (righe 7-19)

Identico a `hm3PlantTest`: una `PathFixture` che aggiunge `HM3/` al path (righe
11-14, con `fileparts(fileparts(...))` che risale da `HM3/tests/` a `HM3/`) e un
caricamento unico dei parametri nominali (righe 16-18). I parametri usati qui
sono `wTVC = 70` rad/s, `zTVC = 0.7`, `tau = 0.020` s, `wBM = 18.9` rad/s.

---

## `testTvcUnityDcGain` (righe 22-26)

```matlab
Wtvc = build_tvc(testCase.p);
testCase.verifyEqual(dcgain(Wtvc), 1, 'AbsTol', 1e-9);
```

- Riga 25: il guadagno statico dell'intera catena TVC deve essere **esattamente 1**.

**Derivazione.** `build_tvc.m` (righe 15-21) costruisce

      W_act(s)   = wTVC^2 / (s^2 + 2*zTVC*wTVC*s + wTVC^2)
      W_delay(s) = pade(tau, n)
      W_tvc(s)   = W_act(s) * W_delay(s)

A s = 0: `W_act(0) = wTVC^2 / wTVC^2 = 1`. Il Pade approssima `exp(-tau*s)`, che
in s = 0 vale 1, e l'approssimante di Pade **preserva esattamente** il valore in
s = 0 (e' costruito facendo coincidere i primi 2n termini dello sviluppo di
Taylor, a partire dal termine di ordine zero). Quindi `W_delay(0) = 1` e il
prodotto vale 1. Verificato in MATLAB: `dcgain(Wtvc) = 1.000000000000`.

**Perche' e' l'invariante giusto.** Un attuatore e' un dispositivo che *segue* il
comando: se gli chiedo 1 grado di deflessione, a regime deve darmi 1 grado, non
0.9. Se `W_tvc(0) != 1`, l'attuatore introdurrebbe un **guadagno nascosto**
nell'anello: tutti i margini calcolati sul Task 2 sarebbero traslati di
`20*log10(W_tvc(0))` dB rispetto a quelli del Task 1, e il confronto rigido/full
del Task 2 misurerebbe in parte un artefatto di normalizzazione invece
dell'effetto fisico del ritardo. E' il tipo di bug che non fa esplodere niente e
falsa tutte le conclusioni.

Nota: il fatto che il guadagno statico sia unitario **non** dice nulla sul
guadagno alla frequenza di crossover -- li' l'attuatore attenua e sfasa, ed e'
esattamente quello che il test successivo va a misurare.

---

## `testTvcOrderIsActuatorPlusPade` (righe 28-32)

```matlab
testCase.verifyEqual(order(build_tvc(testCase.p)),    5);  % default n = 3
testCase.verifyEqual(order(build_tvc(testCase.p, 2)), 4);
```

- Righe 30-31: l'ordine della catena e' `2 + n`, con n = 3 di default (5 stati) e
  n = 2 se passato esplicitamente (4 stati). Verificato: 5 e 4.
- **Che cosa protegge davvero.** Sembra contabilita', ma cattura due errori reali:
  (a) una **cancellazione polo-zero accidentale** (se l'attuatore e il Pade
  avessero un fattore comune, l'ordine scenderebbe e la fase non sarebbe piu'
  quella attesa); (b) un cambio silenzioso dell'ordine di default del Pade. Il
  punto (b) e' importante perche' l'ordine del Pade determina **fino a che
  frequenza l'approssimazione e' fedele**: un Pade di ordine n riproduce bene
  `exp(-j*omega*tau)` per `omega*tau` piccolo, e degrada oltre. Con tau = 0.020 s,
  la scala di frequenza del ritardo e' `1/tau = 50 rad/s`. Il modo di bending sta
  a omega_BM = 18.9 rad/s, cioe' `omega*tau = 0.378`: non e' banalmente piccolo.
  Il commento di `build_tvc.m` (righe 5-6) dice esplicitamente che un ordine piu'
  alto da' una fase migliore vicino a omega_BM. Se qualcuno abbassasse il default
  a 1, il notch verrebbe progettato contro una fase sbagliata: il test lo blocca.
- La riga 31 e' anche l'unico posto della suite che esercita il **secondo
  argomento opzionale** di `build_tvc` con un valore valido (la riga 83 lo passa
  solo per farlo rifiutare dal validatore).

> **Possibile domanda d'esame** -- Perche' Pade di ordine 3 e non un ritardo esatto?
> *Risposta:* Perche' tutta l'analisi di HM3 e' basata su strumenti razionali:
> `allmargin`, `pole`, `isstable`, `minreal`, `connect`. Un ritardo puro
> `exp(-tau*s)` non e' razionale, non ha una realizzazione in spazio di stato
> finita, e la Control System Toolbox lo tratterebbe come `InternalDelay`,
> rendendo impossibile `isstable(T)` sull'anello chiuso e complicando `allmargin`.
> Il Pade lo rende una `tf` a coefficienti reali, al prezzo di un errore che
> cresce con la frequenza. L'ordine 3 e' un compromesso: e' esatto in fase entro
> 1e-6 rad a 2 rad/s (dove leggo il margine di fase rigido) e ancora ragionevole
> a 18.9 rad/s (dove valuto il bending). Nota che il Pade sposta anche il **modulo**
> in modo lieve, mentre il ritardo vero e' all-pass esatto: e' un'approssimazione,
> non una riscrittura.

---

## `testTvcDelayPhaseAtLowFrequency` (righe 34-43)

```matlab
wTest = 2.0;                                   % rad/s, low frequency
h  = freqresp(build_tvc(pp), wTest);
hA = freqresp(tf(pp.wTVC^2, ...
      [1, 2*pp.zTVC*pp.wTVC, pp.wTVC^2]), wTest);
phaseDelay = angle(h) - angle(hA);
testCase.verifyEqual(phaseDelay, -wTest*pp.tau, 'AbsTol', 1e-6);
```

E' il test piu' istruttivo del file: **verifica che il ritardo sia un ritardo**.

- Riga 38: la frequenza di prova e' 2.0 rad/s. **Non e' arbitraria**: e' la banda
  in cui vive il crossover di controllo del Task 1 (2.45 rad/s secondo il README,
  2.455 rad/s verificato). Il test controlla la fedelta' dell'approssimazione
  **proprio dove viene letto il margine di fase**, che e' l'unico posto dove
  conta davvero.
- Righe 39-40: si calcola la risposta in frequenza della catena completa (`h`) e
  quella del **solo attuatore** (`hA`), ricostruito a mano dalla stessa formula
  usata in `build_tvc.m` riga 16.
- Riga 41: la differenza di fase **isola il contributo del Pade**:

      arg(W_tvc) = arg(W_act) + arg(W_delay)
      ->  arg(W_delay) = arg(W_tvc) - arg(W_act)

  E' l'unico modo di testare il fattore di ritardo senza riesporlo dall'interno
  della funzione: si sfrutta il fatto che la fase di un prodotto e' la somma
  delle fasi. Trucco pulito e riusabile.
- Riga 42: il valore atteso e' **la definizione stessa di ritardo di trasporto**:

      exp(-j*omega*tau)   ->   modulo 1,   fase = -omega*tau

  A omega = 2 rad/s e tau = 0.020 s: `-omega*tau = -0.04 rad = -2.29 deg`.
  Verificato in MATLAB: `phaseDelay = -0.04000000 rad`, cioe' l'errore del Pade
  di ordine 3 a questa frequenza e' sotto 1e-8 rad. La `AbsTol` di 1e-6 e' quindi
  onestamente raggiungibile e non e' una tolleranza gonfiata.
- **Perche' proprio questo e' l'invariante giusto.** Il ritardo e' pericoloso non
  perche' attenua (non attenua: `|exp(-j*omega*tau)| = 1`), ma perche' **sfasa in
  modo proporzionale a omega**. In un progetto in frequenza il ritardo e' un
  consumatore puro di margine di fase, e il suo costo cresce linearmente con la
  banda: e' la ragione per cui non si puo' semplicemente "alzare il guadagno"
  per gestire il bending. Il test verifica proprio la costante di proporzionalita'
  di questo consumo. Se il Pade fosse costruito con `tau` in millisecondi invece
  che in secondi, l'anello sembrerebbe avere 1000 volte piu' margine di ritardo e
  il test fallirebbe immediatamente.
- Contesto quantitativo utile all'orale: a omega = 2.45 rad/s (crossover rigido)
  il ritardo toglie `2.45 * 0.020 = 0.049 rad = 2.8 deg`; l'attuatore, misurato,
  ne toglie altri ~2.3 deg a 2 rad/s. Insieme sono ~5 deg su un margine di 30 deg:
  gestibili. E' il **notch** che, aggiungendo il suo ritardo di fase, fa crollare
  il margine da 30 a 14.6 deg (README, tabella Task 2) e costringe al re-tuning
  del PD.
- Piccola nota di onesta': il test ricostruisce `W_act` **riscrivendo la stessa
  formula** del sorgente. La parte "attuatore" del confronto e' quindi
  tautologica; il contenuto informativo sta tutto nella differenza di fase, che
  isola il Pade. Va bene cosi', ma va detto.

> **Possibile domanda d'esame** -- Perche' verifichi la fase del ritardo a 2 rad/s
> e non a omega_BM = 18.9 rad/s, dove il ritardo e' molto piu' grande?
> *Risposta:* Perche' a 2 rad/s l'approssimante di Pade e' praticamente esatto
> (errore < 1e-8 rad), quindi posso scrivere un'asserzione stringente con un
> valore atteso in forma chiusa (`-omega*tau`) e nessun golden value. A 18.9 rad/s
> il Pade di ordine 3 devia visibilmente dal ritardo ideale (`omega*tau = 0.378`
> non e' piccolo), e l'unico valore atteso possibile sarebbe un numero pinnato,
> cioe' un golden value senza contenuto teorico. Ho scelto di testare la
> proprieta' dove e' esatta e dove serve: 2 rad/s e' la banda del crossover
> rigido, che e' esattamente dove leggo il margine di fase.

---

## `testNotchUnityGainFarFromCentre` (righe 45-50)

```matlab
Hx = build_notch_filter(wx, 0.002, 0.7, +1);
testCase.verifyEqual(dcgain(Hx), 1, 'AbsTol', 1e-9);
testCase.verifyEqual(abs(freqresp(Hx, 1e4)), 1, 'AbsTol', 1e-3);
```

- Righe 48-49: il notch ha guadagno unitario **a entrambi gli estremi** dello
  spettro: in continua e a 1e4 rad/s (tre ordini di grandezza sopra omega_BM).

**Derivazione.** Il filtro (righe 20-23 di `build_notch_filter.m`) e'

      H(s) = (s^2 + sgn*2*zeta_N*omega_x*s + omega_x^2)
             --------------------------------------------
             (s^2 +     2*zeta_D*omega_x*s + omega_x^2)

- A s -> 0: numeratore e denominatore tendono entrambi a `omega_x^2` -> `H(0) = 1`
  **esattamente**, per qualunque zeta_N e zeta_D. Da qui la `AbsTol` di 1e-9.
- A s -> infinito: numeratore e denominatore sono entrambi dominati da `s^2` ->
  `|H| -> 1`, ma solo **asintoticamente**. Il residuo si ricava in forma chiusa dal
  modulo quadro:

        |H(j*omega)|^2 = (omega^2 - omega_x^2)^2 + (2*zeta_N*omega_x*omega)^2
                         ----------------------------------------------------
                         (omega^2 - omega_x^2)^2 + (2*zeta_D*omega_x*omega)^2

  Per omega >> omega_x lo scostamento e' **quadratico**, non lineare:
  `1 - |H| ~ (1/2)*(2*zeta_D*omega_x/omega)^2`. A omega = 1e4, con zeta_D = 0.7 e
  omega_x = 18.9: `0.5*(2*0.7*18.9/1e4)^2 = 3.5e-6`, e infatti
  `abs(freqresp(Hx, 1e4)) = 0.9999965`. La `AbsTol` di 1e-3 e' quindi
  **abbondantemente lasca** (basterebbe 1e-5): non e' calibrata sull'asintoto, e'
  solo un margine prudenziale. E' l'unica tolleranza del file che non sia stretta
  al limite numerico, e conviene saperlo dire prima che lo chieda il professore.

**Perche' e' l'invariante giusto.** E' la proprieta' che rende il notch **usabile**:
e' una modifica **locale in frequenza**. Deve scavare un buco attorno a omega_BM
e **non toccare niente altrove** -- in particolare non deve alterare il guadagno
d'anello attorno al crossover rigido (2-3 rad/s), dove vive il margine di fase,
ne' il guadagno in bassa frequenza, dove vive il margine di guadagno aerodinamico.
Se il notch avesse un guadagno statico diverso da 1, tutti i margini di bassa
frequenza del Task 2 sarebbero traslati e il confronto con il Task 1 non
significherebbe piu' nulla.

Attenzione pero': "guadagno unitario **in modulo** lontano dal centro" non vuol
dire "invisibile lontano dal centro". Il notch continua a introdurre **fase**
anche lontano dal buco -- ed e' proprio quella fase residua a far crollare il
margine da 30 a 14.6 deg nel Task 2. Il test verifica il modulo, non la fase, e
va detto esplicitamente: **il test non cattura il costo reale del notch**, che e'
in fase, non in guadagno. Quel costo emerge solo a valle, in `hm3LoopTest` e nei
`main_task2.m`.

---

## `testNotchDepthIsZetaRatio` (righe 52-59)

```matlab
wx = testCase.p.wBM;  zN = 0.002;  zD = 0.7;
Hn  = build_notch_filter(wx, zN, zD, +1);
Hll = build_notch_filter(wx, zN, zD, -1);
testCase.verifyEqual(abs(freqresp(Hn,  wx)), zN/zD, 'AbsTol', 1e-9);
testCase.verifyEqual(abs(freqresp(Hll, wx)), zN/zD, 'AbsTol', 1e-9);
```

E' il test che pinna **il buco**: dove sta, quanto e' profondo, e che non dipende
dal segno.

**Derivazione (da saper fare alla lavagna).** Valutiamo `H(s)` in `s = j*omega_x`,
cioe' esattamente al centro:

    numeratore:    (j*w_x)^2 + sgn*2*zN*w_x*(j*w_x) + w_x^2
                 = -w_x^2 + sgn*2*zN*w_x^2*j + w_x^2
                 = j * sgn * 2 * zN * w_x^2

    denominatore:  (j*w_x)^2 + 2*zD*w_x*(j*w_x) + w_x^2
                 = j * 2 * zD * w_x^2

I due termini reali (`-omega_x^2` e `+omega_x^2`) **si cancellano esattamente** in
entrambi. Resta un rapporto puramente immaginario:

    H(j*w_x) = (j*sgn*2*zN*w_x^2) / (j*2*zD*w_x^2) = sgn * zN/zD

    ->   |H(j*w_x)| = zN / zD          (indipendente da sgn e da w_x)

Ecco perche' i due `verifyEqual` alle righe 57-58 confrontano **entrambe le
varianti di segno** con lo stesso valore: il segno cambia la **fase** (0 oppure
180 deg al centro), non il **modulo**. E' un invariante forte, in forma chiusa,
che giustifica una `AbsTol` di 1e-9.

**La formula della profondita'** e' quindi:

    profondita' [dB] = 20*log10(zeta_N / zeta_D)

Con i valori del progetto ritenuto (zeta_N = 0.002, zeta_D = 0.7):

    zN/zD = 2.857e-3   ->   20*log10(2.857e-3) = -50.9 dB

Verificato in MATLAB: -50.88 dB, che e' il "-51 dB depth" citato nel README. **Il
notch e' un rapporto di smorzamenti, non un guadagno da tarare**: la profondita'
si progetta scegliendo zeta_N molto piccolo (zero quasi sull'asse immaginario,
che quasi cancella il polo di bending) e zeta_D dell'ordine dell'unita' (polo ben
smorzato, che non aggiunge risonanza propria). Questa e' la ragione fisica per cui
si chiama "notch": lo zero poco smorzato scava, il polo smorzato riempie e
richiude.

**Cosa il test NON dice.** Verifica il valore di `|H|` **in** omega_x, non che
omega_x sia il **punto di minimo** di `|H|`. Per zeta_N < zeta_D il minimo cade
(essenzialmente) in omega_x, quindi in pratica coincide, ma formalmente
l'asserzione e' piu' debole di "il buco e' alla frequenza giusta". Inoltre il test
passa `wx = p.wBM` ma **la funzione non sa nulla di omega_BM**: e' il chiamante a
doverla centrare. Che il notch sia effettivamente centrato sul bending **nel loop
reale** e' verificato altrove, in `hm3LoopTest.testDeepNotchStabilisesFullModel`,
dove si controlla `|L(j*omega_BM)| < -10 dB`.

**Il prezzo del notch profondo (onesta').** Il README lo dice e va ripetuto: la
larghezza del buco e' proporzionale a `zeta_D * omega_x`, ma la sua **parte
utile** (dove l'attenuazione e' davvero forte) e' strettissima con zeta_N = 0.002.
Il risultato e' che il progetto tollera -10 % di detuning di omega_BM ma diventa
instabile a +5 %. Un notch a -51 dB e' un notch che **pretende di conoscere
omega_BM**. Nessun test in questo file cattura quella fragilita': e' verificata
solo dai sweep di `main_task2.m` e dal Monte Carlo.

> **Possibile domanda d'esame** -- Come scegli la profondita' di un notch?
> *Risposta:* La profondita' e' fissata analiticamente dal rapporto degli
> smorzamenti, `20*log10(zeta_N/zeta_D)` dB, e va scelta in modo che il picco di
> risonanza del bending finisca sotto il requisito di gain margin flessionale. Nel
> mio caso il picco nudo e' +29 dB e il requisito tipico e' -12 dB (cioe' 12 dB di
> margine sotto lo 0 dB): serve almeno ~41 dB di attenuazione. Con zeta_N = 0.002
> e zeta_D = 0.7 ottengo -51 dB, che porta il lobo a |L(omega_BM)| = -21.9 dB con
> i guadagni del Task 1 e a -18 dB dopo il re-tuning: sopra i 12 dB richiesti. Il
> prezzo e' la sensibilita' al detuning di omega_BM e la fase persa attorno al
> crossover.

---

## `testDefaultVariantIsNonMinimumPhase` (righe 61-65)

```matlab
% sgn = -1 (assignment Eq. 4 as printed): both zeros in the RHP
Hx = build_notch_filter(testCase.p.wBM, 0.2, 0.5);
testCase.verifyTrue(all(real(zero(Hx)) > 0));
```

## `testMinimumPhaseVariantHasLhpZeros` (righe 67-70)

```matlab
Hx = build_notch_filter(testCase.p.wBM, 0.2, 0.5, +1);
testCase.verifyTrue(all(real(zero(Hx)) < 0));
```

Questi due test vanno letti insieme: sono la **distinzione concettuale piu'
importante di tutto il Task 2**.

**Derivazione.** Gli zeri sono le radici di `s^2 + sgn*2*zeta_N*omega_x*s + omega_x^2`:

    s = omega_x * ( -sgn*zeta_N  +/-  j*sqrt(1 - zeta_N^2) )

    ->  Re(s) = -sgn * zeta_N * omega_x

Quindi:
- `sgn = -1` (**default**, l'Eq. 4 **come stampata nella traccia**):
  `Re(s) = +zeta_N*omega_x > 0` -> **entrambi gli zeri nel semipiano DESTRO** ->
  sistema a **fase non minima (NMP)**. Verificato con zeta_N = 0.2, omega_x = 18.9:
  zeri in `+3.78 +/- 18.52j`, e infatti `+0.2*18.9 = 3.78`.
- `sgn = +1`: `Re(s) = -zeta_N*omega_x < 0` -> **zeri nel semipiano SINISTRO** ->
  **fase minima**. Verificato con la parametrizzazione del progetto ritenuto
  (zeta_N = 0.002): zeri in `-0.0378 +/- 18.9j`, e infatti `-0.002*18.9 = -0.0378`.

**Perche' questa distinzione e' il cuore del Task 2.** Un sistema a fase minima ha
la fase **univocamente determinata dal modulo** (relazione di Bode
guadagno-fase): se il modulo scende, la fase scende di conseguenza, ma non c'e'
"lag extra". Un sistema a fase **non** minima aggiunge fase negativa **a parita' di
modulo**: gli zeri RHP si comportano, in fase, come poli. Questo si traduce nelle
due strategie classiche di stabilizzazione dei modi flessionali:

- **Gain stabilisation** (quella ritenuta): si usa il notch **a fase minima**
  (`sgn = +1`) per **abbassare il modulo** del lobo di bending sotto 0 dB con
  ampio margine. Non importa piu' che fase abbia l'anello a omega_BM: se
  `|L| << 1`, la curva di Nichols non puo' circondare il punto critico li'.
  Richiede pero' di **conoscere omega_BM** con precisione (il buco e' stretto).
- **Phase stabilisation**: si lascia il lobo sopra 0 dB e si **modella la fase**
  perche' il lobo passi dalla parte giusta del punto critico. E' quello che fa
  la sezione lead-lag **a fase non minima** dell'Eq. 4 come stampata. E' piu'
  robusta al detuning in frequenza (non c'e' un buco stretto da centrare) ma
  molto piu' delicata: se la fase e' sbagliata di poco, il modo diventa instabile
  invece che stabilizzarsi.

La tabella del Task 2 nel README mostra il trade quantitativo: il lead-lag da
solo (Eq. 4 as printed, migliore di 75 combinazioni provate) porta il bending a
+23 dB e lascia un margine di fase rigido di 11.4 deg -- "barely" stabile. Il
notch profondo lo porta a -21.9 dB. Vince il notch.

**Onesta' -- che cosa i due test coprono davvero.** Coprono la **struttura** della
`tf` (dove finiscono gli zeri in funzione del segno), **non** la parametrizzazione
del filtro effettivamente montato nel progetto. Usano infatti `zeta_N = 0.2,
zeta_D = 0.5`, che sono i valori **di linea guida** citati nel commento di
`build_notch_filter.m` (righe 6-7: "guideline 0.1-0.3" e "0.4-0.6"), non i
`0.002 / 0.7` del notch ritenuto. La scelta e' comunque sensata: con zeta_N = 0.002
gli zeri stanno a `-0.0378` di parte reale, cioe' quasi sull'asse immaginario, e
un test di segno stretto (`> 0` / `< 0`) sarebbe numericamente al limite. Nota
inoltre il caso degenere: con **zeta_N = 0** gli zeri sarebbero **esattamente
sull'asse immaginario** e `real(zero(Hx)) < 0` fallirebbe -- eppure
`build_notch_filter` accetta `zeta_N = 0` (`mustBeNonnegative`, riga 15). E' un
buco nella validazione che nessun test copre; nella pratica non si presenta perche'
il progetto usa 0.002.

> **Possibile domanda d'esame** -- L'Eq. 4 della traccia, cosi' com'e' scritta, e'
> a fase non minima. Perche' l'hai riscritta a fase minima?
> *Risposta:* Non l'ho riscritta: `build_notch_filter` implementa **entrambe** le
> varianti tramite il parametro `numSign`, e il **default** e' proprio -1, cioe'
> l'Eq. 4 come stampata (zeri RHP). Ho poi **confrontato** le due nel trade dei
> filtri del Task 2. La sezione lead-lag NMP e' una tecnica di *phase
> stabilisation*: lascia il lobo di bending sopra 0 dB e cerca di ruotarne la fase
> perche' passi dalla parte sicura del punto critico. Nel mio caso, alla migliore
> delle 75 parametrizzazioni provate, lasciava il bending a +23 dB e collassava il
> margine di fase rigido a 11.4 deg: stabile ma inaccettabile. Il notch a fase
> minima (`numSign = +1`) fa *gain stabilisation*, porta il lobo a -21.9 dB, ed e'
> il progetto che ho ritenuto. Il prezzo, che dichiaro nel report, e' la
> sensibilita' al detuning di omega_BM.

---

## Test negativi sull'`arguments` block (righe 72-85)

```matlab
testCase.verifyError(@() build_notch_filter(-1, 0.2, 0.5), ...
    'MATLAB:validators:mustBePositive');
testCase.verifyError(@() build_notch_filter(18.9, 0.2, 0.5, 0), ...
    'MATLAB:validators:mustBeMember');
testCase.verifyError(@() build_tvc(testCase.p, 2.5), ...
    'MATLAB:validators:mustBeInteger');
```

- Righe 72-75: frequenza di centro negativa -> `mustBePositive` (riga 14 di
  `build_notch_filter.m`). Una `omega_x` negativa produrrebbe una `tf` con
  `omega_x^2` invariato ma il termine lineare col segno sbagliato: silenziosamente
  scambierebbe la variante NMP con quella a fase minima. Rifiutarla e' giusto.
- Righe 77-80: `numSign = 0` -> `mustBeMember(numSign, [-1, 1])` (riga 17).
  Con `numSign = 0` il numeratore diventerebbe `s^2 + omega_x^2`, cioe' zeri
  **esattamente sull'asse immaginario**: una cancellazione ideale del modo di
  bending. E' un progetto matematicamente elegante e **fisicamente indifendibile**
  (richiede conoscenza esatta di omega_BM e sarebbe instabile a qualunque
  detuning). Il validatore lo esclude per costruzione.
- Righe 82-85: ordine di Pade non intero -> `mustBeInteger` (riga 12 di
  `build_tvc.m`). `pade(tau, 2.5)` non ha senso.
- **Osservazione di metodo**: si verifica l'**identificatore** dell'errore
  (`'MATLAB:validators:mustBePositive'`), non il messaggio. E' la pratica corretta:
  i messaggi cambiano con la versione di MATLAB e con la lingua dell'installazione,
  gli ID no. Questi tre test sono anche l'unica documentazione eseguibile del
  **contratto d'uso** delle due funzioni.

---

## Cosa NON e' coperto (limiti noti)

- **La fase del notch attorno al crossover rigido.** E' il vero costo del notch
  (fa crollare il margine da 30 a 14.6 deg), e nessun test di questo file la
  misura. Solo il modulo e' verificato.
- **La larghezza del buco.** Nessun test lega `zeta_D` alla banda di attenuazione,
  quindi nessun test cattura la fragilita' al detuning (-10 % OK, +5 % instabile).
- **`zeta_N = 0`** e' accettato dall'`arguments` block e produrrebbe zeri sull'asse
  immaginario; nessun test lo documenta.
- **Fedelta' del Pade a omega_BM**, dove `omega*tau = 0.378` e l'errore non e'
  trascurabile: testata solo a 2 rad/s.

---

## Possibili domande d'esame

**D: Da dove viene la formula della profondita' del notch, e quanto vale nel tuo
progetto?**
R: Valutando `H(s)` in `s = j*omega_x` i termini reali `-omega_x^2` e `+omega_x^2`
si cancellano sia al numeratore sia al denominatore, e resta il rapporto di due
immaginari puri: `H(j*omega_x) = sgn * zeta_N/zeta_D`. Quindi il modulo al centro
e' esattamente `zeta_N/zeta_D`, indipendente da omega_x e dal segno, e la
profondita' in dB e' `20*log10(zeta_N/zeta_D)`. Con zeta_N = 0.002 e zeta_D = 0.7
il notch ritenuto vale 2.857e-3, cioe' **-50.9 dB**. Il test
`testNotchDepthIsZetaRatio` verifica questa identita' su entrambe le varianti di
segno, con tolleranza 1e-9, perche' e' esatta in forma chiusa.

**D: Che differenza c'e' fra gain stabilisation e phase stabilisation di un modo
flessionale, e quale hai usato?**
R: Nella gain stabilisation si abbassa il modulo del guadagno d'anello alla
frequenza del modo ben sotto 0 dB (tipicamente -12 dB o meglio), cosi' che il lobo
non possa circondare il punto critico qualunque sia la sua fase; serve un notch a
**fase minima** e la conoscenza precisa di omega_BM. Nella phase stabilisation si
lascia il lobo sopra 0 dB e si modella la **fase** perche' passi dalla parte
sicura; e' quello che fa la sezione lead-lag a **fase non minima** dell'Eq. 4 come
stampata nella traccia. Io ho usato la gain stabilisation: il notch a fase minima
porta il bending da +29 dB a -21.9 dB. Il lead-lag da solo, alla migliore delle 75
parametrizzazioni provate, lasciava +23 dB e un margine di fase rigido di 11.4 deg.

**D: Come fai a distinguere, con un test, un filtro a fase minima da uno a fase
non minima?**
R: Guardando i **segni delle parti reali degli zeri**. Per il mio biquad,
`Re(zero) = -numSign * zeta_N * omega_x`: con `numSign = -1` (default, Eq. 4 come
stampata) gli zeri stanno nel semipiano destro e il filtro e' NMP; con
`numSign = +1` stanno nel semipiano sinistro ed e' a fase minima. I due test
`testDefaultVariantIsNonMinimumPhase` e `testMinimumPhaseVariantHasLhpZeros`
verificano esattamente questo con `all(real(zero(Hx)) > 0)` e `< 0`. Fisicamente
conta perche' uno zero RHP aggiunge ritardo di fase a parita' di modulo: si
comporta, in fase, come un polo.

**D: Perche' il ritardo di 20 ms e' un problema, se non attenua nulla?**
R: Proprio perche' non attenua: e' all-pass, `|exp(-j*omega*tau)| = 1`. Il suo
effetto e' interamente in fase, `-omega*tau`, e cresce **linearmente con la
frequenza**. Quindi consuma margine di fase in proporzione alla banda che voglio
usare, e sul lanciatore la banda non e' negoziabile verso il basso (il polo
instabile a +1.84 rad/s impone un crossover ben sopra 1.84). Al mio crossover di
2.45 rad/s il ritardo costa 2.8 deg, che sono gestibili; a 18.9 rad/s costerebbe
21.7 deg, ed e' una delle ragioni per cui il bending si stabilizza in guadagno e
non in fase. Il test `testTvcDelayPhaseAtLowFrequency` verifica che il Pade
riproduca esattamente questa fase, isolandola dalla fase dell'attuatore per
differenza.

**D: Perche' il guadagno statico dell'attuatore e del notch deve essere
esattamente 1?**
R: Perche' altrimenti introdurrebbero un guadagno d'anello nascosto. Tutto il
progetto e' letto sui **margini**, che sono distanze in dB e in gradi rispetto al
punto critico: se un blocco "di compensazione" moltiplicasse `L` per una costante
diversa da 1, tutti i margini di bassa frequenza (in particolare il gain margin
aerodinamico, letto a 0.59 rad/s) si sposterebbero di `20*log10(k)` dB, e il
confronto fra il modello rigido del Task 1 e il modello completo del Task 2
misurerebbe in parte un artefatto di normalizzazione invece dell'effetto fisico
dell'attuatore e del bending. E' un errore che non fa esplodere nulla e falsa
tutte le conclusioni: da qui i due test con tolleranza 1e-9.
