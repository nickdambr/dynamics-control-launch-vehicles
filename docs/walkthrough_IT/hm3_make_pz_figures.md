# HM3/make_pz_figures.m

## Ruolo del file nel progetto

E' uno **script di sola visualizzazione**: non progetta niente, non cambia
nessun guadagno, non produce numeri che finiscono nelle tabelle del report.
Prende gli stessi oggetti costruiti dai `main_task*.m` (plant, attuatore TVC,
notch, PD) e li guarda **nel piano s** invece che sulla carta di Nichols.
Serve a rispondere a una domanda che la Nichols non risponde in modo diretto:
*dove stanno fisicamente i poli, e dove li porta la retroazione?*

La scelta di progetto di HM3 e' interamente in frequenza (Nichols, margini
classificati per banda). La mappa poli-zeri e' il **controllo incrociato**: se
la Nichols dice "stabile con 6 dB e 30 gradi", la mappa deve mostrare tutti i
poli di anello chiuso nel semipiano sinistro. Le quattro figure prodotte sono,
nell'ordine dei banner `%%`:

1. tre **luoghi delle radici** (`intro_rootlocus`, `task2_rootlocus`,
   `task2_rootlocus_full`) -- la stabilita' condizionata;
2. `task1_polemap` -- poli ad anello aperto vs. anello chiuso del loop rigido;
3. `task2_notch_poles` -- la coppia del bending con e senza notch;
4. `task3_pole_migration` -- migrazione dei poli sui vertici del box +/-30%.

Dipende da: `load_hw3_params`, `build_plant_rigid`, `build_plant_full`,
`build_tvc`, `build_notch_filter`, `design_controller`, `assemble_loop`.
Scrive PNG in `HM3/figures/` (che poi il report LaTeX pesca via
`\graphicspath`). Non e' chiamato da nessun altro script: si lancia a mano.

> **Attenzione, e' il punto piu' importante di questa pagina.** Le figure 2, 3
> e 4 **non usano gli stessi guadagni**. La figura 2 usa il PD del Task 1
> (`Kp = 1.7845`, `Kd = 0.4433`, riga 40). Le figure 3 e 4 usano il PD
> **ri-sintonizzato** del Task 2 (`Kp = 1.7318`, `Kd = 0.6867`, righe 55-56).
> Verificato eseguendo il codice. Vedi sotto per le conseguenze.

---

## Intestazione, stile grafico, cartella di uscita (righe 1-35)

- Righe 1-16: commento di testa che dichiara le quattro figure. E' onesto e
  corretto. Le righe 6-7 promettono *"the +sqrt(A6) open-loop pole pulled into
  the LHP by the rigid PD"*, e la figura 2 **lo mostra davvero** -- ma solo
  dopo l'allargamento del `ylim` di riga 73. Fino a poco fa la finestra era
  troppo stretta e la promessa restava sulla carta: la storia e' raccontata
  per esteso nella sezione righe 67-79, e vale la pena saperla.
- Riga 18: `clear; close all; clc;` -- script, non funzione: gira in base
  workspace.
- Riga 19: `warning('off','Control:analysis:MarginUnstable')`. Necessario
  perche' l'anello e' aperto-instabile: ogni chiamata a `margin`/`allmargin`
  dentro `design_controller` emette un warning. Qui viene solo silenziato,
  non aggirato.
- Righe 21-23: `here = fileparts(mfilename('fullpath'))` -- lo script e'
  eseguibile da qualunque cartella corrente, le figure vanno sempre in
  `HM3/figures/`.
- Righe 26-35: palette. `cUnst` (rosso) = instabile / anello aperto /
  senza notch; `cStab` (blu) = stabile / anello chiuso / con notch; `cRHP` =
  ombreggiatura del semipiano destro; `cCorn` = 5 colori per Nominale + V1..V4.

---

## Punto di progetto condiviso: quali guadagni entrano nelle figure (righe 37-41)

```matlab
p0     = load_hw3_params();
Grigid = build_plant_rigid(p0);
K      = design_controller(Grigid, [], 'verbose', false);
zN = 0.002; zD = 0.7; sgn = +1;   % deep min-phase notch
```

- Riga 38: `load_hw3_params()` **senza argomenti** -> caso nominale
  (`mu_alpha_scale = mu_c_scale = 1`), istante di riferimento `t_ref = 72 s`
  (max-qbar). I coefficienti arrivano dal dataset LPV, non dai letterali della
  Tabella 1: eseguendo il codice si ottiene
  `A6 = 3.381828`, `K1 = 4.564724`, `a1 = -0.015423`, `a3 = 20.6090`,
  `a4 = -27.2710`, `V = 937.709 m/s`, `wBM = 18.9 rad/s`, `zBM = 0.005`,
  `qbar = 81129 Pa`, `Alt = 15143 m`.
- Riga 40: `K` e' il PD del **Task 1**, ritarato da `design_controller` sul
  loop rigido con attuatore ideale. Eseguito: `Kp_th = 1.7845`,
  `Kd_th = 0.4433`, `Kp_z = Kd_z = -1e-3`.
- Riga 41: parametri del notch. `zN = 0.002` (numeratore, quasi senza
  smorzamento -> notch profondo), `zD = 0.7` (denominatore, ben smorzato),
  `sgn = +1` -> **notch a fase minima** (zeri nel semipiano sinistro). Il
  default di `build_notch_filter` sarebbe `-1` (Eq. 4 come stampata nella
  traccia, zeri nel semipiano destro): qui viene **deliberatamente scartato**.
  Profondita' analitica del notch alla frequenza centrale:

      |Hx(j*wBM)| = zN/zD = 0.002/0.7 = 0.002857  ->  -50.9 dB

  (verificato numericamente: -50.88 dB).

---

## Anatomia del piano s: da dove viene ogni polo e ogni zero (righe 43-58)

Queste righe costruiscono tutto cio' che poi si vede nelle figure: TVC, plant
flessibile, notch, PD ri-tarato, e i due anelli aperti `Lrig` e `Lful`.

```matlab
Wtvc   = build_tvc(p0,3);
Gfull  = build_plant_full(p0,'ins');
notch  = build_notch_filter(p0.wBM, zN, zD, sgn);
K2 = design_controller(Gfull, Wtvc*notch, 'w_flex',0.6*p0.wBM, ...
                       'w_flex_hi',1.5*p0.wBM, 'w_bending',p0.wBM, ...
                       'verbose',false);
[Lrig,~] = assemble_loop(Grigid, K,  []);
[Lful,~] = assemble_loop(Gfull,  K2, Wtvc*notch);
```

- Riga 51: `build_tvc(p0,3)` -- attuatore del secondo ordine
  (`wTVC = 70 rad/s`, `zTVC = 0.7`) **in serie** con l'approssimante di Pade di
  ordine 3 del ritardo puro `tau = 20 ms`.
- Riga 52: `build_plant_full(p0,'ins')` -- plant a 6 stati con la
  contaminazione INS della misura (Eq. 2): `theta_m = theta + sigma_ins*eta`,
  `z_m = z - phi_ins*eta`. E' proprio questo canale che destabilizza il modo.
- Righe 55-56: `K2` e' il PD **ri-sintonizzato sull'anello completo**
  (plant flessibile + TVC + ritardo + notch). Eseguito: `Kp_th = 1.7318`,
  `Kd_th = 0.6867`. Le opzioni `w_flex = 0.6*wBM = 11.34`,
  `w_flex_hi = 1.5*wBM = 28.35`, `w_bending = wBM` servono a
  `classify_margins` per capire quali crossing sono "rigidi" e quali
  "flessibili".
- Righe 57-58: i due anelli aperti (convenzione `1 + L`, rottura sul segnale
  `delta`).

### Tabella dei poli/zeri (tutti verificati eseguendo il codice)

| oggetto | poli | zeri | origine fisica |
|---|---|---|---|
| `Grigid` (4 stati) | `0`, `+0.0291`, `+1.8165`, `-1.8610` | (verso theta) `-0.0317` | rigido + deriva laterale |
| bending in `Gfull` | `-0.0945 +/- 18.8998i` | -- | `wBM=18.9`, `zBM=0.005` |
| attuatore TVC | `-49.00 +/- 49.99i` | -- | `wTVC=70`, `zTVC=0.7` |
| Pade(tau=0.02, ord. 3) | `-232.2`, `-183.9 +/- 175.4i` | `+232.2`, `+183.9 +/- 175.4i` | ritardo puro -> **fase non minima** |
| notch | `-13.23 +/- 13.50i` | `-0.0378 +/- 18.9i` | `zD=0.7` / `zN=0.002` |
| coppia INS/bending (in `Lful`) | -- | `+15.26`, `-15.08` | `sigma_ins * phi_tvc * Tc` |

### 1. Il polo instabile: da dove esce `+sqrt(A6)`

La dinamica rotazionale disaccoppiata (traccia, Eq. 1, riga 4 della matrice
di stato) e'

    theta_ddot = A_6*theta + (A_6/V)*z_dot + K_1*delta - A_6*alpha_w

Se si ignora l'accoppiamento con la deriva (`z_dot ~ 0`) e il vento, resta

    theta_ddot = A_6*theta + K_1*delta   ->   G(s) = K_1/(s^2 - A_6)

cioe' **due poli reali** in `s = +/- sqrt(A_6) = +/- 1.83897 rad/s`. Il segno
`+` di `A_6` e' il cuore del problema: `A_6 = N_alpha*l_alpha/I_yy` con
`l_alpha > 0` = distanza fra centro di pressione e centro di massa, con il
**CP davanti al CM**. Un'incidenza positiva genera un momento che *aumenta*
l'incidenza -> divergenza aerodinamica. Verifica dei numeri di Tabella 1:
`1.07e6 * 10.39 / 3.28e6 = 3.389 ~ A_6 = 3.3818`.

Il tempo di raddoppio dell'errore d'assetto e'

    t_2 = ln(2)/sqrt(A_6) = 0.693/1.839 = 0.377 s

Ecco perche' il controllo e' obbligatorio: senza retroazione il lanciatore
raddoppia l'errore d'assetto ogni 0.38 s.

### 2. Ma nella mappa il polo NON e' a +1.839: e' a +1.8165 (e ce n'e' un secondo)

E' un dettaglio che vale un punto all'orale. Il plant rigido **completo**
include la deriva laterale, e la sua equazione caratteristica non e'
`s^2 - A_6`. Sviluppando il determinante della `A` di `build_plant_rigid`
(lo stato `z` non retroagisce su nessuno, quindi da' un polo in `0`; restano
`z_dot, theta, theta_dot`):

    s*[ (s - a_1)*(s^2 - A_6) - (a_1*V + a_4)*A_6/V ] = 0

Con i valori nominali `a_1*V + a_4 = -14.46 - 27.27 = -41.73`, il termine di
accoppiamento vale `+0.1505` e il polinomio cubico diventa

    s^3 + 0.0154*s^2 - 3.3818*s + 0.0984 = 0

le cui radici sono `+1.8165`, `-1.8610`, `+0.0291`. Quindi:

- il polo aerodinamico **si sposta** da `+1.839` a `+1.8165` (la deriva lo
  addolcisce di circa l'1%);
- compare un **secondo polo instabile lentissimo** in `+0.0291 rad/s`
  (costante di tempo 34 s), il modo di **deriva**: il veicolo, se non
  controllato, oltre a capottare scivola lateralmente sottovento;
- nel canale `delta -> theta` questo polo e' quasi cancellato da uno zero in
  `-0.0317` (dipolo lento). *Quasi*, non del tutto: il polo e' a destra, lo
  zero a sinistra. E' un dipolo di quelli che non si possono ignorare, perche'
  contiene il segno.

Quindi il plant rigido ha **due poli RHP**, non uno. Il commento del codice e
il README parlano solo di `+sqrt(A6)`: e' la lettura disaccoppiata, comoda per
il progetto, ma la mappa disegnata dallo script mostra la verita' completa.

> **Possibile domanda d'esame** -- *Nella mappa vedo due croci a destra
> dell'asse immaginario. La seconda e' un errore numerico?*
> *Risposta:* No. La prima e' il polo aerodinamico (`+1.82`, ~`+sqrt(A6)`
> perturbato dall'accoppiamento con la deriva). La seconda, in `+0.029`, e' il
> modo di deriva laterale: nasce dai termini `a_1` e `a_4` della riga di
> `z_ddot` combinati con `A_6/V`. E' instabile ma lentissimo (`tau = 34 s`), e
> nel solo canale `delta -> theta` e' quasi cancellato da uno zero in `-0.032`.
> La retroazione debole sulla deriva (`Kp_z = Kd_z = -1e-3`) serve proprio a
> gestirlo.

### 3. La coppia del bending: perche' e' praticamente sull'asse

Il primo modo flessionale entra come oscillatore del secondo ordine
`s^2 + 2*zBM*wBM*s + wBM^2`, con `wBM = 18.9 rad/s` e `zBM = 0.005`. I poli:

    s = -zBM*wBM +/- i*wBM*sqrt(1 - zBM^2) = -0.0945 +/- 18.8998i

La parte reale e' **-0.09**, contro una parte immaginaria di **18.9**: nella
scala della figura sono di fatto *sull'asse immaginario*. Fisicamente:
lo smorzamento strutturale di un serbatoio in alluminio riempito di propellente
e' dello 0.5%, quindi il modo, se eccitato, suona per
`tau = 1/(zBM*wBM) = 10.6 s`. Una struttura di questo tipo non si "smorza da
sola" nei tempi del volo: **o la si evita, o la si spegne attivamente**.

### 4. Gli zeri a fase non minima: due famiglie diverse (attenzione!)

**(a) Pade.** Il ritardo `exp(-tau*s)` non e' razionale; `pade(0.02, 3)` lo
approssima con una funzione **tutta-passa**: gli zeri sono l'immagine
speculare dei poli rispetto all'asse immaginario. Con l'approssimante [3/3]
il polinomio (in `x = tau*s`) e' `x^3 - 12x^2 + 60x - 120`, con radici
`x = 4.644` e `x = 3.678 +/- 3.509i`; dividendo per `tau = 0.02`:

    zeri  : +232.2 ,  +183.9 +/- 175.4i   (RHP, fase non minima)
    poli  : -232.2 ,  -183.9 +/- 175.4i

E' esattamente cio' che dichiara il commento di riga 65 ("+232 and
+184+/-175i"), verificato. Il modulo e' unitario a ogni frequenza (all-pass):
il ritardo **non toglie guadagno, toglie solo fase** -- circa
`-tau*w = -0.02*3.2 = -0.064 rad = -3.7 gradi` al crossover rigido, ma
`-0.02*18.9 = -21.6 gradi` alla frequenza del bending. E' il motivo per cui la
fase precipita ad alta frequenza e per cui il *delay margin* (213 ms nel Task 1,
165 ms nel Task 2) e' una metrica riportata separatamente.

**(b) La coppia +/-15 dalla contaminazione INS.** Questa il commento del codice
la *non* menziona correttamente (riga 63 parla solo di "the -2.5/-15 real
zeros"), ma eseguendo `zero(Lful)` compare uno **zero reale nel semipiano
destro a +15.26** (e il suo compagno a `-15.08`). Da dove viene? La misura e'
`theta_m = theta + sigma_ins*eta`, e il TVC forza il bending con
`-phi_tvc*Tc`. Trascurando la deriva:

    theta/delta ~ K_1/(s^2 - A_6)
    eta/delta   = -c/(s^2 + 2*zBM*wBM*s + wBM^2),   c = phi_tvc*Tc = 65.6

quindi il numeratore di `theta_m/delta` e'

    (K_1 - sigma_ins*c)*s^2 + K_1*2*zBM*wBM*s + (K_1*wBM^2 + sigma_ins*c*A_6)

Il coefficiente di `s^2` vale `4.5647 - 0.178*65.6 = 4.565 - 11.67 = -7.11`:
**e' negativo**. Il segno del termine dominante si ribalta, e il numeratore ha
due radici reali quasi simmetriche rispetto all'asse:

    s ~ +/- sqrt(1670/7.11) = +/- 15.3 rad/s

(numericamente `+15.26` e `-15.08`). Significato fisico: la deflessione
strutturale vista dall'INS e' **in controfase** rispetto alla risposta rigida,
e in modo cosi' forte da dominarla. E la conseguenza per il progetto e'
decisiva: gli zeri stanno **sotto** `wBM` (15.3 < 18.9). Quando gli zeri
precedono il polo del modo, il lobo di bending sulla Nichols entra con la fase
"sbagliata" e **non e' stabilizzabile in fase**: l'unica strada e'
**stabilizzarlo in guadagno**, cioe' schiacciarlo sotto 0 dB con il notch.
Questo e' il motivo strutturale, non estetico, per cui il Task 2 usa un notch
e non un lead-lag.

> **Possibile domanda d'esame** -- *Perche' non stabilizzi il bending in fase
> (phase stabilization) invece di usare un notch profondo?*
> *Risposta:* Perche' la contaminazione INS produce una coppia di zeri reali a
> `+/-15.3 rad/s`, cioe' **sotto** `wBM = 18.9` e uno dei due nel semipiano
> destro. Con zeri sotto il polo, il lobo di bending attraversa la Nichols
> dalla parte sbagliata e nessuna rotazione di fase realizzabile lo porta fuori
> dal punto critico. La stabilizzazione in fase e' possibile solo quando il
> sensore vede il modo con la fase "giusta" (zeri sopra il polo); qui non e' il
> caso, quindi si va per forza in gain stabilization.

---

## `(1)` I tre luoghi delle radici (righe 60-65)

- Righe 60-61: root locus del **loop rigido del Task 1** (attuatore ideale,
  guadagni `K`). Finestra `[-10, 2.5] x [-4, 4]`.
- Righe 62-63: root locus del **loop completo del Task 2**, zoom vicino
  all'origine, `[-16, 3] x [-3, 3]`.
- Righe 64-65: lo stesso loop, vista completa `[-250, 250] x [-200, 200]`, per
  far vedere i rami che corrono verso gli zeri RHP di Pade.

`rlocusplot(L)` traccia le radici di `1 + k*L(s) = 0` al variare di `k` da 0 a
infinito. **Il punto di progetto e' `k = 1`**: il luogo dice quanto sono
"lontani" da instabilita' i poli se il guadagno d'anello viene scalato.

`Lrig` in forma zeri-poli-guadagno (verificato):

    L_rig(s) = 2.003*(s + 4.032)*(s^2 + 0.0562*s + 0.03222)
               ---------------------------------------------
                  s*(s + 1.861)*(s - 1.816)*(s - 0.02909)

Note utili all'orale:
- il guadagno di alta frequenza `2.003` e' essenzialmente `Kd*K1 = 0.4433 *
  4.5647 = 2.024`: a frequenza alta il ramo derivativo domina il PD;
- lo zero in `-4.032` e' il rapporto `Kp/Kd = 1.7845/0.4433 = 4.03` -- e' lo
  zero del PD, non del plant. Il commento di riga 61 lo chiama giustamente
  "the -4.03 zero";
- l'anello ha **4 poli e 3 zeri** -> un solo asintoto (a 180 gradi, verso
  `-inf`);
- il polo in `0` c'e' perche' la deriva e' retroazionata (`Kp_z*z_m`): senza
  la retroazione di deriva quel polo non comparirebbe in `L`.

**Stabilita' condizionata.** Eseguendo `classify_margins(Lrig)` si ottiene
`aeroGM = -6.00 dB @ 0.593 rad/s` (margine di **riduzione** di guadagno) e
nessun margine di aumento (`rigidGM = NaN`, l'attuatore e' ideale). Tradotto:

    k_min = 10^(-6/20) = 0.501    (nessun k_max nel Task 1)

Il sistema e' stabile **solo se il guadagno d'anello supera il 50% del
nominale**. E' l'opposto dell'intuizione da corso base ("abbassa il guadagno e
sei al sicuro"): con un plant aperto-instabile serve *guadagno minimo*, perche'
sotto quella soglia il ramo che parte dal polo `+1.8165` non ha ancora
attraversato l'asse immaginario. E' proprio quello che si vede nel luogo:
per `k` piccolo i poli chiusi coincidono con quelli aperti (due a destra), e
solo crescendo `k` migrano a sinistra.

Nel loop completo (`Lful`) si aggiunge il **limite superiore**: il ritardo e
l'attuatore fanno tornare la fase a `-180` gradi ad alta frequenza, quindi
`k` troppo grande destabilizza. Da qui il "rigid GM" (7.56 dB nel Task 2) --
la banda di guadagno ammissibile e' *chiusa*, ed e' la firma della stabilita'
condizionata sul lanciatore. `Lful` ha 13 poli e 10 zeri -> 3 asintoti, e i
rami "lontani" corrono verso gli zeri RHP di Pade (`+232`, `+184 +/- 175i`):
per questo la figura di riga 64-65 esiste, altrimenti si vedrebbero solo rami
che escono dal quadro.

> **Possibile domanda d'esame** -- *Che differenza c'e' fra il margine di
> guadagno "aerodinamico" e quello "rigido"?*
> *Risposta:* Sono due attraversamenti diversi della fase `-180` gradi dello
> stesso anello. Quello aerodinamico e' a bassa frequenza (0.59 rad/s) ed e' un
> margine di **riduzione**: se il guadagno cala di 6 dB il ramo del polo
> instabile torna a destra. Quello rigido e' al crossover di controllo ed e' un
> margine di **aumento**: se il guadagno cresce, il ritardo TVC porta la fase a
> `-180` e si destabilizza dall'altra parte. Un singolo numero di `margin()`
> non ha senso qui: il codice li separa in `classify_margins.m`.

---

## `(2)` Task 1: mappa dei poli in anello chiuso (righe 67-79)

```matlab
[~,T1] = assemble_loop(Grigid, K, []);
pOL = pole(Grigid);  pCL = pole(T1);
...
xl = [-6 2];  yl = [-2.5 2.5];   % tall enough for the PD-placed pair at -0.95 +/- 1.90i
```

- Riga 68: `assemble_loop(Grigid, K, [])` -- `[]` = attuatore ideale (`Wact = 1`),
  cioe' il Task 1 puro.
- Riga 69: si estraggono i poli aperti (croci rosse) e chiusi (cerchi blu).
- Righe 74-76: `shade_rhp` ombreggia il semipiano destro, poi le due serie.
- Riga 73: **i limiti degli assi**. `yl = [-2.5, 2.5]` e' abbastanza alto da
  contenere la coppia dominante di anello chiuso a `-0.95 +/- 1.90i`, che e'
  **il soggetto stesso della figura**. Su questa riga c'era un difetto, corretto
  di recente: vedi il riquadro in fondo alla sezione.

I poli, calcolati eseguendo il codice:

| | poli |
|---|---|
| aperti (`Grigid`) | `0`, `+0.0291`, `+1.8165`, `-1.8610` (tutti reali) |
| chiusi (`T1`) | `-0.9533 +/- 1.9047i` e `-0.0559 +/- 0.2329i` |

La coppia dominante `-0.9533 +/- 1.9047i` ha `wn = 2.13 rad/s`,
`zeta = 0.448`; la coppia lenta `-0.0559 +/- 0.2329i` ha `wn = 0.239 rad/s`,
`zeta = 0.234`, `tau = 17.9 s` (e' il modo di deriva stabilizzato -- lo stesso
che giustifica l'orizzonte di 80 s scelto in `main_task1.m`).

**Il legame con il progetto analitico.** Sulla dinamica rotazionale
disaccoppiata, chiudere con `delta = -(Kp + Kd*s)*theta` da'

    s^2 + Kd*K_1*s + (Kp*K_1 - A_6) = 0

da cui, immediatamente,

    wn   = sqrt(Kp*K_1 - A_6)  = sqrt(1.7845*4.5647 - 3.3818) = 2.183 rad/s
    zeta = Kd*K_1/(2*wn)       = 0.4433*4.5647/(2*2.183)      = 0.464

che sono esattamente i valori (`w_c = 2.18`, `zeta = 0.46`) riportati nel
README. I poli **veri** dell'anello chiuso completo sono `2.13 / 0.448`: la
differenza (~2%) e' l'effetto della retroazione di deriva, che il calcolo
analitico ignora. Vale la pena saperlo: se all'orale si chiede "coincidono?",
la risposta corretta e' "quasi, e la discrepanza e' il prezzo della retroazione
di deriva".

Nota sui **semi** usati da `design_controller` (D'Antuono Eq. 3.6-3.7):
`Kp0 = 2*A_6/K_1 = 1.4817`, `Kd0 = sqrt(A_6)/K_1 = 0.4029`. Sostituendo nelle
formule sopra: `wn = sqrt(2*A_6 - A_6) = sqrt(A_6) = 1.839` e
`zeta = sqrt(A_6)/(2*sqrt(A_6)) = 0.5` esatti. Cioe' i semi analitici sono la
scelta "specchia il polo instabile a sinistra con `zeta = 0.5`". Poi
`fminsearch` li sposta a `1.7845 / 0.4433` per centrare i target 6 dB / 30
gradi *sull'anello completo*.

> **Difetto corretto, ma da saper raccontare.** Oggi la riga 73 impone
> `ylim = [-2.5, 2.5]` e la figura mostra **tutti e quattro** i poli di anello
> chiuso -- la coppia veloce piazzata dal PD a `-0.95 +/- 1.90i` *e* la coppia
> lenta di deriva a `-0.06 +/- 0.23i` -- insieme al polo instabile ad anello
> aperto a `+1.82` dentro la fascia RHP ombreggiata. E' esattamente cio' che il
> titolo promette.
>
> **Fino a poco fa non era cosi'.** La riga 73 imponeva `ylim = [-0.28, 0.28]`,
> cioe' una finestra alta poco piu' di un quarto di rad/s: dentro ci stava
> **solo** la coppia lenta a `+/-0.233i`. La coppia dominante, a `+/-1.905i`,
> cadeva **fuori quadro**. La figura intitolata *"closed-loop pole placement"*
> tagliava fuori proprio i poli che il PD piazza: si vedeva il problema (il polo
> a `+1.82`) ma non la soluzione (la sua immagine stabilizzata). Il difetto era
> tanto piu' insidioso perche' la figura *sembrava* corretta -- croci rosse,
> cerchi blu, tutto al suo posto -- e solo confrontando i numeri con
> `pole(T1)` si scopriva che meta' dei cerchi mancava.
>
> **Perche' `[-0.28, 0.28]` era la scelta sbagliata anche in linea di
> principio.** Una mappa di *pole placement* deve essere scalata sulla coppia
> **piazzata**, non su quella residua. La coppia veloce ha `wn = 2.13 rad/s`: il
> riquadro deve contenere almeno `+/-wn` in immaginario, altrimenti si sta
> inquadrando la dinamica che il progetto **non** controlla (il modo di deriva)
> e si nasconde quella che controlla. Il difetto e' il tipo di errore che passa
> inosservato in una figura che "gira" senza errori -- e per questo e' una bella
> domanda d'orale: *come fai a sapere che la tua figura mostra cio' che dici che
> mostra?* Risposta: la si valida contro i numeri, non contro l'occhio.

---

## `(3)` Task 2: la coppia del bending con e senza notch (righe 81-96)

```matlab
[~,Tno] = assemble_loop(Gfull, K2, Wtvc);         % senza notch
[~,Tnf] = assemble_loop(Gfull, K2, Wtvc*notch);   % con notch
```

- Riga 82: anello chiuso **senza** notch (solo TVC + ritardo).
- Riga 83: anello chiuso **con** il notch profondo.
- Righe 88: finestra `[-30, 25] x [-45, 45]`.
- Righe 90-91: `yline` a `wBM` per ancorare visivamente la frequenza del modo.

Risultati (verificati eseguendo il codice):

| | poli in zona bending | stabile? |
|---|---|---|
| senza notch | `+2.7582 +/- 16.3135i` -- **RHP** | no (`isstable = 0`) |
| con notch | `-0.1023 +/- 18.9191i` | si' (`isstable = 1`) |

Lettura fisica:
- **Senza notch**, la coppia flessibile (aperta a `-0.0945 +/- 18.90i`, cioe'
  a 0.09 dalla stabilita') viene **spinta a destra dall'anello**: il guadagno
  d'anello a `wBM` e' positivo (+29 dB con i guadagni del Task 1), la fase e'
  quella sbagliata per via degli zeri INS a `+/-15`, e la retroazione **pompa
  energia nel modo** invece di sottrarla. E' la classica instabilita' di
  interazione controllo-struttura: il TVC eccita la flessione, l'INS la legge
  come assetto, il controllore reagisce alla flessione e la rinforza.
- **Con notch**, la coppia torna a `-0.1023 +/- 18.9191i`, cioe'
  **praticamente dove stava ad anello aperto** (`-0.0945 +/- 18.90i`).
  E' *esattamente* il significato di *gain stabilization*: il notch non smorza
  il modo, gli toglie il guadagno d'anello (`|L(wBM)| = -18 dB`) cosi' che la
  retroazione **non lo veda**. Lo smorzamento di anello chiuso resta
  `zeta = 0.1023/18.92 = 0.0054`, contro lo `0.005` strutturale.
- Compare inoltre una coppia nuova, `-7.7342 +/- 10.8265i`: sono i **poli del
  notch stesso** (`-13.23 +/- 13.50i` ad anello aperto), trascinati dalla
  retroazione. Ben smorzati, innocui.

> **Limite onesto della figura.** La comparazione "senza notch" e' fatta con
> `K2`, cioe' i guadagni gia' **ri-sintonizzati sul loop col notch** (il
> commento di riga 82 lo ammette). Nel README/report, invece, la riga di
> tabella "no filter, |L(wBM)| = +29 dB" e' valutata con i guadagni del
> **Task 1**. Non e' quindi la stessa configurazione: la figura risponde alla
> domanda "il PD ri-tarato da solo basterebbe a domare il bending?" (no), non
> alla domanda "cosa succede al progetto del Task 1 se aggiungo il modo
> flessibile?". La conclusione qualitativa (instabile senza notch) non cambia,
> ma le due cifre non sono confrontabili.

> **Limite onesto della finestra.** Con `xlim = [-30, 25]`, i poli
> dell'attuatore TVC (`-49 +/- 50i`) e quelli di Pade (`-232`,
> `-184 +/- 175i`) sono **fuori quadro**. La figura mostra il cluster
> rigido/bending/notch, non l'insieme completo dei 13 poli. Il titolo dice
> "bending poles", quindi e' legittimo, ma non si puo' dire "la figura mostra
> che tutti i poli sono stabili".

> **Possibile domanda d'esame** -- *Se il notch riporta i poli del bending
> praticamente dove stavano ad anello aperto, a che serve? Il modo resta
> smorzato allo 0.5%.*
> *Risposta:* E' proprio l'obiettivo. Il modo flessibile e' una proprieta'
> della struttura: il controllo d'assetto non ha ne' l'autorita' ne' la banda
> per smorzarlo (servirebbe un attuatore molto piu' veloce e un sensore
> pulito). Cio' che si deve garantire e' che il loop **non lo destabilizzi**.
> Senza notch il loop lo porta a `+2.76 +/- 16.3i`, instabile; con il notch lo
> lascia a `-0.10 +/- 18.9i`, cioe' esattamente com'e' in natura. La riserva di
> sicurezza si legge in guadagno (`|L(wBM)| = -18 dB`, sopra i 12 dB
> tipicamente richiesti), non in smorzamento.

---

## `(4)` Task 3: migrazione dei poli sul box +/-30% (righe 98-117)

```matlab
cases = {'Nominal',1.00,1.00; 'V1',0.70,0.70; 'V2',0.70,1.30; ...
         'V3',1.30,0.70; 'V4',1.30,1.30};
...
p  = load_hw3_params('mu_alpha_scale',cases{i,2},'mu_c_scale',cases{i,3});
Gf = build_plant_full(p,'ins');
Wf = build_tvc(p,3) * build_notch_filter(p0.wBM, zN, zD, sgn);
[~,T] = assemble_loop(Gf, K2, Wf);
```

- Righe 99-100: i 4 vertici del box di incertezza su `mu_alpha = A_6`
  (piu' o meno instabile) e `mu_c = K_1` (piu' o meno autorita' di controllo),
  piu' il nominale.
- Riga 109: **solo** `A_6` e `K_1` vengono scalati (lo fa `load_hw3_params`
  alle righe 89-90).
- Riga 111: il notch e' costruito su `p0.wBM` -- "frozen @ nominal". Il
  commento e' corretto in spirito, ma **in Task 3 e' ridondante**: `wBM`,
  `wTVC`, `zTVC`, `tau` non sono fra i parametri scalati, quindi
  `build_tvc(p,3)` e `build_notch_filter(p.wBM,...)` darebbero comunque gli
  stessi identici sistemi. La distinzione conta davvero solo in
  `main_montecarlo.m`, dove `wBM` viene disperso.
- Riga 112: il controllore **non cambia mai** (`K2` fisso). E' la definizione
  stessa dello studio di robustezza: guadagni congelati, plant che si muove.

Massima parte reale dei poli di anello chiuso (verificata eseguendo il codice):

| caso | mu_alpha | mu_c | max Re(polo) | stabile |
|---|---|---|---|---|
| Nominale | 1.00 | 1.00 | `-0.0497` | si' |
| V1 | 0.70 | 0.70 | `-0.0502` | si' |
| V2 | 0.70 | 1.30 | `-0.0297` | si' |
| **V3** | **1.30** | **0.70** | `-0.1023` | si' |
| V4 | 1.30 | 1.30 | `-0.0494` | si' |

Tutti nel semipiano sinistro -> la conclusione del commento di riga 12
("all poles stay in the LHP -> robust") e' **corretta**.

> **Trappola concettuale, da preparare.** Il polo piu' vicino all'asse **non e'
> mai** quello aerodinamico ne' quello del bending: e' il **modo lento di
> deriva** (`Re ~ -0.03 ... -0.10`). Di conseguenza la mappa **non ordina i
> casi come li ordinano i margini**: V3, che il Task 3 identifica come il caso
> peggiore (aero GM = 0.91 dB, PM = 18 gradi), ha in realta' la parte reale
> **piu' negativa** di tutte (`-0.1023`), mentre V2 -- un caso comodo -- e'
> quello con il polo piu' vicino all'asse (`-0.0297`). Non e' una
> contraddizione: la parte reale misura la *velocita' di decadimento*, il
> margine di guadagno misura *quanto si puo' sbagliare il guadagno prima di
> perdere la stabilita'*. Sono robustezze diverse. La mappa poli-zeri **non e'
> in grado di vedere** che V3 e' a 0.91 dB dal disastro. Per quello serve la
> Nichols. Se all'orale si presenta questa figura come "prova di robustezza",
> bisogna sapere che e' una prova di *stabilita' nominale ai vertici*, non di
> *robustezza dei margini*.

---

## Export (righe 119-128)

- Righe 120-127: per ogni figura si forza il tema chiaro (`theme(f,'light')`,
  in `try/catch` perche' `theme` esiste solo da R2025a in poi) e si esporta a
  220 DPI con `exportgraphics`.
- Riga 126: il nome del PNG e' `get(f,'Name')` -- quindi i nomi dei file
  (`task1_polemap.png`, `task2_notch_poles.png`, `task3_pole_migration.png`)
  sono decisi dalla proprieta' `Name` delle `figure` alle righe 71, 86, 103.
  Le tre figure di root locus vengono invece esportate dentro `rl_fig`.

---

## Helper (righe 131-168)

- `shade_rhp` (righe 131-137): disegna una `patch` sul semipiano destro con
  `FaceAlpha = 0.08` e `HandleVisibility='off'` (non entra in legenda). C'e' la
  guardia `if xl(2) > 0`: se la finestra e' tutta a sinistra dell'asse non
  disegna nulla.
- `finish_ax` (righe 139-149): stile comune (font, griglia, assi, etichette in
  `s^-1`). Nota che l'unita' e' corretta: le parti reale e immaginaria di un
  polo sono frequenze, in rad/s = s^-1.
- `rl_fig` (righe 151-168): wrapper attorno a `rlocusplot` -- **non ridisegna
  nulla a mano**, usa il rendering nativo del Control System Toolbox e ne cambia
  solo assi, titolo e font via `setoptions`. Scelta deliberata (dichiarata nel
  commento): le convenzioni di croce/cerchio/direzione dei rami restano quelle
  standard MATLAB, cosi' la figura e' leggibile da chiunque conosca il tool.
  Le tre figure di root locus **non** ricevono l'ombreggiatura RHP.

---

## Possibili domande d'esame

**D: Cosa aggiunge la mappa poli-zeri rispetto alla carta di Nichols, se il
progetto e' fatto tutto in frequenza?**
R: Tre cose che la Nichols non dice esplicitamente. (1) *Dove* stanno i poli
instabili e quanto sono veloci: `+1.84 rad/s` significa raddoppio in 0.38 s, e
questo fissa la banda minima del controllo. (2) *Dove* finisce ogni polo dopo la
chiusura: la Nichols garantisce la stabilita' ma non dice se il modo di deriva
si e' fermato a `tau = 18 s` (informazione che serve per scegliere l'orizzonte
di simulazione). (3) La distinzione fra gain e phase stabilization del bending:
sulla mappa e' evidente che il notch **riporta** i poli flessibili dove stavano,
mentre non li smorza. Di contro, la mappa non vede i margini di guadagno: e' il
motivo per cui le due viste sono complementari e non alternative.

**D: Con quali guadagni sono state generate queste figure?**
R: Dipende dalla figura, ed e' importante saperlo. La figura 2 (`task1_polemap`)
e il primo root locus usano il PD del Task 1, `Kp = 1.7845`, `Kd = 0.4433`,
progettato su plant rigido con attuatore ideale. Le figure 3
(`task2_notch_poles`), 4 (`task3_pole_migration`) e i due root locus del loop
completo usano il PD **ri-sintonizzato** del Task 2 sul loop con TVC + ritardo +
notch, `Kp = 1.7318`, `Kd = 0.6867`. In tutti i casi i guadagni di deriva sono
fissi a `Kp_z = Kd_z = -1e-3`.

**D: Perche' il PD ri-tarato ha un `Kd` molto piu' alto (0.69 contro 0.44) ma
un `Kp` quasi uguale?**
R: Perche' cio' che si e' perso, aggiungendo attuatore + ritardo + notch, e'
**fase** al crossover, non guadagno. La catena TVC (2o ordine a 70 rad/s) + Pade
(20 ms) + notch (denominatore a `zD = 0.7`) introduce insieme una ventina di
gradi di ritardo di fase intorno a 3 rad/s, e il PM del Task 1 crolla da 30 a
14.6 gradi con i guadagni originali. Il termine derivativo e' l'unico che
*aggiunge* fase (`+atan(w*Kd/Kp)`): alzarlo sposta lo zero del PD da
`Kp/Kd = 4.03` a `2.52 rad/s`, cioe' piu' vicino al crossover, recuperando
anticipo di fase. `Kp` invece fissa il margine aerodinamico a bassa frequenza
(6 dB), che non e' cambiato: quindi resta dov'era.

**D: Nel root locus del loop completo alcuni rami vanno verso destra, verso
`+232` e `+184 +/- 175i`. Il sistema e' instabile ad alto guadagno?**
R: Si', ed e' fisicamente corretto. Quegli zeri nel semipiano destro sono
l'approssimante di Pade del ritardo di 20 ms: un ritardo puro e' un all-pass che
a frequenza crescente accumula fase senza limite, quindi qualunque anello con
ritardo diventa instabile se si alza abbastanza il guadagno. E' esattamente la
lettura in piano-s del *delay margin* (165 ms nel Task 2). Va detto pero' che
gli zeri sono **artefatti dell'approssimazione**: il ritardo vero non ha zeri, ha
solo fase. L'approssimante di Pade di ordine 3 e' fedele fino a circa
`3/tau = 150 rad/s`; i rami che corrono verso `+232` sono ben oltre quella banda
e non vanno interpretati quantitativamente.

**D: Il modo di bending, ad anello chiuso, ha smorzamento 0.0054. Non e'
pochissimo?**
R: Si', ed e' voluto. Il modo e' *gain-stabilized*, non *phase-stabilized*: il
notch lo mette a `-18 dB` di guadagno d'anello, quindi il controllore non lo
eccita e non lo smorza. Resta lo smorzamento strutturale (0.5%). La conseguenza
pratica: se qualcosa lo eccita comunque (separazione, transitorio del motore,
turbolenza a banda larga) la struttura suona per `1/(zeta*wBM) ~ 10 s`. Il
requisito che si soddisfa non e' "il modo e' smorzato" ma "il modo non e'
destabilizzato dal controllo, con almeno 12 dB di margine". Il prezzo, dichiarato
nel README, e' che il notch profondo richiede di conoscere `wBM` con precisione:
tollera -10% di detuning ma va instabile a +5%.

**D: Qual e' il limite piu' importante di queste figure?**
R: La figura del Task 3 dimostra la stabilita' ai quattro vertici, ma il polo
piu' vicino all'asse e' **sempre** il modo di deriva, quindi la figura **non
ordina** i casi per criticita': V3, il caso peggiore per margini (0.91 dB),
risulta il piu' "lontano" dall'asse. Non e' un errore di calcolo, e' un limite
intrinseco della lettura in piano-s: la robustezza in guadagno non si legge
dalla posizione dei poli nominali. Serve la Nichols.

C'era anche un **difetto vero**, ora corretto, che vale la pena citare
spontaneamente: `task1_polemap` aveva `ylim = [-0.28, 0.28]`, mentre la coppia
di anello chiuso piazzata dal PD sta a `-0.95 +/- 1.90i` -- **fuori dal
riquadro**. La figura che doveva mostrare il pole placement mostrava solo il
modo lento di deriva. Oggi la riga 73 usa `ylim = [-2.5, 2.5]` e tutti e quattro
i poli chiusi sono in quadro, insieme al polo aperto instabile a `+1.82`.
