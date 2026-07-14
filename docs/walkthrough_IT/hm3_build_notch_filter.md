# HM3/build_notch_filter.m

## Ruolo del file nel progetto

`build_notch_filter.m` costruisce la **sezione di filtro per la stabilizzazione
del modo di bending**: una biquadratica con numeratore e denominatore centrati
sulla stessa frequenza `wx` ma con smorzamenti diversi. E' l'Eq. 4 della
traccia, e in `assemble_loop` viene messa in cascata con il TVC nel ramo diretto
(`main_task2.m` riga 140: `Wfull = Wtvc * Hn`).

Serve perche' il modo flessionale del lanciatore (omega_BM = 18.9 rad/s,
zeta_BM = 0.005 -- praticamente non smorzato) **entra nelle misure** attraverso
l'INS (`build_plant_full.m`, righe 33-36: `theta_m = theta + sigma_ins*eta`) e
viene **eccitato dal TVC** (riga 27: la colonna `delta` forza `etadot` con
`-phi_tvc*Tc`). Chiudere il PD su una misura contaminata da un modo a
smorzamento 0.005 significa avere nel loop aperto un picco di risonanza di
**+29 dB** (il codice calcola `|L(omega_BM)| = 29.01 dB`, che stampa come
`29.0 dB`) che destabilizza l'anello:
`main_task2.m` Step B lo mostra con `isstable(Tb) = 0`. L'attuatore non aiuta
(a 18.9 rad/s attenua -0.01 dB, vedi `hm3_build_tvc.md`), quindi serve un filtro
dedicato.

Il file e' **una sola funzione parametrica** che copre due usi opposti, decisi
dal flag `numSign`:

- `numSign = +1` -> **notch a fase minima**: zeri complessi nel semipiano
  sinistro, buco profondo a `wx`, guadagno unitario lontano. Serve a
  **gain-stabilizzare** il modo (abbassarlo sotto 0 dB).
- `numSign = -1` (default, "Eq. 4 come stampata") -> **zeri nel semipiano
  destro**, filtro a **fase non minima**: stesso modulo del precedente, ma fase
  completamente diversa. Serve a **fase-stabilizzare** il modo (ruotarlo sulla
  Nichols invece di attenuarlo).

`main_task2.m` usa entrambi nel trade study dello Step C (righe 51-113) e alla
fine **ritiene il notch profondo a fase minima** (righe 74-78: `zN = 0.002`,
`zD = 0.7`, `sgn = +1`).

---

## `build_notch_filter` (righe 1-25)

```matlab
function Hx = build_notch_filter(wx, zN, zD, numSign)
num = [1, numSign*2*zN*wx, wx^2];
den = [1,         2*zD*wx, wx^2];
Hx  = tf(num, den);
```

- **Riga 1**: firma `Hx = build_notch_filter(wx, zN, zD, numSign)`. Quattro
  scalari -> una `tf` SISO. Chiamanti: `main_task2.m` (righe 58, 71, 78, 81-83,
  91, 106), `main_task3.m`, `main_montecarlo.m`, `init_simulink_hm3.m`,
  `tests/hm3FilterTest.m`.

- **Righe 13-18**: blocco `arguments`. Da notare i vincoli **effettivi**:
  `zN` e' `mustBeNonnegative` (quindi **zN = 0 e' ammesso**), `zD` e'
  `mustBePositive`, `numSign` e' vincolato a `{-1, +1}`. **Non c'e' nessun
  limite superiore** su `zN`/`zD`: i "guideline 0.1-0.3 / 0.4-0.6" della
  docstring (righe 6-7) sono un commento, non un vincolo. E infatti il design
  ritenuto (`zN = 0.002`) sta **due ordini di grandezza sotto** la guideline.
  Questa e' una scelta consapevole, non un bug, ma va detta all'orale: il filtro
  ritenuto **non e' il filtro dell'Eq. 4 con i valori suggeriti**, e' un notch
  profondo custom.

### Righe 20-21 -- la forma della biquadratica

    Hx(s) = (s^2 + sgn*2*zeta_N*wx*s + wx^2)
            -----------------------------------
            (s^2 +     2*zeta_D*wx*s + wx^2)

**Perche' questa forma e non un'altra.** Numeratore e denominatore hanno lo
**stesso termine costante** `wx^2` e lo **stesso coefficiente di s^2** (= 1).
Le conseguenze sono immediate e sono esattamente quello che si vuole da un
filtro di notch:

- **Guadagno unitario in continua**: `Hx(0) = wx^2 / wx^2 = 1`. Il filtro non
  tocca il guadagno di anello a bassa frequenza, quindi **non erode l'aero gain
  margin** (che vive a ~0.6-0.8 rad/s).
- **Guadagno unitario all'infinito**: per s -> infinito, `Hx -> s^2/s^2 = 1`.
  Il filtro non ha rolloff: **non e' un passa-basso**, non aiuta contro i modi
  flessionali superiori. E' un buco chirurgico e nient'altro.
- **Relativo di grado 0** (proprio ma non strettamente proprio): non aggiunge
  ritardo asintotico.

L'unica differenza fra num e den e' il **coefficiente del termine in s**. Quello
e' l'unico grado di liberta' che il filtro usa.

**La profondita' del buco.** Valutando in s = j*wx, i termini `s^2` e `wx^2` si
cancellano esattamente (`-wx^2 + wx^2 = 0`) e resta solo il termine in s:

    Hx(j*wx) = (sgn*2*zeta_N*wx * j*wx) / (2*zeta_D*wx * j*wx)
             = sgn * zeta_N / zeta_D

    => |Hx(j*wx)| = zeta_N / zeta_D      (profondita' del notch)

**Questa e' l'equazione di progetto.** Con `zN = 0.002` e `zD = 0.7`:

    |Hx(j*wBM)| = 0.002/0.7 = 0.002857  ->  20*log10(0.002857) = -50.88 dB

Verificato numericamente: la funzione restituisce esattamente -50.88 dB a
18.9 rad/s. Con `zN << zD` si ottiene il **buco profondo**, con `zN = zD` il
filtro degenera in `Hx = 1` (nessun effetto), con `zN > zD` diventa un
**amplificatore** risonante. Il caso limite `zN = 0` (ammesso dal validatore!)
mette gli zeri **esattamente sull'asse immaginario** in `+/- j*wx` e da'
profondita' **infinita** (`|Hx(j*wx)| = 0`), un notch ideale -- irrealizzabile in
pratica e numericamente fragile.

**Verifica del bilancio di guadagno sul bending** (il conto che chiude il
progetto di Task 2):

    |L(wBM)| senza filtro  = +29.01 dB
    profondita' del notch  = -50.88 dB
    ------------------------------------
    |L(wBM)| con il notch  =  -21.87 dB     (classify_margins: LwBM_dB = -21.87,
                                             stampato -21.9 nella tabella)

L'attenuazione e' **esattamente additiva in dB**, perche' il filtro e' in serie e
il suo modulo a `wx` non dipende dal resto del loop.

> **Nota di onesta' sui commenti del sorgente.** L'header di `main_task2.m`
> (righe 7-8) dichiara "the +39 dB bending resonance", con il numero scritto
> alla **riga 7**. Il valore vero, che lo script stampa a runtime alla riga 33
> (formato `%.1f`), e' **+29.0 dB** (29.01 dB il valore interno). Il README
> riporta correttamente +29 dB, ed e' consistente con l'aritmetica qui sopra
> (29.0 - 50.9 = -21.9 dB, che e' esattamente il numero nella tabella del
> README). Il **commento a riga 7 di `main_task2.m` e' quindi stale**.

### Righe 20 -- il flag `numSign` e i due modi di stabilizzare

Il numeratore e':

    numSign = +1:  s^2 + 2*zeta_N*wx*s + wx^2
                   zeri in  -zeta_N*wx +/- j*wx*sqrt(1 - zeta_N^2)   -> LHP
                   => FASE MINIMA

    numSign = -1:  s^2 - 2*zeta_N*wx*s + wx^2
                   zeri in  +zeta_N*wx +/- j*wx*sqrt(1 - zeta_N^2)   -> RHP
                   => FASE NON MINIMA

**Il modulo e' identico nei due casi.** Infatti

    num(j*w) con sgn=-1  =  coniugato di [ num(j*w) con sgn=+1 ]

perche' cambia solo il segno della parte immaginaria. Quindi |num| e' lo stesso,
e **le due varianti hanno lo stesso diagramma di Bode del modulo**. Cambia solo
la fase:

    arg(num_NMP) = -arg(num_MP)

    => fase(Hx_NMP) - fase(Hx_MP) = -2*arg(num_MP)

Sopra `wx`, `arg(num_MP)` tende a +180 gradi, quindi la variante a fase non
minima **perde ~360 gradi di fase attraversando la risonanza**. E' una rotazione
enorme e completamente controllata: e' lo strumento con cui si **ruota il lobo di
bending sulla carta di Nichols** attorno al punto critico -- che nella
convenzione degli appunti del corso sta a **(-180 deg, 0 dB)**, da 1 + L = 0
(D'Antuono lo etichetta +180: stesso punto, fase mod 360, ed e' la
rietichettatura che il codice usava fino a poco fa).

**Gain stabilization vs phase stabilization** -- la distinzione centrale
dell'homework:

| | **Gain stabilization** | **Phase stabilization** |
|---|---|---|
| Obiettivo | portare il lobo **sotto 0 dB** | lasciarlo sopra 0 dB ma con la **fase giusta** |
| Requisito tipico | `\|L(wBM)\|` <= -12 dB (spesso -20) | il lobo passa dalla parte "sicura" del punto critico |
| Strumento | notch profondo, `numSign = +1` | sezione NMP / lead-lag, `numSign = -1` |
| Cosa devo sapere bene | **omega_BM** (la frequenza) | **fase** e **segno della forma modale** al sensore |
| Fragilita' | detuning di omega_BM | errore di fase/segno del mode shape |
| Nel codice | `Hn` (righe 74-78 di main_task2) -- **ritenuto** | `Hll` (riga 71) -- **scartato** |

Il progetto di Task 2 sceglie la **gain stabilization**: il notch profondo porta
il bending a -21.9 dB (poi -18 dB dopo il ri-tuning del PD), cioe' un **bending
gain margin di 18 dB**, sopra i 12 dB tipicamente richiesti. La ragione e'
pragmatica: la gain stabilization richiede di conoscere solo *dove* sta il modo,
non *con che fase* arriva al sensore -- e la fase al sensore dipende dalla forma
modale, da `sigma_ins`, da `phi_ins`, cioe' da parametri strutturali molto meno
affidabili della frequenza.

### Il prezzo del notch profondo (due prezzi, entrambi pagati nel codice)

**Prezzo 1 -- il detuning.** Il buco e' profondo *ma strettissimo*. Vicino a `wx`
si puo' linearizzare: posto w = wx + dw con dw << wx,

    wx^2 - w^2 ~= -2*wx*dw
    |num(j*w)| ~= 2*wx*sqrt( dw^2 + (zeta_N*wx)^2 )

Il fondo del buco (dw = 0) vale `2*wx*zeta_N*wx`. Il modulo **raddoppia (+6 dB)**
gia' quando `dw = sqrt(3)*zeta_N*wx`. Con `zN = 0.002` e `wx = 18.9`:

    zeta_N * wx = 0.0378 rad/s      -> il "fondo" del notch e' largo ~0.08 rad/s

Cioe' **0.4% di omega_BM**. Un errore dell'1% sulla conoscenza di omega_BM
(0.19 rad/s) e' gia' 5 volte la mezza-larghezza del null. Verifica numerica
della profondita' effettiva al variare del disaccordo:

| omega_BM vero | `\|Hx\|` [dB] | `\|L(omega_BM)\|` risultante |
|---------------|-------------|----------------------------|
| 0.90 * wx | -16.53 | 29.01 - 16.53 = **+12.5 dB** |
| 0.95 * wx | -22.71 | **+6.3 dB** |
| **1.00 * wx** | **-50.88** | **-21.9 dB** |
| 1.05 * wx | -23.15 | **+5.9 dB** |
| 1.10 * wx | -17.38 | **+11.6 dB** |

Un disaccordo del 5% **butta via 28 dB di attenuazione** e riporta il bending
sopra 0 dB. E' esattamente lo Step D di `main_task2.m` (righe 163-179), che
verifica la stabilita' con i filtri fissi e omega_BM vero perturbato di +/-10%,
e il README riporta il risultato: il notch profondo regge -10% ma **va instabile
a +5%**. Da notare l'onesta' del dettaglio: la tabella qui sopra e' quasi
**simmetrica**, quindi l'asimmetria -10% ok / +5% instabile **non viene dalla
profondita'** ma dalla **fase** del loop alla nuova frequenza del modo -- cioe'
da quale lato del punto critico (-180 deg, 0 dB) passa il lobo. Il codice
non chiarisce questo punto; la spiegazione va cercata sulla Nichols.

**Prezzo 2 -- il ritardo di fase sul crossover rigido.** Sotto `wx` il notch e'
in **ritardo di fase**. Il conto: per w < wx entrambe le parti reali sono
positive, quindi

    fase(Hx) = atan(2*zeta_N*wx*w / (wx^2 - w^2)) - atan(2*zeta_D*wx*w / (wx^2 - w^2))

e con `zN << zD` il primo termine e' trascurabile: **resta solo il ritardo del
denominatore**, che e' governato da `zeta_D`. Al crossover rigido di Task 1
(w = 2.4547 rad/s, `mR.rigidPM_w`; il README arrotonda a 2.45) il notch
(`zN=0.002`, `zD=0.7`, `wx=18.9`) costa **-10.45 gradi**. Confronto con gli
altri contributi alla stessa frequenza:

| contributo | fase a 2.4547 rad/s |
|---|---|
| servo TVC (70 rad/s) | -2.81 gradi |
| ritardo 20 ms (Pade-3) | -2.81 gradi |
| **notch profondo** | **-10.45 gradi** |
| **totale** | **-16.07 gradi** |

**Il notch e' il contributore dominante di ritardo di fase**, piu' dell'attuatore
e del ritardo messi insieme. Ed e' precisamente il conto che spiega lo Step C
decision di `main_task2.m` (righe 131-140): il phase margin rigido di Task 1
(30.0 gradi -- il target su cui `design_controller` converge) meno ~16 gradi da'
~14 gradi -- e infatti il codice misura **14.6 gradi** con i guadagni di Task 1
sul loop completo (righe 144-149, `mB.rigidPM_deg`; e' lo stesso numero della
tabella del README). Da qui l'obbligo di **ri-tarare il PD sul loop pieno**
(righe 152-154), alzando `Kd_th` da 0.44 a 0.69 per recuperare la fase persa.

Nota: `zeta_D` e' il parametro che governa questo prezzo. `zD = 0.7` e' generoso
(largo) e costa i 10 gradi. Un `zD` piu' piccolo restringerebbe il denominatore e
ridurrebbe il ritardo lontano da wx, ma renderebbe il filtro piu' risonante e
piu' fragile. Il codice non esplora questo trade-off: `zD = 0.7` e' fissato a
mano alla riga 76 di `main_task2.m`.

### Riga 24 -- `Hx.Name = 'Notch_Hx'`

Il nome e' hard-coded a `'Notch_Hx'` **anche quando la funzione viene usata come
lead-lag a fase non minima** (`numSign = -1`). E' cosmetico (`Name` non entra in
`connect`, che lavora su `InputName`/`OutputName`), ma e' fuorviante leggendo un
`Hll` che si chiama "Notch".

> **Possibile domanda d'esame** -- Perche' hai scelto la gain stabilization e non
> la phase stabilization, visto che l'Eq. 4 della traccia e' stampata proprio
> nella forma a fase non minima?
> *Risposta:* Perche' `main_task2.m` (Step C, righe 51-70) fa il trade study e lo
> misura: sui 75 candidati lead-lag della guideline (`wx = wBM +/- 4`,
> `zN` in 0.1-0.3, `zD` in 0.4-0.6) **nessuno** stabilizza da solo il loop in modo
> soddisfacente -- il migliore lascia il bending a +23 dB e un phase margin rigido
> di 11 gradi. Il motivo strutturale e' che la guideline da' profondita'
> `zN/zD` fra 0.1/0.6 = -15.6 dB e 0.3/0.4 = -2.5 dB: **troppo poco per un picco
> di +29 dB**. La sezione NMP puo' solo *ruotare* il lobo, e per farlo passare
> dalla parte giusta del punto critico servirebbe conoscere con precisione la fase
> del modo al sensore (`sigma_ins`, `phi_ins`, la forma modale) -- informazione
> molto meno affidabile della sola frequenza. La gain stabilization compra
> robustezza con l'unico parametro di cui mi fido: omega_BM.

---

## Possibili domande d'esame

**D: Ricava la profondita' del notch e spiega perche' dipende solo dal rapporto
degli smorzamenti.**
R: In s = j*wx i termini `s^2` e `wx^2` si cancellano identicamente in numeratore
e denominatore (`-wx^2 + wx^2 = 0`), perche' i due polinomi hanno lo stesso
termine costante e lo stesso coefficiente di s^2. Resta solo il termine lineare:
`Hx(j*wx) = (2*zeta_N*wx*j*wx)/(2*zeta_D*wx*j*wx) = zeta_N/zeta_D`. Tutto il
resto si semplifica, quindi la profondita' e' **puramente il rapporto degli
smorzamenti**, indipendente da wx. Con 0.002/0.7 si ottengono -50.9 dB.
Corollario di progetto: la profondita' e' governata da `zN`, la larghezza (e il
ritardo di fase lontano dal buco) da `zD`. Sono due manopole quasi ortogonali.

**D: Perche' il filtro ha guadagno unitario sia in continua sia all'infinito, e
perche' e' importante?**
R: Perche' num e den condividono il coefficiente di `s^2` (= 1) e il termine
costante (`wx^2`). In DC: `wx^2/wx^2 = 1`; all'infinito: `s^2/s^2 = 1`. E'
importante per due motivi. Primo: **non tocca l'aero gain margin**, che vive a
~0.6-0.8 rad/s, dove il filtro e' trasparente in modulo (il loop e'
condizionalmente stabile e un'attenuazione a bassa frequenza sarebbe
destabilizzante, non conservativa). Secondo: **non e' un passa-basso**. Non da'
alcun rolloff, quindi non fa nulla contro modi flessionali di ordine superiore o
rumore ad alta frequenza -- per quelli servirebbe un filtro diverso.

**D: Gain stabilization e phase stabilization -- cosa cambia operativamente?**
R: Gain stabilization = porto `|L(omega_BM)|` ben sotto 0 dB (requisito tipico
-12 dB o meglio), cosi' il lobo di bending sulla Nichols non puo' avvolgere il
punto critico perche' e' semplicemente troppo in basso. Serve conoscere bene
**omega_BM**. Phase stabilization = lascio il lobo sopra 0 dB ma ne aggiusto la
fase perche' passi dalla parte sicura del punto critico (mantenendo il numero
corretto di avvolgimenti richiesto da Nyquist). Serve conoscere bene la **fase**
del modo al sensore, cioe' segno e ampiezza della forma modale -- un'informazione
strutturale, tipicamente incerta. Qui la scelta e' gain stabilization
(`numSign = +1`, `zN = 0.002`), e la sezione NMP dell'Eq. 4 (`numSign = -1`), che
e' lo strumento per la phase stabilization, viene testata e scartata.

**D: Il notch profondo costa 10 gradi di fase al crossover rigido. Come li
recuperi?**
R: Ri-tarando il PD **sul loop completo**, non su quello rigido. `main_task2.m`
riga 153 richiama `design_controller` passando `Wfull = Wtvc*Hn` come catena
attuatore: il tuner (`fminsearch` sui log dei guadagni) ottimizza `Kp_th`, `Kd_th`
contro gli stessi target (aero |GM| 6 dB, rigid PM 30 gradi) ma **vedendo** il
ritardo di fase di attuatore + delay + notch. Il risultato e' `Kd_th` che sale da
0.44 a 0.69: piu' azione derivativa = piu' anticipo di fase al crossover, che
ricompra i gradi persi. Il concetto e' quello di D'Antuono: le dinamiche del TVC
(e qui anche del filtro) **devono entrare nel piazzamento del PD**, non essere
trattate come una perturbazione a posteriori.

**D: Il tuo notch e' fuori dalla guideline della traccia (zN 0.1-0.3, zD 0.4-0.6).
Come lo giustifichi?**
R: Con l'aritmetica. La guideline da' profondita' `zN/zD` fra -2.5 dB e -15.6 dB.
Il picco di bending nel loop e' **+29 dB**. Nessuna combinazione della guideline
puo' portarlo sotto 0 dB, e infatti lo sweep sui 75 candidati
(`main_task2.m`, righe 54-67) lo conferma. La guideline dell'Eq. 4 descrive una
sezione **lead-lag / phase-shaper** (con il numeratore NMP), non un notch di
gain-stabilization: e' lo strumento per l'altra strategia. Per gain-stabilizzare
serve `zN/zD` dell'ordine di 1e-3, che e' cio' che il codice usa. Il validatore
lo permette (nessun limite superiore su zN/zD, `zN` solo `mustBeNonnegative`), ma
va dichiarato: il filtro ritenuto **non e' l'Eq. 4 con i valori suggeriti**.

**D: Cosa succede se il notch e' disaccordato del 5%?**
R: L'attenuazione crolla da -50.9 dB a circa -23 dB, e il bending nel loop torna a
+5.9 dB, cioe' **sopra 0 dB**: il modo non e' piu' gain-stabilizzato. Il motivo e'
che la mezza-larghezza del fondo del buco vale ~`zeta_N*wx` = 0.038 rad/s, cioe'
lo **0.2% di omega_BM**: e' un buco largo 0.08 rad/s su una frequenza di 18.9
rad/s. Lo Step D di `main_task2.m` verifica questo e trova che il notch profondo
va instabile a +5% di errore su omega_BM. La contromisura di corso (il "triplet"
di notch a 0.9/1.0/1.1*wBM, righe 81-83) copre la banda ma costa ~30 gradi di
ritardo al crossover rigido e destabilizza il loop rigido (PM = -7.3 gradi nella
tabella del README) -- il codice lo prova e lo scarta.
