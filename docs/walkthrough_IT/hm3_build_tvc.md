# HM3/build_tvc.m

## Ruolo del file nel progetto

`build_tvc.m` costruisce il blocco **attuatore** della catena di controllo: la
funzione di trasferimento che sta fra il comando di deflessione calcolato dal
PD (`u_pd`, o `delta_cmd`) e la deflessione fisica dell'ugello (`delta`) che
entra nel piano dell'`build_plant_rigid` / `build_plant_full`. Implementa
l'Eq. 3 della traccia: un **servo di secondo ordine** (omega_TVC = 70 rad/s,
zeta_TVC = 0.7) in serie con un **ritardo puro di trasporto** tau = 20 ms.

Nel flusso dell'homework il file non esiste in Task 1: `main_task1.m` chiama
`design_controller(G, [])` (riga 22) e legge il loop da `m.L` / `m.T` (righe
30-31); l'alias `[]` diventa `Wact = tf(1)` dentro `assemble_loop`
(`assemble_loop.m`, righe 18 e 21), cioe' **attuatore ideale**.
Il file entra in scena in `main_task2.m` (riga 30, `Wtvc = build_tvc(p, 3)`) e
poi in `main_task3.m` / `main_montecarlo.m`. E' esattamente il passaggio in cui
il progetto smette di essere accademico: l'attuatore e il ritardo aggiungono
**fase negativa senza aggiungere attenuazione utile**, e insieme al modo di
bending sono cio' che fa crollare il phase margin rigido da 30 gradi a 14.6
gradi (vedi `main_task2.m`, "Step C decision", righe 131-139 -- il commento nel
codice arrotonda a ~15), obbligando a ri-tarare il PD sul loop completo.

Il punto di sostanza -- ed e' il punto che va difeso all'orale -- e' che il
ritardo puro **non e' una funzione di trasferimento razionale**: exp(-tau*s) e'
trascendente. Tutta la catena a valle (`allmargin`, `pole`, `isstable`,
`minreal`, `fminsearch` nel tuner) lavora su modelli razionali. La scelta di
approssimare con **Pade** non e' cosmetica: e' cio' che rende il resto del
codice eseguibile, e ha una motivazione fisica precisa (l'approssimante di Pade
e' **all-pass**: modulo esattamente unitario, solo fase -- come il ritardo
vero).

Dipendenze: solo `p` (struct da `load_hw3_params`, campi `wTVC`, `zTVC`, `tau`)
e la Control System Toolbox (`tf`, `pade`). Nessuno stato, nessun file su disco.

---

## `build_tvc` (righe 1-23)

```matlab
function Wtvc = build_tvc(p, padeOrder)
s = tf('s');
Wact = p.wTVC^2 / (s^2 + 2*p.zTVC*p.wTVC*s + p.wTVC^2);
[nd, dd] = pade(p.tau, padeOrder);
Wdelay = tf(nd, dd);
Wtvc = Wact * Wdelay;
```

- **Riga 1**: firma `Wtvc = build_tvc(p, padeOrder)`. Restituisce una singola
  `tf` SISO da `u_pd` a `delta`. Il chiamante la passa come terzo argomento di
  `assemble_loop`, dove viene inserita **in cascata nel ramo diretto**, fra il
  guadagno statico del controllore e l'ingresso `delta` del plant
  (`assemble_loop.m`, righe 30-32). Chiamanti reali: `main_task2.m` (riga 30),
  `main_task3.m`, `main_montecarlo.m`, `init_simulink_hm3.m`.

- **Righe 10-13**: blocco `arguments`. `padeOrder` ha default **3** ed e'
  vincolato intero positivo. Nota: questa e' una funzione di *boundary* (viene
  chiamata una volta per run, non dentro un loop di `fsolve`/`fmincon`), quindi
  la validazione qui e' coerente con la convenzione della repo (le funzioni
  hot-loop non ce l'hanno per costo).

### Riga 16 -- il servo di secondo ordine

    Wact(s) = omega_TVC^2 / (s^2 + 2*zeta_TVC*omega_TVC*s + omega_TVC^2)

Forma canonica di un secondo ordine con **guadagno statico unitario**
(`Wact(0) = omega^2/omega^2 = 1`): a regime l'ugello raggiunge esattamente la
deflessione comandata, come deve essere per un servo posizionale. I valori
vengono da `load_hw3_params.m` (righe 55-56): omega_TVC = 70 rad/s, zeta_TVC =
0.7.

Perche' questi numeri contano: la banda dell'attuatore (70 rad/s) sta **quasi un
ordine di grandezza sopra** il crossover rigido di progetto (~2.5-3.2 rad/s) ma
solo un fattore ~3.7 sopra il modo di bending (omega_BM = 18.9 rad/s). Le
conseguenze sono due, e sono opposte a quello che l'intuizione suggerisce:

- **Al crossover rigido** (w = 2.445 rad/s, il crossover di Task 1) l'attuatore
  costa pochissima fase: sfasamento = -atan2(2*zeta*w*omega, omega^2 - w^2) =
  -2.80 gradi. Praticamente trasparente.
- **A omega_BM = 18.9 rad/s** costa -22.18 gradi di fase, ma il suo modulo e'
  |Wact(j*18.9)| = 0.9988, cioe' **-0.01 dB**. Zero attenuazione.

Questo secondo punto e' il piu' importante e va detto esplicitamente: **non si
puo' contare sull'attuatore per attenuare il bending**. Il rolloff del servo
comincia a 70 rad/s, il modo flessionale sta a 18.9 rad/s: l'attuatore lo lascia
passare integro (+29 dB di picco nel loop, verificato: `|L(omega_BM)| = 28.95
dB` senza filtro) e in compenso gli aggiunge 22 gradi di ritardo di fase. E'
esattamente il peggiore dei due mondi, ed e' la ragione per cui serve
`build_notch_filter.m`.

### Righe 18-19 -- l'approssimante di Pade

    [nd, dd] = pade(p.tau, padeOrder);   % tau = 0.020 s, ordine 3
    Wdelay   = tf(nd, dd);

**Da dove viene.** Il ritardo puro e' exp(-tau*s). L'approssimante di Pade
diagonale di ordine (n,n) e'

    exp(-tau*s) ~= N(s)/D(s),   con   N(s) = D(-s)

    D(s) = somma_{k=0..n} c_k * (tau*s)^k,
    c_k  = (2n-k)! * n! / [ (2n)! * k! * (n-k)! ]

Per n = 3 i coefficienti valgono c_0 = 1, c_1 = 1/2, c_2 = 1/10, c_3 = 1/120,
quindi

    D(x) = 1 + x/2 + x^2/10 + x^3/120,    N(x) = D(-x),    x = tau*s

**Perche' e' la forma giusta per un ritardo.** La proprieta' chiave e'
`N(s) = D(-s)`. Valutata sull'asse immaginario:

    N(j*w) = D(-j*w) = coniugato di D(j*w)      (coefficienti reali)

    => |N(j*w)| = |D(j*w)|   per OGNI w
    => |Wdelay(j*w)| = 1     esattamente, a tutte le frequenze

L'approssimante e' **all-pass**: modulo unitario ovunque, esattamente come il
ritardo vero (|exp(-j*w*tau)| = 1). Tutto quello che fa e' fase:

    arg(Wdelay(j*w)) = -2*arg(D(j*w))

che approssima -w*tau. **Questo e' il motivo per cui Pade e' il modo corretto di
modellare un ritardo**, e non un semplice "filtro passa-basso equivalente": un
lag del primo ordine 1/(1+tau*s) darebbe fase negativa *ma anche* attenuazione,
e l'attenuazione fittizia farebbe sembrare il loop piu' robusto di quanto sia
(nasconderebbe il picco di bending). Pade toglie fase e **non regala guadagno**.

**Fase non minima.** Poiche' N(s) = D(-s), gli zeri di N sono l'immagine
speculare dei poli di D rispetto all'asse immaginario. I poli di D sono stabili
(semipiano sinistro), quindi **gli zeri sono tutti nel semipiano destro**: il
blocco e' a **fase non minima**. Verificato numericamente per tau = 0.02 s,
ordine 3:

    zeri  (RHP):  +232.2  ,  +183.9 +/- 175.4j     [rad/s]
    poli  (LHP):  -232.2  ,  -183.9 +/- 175.4j     [rad/s]

Per l'ordine 1 il conto e' a mano: Wdelay = (1 - tau*s/2)/(1 + tau*s/2), zero in
s = +2/tau = +100 rad/s, polo in -100 rad/s. E' la firma classica del ritardo:
uno zero RHP costa fase senza dare attenuazione, ed e' l'unico modo di ottenere
il comportamento "fase che scende linearmente e senza limite" con una funzione
razionale.

**Accuratezza e scelta dell'ordine.** La regola pratica e' che Pade di ordine n
e' fedele fino a circa w*tau ~ n, cioe' w ~ n/tau. Con tau = 20 ms: ordine 1 ->
buono fino a ~50 rad/s, ordine 3 -> fino a ~150 rad/s. Verifica numerica della
fase (Pade-3 vs ritardo esatto -w*tau):

| w [rad/s] | Pade-3 | esatto -w*tau |
|-----------|--------|---------------|
| 2.445 (crossover rigido) | -2.80 gradi | -2.80 gradi |
| 18.9 (omega_BM)          | -21.66 gradi | -21.66 gradi |

**Nota di onesta'**: il commento alle righe 5-6 dice "higher = better phase lag
near wBM". E' vero in generale ma sovrastimato nel caso concreto: gia' Pade-1
sbaglia solo ~0.24 gradi a omega_BM (-21.42 contro -21.66). Il vero guadagno
dell'ordine 3 sta **sopra ~50 rad/s** -- dove Pade-1 sbaglia di 10 gradi al
corner dell'attuatore (70 rad/s: -69.98 contro -80.21 esatti) -- e nel calcolo
del **delay margin**, che `allmargin` legge su tutti i crossover, compresi quelli
ad alta frequenza. Cioe': l'ordine 3 serve, ma per una ragione diversa da quella
scritta nel commento.

### Riga 21 -- la cascata

    Wtvc = Wact * Wdelay;

Prodotto in serie: 2 stati (servo) + 3 stati (Pade-3) = **5 stati**. Verificato:
il loop aperto di Task 2 senza notch ha ordine 11 = 6 (plant full) + 2 + 3, senza
cancellazioni.

Il **bilancio di fase totale** dell'attuatore alle due frequenze che contano:

| w [rad/s] | servo | Pade-3 | totale TVC |
|-----------|-------|--------|------------|
| 2.445 (crossover rigido) | -2.80 | -2.80 | **-5.60 gradi** |
| 18.9 (omega_BM)          | -22.18 | -21.66 | **-43.84 gradi** |

Da leggere cosi': al crossover rigido il TVC costa poco (5.6 gradi su 30 di
phase margin). A omega_BM ruota il lobo di bending di **44 gradi** sulla carta di
Nichols -- e un lobo che sta a +29 dB sopra lo 0 dB, ruotato di 44 gradi, e' cio'
che decide se avvolgi o no il punto critico.

> **Possibile domanda d'esame** -- Perche' il ritardo di 20 ms e' cosi' velenoso
> per questo loop, che e' *condizionalmente stabile*?
> *Risposta:* Perche' un loop condizionalmente stabile deve **infilarsi** fra due
> punti critici sulla Nichols: la fase deve restare in una fascia sia alle basse
> frequenze (dove il polo aerodinamico instabile a +1.84 rad/s impone un guadagno
> minimo, l'aero gain margin) sia alle alte (dove il rolloff e il bending impongono
> un guadagno massimo). Il ritardo aggiunge una fase -w*tau che cresce
> **linearmente e senza limite** con la frequenza, e lo fa **a modulo invariato**.
> Non e' quindi "compensabile" abbassando il guadagno: attenuare non restituisce
> fase. Alle alte frequenze il ritardo trascina l'intera curva verso sinistra
> proprio dove il margine e' piu' stretto, e per questo il **delay margin**
> (`classify_margins.m`, riga 64: `min(am.DelayMargin)` su TUTTI i crossover) e' il
> numero di robustezza piu' stringente del progetto.

---

## Possibili domande d'esame

**D: Perche' approssimare il ritardo con Pade invece di usare `InputDelay` di
MATLAB, che sarebbe esatto?**
R: Perche' con un ritardo esatto il modello ad anello chiuso non e' a dimensione
finita: `pole(T)` e `isstable(T)` non sono definite su un sistema con ritardo
interno, e `main_task2.m` le usa continuamente (riga 34, righe 60-61 e riga 93,
dentro i loop dello Step C per scartare i candidati instabili). Anche il tuner
(`design_controller.m`, riga 89) chiama `isstable` a ogni valutazione di costo. Pade
rende il loop **razionale** e quindi rende eseguibile tutta la catena a valle.
Il prezzo e' che i margini sono accurati solo dove l'approssimante lo e': con
ordine 3 e tau = 20 ms, fino a ~150 rad/s -- ampiamente sopra il crossover
(~3 rad/s) e sopra il bending (18.9 rad/s), quindi il prezzo qui e' trascurabile.

**D: Cosa vuol dire che l'approssimante di Pade e' "all-pass", e perche' e'
esattamente la proprieta' che vogliamo?**
R: All-pass = modulo unitario a tutte le frequenze, |Wdelay(j*w)| = 1. Segue da
N(s) = D(-s): sull'asse immaginario N(j*w) e' il coniugato di D(j*w), quindi i
moduli si cancellano identicamente. E' la proprieta' del ritardo vero
(|exp(-j*w*tau)| = 1): un ritardo **non attenua nulla, sposta solo la fase**. Un
modello che invece attenuasse (per esempio un lag 1/(1+tau*s)) falserebbe il
progetto verso l'ottimismo: farebbe sembrare il picco di bending piu' basso di
quanto sia, e quindi il notch meno necessario.

**D: Perche' l'approssimante di Pade e' a fase non minima, e che effetto ha
sull'anello?**
R: Perche' i suoi zeri sono l'immagine speculare dei suoi poli (N(s) = D(-s)):
se i poli sono stabili, gli zeri sono nel semipiano destro. Per tau = 20 ms e
ordine 3 stanno a +232 rad/s e +184 +/- 175j rad/s. Uno zero RHP contribuisce
fase **negativa** con modulo crescente, cioe' fa esattamente cio' che fa un
ritardo. La conseguenza pratica e' un limite fondamentale di banda: non si puo'
spingere il crossover verso lo zero RHP senza perdere fase in modo irrecuperabile
(regola pratica: crossover ben sotto il modulo dello zero RHP). Qui il crossover
e' ~3 rad/s contro zeri a ~230 rad/s, quindi il vincolo non morde -- il ritardo
ci fa male attraverso il **bending** e il **delay margin**, non attraverso il
crossover rigido.

**D: L'attuatore TVC ha banda 70 rad/s. Non basta il suo rolloff a
gain-stabilizzare il bending a 18.9 rad/s?**
R: No, ed e' un errore comune. A 18.9 rad/s il servo del secondo ordine ha modulo
|Wact| = 0.9988, cioe' **-0.01 dB**: non attenua niente, perche' il suo corner e'
a 70 rad/s, un fattore 3.7 piu' in alto. In compenso gli regala -22 gradi di
fase, e il Pade altri -22. Quindi l'attuatore, sul bending, fa **solo danni**:
lascia passare i +29 dB del picco e ci aggiunge 44 gradi di ritardo. Il modo
va attenuato con un filtro dedicato (`build_notch_filter.m`) oppure
fase-stabilizzato -- l'attuatore non e' un'opzione.

**D: Il ritardo di 20 ms rappresenta solo il trasporto dell'attuatore?**
R: Nella traccia e' dato come ritardo puro dell'Eq. 3, e il codice lo prende cosi'
(`p.tau = 0.020` in `load_hw3_params.m`, riga 57). Fisicamente, in un progetto di
volo, quel numero e' un aggregato: latenza di calcolo del flight computer,
filtraggio dei sensori, campionamento (uno ZOH a periodo T aggiunge un ritardo
efficace ~T/2), e il trasporto vero dell'attuatore. Il codice non modella nulla
di tutto questo separatamente -- e' un singolo `tau` a valle del controllore.

**D: Cosa succede se metto `padeOrder = 1`?**
R: La funzione accetta (il validatore chiede solo intero positivo) e il progetto
non esplode: a omega_BM l'errore di fase e' solo 0.24 gradi. Cambiano invece i
numeri ad alta frequenza -- al corner dell'attuatore l'errore e' 10 gradi -- e
quindi cambiano il **delay margin** (che `allmargin` legge anche sui crossover
alti) e la posizione dei crossover di guadagno sopra i 100 rad/s. In pratica: i
margini pubblicati nel report sono legati alla scelta `padeOrder = 3` fatta a
`main_task2.m` riga 30, ed e' onesto dichiararlo.
