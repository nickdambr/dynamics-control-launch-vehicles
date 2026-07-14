# HM3/LTV_FULL_ASCENT/main_q_scheduling.m

## Ruolo del file nel progetto

`main_full_ascent.m` schedula i guadagni sul **tempo di volo**. E' comodo ma e' barare:
un computer di bordo non "misura" il tempo dall'accensione come parametro fisico -- o
meglio, lo conosce, ma se la missione devia dal profilo nominale (ritardo di accensione,
spinta fuori tolleranza, vento che rallenta la salita) il tempo non e' piu' una buona
etichetta per lo stato del veicolo. La scelta canonica in letteratura LPV per un lanciatore
e' schedulare su una grandezza **misurabile e fisicamente legata al plant**, e il candidato
naturale e' la **pressione dinamica `q_bar`** (che il veicolo stima da pressione statica e
velocita' relativa).

Questo script (ticket T008, Goal 2) prende **gli stessi identici guadagni** progettati da
`design_controller` sulla griglia temporale e li **ri-indicizza su `q`**, poi chiede: `q`
e' una buona variabile di scheduling per questo veicolo? La risposta del codice -- scritta
in chiaro nell'header, righe 7-11 -- e' **no, non pulitamente**, e lo script quantifica il
perche': `A6(t)` (l'instabilita' aerodinamica) ha il massimo a t ~ 72 s mentre `Q(t)` ha il
massimo prima (t = 65-67 s), quindi la mappa guadagno-vs-`q` e' **isteretica**: il ramo
crescente e quello calante di `q` chiedono guadagni diversi allo stesso `q`.

Struttura: costruisce la lookup `Kp(q)` dal solo **ramo ascendente** (dove `q` e'
monotona), misura l'isteresi sul ramo discendente, integra tre risposte LTV (frozen,
t-scheduled, q-scheduled) sullo stesso plant e sullo stesso vento, e fa uno sweep di
margini congelati confrontando t-schedule e q-schedule. Dipende da `init_simulink_lpv`,
`ode_lpv_ascent`, `build_plant_rigid`, `load_hw3_params`, `assemble_loop`. Il plant e'
quello **rigido**: lo scheduling dei guadagni e' indipendente dal bending (che e' il
Goal 1, in `main_flex.m`).

Anche qui vale l'avvertenza generale della cartella: **i numeri del README non si
riproducono con il codice attuale** (l'isteresi dichiarata e' 41%, oggi esce 125%; il
crollo di margine dichiarato "0.5 dB / 2 deg attorno a t = 105 s" oggi non si legge li').
La causa e' la stessa spiegata in `lpv_main_full_ascent.md` (il refactor di
`design_controller`). La **conclusione qualitativa dello studio, invece, regge** -- e la
regge in modo anche piu' netto.

---

## Intestazione e preambolo (righe 1-21)

- Righe 1-13: l'header dichiara gia' il risultato e la raccomandazione ("Mach (monotonic)
  would be the better measurable"). E' un buon esempio di come si scrive un'analisi
  negativa: si dice cosa non funziona e si propone l'alternativa.
- Riga 16: `warning('off', ...MarginUnstable)` -- come in `main_full_ascent`, l'anello e'
  condizionalmente stabile.
- Righe 20-21: `S = init_simulink_lpv()` rifa tutto il setup (dati LPV, 28 design PD sulla
  griglia `5:5:140`, simulazione del generatore di vento). Nota: lo script **ri-esegue** il
  generatore e i 28 tuning, non li riusa da `main_full_ascent` -- ogni script e'
  autoconsistente, al costo di ripetere il lavoro.

---

## `%% Build a q-keyed gain lookup` (righe 23-37)

```matlab
Qs = S.fQ(S.tsched)/1000;          % q ai punti di schedule [kPa]
[~, ipk] = max(Qs);                % il picco separa i due rami
qa  = Qs(1:ipk);  Kpa = S.Kp_sched(1:ipk);  Kda = S.Kd_sched(1:ipk);
Kp_of_q = griddedInterpolant(qa, Kpa, 'linear', 'nearest');
Kd_of_q = griddedInterpolant(qa, Kda, 'linear', 'nearest');
fKp_q = @(t) Kp_of_q(S.fQ(t)/1000);
fKd_q = @(t) Kd_of_q(S.fQ(t)/1000);
```

- Riga 24: valuta `q` **ai punti di design** (28 valori, `t = 5:5:140`).
- Riga 25: trova l'indice del picco. Verificato: `ipk = 13`, cioe' **t = 65 s**, `q = 43.29
  kPa` (il vero massimo del dataset e' a t = 67 s con 43.9 kPa, ma sulla griglia a 5 s il
  massimo cade a 65 s). Le etichette delle legende alle righe 82-83 dicono infatti
  "t <= 65 s" / "t >= 65 s": sono **hardcoded** e coerenti solo perche' il picco cade li'
  con questi parametri di default; se si cambia `tsched_step` diventano sbagliate.
- Riga 26: taglia il **ramo ascendente** (da t = 5 a t = 65 s). Verificato: su quel ramo `q`
  e' strettamente monotona (da 0.39 a 43.29 kPa), che e' la condizione necessaria perche'
  `griddedInterpolant` accetti `qa` come griglia (richiede nodi crescenti). Fuori da li'
  la mappa `t -> q` **non e' invertibile**, ed e' esattamente questo il problema.
- Righe 27-28: la lookup e' costruita **solo** sul ramo ascendente, con estrapolazione
  `'nearest'` (clip agli estremi).
- Righe 30-31: la composizione `fKp_q(t) = Kp_of_q(q(t))` e' la simulazione onesta di cio'
  che farebbe il computer di bordo: **misura q, cerca il guadagno in tabella**. Il fatto
  che qui `q` sia riottenuta da `t` e' solo un artificio di simulazione (il veicolo la
  misurerebbe da sensori); la logica di controllo, pero', vede solo `q`.

**L'isteresi (righe 33-37).**

```matlab
Kp_cmd_dn  = fKp_q(S.tsched(ipk+1:end));    % cio' che la lookup comanda
Kp_need_dn = S.Kp_sched(ipk+1:end);         % cio' che quel ramo richiede
hyst = max(abs(Kp_cmd_dn - Kp_need_dn) ./ Kp_need_dn) * 100;
```

La metrica: sul ramo discendente, confronta il guadagno che la lookup-in-`q` comanda con
il guadagno che il design a quel *tempo* aveva prodotto. Se `q` fosse una buona variabile
di scheduling i due coinciderebbero. Verificato eseguendo: **`hyst` = 125.4 %**, massimo
raggiunto a **t = 140 s**.

---

## Perche' `q` non basta -- la derivazione (il cuore dello studio)

Vale la pena farla per bene, perche' e' *la* domanda d'orale su questo file.

L'argomento **a favore** di `q` e' quello da manuale: il coefficiente destabilizzante e'

    A6 = mu_alpha = N_alpha * l_alpha / Iyy ,   N_alpha ~ q_bar * S_ref * C_N_alpha

quindi il termine aerodinamico "scala con la pressione dinamica" e `q` sembra il parametro
di scheduling naturale. E' inoltre **misurabile a bordo** (pressione statica + velocita'
relativa), a differenza del tempo, che e' legato a un timeline nominale specifico.

L'argomento e' pero' **incompleto su due punti**, ed entrambi si verificano sul dataset.

**Punto 1: `A6` non e' proporzionale a `q`.** In `N_alpha ~ q_bar * S * C_N_alpha(Mach)`
c'e' la dipendenza da Mach del coefficiente aerodinamico, e nel denominatore c'e' `Iyy`,
che **cala** con il consumo di propellente. Verificato sul dataset (rapporto `A6 / q_bar`,
1/(s^2 kPa)):

| t [s] | 20 | 40 | 72 | 95 | 105 |
|---|---|---|---|---|---|
| A6 / q_bar | 0.0438 | 0.0438 | 0.0797 | 0.1045 | 0.1453 |

Il rapporto **triplica**: `A6` e `q` non sono la stessa informazione. E' precisamente per
questo che il massimo di `A6` (t = 72 s) e' 5 s dopo il massimo di `q` (t = 67 s) -- e da
qui l'isteresi.

**Punto 2 (il piu' forte, e il piu' spesso dimenticato): `K1` non c'entra nulla con `q`.**
L'efficacia di controllo e'

    K1 = mu_c = Tc * l_c / Iyy

cioe' **spinta per braccio diviso inerzia**: e' una grandezza **propulsivo-inerziale**, non
aerodinamica. Non dipende da `q_bar` in alcun modo. Verificato: `K1` cresce
**monotonicamente** da 3.21 a 9.64 1/s^2 sull'orizzonte usato (le tabelle sono troncate a
`Tstop` = 140 s; il dataset prosegue fino a t = 150 s, dove `K1` arriva a 10.84, ma quel
punto e' fuori dalla simulazione), mentre `q` sale e poi scende. E il guadagno canonico
**da cui il tuner parte** e'

    Kp = 2*A6/K1        Kd = sqrt(A6)/K1

(`design_controller`, riga 55: `x0 = log([2*A6/K1, sqrt(A6)/K1])`; l'obiettivo effettivo
della cost, riga 88, e' il matching dei margini sull'anello completo -- ma dove la
ritaratura non trova l'attraversamento da agganciare i guadagni **restano** il seed, ed e'
il caso di buona parte della griglia, vedi `lpv_main_full_ascent.md`). Quel guadagno
dipende da **entrambi** i coefficienti. Una lookup in `q` e' cieca su `K1`.

La prova numerica piu' pulita, letta direttamente sui dati:

| t [s] | q_bar [kPa] | A6 | K1 | Kp = 2A6/K1 |
|---|---|---|---|---|
| 30 (salita) | 10.57 | 0.479 | 3.68 | 0.260 |
| 140 (fine) | 10.59 | 0.558 | 9.64 | 0.116 |

**Stessa identica pressione dinamica (10.6 kPa), `A6` quasi uguale, ma il guadagno
richiesto e' 2.25 volte diverso** -- interamente per colpa di `K1`, che `q` non puo'
vedere. Ed e' esattamente il punto dove lo script misura il massimo di isteresi: la lookup
comanda `Kp = 0.2607` mentre servirebbe `Kp = 0.1157` (**+125.4 %**, sovra-guadagno).

Il difetto duale, il **sotto**-guadagno, si legge a t = 105 s: li' `q = 13.53 kPa`, che sul
ramo ascendente corrisponde a t ~ 33 s, dove pero' `A6` valeva 0.589 contro 1.965 a
t = 105 s (3.3 volte meno instabile). La lookup comanda `Kp = 0.3221` quando servirebbe
`Kp = 0.8772`: **63 % di sotto-guadagno** proprio nella discesa ad alto `A6` -- che e' la
frase del README, e questa parte si riproduce.

**Conclusione corretta.** Serve un parametro di scheduling `rho(t)` che sia (a) misurabile,
(b) **monotono** lungo la traiettoria (altrimenti la mappa `rho -> plant` non e' una
funzione), e (c) sufficientemente informativo sul plant. `q` fallisce (b) e in parte (c).
Il **Mach** e' monotono (verificato: da 0.03 a 8.34, `all(diff(Mach)) > 0`) ed e'
misurabile: e' il candidato corretto, ed e' quello che il README indica come follow-up. In
alternativa si schedula su **due** parametri (`q` e Mach, o `q` e massa), che e' cio' che si
fa nei lanciatori veri.

> **Possibile domanda d'esame** -- "ma allora perche' i libri dicono di schedulare su q?"
> *Risposta:* perche' nella maggior parte dei casi si schedula la **parte aerodinamica** del
> problema (per esempio un load-relief o il guadagno di un anello che deve contrastare
> `N_alpha`) e su un intervallo di volo in cui `q` e' monotona. Qui invece si schedula il PD
> completo, il cui guadagno dipende anche da `K1 = Tc*l_c/Iyy`, e su un orizzonte che
> attraversa il massimo di `q`. Le due ipotesi implicite del "manuale" saltano entrambe.

---

## `%% Closed-loop responses` (righe 39-50)

```matlab
runs = struct('name', {'frozen', 't-sched', 'q-sched'});
mdls = {make_model(S, 'frozen', fKp_q, fKd_q), ...
        make_model(S, 't',      fKp_q, fKd_q), ...
        make_model(S, 'q',      fKp_q, fKd_q)};
for k = 1:3
    [~, x] = ode45(@(t,x) ode_lpv_ascent(t,x,mdls{k}), tt, zeros(4,1), odeo);
    runs(k).r = unpack(tt, x, S, mdls{k});
end
```

Tre integrazioni sullo **stesso** plant LTV, **stesso** vento (`S.windfun`, il generatore
del professore), **stesso** stato iniziale nullo a t0 = 5 s, tolleranze 1e-8/1e-10. Cambia
solo la sorgente dei guadagni di beccheggio. I guadagni di deriva `Kp_z = Kd_z = -1e-3`
sono sempre quelli di `S.K0`, in tutti e tre i casi (righe 137-138): **la deriva non e' mai
schedulata**.

Verificato eseguendo (picchi su tutta l'ascesa):

| controllore | picco `\|theta\|` | picco `\|z\|` | picco `\|delta\|` |
|---|---|---|---|
| frozen | 0.971 deg | 29.5 m | 1.06 deg |
| t-sched | 25.36 deg | 36.2 m | 3.44 deg |
| q-sched | **530.9 deg** | **1517 m** | 59.4 deg |

La corsa `q-sched` **diverge**. Attenzione a come si legge: una parte della divergenza e'
ereditata dal difetto del tuner descritto in `lpv_main_full_ascent.md` (lo schedule di
partenza e' gia' cattivo per `t <= 35 s`), e la lookup in `q` lo **amplifica**, perche' a
`q` bassa (inizio volo) restituisce guadagni quasi nulli e li ri-applica anche in coda al
volo. La conclusione qualitativa (`q` non e' una buona variabile di scheduling) resta
valida e anzi e' rafforzata; ma la magnitudine (531 deg) **non e' un risultato pulito** e
non va presentata come tale.

---

## `%% Frozen-time margin sweep` (righe 52-61)

Stesso schema di `main_full_ascent` (e stesso `loop_margin` locale con `margin()` invece di
`classify_margins`, e stesso `abs()` nei plot -- **stessi due difetti**, vedi
`lpv_main_full_ascent.md`): congela il plant ogni 2.5 s e legge GM/PM con i guadagni
`t-sched` e con quelli `q-sched`.

Verificato, valori `|GM|` in dB (t-sched vs q-sched):

| t [s] | 45 | 65 | 75 | 85 | 95 | 105 |
|---|---|---|---|---|---|---|
| t-sched | 6.00 | 6.00 | 6.00 | 6.00 | 6.00 | 6.00 |
| q-sched | 6.00 | 6.00 | **3.72** | **2.75** | **1.71** | 8.03 (*) |

(*) da t ~ 100 s in poi la lettura `margin()` cade su un'attraversata di deriva e il numero
non e' un margine di corpo rigido: e' l'artefatto discusso sopra, che l'`abs()` fa sembrare
ottimo. **Il vero segnale e' la banda 75-95 s**: mentre il t-schedule tiene i 6 dB di
progetto, il q-schedule li perde progressivamente (3.7 -> 2.7 -> 1.7 dB) proprio nella
discesa ad alto `A6`, cioe' dove la lookup sotto-guadagna. Questo e' il grafico che
dimostra la tesi dello studio.

Il README dichiara invece "il margine crolla da 6 dB / 30 deg a 0.5 dB / 2 deg attorno a
t = 105 s". Oggi non si legge cosi': il minimo del q-schedule esce **0.83 dB a t = 35 s**
(dove pero' anche il t-schedule e' rotto, 0.80 dB a t = 120 s), quindi quel minimo **non e'
attribuibile a `q`**. La degradazione attribuibile a `q` in modo pulito e' quella di
75-95 s.

---

## `%% Summary` (righe 63-74)

Stampa il tempo del picco di `Q`, l'isteresi in %, e per ciascuna delle tre corse i picchi
piu' il minimo di `|GM|`/`|PM|` (solo per t-sched e q-sched: alla riga 69 `g` e `p` sono
inizializzati a `NaN` e riempiti solo per `k = 2, 3`, quindi la riga "frozen" mostra NaN
nelle colonne dei margini -- **non e' un bug**, e' voluto: lo sweep dei margini frozen sta
nell'altro script).

---

## `%% Figures` (righe 76-115) e `%% Export` (righe 117-124)

- f1 (righe 78-90): **il grafico chiave**. Plotta `Kp` e `Kd` **contro `q`**, separando in
  blu il ramo `q` crescente e in rosso quello calante. Se `q` fosse una buona variabile di
  scheduling, le due curve si sovrapporrebbero; invece formano un **ciclo di isteresi**. E'
  la rappresentazione visiva di "la mappa `t -> q` non e' invertibile".
- f2 (righe 93-102): risposte t-sched vs q-sched (theta, z, delta). Nota che la cell array
  `sig` (riga 96) usa i nomi dei campi di `r` per generare i tre riquadri in loop: elegante,
  ma implica che `unpack` **deve** produrre campi con quei nomi esatti.
- f3 (righe 105-115): sweep dei margini con i target 6 dB / 30 deg.

---

## `make_model` (righe 127-144)

```matlab
switch mode
    case 'frozen', M.sched = false; M.fKp = S.fKp;  M.fKd = S.fKd;
    case 't',      M.fKp = S.fKp;   M.fKd = S.fKd;
    case 'q',      M.fKp = fKp_q;   M.fKd = fKd_q;
end
```

Il campo `sched` e' inizializzato a `true` alla riga 138 e messo a `false` solo nel caso
`'frozen'`. `ode_lpv_ascent` guarda `M.sched`: se falso usa `Kp_th0/Kd_th0` (i frozen di
max-q), se vero valuta `M.fKp(t)`. Quindi la **stessa** funzione RHS serve i tre casi: nel
caso `'q'` la "funzione del tempo" `fKp_q` e' in realta' una composizione
`Kp_of_q(q(t))` -- il RHS non se ne accorge. E' un'astrazione pulita: il modello di
scheduling e' completamente incapsulato in un handle `t -> guadagno`.

---

## `unpack` (righe 146-161)

Come in `main_full_ascent`, ma piu' snello (niente `alpha` ne' `q_bar*alpha`). Ricostruisce
`delta` a posteriori con i guadagni della corsa (riga 157) -- lecito perche' l'attuatore e'
ideale. Nota alla riga 155: legge `M.sched`, `M.fKp`, `M.fKd` **dal modello**, non da un
flag esterno: cosi' i guadagni usati per ricostruire `delta` sono garantiti identici a
quelli usati dentro l'ODE. E' un dettaglio di igiene che evita l'errore classico di
plottare un `delta` incoerente con la simulazione.

---

## `loop_margin` (righe 163-174)

Identica a quella di `main_full_ascent` (codice duplicato fra i due script -- si sarebbe
potuta estrarre in una funzione condivisa). Vedi li' per la discussione su `margin()` vs
`classify_margins`.

---

## Possibili domande d'esame

**D: Perche' si vorrebbe schedulare su `q_bar` invece che sul tempo?**
R: Perche' il tempo di volo non e' un parametro *fisico* del sistema: e' l'etichetta di un
timeline nominale. Se la missione si discosta dal nominale (spinta fuori tolleranza, vento
che rallenta la salita, accensione ritardata) il veicolo a t = 72 s non e' piu' nello stato
per cui i guadagni erano stati progettati. Un parametro misurato a bordo -- `q_bar`, Mach,
massa stimata -- **segue lo stato reale del veicolo** e rende lo scheduling robusto alle
dispersioni di profilo. E' l'idea base del gain scheduling / LPV: schedulare su `rho(t)`
misurato, non sull'indice temporale.

**D: E allora perche' il tuo studio conclude che `q_bar` qui e' una scelta cattiva?**
R: Per due motivi verificati sui dati. (1) `q(t)` **non e' monotona**: ha il massimo a
t = 67 s, quindi la mappa `t -> q` non e' invertibile e a uno stesso `q` corrispondono due
istanti con plant diversi -- la tabella dei guadagni diventa **isteretica**. (2) Il guadagno
richiesto e' `Kp = 2*A6/K1`, e mentre `A6` e' aerodinamico (legato a `q`, ma non
proporzionale: `A6/q` triplica lungo il volo per effetto del Mach e della caduta di `Iyy`),
`K1 = Tc*l_c/Iyy` e' **propulsivo-inerziale** e cresce monotonamente da 3.2 a 9.6
sull'orizzonte simulato (10.8 all'ultimo punto del dataset, t = 150 s, fuori
dall'orizzonte) -- `q` e' completamente cieca su di esso. Esempio numerico: a t = 30 s e a t = 140 s il veicolo ha la
**stessa** `q_bar` (10.6 kPa) ma richiede guadagni 2.25 volte diversi.

**D: Come misuri l'isteresi nel codice, e quanto vale?**
R: Costruisco la lookup `Kp(q)` sul solo ramo ascendente (dove `q` e' monotona, t = 5-65 s),
poi sul ramo discendente confronto il guadagno che quella lookup comanda con quello che il
design a quel tempo aveva prodotto: `hyst = max(|Kp_cmd - Kp_need| / Kp_need)`. Verificato
oggi: **125.4 %**, con il massimo a t = 140 s. (Il README riporta 41 %: e' un numero
antecedente al refactor di `design_controller`.)

**D: Qual e' la variabile di scheduling corretta, allora?**
R: Il **Mach**: e' misurabile a bordo, e' monotono lungo l'ascesa (verificato sul dataset:
da 0.03 a 8.34, derivata sempre positiva) e cattura la fisica che manca a `q` (la
dipendenza da Mach di `C_N_alpha`). Non risolve pero' da solo il problema di `K1`: la
soluzione industriale e' uno scheduling a **due parametri** (per esempio `q_bar` e Mach, o
`q_bar` e massa), oppure schedulare su una variabile derivata come `A6/K1` stessa se e'
stimabile.

**D: Che differenza c'e' fra gain scheduling e vero controllo LPV?**
R: Il gain scheduling classico e' quello che fa questo codice: si progetta un controllore
LTI su ogni punto congelato e si **interpola** fra i guadagni; la stabilita' del sistema
tempo-variante risultante non e' garantita da nessun teorema, va verificata a posteriori
(ed e' esattamente cio' che fa l'integrazione ode45 dell'LTV). Il controllo LPV vero e
proprio progetta un unico controllore parametrizzato in `rho` con garanzie di stabilita'
sull'intera famiglia (via LMI / funzioni di Lyapunov parameter-dependent), tipicamente
imponendo un limite sulla velocita' di variazione `rho_dot`. Il primo e' quello che si usa
sui lanciatori reali; il secondo e' la copertura teorica.

**D: Se la corsa `q-sched` diverge, come fai a dire che il problema e' l'isteresi e non il
tuner rotto?**
R: Non lo dico dalla risposta temporale, che e' contaminata (lo schedule di base e' gia'
degradato per `t <= 35 s`). Lo dico dallo **sweep di margini congelati nella banda
75-95 s**, dove il t-schedule tiene esattamente i 6 dB di progetto mentre il q-schedule
scende a 3.7 -> 2.7 -> 1.7 dB. In quella banda entrambi gli schedule sono ben tarati, quindi
la differenza e' interamente attribuibile alla **ri-indicizzazione in `q`**: e' li' che si
misura il costo dell'isteresi, non nella divergenza temporale.
