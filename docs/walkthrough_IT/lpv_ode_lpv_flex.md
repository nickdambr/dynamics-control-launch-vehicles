# HM3/LTV_FULL_ASCENT/ode_lpv_flex.m

## Ruolo del file nel progetto

E' la RHS della **versione flessibile** dell'estensione LPV: stesso plant
rigido di `ode_lpv_ascent.m`, ma con **il primo modo di bending riattivato**,
la contaminazione INS delle misure, la catena di attuazione reale (TVC di
secondo ordine + ritardo di trasporto approssimato con Pade) e -- il punto
dello studio -- un **notch che puo' inseguire omega_BM(t)**.

Il problema che risolve. In HM3 il notch profondo del Task 2 e' centrato su
`omega_BM(72) = 18.9 rad/s`, la frequenza del primo modo **al punto di
progetto**. Ma omega_BM non e' una costante: dipende da massa e rigidezza
generalizzate, che cambiano mentre il primo stadio brucia propellente. Sul
dataset la frequenza sale (il README riporta uno sweep 16.5 -> 31.8 rad/s):
un notch fisso a 18.9 rad/s **si detuna**, la risonanza esce dalla buca, il
guadagno d'anello alla frequenza di bending torna sopra 0 dB e l'anello va
instabile. Questo file permette di simulare i due casi -- notch fisso e notch
variabile -- **sullo stesso identico plant LTV**, cambiando una sola cosa: il
handle `M.fwn`.

Dove sta nel flusso: e' il gemello dinamico di `HM3/build_plant_full.m` (che
costruisce lo stesso modello a 6 stati ma congelato), viene chiamato da
`HM3/LTV_FULL_ASCENT/main_flex.m` (righe 40-41, due `ode45`: `S.fomega` per il
notch variabile, `@(t) w72` per quello fisso) e riceve la struttura `M` dalla
helper locale `make_flex` (righe 93-106 dello stesso file), riempita con gli
interpolanti creati da `init_simulink_lpv.m` (righe 104-117). Il modello
Simulink `hm3_full_ascent_flex.slx` (`build_hm3_full_ascent_flex.m`) e' il suo
mirror a blocchi elementari e viene validato contro di lui, non viceversa.

Nota di scoping: qui i **guadagni PD restano congelati** a max-q
(`make_flex` passa `'sched', false`, riga 105 di `main_flex.m`). Lo showcase e'
il notch, non lo scheduling. Il codice *supporta* `sched = true` (riga 26) ma
nessuno script lo esercita in questa configurazione.

---

## Firma, stato e contratto (righe 1-16)

```matlab
function dx = ode_lpv_flex(t, x, M)
% x = [z zdot theta thetadot eta etadot | xn1 xn2 | x_tvc(1:nt)]
% No arguments validation by design: ode45 inner loop.
```

- Riga 1: firma `(t, x, M)`, stessa convenzione di `ode_lpv_ascent`.
- Righe 5-9: lo stato e' **concatenato per blocchi**, e questa e' la chiave di
  lettura di tutto il file:
  - `x(1:6)` = plant flessibile `[z zdot theta thetadot eta etadot]`;
  - `x(7:8)` = stati del **notch** (realizzazione in forma canonica);
  - `x(9:end)` = stati del **TVC** (attuatore + ritardo di Pade).
- La docstring dichiara 13 stati: 6 + 2 + nt, con nt = 5 perche' il TVC e'
  `build_tvc(p0, 3)` = attuatore di 2 ordine per un Pade di ordine 3
  (`HM3/build_tvc.m`, righe 16-21). `main_flex.m` non lo hard-codea: legge
  `nt = size(S.tvc.At, 1)` (riga 37) e dimensiona `x0` di conseguenza (riga 38).
- Riga 13: di nuovo **nessun blocco `arguments`, per design** -- hot loop di
  `ode45`. Vale identico a `ode_lpv_ascent`.
- Riga 15: lo split (`xp`, `xn`, `xt`). Riga 16: destrutturazione in nomi
  parlanti. Costa qualche copia in piu' per chiamata, ma rende leggibili le
  righe 38-40; e' un compromesso deliberato fra velocita' e chiarezza.

---

## Coefficienti tempo-varianti: qui si interpola il **grezzo** (righe 19-21)

```matlab
a1 = M.fa1(t); a3 = M.fa3(t); a4 = M.fa4(t);
A6 = M.fA6(t); K1 = M.fK1(t); V = M.fV(t);
w  = M.fomega(t); aqk = M.faqk(t);
sig = M.fsig(t); phi = M.fphi(t);
aw = M.windfun(t);
```

- Righe 19-20: **dieci** lookup per valutazione (contro le sette di
  `ode_lpv_ascent`). Sono tutti `griddedInterpolant(tg, y, 'linear',
  'nearest')` -- la stessa closure `gi`, riga 101 di `init_simulink_lpv.m` --
  costruiti alle righe 110-111 (`fV` fa eccezione: nasce alla riga 104, fuori
  dal blocco flex, insieme a `fQ`) sulla
  griglia del dataset `GreensiteLPV_DATA.mat`. Interpolazione **lineare**
  dentro la griglia, **clamp all'estremo** fuori ('nearest' come metodo di
  estrapolazione): nessun errore, nessun NaN, nessun warning se si esce.
- **Differenza sostanziale rispetto a `ode_lpv_ascent`**: li' venivano
  interpolati i coefficienti **gia' combinati** (`fc2` = tabella di
  `a1*V + a4`, `fc5` = tabella di `A6/Vsafe`), qui si interpolano i
  **coefficienti grezzi** e le combinazioni si formano **dentro** la RHS
  (righe 38-39). Poiche' l'interpolazione lineare non commuta con la
  moltiplicazione, fra due breakpoint i due file integrano sistemi LTV
  **leggermente diversi** (differiscono di un O(dt^2) proporzionale alle
  variazioni dei fattori; coincidono esattamente sui nodi). Nessuno dei due e'
  sbagliato, ma non ci si deve aspettare che le uscite coincidano bit a bit.
- Righe 19: `V = M.fV(t)` e' l'interpolante della **V grezza**, senza la
  guardia `Vsafe = max(V,1)` che `init_simulink_lpv.m` applica al ramo rigido
  (riga 59). Alla riga 39 si divide per `V`. Al lift-off V(0) = 0, quindi
  **questa RHS non e' integrabile da t = 0**: la divisione esploderebbe.
  Funziona solo perche' `main_flex.m` parte da `t0 = S.t0` (riga 19), e quel
  valore vale 5 s (default dell'`arguments` block di `init_simulink_lpv.m`,
  riga 34). E' una
  fragilita' reale, non protetta da nessun controllo; nel ramo rigido la
  guardia c'e', qui no.
- Riga 21: il vento, `alpha_w(t)` [rad], stesso `griddedInterpolant` del ramo
  rigido (`init_simulink_lpv.m`, riga 107), gia' diviso per V a monte
  (riga 90). Il generatore del professore gira **una volta sola**, quindi le
  due corse (notch fisso e variabile) vedono **la stessa realizzazione di
  vento**: il confronto e' pulito.

---

## Misure INS contaminate dal bending (righe 24-25)

```matlab
theta_m = theta + sig*eta;
thetadot_m = thetadot + sig*etadot;
z_m     = z - phi*eta;
zdot_m  = zdot - phi*etadot;
```

- Righe 24-25: e' letteralmente l'Eq. (2) della traccia, la stessa matrice
  `Cm` di `HM3/build_plant_full.m` (righe 33-36): `sigma_ins` [rad/m] pesa la
  **pendenza** della forma modale nella posizione della IMU (contamina
  l'assetto misurato), `phi_ins` [-] pesa lo **spostamento** modale (contamina
  la posizione misurata). I segni sono opposti (+sig, -phi) e replicano
  esattamente il modello congelato.
- **Qui e' dove nasce tutto il problema del bending**. Il modo elastico non
  entra nella dinamica rigida (righe 38-39: nessun `eta`); entra **solo**
  attraverso la misura. Il PD retroaziona `theta_m`, che contiene `sig*eta`;
  il comando `delta` che ne risulta eccita il modo elastico (riga 40,
  `aqk*delta`), che a sua volta contamina di piu' la misura: e' un anello
  chiuso attorno a un modo con smorzamento `zBM = 0.005`, cioe' praticamente
  nullo. Senza notch il picco di risonanza in `L(j*omega_BM)` sta sopra 0 dB
  e l'anello e' instabile. Da qui il notch.
- Attenzione: `sig` e `phi` sono **anch'essi tempo-varianti** (interpolati
  alla riga 20 dal dataset, campi `sigma_ins` e `phi_ins`). Non e' solo la
  frequenza a muoversi lungo l'ascesa: anche il *guadagno* di contaminazione.
  Il codice ne tiene conto senza dirlo esplicitamente.

---

## Il notch tempo-variante (righe 30-33)

```matlab
wn = M.fwn(t);
v  = 2*(M.zN - M.zD)*wn*xn(2) + u_pd;
xn1d = xn(2);
xn2d = -wn^2*xn(1) - 2*M.zD*wn*xn(2) + u_pd;
```

**Derivazione.** Il notch di HM3 (Eq. 4, variante a fase minima usata come
design ritenuto -- `main_task2.m`, righe 73-78, con `zN = 0.002`,
`zD = 0.7`, `sgn = +1`) e':

    N(s) = (s^2 + 2*zN*wn*s + wn^2) / (s^2 + 2*zD*wn*s + wn^2)

Realizzazione in **forma canonica di controllabilita'**, con denominatore
`s^2 + alpha_1*s + alpha_0` (alpha_1 = 2*zD*wn, alpha_0 = wn^2) e numeratore
`b_2*s^2 + b_1*s + b_0` (b_2 = 1, b_1 = 2*zN*wn, b_0 = wn^2):

    xn1' = xn2
    xn2' = -alpha_0*xn1 - alpha_1*xn2 + u
    y    = (b_0 - alpha_0*b_2)*xn1
           + (b_1 - alpha_1*b_2)*xn2 + b_2*u

Sostituendo: `b_0 - alpha_0 = wn^2 - wn^2 = 0` (il coefficiente di `xn1`
sparisce, ed e' il motivo per cui alla riga 31 `xn(1)` non compare) e
`b_1 - alpha_1 = 2*(zN - zD)*wn`. Il feedthrough e' `b_2 = 1`, quindi il
guadagno in continua e ad alta frequenza vale 1: il notch **non tocca** la
banda utile, agisce solo attorno a `wn`. Il codice alle righe 31-33 e'
esattamente questo, senza semplificazioni nascoste.

- Con `zN = 0.002` e `zD = 0.7`, il coefficiente d'uscita e' circa `-1.4*wn`:
  buca profondissima e strettissima (attenuazione ~ zN/zD, circa -51 dB al
  centro). E' un notch **molto** piu' aggressivo delle linee guida della
  traccia (che nella docstring di `build_notch_filter.m`, righe 6-7,
  suggerisce zN in 0.1-0.3, zD in 0.4-0.6): scelta deliberata di HM3
  ("deep, narrow null"), ereditata qui identica. Un notch cosi' stretto e'
  precisamente quello che **si detuna di piu'**: il fatto che l'esperimento
  T008 mostri la divergenza e' in parte una conseguenza di questa scelta di
  progetto, non solo della fisica.
- Riga 30: `wn = M.fwn(t)`. **E' l'unico bottone dell'esperimento.**
  `main_flex.m` passa `S.fomega` (il `griddedInterpolant` di omega(t)) per il
  notch variabile e `@(t) w72` (una anonymous function costante) per quello
  fisso: la RHS non sa e non deve sapere quale delle due riceve -- duck typing,
  entrambe si chiamano con `f(t)`. Elegante, e rende il confronto rigoroso
  perche' **nient'altro cambia**.
- Solo la variante a **fase minima** e' realizzabile con queste righe: il
  numeratore ha `+2*zN*wn*s`. La variante di `build_notch_filter` con
  `numSign = -1` (Eq. 4 come stampata sulla traccia, zeri in semipiano destro)
  **non e' rappresentabile** senza cambiare il segno alla riga 31. Limite del
  file, da dichiarare.

### Il caveat teorico che vale un punto all'orale

Sostituire `wn -> wn(t)` **dentro una realizzazione** non e' un'operazione
innocua. Una funzione di trasferimento ha infinite realizzazioni equivalenti;
sono equivalenti **solo per parametri costanti**. Se i coefficienti variano nel
tempo, realizzazioni diverse della *stessa* N(s) congelata producono sistemi
LTV **diversi** (i cambi di base diventano `T(t)`, e la trasformazione di
similarita' genera termini in `Tdot`). Inoltre la derivazione "corretta"
dell'equazione di stato di un filtro con parametro variabile fa comparire
termini in `wn_dot`, che qui **non ci sono**: le righe 31-33 sono la forma
"frozen-form" (congela la struttura, sostituisci il parametro), che e' la
prassi ingegneristica ma non un'identita' matematica. Vale finche' `wn_dot` e'
piccola rispetto a `wn^2` -- ipotesi di variazione lenta, di nuovo assunta e
non verificata dal codice. Lo stesso identico caveat vale per la realizzazione
a blocchi elementari del notch in Simulink
(`build_hm3_full_ascent_flex.m`, righe 125-141), che ricalca queste righe.

---

## TVC e ritardo: unico blocco LTI (righe 34-35)

```matlab
delta = M.Ct*xt + M.Dt*v;
xtd   = M.At*xt + M.Bt*v;
```

- Righe 34-35: la catena di attuazione e' una **quadrupla costante**
  `(At, Bt, Ct, Dt)` (uscita alla 34, derivata di stato alla 35), ottenuta in
  `init_simulink_lpv.m` (righe 112-115) da
  `ssdata(ss(build_tvc(p0, 3)))`, con `p0 = load_hw3_params()`, cioe'
  **congelata a t = 72 s**. E' legittimo: `wTVC = 70 rad/s`, `zTVC = 0.7`,
  `tau = 0.020 s` sono costanti di tabella, non funzioni del tempo
  (`load_hw3_params.m`, righe 55-57). Quindi il TVC e' l'**unico** pezzo
  dell'anello che non e' LPV -- e non lo e' per motivi fisici, non per
  approssimazione.
- Il ritardo di trasporto di 20 ms e' un **Pade di ordine 3**
  (`build_tvc.m`, righe 18-19): un'approssimazione razionale, che introduce
  zeri in semipiano destro. E' fedele a fase bassa e sbaglia sempre di piu'
  ad alta frequenza. Poiche' omega_BM arriva fino a ~32 rad/s lungo l'ascesa
  (contro i 18.9 rad/s a cui e' stato scelto l'ordine del Pade), **la qualita'
  del ritardo simulato degrada proprio dove serve**: l'analisi di detuning e'
  quantitativamente ottimistica o pessimistica a seconda del segno dell'errore
  di fase, e il codice non lo quantifica.
- **Ordine dei blocchi**: `u_pd -> notch -> TVC -> delta`. Il notch filtra il
  comando *prima* dell'attuatore, e questo replica l'ordine di `assemble_loop`
  in HM3 (`Wact = Wtvc * notch`, vedi `main_flex.m` righe 27-28). Coerente.
- Non ci sono **saturazioni** ne' limiti di rate sul TVC: `delta` puo' assumere
  qualunque valore. Nel caso "notch fisso" la simulazione diverge, e i valori
  di `delta` a fine corsa sono privi di senso fisico (un attuatore vero
  saturerebbe). La divergenza rimane il risultato corretto sul piano della
  stabilita' lineare, ma le ampiezze non vanno lette come numeri fisici.

---

## Il plant flessibile (righe 38-40)

```matlab
zdd   = a1*zdot + (a1*V + a4)*theta + a3*delta - a1*V*aw;
thdd  = (A6/V)*zdot + A6*theta + K1*delta - A6*aw;
etadd = -w^2*eta - 2*M.zBM*w*etadot + aqk*delta;
```

- Righe 38-39: **identiche** alle prime quattro righe della A di
  `build_plant_full.m` (righe 20-28), con i coefficienti sostituiti dai
  loro interpolanti. Confermano la lettura del segno del vento: raccogliendo,
  `thdd = A6*(theta + zdot/V - alpha_w) + K1*delta`, quindi l'incidenza
  aerodinamica **effettiva** e' `theta + zdot/V - alpha_w` (conta la velocita'
  relativa all'aria: la componente laterale relativa e' `zdot - v_w`).
- **Il bending non retroagisce sulla dinamica rigida**: nelle righe 38-39 non
  compaiono `eta` ne' `etadot`. E' la stessa struttura triangolare a blocchi
  di `build_plant_full.m` (le prime 4 righe della A hanno zeri nelle colonne
  5-6). Fisicamente e' un'approssimazione del modello di corso: la
  deformazione elastica non perturba il moto del corpo rigido, ma sporca le
  misure ed e' eccitata dal TVC.
- Riga 40: oscillatore del primo modo. `w = omega(t)` **varia** (e' questo che
  detuna il notch), `M.zBM = 0.005` e' **costante** -- viene da
  `S.notch.zBM = p0.zBM` (`init_simulink_lpv.m`, riga 116), cioe' dal
  letterale di tabella (`load_hw3_params.m`, riga 49), non dal dataset. Il
  dataset fornisce `omega` ma il codice non cerca uno `zeta` tempo-variante.
- `aqk = -phi_tvc*Tc` [interpolato dal dataset, riga 20] e' il forzamento del
  bending da parte del TVC: la spinta deviata agisce sulla forma modale nel
  punto di applicazione. Segno negativo, quindi `delta > 0` eccita `eta < 0`.
- **Il vento non eccita il bending**: nella riga 40 non c'e' `aw`. Coerente
  con `build_plant_full.m` (riga 28: `Bw` ha 0 nella riga di `etadot`), ma
  fisicamente e' una semplificazione forte -- una raffica reale carica la
  struttura direttamente. Il modo elastico qui puo' essere eccitato **solo**
  attraverso il comando di TVC.
- Nota sull'omega variabile: come per il notch, `etadd` con `w = w(t)` e' la
  forma "frozen-form" dell'oscillatore. L'equazione modale rigorosa di una
  struttura a massa/rigidezza variabili contiene termini aggiuntivi (derivate
  della massa generalizzata e della forma modale). Il codice li ignora --
  scelta standard, ma va detta.

---

## Assemblaggio (riga 42)

```matlab
dx = [zdot; zdd; thetadot; thdd; etadot; etadd;
      xn1d; xn2d; xtd];
```

- Riga 42: concatenazione a colonna nell'**esatto ordine** dello split di
  riga 15. Se qualcuno cambiasse l'ordine di `x0` in `main_flex.m` (riga 38)
  senza cambiare qui, il sistema integrato sarebbe silenziosamente un altro:
  il contratto e' posizionale, non nominale. E' il costo dell'avere una RHS
  veloce senza validazione.
- `xtd` e' un vettore nt x 1: la concatenazione funziona solo perche' `At*xt`
  restituisce una colonna. Nessun `(:)` difensivo -- corretto e voluto in un
  hot loop.

---

## Rigido vs flessibile: cosa cambia davvero

| | `ode_lpv_ascent.m` | `ode_lpv_flex.m` |
|---|---|---|
| stati | 4 | 13 (6 + 2 notch + 5 TVC) |
| coefficienti | combinati (`fc1..fc7`) | grezzi, combinati in linea |
| guardia su V | `Vsafe = max(V,1)` a monte | nessuna: divide per `fV(t)` |
| attuatore | ideale (`delta = u_pd`) | TVC 2 ordine + Pade(3), LTI |
| filtro | nessuno | notch a `wn(t)` (o fisso) |
| misura | stato vero | contaminata dall'INS (Eq. 2) |
| bending | assente | modo 1, `omega(t)`, `zBM` costante |
| guadagni PD | congelati o schedulati | congelati (schedule supportata, mai usata) |

La **variabile aggiuntiva chiave** e' `omega_BM(t)`: la frequenza del primo
modo cambia lungo l'ascesa perche' cambiano massa e rigidezza generalizzate
del veicolo mentre brucia propellente. Nel modello rigido questo non ha alcun
effetto (il modo non c'e'); nel modello flessibile e' la cosa che rompe il
progetto congelato di HM3, perche' il notch e' un filtro **puntato** su una
frequenza e quella frequenza si sposta.

---

## La fallacia del frozen-time, vista da qui

Il ragionamento "ho verificato i margini su ogni plant congelato, dunque il
sistema tempo-variante e' stabile" **non e' valido**. Per `xdot = A(t)*x` la
stabilita' dipende dalla matrice di transizione Phi(t,t0), non dagli autovalori
istantanei: esiste il classico controesempio (Khalil, *Nonlinear Systems*) di
una A(t) i cui autovalori congelati stanno **sempre** in
`(-1 +/- i*sqrt(7))/4`, parte reale -0.25 per ogni t, e la cui soluzione
diverge come `exp(+0.5*t)`. L'implicazione diventa vera solo sotto **ipotesi di
variazione lenta** (`||Adot||` sufficientemente piccola), che qui e' assunta e
mai quantificata.

Va sottolineato che **anche la direzione opposta e' invalida**: dal fatto che
gli LTI congelati siano instabili non segue che l'LTV lo sia. Questo tocca
direttamente lo studio T008, perche' `main_flex.m` (righe 24-33) conclude
"il notch fisso va instabile da t ~ 75 s" leggendo `isstable(Tf)` **sui plant
congelati**. Preso da solo sarebbe un argomento formalmente insufficiente.

Come il codice si difende: non si ferma li'. Alle righe 40-41 di `main_flex.m`
**propaga davvero l'LTV** con `ode45` (`RelTol = 1e-8`, `AbsTol = 1e-10`) sulla
stessa realizzazione di vento, e osserva che il coordinato di bending `eta`
diverge davvero nel caso fisso e resta dell'ordine di 1e-3 nel caso variabile.
I due argomenti -- sweep congelato + integrazione LTV -- si sostengono a
vicenda: nessuno dei due, da solo, sarebbe una prova. E nemmeno insieme lo
sono in senso stretto: la simulazione e' **una** traiettoria, con **un** vento,
da **una** condizione iniziale (`x0 = zeros(6 + 2 + nt, 1)`, cioe' `zeros(13,1)`,
riga 38). Una prova vera
richiederebbe di propagare Phi(t,t0) a ingresso nullo, oppure una funzione di
Lyapunov tempo-variante / una LMI parameter-dependent. Il repo non le fa e non
le rivendica.

> **Possibile domanda d'esame** -- il notch variabile e' un controllore LPV
> legittimo o un trucco?
> *Risposta:* E' gain scheduling nel senso classico: si progetta una famiglia
> di filtri sui plant congelati e si interpola il parametro. E' legittimo sotto
> variazione lenta, ma non porta con se' alcuna garanzia di stabilita' LTV:
> serve la verifica a posteriori (l'integrazione LTV), oppure una sintesi LPV
> vera (LMI con funzione di Lyapunov dipendente dal parametro, che vincola
> anche `wn_dot`). Il codice fa la prima strada e lo dichiara.

---

## Limiti e punti fragili (riepilogo onesto)

- **Nessuna guardia su V**: la riga 39 divide per `fV(t)`; V(0) = 0. La RHS non
  puo' partire da t = 0. Il ramo rigido ha `Vsafe`, questo no.
- **Solo notch a fase minima** realizzabile (riga 31 hard-codea il segno `+`).
- **`zBM` costante** (0.005, letterale di tabella), mentre `omega` e' tabulato.
- **Il vento non forza il bending** e **il bending non forza il corpo rigido**:
  entrambi ereditati da `build_plant_full.m`, entrambi semplificazioni.
- **Nessuna saturazione** su `delta`: nel caso divergente le ampiezze non hanno
  significato fisico.
- **Pade(3) tarato sul punto a 18.9 rad/s** ma usato fino a ~32 rad/s.
- **Forma "frozen-form"** per notch e bending: mancano i termini in `wn_dot` e
  `omega_dot`; valida solo sotto variazione lenta.
- **Contratto posizionale sullo stato** (riga 15): fragile a modifiche di `x0`.
- `M.sched` esiste ma nessuno script lo attiva su questa RHS: e' codice
  supportato ma non esercitato (quindi non testato).

---

## Possibili domande d'esame

**D: Che cosa aggiunge questa RHS rispetto a `ode_lpv_ascent`, e perche' e'
proprio quello a rompere il progetto congelato?**
R: Aggiunge il primo modo di bending (riga 40), la contaminazione INS delle
misure (righe 24-25) e la catena TVC + ritardo con il notch (righe 30-35). Il
punto critico e' che `omega_BM` **varia** lungo l'ascesa, perche' massa e
rigidezza generalizzate cambiano mentre brucia il primo stadio. Il notch
profondo di HM3 e' centrato su omega_BM(72) = 18.9 rad/s: quando la frequenza
si sposta, la risonanza esce dalla buca, il guadagno d'anello alla frequenza
di bending risale sopra 0 dB e l'anello -- chiuso su una misura contaminata da
un modo con smorzamento 0.005 -- diverge.

**D: Perche' l'instabilita' passa dalle misure e non dalla dinamica?**
R: Perche' nel modello del corso il bending **non retroagisce** sul corpo
rigido: nelle righe 38-39 non compaiono `eta` ne' `etadot`. Il modo entra solo
in `theta_m = theta + sigma*eta` e `z_m = z - phi*eta` (Eq. 2). Il PD
retroaziona quella misura sporca, il comando `delta` che ne esce eccita il modo
(`aqk*delta`, riga 40), e il ciclo si chiude. Senza contaminazione INS il
bending sarebbe uno stato non osservato e innocuo.

**D: La stabilita' di ogni plant congelato garantisce la stabilita' del sistema
tempo-variante?**
R: No. Per un LTV la stabilita' dipende dalla matrice di transizione, non dagli
autovalori istantanei; esistono controesempi standard con autovalori congelati
costanti a parte reale negativa e soluzioni divergenti (Khalil). Serve
l'ipotesi di variazione lenta. E vale anche il contrario: da "tutti gli LTI
congelati sono instabili" non segue che l'LTV lo sia -- il che indebolisce, se
preso da solo, lo sweep `isstable` di `main_flex.m` (righe 24-33). Per questo
il codice propaga davvero l'LTV con `ode45` (righe 40-41) e usa i due argomenti
insieme.

**D: Come e' realizzato il notch e perche' `xn(1)` non compare nell'uscita?**
R: E' la forma canonica di controllabilita' di
`N(s) = (s^2 + 2 zN wn s + wn^2)/(s^2 + 2 zD wn s + wn^2)`. Il coefficiente
d'uscita su `xn1` vale `b_0 - alpha_0*d = wn^2 - wn^2 = 0`, quindi sparisce; su
`xn2` vale `b_1 - alpha_1*d = 2*(zN - zD)*wn` (riga 31) e il feedthrough vale 1
(guadagno unitario fuori dalla buca). Con `zN = 0.002` e `zD = 0.7` il notch e'
profondissimo e strettissimo -- che e' anche il motivo per cui si detuna cosi'
male.

**D: Sostituire `wn(t)` in una realizzazione congelata e' matematicamente
corretto?**
R: Non esattamente. Due realizzazioni della stessa funzione di trasferimento
sono equivalenti solo a parametri costanti; con parametri variabili
realizzazioni diverse danno sistemi LTV diversi, e la derivazione rigorosa fa
comparire termini in `wn_dot` che qui non ci sono. E' la prassi standard
("frozen-form"), valida sotto variazione lenta, ma e' un'approssimazione, non
un'identita'. Lo stesso vale per l'oscillatore di bending con `omega(t)`.

**D: Perche' il TVC e' l'unico blocco LTI in un modello LPV?**
R: Perche' i suoi parametri sono davvero costanti: `wTVC = 70 rad/s`,
`zTVC = 0.7`, `tau = 20 ms` sono proprieta' dell'attuatore, non del volo. La
quadrupla `(At, Bt, Ct, Dt)` viene estratta una volta da `build_tvc(p0, 3)` e
riusata. L'unica approssimazione e' il Pade di ordine 3 per il ritardo, tarato
implicitamente attorno alla frequenza di bending del punto di progetto
(18.9 rad/s) ma usato fino a ~32 rad/s.
