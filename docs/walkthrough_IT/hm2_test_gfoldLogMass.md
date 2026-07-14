# HM2_powered_descent/tests/gfoldLogMassTest.m

## Ruolo del file nel progetto

E' la suite `matlab.unittest` che copre i **due kernel numerici della trascrizione
GFOLD log-mass**: `ode_descent_uacc.m` (il right-hand side nonlineare con
l'accelerazione tenuta costante) e `lti_zoh.m` (la discretizzazione ZOH esatta via
matrice esponenziale). Sono i due file che il prototipo `proto_gfold_logmass.m` e la
"variant (d)" dentro `main_task2.m` **condividono**: testarli una volta copre
entrambi i chiamanti.

La suite ha 7 test: 4 sul RHS e 3 sulla discretizzazione. Presi insieme non
verificano una cosa qualunque -- verificano **la tesi centrale della trascrizione**,
cioe' che il cambio di variabile `z = ln(m)`, `u = T/m` rende la dinamica
esattamente LTI e che quindi la mappa discreta e' esatta e non approssimata. Ogni
test isola una delle tre proprieta' che servono a quella tesi:

1. le righe di velocita' **non contengono la massa** (`testAccelerationIsDirect`);
2. la riga di massa e' **lineare in `z`** -- equivalentemente `m' e' proporzionale a
   `m` (`testMassFlowScalesWithMass`);
3. con il cono lossless attivo (`sigma = |u|`), la riga discreta di `z` riproduce
   **esattamente** lo svuotamento nonlineare della massa (`testMassRowConsistency`).

Le proprieta' 1 e 2 sono precisamente le due condizioni sotto cui `A` e `B` di
`lti_zoh` sono matrici *costanti*. La 3 e' il ponte fra il modello convesso e la
realta' nonlineare: e' quella che giustifica il replay `ode45` del prototipo.

Onesta' preliminare: la suite **non tocca l'ottimizzatore**. `solve_gfold_socp`,
`solve_gfold_scvx`, il vincolo di spinta linearizzato, il glide-slope e la
losslessness all'ottimo non sono testati qui (la losslessness e' solo *stampata* a
runtime dal prototipo, `cone_gap`). Vengono testati i mattoni, non l'edificio.

---

## `properties (Constant)` (righe 6-9)

```matlab
Vc = 0.0777;       % V_ref/c
dt = 0.0444;       % tf_nd/(N-1) with N = 50
```

- Righe 7-8: le due costanti del run nominale, **hard-coded e arrotondate a 3 cifre**.
  Non vengono ricalcolate dai dati di Tabella 1: se domani cambiassero `Isp`, `y0` o
  `tf`, il test continuerebbe a girare con i vecchi numeri. Questo **non e' un bug**,
  perche' tutte le identita' verificate sono algebriche e valgono per *qualunque*
  `Vc` e `dt` positivi: i valori servono solo a stare in un regime numerico
  realistico. Ma va saputo: la suite non e' una regression sui dati del problema.
- Per la cronaca, con i dati del prototipo (`y0 = 3000`, `g = 9.81`, `Isp = 225`,
  `tf = 38`, `N = 50`) i valori esatti sono `Vc = V_rif/(Isp*g0) ~ 0.0777` e
  `dt = tf/sqrt(y0/g)/49 ~ 0.0444`: coerenti.

---

## `addHm2ToPath` (righe 11-16)

```matlab
hm2 = fileparts(fileparts(mfilename('fullpath')));
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(hm2));
```

- Riga 13: `mfilename('fullpath')` da' il path di questo file (che sta in
  `tests/`); il doppio `fileparts` risale a `HM2_powered_descent/`.
- Riga 14: `PathFixture` aggiunge quella cartella al path **per la durata della
  classe** e la rimuove alla fine. E' cio' che permette di chiamare
  `ode_descent_uacc` e `lti_zoh` (che vivono nella cartella padre) senza sporcare il
  path globale dell'utente. Fixture di classe (`TestClassSetup`), quindi si paga una
  volta sola.

---

## `testUaccDerivativeDefinition` (righe 20-27)

```matlab
x    = [1; 2; 0.3; -0.4; 0.8];
uacc = [0.6; 0.8];                 % ||u|| = 1
dx   = ode_descent_uacc(x, uacc, testCase.Vc);
expected = [0.3; -0.4; 0.6; 0.8 - 1; -testCase.Vc * 0.8 * 1];
```

- Riga 24: chiama il RHS con stato `[x; y; vx; vy; m] = [1; 2; 0.3; -0.4; 0.8]` e
  accelerazione `u = [0.6; 0.8]`, scelta apposta con **norma unitaria** (terna
  pitagorica 3-4-5 scalata) cosi' che il termine `|u|` nel valore atteso sia
  esattamente 1 e non introduca errore di arrotondamento.
- Riga 25: il valore atteso e' **scritto a mano**, non ricalcolato con la stessa
  formula del codice sotto test -- condizione necessaria perche' il test abbia senso.
  Componente per componente:
  - `dx(1) = vx = 0.3`, `dx(2) = vy = -0.4` -- cinematica, banale;
  - `dx(3) = ux = 0.6` -- **niente `/m`**;
  - `dx(4) = uy - 1 = -0.2` -- gravita' non-dim pari a 1, diretta verso il basso.
    E' qui che si fissa la convenzione di segno: `y` cresce verso l'alto;
  - `dx(5) = -Vc * m * |u| = -Vc * 0.8 * 1` -- la riga di massa **in variabile `m`**
    (non `z`): `m' = -Vc*|T| = -Vc*m*|u|`.
- Riga 26: `AbsTol 1e-15`, cioe' a livello di epsilon macchina. Legittimo: sono somme
  e prodotti esatti in floating point, non c'e' integrazione.

> **Possibile domanda d'esame** -- Perche' il RHS prende `uacc` e non la spinta `T`?
> *Risposta:* Perche' e' la convenzione ZOH nativa del GFOLD: quello che si tiene
> costante sull'intervallo e' **l'accelerazione** `u = T/m`, non la spinta. Con `u`
> costante, `T(t) = m(t)*u` decresce durante l'intervallo insieme alla massa. E'
> esattamente questa scelta che rende `v' = u + g` indipendente dalla massa e quindi
> il sistema LTI. `ode_descent.m` (le altre trascrizioni) tiene costante `T` e quindi
> ha `v' = T/m(t)`, che dentro l'intervallo varia.

---

## `testAccelerationIsDirect` (righe 29-37)

```matlab
dx1 = ode_descent_uacc([0;1;0;0;0.9], uacc, testCase.Vc);
dx2 = ode_descent_uacc([0;1;0;0;0.3], uacc, testCase.Vc);
testCase.verifyEqual(dx1(3:4), [0.5; 0.2], 'AbsTol', 1e-15);
testCase.verifyEqual(dx2(3:4), [0.5; 0.2], 'AbsTol', 1e-15);
```

- Righe 33-36: stesso `uacc`, **due masse molto diverse** (0.9 e 0.3, cioe' 1800 kg e
  600 kg), stesse derivate di velocita'. E' il test **strutturale** della
  trascrizione: dimostra che le righe 3-4 del RHS non dipendono da `x(5) = m`.
- Perche' e' *la* proprieta' giusta: nel sistema `xi' = A*xi + B*w + c`, la colonna
  di `B` che moltiplica `u` e' **costante**. Se `v'` dipendesse dalla massa, `B`
  sarebbe `B(xi)` e il sistema non sarebbe LTI: nessun `expm` una-tantum, nessuna
  discretizzazione esatta, e il SOCP diventerebbe un NLP. Questo test e' quindi il
  guardiano della premessa su cui poggia tutta la variante (d).
- Il commento (righe 30-31) lo dice esplicitamente: "unlike ode_descent.m where the
  acceleration is T/m".

---

## `testMassFlowScalesWithMass` (righe 39-46)

```matlab
uacc = [0; 1];                     % ||u|| = 1
dxA  = ode_descent_uacc([0;1;0;0;1.0], uacc, vc);   % m = 1.0
dxB  = ode_descent_uacc([0;1;0;0;0.5], uacc, vc);   % m = 0.5
% attesi: -Vc*1.0  e  -Vc*0.5
```

- Righe 42-45: la portata di massa **scala linearmente con la massa**:
  `m' = -Vc*m*|u|`. Con `|u| = 1`, dimezzando `m` si dimezza `m'`.
- **Perche' questa e' la proprieta' giusta**: e' l'enunciato "in variabile `m`" del
  fatto che la riga di massa e' **lineare in `z = ln(m)`**. Infatti
  `m'/m = -Vc*|u|` costante equivale a `d(ln m)/dt = -Vc*|u|`, cioe' `z' = -Vc*|u|`:
  un'equazione **senza stato a destra**, che dopo l'introduzione dello slack diventa
  `z' = -Vc*sigma`, la quinta riga di `B` in `lti_zoh` (`B(5,3) = -Vc`). Il commento
  alla riga 40 lo scrive nero su bianco.
- Insieme, `testAccelerationIsDirect` e `testMassFlowScalesWithMass` **sono** la
  verifica che il cambio di variabile linearizza esattamente la dinamica: uno toglie
  la massa dal denominatore delle velocita', l'altro la toglie dalla riga di massa.

> **Possibile domanda d'esame** -- Perche' il logaritmo? Non bastava usare `m` come
> stato?
> *Risposta:* No. Con `m` come stato la riga di massa e' `m' = -Vc*m*|u|`, che e'
> **bilineare** in (stato, controllo): non e' LTI. Passando a `z = ln(m)` la stessa
> equazione diventa `z' = -Vc*|u|`, che non contiene lo stato: la riga di massa e'
> lineare e completamente disaccoppiata. Il logaritmo e' esattamente la trasformazione
> che uccide quella bilinearita' -- e in piu' rende lineare anche il costo, perche'
> massimizzare `z_N` equivale a massimizzare `m_N = exp(z_N)`.

---

## `testBallisticCoast` (righe 48-52)

```matlab
dx = ode_descent_uacc([0.5;1;0.2;-0.1;0.9], [0;0], testCase.Vc);
testCase.verifyEqual(dx, [0.2; -0.1; 0; -1; 0], 'AbsTol', 1e-15);
```

- Riga 50-51: caso limite `u = 0`. Attesi: caduta libera (`v' = [0; -1]`, cioe' solo
  gravita') e **portata di massa nulla** (`dx(5) = 0`).
- Verifica due cose non banali: (i) il termine costante `c = [0;0;0;-1;0]` di
  `lti_zoh` e' coerente con il RHS nonlineare; (ii) `|u| = 0` non produce NaN o
  divisioni per zero nel calcolo di `umag` (`sqrt(0+0) = 0`, poi moltiplicato: nessuna
  patologia). Il coast e' l'arco centrale della soluzione ottima (profilo
  max-coast-max secondo il README), quindi non e' un caso accademico: e' meta' della
  traiettoria.

---

## `testZohClosedForm` (righe 55-69)

```matlab
[Abar, Bbar, cbar] = lti_zoh(h, vc);
A_exp = eye(5);  A_exp(1,3) = h;  A_exp(2,4) = h;
B_exp = [h^2/2, 0,     0;
         0,     h^2/2, 0;
         h,     0,     0;
         0,     h,     0;
         0,     0,    -vc*h];
c_exp = [0; -h^2/2; 0; -h; 0];
```

- Righe 59-65: i valori attesi sono le **formule del moto uniformemente accelerato**,
  scritte a mano. Il test confronta l'output di `expm` (van Loan) con la forma chiusa.
- **Da dove viene la forma chiusa.** La matrice `A` di `lti_zoh` ha solo `A(1,3) = 1`
  e `A(2,4) = 1`, quindi `A^2 = 0` (**nilpotente**: la riga 3 e la riga 4 di `A` sono
  nulle). La serie esponenziale si tronca dopo due termini:

      e^{A*h} = I + A*h

  e l'integrale di convoluzione:

      integrale_0^h e^{A*s} ds = h*I + (h^2/2)*A

  da cui

      Abar = I + A*h
      Bbar = h*B + (h^2/2)*A*B
      cbar = h*c + (h^2/2)*A*c

  Facendo i prodotti: `A*B` ha `1` in posizione (1,1) e (2,2) -- cioe' l'accelerazione
  che si integra due volte in posizione -- e `A*c = [0; -1; 0; 0; 0]`. Si ottiene
  esattamente `B_exp` e `c_exp`: `x_{k+1} = x_k + vx_k*h + ux*h^2/2`,
  `y_{k+1} = y_k + vy_k*h + (uy - 1)*h^2/2`, `vx_{k+1} = vx_k + ux*h`,
  `vy_{k+1} = vy_k + (uy - 1)*h`, `z_{k+1} = z_k - Vc*sigma*h`. Il `-h^2/2` e il `-h`
  in `c_exp` sono il contributo della gravita' unitaria.
- **Perche' e' la proprieta' giusta da testare**: certifica che la mappa discreta non
  e' un'approssimazione ma **la soluzione esatta** dell'intervallo. Se `A` non fosse
  nilpotente, `Abar` conterrebbe tutta la serie e non ci sarebbe forma chiusa con cui
  confrontarsi. E' anche un test dell'**implementazione** del trucco di van Loan
  (indici delle sotto-matrici, `E(1:5,1:5)`, `E(1:5,6:8)`, `E(1:5,9)`, e il blocco
  `zeros(4,9)` che corrisponde alle 3 righe del controllo piu' 1 riga della costante):
  un errore di slicing verrebbe preso qui.
- Riga 66-68: `AbsTol 1e-12`. Piu' lasca di `1e-15` perche' `expm` fa scaling-and-
  squaring e non e' esatta a bit.

> **Possibile domanda d'esame** -- Perche' `lti_zoh` costruisce una matrice 9x9?
> *Risposta:* Perche' il sistema e' affine, non lineare puro: c'e' il termine costante
> di gravita' `c`. Il trucco di van Loan aumenta lo stato con il controllo (3
> componenti) e con la costante (1 componente), 5 + 3 + 1 = 9, e sfrutta l'identita'
> `expm([A B c; 0 0 0]*h) = [e^{A h}, integrale*B, integrale*c; 0, I]`: le
> sotto-matrici in alto a destra **sono** `Bbar` e `cbar`. Cosi' `Abar`, `Bbar`, `cbar`
> escono da **un solo** `expm`, senza integrare nulla numericamente.

---

## `testZohMatchesOde45` (righe 71-83)

```matlab
xi0 = [0.3; 1; -0.2; -0.6; 0];        % [x;y;vx;vy;z]
w   = [0.4; 1.1; 1.3];                % [ux;uy;sigma]
rhs = @(~, xi) [xi(3); xi(4); w(1); w(2)-1; -vc*w(3)];
[~, XI] = ode45(rhs, [0, h], xi0, opts);
testCase.verifyEqual(Abar*xi0 + Bbar*w + cbar, XI(end,:).', 'AbsTol', 1e-9);
```

- Righe 76-82: un passo della mappa discreta contro un'integrazione `ode45` della
  **stessa** dinamica LTI con `w` costante sull'intervallo. Tolleranze
  `RelTol 1e-12 / AbsTol 1e-14` sull'integratore, confronto a `1e-9` -- che e' il floor
  realistico di `ode45`, non un margine generoso.
- Nota che `w(3) = sigma = 1.3` e' **volutamente maggiore** di
  `|u| = sqrt(0.4^2 + 1.1^2) ~ 1.17`: qui il cono e' **inattivo**. Il test verifica la
  mappa discreta del **modello rilassato**, non della fisica. E' corretto cosi': deve
  certificare che `Abar/Bbar/cbar` riproducono la dinamica *che il SOCP crede*, con
  `sigma` come ingresso indipendente.
- **Debolezza dichiarata**: la RHS della riga 78 e' **riscritta a mano dentro il
  test**, non presa da un file sorgente (nessun file della repo espone la dinamica LTI
  in coordinate `xi`; `ode_descent_uacc` lavora in `m`, non in `z`). Se ci fosse un
  errore di segno *identico* in `lti_zoh.m` e in questa lambda, il test passerebbe lo
  stesso. La copertura e' salvata dal test precedente (`testZohClosedForm`), che
  confronta con formule chiuse **indipendenti** dal codice: i due test insieme sono
  ridondanti nel modo giusto.

---

## `testMassRowConsistency` (righe 85-97)

```matlab
uacc = [0.5; 1.0];   umag = norm(uacc);
[~, X] = ode45(@(~,x) ode_descent_uacc(x, uacc, vc), [0 h], [0;1;0;0;1], opts);
z_nl  = log(X(end,5));               % log-massa nonlineare
z_lti = -vc * umag * h;              % predizione LTI con sigma = ||u||
testCase.verifyEqual(z_nl, z_lti, 'AbsTol', 1e-9);
```

E' **il test piu' importante della suite**, e va capito bene perche'.

- Righe 92-94: si integra la dinamica **nonlineare in massa** (`m' = -Vc*m*|u|`,
  quindi `m(h) = exp(-Vc*|u|*h)` partendo da `m(0) = 1`), poi si prende il logaritmo
  della massa finale.
- Riga 95: si calcola la predizione del **modello discreto LTI** per la variazione di
  `z`, ponendo `sigma = |u|`.
- Riga 96: devono coincidere a `1e-9` (floor `ode45`).

**Cosa certifica davvero.** Certifica il **ponte** fra il modello convesso e la
dinamica vera: *se e solo se il cono e' attivo* (`sigma = |u|`), la quinta riga della
mappa LTI riproduce **esattamente** lo svuotamento nonlineare della massa. Il partire
da `m(0) = 1` (quindi `z(0) = 0`) rende il confronto diretto:
`z_nl = ln(m(h)) - 0 = z(h) - z(0)`.

Ed e' esattamente la proprieta' giusta, perche' e' l'unica cerniera fragile
dell'intera trascrizione. Ragionamento: le righe di posizione e velocita' del modello
discreto sono esatte **sempre** (non contengono la massa, e la mappa e' la soluzione
in forma chiusa). Quindi **l'unica** sorgente possibile di discrepanza fra modello
discreto e realta' nonlineare e' la riga di massa, e la riga di massa e' esatta
**se e solo se `sigma = |u|`**. Da qui due conseguenze pratiche, entrambe visibili nel
prototipo:

- il replay `ode45` (`fwd_integrate_uacc`, `proto_gfold_logmass.m` righe 241-254)
  torna al floor dell'integratore **solo** perche' il rilassamento e' lossless;
- il "drift" di massa modello-vs-replay (`dmf`, riga 64 del prototipo) e il `cone_gap`
  (riga 71) misurano **la stessa cosa**.

Se il rilassamento non fosse lossless, il modello brucerebbe piu' propellente del
replay (`sigma > |u|`), il costo sarebbe pessimistico e questo test sarebbe il
sintomo. Il test **non** dimostra che all'ottimo il cono e' attivo (quello e' un
risultato di ottimalita', via PMP): dimostra che **se** lo e', la trascrizione e'
esatta.

> **Possibile domanda d'esame** -- Questo test verifica la losslessness del
> rilassamento?
> *Risposta:* No, e la distinzione e' importante. Verifica la **conseguenza** della
> losslessness, cioe' che *ponendo* `sigma = |u|` la mappa discreta della massa e'
> esatta. Che `sigma = |u|` valga **all'ottimo** e' un risultato di teoria del
> controllo ottimo (PMP: `u` compare in Hamiltoniana solo via `lambda_v . u`, quindi
> il minimo su `{|u| <= sigma}` sta sul bordo), e nel codice viene solo **verificato a
> posteriori** dal `cone_gap` stampato dal prototipo. Nessun test unitario controlla
> questa proprieta'.

---

## Cosa NON e' coperto (onesta')

1. **Nessun test sull'ottimizzatore.** `solve_gfold_socp` e `solve_gfold_scvx` sono
   funzioni *locali* dello script `proto_gfold_logmass.m` (e duplicate in
   `main_task2.m`): non sono raggiungibili dalla suite. Il vincolo di spinta
   linearizzato (`sigma <= Tmax*e^{-z_ref}*(1 - (z - z_ref))`), il glide-slope, la
   trust region e il ratio test **non hanno test**.
2. **Nessun test della losslessness all'ottimo** (`cone_gap`): e' solo una stampa a
   runtime.
3. **Nessun test di conservativita' della tangente**, cioe' che il bound linearizzato
   implichi davvero `|T| <= Tmax` (proprieta' garantita dalla convessita' di
   `exp(-z)`, ma mai verificata numericamente).
4. **Le costanti `Vc`/`dt` sono arrotondate e hard-coded**: la suite non e' una
   regression sui dati di Tabella 1.
5. **La dinamica LTI in coordinate `xi` e' duplicata nel test** (riga 78) invece di
   essere importata da un sorgente: rischio teorico di errori correlati, mitigato da
   `testZohClosedForm`.
6. **Nessun performance test** (`matlab.perftest`) per questi kernel, a differenza di
   quanto prevedono le convenzioni della repo per altre feature.

---

## Possibili domande d'esame

**D: I 7 test cosa dimostrano, messi insieme?**
R: Dimostrano la tesi su cui poggia l'intera variante GFOLD: che il cambio di
variabile rende la dinamica **esattamente LTI** e quindi la discretizzazione ZOH
**esatta**. `testAccelerationIsDirect` dimostra che le righe di velocita' non
contengono la massa (quindi `B` e' costante); `testMassFlowScalesWithMass` dimostra
che `m'/m` non dipende dallo stato (quindi la riga di `z` e' lineare);
`testZohClosedForm` dimostra che `expm` restituisce le formule del moto uniformemente
accelerato (nessun errore di troncamento, perche' `A^2 = 0`); `testZohMatchesOde45`
chiude il cerchio contro un integratore; `testMassRowConsistency` dimostra che con il
cono attivo il modello convesso e la dinamica nonlineare coincidono.

**D: Perche' `testMassRowConsistency` e' il test critico?**
R: Perche' isola **l'unico punto in cui il modello discreto puo' mentire**. Le righe
di posizione/velocita' sono esatte per costruzione (soluzione in forma chiusa di un
moto uniformemente accelerato, con `u` costante per ipotesi ZOH). L'unica
approssimazione concettuale dell'intera trascrizione e' aver sostituito `|u|` con lo
slack `sigma` nell'equazione della massa. Questo test verifica che, quando lo slack e'
saturo (`sigma = |u|`, cioe' quando il rilassamento e' lossless), la sostituzione non
costa nulla: `z(h) - z(0) = -Vc*|u|*h` a precisione di integratore.

**D: Perche' `testZohMatchesOde45` usa `sigma = 1.3 > |u| ~ 1.17`?**
R: Perche' quel test riguarda la **mappa discreta del modello rilassato**, dove
`sigma` e' un ingresso indipendente da `u`. Deve funzionare per qualunque `w`
ammissibile, cono attivo o no. Il caso "cono attivo" e' invece il soggetto specifico
di `testMassRowConsistency`. Sono due domande diverse: "la mia matrice discreta e'
giusta?" e "il mio modello coincide con la fisica quando il rilassamento e' stretto?".

**D: Perche' le tolleranze sono cosi' diverse (1e-15 contro 1e-9)?**
R: I test su `ode_descent_uacc` (righe 20-52) confrontano una manciata di somme e
prodotti in floating point con il loro valore atteso: l'errore ammissibile e' l'epsilon
macchina, `1e-15`. `testZohClosedForm` usa `1e-12` perche' `expm` implementa
scaling-and-squaring e introduce un errore di qualche ulp. I due test che integrano
(`testZohMatchesOde45`, `testMassRowConsistency`) sono limitati dal controllo di passo
di `ode45`: anche con `RelTol 1e-12 / AbsTol 1e-14` il floor realistico e' `1e-9`.
Stringere di piu' produrrebbe test intermittenti.

**D: Se volessi testare davvero la losslessness, come faresti?**
R: Servirebbe estrarre `solve_gfold_socp` in un file autonomo (oggi e' una funzione
locale dello script) e scrivere un test che: (a) risolve il SOCP su un caso piccolo
(N basso), (b) verifica `max_k |sigma_k - |u_k|| < 1e-8`, (c) verifica che la spinta
ricostruita `T = m*u` rispetti `|T| <= Tmax` a tutti i nodi, (d) confronta il costo
ottimo con quello di una delle altre trascrizioni entro qualche kg. Sarebbe un test di
integrazione (richiede YALMIP + ECOS), quindi andrebbe marcato come tale e saltato con
grazia se i pacchetti non ci sono.
