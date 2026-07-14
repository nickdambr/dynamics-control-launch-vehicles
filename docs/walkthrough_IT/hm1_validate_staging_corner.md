# HM1/validate_staging_corner.m

## Ruolo del file nel progetto

Diciamolo subito e senza giri di parole: **questo non e' uno script di produzione,
e' uno script di sanity-check.** Non calcola nessun risultato nuovo per il report,
non genera figure, non entra in nessuna suite `matlab.unittest` (non e' una classe
di test, e' uno script `clear; close all; clc` con `fprintf` finale). Il suo unico
scopo e' **verificare che l'ottimo trovato dallo sweep di `main_task4.m` sia
davvero il punto stazionario che la teoria prevede**, e non un artefatto della
griglia o di un errore di bookkeeping delle masse.

Il meccanismo e' una **cross-validation fra due strade indipendenti** verso lo
stesso numero:

- **strada A (`main_task4.m`)** -- ottimizzazione di *ordine zero*: fissa `ts`,
  risolvi il BVP a 4 incognite, leggi il payload; ripeti su 50 valori di `ts`;
  prendi il `max()`;
- **strada B (questo file)** -- condizione necessaria di *primo ordine*: promuovi
  `ts` a **quinta incognita** dello shooting e aggiungi la **corner condition**
  (Weierstrass-Erdmann all'interior junction) come **quinto residuo**; una sola
  chiamata a `fsolve` restituisce direttamente `ts*`.

Se A e B danno lo stesso `ts`, allora: (i) la corner condition e' stata derivata
giusta; (ii) lo sweep non e' finito su un massimo spurio; (iii) le formule di
payload e di massa jettisonata sono coerenti fra i due file. Se divergessero, uno
dei due sarebbe sbagliato -- e non si saprebbe quale, ma si saprebbe che c'e' un
problema. Questo e' **l'invariante giusto** da controllare, perche' e' l'unico che
mette alla prova la *derivazione analitica* (la parte in cui e' facile sbagliare
segni e normalizzazioni) contro un risultato ottenuto senza usarla.

Il file corrisponde allo **Step 7 dell'Appendice C del report**
(`HM1/report/chapters/AppendixStaging.tex`), come dichiara la riga 1. Riusa
`ode_burn.m` -- stesso RHS, stessa fisica, stessi parametri di `main_task4.m`
(righe 16-19): se cambiassero i parametri, il confronto non sarebbe piu' valido --
anche perche' i valori di riferimento dello sweep sono hardcoded (righe 51-53).

C'e' pero' un **secondo controllo, piu' sottile e piu' interessante**, nascosto
alle righe 62-72: una *cautionary check* che dimostra che la forma della corner
condition deve essere **invariante di gauge**. Ne parliamo alla fine -- e' il pezzo
piu' istruttivo del file.

---

## La condizione che stiamo validando

Riportata verbatim nell'header (righe 7-9):

    eta*[lam_m(tf) - lam_m(ts)] = c*|lam_v(ts)|*(1/m_plus - 1/m_minus)
    con  m_minus = 1 - Q*ts,   m_plus = 1 - (1+eta)*Q*ts

Da dove viene (derivazione completa nella pagina di `main_task4.m`, qui il
succo). Allo staging la mappa di giunzione e' `m(ts+) = m(ts-) - eta*Q*ts`, con
tutto il resto continuo. Siccome il salto **non dipende dallo stato** ma solo dal
*tempo* di giunzione:

- **tutti i costati sono continui** (`dDelta/dx = 0` => la Jacobiana della mappa di
  salto e' l'identita');
- **l'Hamiltoniana salta**, di

      H(ts+) - H(ts-) = T*|lam_v(ts)|*(1/m_plus - 1/m_minus)  > 0

  (l'unico termine discontinuo in `H = lam_y*vy - lam_vy + T*(|lam_v|/m - lam_m/c)`
  e' `T*|lam_v|/m`, perche' `m` scende);
- il coefficiente di `dts` nella variazione fissa **quanto** deve saltare `H`:

      H(ts+) - H(ts-) = eta*Q*[lam_m(tf) - lam_m(ts)]

Uguagliando le due espressioni e dividendo per `Q` (con `T = c*Q`) si ottiene la
corner condition. **Entrambi i membri sono positivi** (`lam_m` e' strettamente
crescente perche' `lam_m_dot = (T/m^2)*|lam_v| > 0`, e `1/m_plus > 1/m_minus`
perche' `m_plus < m_minus`): e' un buon controllo di segno gratuito.

Lettura fisica: a destra il **beneficio propulsivo** dello staging (il salto di
accelerazione `T/m`, prezzato dal primer vector `|lam_v|` = quanto vale un
incremento di `Delta v` in quell'istante); a sinistra la **penalita' strutturale**
(stadiare piu' tardi significa portarsi dietro `eta*Q` di massa inerte in piu' per
unita' di tempo). All'ottimo si pareggiano, che e' esattamente `dtf/dts = 0`.

---

## `%% Parameters` (righe 15-27)

```matlab
c = 0.6; eta = 0.1; yf = 0.04; Q = 2; T = c * Q;
opts_ode = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
opts_fs  = optimoptions('fsolve', 'Display', 'off', ...
    'MaxIterations', 500, 'MaxFunctionEvaluations', 5000, ...
    'FunctionTolerance', 1e-12, 'StepTolerance', 1e-12);
```

- Righe 16-20: **identici a `main_task4.m`** (il commento della riga 15 lo dichiara).
  E' condizione necessaria per il confronto: se i parametri divergessero, i due `ts`
  non avrebbero motivo di coincidere.
- Riga 22: tolleranze dell'integratore (`RelTol 1e-10`, `AbsTol 1e-12`), **identiche**
  a quelle di `main_task4.m` (riga 14).
- Righe 23-25: tolleranze di `fsolve` **piu' strette** di `main_task4.m` (riga 17):
  `1e-12` invece di `1e-10` su funzione e passo (i cap `MaxIterations = 500` e
  `MaxFunctionEvaluations = 5000` sono invece gli stessi). Sensato: qui la corner
  condition e' un residuo *aggiuntivo* e piu' delicato dei tre terminali, e vogliamo
  che `ts` esca con piu' cifre significative dello sweep con cui lo confrontiamo.

---

## `%% Step 1: warm-start from a fixed-ts inner BVP` (righe 29-36)

```matlab
ts0 = 0.33;
p_inner = p; p_inner.ts = ts0;
z4_guess = [0.6; 3.8; 14; 0.42];
[z4, ~, ef1] = fsolve(@(z) shooting_inner(z, p_inner, ...
                      opts_ode), z4_guess, opts_fs);
```

- Riga 30: `ts0 = 0.33`, cioe' **gia' vicinissimo** alla risposta attesa (0.336).
  Va detto con onesta': il warm start e' *informato dalla risposta*. Non e' una
  circolarita' logica (il solve a 5 incognite deve comunque far annullare la corner
  condition, che a `ts=0.33` **non** e' soddisfatta), ma e' una scorciatoia: se si
  partisse da `ts0 = 0.10` il bacino di convergenza dello shooting indiretto
  potrebbe non bastare. Lo script non lo dice.
- Riga 32: guess `[lam_vx0; lam_vy0; lam_y; tf] = [0.6; 3.8; 14; 0.42]`. I primi
  tre sono gli stessi di `main_task4.m`; `tf = 0.42` e' pero' gia' il valore del
  caso *a due stadi* (in `main_task4.m` il guess monostadio era `0.30`).
- Righe 33-36: si risolve prima il BVP interno a `ts` **fissato** (4 incognite,
  4 residui -- esattamente la `shooting_twostage` di `main_task4.m`). Serve solo a
  produrre un punto di partenza consistente per il solve aumentato. Se non converge,
  `error()` immediato.

---

## `%% Step 2: augmented solve` (righe 38-41)

```matlab
w_guess = [z4; ts0];
[w, res, ef] = fsolve(@(w) shooting_corner(w, p, opts_ode), ...
                      w_guess, opts_fs);
```

- Riga 40: si appende `ts0` al vettore delle incognite: da 4 a **5 incognite**
  `[lam_vx0; lam_vy0; lam_y; tf; ts]` (la riga 39 e' il commento che elenca i cinque
  residui).
- Riga 41: un solo `fsolve` su 5 residui. **Struttura del sistema 5x5** (vale la
  pena capirla, perche' spiega perche' e' ben posto):
  - i **3 residui terminali** (`y(tf)=yf`, `vx(tf)=1`, `vy(tf)=0`) dipendono dai
    costati **solo attraverso i rapporti** `lam_vy0/lam_vx0` e `lam_y/lam_vx0`
    (l'angolo `phi = atan2(lam_vy, lam_vx)` e' invariante per riscalamento
    positivo). Insieme a `tf`, quei tre residui determinano **la traiettoria**, dato
    `ts`;
  - il residuo **corner** determina `ts` -- ed e' l'unico che porta informazione di
    ottimalita' sulla giunzione;
  - il residuo **`H0 = 0`** fissa la **scala** dei costati (che altrimenti resterebbe
    libera per omogeneita': le condizioni necessarie sono omogenee di grado 1 in
    `lambda`). Non e' -- nel caso a due stadi -- la condizione di trasversalita' per
    `tf` libero, perche' `H` salta a `ts` e quindi `H(0)=0` non implica `H(tf)=0`.
    E' un **gauge fixing**, e la sua innocuita' e' garantita dal fatto che la corner
    condition e' scritta in forma invariante di gauge (vedi sotto).

---

## `%% Extract and report` (righe 43-60)

```matlab
payload = mf - eta*Q*(tf - ts);
...
if ef > 0 && abs(ts - 0.336) < 5e-3
    fprintf('\nPASS: corner solve matches the swept optimum.\n');
```

- Riga 47: stessa formula di payload di `main_task4.m` (riga 106): `mf` meno la
  struttura del **secondo** stadio, `ms2 = eta*Q*(tf-ts)`. Coerente per costruzione.
- Righe 51-53: i valori di riferimento dello sweep sono **hardcoded nelle stringhe
  di formato**: `(sweep: 0.336)`, `(sweep: 0.424)`, `(sweep: 0.068)`.
- Riga 56: il **criterio di PASS**: `ef > 0 && abs(ts - 0.336) < 5e-3`.

Qui vanno fatte tre osservazioni oneste.

1. **Lo script non ri-esegue `main_task4.m`.** Confronta contro tre *costanti
   letterali* trascritte a mano da una run precedente. Se cambiassero `c`, `eta`,
   `yf` o `Q`, quelle costanti diventerebbero stantie e il test darebbe FAIL pur
   essendo il codice corretto. E' un **oracolo registrato**, non un oracolo
   ricalcolato. E' una scelta legittima per uno script di verifica una-tantum, ma
   e' una fragilita' da conoscere.
2. **La tolleranza `5e-3` non e' arbitraria, ed e' la scelta giusta.** Lo sweep di
   `main_task4.m` usa 50 punti su un intervallo di circa `[0.01, 0.42]`, quindi ha
   un passo di griglia di circa `0.008`: il suo `ts_opt` e' **agganciato a un nodo**
   e non puo' essere piu' preciso di ~mezzo passo. Una tolleranza di `5e-3` e'
   esattamente dell'ordine di mezzo passo di griglia. Chiedere `1e-6` sarebbe stato
   scorretto: si starebbe pretendendo dallo sweep una precisione che non ha.
   (Corollario onesto: e' il **solve di questo file** ad essere accurato, non lo
   sweep. Il numero "buono" di `ts*` e' quello che esce da qui.)
3. **Il `norm(res)` viene stampato ma non usato come gate** (riga 54). Un test piu'
   severo avrebbe messo anche `norm(res) < 1e-8` nella condizione di PASS: `fsolve`
   puo' restituire `ef > 0` con criteri di arresto (es. `StepTolerance`) che non
   garantiscono un residuo piccolo. Nella pratica il residuo esce a `1e-10` o meglio
   (cosi' riporta il report), ma il gate non lo verifica.

Da segnalare anche una piccola incoerenza documentale sul valore attribuito allo
**sweep**: il README di HM1 e la tabella dei risultati di Task 4 nel report
(`Task4.tex`, `ts* = 0.337`) danno `0.337`, mentre le stringhe hardcoded di questo
script (righe 51-53, etichettate proprio `sweep:`), il testo di `Task4.tex` e lo
Step 7 dell'Appendice C danno `0.336`, che e' il numero del **corner solve**. E'
una differenza nel terzo decimale, del tutto compatibile con l'aggancio alla griglia
dello sweep (passo ~`0.008`), ma la provenienza dei due numeri va tenuta distinta:
`0.337` = nodo della griglia, `0.336` = zero della corner condition.

---

## `%% Check: the burnout reference lam_m(tf) is essential` (righe 62-72)

```matlab
[w_bad, ~, ef_bad] = fsolve(@(w) shooting_corner_wrong(w, p, ...
                            opts_ode), [z4; ts0], opts_fs);
fprintf('  eta*lam_m(ts) alone => ts = %.6f (spurious)\n', ...
        w_bad(5));
```

**Questa e' la parte piu' importante del file**, e non e' un semplice
"controesempio decorativo": e' una dimostrazione numerica di **invarianza di
gauge**.

Il punto. Il BVP integra `lam_m` a partire da `lam_m(0) = 1` (riga 87). Quel `1`
**non e' fisica, e' una convenzione**. Le condizioni necessarie hanno una liberta'
a **due parametri**:

    lam -> k*lam        (k > 0: scala globale di tutti i costati)
    lam_m -> lam_m + b  (offset additivo del solo lam_m)

L'offset `b` e' ammissibile perche' `lam_m` entra in `H` solo tramite `-Q*lam_m`
(un termine costante nel tempo... a meno del salto, che pero' e' preservato) e
perche' `lam_m_dot = (T/m^2)*|lam_v|` **non dipende da `lam_m` stesso**: la sua ODE
e' un'integrazione pura, quindi cambiare la costante d'integrazione trasla l'intera
curva `lam_m(t)` rigidamente. Nel report, la formulazione a minimo `tf` porta a
`lam_m(tf) = 0`; nella formulazione a massimo payload porta a `lam_m(tf) = 1`; il
codice usa `lam_m(0) = 1`. **Sono tre gauge dello stesso problema**, e la
traiettoria e' la stessa in tutti e tre.

Conseguenza: qualunque condizione di ottimalita' scritta in termini di `lam_m` deve
essere **invariante** sotto `(k, b)`, altrimenti il suo zero dipende da una
convenzione arbitraria. La forma

    eta*[lam_m(tf) - lam_m(ts)] = c*|lam_v(ts)|*(1/m_plus - 1/m_minus)

lo e':
- la **differenza** `lam_m(tf) - lam_m(ts)` **cancella l'offset `b`**;
- entrambi i membri sono **omogenei di grado 1 in `k`** (`lam_m(tf)-lam_m(ts)` e'
  l'integrale di `(T/m^2)*|lam_v|`, che scala con `k`, e a destra `|lam_v(ts)|`
  scala con `k`).

La forma

    eta*lam_m(ts) = c*|lam_v(ts)|*(1/m_plus - 1/m_minus)

(riga 185, `shooting_corner_wrong`) **non lo e'**: `lam_m(ts)` da solo si porta
dietro l'offset arbitrario `b = 1` ereditato da `lam_m0 = 1`. Il suo zero e' quindi
**una funzione della convenzione**, e infatti cade altrove. L'appendice del report
riporta il valore spurio: `ts ~ 0.225` invece di `0.336`.

- Righe 65-72: il solve "sbagliato" viene lanciato con lo **stesso warm start**
  `[z4; ts0]` del solve corretto, e converge (`ef_bad > 0`) -- ecco perche' l'errore
  e' insidioso: **non esplode**, ti da' un numero perfettamente plausibile.
- Riga 71: se non converge, lo script lo dice e basta. Non c'e' assertion su
  `w_bad(5)`: il check e' **informativo**, stampa e non fallisce mai. Onesta': non e'
  un test, e' una demo.

> **Possibile domanda d'esame** -- perche' nella corner condition compare la
> *differenza* `lam_m(tf) - lam_m(ts)` e non semplicemente `lam_m(ts)`?
> *Risposta:* perche' la condizione deve essere indipendente dalla normalizzazione
> dei costati. `lam_m` obbedisce a `lam_m_dot = (T/m^2)*|lam_v|`, che non contiene
> `lam_m`: cambiare la sua costante d'integrazione trasla rigidamente tutta la
> curva. Il mio shooting parte da `lam_m0 = 1`, la teoria a minimo tempo vuole
> `lam_m(tf) = 0`, quella a massimo payload vuole `lam_m(tf) = 1`. Scrivendo la
> condizione come differenza, l'offset si cancella e posso valutarla direttamente
> sui costati che il BVP integra. `shooting_corner_wrong` dimostra numericamente
> cosa succede se lo si dimentica: converge lo stesso, ma su un `ts` spurio
> (~0.225 invece di 0.336).

---

## `propagate` (righe 75-94)

```matlab
ic = [0;0;0;0;1;1];                       % lam_m0 = 1
[~, Z1] = ode45(@(t,z) ode_burn(t,z,pp), [0 ts], ic, opts_ode);
zs_minus = Z1(end,:)';
z_plus = zs_minus;
z_plus(5) = z_plus(5) - p.eta*p.Q*ts;     % jettison
[~, Z2] = ode45(@(t,z) ode_burn(t,z,pp), [ts tf], z_plus, opts_ode);
zf = Z2(end,:)';
```

- Riga 75: firma. Restituisce **due** cose: lo stato finale `zf` e lo stato a
  `ts^-` (`zs_minus`), *prima* del jettison. Serve entrambi: `zf(6) = lam_m(tf)` e
  `zs(6) = lam_m(ts)` per la corner condition, `zs(5) = m_minus` per il salto di
  massa. E' il **refactoring** che `main_task4.m` non ha fatto (la' la propagazione
  e' duplicata a mano dentro `shooting_twostage` e di nuovo nel corpo dello script).
- Righe 84-85 (dove `pp` viene costruita) + righe 88 e 92 (le due chiamate a
  `ode45`): la **stessa struct `pp`** viene passata a `ode_burn` su entrambi gli
  archi. E' cosi' che i costati `lam_vx`, `lam_vy`, `lam_y` restano continui: sono
  formule analitiche nei parametri, non stati integrati, quindi non hanno alcun
  kink a `ts`.
- Riga 87: `lam_m(0) = 1`. Il gauge.
- Riga 91: **il salto**. Solo la componente 5 (massa). `z_plus(6) = lam_m` **non
  viene toccato**: e' la continuita' del costato di massa, che vale perche' la massa
  jettisonata `eta*Q*ts` dipende solo da `ts` e non da `m(ts)`.
- Riga 92: arco 2 con lo **stesso** RHS: stesso `T`, stesso `Q`, stesso `c`. L'unica
  cosa che cambia e' `m`, quindi l'accelerazione `T/m` fa un gradino verso l'alto.

> **Possibile domanda d'esame** -- se avessi modellato il jettison come "butto una
> frazione `k` della massa corrente", cosa sarebbe cambiato?
> *Risposta:* tutto. Con `m(ts+) = (1-k)*m(ts-)` la mappa di salto dipende dallo
> stato, la sua Jacobiana rispetto a `m(ts-)` vale `1-k` e non `1`, e la condizione
> di corner `lam(ts-)' = lam(ts+)' * dg/dx-` **non** collassa piu' in continuita':
> si otterrebbe `lam_m(ts-) = (1-k)*lam_m(ts+)`, cioe' **il costato di massa
> salterebbe**. E' proprio il fatto che la massa buttata sia *assoluta* (funzione del
> solo `ts`) a lasciare passare `lam_m` intatto.

---

## `shooting_inner` (righe 96-121)

- Riga 96: 4 incognite, 4 residui, `ts` fisso in `p.ts`. E' la **stessa formulazione**
  di `shooting_twostage` in `main_task4.m`, riscritta sopra `propagate`.
- Righe 106-108: guardia di bracket (`tf <= ts || tf > 2 || ts <= 0`) che restituisce
  `1e6*ones(4,1)`.
- Righe 110-117: **qui c'e' un difetto di ordinamento**. La guardia su `m_plus > 0`
  (riga 115) e' eseguita **dopo** `propagate` (riga 111), cioe' dopo che l'arco 2 e'
  gia' stato integrato -- potenzialmente con `m <= 0`, che in `ode_burn` finisce nei
  denominatori `T/m` e `T/m^2`. In `main_task4.m` (`shooting_twostage`, righe 282-284)
  la stessa guardia e' correttamente **prima** dell'integrazione dell'arco 2. In
  pratica non morde, perche' il `try/catch` cattura gli errori e i guess restano
  lontani dal bordo, ma se `ode45` restituisse `NaN` senza sollevare eccezione il
  `catch` non scatterebbe e il residuo sarebbe `NaN`.
- Righe 118-120: residui. `H0` calcolato dalla forma algebrica a `t=0`
  (`vx=vy=0, m=1, lam_m0=1`), e i tre terminali `y(tf)-yf`, `vx(tf)-1`, `vy(tf)`.

---

## `shooting_corner` (righe 123-154)

```matlab
lam_vy_ts = lam_vy0 - lam_y*ts;
lam_v_ts  = sqrt(lam_vx0^2 + lam_vy_ts^2);
corner = p.eta*(zf(6) - zs(6)) ...
         - p.c*lam_v_ts*(1/m_plus - 1/m_minus);
res = [zf(2)-p.yf; zf(3)-1; zf(4); H0; corner];
```

- Riga 123: **il cuore del file**: 5 incognite `[lam_vx0; lam_vy0; lam_y; tf; ts]`,
  5 residui.
- Riga 133: guardia estesa. Oltre a `tf > ts > 0` e `tf < 2`, c'e'
  `ts > 1/(Q*(1+eta)) - 1e-3` -> reject. E' il vincolo fisico `m_plus > 0`:
  `m_plus = 1 - (1+eta)*Q*ts` si annulla a `ts = 1/(Q*(1+eta)) = 0.4545`. Lo stesso
  bound compare in `main_task4.m` alla riga 45.
- Righe 141-142: `m_minus = zs(5)` (letto dall'integrazione) e
  `m_plus = m_minus - eta*Q*ts`. Nota: `m_minus` **potrebbe** essere calcolato in
  forma chiusa (`1 - Q*ts`, perche' `mdot = -Q` costante e `m0 = 1`) -- infatti
  l'header lo fa (riga 9). Leggerlo dall'integrazione e' piu' robusto e serve anche
  come controllo incrociato implicito: se l'integratore fosse fuori tolleranza,
  `zs(5)` e `1 - Q*ts` divergerebbero.
- Righe 150-151: **`lam_v(ts)` viene ricostruito analiticamente**, non letto dallo
  stato. E' obbligatorio: `ode_burn` **non integra** `lam_vx` e `lam_vy` (li ricava
  dai parametri), quindi non esistono come componenti di `z`. La formula e' la
  soluzione chiusa delle equazioni di Eulero-Lagrange:

      lam_vx(t)  = lam_vx0            (costante)
      lam_vy(t)  = lam_vy0 - lam_y*t  (rampa lineare)
      |lam_v(ts)| = sqrt(lam_vx0^2 + (lam_vy0 - lam_y*ts)^2)

- Riga 152: il residuo di corner, nella forma **burnout-referenced**
  `eta*(lam_m(tf) - lam_m(ts))`, con `zf(6) = lam_m(tf)` e `zs(6) = lam_m(ts)`.
  Osservare che `zs` e' lo stato a `ts^-`, ma `lam_m` e' continuo, quindi
  `lam_m(ts^-) = lam_m(ts^+)` e la scelta e' irrilevante -- cosa che, dovendo scriverlo
  in un report, va detta esplicitamente.
- Riga 153: l'ordinamento dei residui. I primi 4 sono identici a quelli del BVP
  interno; il quinto e' la novita'.

> **Possibile domanda d'esame** -- nel solve a 5 incognite, chi determina cosa?
> *Risposta:* i 3 residui terminali determinano la traiettoria (i due rapporti dei
> costati e `tf`), dato `ts`; il residuo di corner determina `ts`; il residuo
> `H0 = 0` fissa la scala globale dei costati, che i residui terminali non vedono
> (dipendono solo dalla *direzione* del primer vector). Il sistema e' ben posto
> proprio perche' i tre gruppi sono indipendenti. Il residuo di corner puo' essere
> valutato con i costati "sporchi" (`lam_m0 = 1`, `H(0) = 0`) perche' e' scritto in
> forma invariante di gauge.

---

## `shooting_corner_wrong` (righe 156-187)

- Riga 156: **copia deliberatamente corrotta** di `shooting_corner`. Unica differenza:
  la riga 185, dove `eta*(zf(6) - zs(6))` diventa `eta*zs(6)`, cioe'
  `eta*lam_m(ts)` senza il riferimento a burnout.
- Onesta' massima: **non e' un "quasi-errore" che qualcuno commetterebbe per
  distrazione derivando la teoria** -- e' un residuo costruito apposta come
  *esperimento di controllo*. Serve a rispondere alla domanda "e se avessi
  dimenticato il termine `lam_m(tf)`?", e la risposta e' "converge lo stesso, e ti
  da' un `ts` sbagliato". Il valore didattico e' esattamente questo: **l'errore non
  si manifesta come un fallimento numerico**, si manifesta come un risultato
  plausibile ma falso. Senza il confronto con lo sweep non te ne accorgeresti.
- Il resto della funzione (righe 167-184) e' identico a `shooting_corner`,
  duplicazione di codice inclusa. Accettabile per uno script di validazione, ma e'
  copia-incolla.

---

## Cosa questo script valida davvero, e cosa no

**Valida:**
- che la corner condition derivata nel report abbia lo **stesso zero** del massimo
  numerico trovato dallo sweep di `main_task4.m` (con `ts` a `5e-3`, cioe' entro la
  risoluzione dello sweep);
- che segni, normalizzazioni e bookkeeping delle masse (`m_minus`, `m_plus`,
  `ms1 = eta*Q*ts`) siano coerenti fra derivazione analitica e codice;
- che la forma **burnout-referenced** sia necessaria (righe 62-72).

**Non valida:**
- che l'ottimo sia un **massimo** e non un minimo o un flesso -- la corner condition
  e' una condizione **necessaria del primo ordine**, non sufficiente. E' lo sweep di
  `main_task4.m` a mostrare che la curva payload-vs-`ts` ha un massimo li'. (I due
  file quindi si sostengono a vicenda: lo sweep da' l'ordine zero e la natura del
  punto, il corner solve la precisione. Nessuno dei due basta da solo.)
- che la soluzione sia **globale**: `fsolve` e' un metodo locale, e il warm start e'
  a `ts0 = 0.33`.
- niente **al di fuori del singolo set di parametri** `c=0.6, eta=0.1, yf=0.04, Q=2`.
  Non c'e' nessuno sweep su `Q` o `yf` che verifichi la robustezza della condizione.
- il **residuo finale**, che viene stampato ma non incluso nel criterio di PASS.

---

## Possibili domande d'esame

**D: A cosa serve `validate_staging_corner.m`? Non e' ridondante rispetto a
`main_task4.m`?**
R: No, e' proprio il contrario: e' l'unico controllo indipendente che ho.
`main_task4.m` trova `ts*` con una ricerca esaustiva su griglia -- un metodo di
ordine zero che non usa mai la teoria di Pontryagin sulla giunzione. Questo script
la usa e basta: promuove `ts` a quinta incognita e impone la corner condition come
quinto residuo. Il fatto che le due strade cadano sullo stesso `ts` (0.336 dal corner
solve contro 0.336-0.337 dello sweep, cioe' entro la risoluzione della griglia) e'
l'invariante che verifico. Se avessi
sbagliato un segno o una normalizzazione nella derivazione, i due numeri
divergerebbero.

**D: Perche' la tolleranza del PASS e' `5e-3` e non `1e-8`?**
R: Perche' il termine di paragone non e' accurato. Lo sweep di `main_task4.m` usa
50 punti su circa `[0.01, 0.42]`, quindi un passo di ~`0.008`: il suo `ts_opt` e'
un nodo della griglia e non puo' essere piu' preciso di mezzo passo. `5e-3` e'
proprio quell'ordine di grandezza. Chiedere `1e-8` significherebbe pretendere dallo
sweep una precisione che per costruzione non ha. Nota il rovescio: e' il corner
solve ad essere il numero *buono*; lo sweep e' il numero *approssimato*.

**D: I costati sono continui allo staging? E l'Hamiltoniana?**
R: I costati sono **tutti continui**, perche' la massa jettisonata
`ms1 = eta*Q*ts` dipende solo dal tempo di giunzione e non dallo stato: la
Jacobiana della mappa di salto rispetto a `x(ts-)` e' l'identita' e la condizione
`lam(ts-)' = lam(ts+)' * dg/dx-` diventa `lam(ts-) = lam(ts+)`. Nel codice si vede
alla riga 91: si tocca **solo** `z_plus(5)`, la massa; `z_plus(6) = lam_m` passa
intatto. L'Hamiltoniana invece **salta**, di
`T*|lam_v(ts)|*(1/m_plus - 1/m_minus) > 0`, perche' la massa a denominatore in
`T*|lam_v|/m` scende di colpo. Ed e' proprio l'entita' di quel salto che la corner
condition prescrive.

**D: Perche' `lam_v(ts)` e' ricalcolato a mano (righe 150-151) invece di leggerlo
dallo stato integrato?**
R: Perche' non e' nello stato. `ode_burn` integra `[x; y; vx; vy; m; lam_m]`: gli
altri costati hanno **soluzione in forma chiusa** (`lam_vx = lam_vx0` costante,
`lam_vy = lam_vy0 - lam_y*t` rampa lineare, `lam_x = 0`, `lam_y` costante) e non
vengono propagati numericamente. `lam_m` e' l'unico che richiede un'integrazione,
perche' la sua ODE `lam_m_dot = (T/m^2)*|lam_v|` dipende dalla massa. Quindi
`|lam_v(ts)|` si valuta direttamente dalla formula, il che e' anche piu' accurato
che leggerla da un integratore.

**D: Cosa dimostra il "cautionary check" delle righe 62-72?**
R: Che la corner condition deve essere scritta in forma **invariante di gauge**. Il
mio BVP integra `lam_m` da `lam_m0 = 1`, che e' una convenzione arbitraria: siccome
`lam_m_dot` non dipende da `lam_m`, cambiare la costante d'integrazione trasla
rigidamente tutta la curva. Se scrivo la condizione con `lam_m(ts)` da solo, il suo
zero dipende da quella convenzione ed e' privo di significato; scrivendola con la
**differenza** `lam_m(tf) - lam_m(ts)` l'offset si cancella. Il check lo mostra
numericamente: la versione senza riferimento converge senza errori, ma su un `ts`
spurio (~0.225 secondo l'appendice del report) invece di 0.336. E' un errore che
non fa rumore, ed e' per questo che vale la pena averlo messo nero su bianco.

**D: Il solve a 5 incognite garantisce che quello sia un massimo del payload?**
R: No. La corner condition e' una condizione **necessaria del primo ordine**
(`dtf/dts = 0`): identifica un punto stazionario, non ne certifica la natura. E' lo
sweep di `main_task4.m` a mostrare che la curva payload-vs-`ts` ha li' un massimo
interno e non un minimo. I due file sono complementari: lo sweep da' la *forma*
della curva, il corner solve la *posizione precisa* dell'estremo.
