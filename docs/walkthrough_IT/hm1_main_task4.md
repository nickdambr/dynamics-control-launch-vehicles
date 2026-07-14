# HM1/main_task4.m

## Ruolo del file nel progetto

`main_task4.m` e' lo script del **Task 4 di HM1: staging ottimo**. La variante e'
un lanciatore a **due stadi** che parte da fermo all'origine (nessuna salita
verticale, a differenza di Task 2 e Task 3), accende il motore e lo tiene acceso
per tutto il volo. A un istante `ts` -- l'unica novita' rispetto al Task 1 -- la
struttura del primo stadio viene **buttata via** (jettison): la massa cade di
colpo di `ms1 = eta*Q*ts`, e il volo prosegue con lo stesso motore (stesso `T`,
stesso `Q`, stesso `c`) fino all'iniezione.

Il file esiste per rispondere a una domanda quantitativa: **quanto guadagna il
payload lo staging, e quando conviene farlo?** La struttura e' a **due livelli**:

- **livello interno (BVP)** -- per un `ts` *fissato* si risolve lo stesso shooting
  a 4 incognite del Task 1, ma su due archi separati dal salto di massa
  (`shooting_twostage`, righe 248-298);
- **livello esterno (sweep)** -- `ts` viene spazzato su una griglia di 50 punti e
  per ognuno si legge il payload risultante; il massimo della curva e' lo staging
  ottimo (righe 38-131).

Cioe': lo script **non impone la condizione di ottimalita' su `ts`**. La cerca
numericamente. La condizione analitica corrispondente (la *corner condition*
all'istante di giunzione) esiste e viene derivata piu' sotto, ma e' implementata
in un file separato, `validate_staging_corner.m`, che serve appunto a verificare
che lo sweep di questo script cada dove la teoria dice.

Dipende da `ode_burn.m` (RHS condiviso: dinamica planare + costati + legge di
controllo linear-tangent gia' sostituita) e da nient'altro. Non c'e' atmosfera,
non c'e' drag, la spinta e' sempre accesa: `T = c*Q`, `g = 1` nondimensionale.

---

## `%% Parameters` (righe 7-17)

```matlab
c   = 0.6;
eta = 0.1;
yf  = 0.04;
Q   = 2;         % same Q for both stages
T   = c * Q;
```

- Righe 8-12: i dati della traccia in unita' nondimensionali. `c = 0.6` e'
  la velocita' efficace di scarico (Isp*g0 / Vrif), `eta = 0.1` il coefficiente
  strutturale `ms/mp`, `yf = 0.04` la quota target, `Q = 2` la portata di massa.
  `T = c*Q = 1.2` e' la spinta, costante e sempre accesa.
- Riga 11: il commento `same Q for both stages` va preso sul serio: e' una
  **semplificazione di modello**. Un lanciatore reale ha un secondo stadio con
  spinta e portata molto diverse (e spesso Isp maggiore, motore da vuoto). Qui i
  due stadi condividono `T`, `Q`, `c`: l'unica cosa che cambia allo staging e' la
  massa. E' anche il motivo per cui il guadagno di payload che esce (+182% nel
  README) e' cosi' generoso rispetto ai lanciatori veri.
- Righe 14-17: tolleranze. `ode45` a `RelTol=1e-10, AbsTol=1e-12`, `fsolve` a
  `1e-10` su funzione e passo. Serve: i residui dello shooting sono differenze
  fra numeri O(1) e il Jacobiano di `fsolve` e' calcolato per differenze finite
  sul risultato di un'integrazione -- se l'integratore e' lasco, il Jacobiano e'
  rumore.

---

## `%% First solve single-stage (Task 1) as reference` (righe 19-36)

```matlab
z_guess_1 = [0.6; 3.8; 14; 0.30];
[z_ref, ~, ef] = fsolve(@(z) shooting_single(z, p1, ...
                        opts_ode), z_guess_1, opts_fs);
...
payload_single = mf_single*(1+eta) - eta;
```

- Righe 21-24: si risolve prima il **caso monostadio** (identico al Task 1). Serve
  a due cose: (a) e' il termine di paragone per il guadagno di staging; (b) la sua
  soluzione `z_ref` e' il **warm start** dello sweep a due stadi (riga 77).
- Riga 23: il guess `[lam_vx0; lam_vy0; lam_y; tf] = [0.6; 3.8; 14; 0.30]`. Con
  `lam_vy0/lam_vx0 = 3.8/0.6` l'angolo di spinta iniziale e'
  `phi(0) = atan2(3.8, 0.6) ~ 81 gradi`: quasi verticale, come dev'essere.
- Riga 28: `ic = [0;0;0;0;1;1]` -- stato iniziale `[x; y; vx; vy; m; lam_m]`, cioe'
  origine, fermo, `m0 = 1` (nondimensionalizzato su se stesso) e `lam_m(0) = 1`.
  Quest'ultimo **non e' un dato fisico**: e' la normalizzazione dei costati (vedi
  la sezione sulla matematica).
- Riga 29: la soluzione convergente viene **ri-integrata** per estrarre `mf`. E'
  lavoro duplicato (`shooting_single` aveva gia' integrato la stessa traiettoria),
  ma costa una frazione di secondo ed evita di far restituire lo stato al residuo.
- Riga 32: `payload_single = mf*(1+eta) - eta`. Derivazione: la massa iniziale si
  ripartisce in `m0 = 1 = mu + ms + mp` (payload + struttura + propellente). Il
  propellente bruciato e' `mp = Q*tf = 1 - mf`, la struttura e' `ms = eta*mp`.
  Quindi

      mu = mf - ms = mf - eta*(1 - mf) = mf*(1 + eta) - eta

  Il payload e' cio' che resta **dopo aver scartato anche la struttura**: e' una
  *differenza di numeri grandi* (qui `mf ~ 0.11`, `mu ~ 0.024`), e questo spiega
  perche' guadagni percentualmente enormi sul payload corrispondano a variazioni
  minuscole di `tf`.

> **Possibile domanda d'esame** -- perche' `payload = mf*(1+eta) - eta` e non
> semplicemente `mf`?
> *Risposta:* perche' `mf` include ancora la struttura (serbatoi, motore) del
> lanciatore. Il modello a coefficiente strutturale dice che la struttura e'
> proporzionale al propellente imbarcato, `ms = eta*mp`; siccome `mp = 1 - mf`,
> si ricava `ms = eta*(1 - mf)` e il payload utile e' `mf - ms`. Se si riportasse
> `mf` come figura di merito si premierebbero i lanciatori che bruciano poco
> propellente ma sono tutti serbatoio.

---

## La matematica dello staging: salto di massa, costati, Hamiltoniana

Questa e' la parte che il codice usa ma non scrive. Vale la pena derivarla per
intero, perche' e' il cuore della domanda "i costati sono continui allo staging?".

### Il problema di controllo ottimo di base

Con lo stato `x = [x; y; vx; vy; m]` e la spinta sempre accesa, la dinamica
(in `ode_burn.m`) e'

    xdot  = vx
    ydot  = vy
    vxdot = (T/m)*cos(phi)
    vydot = (T/m)*sin(phi) - 1        (g = 1 nondim)
    mdot  = -Q

L'Hamiltoniana (convenzione di massimizzazione, quella usata da `ode_burn`) e'

    H = lam_x*vx + lam_y*vy
        + lam_vx*(T/m)*cos(phi)
        + lam_vy*((T/m)*sin(phi) - 1)
        - lam_m*Q

Le equazioni di Eulero-Lagrange (`lam_dot = -dH/dx`) danno:

    lam_x_dot  = 0                 ->  lam_x  = const = 0  (x(tf) libero)
    lam_y_dot  = 0                 ->  lam_y  = const
    lam_vx_dot = -lam_x = 0        ->  lam_vx = lam_vx0
    lam_vy_dot = -lam_y            ->  lam_vy = lam_vy0 - lam_y*t
    lam_m_dot  = (T/m^2)*|lam_v|   ->  lam_m crescente (unica ODE da integrare)

La massimizzazione su `phi` allinea la spinta al **primer vector** `lam_v`:
`phi = atan2(lam_vy, lam_vx)`, da cui `tan(phi)` lineare in `t` (legge
*linear-tangent*). Sostituendo il `phi` ottimo si ottiene la forma compatta

    H = lam_y*vy - lam_vy + T*( |lam_v|/m - lam_m/c )

usando `Q = T/c`. A `t = 0` (dove `vy = 0`, `m = 1`, `lam_m0 = 1`) diventa
**algebrica**:

    H(0) = -lam_vy0 + T*( |lam_v(0)| - 1/c )

che e' esattamente la riga 243 (e 295) del codice.

### Il salto di massa

Allo staging la mappa di giunzione e'

    m(ts+) = m(ts-) - eta*Q*ts,        x, y, vx, vy continui

Perche' `eta*Q*ts`? Il propellente bruciato dal primo stadio e' `mp1 = Q*ts`; il
modello a coefficiente strutturale dice `ms1 = eta*mp1 = eta*Q*ts`. **Il salto
dipende esplicitamente da `ts`**: questo e' il punto tecnico chiave, e' cio' che
rende la giunzione diversa da un semplice cambio di fase.

Siccome `mdot = -Q` e `m0 = 1`:

    m_minus = 1 - Q*ts
    m_plus  = 1 - (1+eta)*Q*ts

(righe 43-44 del codice, e la condizione `m_plus > 0` da' il vincolo
`ts < 1/(Q*(1+eta))` della riga 45).

### I costati sono continui?

**Si', tutti e cinque.** Ecco perche'. Si adiunge il vincolo di giunzione
`x(ts+) = x(ts-) + Delta(ts)` con un moltiplicatore `pi`:

    Jaug = ... + pi'*[ x(ts-) + Delta(ts) - x(ts+) ] + integrali dei due archi

Il coefficiente della variazione `delta x(ts-)` da' `lam(ts-) = pi`, quello di
`delta x(ts+)` da' `lam(ts+) = pi`. Quindi

    lam(ts-) = lam(ts+)

La ragione strutturale: `Delta` **non dipende dallo stato** `x(ts-)`, dipende solo
da `ts`. In generale la condizione e' `lam(ts-)' = lam(ts+)' * dDelta/dx-`, e qui
`dDelta/dx- = I` (matrice identita'), quindi i costati passano lisci. Se il modello
avesse `ms1 = eta*(m0 - m(ts))` scritto in funzione della *massa* (invece che del
tempo), la Jacobiana non sarebbe l'identita' e `lam_m` **salterebbe**.

Nel codice questa continuita' e' imposta implicitamente e in modo pulito:
`lam_vx`, `lam_vy`, `lam_y` sono formule analitiche nei *parametri* `pp.lam_vx0`,
`pp.lam_vy0`, `pp.lam_y`, e la **stessa** struct `pp` viene passata a `ode_burn`
sui due archi (righe 90 e 98) -- quindi sono automaticamente continui. `lam_m` e'
l'unico costato integrato, e la riga 280 tocca **solo** `z_s(5)` (la massa),
lasciando `z_s(6) = lam_m` intatto (commento esplicito, riga 281).

### L'Hamiltoniana e' continua? NO.

Attraverso `ts` tutto e' continuo tranne `m`, che scende. Nella forma compatta
`H = lam_y*vy - lam_vy + T*(|lam_v|/m - lam_m/c)`, l'unico termine che cambia e'
`T*|lam_v|/m`:

    H(ts+) - H(ts-) = T*|lam_v(ts)|*( 1/m_plus - 1/m_minus )  > 0

**L'Hamiltoniana salta verso l'alto.** Fisicamente e' ovvio: buttando massa,
l'accelerazione `T/m` aumenta di colpo, e `H` (il "tasso istantaneo di accumulo
di merito") sale. Questo e' il caso classico in cui **H non e' continua** perche'
la mappa di giunzione dipende dal tempo di giunzione.

### La condizione di ottimalita' su `ts` (corner condition)

Il coefficiente di `dts` nella variazione da' la condizione mancante:

    dphi/dts + H(ts-) - H(ts+) + lam(ts)'*(dDelta/dts) = 0

dove `phi` e' il payoff terminale. Qui il payoff e' il **payload**:

    mu = m(tf) - eta*Q*(tf - ts)      ->   dphi/dts = +eta*Q

e `dDelta/dts = (0,0,0,0, -eta*Q)'` da cui `lam'*dDelta/dts = -eta*Q*lam_m(ts)`.
Sostituendo il salto di `H` calcolato sopra:

    eta*Q - T*|lam_v(ts)|*(1/m_plus - 1/m_minus) - eta*Q*lam_m(ts) = 0

e, dividendo per `Q` e usando `T = c*Q`:

    eta*( 1 - lam_m(ts) ) = c*|lam_v(ts)|*( 1/m_plus - 1/m_minus )

Con la trasversalita' `lam_m(tf) = 1` (una massa in piu' a burnout e' una massa
in piu' di payload) si ottiene la **forma usata in `validate_staging_corner.m`**:

    eta*[ lam_m(tf) - lam_m(ts) ] = c*|lam_v(ts)|*(1/m_plus - 1/m_minus)

Lettura fisica del bilancio (ritardare lo staging di `dts`):
- **guadagno** `eta*Q*(1 - lam_m(ts))*dts`: se stadi piu' tardi, il primo stadio
  ha bruciato piu' propellente, quindi la sua struttura e' piu' pesante e viene
  buttata; di conseguenza la struttura del *secondo* stadio, `ms2 = eta*Q*(tf-ts)`,
  e' piu' leggera e resta piu' payload;
- **perdita** `T*|lam_v|*(1/m_plus - 1/m_minus)*dts`: hai volato un `dts` in piu'
  con il veicolo *pesante* invece che con quello *leggero*.

L'ottimo e' dove i due si pareggiano. `main_task4.m` **non risolve questa
equazione**: la cerca facendo `max()` su una griglia. Il file
`validate_staging_corner.m` la risolve.

**Nota sulle due formulazioni equivalenti.** L'appendice del report
(`HM1/report/chapters/AppendixStaging.tex`) deriva la stessa identica condizione
partendo dal payoff `phi = tf` (**minimo tempo di combustione**), per cui
`lam_m(tf) = dphi/dm = 0` e la giunzione non ha termine `dphi/dts`. Il risultato
finale coincide, perche' le due formulazioni differiscono solo per un **gauge**:
`lam_m` entra in `H` solo tramite `-Q*lam_m`, e `lam_m_dot = (T/m^2)*|lam_v|` non
dipende da `lam_m` stesso. Le condizioni necessarie sono quindi invarianti sotto

    lam -> k*lam   (k > 0, scala globale)
    lam_m -> lam_m + b   (offset additivo del solo lam_m)

Passare dal payoff "payload" (`lam_m(tf)=1`) a quello "min tf" (`lam_m(tf)=0`) e'
esattamente `k = (1+eta)*Q` e `b = 1`. La corner condition, scritta come
**differenza** `lam_m(tf) - lam_m(ts)`, e' invariante sotto entrambi: la differenza
uccide `b`, e i due membri sono omogenei di grado 1 in `k`. Ed e' proprio per
questo che si puo' valutarla sui costati "sporchi" che il BVP integra con
`lam_m0 = 1` e `H(0) = 0`.

> **Possibile domanda d'esame** -- allo staging i costati saltano?
> *Risposta:* no, sono tutti continui, perche' la massa jettisonata
> `ms1 = eta*Q*ts` dipende solo dal tempo di giunzione e non dallo stato: la
> Jacobiana della mappa di salto rispetto a `x(ts-)` e' l'identita', e la
> condizione di corner `lam(ts-)' = lam(ts+)' * dDelta/dx-` collassa in
> `lam(ts-) = lam(ts+)`. Quello che *salta* e' l'Hamiltoniana, di
> `T*|lam_v|*(1/m+ - 1/m-) > 0`, perche' il salto di massa dipende esplicitamente
> da `ts`.

---

## `%% Sweep staging time` (righe 38-112)

```matlab
ts_max = min(tf_single*0.95, 1/(Q*(1+eta)) - 0.01);
ts_vec = linspace(0.01, ts_max, 50);
...
[~, idx_mid] = min(abs(ts_vec - tf_single*0.4));
for pass = 1:2
```

- Righe 41-45: il range ammissibile di `ts`. Due vincoli: `ts < tf` (approssimato
  con `0.95*tf_single`, perche' `tf` del caso a due stadi non e' ancora noto) e
  `m_plus > 0`, cioe' `ts < 1/(Q*(1+eta)) = 0.4545` con i dati del problema, con
  un margine di `0.01`. Con `Q=2, eta=0.1` vince quasi sempre il primo vincolo.
- Riga 46: **50 punti di griglia**. E' il limite di risoluzione dell'intero Task 4:
  il passo e' circa `0.008` in `ts`, quindi il `ts_opt` riportato e' **agganciato
  alla griglia** e non e' il vero punto stazionario. Non c'e' nessun raffinamento
  (`fminbnd`, interpolazione parabolica del massimo, golden section) -- e' un
  limite reale dello script, non un difetto nascosto.
- Righe 48-51: pre-allocazione a `NaN`. I punti che non convergono restano `NaN`
  e vengono filtrati alla riga 115 con `valid = ~isnan(pay_two)`. E' il pattern
  giusto: un punto fallito non sporca il massimo.
- Righe 57-64: **continuation a due passate**. Si parte da `idx_mid` (il punto piu'
  vicino a `0.4*tf_single`), si sale fino in fondo (pass 1), poi si torna indietro
  da `idx_mid-1` a 1 (pass 2). Ogni solve usa come guess la soluzione del `ts`
  precedente (`z_prev_loc`, righe 74-78 e 109). Questa e' la ricetta standard della
  repo per far convergere lo shooting su uno sweep: il BVP indiretto ha un bacino
  di convergenza stretto, e partire dal vicino gia' risolto e' quasi sempre
  sufficiente.
- Riga 66: `z_prev_loc = []` **dentro** il loop `pass`. Conseguenza: la seconda
  passata riparte dal warm start monostadio `z_ref` invece che dalla soluzione
  gia' convergente a `idx_mid` (che e' in `sol_two{idx_mid}`). Funziona, ma e' un
  warm start sprecato -- la soluzione a `ts` adiacente sarebbe stata un guess
  migliore.
- Riga 109: `z_prev_loc` viene aggiornato **solo in caso di successo** (`ef > 0`).
  Corretto: se un punto fallisce, il vicino successivo riparte comunque dall'ultima
  soluzione buona invece che da spazzatura.
- Righe 88-98: dopo la convergenza la traiettoria viene **ri-integrata a mano**:
  arco 1 su `[0, ts]`, poi il salto `z_s(5) = z_s(5) - ms1` (riga 95), poi arco 2
  su `[ts, tf]`. Osservare che `z_s(6)` (cioe' `lam_m`) **non viene toccato**: e'
  la continuita' dei costati, fatta operativamente.
- Righe 103-106: il payload a due stadi.

      mp2 = Q*(tf - ts)         propellente del secondo stadio
      ms2 = eta*mp2             sua struttura
      mu  = mf - ms2

  Nota che qui **non** si usa la formula monostadio `mf*(1+eta) - eta`: sarebbe
  sbagliata, perche' la struttura del primo stadio e' gia' stata buttata e non
  va sottratta una seconda volta.

### L'identita' nascosta: il payload dipende solo da `tf`

Vale la pena esplicitarla, perche' il codice non la scrive e cambia
completamente la lettura del Task 4. Sostituendo
`mf = 1 - (1+eta)*Q*ts - Q*(tf - ts)` in `mu = mf - eta*Q*(tf - ts)`:

    mu = 1 - (1+eta)*Q*ts - (1+eta)*Q*(tf - ts)
       = 1 - (1+eta)*Q*tf

**Il payload non dipende esplicitamente da `ts`: dipende solo dal tempo di volo
totale `tf`.** (La stessa formula vale nel caso monostadio: `mf*(1+eta) - eta`
con `mf = 1 - Q*tf` da' anch'essa `1 - (1+eta)*Q*tf`.) Quindi:

- massimizzare il payload equivale a **minimizzare `tf`**;
- lo staging aiuta solo perche', alleggerendo il veicolo a `ts`, permette di
  raggiungere le stesse condizioni terminali in **meno tempo di combustione**;
- per `ts` fissato, il BVP interno massimizza `mf`, e siccome
  `mf = m(ts+) - Q*(tf - ts)` con `m(ts+)` gia' determinato da `ts`, massimizzare
  `mf` equivale a minimizzare `tf`, che equivale a massimizzare il payload.
  **L'obiettivo del livello interno e' quindi coerente con quello esterno** -- cosa
  tutt'altro che scontata in uno schema a due livelli.

Verifica numerica con i numeri del README: monostadio `mu = 0.0241` da'
`tf = (1 - 0.0241)/2.2 = 0.4436`; due stadi `mu = 0.0680` da' `tf = 0.4236`. Il
`validate_staging_corner.m` stampa proprio `tf` di sweep `= 0.424` (riga 52). Quindi
il famoso **"+182% di payload" e' in realta' un taglio del 4.5% sul tempo di
combustione**: il payload e' una piccola differenza di numeri grandi e amplifica
tutto.

> **Possibile domanda d'esame** -- il livello interno massimizza `mf`, ma
> l'obiettivo vero e' il payload. Non e' un'incoerenza?
> *Risposta:* no, e si dimostra. Per `ts` fissato, `m(ts+) = 1-(1+eta)*Q*ts` e'
> gia' determinato, quindi `mf` e il payload sono entrambi funzioni affini
> strettamente decrescenti di `tf`: massimizzare l'uno o l'altro equivale a
> minimizzare `tf` e produce la stessa traiettoria. L'unica differenza fra le due
> formulazioni sta nella *normalizzazione dei costati* (in particolare nel valore
> di `lam_m(tf)` e in `H(tf)`), non nella soluzione fisica.

---

## `%% Find optimal staging` + plot (righe 114-191)

- Righe 115-119: `max()` sul vettore `pay_two` filtrato, poi rimappatura
  dell'indice sul vettore completo. E' l'ottimizzatore esterno: una **ricerca
  esaustiva su griglia**, niente di piu'.
- Righe 124-131: report, incluso il guadagno percentuale sul payload. Come detto
  sopra, quel `%` va letto con cautela: e' l'amplificazione di una differenza.
- Righe 134-144: payload vs `ts`, con la retta orizzontale del monostadio. La riga
  144 (`lg.Position(2) = lg.Position(2) + 0.12`) e' un aggiustamento cosmetico per
  non sovrapporre la legenda alla `yline`.
- Righe 156-188: ri-integrazione e plot della traiettoria ottima e del profilo di
  massa. La riga 183 disegna il salto verticale della massa a `ts` come segmento
  tratteggiato -- e' l'unica firma visibile dello staging.
- Nota: nessuna di queste righe verifica la corner condition. Lo script si fida
  del massimo di griglia.

---

## `%% EXPORT FIGURES` (righe 195-211)

- Righe 196-199: cartella `figures/` accanto allo script (creata se manca) e
  `slugify`, che genera i nomi dallo `Name` della figura; il prefisso `task4_`
  viene aggiunto alla riga 210, nella chiamata a `exportgraphics`.
- Righe 204-208: `theme(fig, 'light')` con `try/catch` che ricade su
  `fig.Color = 'w'` per MATLAB pre-R2025a. Serve a non esportare figure con lo
  sfondo scuro del desktop dentro il report LaTeX.

---

## `shooting_single` (righe 215-246)

```matlab
Lam0 = sqrt(lam_vx0^2 + lam_vy0^2);
H0 = -lam_vy0 + p.T*(Lam0 - 1/p.c);
res = [zf(2)-p.yf; zf(3)-1; zf(4); H0];
```

- Riga 215: firma. Incognite `z0 = [lam_vx0; lam_vy0; lam_y; tf]`, 4 residui.
  Chiamata solo da `fsolve` alla riga 24. **Nessun blocco `arguments`** (riga 224):
  scelta di design, sta nel loop interno di `fsolve`/`ode45`.
- Righe 227-229: guardia `tf <= 0 || tf > 2` che restituisce `1e6*ones(4,1)`. E' il
  trucco standard per tenere `fsolve` dentro un bracket fisico, ma va detto: rende
  il residuo **discontinuo**, e se `fsolve` calcola il Jacobiano per differenze
  finite a cavallo del bordo ottiene derivate mostruose. Funziona perche' l'ottimo
  e' lontano dal bordo.
- Riga 234: `ic = [0;0;0;0;1;1]`, con `lam_m0 = 1`.
- Righe 242-243: il residuo di trasversalita'. Deriva dalla forma compatta di `H`
  valutata a `t=0` con `vy=0, m=1, lam_m0=1`.
- Riga 245: i tre residui terminali. `zf(2) = y(tf)`, `zf(3) = vx(tf)`,
  `zf(4) = vy(tf)`: quota target, velocita' orizzontale unitaria (nondim), velocita'
  verticale nulla. `x(tf)` e' **libero** (nessun vincolo di downrange) -- ed e'
  esattamente da li' che viene `lam_x(tf) = 0`, quindi `lam_x = 0` ovunque, quindi
  la struttura dei costati che `ode_burn` da' per scontata.

---

## `shooting_twostage` (righe 248-298)

```matlab
[~,Z1] = ode45(@(t,z) ode_burn(t,z,pp), [0 ts], ic, opts_ode);
z_s = Z1(end,:)';
ms1 = p.eta * p.Q * ts;
z_s(5) = z_s(5) - ms1;
% lam_m is continuous across staging
[~,Z2] = ode45(@(t,z) ode_burn(t,z,pp), [ts tf], z_s, opts_ode);
```

- Riga 248: firma. Stesse 4 incognite di `shooting_single`; `ts` **non** e' incognita,
  arriva come parametro `p.ts` fissato dallo sweep.
- Righe 262-264: guardia, ora anche su `tf <= ts` e `ts <= 0`.
- Righe 266-267: **la stessa `pp` viene usata su entrambi gli archi**. Questa singola
  riga e' la continuita' dei costati `lam_vx`, `lam_vy`, `lam_y`: `ode_burn` li
  ricostruisce dai parametri, quindi le loro funzioni del tempo non hanno alcun
  kink a `ts`. La rampa `lam_vy(t) = lam_vy0 - lam_y*t` prosegue liscia attraverso
  lo staging.
- Righe 271-276: arco 1 su `[0, ts]`.
- Righe 279-281: **il salto**. Solo `z_s(5)` (massa) e' modificato; `z_s(1:4)` (stato
  cinematico) e `z_s(6)` (`lam_m`) passano intatti. Il commento della riga 281 dice
  esplicitamente cio' che abbiamo derivato sopra.
- Righe 282-284: guardia `m_plus > 0`, correttamente **prima** di integrare l'arco 2
  (in `validate_staging_corner.m` la stessa guardia arriva *dopo* la propagazione --
  vedi quella pagina).
- Righe 287-292: arco 2 su `[ts, tf]`, con lo stesso RHS: stesso motore, stesso `Q`.
  L'unica cosa cambiata e' `m`, quindi l'accelerazione `T/m` fa un gradino verso
  l'alto.
- Righe 294-297: **residui identici a quelli del monostadio**. E qui c'e' il punto
  onesto da segnalare: il quarto residuo `H0 = 0` a `t=0` (riga 295) e' ereditato
  pari pari dal Task 1 -- il commento della riga 256 lo enuncia soltanto come
  `H=0 imposed at t0`, senza qualificarlo. Nel monostadio quella riga *e'* davvero
  la trasversalita' per `tf` libero -- `H` e' costante lungo l'arco, quindi `H(0)=0`
  e `H(tf)=0` sono la stessa cosa. **Nel caso a due stadi non lo e' piu'**, perche'
  `H` salta a `ts`. In realta' qui `H0 = 0` funziona come **condizione di gauge**
  (fissa la scala dei costati), non come trasversalita'. E' innocuo, e il motivo e'
  sottile:
  - i 3 residui terminali dipendono dai costati **solo attraverso i rapporti**
    `lam_vy0/lam_vx0` e `lam_y/lam_vx0` (perche' `phi = atan2(...)` e' invariante
    per riscalamento positivo), quindi determinano da soli i due rapporti e `tf`,
    cioe' **tutta la traiettoria**;
  - la quarta incognita e' di fatto la *scala* di `lam_v`, che non tocca la
    traiettoria; `H0 = 0` la fissa e rende il Jacobiano 4x4 non singolare.

  Conseguenza pratica: traiettoria, `tf`, `mf` e payload sono **corretti**, ma i
  `lam_m` che escono da questo script **non** soddisfano `lam_m(tf) = 1` e non sono
  i costati "veri" a meno di un riscalamento e di un offset additivo. Chi vuole
  usarli in una condizione di corner deve usarne una forma **invariante di gauge** --
  ed e' precisamente quello che fa `validate_staging_corner.m` con la differenza
  `lam_m(tf) - lam_m(ts)`.

> **Possibile domanda d'esame** -- nel tuo Task 4 imponi `H(0) = 0`. Ma se
> l'Hamiltoniana salta a `ts`, come puo' essere la condizione di trasversalita' per
> `tf` libero?
> *Risposta:* non lo e'. Nel monostadio `H` e' costante e `H(0)=0` coincide con
> `H(tf)=0`; con lo staging `H` salta di `T*|lam_v|*(1/m+ - 1/m-)`, quindi
> `H(0)=0` e `H(tf)=0` sono condizioni diverse. Nella mia formulazione `H(0)=0`
> serve a fissare la scala dei costati (che altrimenti resterebbe indeterminata per
> omogeneita'). La traiettoria non se ne accorge, perche' dipende solo dalla
> *direzione* del primer vector e i tre residui terminali la determinano gia' da
> soli. Nella normalizzazione "vera" del payoff-payload (`lam_m(tf)=1`) si avrebbe
> invece `H(tf) = eta*Q`, perche' `mu = m(tf) - eta*Q*(tf-ts)` dipende
> esplicitamente da `tf`; nella formulazione equivalente a minimo `tf` usata
> nell'appendice del report si ha `lam_m(tf)=0` e `H(tf)=1`. Tutte e tre sono lo
> stesso problema, in tre gauge diversi.

---

## Limiti e onesta' sul file

1. **Lo sweep e' a 50 punti, senza raffinamento.** `ts_opt` e' un nodo della
   griglia (passo ~0.008); il vero punto stazionario sta fra due nodi. Nessuna
   parabola sul massimo, nessun `fminbnd`.
2. **La corner condition non e' mai verificata qui.** Lo script trova un massimo
   numerico e si ferma. E' `validate_staging_corner.m` a chiudere il cerchio.
3. **Stesso `Q`, stesso `c`, stesso `T` per i due stadi.** Semplificazione forte:
   niente motore da vuoto, niente Isp diverso, niente coefficienti strutturali
   diversi per stadio. Il "+182%" ne e' figlio.
4. **Nessuna salita verticale.** Il veicolo parte da fermo e ruota subito: a `t=0`
   la spinta e' gia' inclinata di ~81 gradi ma non 90. E' la variante richiesta dal
   Task 4, non un errore, ma non e' un profilo di lancio realistico.
5. **Nessun coast allo staging.** Jettison istantaneo e riaccensione immediata.
6. **`z_prev_loc` azzerato a ogni passata** (riga 66): la seconda passata butta via
   un warm start migliore gia' disponibile.
7. **Le guardie a `1e6`** rendono i residui discontinui: se il guess iniziale
   cadesse fuori bracket, `fsolve` avrebbe un Jacobiano privo di senso.

---

## Possibili domande d'esame

**D: Allo stage separation la massa salta. I costati saltano?**
R: No. La condizione di corner generale per una mappa di giunzione
`x(ts+) = g(x(ts-), ts)` e' `lam(ts-)' = lam(ts+)' * dg/dx-`. Qui
`g(x, ts) = x + Delta(ts)` con `Delta = (0,0,0,0,-eta*Q*ts)'`: `Delta` dipende solo
dal *tempo* di giunzione, non dallo stato, quindi `dg/dx- = I` e tutti i costati
sono continui. In `main_task4.m` questo si vede alla riga 280: si sottrae la massa
`z_s(5)` e basta, `z_s(6) = lam_m` resta com'e'. Se invece la massa jettisonata
fosse stata scritta come funzione di `m(ts)`, `lam_m` avrebbe avuto un salto.

**D: L'Hamiltoniana e' continua allo staging?**
R: No, salta verso l'alto di `T*|lam_v(ts)|*(1/m_plus - 1/m_minus) > 0`. Tutti i
termini di `H = lam_y*vy - lam_vy + T*(|lam_v|/m - lam_m/c)` sono continui tranne
`1/m`, che aumenta di colpo. Fisicamente: buttata la struttura, l'accelerazione
`T/m` cresce e il veicolo "guadagna merito" piu' in fretta. L'Hamiltoniana si
conserva solo *dentro* ciascun arco (il sistema e' autonomo), non attraverso la
giunzione.

**D: Qual e' allora la condizione di ottimalita' sull'istante di staging?**
R: Il coefficiente di `dts` nella variazione:
`dphi/dts + H(ts-) - H(ts+) + lam(ts)'*dDelta/dts = 0`. Con
`phi = m(tf) - eta*Q*(tf-ts)` (payload) si ottiene

    eta*[ lam_m(tf) - lam_m(ts) ] = c*|lam_v(ts)|*( 1/m_plus - 1/m_minus )

E' un bilancio fra due effetti opposti del ritardare lo staging: a sinistra il
guadagno (piu' struttura buttata al primo stadio => meno struttura al secondo, piu'
payload), a destra la perdita (hai volato un `dt` in piu' con il veicolo pesante).
`main_task4.m` **non** impone questa equazione: la sostituisce con uno sweep su
`ts` e un `max()`.

**D: Perche' non risolvi direttamente il problema con `ts` come quinta incognita?**
R: Si puo', ed e' esattamente quello che fa `validate_staging_corner.m` (5 incognite,
5 residui, il quinto e' la corner condition). Lo sweep di `main_task4.m` ha due
vantaggi didattici: (a) produce la **curva** payload-vs-`ts`, che mostra dove sta
l'ottimo e quanto e' piatto attorno; (b) e' robusto -- non serve indovinare `ts`,
basta il warm start del monostadio e la continuation. Lo svantaggio e' la
risoluzione: l'ottimo esce agganciato alla griglia.

**D: Il livello interno massimizza `mf`, ma il payload e' `mf - eta*Q*(tf - ts)`.
Non stai ottimizzando la cosa sbagliata?**
R: No. Con `mdot = -Q` costante si ha `mf = 1 - (1+eta)*Q*ts - Q*(tf - ts)`, quindi
sostituendo nel payload: `mu = 1 - (1+eta)*Q*tf`. Il payload dipende **solo da
`tf`**. A `ts` fissato, sia `mf` sia `mu` sono affini strettamente decrescenti in
`tf`: massimizzare l'uno o l'altro da' la **stessa traiettoria**. Le due
formulazioni differiscono solo nella normalizzazione dei costati.

**D: Come mai un guadagno del +182% di payload?**
R: Perche' `mu = 1 - (1+eta)*Q*tf` e' una piccola differenza di numeri grandi.
Passando dal monostadio (`tf ~ 0.444`) ai due stadi (`tf ~ 0.424`) il tempo di
combustione cala del ~4.5%, ma `mu` passa da ~0.024 a ~0.068. Il fattore
amplificatore e' `(1+eta)*Q/mu ~ 32`. La cifra e' inoltre gonfiata dalle
semplificazioni: stesso motore e stesso `eta` per i due stadi, niente atmosfera,
niente drag, niente coast di separazione.

**D: Cosa fissa il vincolo superiore su `ts` alla riga 45?**
R: Due cose. `m_plus = 1 - (1+eta)*Q*ts > 0` da' `ts < 1/(Q*(1+eta))`: oltre quel
punto il veicolo avrebbe buttato piu' massa di quanta ne abbia (propellente bruciato
+ struttura corrispondente supererebbero `m0`). E `ts < tf`, approssimato con
`0.95*tf_single` perche' il `tf` a due stadi non e' ancora noto quando si costruisce
la griglia. Con `Q=2, eta=0.1` il primo bound vale `0.4545` e il secondo, piu'
stringente, ~`0.421`.
