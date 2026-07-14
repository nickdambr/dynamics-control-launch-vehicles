# HM1/main_task3.m

## Ruolo del file nel progetto

`main_task3.m` risolve la **terza variante** della missione di HM1: alla sequenza
del Task 2 (`salita verticale -> arco propulso`) si aggiunge un **arco balistico
di coasting** finale. Il motore viene spento **prima** dell'iniezione, e il
veicolo raggiunge la quota target `yf` in volo libero, arrivando all'apogeo
esattamente con `vy = 0` e `vx = 1`. Se questo e' possibile, si risparmia
propellente: la spinta non deve piu' combattere la gravita' negli ultimi
istanti.

La differenza rispetto al Task 2 e' **strutturale**, non solo di condizione
iniziale: adesso il problema e' genuinamente **bang-bang** in spinta
(motore acceso -> motore spento), e il punto di spegnimento non e' un dato ma
un'**incognita determinata dall'ottimalita'**. La condizione che lo determina e'
l'annullarsi della **funzione di switching** `S = |lam_v|/m - lam_m/c`.
Conseguenza diretta: `lam_m` **entra nei residui** (compare in `S`), quindi il
trucco di normalizzazione del Task 1/2 (`lam_m0 = 1`) non e' piu' comodo, e
`lam_m0` torna a essere un'incognita. Lo shooting passa da **4x4 a 5x5**.

La seconda idea chiave del file e' che **l'arco di coasting non viene mai
integrato numericamente dentro lo shooting**. Con `T = 0` la dinamica e' un moto
parabolico esatto, quindi le tre condizioni terminali del Task 2
(`y(tf)=yf, vx(tf)=1, vy(tf)=0`) vengono **trasportate all'indietro in forma
chiusa** fino all'istante di cutoff, dove diventano due sole condizioni
algebriche (`vx(tc)=1` e `yc + 0.5*vyc^2 = yf`). La terza (`vy(tf)=0`) e'
soddisfatta per costruzione, perche' il tempo di coasting viene definito come
il tempo di salita all'apogeo, `t_coast = vyc`.

Il file dipende da `ode_burn.m` (RHS condiviso stato + `lam_m`) e definisce tre
funzioni locali: `ode_vertical` e `event_altitude` (copiate identiche da
`main_task2.m`) e `shooting3`. Esporta le figure in `HM1/figures/` con prefisso
`task3_`.

---

## Richiamo: la funzione di switching e il problema bang-bang

Nel Task 1/2 la spinta e' sempre accesa e `Q` e' costante. Per far comparire un
coast bisogna ammettere una manetta `u` in `[0,1]`, con `T = u*c*Q` e
`m_dot = -u*Q`. L'Hamiltoniana (convenzione del massimo, `lam_x = 0`, `g = 1`)
diventa

    H = lam_y*vy - lam_vy
        + u*c*Q*( |lam_v|/m - lam_m/c )

perche' il massimo su `phi` di `lam_vx*cos(phi) + lam_vy*sin(phi)` e' `|lam_v|`,
e i termini in `u` si raccolgono. La quantita' tra parentesi e' la **funzione di
switching di Lawden**:

    S = |lam_v|/m - lam_m/c

`H` e' **lineare in `u`**, quindi il massimo si ottiene agli estremi:

    u = 1  se  S > 0      (motore acceso)
    u = 0  se  S < 0      (coast)

e lo **spegnimento avviene dove `S = 0`**. E' questa l'equazione che il codice
impone come quinto residuo (`S` e' definita alla riga 298 e inserita nel vettore
`res` alla riga 304), ed e' l'unica ragione per cui il cutoff non e'
una variabile arbitraria da spazzare: e' fissato dalle condizioni necessarie.

> **Possibile domanda d'esame** -- perche' la spinta ottima e' bang-bang e non
> "a manetta intermedia"?
> *Risposta:* perche' `H` dipende linearmente dalla manetta `u`. Il PMP richiede
> di massimizzare `H` rispetto a `u` in `[0,1]`: il massimo di una funzione lineare
> su un intervallo sta sempre a un estremo, salvo il caso singolare `S == 0` su un
> intervallo di tempo di misura non nulla (arco singolare), che qui non si presenta.
> Quindi `u` vale 1 o 0 e commuta quando `S` cambia segno.

---

## Intestazione e parametri (righe 1-18)

- Righe 8-13: stessi dati del Task 2 -- `c = 0.6`, `eta = 0.1`, `y1 = 1e-4`,
  `yf = 0.04`, `Q = 2`, `T = c*Q = 1.2`.
- Righe 15-18: stesse tolleranze strette (`ode45` a `1e-10/1e-12`, `fsolve` a
  `1e-10`).
- **Attenzione al commento di riga 3**: *"Search for optimal engine cutoff time
  to maximize payload"*. E' **fuorviante**: nel codice non esiste nessuno sweep
  o ricerca sul tempo di cutoff. Il cutoff `t_burn` e' una delle 5 incognite dello
  shooting ed e' determinato **implicitamente** dalla condizione di ottimalita'
  `S(tc) = 0`. E' esattamente il punto di forza del metodo indiretto: la struttura
  ottima esce dalle condizioni necessarie, non da un'ottimizzazione esterna.

---

## Fase 1 -- salita verticale (righe 20-31)

```matlab
ic_vert = [0; 0; 1];
opts_vert = odeset('RelTol', 1e-12, 'AbsTol', 1e-14, ...
    'Events', @(t,z) event_altitude(t, z, y1));
[T_vert, Z_vert] = ode45(@(t,z) ode_vertical(t,z,T,Q), ...
                         [0 1], ic_vert, opts_vert);
t1  = T_vert(end);
vy1 = Z_vert(end,2);
m1  = Z_vert(end,3);
```

- Identica al Task 2 (righe 19-37 di `main_task2.m`), tolleranze comprese.
- Righe 27-29: si estraggono **solo** `t1`, `vy1`, `m1`. La quota raggiunta
  `Z_vert(end,1)` **non** viene letta: alla riga 48 il codice usa il valore
  nominale `p.y0 = y1`. Nel Task 2 invece si usa `y_1 = Z_vert(end,1)`.
  Incoerenza reale tra i due script; l'impatto e' dell'ordine della tolleranza
  dell'evento (trascurabile), ma va dichiarata.
- Stessa fragilita' del Task 2: se `Q < 1/c = 1.667` il veicolo non decolla,
  l'evento non scatta e lo script prosegue silenziosamente con `t1 = 1`.

---

## Fasi 2+3 -- impostazione del BVP burn + coast (righe 33-61)

Il commento alle righe 34-41 e' la specifica del BVP, e vale la pena leggerlo
come mappa:

```matlab
% BVP unknowns: [lam_vx0, lam_vy0, lam_y, lam_m0, t_burn]
%   vx(tc) = 1
%   y(tc) + 0.5*vy(tc)^2 = yf
%   lam_m(tc) = 1
%   lam_vy(tc) = lam_y*vy(tc)
%   S(tc) = |lam_v|/m - lam_m/c = 0
```

- Riga 34: **il vettore delle incognite cresce a 5**. Rispetto al Task 2 si
  aggiunge `lam_m0`. Motivo: `lam_m` compare esplicitamente nella funzione di
  switching, quindi la sua **scala assoluta rispetto a `|lam_v|` conta**. Nel
  Task 1/2 `lam_m` era invece un puro quadratura, disaccoppiato da tutto, e lo si
  poteva normalizzare a 1 a costo zero.
- Righe 43-51: la struct `p` con parametri e stato iniziale dell'arco propulso --
  come nel Task 2, la giunzione salita/burn trasferisce **solo lo stato**, senza
  condizioni sui costati (la fase verticale non ha gradi di liberta').
- Righe 53-61: quattro guess iniziali in una cell array. Tutti hanno
  `lam_m0 < 0` (`-0.5, -0.6, -0.4, -0.5`) e `t_burn` corto (`0.10 - 0.18`).
  I guess servono a due scopi (vedi sotto): entrare nel bacino della radice
  **fisica** e non della radice **spuria**.

### La radice spuria (commento righe 53-57)

La condizione `res(2)`, `yc + 0.5*vyc^2 = yf`, e' **quadratica in `vyc`**:
ammette due rami.

- Ramo **fisico**: `vyc > 0`. Il motore si spegne **sotto** `yf` e il veicolo
  sale balisticamente fino all'apogeo, che cade esattamente su `yf`. Il tempo di
  coast e' `t_coast = vyc > 0`.
- Ramo **spurio**: `vyc < 0`. L'equazione algebrica e' soddisfatta lo stesso (c'e'
  `vyc^2`), ma il "coast" avrebbe durata negativa: `t_coast = vyc < 0`. Il veicolo
  sarebbe gia' oltre l'apogeo e starebbe scendendo. Non e' una soluzione fisica,
  ma **fsolve non lo sa**: e' solo un'altra radice del sistema.

Il codice non aggiunge un vincolo `vyc > 0` (fsolve e' unconstrained): lo gestisce
con un **filtro a posteriori** (righe 65-76). Vedi sotto.

> **Possibile domanda d'esame** -- perche' la condizione di coast e'
> `yc + 0.5*vyc^2 = yf` e cosa rappresenta?
> *Risposta:* durante il coast `T = 0`, quindi `vy_dot = -g = -1` (adimensionale)
> e `vy(s) = vyc - s`: l'apogeo (`vy = 0`) e' raggiunto dopo `s = vyc`. La quota
> guadagnata e' `vyc*s - 0.5*s^2` valutata in `s = vyc`, cioe' `0.5*vyc^2`. Quindi
> la quota all'apogeo e' `yc + 0.5*vyc^2` e deve valere `yf`. E' la conservazione
> dell'energia (`v^2/2g` in dimensionale) e sostituisce **due** condizioni
> terminali del Task 2 (`y(tf)=yf` insieme a `vy(tf)=0`), perche' `vy(tf)=0` e' vera
> per costruzione all'apogeo.

---

## Loop dei guess e filtro sulla radice fisica (righe 63-77)

```matlab
for gg = 1:numel(guess_list)
    [z_try, ~, ef_try] = fsolve(@(z) shooting3(z,p,opts_ode), ...
                                guess_list{gg}, opts_fs);
    if ef_try > 0
        ... reintegra l'arco con z_try ...
        if Zc(end,4) > 1e-6      % vyc > 0: accept the physical coast
            z_sol = z_try; ef = ef_try; break;
        end
    end
end
```

- Righe 65-76: **multi-start, non continuazione**. Si prova ogni guess in
  sequenza; se `fsolve` converge (`ef_try > 0`) si **reintegra** l'arco propulso
  (righe 70-71) solo per leggere `vyc = Zc(end,4)` e si accetta la radice **solo
  se `vyc > 1e-6`** (riga 72). E' il filtro che scarta il ramo spurio.
- Da dichiarare con onesta': **in questo script non c'e' nessuna continuazione**
  (nessun warm start dalla soluzione di un problema vicino). Il Task 1 usa la
  continuazione sullo sweep in `Q`; il Task 2 riusa il guess cablato del Task 1;
  il Task 3 usa una lista di 4 guess tarati a mano. Una continuazione naturale
  qui sarebbe: partire dalla soluzione del Task 2 (coast di durata nulla) e
  "aprire" gradualmente il coast, oppure fare continuazione su `yf`.
- Costo: la reintegrazione della riga 70 e' ridondante (l'informazione era gia'
  dentro `fsolve`), ma e' eseguita al massimo 4 volte, quindi e' irrilevante.
- Riga 77: se nessun guess produce una radice fisica, `ef = -1` e lo script stampa
  l'errore alla riga 198 senza produrre nulla.

---

## Post-processing: ricostruzione del coast in forma chiusa (righe 79-119)

```matlab
t_coast = vyc;               % time to reach vy = 0
x_coast = xc + vxc * t_c;
y_coast = yc + vyc * t_c - 0.5 * t_c.^2;
vx_coast = vxc * ones(size(t_c));
vy_coast = vyc - t_c;
mf_task3 = mc;  % mass doesn't change during coast
```

- Righe 87-88: reintegrazione dell'arco propulso su griglia densa (500 punti) per
  i grafici, con `ic2(6) = z_sol(4) = lam_m0`.
- Righe 91-95: stato al cutoff `(xc, yc, vxc, vyc, mc)`.
- Riga 98: `t_coast = vyc`. In adimensionale `g = 1`, quindi il tempo per portare
  `vy` da `vyc` a zero e' esattamente `vyc`. E' il motivo per cui la condizione
  `vy(tf) = 0` **non compare tra i residui**: e' imposta per costruzione.
- Righe 102-106: il coast e' **analitico** -- moto parabolico con `vx` costante
  (nessuna forza orizzontale) e `y` parabola. Non si chiama `ode45`: cio' elimina
  ogni errore di integrazione dall'arco balistico e dimezza il costo dello
  shooting.
- Riga 108: `mf_task3 = mc` -- durante il coast la massa non cambia (`Q = 0`).
  Quindi **massimizzare `m(tf)` equivale a massimizzare `m(tc)`**: e' la massa al
  cutoff che conta.
- Righe 116-118: i print di verifica confrontano `y_coast(end)`, `vx_coast(end)`,
  `vy_coast(end)` con i target `yf`, `1`, `0`. Sono **identita' per costruzione**
  (discendono direttamente dai residui e dalla definizione di `t_coast`), non un
  test indipendente: se lo shooting e' convergito, non possono che tornare.
- Riga 119: payload `mf*(1+eta) - eta`, stessa formula del Task 1/2.

---

## Plots (righe 121-195)

- Righe 122-133: traiettoria a tre colori -- rosso verticale, blu burn, verde
  tratteggiato coast. Il coast e' visivamente l'arco che si "appiattisce"
  raggiungendo l'apogeo con tangente orizzontale.
- Righe 135-153: inset di zoom sulla salita verticale (identico al Task 2, e per
  la stessa ragione: `y1 = 1e-4` e' invisibile a scala piena).
- Righe 156-167: angoli lungo il burn. `lam_vy_k = z_sol(2) - z_sol(3)*T2(kk)` e'
  ricalcolato **analiticamente** (il costato e' lineare), poi
  `phi = atan2(lam_vy_k, lam_vx0)` e `psi = atan2(vy, vx)`.
- Righe 170 e 180: `psi_coast` viene calcolato alla riga 170 e tracciato alla
  riga 180. Durante il coast si traccia **solo `psi`** -- non c'e' spinta, quindi
  `phi` non e' definito. E' la firma grafica del bang-bang.
- Righe 181-182: `xline` marcano `t1` (fine verticale) e `t1 + t_burn` (engine
  cutoff).
- Righe 188-195: profilo di massa. Tre tratti: discesa lineare nella verticale
  (pendenza `-Q`), discesa lineare nel burn (stessa pendenza), **piatto** nel
  coast. Il plateau finale e' la firma visiva del risparmio di propellente.

---

## Export figure (righe 201-217)

- Identico ai Task 1/2: `theme(fig,'light')` (evita PNG a sfondo nero se il
  desktop MATLAB e' in dark mode), `exportgraphics` a 200 dpi, prefisso `task3_`.

---

## `ode_vertical` (righe 221-230) e `event_altitude` (righe 232-242)

- Copie **identiche** delle omonime in `main_task2.m` (righe 178-199). Nessuna
  fattorizzazione in file condiviso (a differenza di `ode_burn.m`).
- `ode_vertical`: dinamica con `phi = pi/2` sostituito a mano,
  `dz = [vy; T/m - 1; -Q]` (con `g = 1`).
- `event_altitude`: `value = y - y1`, `isterminal = 1`, `direction = 1`
  (attraversamento in salita).

---

## `shooting3` (righe 244-305)

Il cuore del file. Cinque incognite, cinque residui.

- Riga 244: firma `res = shooting3(z0, p, opts_ode)`. Nessun blocco `arguments`
  (dichiarato alla riga 257): gira dentro il loop di `fsolve`.
- Righe 259-263: unpack `[lam_vx0; lam_vy0; lam_y; lam_m0; t_burn]`.
- Righe 265-268: stessa box constraint via penalita' del Task 2
  (`t_burn <= 0 || t_burn > 2` -> `res = 1e6`). Stesso caveat: rende il residuo
  discontinuo sulla frontiera e la Jacobiana per differenze finite ivi calcolata
  e' priva di senso.
- Riga 275: `ic = [x0; y0; vx0; vy0; m0; lam_m0]` -- il sesto stato ora parte
  dall'**incognita** `lam_m0`, non piu' da 1.
- Righe 277-283: `try/catch` attorno a `ode45`.
- Righe 285-292: si legge lo stato al cutoff e si ricostruisce **analiticamente**
  il costato:

      lam_vy_c   = lam_vy0 - lam_y*t_burn
      lam_v_norm = sqrt(lam_vx0^2 + lam_vy_c^2)

  Legittimo: `lam_vy` e' esattamente lineare, quindi non serve integrarlo. Nota
  che `lam_mc` viene invece **letto dall'integrazione** (`zf(6)`), perche' la sua
  evoluzione `lam_m_dot = (T/m^2)*|lam_v|` non ha forma chiusa elementare.

### I cinque residui (righe 300-304)

```matlab
res = [vxc - 1;
       yc + 0.5 * vyc^2 - p.yf;
       lam_mc - 1;
       lam_vy_c - lam_y * vyc;
       S];
```

**`res(1)` -- `vx(tc) = 1`.** Durante il coast non c'e' forza orizzontale
(`T = 0`, gravita' verticale, niente drag), quindi `vx` e' **costante**:
`vx(tf) = vx(tc)`. La condizione terminale `vx(tf) = 1` viene percio' trasportata
all'indietro senza approssimazioni fino al cutoff.

**`res(2)` -- `yc + 0.5*vyc^2 = yf`.** La condizione balistica derivata sopra.
Assorbe in una sola equazione le **due** condizioni terminali `y(tf) = yf` e
`vy(tf) = 0` del Task 2.

**`res(3)` -- `lam_m(tc) = 1`.** E' la **trasversalita' vera**, `lam_m(tf) = 1`
(costo = `m(tf)`), riportata al cutoff. Legittimo perche' durante il coast
`T = 0` e quindi

    lam_m_dot = (T/m^2)*|lam_v| = 0   ->   lam_m costante nel coast

quindi `lam_m(tc) = lam_m(tf) = 1`. Questa e' anche la condizione che **fissa la
scala** dei costati (nel Task 1/2 la scala era fissata dalla normalizzazione
`lam_m0 = 1`).

**`res(4)` -- `lam_vy(tc) = lam_y*vyc`.** E' la condizione `H = 0` **scritta
sull'arco di coast**. Derivazione: con `T = 0` e `lam_x = 0`,

    H_coast = lam_y*vy + lam_vy*(-g) = lam_y*vy - lam_vy

(gli altri termini contengono `T` o `Q` e si annullano). Il tempo finale e'
libero e il sistema e' autonomo, quindi `H = 0`; valutando in `tc`:

    lam_y*vyc - lam_vy(tc) = 0

Si verifica facilmente che `H_coast` e' costante lungo il coast:
`d/dt (lam_y*vy - lam_vy) = lam_y*(-1) - (-lam_y) = 0`. Coerente.

**`res(5)` -- `S(tc) = 0`.** La condizione di switching. E' l'equazione che
**determina la durata del burn** `t_burn`: senza di essa il cutoff sarebbe
arbitrario.

### Perche' non c'e' un residuo `H(0) = 0` come nei Task 1/2?

Perche' e' **implicato** da `res(4)` e `res(5)` insieme. Valutiamo `H` sull'arco
propulso all'istante di cutoff (`u = 1`, `T = c*Q`, `lam_x = 0`, `g = 1`):

    H_burn(tc) = lam_y*vyc + (T/mc)*|lam_v_c| - lam_vy_c - lam_mc*Q

Raccogliendo (con `T = c*Q`):

    H_burn(tc) = ( lam_y*vyc - lam_vy_c )
                 + c*Q*( |lam_v_c|/mc - lam_mc/c )
               = res(4) + c*Q * S(tc)

Se `res(4) = 0` e `S(tc) = 0`, allora `H_burn(tc) = 0`; e poiche' il sistema e'
autonomo, `H` e' costante sull'arco propulso, quindi `H = 0` **su tutto il burn**,
inclusi `t = 0`. La condizione `H(0) = 0` che nei Task 1/2 era un residuo esplicito
qui e' **automaticamente soddisfatta**. Questo e' anche l'enunciato della
condizione angolare di Weierstrass-Erdmann: `H` e' **continua** attraverso il
corner, e la continuita' e' garantita proprio da `S(tc) = 0`.

### Le condizioni di corner (Weierstrass-Erdmann) usate implicitamente

Alla giunzione burn -> coast **non c'e' salto di stato ne' vincolo interno**,
quindi:

- i costati sono **continui**: `lam(tc-) = lam(tc+)`. Il codice lo usa senza dirlo:
  `lam_vy_c` calcolato **sul burn** viene inserito nella condizione `H_coast = 0`
  (`res(4)`), che vive **sul coast**; e `lam_mc` calcolato sul burn viene
  eguagliato al suo valore terminale sul coast (`res(3)`).
- `H` e' **continua**: garantita da `S(tc) = 0`, come mostrato sopra.

Non c'e' quindi nessuna condizione aggiuntiva da imporre: la continuita' e'
"cablata" nel modo in cui i residui sono scritti.

### L'hack numerico sulla funzione di switching (righe 294-298)

```matlab
% Using the 1/c form (not the literal lam_mc/c) keeps fsolve in the
% basin of the physical vyc>0 root; the lam_mc/c form drifts to the
% spurious vyc<0 branch. The two coincide at cutoff.
S = lam_v_norm / mc - 1 / p.c;
```

Questo va dichiarato per intero, perche' e' una **deviazione consapevole dalla
forma letterale**:

- La funzione di switching e' `S = |lam_v|/m - lam_m/c`. Alla soluzione,
  `res(3)` impone `lam_mc = 1`, quindi `lam_mc/c = 1/c`: **le due forme
  coincidono nella radice**, e i due sistemi 5x5 hanno **esattamente lo stesso
  insieme di soluzioni**.
- Ma **fuori** dalla soluzione (cioe' durante le iterazioni di Newton) sono
  diverse: la forma `1/c` **rimuove la dipendenza di `res(5)` da `lam_m0`**, cioe'
  azzera una colonna della Jacobiana. Questo cambia il campo di Newton e, secondo
  il commento dell'autore, mantiene `fsolve` nel bacino della radice fisica
  (`vyc > 0`) invece di farlo scivolare su quella spuria.
- Onesta': **e' una scelta di condizionamento numerico, non teorica**. La forma
  teoricamente corretta e' `lam_mc/c`. Le due sono equivalenti alla convergenza;
  se all'orale viene chiesta la funzione di switching, la risposta e'
  `|lam_v|/m - lam_m/c`.

> **Possibile domanda d'esame** -- perche' nel Task 3 non si puo' usare la
> normalizzazione `lam_m0 = 1` del Task 1/2?
> *Risposta:* nei Task 1/2 `lam_m` e' completamente disaccoppiato: non entra ne'
> nella dinamica dello stato, ne' nella legge del controllo, ne' nei residui, e
> puo' quindi essere scalato liberamente (le condizioni necessarie sono omogenee
> di grado 1 nei costati). Nel Task 3 `lam_m` entra nella funzione di switching
> `S = |lam_v|/m - lam_m/c`, che confronta la sua **scala assoluta** con quella di
> `|lam_v|`: fissarlo arbitrariamente a 1 all'istante iniziale falserebbe il punto
> di switching. La scala viene percio' fissata dove la teoria la fissa davvero,
> cioe' con la trasversalita' `lam_m(tf) = 1`, riportata al cutoff come `res(3)`;
> `lam_m0` diventa la quinta incognita. Coerentemente, i guess usati nel codice
> hanno `lam_m0 < 0`: `lam_m` cresce lungo il burn
> (`lam_m_dot = (T/m^2)*|lam_v| > 0`) e deve arrivare a 1 al cutoff, quindi parte
> piu' in basso. Il codice non stampa mai il `lam_m0` convergito, quindi dal solo
> sorgente non se ne puo' confermare il segno alla soluzione.

---

## Limiti noti / punti onesti da dichiarare

- **Il commento di riga 3 e' fuorviante**: non c'e' nessuna "ricerca" del tempo di
  cutoff; `t_burn` e' un'incognita dello shooting fissata da `S(tc) = 0`.
- **`S` e' imposta solo al cutoff**, non verificata sull'arco. Il codice **non
  controlla** che `S > 0` per tutto il burn (cioe' che la struttura ottima sia
  davvero `burn -> coast` e non, per esempio, `coast -> burn -> coast`), ne' che
  `S` attraversi lo zero **in discesa** in `tc`. E' un controllo di ottimalita'
  della struttura che manca.
- **La forma di `S` usata (`1/c` invece di `lam_mc/c`)** e' un espediente di
  convergenza (vedi sopra): corretta alla soluzione, non alla lettera durante le
  iterazioni.
- **Nessuna continuazione**: solo 4 guess cablati + filtro sul segno di `vyc`.
  Se nessuno converge, lo script si arrende.
- **Incoerenza con il Task 2**: qui `p.y0 = y1` (nominale) invece di
  `Z_vert(end,1)` (raggiunta).
- **Nessuna guardia sull'evento** della salita verticale (come nel Task 2).
- **Codice duplicato**: `ode_vertical` e `event_altitude` copiate dal Task 2.
- Il modello resta senza drag e a Terra piatta: un coast "vero" a 620 m -
  240 km di quota adimensionale attraverserebbe l'atmosfera, ma qui non c'e'.

---

## Possibili domande d'esame

**D: Qual e' esattamente la variante del Task 3 rispetto al Task 2, letta dal
codice?**
R: Si aggiunge un arco di coasting finale. In termini di BVP: (1) il vettore delle
incognite passa da 4 a 5, aggiungendo `lam_m0`; (2) i tre residui terminali del
Task 2 (`y(tf)-yf`, `vx(tf)-1`, `vy(tf)`) vengono trasportati analiticamente al
cutoff e diventano due (`vxc - 1` e `yc + 0.5*vyc^2 - yf`), perche' `vy(tf) = 0`
e' vera per costruzione all'apogeo; (3) il residuo `H(0) = 0` sparisce e viene
sostituito da due condizioni nuove: `H_coast = 0` scritta al cutoff
(`lam_vy_c = lam_y*vyc`) e la condizione di switching `S(tc) = 0`; (4) si aggiunge
la trasversalita' `lam_m(tc) = 1`. Totale: 5 incognite, 5 residui.

**D: Come si ricava la funzione di switching e cosa impone?**
R: Si ammette una manetta `u` in `[0,1]` con `T = u*c*Q` e `m_dot = -u*Q`.
Raccogliendo i termini in `u` nell'Hamiltoniana si ottiene
`H = lam_y*vy - lam_vy + u*c*Q*S` con `S = |lam_v|/m - lam_m/c`. Poiche' `H` e'
lineare in `u`, il PMP da' `u = 1` se `S > 0` e `u = 0` se `S < 0`; il passaggio
burn -> coast avviene dove `S = 0`. Imporre `S(tc) = 0` e' l'equazione che
determina la durata del burn: senza di essa il cutoff sarebbe arbitrario.

**D: Perche' l'arco di coasting non viene integrato numericamente?**
R: Perche' con `T = 0` la dinamica e' esatta in forma chiusa: `vx` costante,
`vy(s) = vyc - s`, `y(s) = yc + vyc*s - 0.5*s^2`, `m` costante, `lam_m` costante,
`lam_vy` lineare. Tutte le condizioni terminali si possono quindi riportare
all'indietro fino al cutoff in modo esatto. Questo (a) elimina l'errore di
integrazione sull'arco balistico, (b) elimina il tempo di coast dal vettore delle
incognite (`t_coast = vyc` per definizione di apogeo), (c) dimezza il costo di ogni
valutazione del residuo dentro `fsolve`.

**D: Quali condizioni di giunzione (corner conditions) valgono tra burn e coast, e
dove sono nel codice?**
R: Non c'e' salto di stato ne' vincolo di interior-point, quindi valgono le
condizioni di Weierstrass-Erdmann nella forma "liscia": **costati continui**
`lam(tc-) = lam(tc+)` e **`H` continua**. Il codice le usa implicitamente: `lam_vy_c`
e `lam_mc`, calcolati sull'**arco propulso**, vengono usati in condizioni che
vivono sul **coast** (`res(4)` e `res(3)`) -- cioe' si sta assumendo la continuita'.
La continuita' di `H` e' invece garantita esattamente da `S(tc) = 0`: si dimostra
che `H_burn(tc) = (lam_y*vyc - lam_vy_c) + c*Q*S(tc)`, quindi con `res(4) = 0` e
`res(5) = 0` si ha `H_burn(tc) = H_coast(tc) = 0`.

**D: Perche' esiste una radice spuria e come la scarta il codice?**
R: Perche' il residuo `yc + 0.5*vyc^2 = yf` e' quadratico in `vyc` e ammette anche
`vyc < 0` (veicolo gia' in discesa dopo l'apogeo). Formalmente e' una radice del
sistema, ma corrisponde a `t_coast = vyc < 0`, cioe' un coast di durata negativa:
e' non fisica. `fsolve` non ha vincoli e ci puo' cadere. Il codice non impone
`vyc > 0` come vincolo: prova quattro guess con burn corto (`t_burn = 0.10-0.18`,
che cadono nel bacino della radice fisica), reintegra l'arco e **accetta la prima
radice con `vyc > 1e-6`** (riga 72). E' un filtro a posteriori, non un vincolo.

**D: Il Task 3 usa la continuazione?**
R: No. Usa un **multi-start** su una lista di 4 guess cablati, con un filtro sul
segno di `vyc`. La continuazione (warm start a catena da una soluzione vicina) e'
usata solo nel Task 1, per lo sweep su `Q`. Una continuazione ragionevole per il
Task 3 sarebbe partire dalla soluzione del Task 2 (coast di durata nulla) e
aumentare gradualmente il parametro che apre il coast, oppure fare continuazione
su `yf`.

**D: Perche' il coasting fa risparmiare propellente?**
R: Perche' la massa non cambia durante il coast (`Q = 0`), quindi
`m(tf) = m(tc)`: massimizzare la massa finale equivale a massimizzare la massa al
cutoff. Negli ultimi istanti di un burn "fino all'iniezione" (Task 2) la spinta
deve anche sostenere il veicolo contro gravita' mentre porta `vy` a zero: e' una
perdita gravitazionale. Se invece si spegne prima, e' la **gravita' stessa** a
portare `vy` a zero gratuitamente durante la salita balistica, mentre `vx` e' gia'
al valore richiesto e resta tale. La condizione `S(tc) = 0` individua esattamente
il punto in cui questo scambio diventa conveniente.

**D: Perche' `t_coast = vyc` senza costanti?**
R: Perche' il problema e' adimensionalizzato con `a_rif = g`, quindi `g = 1` in
unita' adimensionali. Durante il coast `vy_dot = -g = -1`, cioe' `vy(s) = vyc - s`,
e l'apogeo (`vy = 0`) e' raggiunto in `s = vyc`. In dimensionale sarebbe
`t_coast = vyc/g`.
