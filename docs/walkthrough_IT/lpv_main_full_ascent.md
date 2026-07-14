# HM3/LTV_FULL_ASCENT/main_full_ascent.m

## Ruolo del file nel progetto

E' lo script principale dello studio "oltre la traccia" (ticket T007). HM3 congela
tutto all'istante di massima pressione dinamica (t_ref = 72 s): matrici del plant,
guadagni PD, notch e TVC sono valutati li' e l'anello chiuso e' LTI. Questo script
solleva quel design all'**intera ascesa (0-140 s)**: il plant rigido di piano di
beccheggio diventa **tempo-variante (LTV / LPV)**, i coefficienti sono letti istante
per istante da `GreensiteLPV_DATA.mat`, e il generatore di vento del professore
(`strong_wind.slx`, che copre tutta l'ascesa) viene messo **direttamente nell'anello**
invece di essere finestrato attorno a max-q come fa `load_wind_profile`.

Lo scopo dichiarato e' rispondere a una domanda che la traccia non pone: *quanto lontano
arriva un design a punto singolo?* Per farlo confronta, sullo stesso plant LTV e con lo
stesso vento, due controllori: quello **frozen** (l'unica coppia PD progettata a 72 s,
tenuta per tutto il volo) e quello **gain-scheduled** (una coppia PD per ogni punto di
una griglia temporale, interpolata in volo). Le metriche sono il picco di assetto
`theta`, il picco di deriva laterale `z`, il picco di comando TVC `delta` e l'indicatore
di carico strutturale `q_bar * alpha`.

Dipendenze: `init_simulink_lpv` (setup dati + schedule + vento, cartella
`LTV_FULL_ASCENT/`), `ode_lpv_ascent` (RHS LTV, 4 stati), e i moduli HM3 riusati
`load_hw3_params`, `build_plant_rigid`, `design_controller`, `assemble_loop`,
`load_wind_profile`, `simulate_gust_response`. Lo script e' la **sorgente di verita'**
numerica: `hm3_full_ascent.slx` (autorato da `build_hm3_full_ascent.m`) deve
riprodurre la stessa risposta, e `run_full_ascent_simulink.m` sovrappone le due.

Nota importante da tenere presente leggendo tutta la pagina: **i numeri stampati oggi
dal codice non coincidono con quelli scritti nel README della cartella** (che promette
picco di deriva 33 m -> 27 m a favore dello scheduling). Ho eseguito lo script e la
conclusione si e' *ribaltata*. La sezione "Cosa produce oggi il codice" spiega esattamente
perche', con la causa individuata nel sorgente. Vale la pena saperlo prima dell'orale.

---

## Intestazione e preambolo (righe 1-20)

- Righe 1-15: banner di documentazione. Dichiara onestamente lo scopo ("Portfolio
  showcase, NOT part of the HM3 deliverable"), i due controllori confrontati e il fatto
  che ode45 e' la sorgente di verita' e Simulink il replay.
- Riga 18: `warning('off', 'Control:analysis:MarginUnstable')`. E' necessario perche'
  l'anello del lanciatore e' **condizionalmente stabile**: `margin`/`allmargin`
  emettono un warning ad ogni valutazione. Spegnerlo qui e' innocuo, ma va ricordato
  che il warning che si sta zittendo e' proprio quello che segnala la natura
  condizionale dell'anello.
- Righe 19-20: aggiunge al path la cartella padre `HM3/` per riusare gli helper.

---

## `%% Setup` -- coefficienti LPV, gain schedule, vento (righe 22-24)

```matlab
S  = init_simulink_lpv();
t0 = S.t0;  Tend = S.Tstop;
```

Una riga sola, ma e' il 90% del lavoro. `init_simulink_lpv` (default `t0 = 5 s`,
`Tstop = 140 s`, `tsched_step = 5 s`) fa quattro cose:

**1. Costruisce i coefficienti efficaci del plant LTV.** Il modello rigido di HM3
(`build_plant_rigid`) e'

    zddot     = a1*zdot + (a1*V + a4)*theta + a3*delta - a1*V*alpha_w
    thetaddot = (A6/V)*zdot + A6*theta      + K1*delta - A6*alpha_w

con `A6` = coefficiente di momento aerodinamico (mu_alpha, l'**instabilita'**: il polo
rigido e' +sqrt(A6)) e `K1` = efficacia di controllo (mu_c). `init_simulink_lpv`
(righe 61-68) raggruppa i prodotti in **sette coefficienti efficaci**, uno per termine:

    c1 = a1        (* zdot)      c5 = A6/V   (* zdot)
    c2 = a1*V + a4 (* theta)     c6 = A6     (* theta e * (-alpha_w))
    c3 = a3        (* delta)     c7 = K1     (* delta)
    c4 = a1*V      (* alpha_w)   invV = 1/V  (alpha_w = v_w * invV)

Il motivo del raggruppamento e' Simulink: cosi' ogni termine dell'equazione e' *una*
lookup 1-D per *un* segnale, e il modello resta ispezionabile blocco per blocco. Ogni
`c_i` diventa un `griddedInterpolant` lineare in `t` (righe 101-104 di
`init_simulink_lpv`), che e' esattamente cio' che `ode_lpv_ascent` valuta ad ogni passo.

**2. Progetta i guadagni frozen `S.K0`** chiamando `design_controller` sul plant
congelato a 72 s con attuatore ideale (Task 1 di HM3). Verificato eseguendo:
`Kp_th = 1.7845`, `Kd_th = 0.4433`, `Kp_z = Kd_z = -1e-3`. I due guadagni di deriva
`Kp_z, Kd_z` sono **costanti e mai schedulati** in tutto lo studio.

**3. Costruisce la gain schedule**: 28 punti su `tsched = 5:5:140 s`, una chiamata a
`design_controller(build_plant_rigid(load_hw3_params('t_ref', tsched(i))))` per punto.

**4. Simula il vento** una volta sola su tutto l'orizzonte (funzione locale
`run_wind_generator`, righe 142-178) e ne ricava `alpha_w(t) = v_w(t)/V(t)` come
interpolante `S.windfun`. Verificato: 1408 campioni su [0, 140] s, picco `|v_w|`
= 22.60 m/s, picco `|alpha_w|` = 2.42 deg.

> **Possibile domanda d'esame** -- perche' i coefficienti sono raggruppati in c1..c7
> invece di tenere separati a1, V, a4?
> *Risposta:* per il solo ode45 sarebbe indifferente. Serve al gemello Simulink: con
> `c2 = a1*V + a4` il termine in theta e' *un* prodotto (una lookup per il segnale
> theta), mentre tenendo a1, V, a4 separati servirebbero due prodotti e una somma in
> piu' per lo stesso termine. Il prezzo e' che i coefficienti fisici non sono piu'
> leggibili singolarmente nel modello: si vede `c2(t)`, non `a1(t)` e `V(t)`.

---

## `%% LTV closed-loop integration` -- frozen vs scheduled (righe 26-35)

```matlab
tt   = (t0:0.02:Tend).';
x0   = zeros(4, 1);              % [z, zdot, theta, thetadot]
odeo = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
[~, xF] = ode45(@(t,x) ode_lpv_ascent(t, x, make_model(S,0)), tt, x0, odeo);
[~, xS] = ode45(@(t,x) ode_lpv_ascent(t, x, make_model(S,1)), tt, x0, odeo);
```

- Riga 27: griglia di uscita a 20 ms. Non e' il passo del solver (ode45 e' a passo
  variabile): e' **solo** la griglia su cui `ode45` interpola l'uscita densa, e non
  influenza ne' i passi interni ne' la frequenza con cui la RHS campiona il vento
  (`ode_lpv_ascent`, riga 14: `aw = M.windfun(t)` viene valutata ai passi che il
  solver sceglie in base a `RelTol`/`AbsTol`, riga 29). Serve al confronto e alla
  sovrapposizione con Simulink. Il vincolo di risolvere le innovazioni di turbolenza
  del generatore (rumore a 0.1 s) vive invece nel gemello Simulink, dove e' imposto
  da `MaxStep = 0.02` (`build_hm3_full_ascent.m`, righe 33-36).
- Riga 28: **stato iniziale nullo a t0 = 5 s**. Approssimazione onesta da segnalare: il
  generatore di vento gira da t = 0, quindi fra 0 e 5 s il veicolo sarebbe gia' stato
  perturbato; qui si riparte da zero. L'errore e' piccolo (il vento e' debole a bassa
  quota) ma esiste.
- Riga 29: tolleranze strette. Non e' pedanteria: la differenza fra i due controllori
  si legge su picchi di deriva di qualche metro, e il confronto con Simulink dichiara
  accordi a ~1e-7 rad, quindi il baseline deve essere molto piu' accurato di cosi'.
- Righe 31-32: **stesso plant, stesso vento, stesso stato iniziale**; cambia solo il flag
  `sched` dentro `make_model`. E' il confronto controllato che rende il risultato
  interpretabile.
- Righe 34-35: `unpack` ricostruisce `delta`, `alpha` e `q_bar*alpha` a posteriori, con i
  guadagni *effettivamente usati* in quella corsa (vedi sezione `unpack`).

Il modello di controllo dentro `ode_lpv_ascent` (riga 22 di quel file) e'

    delta = -(Kp*theta + Kd*thetadot + Kp_z*z + Kd_z*zdot)

cioe' **attuatore ideale** (nessun TVC, nessun ritardo, nessun notch) e `theta_ref = 0`:
e' un problema di *regolazione* dell'assetto attorno alla traiettoria nominale, non di
inseguimento di un pitch program. Coerente con il Task 1 di HM3.

---

## `%% Frozen-time margin sweep` (righe 37-47)

```matlab
tm = (t0:2.5:Tend).';
for i = 1:numel(tm)
    Gi = build_plant_rigid(load_hw3_params('t_ref', tm(i)));
    [gmF(i), pmF(i)] = loop_margin(Gi, S.K0);          % guadagni frozen
    Ksi = S.K0; Ksi.Kp_th = S.fKp(tm(i));
    Ksi.Kd_th = S.fKd(tm(i));
    [gmS(i), pmS(i)] = loop_margin(Gi, Ksi);           % guadagni schedulati
end
```

E' il cuore concettuale dello script. Ad ogni istante `tm(i)` **congela** il plant
(ricostruisce un `ss` LTI con i coefficienti a quell'istante) e ne legge i margini di
guadagno e fase con (a) i guadagni fissi di max-q e (b) i guadagni schedulati.

La lettura "da manuale" del risultato e': il design frozen raggiunge il **minimo** dei
suoi margini **esattamente a max-q** e ha piu' margine ovunque altrove, perche' max-q e'
l'istante piu' instabile (`A6` massimo -> polo +sqrt(A6) massimo). Questo pezzo
**si riproduce**: eseguendo, `|GM|` frozen varia in [6.10, 37.76] dB e `|PM|` in
[30.2, 52.0] deg, con **entrambi i minimi a t = 72.5 s** (il punto della griglia
adiacente a 72 s). E' la giustificazione quantitativa della scelta della traccia.

**Due difetti reali di questo sweep, da conoscere.**

1. `loop_margin` (righe 165-176) usa `margin(L)`, cioe' il margine **di default** di
   MATLAB. Ma HM3 ha una funzione apposta, `classify_margins`, il cui commento dice
   testualmente che le attraversate sotto `w_drift` sono artefatti della deriva laterale
   e che *"Taking margin()'s default instead would pick one of these"*. Lo sweep qui usa
   proprio quel default. Al punto di design i due coincidono (verificato: a 72 s
   `margin()` da' 6.00 dB / 30.0 deg, identico all'aero-GM e al rigid-PM classificati),
   ma lontano da max-q **non coincidono piu'**.
2. Le figure plottano `abs(gm)` e `abs(pm)` (righe 104 e 108). Il valore assoluto
   **nasconde il segno**: un margine di fase negativo (anello instabile) appare come un
   margine grande e positivo. Verificato: a t = 25 s con i guadagni schedulati l'anello
   chiuso e' **instabile** (`isstable(T) = 0`, aero-GM classificato = NaN, PM = -163.8
   deg) ma `margin()` restituisce GM = +10.22 dB e PM = -163.8 deg, che dopo `abs()`
   diventano "10.2 dB / 163.8 deg", cioe' hanno l'aria di margini eccellenti.

> **Possibile domanda d'esame** -- perche' non basta `margin()` su questo anello?
> *Risposta:* perche' l'anello e' condizionalmente stabile e ha un integratore di deriva
> (la posizione z e' l'integrale di zdot): il diagramma di Nichols entra "dall'alto" e
> produce attraversate a bassissima frequenza che non sono margini di corpo rigido.
> `margin()` restituisce il primo/minimo margine che trova, che puo' essere uno di quegli
> artefatti. HM3 risolve con `classify_margins`, che separa le bande (drift / rigido /
> flessibile) e legge un margine per banda; questo script pero' non la usa.

---

## `%% Consistency check at t_ref` (righe 49-61)

```matlab
p72 = load_hw3_params();                       % t_ref = 72 s
[~, T72] = assemble_loop(build_plant_rigid(p72), S.K0, []);
wg72 = load_wind_profile(p72);                 % raffica 1-cos, 12 s
rHM3 = simulate_gust_response(T72, wg72);
Mf = make_model(S, 0);
cc = @(v) griddedInterpolant([0 200], [v v], 'linear', 'nearest');
Mf.fc1 = cc(S.fc1(72)); ...                    % coefficienti congelati
Mf.windfun = griddedInterpolant(wg72.t(:), wg72.alphaw(:), ...);
[~, xf72] = ode45(@(t,x) ode_lpv_ascent(t, x, Mf), wg72.t, x0, odeo);
err_consistency = max(abs(xf72(:,3) - rHM3.theta));
```

E' il test di non-regressione piu' importante della cartella: **se congelo il modello LPV
a 72 s, deve tornare esattamente il modello di HM3**. Il trucco della riga 55 e' elegante:
`cc(v)` costruisce un `griddedInterpolant` su due soli nodi con lo stesso valore, cioe' una
funzione costante; sostituendo `Mf.fc1..fc7` con queste costanti, `ode_lpv_ascent` gira
"in modalita' LTI" senza modificare una riga di `ode_lpv_ascent`.

Le due strade sono radicalmente diverse dal punto di vista numerico -- da un lato `lsim`
su un `ss` LTI costruito da `connect`/`getLoopTransfer`, dall'altro `ode45` sul RHS LTV --
quindi l'accordo prova che il *plant*, la *convenzione dei segni* e il *cablaggio del
controllore* sono gli stessi. Il vento usato e' la **raffica 1-cosine** di default di
`load_wind_profile` (`profile = 'gust'`, ampiezza da `drywind.mat` alla quota di max-q:
verificato `Vg = 6.38 m/s`, `Tg = 3 s`, orizzonte 12 s), non la finestra di strong wind.

Verificato eseguendo: `err_consistency = 1.68e-9 rad` (il README dichiara ~7e-10; stesso
ordine, la differenza dipende dal fatto che `S.K0` non e' piu' lo stesso di quando il
README e' stato scritto -- vedi sotto). In entrambi i casi e' rumore di tolleranza.

---

## `%% Summary` (righe 63-73)

Stampa la tabella dei picchi per i due controllori, l'errore di consistenza e -- riga 71 --
`[qmax, iqm] = max(S.Q)`, cioe' il **massimo della pressione dinamica del dataset**, che
viene stampato accanto al punto di design.

E' un dettaglio che vale oro all'orale: `max(S.Q)` cade a **t = 67 s** (43.9 kPa) mentre
il punto di design della traccia e' **t = 72 s**. Non e' un errore: il punto di design non
e' il massimo di `q_bar`, e' il massimo di `A6` (verificato sui dati: `A6` ha il suo
massimo esattamente a t = 72 s, valore 3.382 1/s^2). Cioe' il "max-q" del titolo e' in
realta' il **massimo dell'instabilita' aerodinamica**, che per questo veicolo e' 5 s dopo
il massimo di pressione dinamica -- perche' `A6 = N_alpha * l_alpha / Iyy` dipende anche da
Mach (tramite C_N_alpha) e dall'inerzia, che cala per consumo di propellente. Tutta la
sezione sullo scheduling in `q` (`main_q_scheduling.m`) nasce da questo sfasamento.

---

## `%% Figures` (righe 75-111) e `%% Export` (righe 113-120)

- f1 (righe 77-86): tre riquadri theta / z / delta, frozen (blu) vs scheduled (rosso),
  con `xline(72)` a marcare il punto di design.
- f2 (righe 89-98): l'indicatore di carico `q_bar * alpha` con la regione di max-q
  evidenziata da un `patch`. E' il *proxy* del carico strutturale (il momento flettente
  e' proporzionale a `q_bar * alpha` a meno di fattori geometrici).
- f3 (righe 101-111): il margin sweep, con le linee di target 6 dB e 30 deg. Attenzione
  agli `abs()` di cui sopra.
- Righe 116-119: `theme(f,'light')` in `try/catch` (la funzione `theme` esiste solo su
  MATLAB recenti) e export PNG a 200 dpi in `figures/`.

---

## `make_model` (righe 123-134)

Impacchetta `S` nella struct piatta che `ode_lpv_ascent` si aspetta. L'unico parametro e'
`sched` (0/1), che diventa `M.sched` logico. Nota che passa **sempre** `fKp`/`fKd` *e* i
guadagni frozen `Kp_th0`/`Kd_th0`: e' `ode_lpv_ascent` (righe 17-21) a scegliere quale
usare. I guadagni di deriva `Kp_z`, `Kd_z` sono sempre quelli di `S.K0`: **la deriva non
e' mai schedulata**.

---

## `unpack` (righe 136-163)

- Righe 146-150: ricostruisce i guadagni usati nella corsa (schedulati o frozen).
- Riga 151: **ricostruisce `delta` a posteriori** dalla stessa legge di controllo del RHS.
  E' corretto solo perche' l'attuatore e' ideale: `delta` e' una funzione algebrica dello
  stato, non ha dinamica propria. Nel caso flessibile (`ode_lpv_flex`) non si potrebbe
  fare, e infatti li' `delta` e' uscita di uno stato TVC.
- Riga 153: `alpha = theta + zdot/V - alpha_w` -- **angolo di attacco totale**, con il
  **segno meno** sul vento (le righe 154-155 sono il commento che lo giustifica in
  sorgente). La derivazione: `theta` e' l'assetto perturbato, `zdot/V` e' l'angolo di
  traiettoria perturbato (per piccole perturbazioni normali alla traiettoria di
  riferimento, con V la velocita' nominale), e `alpha_w = v_w/V` e' la rotazione del
  **vettore velocita' relativa** prodotta dal vento laterale. `alpha` e' l'incidenza
  rispetto alla velocita' **relativa all'aria**, non rispetto al suolo: un vento che
  spinge nella direzione in cui il muso e' gia' inclinato **riduce** l'incidenza vista dal
  veicolo, da cui la sottrazione.
  - **E' la stessa convenzione del RHS che si sta integrando**, e questa e' la ragione
    decisiva: `ode_lpv_ascent.m` (righe 26 e 28) forza con `- M.fc4(t)*aw` e
    `- M.fc6(t)*aw`, cioe' su `alpha = theta + zdot/V - alpha_w`. Un `+` nel
    post-processing descriverebbe un'incidenza **diversa da quella che guida la
    dinamica appena integrata**.
  - *Nota onesta:* fino a poco fa questa riga aveva un `+`. Era un bug **di
    post-processing** -- il plant LPV ha sempre avuto il meno -- quindi `theta`, `z` e
    `delta` non ne erano toccati (e infatti non cambiano); a sbagliare era **solo** il
    canale diagnostico `alpha`, e con esso `qbar*alpha`. Ora e' allineato. L'effetto sul
    numero e' visibile: il picco di `qbar*alpha` frozen passa da 73.3 a **81.3
    kPa*deg**.
- Riga 156: `qa = (Q/1000) * alpha_deg` in kPa*deg. Usa `S.fQ`, cioe' la **Q del dataset**
  (picco 43.9 kPa), non la `p.qbar` di HM3 (circa 81 kPa, da atmosfera esponenziale a
  15.1 km). Le due non coincidono: verificato `Q(72) = 42.4 kPa` dal dataset contro
  `0.5*1.225*exp(-15143/8000)*937.7^2 = 81.1 kPa` dalla formula di `load_hw3_params`. Il
  README lo dichiara. Onestamente: le due grandezze **non sono la stessa cosa** e il
  dataset non e' internamente coerente con un'atmosfera standard (anche il Mach del
  dataset a 72 s, 2.07, con V = 937.7 m/s implicherebbe una velocita' del suono di 453
  m/s, che non e' quella a 15 km). Quindi `q_bar*alpha` qui va letto come **indicatore
  relativo**, buono per confrontare due controllori sullo stesso dataset, non come carico
  in unita' assolute.

---

## `loop_margin` (righe 165-176)

Tre righe: assembla l'anello aperto con `assemble_loop(G, K, [])` (attuatore ideale) e
chiama `margin`. `assemble_loop` rompe l'anello sul segnale `delta` con
`getLoopTransfer(T, 'delta', -1)`, convenzione `1 + L`. Vedi i due difetti discussi nella
sezione sul margin sweep: e' qui che vivono.

---

## Cosa produce oggi il codice (verifica numerica) e divergenze dal README

Ho eseguito il percorso di calcolo dello script (stessa `init_simulink_lpv`, stessi ode45,
stesso sweep) sul commit corrente. Risultati **reali**:

| grandezza | frozen | scheduled |
|---|---|---|
| picco abs(theta) | 0.971 deg | **25.36 deg** |
| picco abs(z) | 29.53 m | **36.19 m** |
| picco abs(delta) | 1.056 deg | 3.44 deg |
| picco q_bar*alpha | 81.3 kPa*deg | 263.9 kPa*deg |
| min GM (sweep) | 6.10 dB (a t = 72.5 s) | 0.80 dB (a t = 120 s) |
| min PM (sweep) | 30.2 deg (a t = 72.5 s) | 4.5 deg (a t = 35 s) |

> I valori di `q_bar*alpha` sono quelli **dopo** la correzione del segno di `alpha`
> (vedi la sezione su `unpack`): col vecchio `+` erano 73.3 e 263.8. `theta`, `z` e
> `delta` **non cambiano** -- il bug era solo nel post-processing, il plant integrato
> aveva gia' il segno giusto -- quindi l'analisi qui sotto sullo scheduling non e'
> toccata.

Il README della cartella dichiara invece "entrambi limitati, picco theta ~0.95 deg" e
"la deriva cala da 33 m a 27 m grazie allo scheduling", e "lo schedule tiene 6 dB / 30 deg
piatti su tutto il volo". **Oggi non e' cosi': lo schedule e' peggiore del design frozen.**

**Causa, individuata nel sorgente.** `init_simulink_lpv` (righe 74-85) costruisce lo
schedule con una *continuation*: passa il risultato precedente come warm start,
`design_controller(..., 'K0', Kprev, ...)`. Ma il `design_controller` attuale
(`HM3/design_controller.m`) documenta esplicitamente `K0  ignored (kept for call
compatibility)` -- cosi' alla riga 19 dell'help, ribadito dall'`arguments` alla riga 32
(`o.K0 (1,2) ... = [0 0]   % accepted, unused`) -- e parte **sempre** (riga 55) dal
punto canonico di D'Antuono

    Kp = 2*A6/K1        Kd = sqrt(A6)/K1

**Quindi la continuation e' un no-op**: il warm start viene calcolato, passato, e buttato.
Verificato con git: al commit in cui il README e' stato scritto (`036cef6`),
`design_controller` usava davvero `x0 = log(o.K0)` e una cost basata su `margin()`; il
refactor successivo (`b43c5e9`) lo ha riscritto nella forma D'Antuono + `classify_margins`,
con `cost = 1e6` quando una banda di margine non ha attraversamento. Il README e le figure
in `figures/` sono **antecedenti** a quel refactor.

Conseguenza fisica: agli istanti in cui l'instabilita' aerodinamica e' piccola (inizio
ascesa, `A6(5 s) = 0.010`; e fine ascesa, `A6(140) = 0.558`) il punto di partenza canonico
`2*A6/K1` **tende a zero**, e la cost restituisce 1e6 costante (nessun attraversamento
aero-GM da agganciare), quindi `fminsearch` non si muove e restituisce il punto iniziale.
Verificato: per `t <= 35 s` e `t >= 120 s` i guadagni schedulati coincidono **esattamente**
con `2*A6/K1` e `sqrt(A6)/K1`; solo nella finestra `t = 40..115 s` la ritaratura riesce e
centra 6.00 dB / 30.0 deg. Con `Kp_th ~ 0.006` a t = 5 s l'anello di assetto e' di fatto
aperto e i guadagni di deriva fissi (`-1e-3`) lo destabilizzano: l'LTI congelato a t = 25 s
con i guadagni schedulati e' **instabile** (`isstable(T) = 0`). Da li' i 25 deg di
escursione di theta.

Come si presenta all'orale: non nascondere il fatto. La lettura corretta e' che lo studio
mette in luce un limite del **tuner**, non dello scheduling in se': un gain schedule che
insegue guadagni proporzionali ad `A6` collassa dove `A6 -> 0`, e serve un vincolo di
banda minima (o mantenere i guadagni dell'ultimo punto valido) fuori dalla regione
aerodinamicamente significativa. La parte che *si riproduce* -- e che e' il punto piu'
importante -- e' che **il design frozen a max-q ha il suo minimo di margine esattamente a
max-q** (6.10 dB / 30.2 deg a t = 72.5 s) e resta adeguato ovunque: la scelta di progetto
della traccia e' quantitativamente giustificata.

---

## Il generatore di vento del professore: come e' cablato e perche' usarlo

Cablaggio nel percorso ode45 (funzione locale `run_wind_generator` in
`init_simulink_lpv`, righe 142-178):

1. `load_system` di `General/hw3-v3/strong_wind.slx`, con `onCleanup` che chiude il
   modello **senza salvare** (riga 154): il file del professore non viene mai modificato.
2. Il logging viene attivato sulle **porte di uscita** del Subsystem (righe 157-160), non
   sulle linee -- e' l'unico modo che funziona quando le linee non hanno nome.
3. `Simulink.SimulationInput` inietta le variabili che il generatore si aspetta
   (`drywind`, `GreensiteLPV`) e fissa lo `StopTime` (righe 162-167).
4. Le due uscite sono il **profilo medio di vento** `v_wp` (schedulato in quota) e la
   **turbolenza di Dryden** (con sigma schedulata in quota, letta da `drywind.mat`).
   Vengono sommate (righe 175-176) su una griglia unione, con `unique` per proteggersi
   dai tempi ripetuti dei log a passo variabile.
5. `alpha_w(t) = v_w(t) / V(t)` con V dal dataset (riga 90), e diventa
   `S.windfun`, un `griddedInterpolant` che `ode_lpv_ascent` valuta ad ogni passo.

Nel `.slx` (`build_hm3_full_ascent.m`, righe 57-64) il generatore e' **copiato una volta
sola** come blocco (`add_block('strong_wind/Subsystem', ...)`) e alimentato dal Clock;
`alpha_w` e' formato in-linea come `(v_wp + turb) * invV(t)`. Stessa formula, stessa quota,
stesso seme.

**Perche' il generatore del professore e non un vento sintetico.** Tre motivi solidi:
(a) il profilo medio e' funzione della **quota**, quindi ha senso solo se il veicolo vola
per davvero da 0 a 140 s -- una raffica 1-cosine di 3 s non puo' rappresentare il
**plateau** di vento sostenuto nella regione di max-q, che e' proprio la sollecitazione
dimensionante; (b) i semi sono fissi, quindi il confronto frozen/scheduled e'
**riproducibile** e le differenze sono attribuibili al controllore e non al vento;
(c) e' lo stesso disturbo su cui e' costruito il caso "strong wind" di HM3, quindi i due
studi sono confrontabili (il picco di theta frozen, 0.97 deg, e' in effetti coerente con
il caso strong-wind finestrato di HM3).

Limite onesto: e' **una sola realizzazione** di un processo stocastico. Nessuna
affermazione statistica e' possibile da questa corsa; per quello serve il Monte Carlo
(che HM3 ha, in `main_montecarlo.m`, ma sul modello congelato).

---

## La fallacia del frozen-time e come questo script la mette alla prova

La progettazione a punti congelati poggia su un'ipotesi che **non e' un teorema**: che se
ogni sistema LTI congelato `A(t_i)` e' stabile, allora il sistema tempo-variante
`xdot = A(t) x` e' stabile. **E' falso in generale**: esistono matrici `A(t)` i cui
autovalori stanno sempre in `Re = -1` e le cui soluzioni divergono (il controesempio
classico e' `A(t)` con autovettori che ruotano abbastanza in fretta). Vale invece un
risultato di **variazione lenta**: se `||dA/dt||` e' sufficientemente piccola rispetto al
margine di stabilita' dei congelati, la stabilita' congelata implica quella LTV.

Questo script mette alla prova l'ipotesi nel modo empirico giusto: **calcola i margini
congelati** (righe 37-47) *e* **integra il vero LTV** (righe 26-35), poi confronta. Ed e'
proprio la separazione di scale a salvare il caso rigido: la banda dell'anello e'
dell'ordine di `sqrt(A6) ~ 1.8 rad/s` a max-q (periodo ~3.4 s), mentre i coefficienti
variano su decine di secondi -- circa un ordine di grandezza di separazione. Coerentemente,
il controllore frozen, che e' congelato-stabile ovunque con margine >= 6 dB, produce una
risposta LTV limitata (picco theta 0.97 deg).

Il caso in cui la fallacia si vede davvero e' nel gemello flessibile (`main_flex.m`): li'
il notch fisso rende l'anello **congelato-instabile da t = 75 s**, eppure la risposta LTV
resta piccola ancora per una decina di secondi (a t = 85 s la coordinata di bending vale
3.9e-4, come nel caso stabile) e diverge solo dopo ~90-100 s. Instabilita' congelata non
significa divergenza immediata: significa **crescita esponenziale a partire da quel
punto**, e su un orizzonte finito la crescita puo' non fare in tempo a manifestarsi.

---

## Possibili domande d'esame

**D: Perche' progettare a max-q e non, per esempio, a meta' volo?**
R: Perche' max-q e' l'istante **dimensionante**. Il polo instabile del corpo rigido e'
`+sqrt(A6)` e `A6` (il momento aerodinamico destabilizzante) ha il suo massimo li': il
plant e' il piu' difficile da stabilizzare e il vento produce il momento maggiore. Lo
sweep di margini congelati di questo script lo dimostra: con i guadagni fissi di max-q, i
margini toccano il **minimo esattamente a t = 72.5 s** (6.10 dB, 30.2 deg -- cioe' i target
di progetto) e sono piu' ampi ovunque altrove. Un design fatto altrove non garantirebbe
nulla a max-q; un design fatto a max-q e' automaticamente conservativo altrove.

**D: Nel dataset la pressione dinamica ha il massimo a 67 s, ma il punto di design e' 72 s.
Contraddizione?**
R: No, e' un punto sottile. Il design point e' il massimo di `A6`, non di `q_bar`.
`A6 = N_alpha * l_alpha / Iyy` e `N_alpha ~ q_bar * S * C_N_alpha(Mach)`: oltre a `q_bar`
ci sono la dipendenza da Mach di `C_N_alpha` e la caduta dell'inerzia `Iyy` per consumo di
propellente. Verificato sui dati: `A6/q_bar` **non e' costante** (0.044 a t = 40 s, 0.080 a
t = 72 s, 0.145 a t = 105 s). Il massimo di `A6` cade quindi 5 s dopo il massimo di
`q_bar`. Questo sfasamento e' esattamente cio' che rompe lo scheduling in `q`.

**D: Che cosa dimostra il consistency check a 72 s e perche' non e' banale?**
R: Dimostra che il modello LPV, congelato al punto di design, **coincide con il modello di
HM3**. Non e' banale perche' le due strade sono numericamente indipendenti: da un lato
`lsim` su una realizzazione `ss` ottenuta da `connect` + `getLoopTransfer`, dall'altro
`ode45` sul RHS scritto a mano con interpolanti. Un errore di segno su `alpha_w` **dentro
il plant**, uno scambio di stati, un `A6/V` scritto come `A6*V` -- si vedrebbero subito.
L'errore misurato e' ~1.7e-9 rad, cioe' rumore di tolleranza dell'integratore.

**Ma va detto cosa il check NON copre**, perche' e' istruttivo: confronta solo `theta`.
Un errore di segno nel **post-processing** di `alpha` (che non entra nella dinamica, ma
solo nel calcolo diagnostico a valle) passerebbe indisturbato -- ed e' esattamente quello
che e' successo: `unpack` ha avuto per un po' `alpha = theta + zdot/V + alpha_w` invece del
meno, e il consistency check restava verde perche' `theta` era ed e' corretto in entrambe
le strade. La lezione: un check di consistenza **verifica cio' che confronta**, e un canale
diagnostico che non retroagisce sulla dinamica non e' coperto da nessun confronto di stato.
Per coprirlo servirebbe confrontare anche `alpha` fra le due strade, o -- meglio -- ricavare
`alpha` dalla dinamica stessa (`theta_ddot - K1*delta = A6*alpha`) invece di riscriverne la
formula a mano.

**D: Lo scheduling nel tuo codice peggiora la risposta. Perche'?**
R: E' un difetto del *tuner*, non del principio. `design_controller` riparte, ad ogni punto
della griglia, dal valore canonico `Kp = 2*A6/K1`, `Kd = sqrt(A6)/K1`, e la sua cost vale
1e6 (costante) quando l'anello congelato non ha l'attraversamento aero-GM da agganciare.
Dove `A6 -> 0` (t <= 35 s e t >= 120 s) succedono entrambe le cose: i guadagni canonici
tendono a zero e `fminsearch` non ha gradiente da seguire, quindi resta li'. Con
`Kp_th ~ 0.006` l'anello di assetto e' praticamente aperto mentre i guadagni di deriva
restano fissi a -1e-3, e l'LTI congelato risulta **instabile** (verificato a t = 25 s). La
correzione naturale e' saturare lo schedule ai guadagni dell'ultimo punto in cui il tuner
converge, oppure imporre una banda minima. Da segnalare anche che la *continuation*
dichiarata (warm start `K0`) e' oggi un **no-op**: `design_controller` accetta ma ignora
`K0`.

**D: Che cos'e' `q_bar * alpha` e perche' lo si guarda?**
R: E' l'indicatore di **carico strutturale**: il momento flettente sul corpo del lanciatore
e' proporzionale alla forza normale, che scala con `q_bar * alpha`, dove alpha e' l'angolo
di attacco totale

    alpha = theta + zdot/V - alpha_w

cioe' l'incidenza rispetto alla velocita' **relativa all'aria** (il meno lo impone il
plant: `Bw = [0; -a1*V; 0; -A6]`). Un controllore che tiene l'assetto rigidamente
(attitude hold) **non** si limita a subire `alpha ~ -alpha_w`: per opporsi al momento
aerodinamico becca il muso **dentro** il vento relativo, e il termine d'assetto si
**somma** a quello di vento. Misurato su HM3 a max-q: picco `|alpha|` = 0.577 deg contro i
0.390 deg del vento da solo -- **l'attitude hold e' load-aggravating**. Un controllore di
**load relief** lascerebbe invece che il veicolo si orienti nel vento riducendo alpha, a
costo di deriva; ma i guadagni di deriva della traccia (-1e-3) sono troppo piccoli per
farlo davvero. E' un trade-off, e questo script lo mostra apertamente invece di
nasconderlo.

**D: Perche' il vento entra come `alpha_w = v_w / V` e non come forza?**
R: Perche' il modello e' **linearizzato**: la forza normale e' `N = N_alpha * alpha` e un
vento laterale `v_w` produce, al primo ordine, una rotazione del vettore velocita'
relativa pari a `alpha_w = v_w / V` (per piccoli angoli, con V la velocita' relativa
nominale). Da qui i termini `-a1*V*alpha_w` e `-A6*alpha_w` nelle due equazioni: sono le
stesse colonne di `alpha_w` di `build_plant_rigid` (`Bw = [0; -a1*V; 0; -A6]`), e le
stesse dei `- M.fc4(t)*aw` / `- M.fc6(t)*aw` di `ode_lpv_ascent` (righe 26 e 28). Il segno
meno viene dal fatto che l'incidenza **aerodinamica** e'
`alpha = theta + zdot/V - alpha_w`: raccogliendo, le due righe forzano su `a1*V*alpha` e
`A6*alpha`, e la colonna di vento eredita il meno. In altre parole: il vento non aggiunge
incidenza al veicolo, **ruota l'aria sotto di lui**. Un limite del modello: `V` e' quella
**nominale**, non la velocita' perturbata.
