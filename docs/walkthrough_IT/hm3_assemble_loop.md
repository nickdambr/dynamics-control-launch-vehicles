# HM3/assemble_loop.m

## Ruolo del file nel progetto

`assemble_loop.m` e' il **nodo centrale** di HM3: prende il plant (rigido o
completo), la struct dei guadagni e la catena attuatore/filtro, e restituisce i
due oggetti su cui vive tutto il resto dell'homework:

- **`L`** -- il **loop aperto SISO**, rotto in corrispondenza di `delta`, nella
  convenzione `1 + L`. E' l'oggetto su cui si legge la carta di Nichols e su cui
  girano `allmargin` / `classify_margins` / `margin`.
- **`T`** -- il **loop chiuso MIMO**, da `{alpha_w, theta_ref}` a
  `{theta, z, zdot, delta}`. E' l'oggetto che si simula (`lsim` in
  `simulate_gust_response`) e su cui si testa la stabilita' (`isstable`,
  `pole`).

E' chiamato praticamente da ovunque: `main_task1.m` (indirettamente via
`design_controller`), `main_task2.m` (righe 26, 31, 59, 92, 123, 144, 155, 175),
`main_task3.m`, `main_montecarlo.m`, e -- criticamente -- **dentro la funzione di
costo del tuner** (`design_controller.m`, riga 82), quindi viene rieseguito
centinaia di volte per ogni run di `fminsearch`.

La legge di controllo che implementa e' la Eq. dell'assignment:

    delta_cmd = Kp_th*(theta_ref - theta_m) - Kd_th*thetadot_m
                - Kp_z*z_m - Kd_z*zdot_m

cioe' un **PD di assetto** piu' un **feedback debole di deriva laterale**
(`Kp_z = Kd_z = -1e-3`, fissati in `design_controller.m` righe 28-29).

**Il punto piu' delicato del file -- e probabilmente di tutto HM3 -- e' che quel
feedback di posizione laterale `Kp_z*z_m` porta nel loop un integratore libero.**
Il plant ha `z' = zdot` con la prima colonna della matrice A tutta nulla
(`build_plant_rigid.m`, righe 11-14): `z` non retroaziona su nulla, quindi
**s = 0 e' un autovalore del plant** e la funzione `delta -> z` ha un polo
nell'origine. Comunque piccolo sia `Kp_z`, il ramo `Kp_z*z_m` inietta quel `1/s`
in `L`: **|L| -> infinito per w -> 0**. E' questo che fa "venire la curva
dall'alto" sulla Nichols, che produce **crossing multipli** a bassa frequenza, e
che rende un singolo numero di `margin()` privo di significato. E' la ragione per
cui esiste `classify_margins.m`.

---

## Firma e contratto (righe 1-21)

```matlab
function [L, T, info] = assemble_loop(G, K, Wact)
arguments
    G {mustBeA(G, 'lti')}
    K (1,1) struct
    Wact = tf(1)
end
if isempty(Wact), Wact = tf(1); end    % [] => ideal actuator
```

- **Riga 1**: tre output. `info` (riga 41) e' una struct con i due blocchi IO
  costruiti, nei campi `Kc` e `Wact`. **Nessun chiamante lo richiede**: in tutto
  il repo `assemble_loop` e' sempre invocato con due output al massimo
  (`[L, T]`, `[~, T]` o solo `L`), e `init_simulink_hm3.m` non la chiama affatto
  -- ricostruisce i blocchi da se' (`build_plant_*`, `build_tvc`,
  `build_notch_filter`). E' quindi un **dead output**: utile solo per debug
  interattivo.

- **Righe 15-19**: `arguments`. `G` deve essere un `lti` (accetta `ss`, `tf`,
  `zpk`, `genss`). `K` una struct scalare -- **nessun controllo sui campi**: se
  manca `Kd_z` la funzione crasha alla riga 25 con un errore poco parlante. E'
  una fragilita' nota e accettabile (chiamata sempre con struct costruite in
  casa).

- **Riga 21**: `if isempty(Wact), Wact = tf(1); end`. Alias documentato: `[]` e
  argomento assente significano **attuatore ideale**. E' cosi' che
  `main_task1.m` ottiene il caso "Task 1": `design_controller(G, [])`. Con
  `Wact = 1` il comando del PD **e' istantaneamente** la deflessione fisica
  (`delta == u_pd`), che e' l'idealizzazione di Task 1.

## Il controllore come guadagno statico (righe 23-27)

```matlab
Kc = ss([K.Kp_th, -K.Kp_th, -K.Kd_th, -K.Kp_z, -K.Kd_z]);
Kc.InputName  = {'theta_ref','theta_m','thetadot_m','z_m','zdot_m'};
Kc.OutputName = {'u_pd'};
```

- **Riga 25**: il controllore e' una **riga 1x5** promossa a `ss` (sistema con
  zero stati, solo matrice D). Applicata al vettore di ingressi:

      u_pd = [ Kp_th, -Kp_th, -Kd_th, -Kp_z, -Kd_z ] *
             [ theta_ref; theta_m; thetadot_m; z_m; zdot_m ]

           = Kp_th*(theta_ref - theta_m) - Kd_th*thetadot_m
             - Kp_z*z_m - Kd_z*zdot_m

  che e' esattamente la docstring (riga 3) e la legge della traccia. Il primo e
  il secondo coefficiente sono `+Kp_th` e `-Kp_th`: l'**errore** di assetto.

  **Attenzione al segno effettivo del ramo di deriva.** I guadagni sono
  `Kp_z = Kd_z = -1e-3` (negativi), quindi le entrate della riga valgono
  `-Kp_z = +1e-3` e `-Kd_z = +1e-3`: sul canale `z_m`/`zdot_m` il feedback entra
  in `u_pd` con **segno positivo**. Il doppio segno negativo (gain negativo x
  segno negativo nella legge) e' esattamente cio' che produce l'azione di
  **load relief** descritta nel README: il veicolo, invece di tenere l'assetto
  rigidamente, si lascia beccheggiare **verso** il vento e riduce l'incidenza
  aerodinamica, pagando in deriva laterale. E' la contrapposizione classica
  attitude-hold vs load-relief/drift-minimum.

  **Nessuna azione integrale.** Il PD e' puro. Non e' una dimenticanza: un
  integratore nel controllore aggiungerebbe -90 gradi di fase a bassa frequenza,
  proprio dove l'anello ha bisogno di **guadagno minimo** per stabilizzare il polo
  aerodinamico RHP (l'aero gain margin). Sommato all'integratore gia' presente
  nel plant (`z`), darebbe un doppio integratore in DC. E non c'e' errore
  stazionario da azzerare: la raffica e' un transitorio.

- **Righe 26-27**: nomi di IO. Sono **la chiave del cablaggio**: `connect` (riga
  35) non usa indici, usa i **nomi dei segnali**. `theta_m`, `thetadot_m`, `z_m`,
  `zdot_m` sono esattamente i primi quattro `OutputName` del plant
  (`build_plant_rigid.m` riga 33 / `build_plant_full.m` riga 54), quindi
  `connect` li collega da solo.

## La catena attuatore (righe 29-32)

```matlab
Wa = ss(Wact);
Wa.InputName  = {'u_pd'};
Wa.OutputName = {'delta'};
```

- **Righe 30-32**: la catena `Wact` (TVC * Pade * notch, oppure `1`) viene
  inserita **nel ramo diretto**, fra l'uscita del controllore e l'ingresso
  `delta` del plant. Non e' in retroazione, non e' sul sensore: e' un blocco in
  serie. Conseguenza: **tutta la fase che TVC + ritardo + notch tolgono finisce
  dentro `L`**, e quindi dentro i margini.

  Nota: il ramo `u_pd -> delta` chiude il cerchio dei nomi. `Wa.OutputName =
  'delta'` coincide con `G.InputName{1} = 'delta'`, quindi `connect` chiude
  l'anello automaticamente.

## Chiusura dell'anello (righe 34-35)

```matlab
T = connect(G, Kc, Wa, {'alpha_w','theta_ref'}, ...
            {'theta','z','zdot','delta'}, {'delta'});
```

- **Riga 35**: `connect` con quattro gruppi di argomenti.
  1. I **blocchi** (`G`, `Kc`, `Wa`), cablati per nome:
     `G.theta_m -> Kc.theta_m`, ..., `Kc.u_pd -> Wa.u_pd`,
     `Wa.delta -> G.delta`. L'anello e' chiuso.
  2. Gli **ingressi esterni** mantenuti: `{'alpha_w','theta_ref'}`. `alpha_w` e'
     il **disturbo** (entra nel plant sulla seconda colonna di B), `theta_ref` il
     riferimento di assetto.
  3. Le **uscite** mantenute: `{'theta','z','zdot','delta'}`. Nota bene: `theta`,
     `z`, `zdot` sono le uscite di **plotting** del plant (righe 45-47 di
     `build_plant_full.m`), cioe' gli stati **veri**, NON le misure INS
     contaminate `theta_m`/`z_m`. Scelta corretta e importante: **l'anello si
     chiude sulla misura sporca, ma si riporta cio' che il veicolo fa davvero**.
     `delta` e' l'uscita di `Wa`, quindi la deflessione **post-attuatore** (in
     Task 1, con `Wact = 1`, coincide con `u_pd`).
  4. Il **sesto** argomento della chiamata (l'ultimo), `{'delta'}`: sono gli
     **analysis point** (APs).
     `connect` inserisce un punto di rottura *virtuale* sul segnale `delta`: `T`
     non e' piu' un `ss` ma un **`genss`**, con dentro l'informazione di dove
     l'anello puo' essere aperto. **Senza questo, la riga 38 non funzionerebbe.**
     L'AP e' trasparente in simulazione (chiuso al valore nominale), quindi `T`
     resta il vero loop chiuso.

## Loop aperto e riduzione (righe 37-39)

```matlab
L = getLoopTransfer(T, 'delta', -1);
L = minreal(tf(L), 1e-6);
```

- **Riga 38**: `getLoopTransfer(T, 'delta', -1)` apre l'anello all'AP `delta` e
  restituisce la funzione di trasferimento di anello. Il terzo argomento `-1` e'
  il **segno del loop**: chiede la `L` nella convenzione **negativa**, cioe'
  quella per cui il loop chiuso e' `T = L/(1+L)`. E' la convenzione che
  `margin`, `allmargin`, `nichols` e il criterio di Nyquist si aspettano
  (punto critico in -1, ovvero -180 gradi / 0 dB). Senza il `-1` si otterrebbe
  la `L` col segno fisico dell'anello (che, essendo in retroazione negativa,
  vale `-L`) e **tutti i margini uscirebbero ruotati di 180 gradi**.

  Nota tecnica: per un anello SISO singolo la `L` **non dipende** da dove si
  rompe (il prodotto scalare `Kc*Wa*G` e' invariante per rotazione ciclica).
  Rompere in `delta` e' una scelta **fisica**: e' il punto in cui, su un banco,
  si inietterebbe uno stimolo nell'attuatore, ed e' il punto rispetto al quale il
  **delay margin** ha il significato di "quanto ritardo extra sopporta la catena
  di comando".

- **Riga 39**: `minreal(tf(L), 1e-6)`. Due operazioni. `tf(L)` converte il
  `genss` in funzione di trasferimento razionale; `minreal` con tolleranza 1e-6
  cancella coppie polo-zero quasi coincidenti.

  **Nota di onesta'**: nei casi nominali `minreal` **non cancella nulla**.
  Verificato: Task 1 -> `order(L) = 4` = i 4 stati del plant rigido (`Kc` e `Wa`
  sono statici); Task 2 senza notch -> `order(L) = 11` = 6 (plant full) + 2 (servo
  TVC) + 3 (Pade-3); Task 2 con notch -> `order(L) = 13`. Nessuna riduzione. La
  chiamata e' quindi **difensiva**, non necessaria -- utile perche' il tuner la
  ripete centinaia di volte e i modi non-minimi possono comparire con guadagni
  esotici, ma va detto che nel flusso nominale e' un no-op. Va anche detto il
  rischio: `tf()` su un sistema di ordine 13 con poli che spaziano da 0 a ~254
  rad/s (il polo piu' veloce e' quello del Pade-3 sul ritardo di 20 ms)
  e' la rappresentazione **numericamente meno condizionata** (i coefficienti
  del polinomio spaziano molti ordini di grandezza); `ss` o `zpk` sarebbero piu'
  sicuri. Con ordine 13 il problema non morde, ma e' un punto molle noto.

---

## PUNTO CRITICO -- perche' `margin()` da solo non significa niente

Questo e' il cuore concettuale del file, e la domanda su cui l'orale si gioca.

### 1. Il plant ha un integratore libero sulla deriva

Nella matrice A del plant (`build_plant_rigid.m`, righe 11-14) la **prima colonna
e' identicamente nulla**: lo stato `z` non compare in nessuna equazione. E' un
puro integratore di `zdot`. Quindi `s = 0` e' un autovalore. Autovalori del plant
rigido, verificati:

    eig(A) = [ 0 ,  -1.861 ,  +0.0291 ,  +1.8165 ]   [rad/s]

Il polo in **0** e' l'integratore di posizione laterale. Il polo in **+1.8165**
e' il polo aerodinamico instabile (la coppia disaccoppiata sarebbe
`+/- sqrt(A6) = +/- 1.839`; l'accoppiamento con la deriva la sposta leggermente e
produce anche il polo lento **+0.0291**). Quindi il loop aperto ha **due poli nel
semipiano destro**, non uno.

### 2. Quell'integratore entra in `L` attraverso `Kp_z*z_m`

Il ramo `z_m` del controllore (`-Kp_z` nella riga 25) moltiplica una funzione
`delta -> z` che ha un `1/s`. Di conseguenza:

    |L(j*w)| -> infinito   per  w -> 0     (pendenza -20 dB/dec)

Sulla carta di Nichols la curva **entra dall'alto**. E' esattamente cio' che dice
il commento di `classify_margins.m` (righe 13-14): *"the drift position
integrator makes the Nichols come from the top"*.

**E succede comunque, per qualunque `Kp_z` diverso da zero.** `Kp_z = -1e-3` e'
minuscolo, ma non cambia il fatto: cambia solo **a che frequenza** il lobo di
deriva attraversa lo 0 dB, non **se** lo attraversa. Non e' un effetto che si
puo' rendere trascurabile abbassando il guadagno -- si puo' solo spostare.

### 3. Il risultato: crossing multipli, e `allmargin` restituisce vettori

`allmargin(L)` sul loop rigido di Task 1 (guadagni Kp = 1.78, Kd = 0.44):

    GainMargin   :  -6.00 dB              @  0.593 rad/s
    PhaseMargin  : -133.12 , -40.53 , +30.00  gradi
                   @  0.161 ,  0.222 ,   2.455  rad/s

**Tre** crossing di fase. I primi due (0.161 e 0.222 rad/s) sono **artefatti del
lobo di deriva**: `classify_margins.m` li scarta con la maschera
`pf > opts.w_drift` (riga 49), con `w_drift = 0.3*sqrt(A6) = 0.55 rad/s`, e li
riporta separatamente come `mm.drift_w` (riga 62) marcandoli esplicitamente come
"NOT rigid-body margins". Il crossing utile e' il terzo: **PM = 30.00 gradi a
2.455 rad/s**, il crossover rigido. (I due numeri cadono esattamente sui target
di progetto -- 6 dB e 30 gradi -- perche' il `fminsearch` di `design_controller.m`
minimizza proprio `(|aeroGM|-6)^2 + (rigidPM-30)^2`, riga 88, con due gradi di
liberta' per due target.)

Sul loop completo di Task 2 (con TVC + Pade + notch, ordine 13, PD ri-tarato --
il design effettivamente consegnato, `main_task2.m` righe 153-155) la situazione
e' peggiore: **5 crossing di guadagno e 3 di fase**.

    GM dB : -283.43   -6.00    +7.56   +39.18   +87.96
    @ w   :   0.000    0.542   11.110  129.280  846.680

La prima riga, a **w = 0 esatta**, e' l'artefatto numerico dell'integratore --
ed e' precisamente il motivo della maschera `gf > 0` in `classify_margins.m`
(maschere alle righe 44-45, commento a riga 43: *"exclude the DC / integrator
entry at gf == 0"*).

### 4. E c'e' anche la stabilita' condizionale

Con **due poli RHP** in anello aperto (P = 2), il criterio di Nyquist richiede
esattamente **2 avvolgimenti antiorari** del punto critico. Conseguenza: il loop
e' **condizionalmente stabile**, cioe' stabile solo per guadagni in una **fascia**
`[K_min, K_max]`, non in `(0, K_max]`. Se il guadagno **scende** troppo il polo
aerodinamico non viene piu' stabilizzato; se **sale** troppo si destabilizza dal
lato alto.

Questo e' visibile nei segni dei gain margin:

- **Aero GM** = margine di **riduzione** di guadagno, con `gmdb < 0`
  (`classify_margins.m`, riga 44: maschera `gmdb < 0`). Qui **-6.00 dB**: il loop
  tollera al massimo 6 dB **in meno** di guadagno.
- **Rigid GM** = margine di **aumento**, con `gmdb > 0` (riga 45). Nel loop
  rigido di Task 1 **non esiste** (attuatore ideale, nessun rolloff -> nessun
  crossing di guadagno in aumento). Compare in Task 2, a +7.56 dB @ 11.1 rad/s.

Un "gain margin" letto ingenuamente come "quanto guadagno posso aggiungere" e'
quindi **letteralmente al contrario** del numero che conta di piu' in questo
progetto.

### 5. Cosa restituisce `margin()` davvero

Verificato sul loop di Task 1: `margin(L)` restituisce `GM = -6.00 dB @ 0.593`,
`PM = 30.00 deg @ 2.455`. Per **onesta'**: in questo caso specifico il PM che
`margin` sceglie e' proprio quello rigido, quindi non e' "sbagliato" -- ma:

- **nasconde** gli altri due crossing di fase (-133 e -40 gradi), che sono la
  firma del lobo di deriva e della stabilita' condizionale;
- restituisce un GM **negativo** (-6.00 dB), che un lettore distratto
  interpreterebbe come "instabile", mentre e' il margine di **riduzione** ed e'
  esattamente il target di progetto (6 dB);
- sul loop completo di Task 2 **nasconde il Rigid GM** (+7.56 dB @ 11.1 rad/s),
  che e' il margine che il ri-tuning del PD deve preservare.

Cioe': il problema non e' che `margin` restituisce un numero *falso*, e' che
restituisce **un** numero dove ne servono **cinque**, tutti con un significato
fisico diverso, e senza dire a quale banda appartiene. Da qui
`classify_margins.m`, che li **bina per frequenza** (deriva / rigido / flex) e li
riporta con il nome del fenomeno che rappresentano.

> **Possibile domanda d'esame** -- Perche' il loop include il feedback di
> posizione laterale, visto che complica tutta la lettura dei margini? Non potevi
> mettere `Kp_z = 0`?
> *Risposta:* Perche' senza feedback di deriva il veicolo tiene l'assetto e
> incassa **tutta** l'incidenza del vento: a max-q e' la condizione di carico
> peggiore. Il termine `Kp_z*z` (con `Kp_z` negativo) e' l'azione di **load
> relief / drift-minimum**: lascia beccheggiare il veicolo verso il vento
> riducendo `q_bar*alpha`, al prezzo di deriva laterale. Il costo di progetto e'
> quello descritto sopra -- e c'e' anche un costo quantitativo, dichiarato dalla
> docstring di `design_controller.m` (righe 8-10): *"the lateral-drift feedback
> erodes the aerodynamic gain margin (the canonical decoupled 6 dB drops to ~4 dB
> on the full loop)"*. E' proprio per questo che il tuner ri-tara i guadagni **sul
> loop pieno** invece di fidarsi della formula di piazzamento disaccoppiata
> `Kp = 2*A6/K1`, `Kd = sqrt(A6)/K1`.

---

## Possibili domande d'esame

**D: Perche' la funzione restituisce sia `L` sia `T`? Non basta uno dei due?**
R: No, servono entrambi e servono a cose diverse. `L` e' **SISO** e vive nel
dominio della frequenza: e' l'oggetto su cui si legge Nichols, si calcolano
`allmargin`/`classify_margins`, si progettano guadagno e fase. `T` e' **MIMO**
(2 ingressi, 4 uscite) e vive nel tempo: e' cio' che si simula con `lsim` e su
cui si chiama `isstable`/`pole`. Il punto decisivo e' che **`alpha_w` non compare
in `L`**: la raffica e' un disturbo che entra nel plant in un punto diverso dalla
rottura dell'anello, quindi da `L` non si puo' ricavare la risposta alla raffica.
Viceversa da `T` non si leggono i margini. E siccome i due escono dallo **stesso**
`connect`, non possono divergere fra loro (nessun rischio di simulare un modello
diverso da quello su cui si sono letti i margini).

**D: Cos'e' l'analysis point `{'delta'}` nel `connect`, e cosa succede se lo
tolgo?**
R: E' un punto di rottura virtuale che `connect` inserisce sul segnale `delta`.
Trasforma `T` da `ss` a `genss` e memorizza **dove** l'anello puo' essere aperto.
Senza di esso `getLoopTransfer(T,'delta',-1)` (riga 38) non avrebbe un punto a cui
riferirsi e fallirebbe: dopo la chiusura, `delta` sarebbe solo un filo interno.
L'AP e' trasparente in simulazione (viene chiuso al valore nominale), quindi `T`
resta il vero anello chiuso e `lsim` da' la risposta corretta.

**D: Perche' `getLoopTransfer` con il segno `-1`?**
R: Perche' `margin`, `allmargin`, `nichols` e il criterio di Nyquist lavorano
nella convenzione `1 + L`, con punto critico in -1 (equivalentemente -180 gradi
/ 0 dB sulla Nichols). L'anello fisico e' in retroazione negativa, quindi il
guadagno d'anello "come cablato" vale `-L`. Il flag `-1` chiede a MATLAB di
restituire la `L` **gia' nella convenzione positiva**, quella per cui il loop
chiuso e' `L/(1+L)`. Senza il flag tutti i margini uscirebbero sfasati di 180
gradi e i grafici di Nichols sarebbero letti sul punto critico sbagliato.

**D: L'integratore libero della deriva -- da dove viene esattamente, e perche'
non lo puoi eliminare abbassando `Kp_z`?**
R: Viene dal plant, non dal controllore: nella matrice A la prima colonna e' nulla
(`build_plant_rigid.m` righe 11-14), cioe' lo stato `z` non retroaziona su nessuna
equazione -- e' un puro integratore di `zdot`. Quindi `s = 0` e' un autovalore del
plant e la funzione `delta -> z` ha un polo nell'origine. Il ramo `Kp_z*z_m` lo
porta dentro `L`. Abbassare `Kp_z` **non lo elimina**: il termine `Kp_z/s` domina
comunque per `w -> 0`, qualunque sia `Kp_z != 0`. Cambia solo la frequenza a cui
il lobo di deriva attraversa lo 0 dB (piu' `Kp_z` e' piccolo, piu' bassa). L'unico
modo di eliminarlo e' `Kp_z = 0` esatto -- cioe' rinunciare al load relief.

**D: Il loop e' condizionalmente stabile. Cosa significa operativamente sui
margini?**
R: Che il guadagno d'anello ammissibile e' una **fascia** e non una semiretta.
Il plant ha due poli RHP (verificato: +1.8165 e +0.0291 rad/s), quindi Nyquist
richiede 2 avvolgimenti antiorari del punto critico: se il guadagno **scende**
sotto una soglia perdo gli avvolgimenti e il polo aerodinamico non e' piu'
stabilizzato; se **sale** troppo destabilizzo dal lato alto. Operativamente:
esistono **due** gain margin, uno di riduzione (l'**Aero GM**, negativo in dB, il
target da 6 dB dell'assignment) e uno di aumento (il **Rigid GM**, positivo, che
compare solo quando c'e' rolloff, cioe' da Task 2 in poi). Sulla Nichols la curva
si **infila fra i due punti critici** invece di stare tutta da una parte -- e'
proprio la firma grafica della stabilita' condizionale, visibile nelle figure di
Task 1 e Task 2.

**D: Le uscite di `T` sono `theta`, `z`, `zdot` -- perche' non `theta_m`, `z_m`?**
R: Perche' l'anello si chiude sulle **misure** (contaminate dal bending via INS,
`build_plant_full.m` righe 33-36: `theta_m = theta + sigma_ins*eta`) ma si vuole
riportare cio' che il **veicolo fa davvero**. Il plant espone entrambe le famiglie
di uscite apposta: le prime quattro (`theta_m`...`zdot_m`) vanno al controllore, le
ultime tre (`theta`, `z`, `zdot`, righe 45-47) sono gli stati veri per il
plotting. `connect` prende le prime per cablare e le seconde come uscite esterne.
E' la distinzione fra "cosa vede il sensore" e "cosa succede al razzo", ed e'
esattamente il motivo per cui il bending e' pericoloso: **lo vedi anche se non
c'e'**, e il controllore reagisce a un fantasma.

**D: `minreal(tf(L), 1e-6)` -- serve davvero?**
R: Nei casi nominali no. Verificato: `order(L)` vale 4 in Task 1 (= i 4 stati del
plant, dato che `Kc` e `Wa` sono statici), 11 in Task 2 senza notch (6+2+3) e 13
con il notch. **Nessuna cancellazione avviene.** La chiamata e' difensiva:
`assemble_loop` gira dentro `fminsearch` (`design_controller.m` riga 82) centinaia
di volte con guadagni arbitrari, e li' una quasi-cancellazione puo' comparire.
Il rovescio della medaglia e' che `tf()` su un sistema di ordine 13 con poli da 0
a ~254 rad/s e' la rappresentazione numericamente peggio condizionata; `zpk` o `ss`
sarebbero piu' robusti. E' un punto molle noto, non un bug attivo.
