# HM3/LTV_FULL_ASCENT/build_hm3_full_ascent.m

## Ruolo del file nel progetto

Questo file **e' il modello Simulink**. Non e' uno script che carica un `.slx`
gia' esistente e nemmeno un generatore di griglie di plant o di array di `ss`:
e' un *authoring script* che, usando le API `new_system` / `add_block` /
`add_line` / `set_param` / `save_system`, costruisce da zero il file
`HM3/LTV_FULL_ASCENT/hm3_full_ascent.slx` blocco per blocco e filo per filo. Il
`.slx` e' un artefatto **derivato**; la sorgente versionabile e' questo `.m`.

Il modello che costruisce e' la controparte Simulink della baseline MATLAB
`ode_lpv_ascent.m` + `main_full_ascent.m`: il plant rigido di beccheggio di HM3
(`build_plant_rigid`), che nella traccia e' congelato a max-q (t = 72 s), qui
diventa **tempo-variante (LPV)** sull'intera ascesa 0-140 s, e il generatore di
vento del professore (`General/hw3-v3/strong_wind.slx`) viene messo
**dentro l'anello** invece che pre-campionato su una finestra di 12 s attorno a
max-q. Il punto della dimostrazione e' confrontare, sulla stessa tela, i guadagni
**frozen** (l'unico design a max-q tenuto per tutto il volo) contro i guadagni
**gain-scheduled** in tempo.

Dipendenze: nel modello non c'e' un solo numero di **dati o di design**. Tabelle,
breakpoint, guadagni e perfino `StopTime` sono **nomi di variabili** del base
workspace (`lpv_t`, `lpv_c1`...`lpv_c7`, `lpv_invV`, `Kp_th0`, `tsched`,
`Kp_sched`, `sched`, `Tstop`, ...) che `init_simulink_lpv.m` deve avere gia'
pushato. Gli unici numeri letterali cablati nel `.slx` sono costanti
*strutturali*, non parametri del problema: `StartTime = 0` e `MaxStep = 0.02` del
solver (riga 36), `Threshold = 0.5` degli Switch (righe 111-112),
`InitialCondition = 0` degli integratori (righe 73-74, 81-82), `Inputs = 2` dei
Product e `NumberOfTableDimensions = 1` delle lookup. Chi lo chiama:
`run_full_ascent_simulink.m` (con `o.rebuild=true`, oppure automaticamente se il
`.slx` non esiste). Chi consuma il suo output: sempre
`run_full_ascent_simulink.m`, che simula il modello e sovrappone il risultato
alla replica `ode45`.

Le equazioni realizzate sono quelle scritte nel README della cartella e nella RHS
`ode_lpv_ascent.m`:

    zddot     = a1*zdot + (a1*V + a4)*theta + a3*delta - a1*V*alpha_w
    thetaddot = (A6/V)*zdot + A6*theta      + K1*delta - A6*alpha_w

con tutti i coefficienti funzione del tempo di volo.

---

## Firma, header e `arguments` (righe 1-18)

```matlab
function build_hm3_full_ascent(o)
arguments
    o.open (1,1) logical = false
end
```

- Righe 1 e 16-18: la firma prende il solo struct name-value `o`; il blocco
  `arguments` dichiara l'unica opzione, `o.open (1,1) logical = false` (riga 17),
  che decide se lasciare il modello aperto a video dopo la build. Nessun
  argomento di dati: **tutti i dati stanno nel base workspace**, non passano da
  qui.
- Righe 2-14: l'header dichiara il contratto: "Simulink mirror of
  ODE_LPV_ASCENT / MAIN_FULL_ASCENT (source of truth). Run INIT_SIMULINK_LPV
  first so the referenced base variables exist." La gerarchia e' esplicita: la
  verita' e' la RHS MATLAB, il `.slx` la specchia. Questo va detto all'orale
  cosi' com'e': il modello Simulink **non** e' usato per validare la teoria, e'
  usato per mostrare che la stessa matematica si puo' realizzare a blocchi.

---

## Percorsi, generatore del professore, modello nuovo e solver (righe 20-37)

```matlab
load_system(fullfile(gdir, 'strong_wind.slx'));
cleanupSW = onCleanup(@() close_system('strong_wind', 0));
```

- Righe 20-24: `mdl = 'hm3_full_ascent'`; `gdir` punta a
  `HM3/General/hw3-v3`, la cartella del materiale del professore. `addpath(here)`
  e `addpath(hm3)` sono effetti collaterali sul path MATLAB **mai ripristinati**:
  igiene discutibile, ma innocua in una sessione interattiva.
- Righe 27-28: carica `strong_wind.slx` (il generatore di vento del professore)
  solo per poterne **copiare** il sottosistema. L'`onCleanup` lo chiude con il
  flag `0` = *discard changes*: qualunque cosa succeda, anche in caso di errore,
  il file del professore non viene mai salvato. E' la traduzione in codice della
  regola "non si tocca il materiale del corso".
- Riga 31: `close_system(mdl, 0)` se il modello era gia' caricato, poi riga 32
  `new_system(mdl)`. La build e' **distruttiva e idempotente**: si riparte
  sempre da zero, non si fa merge con un modello esistente. Questa e' la ragione
  per cui non ha senso modificare a mano il `.slx`: la modifica sparisce alla
  prossima esecuzione.
- Righe 35-37: solver `ode45` a passo variabile, `StartTime = 0`,
  `StopTime = 'Tstop'` (una **stringa**: e' il nome di una variabile base, non un
  numero), `MaxStep = 0.02`.
- Righe 33-34, il commento: "MaxStep bounded to resolve the generator's 0.1 s
  noise innovations". Il generatore di turbolenza Dryden aggiorna il rumore ogni
  0.1 s; con `MaxStep = 0.02` si hanno 5 passi per innovazione, altrimenti il
  solver potrebbe scavalcare le discontinuita' del rumore e la replica `ode45`
  non riuscirebbe piu' a ricostruire lo stesso ingresso.

> **Possibile domanda d'esame** -- il modello raggiunge ~1e-7 rad di accordo con
> `ode45`: e' merito delle tolleranze?
> *Risposta:* No, ed e' un punto onesto da ammettere. Il builder **non imposta
> mai** `RelTol` / `AbsTol` sul modello, quindi Simulink gira con il suo default
> (`RelTol = 1e-3`), mentre la replica in `run_full_ascent_simulink.m` (riga 29)
> usa `odeset('RelTol',1e-9,'AbsTol',1e-11)`. Le due integrazioni hanno
> tolleranze diversissime; l'accordo si regge su `MaxStep = 0.02`, che inchioda
> il passo di Simulink a un valore cosi' piccolo da rendere l'errore locale molto
> piu' fine di quanto la tolleranza richiederebbe. Se si rilassa `MaxStep`,
> l'accordo degrada.

---

## Le scorciatoie di authoring: `A`, `W`, `LK` (righe 39-42)

```matlab
A  = @(name, src, pos, varargin) add_block(src, ...
        [mdl '/' name], 'Position', pos, varargin{:});
W  = @(s, d) add_line(mdl, s, d, 'autorouting', 'on');
LK = 'simulink/Lookup Tables/n-D Lookup Table';
```

- Riga 40: `A` = "aggiungi blocco". `src` e' il path di libreria
  (`built-in/Product`, `built-in/Integrator`, ...), `pos` il rettangolo
  `[x1 y1 x2 y2]` sul canvas, `varargin` le coppie parametro/valore.
- Riga 41: `W` = "collega". La sintassi `'blocco/porta'` (es. `'c1/1'`) e' la
  forma testuale di `add_line`; `autorouting` fa instradare i fili a Simulink.
- Riga 42: la libreria del blocco lookup. Si usa la **n-D Lookup Table**
  forzandone la dimensione a 1 (`NumberOfTableDimensions = '1'`) invece del
  vecchio blocco `Lookup Table` 1-D, che e' deprecato.

Nota onesta: buona parte del file e' **layout** (le coordinate `[430 60 460 90]`
e simili sono numeri magici scelti a occhio). Non hanno alcun effetto sulla
matematica, ma occupano molte righe e rendono il diff rumoroso.

---

## Clock e le lookup dei coefficienti effettivi (righe 44-54)

```matlab
cs = {'c1','c2','c3','c4','c5','c6','c7','invV'};
tb = {'lpv_c1','lpv_c2', ... ,'lpv_c7','lpv_invV'};
for k = 1:numel(cs)
    A(cs{k}, LK, [...], ...
      'NumberOfTableDimensions','1', ...
      'BreakpointsForDimension1','lpv_t', 'Table',tb{k});
    W('Clock/1', [cs{k} '/1']);
end
```

- Riga 45: un unico blocco `Clock`. E' **il** parametro di scheduling: il modello
  e' un LPV schedulato sul tempo di volo, e il Clock e' l'unica sorgente di
  quel parametro.
- Righe 48-54: otto lookup 1-D, tutte con gli **stessi breakpoint** `lpv_t`,
  ognuna con la propria tabella. Il Clock alimenta tutte e otto. Interpolazione:
  il parametro `InterpMethod` non viene toccato, quindi resta il default
  **`Linear`**.

### Il meccanismo LPV: perche' non un blocco State-Space

Un blocco `State-Space` di Simulink valuta `A, B, C, D` **una volta**, in fase di
compilazione: sono costanti. Non esiste modo di farle variare nel tempo. Per un
plant LPV le strade sono tre:

1. **Blocchi elementari** (la scelta di questo file): una lookup per
   coefficiente, `Product` per moltiplicarla per il segnale, `Sum` per sommare i
   termini, `Integrator` per integrare. Ogni coefficiente resta un segnale
   visibile e sondabile a video.
2. **MATLAB Function block**: si riscriverebbe `ode_lpv_ascent` dentro un blocco.
   Funziona, ma il modello diventa una scatola nera: dal canvas non si vede piu'
   nessun coefficiente, e non c'e' piu' niente da "mostrare".
3. **LPV System block** (Control System Toolbox): richiede un array di `ss`
   costruito su una griglia di scheduling. Va detto onestamente che, a parita' di
   griglia e di interpolazione lineare, sarebbe **matematicamente equivalente**
   alla soluzione 1 (interpolare linearmente ogni entrata della matrice `A` = 
   interpolare linearmente ogni coefficiente). Le ragioni per cui e' scartato
   sono pratiche: dipendenza da un toolbox, opacita', e un array di 141 modelli
   `ss` 4x4 solo per rappresentare 7 scalari.

Il codice non motiva mai esplicitamente la scelta; il README della cartella dice
solo "so every coefficient stays inspectable, mirroring the style of the
professor's own generator". Il resto della motivazione va ricostruito, ed e'
corretto dirlo all'orale come ricostruzione.

### Il trucco dei coefficienti effettivi (il cuore della pagina)

`init_simulink_lpv.m` (righe 61-68) **non** tabula i fattori grezzi `a1, V, a4,
A6, K1`. Tabula i **coefficienti gia' combinati**, uno per termine
dell'equazione:

    c1 = a1          (moltiplica zdot)
    c2 = a1*V + a4   (moltiplica theta)
    c3 = a3          (moltiplica delta)
    c4 = a1*V        (moltiplica alpha_w)
    c5 = A6/V        (moltiplica zdot)
    c6 = A6          (moltiplica theta e -alpha_w)
    c7 = K1          (moltiplica delta)
    invV = 1/V       (alpha_w = v_w * invV)

Perche' e' essenziale e non solo elegante: **l'interpolazione lineare non commuta
con la moltiplicazione**. Sia `L{.}` l'interpolatore lineare su una cella di
ampiezza `h`, con `u = s/h` in `[0,1]`, `Df = f1 - f0`, `Dg = g1 - g0`. Allora

    L{f*g}(u) - L{f}(u)*L{g}(u) = Df*Dg * u*(1-u)

che vale zero **solo** sui breakpoint e raggiunge il massimo `Df*Dg/4` a meta'
cella. E' un errore di ordine `h^2`, ma e' un errore di **modello**, non di
solver: non si riduce stringendo il passo di integrazione.

Conseguenza operativa: il `.slx` e la RHS `ode_lpv_ascent.m` integrano davvero la
**stessa** funzione **soltanto se entrambi interpolano le stesse tabelle
combinate**. Ed e' esattamente quello che accade nel track rigido:

- il `.slx` legge `lpv_c1`...`lpv_c7` (righe 48-49 di questo file);
- `init_simulink_lpv.m` (righe 102-103) costruisce `S.fc1 = gi(c1)` ...
  `S.fc7 = gi(c7)` con `griddedInterpolant(tg, y, 'linear', 'nearest')` sulle
  **stesse identiche** array;
- `ode_lpv_ascent.m` (righe 26-28) usa `M.fc1(t)` ... `M.fc7(t)` e non ricombina
  mai nulla.

Le due catene interpolano gli stessi numeri con la stessa regola. Questo, e non
la fortuna, e' cio' che rende possibile l'accordo a ~1e-7 rad. (Nel modello
**flessibile** questa disciplina viene rotta: vedi la pagina
`lpv_build_full_ascent_flex.md`.)

### Numeri reali della griglia

Verificato leggendo `GreensiteLPV_DATA.mat`: il tempo del dataset e' uniforme,
`t = 0:1:150 s` (151 punti), quindi `lpv_t = 0:1:140` -> **141 breakpoint,
passo h = 1 s**. E' una griglia grossolana: l'errore di non-commutazione, se ci
si cadesse dentro, sarebbe tutt'altro che trascurabile.

> **Possibile domanda d'esame** -- se avessi tabulato `a1`, `V`, `a4`
> separatamente e ricostruito `c2 = a1*V + a4` a valle con un `Product` e un
> `Sum`, il modello sarebbe stato sbagliato?
> *Risposta:* No, sarebbe stato un modello LPV perfettamente legittimo -- ma un
> modello **diverso**, che integra una funzione diversa fra i breakpoint. Il
> confronto con la baseline `ode45` avrebbe avuto un pavimento di errore fissato
> dal termine `Df*Dg*u*(1-u)`, indipendente dalle tolleranze. Su questa griglia a
> 1 s, `L{a1*V+a4} - (L{a1}*L{V} + L{a4})` vale fino a 2.2e-3 in valore assoluto
> (circa 0.003% del coefficiente): abbastanza da mascherare completamente il
> livello 1e-7 che si vuole dimostrare.

### Estrapolazione: una asimmetria che vale la pena conoscere

Le otto lookup dei coefficienti (righe 51-52) **non** impostano `ExtrapMethod`,
quindi restano al default `Linear`. La `griddedInterpolant` della baseline usa
invece `'nearest'` come metodo di estrapolazione (= clip, tiene il valore
d'estremo). Le due regole **non coincidono**. In pratica non morde, perche'
`lpv_t` copre `[0, 140]` e la simulazione gira su `[0, 140]`: non si estrapola
mai. Ma e' una incoerenza latente. Nota che l'autore la conosce, perche' per le
lookup dei guadagni (righe 108-109) mette `ExtrapMethod = 'Clip'` esplicitamente,
e nel modello flessibile lo mette su **tutte** le lookup.

---

## Il generatore di vento e alpha_w (righe 56-64)

```matlab
add_block('strong_wind/Subsystem', [mdl '/WindGen'], ...);
A('vw_sum','built-in/Sum',[...],'Inputs','++', ...);
A('aw_prod','built-in/Product',[...],'Inputs','2');
```

- Riga 57: il sottosistema del professore viene **copiato** dentro il nuovo
  modello. E' una copia, non un link a libreria: se il professore cambiasse
  `strong_wind.slx`, il `.slx` generato non se ne accorgerebbe finche' non si
  rilancia il builder. Il sorgente resta immutato grazie all'`onCleanup` di riga
  28.
- Riga 60: il Clock alimenta l'ingresso del generatore. Vento e veicolo
  condividono **lo stesso orologio**: e' proprio questa la differenza rispetto a
  HM3, dove `load_wind_profile` ritagliava 12 s di vento attorno a max-q.
- Riga 58 + righe 61-62: `vw_sum` (`Inputs = '++'`) somma le due uscite del
  generatore (profilo medio `v_wp` e turbolenza Dryden), coerentemente con
  `init_simulink_lpv/run_wind_generator`, che legge `Outport(1) = sw_vwp` e
  `Outport(2) = sw_turb`.
- Riga 59 + righe 63-64: `aw_prod = vw_sum * invV`, cioe'

      alpha_w(t) = ( v_wp(t) + turbolenza(t) ) / V(t)

  L'angolo d'attacco indotto dal vento e' il rapporto fra velocita' laterale del
  vento e velocita' del veicolo: e' la linearizzazione `alpha_w ~ tan(alpha_w)`,
  valida perche' `v_w << V` (V va da 410 a 3116 m/s nel dataset).

Nota su `invV`: `init_simulink_lpv` calcola `Vsafe = max(V,1)` con il commento
"V(0)=0: guard A6/V and 1/V at lift-off". **Il commento e' falso per questo
dataset**: `min(V) = 410.4 m/s`, quindi `max(V,1)` e' bit-per-bit uguale a `V` e
la guardia non si attiva mai. E' codice morto piu' un commento che contraddice i
dati. Innocuo, ma da non raccontare come se fosse una protezione attiva.

---

## Il plant: prodotti -> somme -> integratori (righe 66-100)

```matlab
A('zdd','built-in/Sum',[...],'Inputs','+++-', ...);
A('int_zd','built-in/Integrator',[...],'InitialCondition','0');
A('int_z', 'built-in/Integrator',[...],'InitialCondition','0');
```

- Righe 68-71: `P1..P4`, quattro `Product` a due ingressi. Ognuno realizza un
  termine: `c1*zdot`, `c2*theta`, `c3*delta`, `c4*alpha_w`.
- Riga 72: `zdd` e' un `Sum` con `Inputs = '+++-'`. **L'ordine delle porte e'
  semantico**: `+P1 +P2 +P3 -P4`, cioe'

      zddot = c1*zdot + c2*theta + c3*delta - c4*alpha_w

  Il segno meno su `alpha_w` e' fisico: il vento riduce l'angolo d'attacco totale
  visto dal corpo rispetto a quello dovuto all'assetto.
- Righe 73-74: catena `int_zd -> int_z`. Il primo integratore trasforma `zddot`
  in `zdot`, il secondo `zdot` in `z`. Condizioni iniziali entrambe a zero,
  coerenti con `zeros(4,1)` usato dalla replica `ode45`
  (`run_full_ascent_simulink.m` riga 49).
- Righe 76-82: identico per il canale di beccheggio, con
  `thetaddot = c5*zdot + c6*theta + c7*delta - c6*alpha_w`.
- Riga 86: `W('c6/1','P6/1'); W('c6/1','P8/1');` -- **una sola lookup** `c6`
  alimenta due prodotti. E' il motivo per cui `A6` compare due volte
  nell'equazione di `thetaddot` (moltiplica `theta` e `-alpha_w`): il termine
  aerodinamico dipende dall'angolo d'attacco totale `theta - alpha_w`, e la
  scrittura a due termini con lo stesso coefficiente lo rende esplicito.
- Righe 95-100: retroazione degli stati sulla **porta 2** dei prodotti. I segnali
  arrivano dalle uscite degli integratori, quindi **non c'e' loop algebrico**:
  ogni anello di feedback attraversa almeno un integratore.

> **Possibile domanda d'esame** -- da dove viene `c5 = A6/V` e non `A6`?
> *Risposta:* Il momento aerodinamico e' proporzionale all'angolo d'attacco
> totale `alpha = theta + zdot/V - alpha_w` (deriva laterale + assetto - vento).
> Sostituendo in `thetaddot = A6*alpha + K1*delta` si ottiene
> `thetaddot = (A6/V)*zdot + A6*theta + K1*delta - A6*alpha_w`: il `1/V`
> nasce dalla conversione della velocita' laterale `zdot` in angolo. Lo stesso
> `V` compare in `c2 = a1*V + a4` per la forza normale.

---

## Il controllore: frozen contro scheduled (righe 102-133)

```matlab
A('sched','built-in/Constant',[...],'Value','sched');
A('Kp_f','built-in/Constant',[...],'Value','Kp_th0');
A('Kp_s',LK,[...],'BreakpointsForDimension1','tsched', ...
   'Table','Kp_sched','ExtrapMethod','Clip');
A('Kp_sw','built-in/Switch',[...], ...
   'Criteria','u2 >= Threshold','Threshold','0.5');
```

- Riga 103: `sched` e' una `Constant` che legge la variabile base omonima.
  `run_full_ascent_simulink.m` (riga 41) fa `assignin('base','sched',sc)` prima
  di ogni `sim`, quindi lo stesso `.slx` serve entrambi gli esperimenti.
- Righe 104-105: `Kp_f`, `Kd_f` sono i guadagni **congelati**, cioe' il design
  HM3 Task-1 a max-q (`Kp_th0`, `Kd_th0` prodotti da `design_controller` su
  `build_plant_rigid(load_hw3_params())`, `t_ref = 72 s`).
- Righe 108-109: `Kp_s`, `Kd_s` sono lookup **1-D su `tsched`**, con tabelle
  `Kp_sched` / `Kd_sched`. Interpolazione lineare (default). I punti di design
  sono `tsched = 5:5:140` -> **28 punti**, uno ogni 5 s: `init_simulink_lpv`
  (righe 75-85) risolve 28 volte `design_controller` su altrettanti plant
  congelati, in **continuazione** (warm start `Kprev` dal punto precedente),
  esattamente la ricetta di continuazione usata altrove nella repo.
- Righe 106-107 e `ExtrapMethod='Clip'`: qui la clip **serve davvero**. La
  griglia dei guadagni parte da `t = 5 s` mentre la simulazione parte da
  `t = 0`, quindi fra 0 e 5 s si estrapola per forza. `Clip` tiene il valore
  d'estremo, che e' esattamente cio' che fa `griddedInterpolant(...,'nearest')`
  nella baseline. Il commento del codice lo dice ed e' corretto.
- Righe 111-112: `Switch` con `Criteria = 'u2 >= Threshold'`, `Threshold = 0.5`.
  Semantica: se `sched >= 0.5` passa la **porta 1** (schedulato), altrimenti la
  **porta 3** (congelato). Il ramo non selezionato **non** viene eseguito: il
  builder non tocca il parametro di configurazione `ConditionallyExecuteInputs`,
  che di default e' `'on'`, e la *conditional input branch execution* di Simulink
  calcola solo il ramo effettivamente selezionato quando quel ramo alimenta in
  esclusiva lo Switch -- ed e' il caso di `Kp_s`/`Kp_f` verso `Kp_sw` e di
  `Kd_s`/`Kd_f` verso `Kd_sw`. Quindi non si paga nemmeno il costo di calcolo del
  ramo inattivo.
- Righe 117-120: qui c'e' una asimmetria da dichiarare. `Kp` e `Kd` di beccheggio
  passano da `Product` (guadagno = segnale variabile). I guadagni di **deriva**
  `Kp_z0` e `Kd_z0` sono invece blocchi `Gain` con valore costante:
  **non sono mai schedulati**, nemmeno con `sched = 1`. Il "gain scheduling" di
  questo showcase riguarda solo la coppia PD di assetto. Coerente con
  `ode_lpv_ascent.m` (riga 22, usa sempre `M.Kp_z`, `M.Kd_z`) e con la tabella
  del README, che infatti dice "pitch gains".
- Riga 121: `delta_sum` con `Inputs = '----'`, cioe'

      delta = -( Kp*theta + Kd*thetadot + Kp_z*z + Kd_z*zdot )

  Legge di controllo PD puro con riferimento `theta_ref = 0`, identica alla riga
  22 di `ode_lpv_ascent.m`.
- Righe 132-133: `delta` rientra nel plant sulle porte 2 di `P3` (`c3*delta`) e
  `P7` (`c7*delta`).

**Limite strutturale importante**: in questo modello rigido `delta` va dal
sommatore direttamente al plant. Non c'e' **nessun attuatore TVC, nessun ritardo
di trasporto, nessuna saturazione, nessun rate limit, nessun modo flessibile,
nessun notch**. E' l'ipotesi di "attuatore ideale" dichiarata nell'header di
`ode_lpv_ascent.m`. Qualunque affermazione sui margini letta su questo anello va
qualificata: e' l'anello rigido idealizzato, non quello di HM3 Task 2/3.

---

## Logging (righe 135-146)

```matlab
tw = {'-1','MaxDataPoints','inf'};
A('log_theta','built-in/ToWorkspace',[...], ...
  'VariableName','theta_sl','SaveFormat','Timeseries', ...
  'SampleTime',tw{:});
```

- Riga 136: `SampleTime = '-1'` (ereditato: su un segnale continuo diventa "ogni
  passo maggiore del solver") e `MaxDataPoints = 'inf'` per togliere il tetto di
  default a 1000 punti. Senza quest'ultimo, il confronto con `ode45` sarebbe
  fatto su una manciata di campioni.
- Righe 137-146: cinque `To Workspace` in formato `Timeseries`: `theta_sl`,
  `z_sl`, `zdot_sl`, `delta_sl`, `alpha_w_sl`. E' il **contratto** con
  `run_full_ascent_simulink.m`, che li legge come `so.theta_sl` ecc.

Osservazioni oneste sul logging:

- `zdot_sl` viene loggato ma **`run_full_ascent_simulink.m` non lo legge mai**
  (riga 44 preleva solo `theta_sl`, `z_sl`, `delta_sl`, `alpha_w_sl`). Log
  inutilizzato.
- Riga 37: il modello abilita `SignalLogging` con `SignalLoggingName = 'logsout'`,
  ma **nessun segnale viene marcato per il logging** in tutta la build. `logsout`
  resta vuoto e non viene mai letto. E' configurazione morta, probabilmente
  trascinata da `run_wind_generator` (dove `logsout` serve davvero).

### Perche' `alpha_w_sl` e' loggato: il punto piu' delicato dell'intera validazione

`run_full_ascent_simulink.m` (riga 92) costruisce la `windfun` della replica
`ode45` come `griddedInterpolant` sull'`alpha_w` **che il modello Simulink ha
appena prodotto**. Il vento **non viene rigenerato** in MATLAB.

Va detto chiaramente: questa non e' una validazione indipendente, e' un
**controllo di autoconsistenza** con il ramo del vento cortocircuitato per
costruzione. Cio' che viene davvero confrontato e': tabelle dei coefficienti,
cablaggio del plant, legge di controllo, integrazione. Il generatore di vento e'
comune ai due lati e quindi non e' sotto test. L'header di
`run_full_ascent_simulink.m` (righe 3-5) lo ammette: "the residual is just wind
interpolation between solver steps".

---

## Salvataggio (righe 148-153)

- Righe 149-150: `save_system` scrive `hm3_full_ascent.slx` **accanto allo
  script**. Il `.slx` e' quindi un artefatto derivato che vive nella stessa
  cartella del suo generatore.
- Riga 152: apre o chiude senza salvare a seconda di `o.open`.

### Perche' costruire il modello da script

1. **Riproducibilita'.** Il `.slx` e' binario: non si diffa, non si revisiona,
   non si merge. Il `.m` si'. Chi legge la repo vede *tutto* il modello --
   blocchi, parametri, fili -- in 153 righe di testo.
2. **La sorgente resta unica.** Se cambia il design (un guadagno, una tabella,
   il passo della schedula), si modifica `init_simulink_lpv.m` o questo file e si
   rilancia: il modello si rigenera. Non esiste il rischio classico del
   "modello Simulink che ha divergito dallo script".
3. **Contrasto voluto con HM3.** Il README della cartella lo dice esplicitamente
   (righe 54-57): `HM3/models/hm3_closed_loop.slx` e' costruito a mano seguendo
   una guida; qui il modello **e'** il codice.

Il prezzo: il modello **non e' modificabile a mano** (ogni edit si perde alla
prossima build) e **non e' autosufficiente** (tutti i parametri di dati sono nomi
di variabili base; senza `init_simulink_lpv` il `.slx` non compila nemmeno,
perche' persino `StopTime` e' la stringa `'Tstop'`). Una scelta piu' robusta sarebbe
stata mettere i dati nel Model Workspace o in un data dictionary.

---

## Limiti, hack e codice stale (riepilogo onesto)

- **Nessuna tolleranza impostata sul modello**: Simulink gira a `RelTol` di
  default (1e-3) contro `1e-9` della replica. L'accordo dipende da `MaxStep`.
- **`SignalLogging`/`logsout` attivi ma vuoti** (riga 37): configurazione morta.
- **`zdot_sl` loggato e mai letto.**
- **`ExtrapMethod` non impostato sulle lookup `c1..c7`** (default `Linear`),
  mentre la baseline usa estrapolazione `nearest`. Innocuo solo perche' la
  griglia copre l'orizzonte; incoerente con il modello flessibile, che mette
  `Clip` ovunque.
- **`Vsafe = max(V,1)` in `init_simulink_lpv` e' inerte** e il commento
  "V(0)=0 ... at lift-off" e' smentito dai dati (`min(V) = 410.4 m/s`).
- **Attuatore ideale**: niente TVC, ritardo, saturazione, rate limit, bending.
- **Guadagni di deriva mai schedulati**, nemmeno con `sched = 1`.
- **`addpath` non ripristinato.**
- **La validazione e' autoconsistente**, non indipendente: il vento e' preso dal
  modello stesso.
- Le variabili base `lpv_Q`, `lpv_V`, `lpv_h`, `Tstart` sono pushate da
  `init_simulink_lpv` ma **nessun blocco di questo modello le usa**.
- Il modello copia il sottosistema del vento: se il file del professore cambiasse,
  serve un rebuild manuale, senza alcun avviso automatico.

---

## Possibili domande d'esame

**D: Come si realizza un plant LPV in Simulink, visto che il blocco State-Space
prende matrici costanti?**
R: Si scende al livello dei blocchi elementari. Il Clock genera il parametro di
scheduling (qui il tempo di volo), una lookup table 1-D per coefficiente
restituisce `c_i(t)`, un `Product` moltiplica il coefficiente per lo stato o
l'ingresso corrispondente, un `Sum` compone la derivata e una catena di
`Integrator` chiude lo stato. Il blocco `State-Space` e' escluso perche' valuta
`A,B,C,D` una sola volta a compile time. Le alternative sono un `MATLAB Function`
block (che pero' rende il modello una scatola nera) o l'`LPV System` block
(matematicamente equivalente a parita' di griglia, ma dipendente da un toolbox e
molto piu' pesante: servirebbe un array di 141 modelli `ss` per rappresentare
7 scalari).

**D: Cos'e' il "trucco dei coefficienti effettivi" e perche' e' indispensabile?**
R: Invece di tabulare i fattori grezzi `a1, V, a4, ...` e ricombinarli con
blocchi, si tabulano i coefficienti **gia' combinati** (`c2 = a1*V + a4`,
`c5 = A6/V`, ...) e sia il `.slx` sia la RHS `ode45` interpolano le stesse
tabelle. Serve perche' l'interpolazione lineare non commuta con il prodotto: su
una cella, `L{f*g} - L{f}*L{g} = Df*Dg*u*(1-u)`, che si annulla solo sui
breakpoint e vale fino a `Df*Dg/4` a meta' cella. Se i due lati combinassero i
fattori in momenti diversi, integrerebbero **funzioni diverse** e il residuo
avrebbe un pavimento di modello indipendente dalle tolleranze. E' proprio la
condivisione delle tabelle combinate che permette l'accordo a ~1e-7 rad.

**D: Perche' `MaxStep = 0.02` e non un valore piu' comodo?**
R: Il generatore di turbolenza del professore aggiorna il rumore ogni 0.1 s.
Con `MaxStep = 0.02` si hanno 5 passi per innovazione: il solver non scavalca le
discontinuita' e i campioni di `alpha_w` loggati bastano a ricostruire lo stesso
ingresso nella replica `ode45`. E' anche cio' che, di fatto, garantisce
l'accuratezza del modello, dato che le tolleranze di Simulink restano ai valori
di default.

**D: Il confronto Simulink vs `ode45` a 1e-7 rad dimostra che il modello e'
corretto?**
R: Dimostra che le due implementazioni della **stessa** matematica coincidono. Non
e' una validazione indipendente: la replica `ode45` viene alimentata con
l'`alpha_w` che il modello Simulink ha appena generato e loggato, quindi il ramo
del vento e' comune ai due lati ed esce dal confronto. Cio' che il test copre e'
il cablaggio del plant, le tabelle dei coefficienti e la legge di controllo; cio'
che non copre e' il generatore di vento e la correttezza fisica del modello LPV.

**D: In cosa consiste esattamente il gain scheduling di questo modello?**
R: `init_simulink_lpv` risolve `design_controller` su 28 plant congelati
(`tsched = 5:5:140`, uno ogni 5 s), in continuazione con warm start dal punto
precedente, e ne salva `Kp_sched`, `Kd_sched`. Nel modello, due lookup 1-D su
`tsched` con estrapolazione `Clip` (necessaria perche' la simulazione parte a
t=0 mentre la griglia parte a t=5 s) restituiscono `Kp(t)`, `Kd(t)`, e due
`Switch` pilotati dalla costante `sched` scelgono fra questi e la coppia
congelata a max-q. **Solo** i guadagni di beccheggio sono schedulati: quelli di
deriva restano bloccati al design di max-q in entrambe le configurazioni.

**D: Perche' costruire il modello da script invece che disegnarlo?**
R: Perche' cosi' il modello e' testo: diffabile, revisionabile, riproducibile.
La build e' distruttiva e idempotente (`new_system` da zero), quindi il `.m` resta
l'unica fonte di verita' e il `.slx` un artefatto derivato che si rigenera quando
il design cambia; sparisce il rischio che il modello Simulink diverga
silenziosamente dallo script. Il costo e' che il modello non si puo' piu'
ritoccare a mano (ogni edit si perde al rebuild) e che non e' autosufficiente:
ogni parametro di dati dei blocchi e' un nome di variabile del base workspace,
quindi senza `init_simulink_lpv` il `.slx` non compila.
