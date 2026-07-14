# HM3/load_wind_profile.m

## Ruolo del file nel progetto

E' il **generatore del disturbo** di tutto HM3. Produce una storia temporale di
vento e la converte in un **angolo d'attacco di disturbo** `alpha_w(t)`, che e'
l'unico ingresso esogeno del plant (`build_plant_rigid` / `build_plant_full`
hanno ingressi `[delta, alpha_w]`). Restituisce una struct `w` che
`simulate_gust_response` passa direttamente a `lsim`.

Chi lo chiama (verificato con grep): `main_task1.m` (riga 44),
`main_task2.m` (riga 182), `main_task3.m` (riga 43), `main_montecarlo.m`
(riga 45), `init_simulink_hm3.m` (riga 60) e
`LTV_FULL_ASCENT/main_full_ascent.m` (riga 52). Tutti i `main_task*` passano
`Tend = 80` invece del default 12, per lasciar decadere il modo lento di deriva
(`tau ~ 18 s` ad anello chiuso).

Il file fa due cose ben distinte:

1. **profili analitici** (`gust`, `step`, `doublet`) -- deterministici,
   calcolati a mano, ampiezza presa dalla dispersione di `drywind.mat`;
2. **profilo `strongwind`** -- simula il modello Simulink del professore
   (`strong_wind.slx`: vento medio + turbolenza di Dryden schedulata in quota),
   e ne ritaglia una finestra intorno al max-qbar.

Il profilo di default (`gust`, severita' `severe`) e' quello della traccia:
raffica 1-coseno con picco `V_g = 6.38 m/s` -> `alpha_w = 0.39 gradi`.

---

## Firma, docstring e blocco `arguments` (righe 1-28)

```matlab
function w = load_wind_profile(p, o)
% Wind angle-of-attack disturbance for the gust simulation.
%   alpha_w = v_w/V.
```

- Riga 1: `p` e' la struct dei parametri (`load_hw3_params`); di essa usa solo
  tre campi: `p.V` (velocita' relativa), `p.Alt` (quota, per pescare la
  dispersione del vento) e `p.t_ref` (solo nel ramo `strongwind`).
- Righe 19-28: `arguments` con name-value. Default: `profile = 'gust'`,
  `severity = 'severe'`, `Vg = []` (= "calcolamela tu"), `Tg = 3 s`,
  `Tend = 12 s`, `dt = 0.005 s`, `t0 = 1 s`.
- Riga 22: `severity` e' validata con `mustBeMember` -> errore immediato se si
  scrive male. `profile` invece e' solo `mustBeTextScalar`: il controllo vero
  avviene nel `switch` (riga 63-64), che lancia `load_wind_profile:profile`.
  E' testato in `tests/hm3LoopTest.m` (righe 62-64).
- Riga 23: `Vg` e' validato `mustBeScalarOrEmpty` ma **non** `mustBePositive`:
  un `Vg` negativo passerebbe silenziosamente (produce semplicemente una
  raffica di segno opposto).
- `dt = 0.005 s` = 200 Hz. La frequenza di Nyquist e' `pi/dt = 628 rad/s`,
  ampiamente sopra il bending (`wBM = 18.9`) e sopra l'attuatore TVC (70 rad/s):
  il campionamento non introduce aliasing nella `lsim`.

---

## Il ramo `strongwind`: uscita anticipata (righe 30-34)

```matlab
if strcmpi(o.profile, 'strongwind')
    w = local_strong_wind(p, o);
    return;
end
```

Dispatch precoce: il ramo Simulink non condivide nulla con i profili analitici
(nemmeno il calcolo di `Vg`). Vedi sotto, righe 71-112.

---

## Ampiezza della raffica dalla dispersione di `drywind.mat` (righe 36-49)

```matlab
Vg = o.Vg;
if isempty(Vg)
    dwfile = fullfile(fileparts(mfilename('fullpath')), ...
                      'General','hw3-v3','drywind.mat');
    if isfile(dwfile)
        S = load(dwfile); dw = S.drywind;
        alt_km = p.Alt/1000;             % drywind.alt is in km
        sig = dw.sigma.(o.severity);
        Vg = interp1(dw.alt, sig, alt_km, 'linear', 'extrap');
    else
        Vg = 8.0;                        % fallback [m/s]
    end
end
```

- Righe 39-40: il dataset del professore, `HM3/General/hw3-v3/drywind.mat`.
  Ispezionandolo (verificato) contiene:
  - `dw.alt = [1 2 4 6 8 10 12 14 18 20 25 30]` -- quote **in km**;
  - `dw.sigma.light` / `.moderate` / `.severe` -- intensita' di turbolenza
    (deviazione standard, m/s) alle stesse quote;
  - `dw.Lh` -- scale di lunghezza della turbolenza (usate dal modello Dryden,
    non da questa funzione).
- Riga 43: `alt_km = p.Alt/1000` -- la conversione **e' necessaria** perche' la
  griglia e' in km e `p.Alt` in m. Il commento nel codice lo dice
  esplicitamente. E' l'errore classico che qui e' stato evitato.
- Riga 45: interpolazione lineare a `Alt = 15143 m -> 15.143 km`. Il punto cade
  fra 14 e 18 km, quindi l'`'extrap'` non entra mai in gioco al max-qbar (ma
  proteggerebbe una chiamata a quota fuori griglia).

Valori ottenuti eseguendo il codice a `t = 72 s` (`Alt = 15143 m`,
`V = 937.709 m/s`):

| severity | `V_g` [m/s] | `alpha_w` di picco |
|---|---|---|
| `light` | 0.2543 | 0.0155 gradi |
| `moderate` | 2.5686 | 0.1569 gradi |
| **`severe`** (default) | **6.3785** | **0.3897 gradi** |

Il `V_g = 6.4 m/s` e il picco `alpha_w = 0.39 gradi` citati nel README e nel
report vengono **da qui**, non da un numero scritto a mano.

- Righe 46-48: se il `.mat` non c'e', fallback a `Vg = 8.0 m/s`. E' un valore
  arbitrario, piu' severo del dato reale (6.38): non rompe nulla ma **cambia
  silenziosamente il risultato**. Vale la pena saperlo se una simulazione
  restituisce numeri che non tornano con il report.

> **Approssimazione da dichiarare all'orale.** `dw.sigma.severe` e' una
> **deviazione standard di turbolenza** (l'intensita' RMS che il modello di
> Dryden usa per generare un processo stocastico). Qui viene riutilizzata come
> **ampiezza di picco** di una raffica discreta deterministica. Non e' la
> definizione standard di una discrete gust (nelle norme aeronautiche
> l'ampiezza di progetto e' legata alla scala della raffica e a un fattore di
> attenuazione, non a sigma). E' una scorciatoia ingegneristica legittima --
> "prendo l'intensita' severa a questa quota e la uso come picco" -- ma va
> chiamata con il suo nome: `V_g` **non** e' un picco a `1 sigma` di un
> processo, e' un'ampiezza scelta.

---

## I tre profili analitici (righe 51-65)

```matlab
t  = 0:o.dt:o.Tend;
vw = zeros(size(t));
switch lower(o.profile)
    case 'gust'   % 1-cosine discrete gust
        idx = t >= o.t0 & t <= o.t0 + o.Tg;
        vw(idx) = 0.5*Vg*(1 - cos(2*pi*(t(idx)-o.t0)/o.Tg));
```

### `gust` -- la raffica 1-coseno (righe 54-56)

    v_w(t) = 0.5*V_g*[1 - cos(2*pi*(t - t_0)/T_g)]   per t_0 <= t <= t_0 + T_g
    v_w(t) = 0                                        altrove

E' la forma **"1-coseno a periodo pieno"**: l'argomento del coseno percorre
`0 -> 2*pi` sull'intervallo `[t_0, t_0 + T_g]`, quindi la raffica **parte da
zero, sale al picco e torna a zero**. Con i default (`t_0 = 1`, `T_g = 3`):

- `v_w(1) = 0`, `v_w(2.5) = V_g = 6.38 m/s` (picco), `v_w(4) = 0`;
- anche la **derivata** e' nulla agli estremi (`sin(0) = sin(2*pi) = 0`), quindi
  il profilo e' `C^1`: non ci sono gradini ne' spigoli.

La continuita' `C^1` non e' un dettaglio estetico. Lo spettro di un impulso
1-coseno (finestra di Hann) decade asintoticamente come `1/w^3`: alla frequenza
del bending (`wBM = 18.9 rad/s`, cioe' circa 9 volte `2*pi/T_g = 2.09 rad/s`)
l'energia della raffica e' ormai decine di dB sotto il picco. Conseguenza
pratica e **limite onesto**: la risposta al gusto **non e' un test del modo
flessibile**. La raffica eccita la banda rigida (crossover 2.45-3.2 rad/s), che
e' esattamente dove serve, ma la stabilizzazione del bending e' giustificata in
frequenza (Nichols, `|L(wBM)| = -18 dB`), non da questa simulazione temporale.

### `step` (righe 57-58) e `doublet` (righe 59-62)

- `step`: `v_w = V_g` per `t >= t_0`. Discontinuo -> spettro `1/w`, contenuto
  a banda larga. Serve a leggere il comportamento **a regime**.
- `doublet`: `+V_g` sulla prima meta' di `T_g`, `-V_g` sulla seconda. Due
  gradini -> eccita bene, ma e' un'onda quadra: contiene energia a tutte le
  frequenze, bending compreso. E' il profilo giusto se si vuole *provocare* il
  modo flessibile, non quello della traccia.

Nessuno dei due e' usato dai `main_task*`; esistono per i test
(`tests/hm3LoopTest.m`, riga 50) e per esplorazione.

---

## Dal vento all'angolo d'attacco: `alphaw = vw/p.V` (righe 67-69)

```matlab
w = struct('t',t,'vw',vw,'alphaw',vw/p.V, 'V',p.V, ...
           'Vg',Vg,'Tg',o.Tg,'profile',o.profile,'severity',o.severity);
```

E' **una riga sola**, ed e' il cuore fisico del file.

### La derivazione

Il veicolo vola con velocita' relativa `V` lungo il suo asse. Un vento
**laterale** `v_w` non cambia (al primo ordine) il modulo della velocita'
relativa, ma ne **ruota la direzione**: il vettore velocita' dell'aria rispetto
al veicolo acquista una componente trasversale `-v_w`. L'incidenza indotta e'
quindi l'angolo fra l'asse e il nuovo vento relativo:

    alpha_w = arctan(v_w / V) ~ v_w / V        (piccoli angoli)

Numeri: `alpha_w = 6.3785/937.709 = 0.006802 rad = 0.3897 gradi`. L'errore
dell'approssimazione ai piccoli angoli e' `alpha^2/3 ~ 1.5e-5` in relativo:
del tutto trascurabile. Ha senso ricordarlo: **`V` e' 150 volte piu' grande di
`v_w`**, quindi anche una raffica "severa" di 6.4 m/s produce meno di mezzo
grado di incidenza. Non e' la raffica a essere grande, e' il **`qbar` a essere
enorme**.

### Perche' il vento entra come disturbo su `alpha` e non come forza

Perche' nel modello linearizzato **tutta** la forza e il momento aerodinamici
sono proporzionali all'incidenza. La riga di `theta_ddot` del plant e'

    theta_ddot = A_6*theta + (A_6/V)*z_dot + K_1*delta - A_6*alpha_w
               = A_6*(theta + z_dot/V - alpha_w) + K_1*delta
               = A_6*alpha + K_1*delta

e la riga di `z_ddot`, raccogliendo `a_1`:

    z_ddot = a_1*z_dot + (a_1*V + a_4)*theta + a_3*delta - a_1*V*alpha_w
           = a_1*V*(theta + z_dot/V - alpha_w) + a_4*theta + a_3*delta
           = a_1*V*alpha + a_4*theta + a_3*delta

**Entrambe** le righe dipendono dal vento **solo attraverso la combinazione**
`alpha = theta + z_dot/V - alpha_w`. Di qui il fatto -- verificabile leggendo
`build_plant_rigid` -- che la colonna del vento e'

    B_w = [0; -a_1*V; 0; -A_6]

cioe' *la colonna di `theta` privata dei termini non aerodinamici e cambiata di
segno*: `A_6` (momento) e `a_1*V` (forza normale) rispondono all'incidenza,
`a_4` no (e' la componente di spinta/resistenza, che dipende dall'assetto, non
dall'incidenza). Questo e' il motivo strutturale per cui **non serve un ingresso
di forza separato**: modellare il vento come `alpha_w` e' *esatto* dentro
l'ipotesi di linearizzazione, non e' una comodita'.

Ordine di grandezza della forza corrispondente:
`N_alpha * alpha_w = 1.07e6 * 0.006802 = 7.3 kN`, che su una massa di 73.8 t
vale `7278/73800 = 0.099 m/s^2`, cioe' appena **0.01 g** di accelerazione
laterale; e momento `N_alpha*l_alpha*alpha_w = 75.6 kN*m`. Il TVC deve
compensarlo con `delta ~ (A_6/K_1)*alpha_w = 0.741*0.39 = 0.29 gradi` a regime
-- coerente con il picco simulato di 0.53 gradi (che include il transitorio).

> **IL SEGNO DI `alpha_w` -- il punto d'orale piu' redditizio di questa pagina.**
> Codice e report usano oggi, coerentemente, il segno **meno**:
>
>     alpha = theta + z_dot/V - alpha_w
>
> Lo usa il **plant** (colonna `B_w = [0; -a_1*V; 0; -A_6]`, `build_plant_rigid.m`
> riga 17) e lo usa il **post-processing** (`simulate_gust_response.m`, riga 29).
> Il report lo deriva in `Introduction.tex` (righe 61-64) e lo riusa
> nell'indicatore di carico (riga 149, e `Task1.tex` riga 240).
>
> **Perche' il meno e' l'unico segno possibile.** `alpha_w = v_w/V` con `v_w` la
> velocita' *dell'aria*. L'incidenza e' l'angolo fra l'asse del veicolo e il
> vento **relativo**, quindi dipende dalla velocita' laterale del veicolo
> **rispetto all'aria**, `z_dot - v_w`, e non da `z_dot` in assoluto:
>
>     alpha = theta + (z_dot - v_w)/V = theta + z_dot/V - alpha_w
>
> Il segno non e' una convenzione libera: e' **imposto dall'Eq. (1)**. La colonna
> del vento ha `-a_1*V` e `-A_6`, cioe' esattamente i coefficienti con cui
> `theta` entra tramite `alpha`, **cambiati di segno**. Se l'incidenza fosse
> `theta + z_dot/V + alpha_w`, la colonna dovrebbe essere `[0; +a_1*V; 0; +A_6]`.
> Un `+` nel post-processing contraddirebbe il plant che si sta simulando.
>
> **Prova indipendente (la piu' bella da citare all'orale).** Con la retroazione
> di deriva disattivata (`Kp_z = Kd_z = 0`) e una raffica a gradino, la
> simulazione converge a `z_dot = +6.379 m/s = +V*alpha_w` -- esattamente il
> valore che **annulla** `alpha = theta + z_dot/V - alpha_w` con `theta = 0`.
> Cioe': il veicolo si lascia trascinare sottovento finche' la sua velocita'
> laterale eguaglia quella del vento, e l'incidenza relativa sparisce. Se il
> segno fosse `+`, l'incidenza a regime raddoppierebbe invece di annullarsi --
> fisicamente assurdo. (Coerentemente, `Task1.tex` riga 261 riporta
> `z_dot_ss = +V*alpha_w`.)
>
> **Nota storica, ed e' esattamente la domanda che un esaminatore ama.** Fino a
> poco fa `simulate_gust_response.m` calcolava l'indicatore di carico con il
> segno **piu'**, mentre il plant aveva (correttamente) il meno. Il plant non e'
> mai stato sbagliato: il bug era **solo** nel post-processing, cioe' nella
> ricostruzione a valle di `alpha`. Conseguenze numeriche, sul Task 1:
>
> | formula | picco \|alpha\| | `qbar*alpha` |
> |---|---|---|
> | `+alpha_w` (vecchio, sbagliato) | 0.255 gradi | 20.7 kPa*grado |
> | **`-alpha_w` (attuale, corretto)** | **0.577 gradi** | **46.8 kPa*grado** |
>
> Cioe' **piu' del doppio**, e -- soprattutto -- si **ribalta la conclusione
> fisica**. Con il `+`, il picco di incidenza (0.255 gradi) risultava *sotto* il
> contributo del solo vento (0.390 gradi), e sembrava che l'anello facesse
> **load relief**. Con il segno giusto il picco (0.577 gradi) lo **supera**: per
> tenere l'assetto il loop becca il muso **dentro** il vento relativo, e quel
> contributo si **somma** a quello del vento. Una legge di puro attitude-hold e'
> quindi **load-aggravating**, non load-relieving.
>
> **Cosa non e' cambiato:** stabilita', margini, Nichols, notch, Task 3, Monte
> Carlo. `alpha_w` e' un ingresso di **disturbo** e non entra in `L(s)`, quindi
> il segno non tocca la dinamica di anello chiuso. Non cambiano nemmeno i picchi
> di `theta`, `z` e `delta`: cambia **solo** la metrica di carico, cioe'
> l'unica quantita' ricostruita a valle dall'`alpha`.

### Perche' il max-qbar e' la condizione critica

`qbar = 0.5*rho*V^2`. Salendo, `V` cresce e `rho` cala esponenzialmente
(`rho = 1.225*exp(-h/8000)` in `load_hw3_params`, riga 83): il prodotto ha un
**massimo**, che per questo lanciatore cade a `t = 72 s`, quota 15143 m, con
`rho = 0.1845 kg/m^3`, `V = 937.7 m/s` -> `qbar = 81.1 kPa`.

Tre cose diventano critiche **contemporaneamente** in quel punto:

1. **Carico strutturale.** Il momento flettente su un lanciatore snello e'
   proporzionale alla forza normale distribuita, cioe' a `qbar * alpha`. Il
   prodotto `qbar*alpha` (in Pa*grado) e' l'**indicatore di carico** standard, ed
   e' massimo dove `qbar` e' massimo -- tanto piu' che le quote 10-15 km sono la
   fascia del jet stream, cioe' dove i venti *sono* piu' forti. Massimo `qbar` e
   massimo vento nello stesso punto: e' li' che si dimensiona la struttura.
2. **Instabilita' aerodinamica massima.** `A_6 = mu_alpha` e' proporzionale a
   `qbar` (`A_6 = N_alpha*l_alpha/I_yy` con `N_alpha ~ qbar*S*C_N_alpha`): il
   polo instabile `+sqrt(A_6) = +1.84 rad/s` e' **piu' a destra che mai** proprio
   a max-qbar. Il punto piu' difficile da controllare e il punto piu' caricato
   coincidono.
3. **Conflitto di progetto, e con esso un mito da smontare.** Tenere l'assetto
   (`theta -> 0`) contro una raffica significa **mantenere l'incidenza**, cioe'
   **caricare** la struttura. La contromossa istintiva sarebbe "lascia che il
   veicolo si allinei al vento" (*weathervaning*): ma **su questo lanciatore la
   banderuola non esiste**. Con `A_6 > 0` il centro di pressione sta **davanti**
   al baricentro, quindi il momento aerodinamico e' **divergente**, non
   raddrizzante: abbandonato a se' stesso il veicolo non si allinea al vento
   relativo, ci **diverge contro**. Non c'e' nessuna stabilita' a banderuola da
   cui estrarre load relief gratis. E' il motivo per cui il load relief, quando
   serve davvero, va **costruito** (guadagni piccoli e negativi sulla deriva,
   `Kp_z = Kd_z = -1e-3`, o -- in un progetto di volo vero -- una retroazione
   esplicita su un accelerometro o su una stima di `alpha`).

Ecco perche' HM3 progetta a tempo congelato **proprio** a `t = 72 s`: e' il caso
peggiore su tutti e tre i fronti. Con la convenzione corretta il picco vale
`qbar*alpha = 46.8 kPa*grado` (README e report riportano ~47): resta un ordine
di grandezza sotto i limiti strutturali tipici, ma il valore di soglia di
riferimento **non e' nel codice** (non e' calcolato, e' un'affermazione del
testo).

---

## `local_strong_wind` (righe 71-112)

Il ramo "professore": invece di una raffica analitica, simula il modello
`strong_wind.slx` (vento medio + turbolenza di Dryden schedulata in quota) e ne
ritaglia una finestra intorno al max-qbar.

```matlab
load_system(fullfile(gdir,'strong_wind.slx'));
cleanup = onCleanup(@() close_system('strong_wind', 0));   % discard edits
```

- Righe 77-80: carica `drywind.mat` e `GreensiteLPV_DATA.mat` (il modello ne ha
  bisogno come variabili di workspace).
- Righe 82-83: `load_system` (non `open_system`: nessuna finestra) e
  `onCleanup` con `close_system(...,0)` -> **il `.slx` non viene mai
  modificato su disco**, anche se la funzione esce per errore. Il secondo
  argomento `0` significa "chiudi senza salvare". E' una precauzione corretta e
  deliberata (il file e' materiale del corso, non va toccato).
- Righe 85-89: il logging viene attivato **sulle porte di uscita del
  sottosistema**, non sulle linee di segnale -- il commento di riga 87 lo dice.
  E' la ragione per cui si passa da `get_param(...,'PortHandles')` e si fa
  `set_param(ph.Outport(k), 'Name', ..., 'DataLogging','on')`. Questo *e'* un
  edit al modello in memoria: motivo in piu' per l'`onCleanup`.
- Riga 91: `t1 = p.t_ref - o.t0 = 72 - 1 = 71 s`. E' l'**allineamento
  temporale**: la finestra estratta parte 1 s prima del max-qbar, cosi' che
  l'istante locale `t = 1 s` (che nei profili analitici e' l'inizio della
  raffica) corrisponda esattamente a `t = 72 s` di volo. Elegante: le due
  famiglie di profili sono confrontabili nello stesso grafico.
- Righe 92-97: `Simulink.SimulationInput`, con le due variabili iniettate e
  `StopTime = t1 + Tend = 83 s`.
- Righe 99-107: si estraggono i due segnali (`sw_vwp` = vento medio,
  `sw_turb` = turbolenza), si eliminano i tempi duplicati con `unique` (i log a
  passo variabile possono ripetere lo stesso istante -- ed `interp1` con ascisse
  ripetute fallisce), e si **ricampiona a passo fisso** `dt` sulla finestra
  `[t1, t1+Tend]`.
- Riga 106-107: `vw = vento_medio + turbolenza`. Il totale, non solo la
  turbolenza.
- Righe 109-111: stessa struct dei profili analitici, cosi' che
  `simulate_gust_response` non debba sapere da dove viene il vento. Ma:
  `Vg = max(abs(vw))` (non e' piu' un'ampiezza di progetto, e' un **massimo a
  posteriori**) e `Tg = o.Tend` (l'intera finestra: non c'e' una durata di
  raffica).

> **Limiti onesti di questo ramo.**
> - La docstring (riga 76) dice *"Seeds fixed -> reproducible"*, ma **la
>   funzione non fissa nessun seme**: i semi sono dentro `strong_wind.slx`. La
>   riproducibilita' e' quindi una proprieta' del modello, non di questo codice.
> - Le opzioni `Vg`, `Tg` e **`severity` sono ignorate** in questo ramo (la
>   severita' e' scelta dentro il `.slx`), eppure la struct restituita dichiara
>   `'severity','severe'` (riga 111) qualunque cosa l'utente abbia passato. Se
>   si chiama `load_wind_profile(p,'profile','strongwind','severity','light')`
>   non succede nulla e la struct mente.
> - `alphaw = vw/p.V` usa il **`V` congelato a `t = 72 s`** anche se la finestra
>   copre 71-83 s, durante i quali `V` cambia di parecchio. E' coerente con
>   l'ipotesi frozen-time di tutto HM3, ma e' un'approssimazione: nel modello
>   LTV completo (`LTV_FULL_ASCENT/`) il vento va diviso per `V(t)`.
> - `interp1(..., 'linear')` senza `'extrap'`: se il log Simulink non copre
>   tutta la finestra richiesta escono `NaN` silenziosi, che poi propagano
>   dentro `lsim`.

---

## Possibili domande d'esame

**D: Perche' il vento entra nel modello come angolo d'attacco e non come una
forza laterale applicata?**
R: Perche' nel modello linearizzato la forza normale e il momento aerodinamico
sono *entrambi* proporzionali alla sola incidenza `alpha`, e un vento laterale
non fa altro che modificare `alpha` (ruota il vento relativo). Raccogliendo, sia
la riga di `z_ddot` sia quella di `theta_ddot` dipendono dal vento solo tramite
`alpha = theta + z_dot/V - alpha_w`. La colonna del disturbo del plant e'
infatti `B_w = [0; -a_1*V; 0; -A_6]`, cioe' la colonna di `theta` privata del
termine di spinta `a_4` (che non e' aerodinamico) e cambiata di segno. Aggiungere
un ingresso di forza separato sarebbe ridondante: `alpha_w` e' *il* canale del
vento, esattamente, dentro l'ipotesi di linearita'.

**D: Da dove vengono i 6.4 m/s e i 0.39 gradi?**
R: `V_g` non e' un numero scritto a mano: e' `dw.sigma.severe` interpolato
linearmente in quota nel dataset `drywind.mat` del corso, alla quota di max-qbar
(15143 m = 15.143 km, fra i nodi di 14 e 18 km) -> 6.3785 m/s. L'angolo e' la
conversione ai piccoli angoli `alpha_w = V_g/V = 6.3785/937.709 = 0.006802 rad =
0.3897 gradi`. Va detto che `sigma` e' formalmente una **RMS di turbolenza** del
modello di Dryden, riusata qui come ampiezza di picco di una raffica discreta:
e' una scelta, non una definizione da norma.

**D: Perche' una raffica 1-coseno e non un gradino?**
R: Perche' un gradino ha una discontinuita' e quindi contenuto spettrale a tutte
le frequenze: eccita l'attuatore e il modo flessibile per ragioni puramente
matematiche, non fisiche. La 1-coseno parte e finisce con valore **e derivata**
nulli (`C^1`), quindi il suo spettro decade come `1/w^3` e l'energia resta
concentrata attorno a `2*pi/T_g = 2.1 rad/s`, cioe' proprio nella banda del
controllo rigido (crossover 2.4-3.2 rad/s). E' un carico realistico: nell'aria
vera una raffica ha una scala spaziale finita, e attraversarla a `V` produce
esattamente una salita e una discesa. Il rovescio della medaglia (limite da
dichiarare): a `wBM = 18.9 rad/s` la raffica non ha praticamente energia, quindi
**la risposta al gusto non prova nulla sulla stabilita' del bending** -- quella
si dimostra solo in frequenza.

**D: Perche' progettare a max-qbar e non, che so, a decollo o a fine primo
stadio?**
R: Perche' a max-qbar tre criticita' coincidono. (1) `qbar = 81 kPa` e' massimo,
quindi il carico strutturale `qbar*alpha` e' massimo; e la fascia 10-15 km e'
anche quella del jet stream, dove il vento e' piu' forte. (2) `mu_alpha = A_6` e'
proporzionale a `qbar`, quindi l'instabilita' aerodinamica e' massima: il polo a
`+1.84 rad/s` (raddoppio dell'errore in 0.38 s) e' il piu' a destra dell'intero
volo. (3) Il conflitto attitude-hold / load-relief e' massimo. Progettare li'
significa progettare nel caso peggiore; se il progetto tiene a max-qbar, per il
resto dell'ascesa e' conservativo (cosa che l'estensione LPV in
`LTV_FULL_ASCENT/` verifica esplicitamente).

**D: `Tend` di default e' 12 s, ma tutti i `main_task*` passano 80 s. Perche'?**
R: Perche' 12 s bastano a vedere la raffica (1-4 s) e il transitorio d'assetto
(la coppia dominante di anello chiuso ha `tau = 1.05 s`), ma **non** a vedere
esaurirsi il modo lento di deriva laterale, che ad anello chiuso ha
`wn = 0.24 rad/s`, `zeta = 0.23`, cioe' `tau = 17.9 s`. Servono circa `5*tau`
per vedere `z` e `theta` tornare a zero: da qui gli 80 s (il commento di riga 44
di `main_task1.m` lo dice esattamente).

**D: Con che segno entra `alpha_w` nell'incidenza, e come lo giustifichi?**
R: Con il **meno**: `alpha = theta + z_dot/V - alpha_w`. Non e' una convenzione
scelta, e' imposta dall'Eq. (1): la colonna del vento e'
`B_w = [0; -a_1*V; 0; -A_6]`, cioe' i coefficienti con cui `theta` entra tramite
`alpha`, cambiati di segno. Fisicamente, l'incidenza dipende dalla velocita'
laterale **relativa all'aria**, `(z_dot - v_w)/V`. Prova indipendente: togliendo
la retroazione di deriva, sotto raffica a gradino la simulazione converge a
`z_dot = +V*alpha_w`, cioe' esattamente al valore che **annulla** `alpha` --
il veicolo si fa trascinare finche' il vento relativo sparisce. Con il `+`
l'incidenza a regime raddoppierebbe: assurdo.

**D: Quindi il loop d'assetto allevia o aggrava il carico?**
R: Lo **aggrava**. Il picco di incidenza del Task 1 e' `0.577 gradi`, cioe'
**sopra** il contributo del solo vento (`0.390 gradi`): per tenere `theta -> 0`
il controllore becca il muso **dentro** il vento relativo, e il suo contributo
si **somma** a quello del vento invece di elidersi. Una legge di puro
attitude-hold e' quindi **load-aggravating**. Vale la pena aggiungere il perche'
strutturale: con `A_6 > 0` il centro di pressione e' davanti al baricentro,
il momento aerodinamico e' divergente e **non esiste stabilita' a banderuola**
da cui ricavare load relief gratuito. Con `qbar = 81 kPa` il picco
dell'indicatore di carico vale `46.8 kPa*grado`.

*Nota storica (spendibile all'orale).* Fino a poco fa il post-processing
(`simulate_gust_response.m`) ricostruiva `alpha` con il segno **piu'**, mentre il
plant aveva gia' -- correttamente -- il meno. Il difetto dava un picco di
`0.255 gradi` (`20.7 kPa*grado`), **sotto** il vento da solo, e suggeriva una
falsa azione di load relief. Il modello non e' mai stato sbagliato: lo era la
formula a valle. La correzione **non tocca margini, stabilita', Nichols o Task 3**
(`alpha_w` e' un disturbo, non entra in `L(s)`) e nemmeno i picchi di `theta`,
`z`, `delta`: cambia solo la metrica di carico, di un fattore 2.3.
