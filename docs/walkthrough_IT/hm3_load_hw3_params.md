# HM3/load_hw3_params.m

## Ruolo del file nel progetto

E' il **punto unico di verita' numerica** di tutto HM3. Restituisce una struct `p`
con i coefficienti del modello pitch-plane del lanciatore fittizio *Greensite*
congelati all'istante di massima pressione dinamica (`t = 72 s`, il punto
strutturalmente critico dell'ascesa). Ogni altro file di HM3 parte da qui:
`build_plant_rigid(p)`, `build_plant_full(p)`, `build_tvc(p)`,
`build_notch_filter(p)`, `load_wind_profile(p)`, i tre `main_task*.m`, il Monte
Carlo e l'estensione LPV in `LTV_FULL_ASCENT/`.

Il file fa tre cose in cascata:

1. scrive i **letterali della Table 1** della traccia (geometria, masse, forze,
   coefficienti aero/controllo, modo di bending, sensore INS, attuatore TVC);
2. se il dataset di riferimento `General/hw3-v3/GreensiteLPV_DATA.mat` e'
   presente, **sovrascrive i coefficienti tempo-varianti** interpolandoli a
   `t_ref` -- la Table 1 resta come fallback e come cross-check;
3. applica gli **scaling di incertezza** `mu_alpha_scale` / `mu_c_scale` usati
   nelle quattro corner di Task 3 (+/-30 % su `mu_alpha` e `mu_c`).

Il parametro `t_ref` (default 72) e' quello che trasforma questa funzione da
"caricatore di una tabella" a **griglia LPV**: `LTV_FULL_ASCENT/main_full_ascent.m`
la chiama in ciclo su `t_ref = 5...140 s` (griglia `tm = t0:2.5:Tstop` con
`t0 = 5 s`, `Tstop = 140 s`, default di `init_simulink_lpv.m`) per costruire il
modello tempo-variante di tutta l'ascesa. La stessa funzione serve quindi sia il
design frozen-time (HM3 propriamente detto) sia l'estensione fuori-traccia.

I valori interpolati dal dataset LPV a `t = 72 s` (verificati eseguendo il
codice) sono:

    A6 = 3.381828   K1 = 4.564724   a1 = -0.015423
    a3 = 20.609026  a4 = -27.270977 V  = 937.7087 m/s
    Tc = 1.5213e6 N wBM = 18.9 rad/s sigma_ins = 0.178  phi_ins = 0.8
    phi_tvc = 4.309544e-05 1/kg      rho = 0.18453      qbar = 81129 Pa

cioe' coincidono con la Table 1 fino alla 4a-5a cifra: e' esattamente il
cross-check che il test `tests/hm3PlantTest.m` (righe 23-31) pinna, con
tolleranze diverse per grandezza -- `AbsTol = 5e-3` su `A6` e `K1` (righe 26-27),
`0.5` su `V`, `0.05` su `wBM`, `5e-2` su `a4` (righe 28-30).

---

## `load_hw3_params` -- firma e blocco `arguments` (righe 1-23)

```matlab
function p = load_hw3_params(opt)
arguments
    opt.mu_alpha_scale (1,1) {mustBeNumeric, mustBeReal} = 1.0
    opt.mu_c_scale     (1,1) {mustBeNumeric, mustBeReal} = 1.0
    opt.t_ref          (1,1) {mustBeNumeric, mustBeReal} = 72
end
```

- Riga 1: firma con **soli argomenti name-value** (`opt`). Chiamata nominale:
  `p = load_hw3_params()`. Chiamata di corner: `load_hw3_params('mu_alpha_scale',1.3,'mu_c_scale',0.7)`.
- Righe 19-23: la validazione qui c'e' (e' un helper di boundary, chiamato una
  volta per run e non dentro un loop di ottimizzazione), coerente con la
  convenzione del repo.
- Righe 3-16: l'header documenta i campi. Attenzione: la riga 9 dichiara
  `A6,K1,a1,a3,a4` tutti in `[1/s^2]`. **Non e' vero dimensionalmente**: vedi il
  blocco successivo. E' un'imprecisione del commento, non del codice.

---

## Letterali Table 1: geometria, masse, forze (righe 25-35)

```matlab
p.m          = 7.38e4;   % kg
p.l_alpha    = 10.39;    % m
p.l_c        = 9.84;     % m
p.Iyy        = 3.28e6;   % kg m^2
p.Alt        = 15143;    % m
p.Tt_minus_D = 1.71e6;   % N
p.N_alpha    = 1.07e6;   % N/rad
```

- Righe 30-31: `l_alpha` e' la **distanza fra centro di pressione e baricentro**
  (braccio della forza normale aerodinamica), `l_c` la distanza fra il **punto di
  gimbal** del motore e il baricentro (braccio del controllo). Sono i due bracci
  che generano i due momenti in competizione.
- Riga 35: `N_alpha` e' la derivata della **forza normale** rispetto all'angolo
  d'attacco, `dN/dalpha` [N/rad]. E' il tipico `qbar*S*CN_alpha`.

Da questi numeri si ricostruiscono i coefficienti della Table 1 -- ed e' la
derivazione da saper rifare a mente:

    A6 = mu_alpha = N_alpha * l_alpha / Iyy
                  = 1.07e6 * 10.39 / 3.28e6 = 3.39 1/s^2      (tab: 3.3818)

    K1 = mu_c     = Tc * l_c / Iyy
                  = 1.52e6 * 9.84 / 3.28e6 = 4.56 1/s^2       (tab: 4.5647)

    a1 = -N_alpha / (m*V)
       = -1.07e6 / (7.38e4 * 937.7) = -0.01546 1/s            (tab: -0.0154)

    a3 = Tc / m = 1.52e6 / 7.38e4 = 20.60 m/s^2 per rad       (tab: 20.6090)

Il significato: `A6` (= `mu_alpha`) e' l'**accelerazione angolare per unita' di
angolo d'attacco** prodotta dalla forza normale aerodinamica applicata **davanti**
al baricentro; `K1` (= `mu_c`) e' l'**efficacia di controllo** del TVC, cioe'
l'accelerazione angolare per unita' di deflessione dell'ugello. Il rapporto
`A6/K1 = 0.74` e' la misura di quanto controllo serve per contrastare
l'instabilita' aerodinamica: sono le due sole quantita' che `design_controller.m`
legge dal plant (righe 44-45 di quel file) per calcolare i guadagni di partenza
`Kp = 2*A6/K1`, `Kd = sqrt(A6)/K1`.

- **Nota di unita' (il commento della riga 9 sbaglia).** Solo `A6` e `K1` sono in
  `1/s^2` (per radiante). `a1 = -N_alpha/(m*V)` ha unita' `1/s` (moltiplica
  `zdot` [m/s] e deve dare `m/s^2`); `a3` e `a4` sono in `m/s^2` per radiante
  (moltiplicano `delta` e `theta`). Il fatto che il codice funzioni comunque
  dipende dal fatto che i numeri entrano nelle matrici A/B alle righe giuste, non
  dal commento.

> **Possibile domanda d'esame** -- Perche' `A6 > 0` implica instabilita' statica?
> *Risposta:* `A6 = N_alpha*l_alpha/Iyy` e' positivo solo se la forza normale
> aerodinamica agisce **davanti** al baricentro (centro di pressione a monte del
> CG, `l_alpha > 0` col segno della traccia). In quella configurazione un angolo
> d'attacco positivo genera un momento che **aumenta** l'angolo d'attacco:
> feedback positivo. La dinamica rotazionale isolata diventa
> `theta_ddot = A6*theta`, cioe' `s^2 - A6 = 0`, con un polo reale a
> `+sqrt(A6) = +1.84 rad/s`. Nessun lanciatore snello e' staticamente stabile in
> volo atmosferico: la stabilita' e' un prodotto del controllo, non della
> configurazione.

---

## L'incongruenza su `a4` (righe 36-38, 44)

```matlab
% a4 inconsistency: -(Tt-D)/m = -23.17 1/s^2 from the values above, but the
% table lists a4 = -27.2710. The LPV set agrees with -27.2710, so that is
% what enters the dynamics.
p.a4 = -27.2710;
```

- Righe 36-38: il commento e' onesto ed e' **corretto**. La formula strutturale
  del coefficiente e' `a4 = -(T - D)/m`; con i valori di tabella
  (`Tt_minus_D = 1.71e6 N`, `m = 7.38e4 kg`) si ottiene `-23.17`, mentre la
  Table 1 dichiara `-27.2710`. Per ottenere `-27.2710` servirebbe
  `T - D = 27.271 * 7.38e4 = 2.01e6 N`, cioe' ~18 % in piu' della spinta netta
  dichiarata.
- Il codice **non risolve** l'incongruenza: usa `-27.2710` perche' e' quello che
  il dataset LPV conferma (`a4 = -27.270977` a `t = 72 s`). E' la scelta giusta
  (il dataset e' la fonte del professore), ma va detto all'orale: la Table 1 non
  e' internamente consistente sulla riga di `Tt - D`, oppure `Tt_minus_D` e'
  tabulato al netto di qualcos'altro (per es. la sola spinta al vuoto contro la
  spinta reale, o una D valutata diversamente).
- La combinazione che entra davvero nella dinamica laterale e'
  `a1*V + a4 = -14.46 - 27.27 = -41.73 m/s^2 per rad`, che con la formula
  strutturale varrebbe `-[(T-D) + N_alpha]/m = -(1.71e6 + 1.07e6)/7.38e4 = -37.7`.
  La differenza (~10 %) sposta di poco i poli ma **non** cambia la fisica.

---

## Bending, sensore INS, attuatore TVC (righe 47-57)

```matlab
p.wBM      = 18.9;     % rad/s
p.zBM      = 0.005;    % -
p.phi_ins  = 0.8;      % -
p.sigma_ins= 0.178;    % rad/m
p.phi_tvc  = 4.31e-5;  % 1/kg
p.wTVC = 70;  p.zTVC = 0.7;  p.tau = 0.020;
```

- Righe 48-49: primo modo di **flessione** del corpo. `zBM = 0.005` e' uno
  smorzamento strutturale bassissimo (0.5 %): il fattore di qualita' e'
  `Q = 1/(2*zBM) = 100`, cioe' **+40 dB** di amplificazione alla risonanza
  rispetto al guadagno statico. E' la ragione per cui il modo di bending non si
  puo' ignorare anche se il suo guadagno DC e' piccolo.
- Righe 51-52: coefficienti della **Eq. (2)** della traccia, l'osservazione INS.
  `sigma_ins` [rad/m] e' la **pendenza** della forma modale alla stazione
  dell'INS (contamina `theta_m`), `phi_ins` [-] e' lo **spostamento** modale alla
  stessa stazione (contamina `z_m`). Sono i due numeri che rendono il modo di
  bending **visibile in retroazione**: senza di loro il loop non lo vedrebbe
  affatto (vedi la pagina di `build_plant_full.m`).
- Riga 54: `phi_tvc` [1/kg] e' il **forcing** del bending da parte del TVC:
  moltiplicato per la spinta di controllo `Tc` da' l'accelerazione modale per
  radiante di deflessione, `phi_tvc*Tc = 65.56 m/s^2 per rad`.
- Righe 55-57: attuatore TVC di **Eq. (3)** -- secondo ordine con `wTVC = 70 rad/s`
  (ben oltre la banda di controllo ~2-3 rad/s e anche **sopra** il triplo di `wBM`:
  `3*18.9 = 56.7 rad/s < 70`, quindi l'attuatore da solo non attenua la risonanza
  strutturale -- il roll-off utile arriva dal notch), `zTVC = 0.7`, piu' un ritardo
  puro `tau = 20 ms` che in `build_tvc.m` viene approssimato con un Pade.

---

## Override dal dataset LPV (righe 59-80)

```matlab
datafile = fullfile(fileparts(mfilename('fullpath')), ...
                    'General', 'hw3-v3', 'GreensiteLPV_DATA.mat');
if isfile(datafile)
    S = load(datafile);  L = S.GreensiteLPV;
    at = @(ts) interp1(ts.Time, squeeze(ts.Data), opt.t_ref);
    p.A6 = at(L.A6);  ...  p.wBM = at(L.omega);
    p.src = 'GreensiteLPV_DATA.mat @ t=72 s';
else
    p.src = 'Table 1 literals (data file not found)';
end
```

- Righe 60-61: il path e' costruito con `fileparts(mfilename('fullpath'))`, quindi
  **relativo al file, non alla cwd**: la funzione si puo' chiamare da qualunque
  directory (i main di `LTV_FULL_ASCENT/` lo fanno).
- Riga 65: `at` e' la closure di interpolazione. `squeeze(ts.Data)` e' **difensivo**:
  nel dataset i campi `timeseries` di `GreensiteLPV` hanno `Data` gia' come vettore
  colonna `Nt x 1`, quindi qui `squeeze` e' un no-op -- servirebbe solo se `Data`
  arrivasse nel layout `1 x 1 x Nt` (quello tipico dei log Simulink, dove infatti
  la stessa chiamata ricorre in `run_full_ascent_simulink.m`).
  **`interp1` di default e' lineare e senza extrapolation**: se
  `t_ref` cade fuori dalla griglia temporale del dataset, i coefficienti tornano
  `NaN` silenziosamente (nessun errore). E' un limite reale della funzione quando
  la si usa come griglia LPV.
- Righe 66-76: cosa viene **sovrascritto** -- `A6, K1, a1, a3, a4, V, Tc,
  sigma_ins, phi_ins, phi_tvc, wBM` (da `L.omega`). Cosa **resta letterale** --
  `zBM`, `wTVC`, `zTVC`, `tau`, `m`, `Iyy`, `Alt`, `l_alpha`, `l_c`, `N_alpha`,
  `Tt_minus_D`. Quindi in un sweep su `t_ref` la massa, l'inerzia e soprattutto
  l'**altitudine restano congelate a quelle di t = 72 s**, mentre `V` varia: la
  `qbar` calcolata sotto e' corretta solo al punto di progetto.
- Righe 77-79: `p.src` e' un flag di provenienza stampato dai main -- utile
  all'orale per dire "questi numeri vengono dal dataset del professore, non li ho
  ricopiati a mano".

> **Possibile domanda d'esame** -- Se hai gia' la Table 1, perche' leggere il
> `.mat`?
> *Risposta:* per tre motivi. (1) Cross-check: se i due insiemi coincidono a 4-5
> cifre (e coincidono, il test lo pinna), la trascrizione della tabella e' giusta.
> (2) Precisione: la tabella e' arrotondata a 4 cifre, il dataset no. (3)
> Riusabilita': la stessa funzione con `t_ref` variabile diventa la griglia LPV di
> tutta l'ascesa, che e' quello che serve all'estensione `LTV_FULL_ASCENT/`. La
> Table 1 rimane come fallback se il file non c'e'.

---

## Pressione dinamica derivata (righe 82-84)

```matlab
p.rho  = 1.225 * exp(-p.Alt/8000);  % kg/m^3
p.qbar = 0.5 * p.rho * p.V^2;       % Pa
```

- Riga 83: **atmosfera esponenziale**, la stessa convenzione usata in HM0
  (`rho0 = 1.225 kg/m^3`, `Hscale = 8000 m`). A `Alt = 15143 m` da'
  `rho = 0.1845 kg/m^3`.
- Riga 84: `qbar = 0.5*0.1845*937.71^2 = 81.1 kPa`. E' il numero con cui i main
  convertono l'angolo d'attacco nell'indicatore di carico `qbar*alpha`
  [kPa*deg], la quantita' dimensionante al max-q.
- **Caveat**: `qbar` e' calcolato con `p.Alt` letterale (mai aggiornato dal
  dataset LPV) e `p.V` interpolato. Per `t_ref = 72` e' esatto; per qualunque
  altro `t_ref` la densita' e' quella sbagliata. Lo stesso vale per la dispersione
  del vento in `load_wind_profile.m` (riga 43), che cerca `sigma` nel dataset
  `drywind` all'altitudine `p.Alt/1000` = 15.1 km sempre.

---

## Scaling di incertezza Task 3 (righe 86-92)

```matlab
p.A6 = p.A6 * opt.mu_alpha_scale;   % mu_alpha
p.K1 = p.K1 * opt.mu_c_scale;       % mu_c
```

- Righe 89-90: le uniche due quantita' scalate. **`a1`, `a3`, `a4`, `V`, `Tc` non
  vengono toccati** -- il test `testUncertaintyScalingAppliesToA6K1` (righe 40-47)
  verifica esplicitamente che `a3` resti invariato.
- Le corner di Task 3 (`+/-30 %` su ciascuno, 4 vertici: V1-V4 in `main_task3.m`)
  agiscono quindi solo
  sulla **equazione dei momenti**: si perturba il momento aerodinamico
  destabilizzante e l'efficacia del controllo, non le forze. Fisicamente equivale
  a un'incertezza sul **braccio** (spostamento del centro di pressione, errore
  sulla posizione del gimbal / sull'allineamento dell'ugello), non sul modulo di
  `N_alpha` o della spinta. Una perturbazione "vera" del `+30 %` su `N_alpha`
  cambierebbe anche `a1` (che ne dipende linearmente) e quindi il termine di
  drift.
- E' comunque **quello che chiede la traccia**, che definisce l'incertezza su
  `mu_alpha` e `mu_c`; ma conviene saperlo dire, perche' e' esattamente la domanda
  "e se ti muovo anche `a1`?".
- **Ordine delle operazioni**: lo scaling avviene **dopo** il calcolo di `qbar`
  (righe 82-84), quindi `p.qbar` resta il valore nominale anche nelle corner --
  il che e' corretto (la pressione dinamica non dipende da `mu_alpha`) e rende
  confrontabili i carichi `qbar*alpha` fra i vertici.
- Righe 87-88: i due fattori vengono anche **memorizzati** nella struct, cosi' i
  main possono etichettare le figure con la corner corrispondente.

---

## Possibili domande d'esame

**D: Cosa sono fisicamente `mu_alpha` e `mu_c`, e perche' HM3 li tratta come le due
uniche grandezze incerte?**
R: `mu_alpha = A6 = N_alpha*l_alpha/Iyy` e' il momento aerodinamico destabilizzante
per unita' di angolo d'attacco; `mu_c = K1 = Tc*l_c/Iyy` e' il momento di controllo
per unita' di deflessione dell'ugello. Sono i due coefficienti dell'equazione di
beccheggio, cioe' i due che decidono da soli la stabilita' del loop di assetto: il
polo instabile e' `+sqrt(mu_alpha)` e l'autorita' per contrastarlo e' `mu_c`. Sono
anche i due meno noti a priori in un progetto reale: `mu_alpha` dipende dalla
posizione del centro di pressione (fortemente dipendente dal Mach e mal predetta
dalla CFD in transonico), `mu_c` dalla spinta effettiva e dagli allineamenti.
Perturbarli del +/-30 % e' la corner analysis della traccia.

**D: Il file dichiara `Tt_minus_D = 1.71e6 N` ma usa `a4 = -27.271`, che
corrisponde a `2.01e6 N`. Come lo giustifichi?**
R: E' un'incongruenza della Table 1, e il codice la documenta nel commento alle
righe 36-38 invece di nasconderla. Ho scelto di usare il valore della tabella /
del dataset LPV (che concordano su `-27.2710`) perche' il dataset e' la fonte
primaria del professore e perche' i coefficienti `a*` sono quelli effettivamente
tabulati; `Tt_minus_D` e' un dato accessorio che nel codice non entra in nessuna
matrice. L'effetto sulla dinamica e' comunque piccolo: `a4` compare solo nella
riga di `z_ddot` e sposta il polo lento di drift da `(T-D)/(m*V) = 0.0247` a
`0.0291 rad/s`.

**D: Perche' `zBM = 0.005` e' cosi' critico?**
R: Perche' fissa il picco di risonanza. Un oscillatore del secondo ordine ha
amplificazione `1/(2*zeta)` alla risonanza: con `zeta = 0.005` sono 100, cioe'
+40 dB sopra il guadagno statico. E' quello che porta il guadagno di anello a
`+29 dB` a 18.9 rad/s nel modello completo (vedi `build_plant_full.m`), cioe'
ben oltre il punto critico: senza filtro il loop e' instabile. Uno smorzamento
strutturale del 5 % invece che dello 0.5 % renderebbe il problema quasi
irrilevante.

**D: Come diventa una griglia LPV questa funzione?**
R: Tramite `t_ref`. `LTV_FULL_ASCENT/main_full_ascent.m` chiama
`build_plant_rigid(load_hw3_params('t_ref', tm(i)))` su una griglia
`tm = 5:2.5:140 s`, ottenendo un plant congelato per ciascun istante e quindi un
modello LTV per tutta l'ascesa. Il limite noto e' che `Alt`, `m` e `Iyy` non sono
interpolati dal dataset (restano ai valori di 72 s), quindi la `qbar` derivata e'
valida solo al punto di progetto; l'estensione LPV infatti si ricalcola la propria
`qbar` per lo scheduling.

**D: `interp1` senza extrapolation: e' un bug?**
R: Non nel caso d'uso nominale (`t_ref = 72` e' ben dentro la griglia), ma e' una
fragilita': se qualcuno chiedesse `t_ref = 200 s` la funzione restituirebbe una
struct piena di `NaN` senza sollevare alcun errore, e il primo sintomo sarebbe un
plant `ss` con matrici `NaN`. Un `mustBeInRange` sulla griglia del dataset (o un
`'extrap'` esplicito con warning) sarebbe la correzione onesta.
