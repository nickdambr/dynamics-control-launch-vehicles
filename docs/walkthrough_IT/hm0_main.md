# HM0_falcon9_ascent/main.m

## Ruolo del file nel progetto

`main.m` e' l'entry point canonico di HM0 e l'unico file della repo che
implementa per intero la simulazione richiesta dalla traccia
(`LVdynamics_Homework_0.pdf`): propagare in avanti l'ascesa del primo stadio del
Falcon 9 dal Kennedy Space Center fino al MECO, con un modello 3-DoF a punto
materiale, atmosfera esponenziale, Terra sferica rotante e gravita' newtoniana.

E' importante inquadrare bene cosa NON e': **non c'e' ottimizzazione, non c'e'
controllo, non c'e' un problema ai limiti**. Il profilo di assetto e' *assegnato*
(verticale -> pitchover -> gravity turn) e il codice si limita a integrare le
equazioni del moto con `ode45` a partire dalle condizioni iniziali. Tutta la
"intelligenza" del file sta (a) nella derivazione corretta delle equazioni in
coordinate sferiche con velocita' in terna UEN, (b) nella modellazione delle
forze (spinta corretta per la contropressione, drag, gravita'), (c) nel
post-processing che ricostruisce Mach, pressione dinamica, angoli e ground track.

Il file e' **autocontenuto**: costanti, condizioni iniziali, integrazione,
post-processing, otto figure e l'export PNG stanno tutti qui; le equazioni del
moto sono una *local function* `eom` in coda al file (righe 375-466). Non chiama
nessuna funzione della repo. E' chiamato da `tests/falcon9AscentTest.m` tramite
il pattern `run()` + harvest del workspace, e le sue figure sono quelle
referenziate dal README di HM0.

`main2.m` risolve lo stesso identico problema fisico in forma adimensionale con
una parametrizzazione a tre archi: e' una variante preparatoria per gli homework
di ottimizzazione, non un'alternativa a `main.m`.

---

## Intestazione e `clear` (righe 1-6)

```matlab
%% main.m - Falcon 9 First Stage Trajectory Simulation
%  3-DoF point-mass model in spherical coordinates
%  Velocity in Up-East-North (UEN) frame

clear; close all; clc;
```

- Righe 1-4: la testata dichiara le due scelte di modellazione che governano
  tutto il resto: **stato in coordinate sferiche** `(r, theta, phi)` e
  **velocita' proiettata in terna UEN** (Up-East-North). Sono due scelte
  indipendenti: si potrebbe benissimo avere posizione sferica e velocita' in
  componenti cartesiane inerziali. La combinazione scelta e' quella che rende le
  equazioni "leggibili" per un lanciatore (u = rateo di salita, v = velocita'
  verso Est, w = verso Nord).
- Riga 6: `clear` all'inizio. E' un dettaglio innocuo per l'uso interattivo, ma
  ha una conseguenza pesante sulla testabilita': un test che chiami `run('main.m')`
  dal proprio workspace se lo vede cancellare, `testCase` incluso. E' esattamente
  il motivo per cui `falcon9AscentTest.m` deve isolare la `run()` dentro una
  funzione locale usa-e-getta (vedi `hm0_test_falcon9Ascent.md`).

---

## `CONSTANTS AND VEHICLE PARAMETERS` (righe 8-54)

Blocco puramente dichiarativo, ma diversi valori meritano una difesa all'orale.

- Righe 13-14: `mu = 3.986004418e14` m^3/s^2 e `RE = 6378137` m. Sono i valori
  WGS-84. `RE` e' il raggio **equatoriale**: usarlo come raggio di una Terra
  sferica sovrastima leggermente il raggio locale al KSC (lat 28.57 deg): il
  raggio geocentrico dell'ellissoide WGS-84 a quella latitudine vale ~6373 km,
  quindi l'errore e' di ~5 km su 6378, irrilevante rispetto alle altre
  approssimazioni del modello.
- Righe 17-20: modello atmosferico. `rho0 = 1.225` kg/m^3, `p0 = 101325` Pa,
  `Hscale = 8000` m, `Tamb = 288.15` K. **Onesta': queste quattro costanti non
  sono mutuamente consistenti.** Per un'atmosfera isoterma in equilibrio
  idrostatico con gas ideale vale H = R*T/g, che con R = 287.058 e T = 288.15 da'
  H = 8435 m, non 8000 m. `Hscale = 8000` e' un valore ingegneristico
  "arrotondato" molto usato in letteratura. La conseguenza pratica e' che
  l'atmosfera decade un po' piu' in fretta del previsto dal gas ideale isotermo;
  l'effetto e' un max-Q leggermente anticipato e piu' basso. Va detto, non
  nascosto.
- Riga 24: `omegaE = 2*pi/Tsid` con `Tsid = 86136` s (riga 23). **Il giorno
  sidereo standard e' 86164.09 s**, non 86136. La differenza (28 s) porta a
  `omegaE = 7.29449e-5` rad/s invece di 7.29212e-5: un errore relativo di
  3.2e-4. E' del tutto trascurabile sui risultati (la velocita' iniziale Est
  cambia di ~0.1 m/s su 408.6), ma e' un errore reale e conviene saperlo
  ammettere. Il valore e' ripetuto identico in `main2.m` riga 52.
- Riga 29: `a_sound = sqrt(gamma_air * Rgas * Tamb)` = **340.30 m/s**, calcolato
  **una volta sola a T = 288.15 K** e poi usato a tutte le quote. Questa e' una
  conseguenza diretta dell'ipotesi di atmosfera isoterma: se T e' costante, anche
  a lo e'. E' una semplificazione forte: nell'atmosfera standard reale a scende a
  ~295 m/s a 11 km, quindi **questo modello sottostima il numero di Mach in
  quota** (di circa il 13-15% attorno alla tropopausa). L'attraversamento di Mach
  1 riportato (t = 61.8 s, h = 8.04 km) e' quindi un po' *tardivo* rispetto alla
  realta'.
- Righe 32-38: primo stadio. `Tvac1 = 8227 kN`, `cvac1 = 3244` m/s,
  `Aex1 = 11.039` m^2, `tb1 = 162` s. La portata di massa e' **derivata**, non
  assegnata: `Qdot1 = Tvac1/cvac1` = **2536.07 kg/s** (riga 38). Questa e' la
  definizione di velocita' di efflusso *efficace* in vuoto: c_vac = T_vac/mdot,
  cioe' c_vac ingloba gia' il termine di spinta di pressione p_e*A_e. Se si
  usasse la velocita' di efflusso "vera" u_e si otterrebbe una portata diversa.
- **Onesta' sul bilancio di massa**: `Qdot1 * tb1` = 410 843 kg, mentre il
  propellente dichiarato e' `mp1 = 410900` kg (riga 33). I due valori differiscono
  di ~57 kg. La terna (T_vac, c_vac, t_b) non e' perfettamente consistente con
  mp1; il codice non se ne accorge perche' `mp1` entra solo nella massa iniziale
  (riga 69) e mai nella dinamica, che usa esclusivamente `Qdot1`. In pratica il
  veicolo arriva a MECO con 57 kg di propellente non bruciato. Irrilevante, ma se
  all'orale viene chiesto "perche' m(t_b) non e' esattamente m0 - mp1?", questa
  e' la risposta.
- Righe 45-46: `CD = 0.329` **costante** e `Sref = 10.52` m^2. Nessuna dipendenza
  da Mach. Nella realta' il C_D di un lanciatore raddoppia abbondantemente nel
  transonico. Questo modello quindi **sottostima le perdite di drag proprio nella
  regione di max-Q**, che e' dove contano di piu'.
- Righe 53-54: `t_end1 = 5` s, `t_end2 = 15` s. Sono i due switch di fase,
  deterministici e noti a priori: e' il motivo per cui non serve una event
  function nell'integratore (vedi la sezione sull'integrazione).

> **Possibile domanda d'esame** -- Perche' la portata `Qdot` si calcola come
> `Tvac/cvac` e non come `Tvac/(Isp*g0)` o con la velocita' di efflusso vera?
> *Risposta:* Perche' `cvac` e' definita come velocita' di efflusso *efficace* in
> vuoto, cioe' l'impulso specifico in vuoto moltiplicato per g0: per definizione
> T_vac = mdot * c_vac. Le tre scritture sono la stessa cosa. Usare invece la
> velocita' di efflusso *fisica* u_e sarebbe sbagliato, perche' T = mdot*u_e +
> (p_e - p_amb)*A_e: u_e non ingloba il termine di pressione, c_vac si'. Il codice
> e' coerente: usa c_vac per la portata e sottrae separatamente p_amb*A_ex per la
> correzione di contropressione.

---

## `INITIAL CONDITIONS - Launch from Kennedy Space Center` (righe 56-71)

```matlab
r0     = RE;                          % Radius [m]
theta0 = lon0;                        % Right ascension [rad]
phi0   = lat0;                        % Declination [rad]
u0     = 0;                           % Up velocity [m/s]
v0     = omegaE * RE * cos(lat0);     % East velocity [m/s]
w0     = 0;                           % North velocity [m/s]
```

- Righe 60-61: KSC, lat 28.573469 deg N, lon -80.651070 deg.
- Riga 64: `theta0 = lon0`. Questo e' il punto piu' sottile di tutto il file.
  `theta` e' etichettata **"right ascension"**, cioe' un angolo **inerziale**, non
  la longitudine terrestre. Al tempo t = 0 le due coincidono per convenzione
  (si sceglie l'origine dell'ascensione retta coincidente con il meridiano di
  Greenwich a t = 0), ma da li' in poi divergono: la Terra ruota sotto il veicolo.
  La prova che `theta` sia davvero inerziale sta a riga 187, dove la longitudine
  di ground track viene recuperata come `theta - omegaE*t` -- cioe' sottraendo
  esplicitamente la rotazione terrestre.
- Riga 67: `v0 = omegaE * RE * cos(lat0)` = **408.6 m/s**. Il veicolo, fermo sulla
  rampa, ha gia' una velocita' **inerziale** verso Est pari a quella del punto di
  lancio trascinato dalla rotazione terrestre. Questa e' l'unica porta d'ingresso
  della rotazione terrestre nello *stato* del sistema: nelle equazioni non ci sara'
  nessun termine di Coriolis esplicito (vedi la derivazione piu' sotto).
  Conseguenza fisica immediata: nel gravity turn la spinta si allinea alla
  velocita' *relativa*, quindi il veicolo "eredita" e sfrutta questi 408.6 m/s, che
  e' il regalo del lancio verso Est da bassa latitudine.
- Riga 69: `m0 = mdry1 + mp1 + mdry2 + mp2 + mfair + mpay` = **569 100 kg**. Il
  primo stadio porta in volo *tutta* la pila: secondo stadio pieno, fairing e
  payload. Corretto, perche' fino al MECO nulla e' stato sganciato.

Un numero utile da avere in tasca per l'orale, ricavabile dalle righe 60-69 piu'
la `eom`:

    du(0) = v0^2/RE - mu/RE^2 + T_SL/m0
          = 0.026 - 9.798 + 12.491  =  2.72 m/s^2

con `T_SL = Tvac - p0*Aex` = 7.108 MN. Il rapporto spinta/peso al decollo e'
**1.275**: il veicolo stacca, ma senza abbondanza. Il termine centrifugo
(0.026 m/s^2) e' il "regalo" della rotazione terrestre ed e' visibilmente
trascurabile al lift-off.

> **Possibile domanda d'esame** -- Se `theta` e' inerziale, perche' la si
> inizializza con la longitudine geografica del KSC?
> *Risposta:* Perche' e' una scelta di origine dell'asse x inerziale: si fissa il
> frame inerziale in modo che a t = 0 coincida con l'ECEF. Da quel momento
> `theta(t)` continua a essere ascensione retta e la longitudine terrestre si
> recupera come `theta(t) - omegaE*t`. Se si volesse un frame inerziale allineato
> all'equinozio vernale basterebbe aggiungere il GMST iniziale: la dinamica non
> cambierebbe di una virgola, cambierebbe solo l'offset del ground track.

---

## `PARAMETER STRUCTURE` (righe 73-90)

- Righe 77-90: tutti i parametri fisici confluiscono in una `struct par` passata
  alla `eom` per riferimento. Nessun uso di variabili globali, nessun
  `evalin`/`assignin`: la `eom` e' una funzione pura di `(t, y, par)`.
- Nota che la struct **non** contiene `a_sound`, `Tamb`, `mp1`, `tb1`: queste
  servono solo prima o dopo l'integrazione, non dentro la dinamica. E' una
  separazione pulita fra "cosa serve al modello" e "cosa serve al post-processing".
- `par.t1` e `par.t2` (righe 89-90) sono i tempi di switch: la logica di fase e'
  quindi dentro la `eom` e non fuori.

---

## `NUMERICAL INTEGRATION (ode45)` (righe 92-98)

```matlab
opts  = odeset('RelTol', 1e-10, 'AbsTol', 1e-12, 'MaxStep', 1);
tspan = [0, tb1];
[t, Y] = ode45(@(t, y) eom(t, y, par), tspan, y0, opts);
```

- Riga 98: `ode45` = Dormand-Prince esplicito 5(4) a passo adattivo. Scelta
  corretta: il problema e' **non-stiff** (le scale temporali sono tutte
  dell'ordine dei secondi, non ci sono modi veloci) e il RHS e' economico da
  valutare, quindi un metodo esplicito e' l'ottimo.
- Riga 96 (`odeset`): `RelTol = 1e-10`, `AbsTol = 1e-12`: tolleranze molto piu' strette del necessario
  per una semplice propagazione (per un grafico basterebbe 1e-6). Sono la
  convenzione della repo (`CLAUDE.md`) ereditata dal lavoro di shooting/indiretto
  degli homework successivi, dove l'accuratezza della propagazione entra
  direttamente nel residuo del BVP. Qui il costo e' solo tempo di calcolo. Il
  vantaggio collaterale e' che rende sensato il test di bookkeeping del
  propellente a `RelTol 1e-8` (vedi `hm0_test_falcon9Ascent.md`).
- `MaxStep = 1` s: **non** e' una tolleranza, e' un vincolo geometrico sul passo.
  Serve a garantire che la griglia di output sia abbastanza fitta da (a) risolvere
  bene il picco di max-Q e l'attraversamento di Mach 1, che vengono localizzati in
  post-processing con `max`/`find` *sui campioni*, e (b) impedire che un singolo
  passo salti a pie' pari uno degli intervalli di fase -- in particolare la fase 2
  dura solo 10 s.
- Riga 97: `tspan = [0, tb1]`. Il MECO e' **imposto per tempo**, non rilevato: si
  integra fino a t = 162 s e basta. Non ci sono event function in `main.m`.

**Perche' nessuna event function?** Perche' i due switch di fase (t = 5 s,
t = 15 s) sono *deterministici e noti a priori*: non c'e' nulla da localizzare. Le
event function servono quando l'istante dell'evento e' esso stesso incognito
(impatto al suolo, esaurimento propellente, attraversamento di una quota).

**Cosa succede allora agli switch?** Va guardato in faccia:
- A t = 5 s (fine ascesa verticale) il RHS e' in realta' **continuo**: la legge di
  pitchover a riga 435 vale `gT = 90 - 0.05*(t - 5)`, che per t = 5 da'
  esattamente 90 deg, quindi `Tu = Tmag*sin(90) = Tmag`, `Tv = 0` -- cioe' la
  stessa cosa della fase 1. Lo switch e' invisibile all'integratore.
- A t = 15 s (fine pitchover) il RHS e' **discontinuo**: la direzione di spinta
  passa di colpo da elevazione 89.5 deg / azimuth Est alla direzione della
  velocita' relativa, che in quell'istante ha un'elevazione simile ma non
  identica. `ode45` non sa che c'e' una discontinuita': semplicemente il suo
  stimatore d'errore locale esplode sul passo che la attraversa, il passo viene
  rifiutato e dimezzato ripetutamente finche' non rientra in tolleranza. Il
  risultato e' corretto ma pagato con qualche decina di passi buttati. La cura
  "pulita" sarebbe spezzare l'integrazione in due chiamate `ode45` con restart a
  t = 15 s.

> **Possibile domanda d'esame** -- `ode45` e' a passo adattivo: perche' allora
> imporre `MaxStep = 1`?
> *Risposta:* Perche' l'adattivita' controlla l'*errore*, non la *risoluzione
> dell'output*. Con tolleranze cosi' strette e una dinamica dolce, `ode45`
> potrebbe comunque scegliere passi lunghi in tratti "facili", producendo una
> griglia rada. Ma max-Q e Mach 1 non vengono localizzati con event function:
> sono estratti con `max(qdyn)` e `find(Mach>=1,1)` *sui campioni*. Una griglia
> rada sposterebbe l'istante riportato di max-Q anche di parecchi secondi.
> `MaxStep = 1` e' quindi un vincolo di *post-processing*, non di accuratezza.

---

## `POST-PROCESSING` (righe 100-188)

### Estrazione e quantita' derivate (righe 104-132)

- Righe 105-111: unpacking delle sette colonne di `Y`.
- Riga 114: `h = r - RE`. Quota su Terra sferica.
- Riga 117: `Vmag` = modulo della velocita' **inerziale**.
- Righe 120-123: velocita' **relativa** all'atmosfera. Qui c'e' della fisica:

      urel = u
      vrel = v - omegaE*r*cos(phi)
      wrel = w

  L'atmosfera e' assunta **co-rotante rigidamente con la Terra** (nessun vento).
  Un punto solidale alla Terra a raggio r e latitudine phi ha velocita' inerziale
  puramente verso Est di modulo `omegaE*r*cos(phi)` (e' la velocita' del moto
  circolare attorno all'asse polare, il cui raggio e' `r*cos(phi)`). Sottraendola
  dalla componente Est della velocita' inerziale si ottiene la velocita' rispetto
  all'aria. Le componenti Up e Nord **non** sono toccate, perche' `omega x r` non
  ha componenti in quelle direzioni.
  Verifica di coerenza: a t = 0, `v = omegaE*RE*cos(lat0)` e `r = RE`,
  `phi = lat0`, quindi `vrel(0) = 0` **esattamente**. Il razzo fermo sulla rampa ha
  velocita' relativa nulla. E' proprio questo che rende necessarie le guardie
  `if Vrel > 1e-10` nella `eom`.
- Riga 126: `rho_traj = rho0*exp(-h/Hscale)`. Atmosfera esponenziale isoterma.
  (In post-processing la densita' si chiama `rho_traj`; il nome `rho` esiste solo
  come variabile locale dentro `eom`, riga 399.)
- Riga 129: `qdyn = 0.5*rho_traj.*Vrel.^2`. **Con la velocita' relativa, non con quella
  inerziale.** E' l'unica scelta fisicamente sensata: la pressione dinamica e'
  quella che il veicolo *sente*, cioe' rispetto all'aria. Usare `Vmag` gonfierebbe
  q al lift-off (dove Vmag = 408 m/s ma il razzo e' fermo rispetto all'aria) --
  errore da matita blu.
- Riga 132: `Mach = Vrel / a_sound`, con `a_sound` costante (vedi sopra).

**La struttura di max-Q.** q = 0.5*rho0*exp(-h/H)*Vrel^2 e' il prodotto di un
fattore che *crolla* esponenzialmente con la quota e di uno che *cresce* col
quadrato della velocita'. Derivando rispetto al tempo e ponendo a zero:

    d/dt[ln q] = -(1/H)*dh/dt + 2*(1/Vrel)*dVrel/dt = 0
    =>   dVrel/dt / Vrel  =  u / (2*H)

cioe' il max-Q cade nell'istante in cui il **rateo relativo di accelerazione**
uguaglia **meta' del rateo relativo di caduta della densita'**. Prima di quel
punto vince V^2, dopo vince exp(-h/H). Nella run nominale (README) esce
t = 74.9 s, h = 12.78 km, q = 29.5 kPa.

**Perche' Mach 1 (61.8 s) viene PRIMA di max-Q (74.9 s)?** Sono due condizioni
diverse: Mach 1 e' una soglia su Vrel (una velocita' fissata, 340.3 m/s), max-Q e'
un massimo di un prodotto. Il veicolo raggiunge 340 m/s molto prima che il
prodotto rho*V^2 giri. Non c'e' nessun legame causale fra i due eventi -- capita in
quest'ordine per qualunque lanciatore con questo profilo, ed e' il motivo per cui
il transonico e il max-Q sono due criticita' *distinte* nel progetto strutturale.

### Angoli di spinta e aerodinamico (righe 134-155)

```matlab
Vh = sqrt(vrel(k)^2 + wrel(k)^2);
if Vrel(k) > 1e-6
    gammaA(k) = atan2d(urel(k), Vh);
else
    gammaA(k) = 90;
end
```

- Righe 140-145: `gammaA` e' l'**elevazione della velocita' relativa** sul piano
  orizzontale locale: `atan2d(componente Up, modulo della componente
  orizzontale)`. 90 deg = salita verticale, 0 deg = volo orizzontale.
  L'uso di `atan2d` invece di `atand(u/Vh)` e' corretto perche' gestisce
  `Vh -> 0` (salita verticale pura) senza divisione per zero.
- Righe 141-145: la guardia sul caso `Vrel ~ 0` (che, come visto, e' esattamente
  il caso a t = 0) assegna convenzionalmente 90 deg. E' una convenzione, non un
  risultato: a velocita' relativa nulla l'angolo di traiettoria e' indefinito.
- Righe 148-154: `gammaT` viene **ricostruito replicando la logica di fase della
  `eom`**, non estratto dall'integrazione.

**Attenzione, punto onesto e importante.** Alla riga 153, in fase 3,
`gammaT(k) = gammaA(k)` -- cioe' l'angolo di spinta e' *assegnato uguale* a quello
aerodinamico. Nella Figura 7 le due curve coincidono perfettamente dopo t = 15 s,
ma **questo non e' una verifica del gravity turn: e' una tautologia**. Le curve
coincidono perche' il codice le ha rese uguali per costruzione in
post-processing. E' vero che la `eom` implementa davvero la spinta lungo Vrel in
fase 3 (righe 443-449), quindi il grafico non *mente*; semplicemente non
dimostra nulla. Se si volesse una vera verifica bisognerebbe ricalcolare
l'elevazione del vettore spinta dalla `eom` in modo indipendente e confrontarla.

C'e' un secondo aspetto notevole: fra t = 0 e t = 15 s, `gammaT != gammaA`. Il
veicolo vola quindi con un **angolo d'attacco non nullo** durante l'ascesa
verticale e il pitchover (banale al decollo, dove Vrel ~ 0 e alpha e'
indefinito, ma reale durante il pitchover). Il modello pero' **non produce
nessuna forza di portanza o laterale**: c'e' solo drag lungo -Vrel. E' una
semplificazione consapevole del 3-DoF a punto materiale. Il "zero-lift gravity
turn" e' letteralmente a alpha = 0 solo in fase 3.

### Eventi di missione (righe 157-161)

- Righe 158-159: `ip1`, `ip2` = indici dei campioni *piu' vicini* a t = 5 e
  t = 15 s. Non esatti: e' proprio per questo che `MaxStep = 1` serve.
- Riga 160: `im1 = find(Mach >= 1, 1, 'first')` -- **primo campione** oltre Mach 1,
  non interpolazione. L'istante riportato ha quindi un'incertezza pari al passo
  locale (<= 1 s). Restituisce `[]` se il veicolo non supera Mach 1: da qui tutte
  le guardie `if ~isempty(im1)` nei plot.
- Riga 161: `[qmax, imQ] = max(qdyn)` -- massimo **sui campioni**, non del
  polinomio interpolante. Il valore di picco e' quindi leggermente sottostimato
  (il vero massimo cade fra due campioni). Con `MaxStep = 1` s e una curva q(t)
  molto piatta attorno al picco, l'errore e' minimo, ma esiste.

### Da ECI a ENU locale (righe 163-184)

- Righe 165-167: coordinate cartesiane inerziali dalla terna sferica:

      X = r*cos(phi)*cos(theta)
      Y = r*cos(phi)*sin(theta)
      Z = r*sin(phi)

  Notare che con `phi` = declinazione (non colatitudine) compare `cos(phi)` nelle
  componenti equatoriali e `sin(phi)` in quella polare. Se `phi` fosse la
  colatitudine sarebbero scambiati: e' l'errore di segno/funzione piu' comune in
  questi codici.
- Righe 175-177: la matrice `R_enu` ha come **righe** i versori East, North, Up
  del sito di lancio espressi in ECI:

      e_East  = [-sin(theta0),            cos(theta0),           0          ]
      e_North = [-sin(phi0)*cos(theta0), -sin(phi0)*sin(theta0), cos(phi0)  ]
      e_Up    = [ cos(phi0)*cos(theta0),  cos(phi0)*sin(theta0), sin(phi0)  ]

  Moltiplicare per `R_enu` proietta un vettore ECI sulla terna locale: e'
  esattamente la definizione di matrice di rotazione costruita per righe di
  versori. Si verifica facilmente che e_East x e_North = e_Up (terna destrorsa).
- Righe 179-180: si sottrae la posizione del sito di lancio **prima** di ruotare --
  cioe' e' una trasformazione affine (traslazione + rotazione), non una semplice
  rotazione.
- **Sottigliezza**: `R_enu` e' costruita con `theta0`, cioe' con la terna locale
  **congelata all'istante del lancio nel frame inerziale**. Il sito di lancio,
  pero', *ruota* con la Terra. Quindi la Figura 1 non e' "la traiettoria vista dal
  sito di lancio" in senso ECEF: e' la traiettoria inerziale proiettata su una
  terna inerziale fissa che *coincideva* con la terna locale a t = 0. Su 162 s la
  Terra ruota di ~0.68 deg, quindi la differenza e' piccola ma non nulla. Il
  codice non esplicita questa scelta.

### Ground track (righe 186-188)

- Riga 187: `lon_ground = rad2deg(theta - omegaE*t)`. Qui invece la rotazione
  terrestre e' tolta correttamente: si passa da ascensione retta (inerziale) a
  longitudine (ECEF). E' la conferma definitiva che `theta` e' inerziale.
- Riga 188: la latitudine `phi` non ha bisogno di correzioni (la rotazione e'
  attorno all'asse polare, non cambia la latitudine).

---

## `CONSOLE SUMMARY` (righe 190-212)

Blocco di `fprintf`. Non fa fisica, ma i suoi output sono gli stessi numeri
riportati nel README e su cui si appoggiano le asserzioni "di banda" del test
(quota finale fra 50 e 120 km). Riga 207: la guardia `if ~isempty(im1)` evita di
stampare Mach 1 se non c'e' stato.

---

## `PLOTS` (righe 214-350)

Otto figure. Meritano un commento solo tre cose.

- Righe 219-222: la palette dei quattro eventi (fine verticale / fine pitchover /
  Mach 1 / max-Q) e' definita una volta e riusata su tutte le figure: gli eventi
  hanno sempre lo stesso colore e lo stesso marker in ogni grafico. E' una buona
  pratica di leggibilita' e la traccia lo richiede esplicitamente ("mark the key
  mission events").
- `HandleVisibility','off'` sulle `xline` degli switch di fase: le linee
  tratteggiate verticali compaiono nel grafico ma non nella legenda. Senza questo,
  ogni legenda avrebbe due voci inutili.
- Figura 7 (righe 326-333): come discusso sopra, la sovrapposizione di gammaT e
  gammaA dopo t = 15 s e' costruita, non misurata.

---

## `EXPORT FIGURES TO figures/` (righe 352-369)

```matlab
slugify = @(s) lower(regexprep(s, '[^a-zA-Z0-9]+', '_'));
fig_handles = findobj(groot, 'Type', 'figure');
for kk = 1:numel(fig_handles)
    try
        theme(fig_handles(kk), 'light');
        drawnow;
    catch
        fig_handles(kk).Color = 'w';
    end
    ...
end
```

- Riga 355: la cartella di destinazione e' risolta da `mfilename('fullpath')`, non
  dal working directory: lo script si puo' lanciare da qualunque cartella.
- Riga 358: `slugify` trasforma il `Name` della figura in un nome file
  (`'3D Trajectory'` -> `3d_trajectory.png`). E' il motivo per cui i nomi delle
  figure alle righe 225, 249, ... non sono cosmetici: **sono i nomi dei file** e
  quindi i link del README.
- Righe 361-366: `theme(fig, 'light')` forza il tema chiaro anche se il desktop
  MATLAB e' in dark mode -- altrimenti i PNG uscirebbero con sfondo nero e
  finirebbero cosi' nel report LaTeX. Il `try/catch` degrada su
  `fig.Color = 'w'` per versioni pre-R2025a, dove `theme` non esiste.
- **Onesta'**: `findobj(groot, 'Type','figure')` esporta **tutte** le figure aperte
  nella sessione, non solo le otto di questo script. In pratica non e' un problema
  perche' la riga 6 fa `close all`, ma e' una dipendenza implicita dallo stato
  della sessione. Se qualcuno rimuovesse il `close all`, l'export sputerebbe fuori
  anche figure estranee.

---

## `eom` -- equazioni del moto (righe 375-466)

E' il cuore del file. La firma:

```matlab
function dydt = eom(t, y, par)
% State vector: y = [r, theta, phi, u, v, w, m]
```

- Riga 375: funzione pura, chiamata ~10^4-10^5 volte da `ode45`. Coerentemente con
  la convenzione della repo (`CLAUDE.md`), **non ha blocco `arguments`**: e' una
  hot loop, la validazione costerebbe piu' del calcolo. Non e' una svista.

### Derivazione delle equazioni (righe 452-460)

Questo e' il punto su cui vale la pena spendere tempo, perche' e' *l'unica vera
matematica del file* e all'orale e' la domanda quasi certa.

**Setup.** La posizione e' `R = r * e_r`, dove la terna locale (in **frame
inerziale**, con `theta` = ascensione retta e `phi` = declinazione) e':

    e_r = ( cos(phi)*cos(theta),  cos(phi)*sin(theta),  sin(phi) )   [Up]
    e_e = (-sin(theta),           cos(theta),           0        )   [East]
    e_n = (-sin(phi)*cos(theta), -sin(phi)*sin(theta),  cos(phi) )   [North]

Si verifica che e_e x e_n = e_r: la terna (East, North, Up) e' destrorsa.

**Cinematica (righe 452-455).** Derivando `R = r*e_r` nel tempo e usando

    de_r/dt = cos(phi)*theta_dot * e_e  +  phi_dot * e_n

si ottiene

    V = r_dot * e_r  +  r*cos(phi)*theta_dot * e_e  +  r*phi_dot * e_n

Confrontando con `V = u*e_r + v*e_e + w*e_n` si legge subito l'identificazione
delle componenti UEN, e quindi le tre equazioni cinematiche del codice:

    u = r_dot                    ->  dr     = u
    v = r*cos(phi)*theta_dot     ->  dtheta = v / (r*cos(phi))
    w = r*phi_dot                ->  dphi   = w / r

Sono **esattamente** le righe 453-455. Notare che `dtheta` diverge ai poli
(`cos(phi) -> 0`): e' una singolarita' **della parametrizzazione**, non della
fisica. Per un lancio da 28.57 deg non e' un problema, ma per una missione polare
lo sarebbe (e sarebbe uno degli argomenti per passare a una formulazione
cartesiana o a quaternioni).

**Dinamica (righe 457-460).** Qui viene il bello. Derivando
`V = u*e_r + v*e_e + w*e_n` bisogna derivare **anche i versori**, che ruotano
perche' il veicolo si muove. Le derivate della terna sono:

    de_r/dt =  cos(phi)*theta_dot * e_e  +  phi_dot * e_n
    de_e/dt = -cos(phi)*theta_dot * e_r  +  sin(phi)*theta_dot * e_n
    de_n/dt = -phi_dot * e_r             -  sin(phi)*theta_dot * e_e

Quindi

    a = u_dot*e_r + v_dot*e_e + w_dot*e_n
        + u*de_r/dt + v*de_e/dt + w*de_n/dt

Raccogliendo componente per componente e sostituendo `theta_dot = v/(r*cos(phi))`
e `phi_dot = w/r`:

    a_Up    = u_dot - (v^2 + w^2)/r
    a_East  = v_dot + ( u*v - v*w*tan(phi) ) / r
    a_North = w_dot + ( u*w + v^2*tan(phi) ) / r

Invertendo (a = F/m) si ottengono **esattamente** le righe 458-460:

    du = (v^2 + w^2)/r         + g_u + (Tu + Du)/m
    dv = (-u*v + v*w*tan(phi))/r     + (Tv + Dv)/m
    dw = (-u*w - v^2*tan(phi))/r     + (Tw + Dw)/m

**Interpretazione dei termini di trasporto** (la domanda d'esame):

- `(v^2 + w^2)/r` in `du`: e' il termine **centrifugo/di curvatura**. Il veicolo
  si muove tangenzialmente con velocita' orizzontale sqrt(v^2+w^2) lungo una
  superficie di raggio r: per restare su quella superficie servirebbe
  un'accelerazione centripeta `V_h^2/r` verso il basso. Nella forma risolta per
  `u_dot`, compare col segno + come *accelerazione apparente verso l'alto*.
  Fisicamente e' **il termine che alleggerisce il veicolo man mano che accelera in
  orizzontale** -- al limite `V_h^2/r = mu/r^2` si ha l'orbita circolare e
  `du = 0` anche a spinta nulla. Se lo si togliesse, il razzo non potrebbe *mai*
  andare in orbita: ricadrebbe sempre.
- `-u*v/r` in `dv` e `-u*w/r` in `dw`: conservazione del **momento angolare**. Se
  il veicolo sale (u > 0) e non ci sono coppie, la velocita' tangenziale deve
  diminuire, perche' r*v ~ cost. E' lo stesso effetto per cui una pattinatrice
  rallenta aprendo le braccia.
- `+v*w*tan(phi)/r` in `dv` e `-v^2*tan(phi)/r` in `dw`: accoppiamento
  **Est-Nord dovuto alla convergenza dei meridiani**. Andando verso i poli i
  meridiani si stringono, quindi muoversi verso Nord con una componente Est fa
  ruotare il vettore velocita'. Sono i termini che degenerano ai poli
  (`tan(phi) -> inf`), stessa singolarita' di prima.

**Il punto piu' importante di tutti: dov'e' Coriolis?** Non c'e', e **non deve
esserci**. Le equazioni sono scritte nel **frame inerziale**: `theta` e'
ascensione retta, `u, v, w` sono componenti della velocita' *inerziale*.
Coriolis e la centrifuga *planetaria* sono accelerazioni **apparenti**, che
compaiono solo se si scrive F = m*a in un frame **rotante** (ECEF). Qui non lo si
fa. I termini di trasporto sopra **non sono** Coriolis: nascono dalla derivazione
dei versori di una terna che ruota perche' *il veicolo si sposta*, non perche' *la
Terra gira* -- infatti `omegaE` non compare in nessuno di essi.

La rotazione terrestre entra in `main.m` in **tre** punti soltanto, tutti espliciti:
1. la condizione iniziale `v0 = omegaE*RE*cos(lat0)` (riga 67);
2. la velocita' relativa usata per drag e direzione di spinta (righe 403-407);
3. la conversione a longitudine di ground track (riga 187).

Se si riscrivesse tutto in ECEF si otterrebbero le *stesse traiettorie*, ma le
equazioni avrebbero termini `-2*omega x V_rel` (Coriolis) e
`-omega x (omega x R)` (centrifuga planetaria) espliciti, e le condizioni
iniziali sarebbero a velocita' nulla. E' pura contabilita': la fisica e' identica.

### Atmosfera e spinta (righe 398-423)

```matlab
rho  = par.rho0 * exp(-alt / par.H);
patm = par.p0   * exp(-alt / par.H);
...
Tmag = par.Tvac - patm * par.Aex;
```

- Righe 399-400: densita' **e** pressione decadono con **la stessa** scala H. Non
  e' un'assunzione gratuita: con gas ideale p = rho*R*T, se p e rho hanno lo
  stesso H allora T = cost. E' precisamente l'ipotesi isoterma. E' anche cio' che
  giustifica `a_sound` costante.
  La derivazione: equilibrio idrostatico dp/dh = -rho*g, con p = rho*R*T e T, g
  costanti, da' dp/p = -dh/(R*T/g) = -dh/H, quindi p = p0*exp(-h/H).
- Riga 423: **la correzione di contropressione**, l'unico pezzo di
  propulsione "vera" del file. L'equazione della spinta di un ugello e':

      T = mdot * u_e + (p_e - p_amb) * A_e

  dove u_e e' la velocita' di efflusso, p_e la pressione al piano d'uscita, A_e
  l'area d'uscita. In vuoto (p_amb = 0):

      T_vac = mdot * u_e + p_e * A_e

  Sottraendo membro a membro:

      T(h) = T_vac - p_amb(h) * A_e

  che e' esattamente la riga 423. Il fatto notevole e' che **non serve conoscere
  ne' u_e ne' p_e separatamente**: bastano T_vac, A_e e la pressione ambiente.
  Il termine `mdot*u_e` (spinta di quantita' di moto) e' invariante con la quota,
  quindi tutta la dipendenza da h e' nel termine di pressione.
  Numericamente: al livello del mare `T_SL = 8.227 - 0.10133*11.039 = 7.108` MN,
  cioe' **il 13.6% di spinta in meno** rispetto al vuoto. La spinta cresce poi
  monotonamente durante la salita fino a T_vac.
- Riga 463: `dm = -par.Qdot`, **costante**. La portata non dipende dalla quota:
  e' l'ugello che e' "choked", la portata e' fissata dalla gola e dalla pressione
  in camera, non dall'ambiente. Solo la *spinta* varia con h, non il consumo.
  Conseguenza diretta e verificabile: `m(t) = m0 - Qdot*t` esattamente, che e'
  proprio cio' che il test `testPropellantBookkeeping` sfrutta come oracolo
  analitico.

### Drag (righe 412-420)

```matlab
if Vrel > 1e-10
    Dcoeff = -0.5 * rho * par.Sref * par.CD * Vrel;
    Du = Dcoeff * ur;
    ...
```

- Il drag e' scritto in **forma vettoriale**:

      D_vec = -0.5*rho*S*CD*|V_rel| * V_rel_vec

  Il modulo e' `0.5*rho*S*CD*|V_rel|^2` (la formula classica) e la direzione e'
  automaticamente `-V_rel/|V_rel|`. Scriverlo cosi' evita di calcolare
  esplicitamente il versore (e quindi una divisione) ed e' numericamente piu'
  robusto: `Dcoeff` va a zero linearmente con Vrel, quindi anche senza la guardia
  il drag sarebbe corretto a Vrel = 0.
- La guardia `Vrel > 1e-10` e' quindi **ridondante per il drag** ma indispensabile
  per la fase 3 (righe 443-449), dove si divide per `Vrel`. Ha senso averla lo
  stesso, per simmetria.
- **Non c'e' portanza**, non c'e' forza laterale, non c'e' dipendenza di CD da
  Mach. Il modello aerodinamico e' "drag-only, CD costante".

### Logica di fase e legge di pitchover (righe 425-450)

- **Fase 1** (righe 425-429), t <= 5 s: `Tu = Tmag`, `Tv = Tw = 0`. Spinta
  puramente radiale. Serve a guadagnare quota e velocita' prima di iniziare a
  ruotare, in una regione dove la pressione dinamica e' bassa e il veicolo e' ancora
  lento (quindi poco controllabile aerodinamicamente).
- **Fase 2** (righe 431-439), pitchover:

  ```matlab
  gT  = deg2rad(90 - 0.05 * (t - par.t1));
  psi = deg2rad(90);
  Tu = Tmag * sin(gT);
  Tv = Tmag * cos(gT) * sin(psi);
  Tw = Tmag * cos(gT) * cos(psi);
  ```

  L'elevazione scende **linearmente** da 90 deg a 89.5 deg in 10 s: un rateo di
  **0.05 deg/s**, per un totale di mezzo grado. E' un "kick" minuscolo, ed e'
  intenzionale: serve solo a *innescare* il gravity turn, non a eseguirlo.
  L'azimuth `psi = 90 deg` significa **Est**, il che spiega perche' il ground track
  si incurva verso Est e perche' il veicolo sfrutta la rotazione terrestre.
  Le formule `Tu = T*sin(gT)`, `Tv = T*cos(gT)*sin(psi)`, `Tw = T*cos(gT)*cos(psi)`
  sono la decomposizione standard di un versore dato elevazione (dal piano
  orizzontale) e azimuth (misurato **da Nord verso Est**, visto che `psi = 0` da'
  spinta puramente Nord e `psi = 90` puramente Est).

  **Onesta' numerica**: `cos(deg2rad(90))` in floating point non fa 0 ma
  6.1e-17. Quindi `Tw` in fase 2 non e' esattamente zero: vale
  ~6e-17*Tmag*cos(gT), cioe' al massimo ~4e-12 N (caso peggiore a fine
  pitchover, dove cos(gT) = cos(89.5 deg) = 8.7e-3 e Tmag ~ 7.1 MN; a t = 5 s,
  con gT = 90 deg, il fattore 6.1e-17 compare due volte e `Tw` e' ancora piu'
  piccolo). Fisicamente irrilevante, ma e'
  una differenza *reale* rispetto a `main2.m`, che pone `Tw = 0` esattamente
  (riga 550 di quel file). Se qualcuno chiedesse perche' i due script non danno
  bit-per-bit lo stesso risultato, questo e' uno dei motivi (minore).

- **Fase 3** (righe 441-450), gravity turn:

  ```matlab
  Tu = Tmag * ur / Vrel;
  Tv = Tmag * vr / Vrel;
  Tw = Tmag * wr / Vrel;
  ```

  La spinta e' allineata al **versore della velocita' relativa**. Questo e' il
  gravity turn (o "zero-lift turn"): dato che l'assetto segue la velocita'
  relativa, l'**angolo d'attacco e' identicamente nullo**, e quindi (in un modello
  con portanza) non ci sarebbe nessuna forza aerodinamica normale. La rotazione
  della traiettoria e' allora prodotta **unicamente dalla componente di gravita'
  normale alla velocita'** -- da cui il nome.

  La ragione ingegneristica e' strutturale, non di prestazione: volare ad alpha = 0
  attraverso la regione di alta pressione dinamica azzera i carichi laterali
  (`q*alpha`), che sono il vincolo dimensionante del corpo del lanciatore. Un
  gravity turn e' *subottimale* dal punto di vista della perdita di gravita' -- un
  profilo di pitch ottimo (HM1) fa meglio -- ma e' quello che la struttura
  sopporta.

  Nota: dopo l'ingresso in fase 3 il veicolo continua a inclinarsi da solo perche'
  la gravita' ha una componente normale a V; ma il **pitchover di mezzo grado
  della fase 2 e' cio' che rende non nulla quella componente**. Senza il kick, la
  velocita' resterebbe esattamente radiale, la gravita' sarebbe esattamente
  antiparallela a V, la componente normale sarebbe nulla e il veicolo salirebbe
  verticalmente per sempre. **Il gravity turn e' un equilibrio instabile che va
  innescato.**

> **Possibile domanda d'esame** -- Il pitchover dura 10 s e ruota il vettore spinta
> di appena 0.5 deg. Come fa un'inclinazione cosi' piccola a produrre una
> traiettoria che a MECO e' quasi orizzontale?
> *Risposta:* Perche' il pitchover non deve *fare* la rotazione, deve solo
> *innescarla*. Dopo il kick, la velocita' non e' piu' esattamente radiale, quindi
> la gravita' acquista una componente normale a V. Da quel momento, in fase 3, la
> spinta segue V e la gravita' normale ruota V continuamente e in modo
> auto-sostenuto per i restanti 147 s. La rotazione totale e' l'integrale di
> -g*cos(gamma)/V nel tempo, che e' grande; il kick iniziale ne fissa solo il
> segno e l'innesco. Con 0 deg di kick, la componente normale di gravita' sarebbe
> identicamente nulla e il razzo salirebbe verticalmente per sempre.

> **Possibile domanda d'esame** -- Perche' in fase 3 la spinta segue la velocita'
> **relativa** e non quella inerziale?
> *Risposta:* Perche' l'obiettivo del gravity turn e' azzerare l'**angolo
> d'attacco**, che e' un angolo *aerodinamico*: e' definito rispetto al flusso
> d'aria, cioe' rispetto alla velocita' relativa all'atmosfera. Allineare la
> spinta a V_inerziale lascerebbe un alpha residuo pari all'angolo fra V_in e
> V_rel -- che al lift-off e' 90 deg (V_rel = 0!) e resta di parecchi gradi nella
> fase densa. Sarebbe autolesionistico proprio dove i carichi q*alpha contano.

---

## Possibili domande d'esame

**D: Le equazioni del moto sono scritte in un frame inerziale o rotante? Come si
riconosce dal codice, e dove finisce la rotazione terrestre?**
R: Inerziale. Si riconosce da tre indizi: (1) nelle equazioni dinamiche (righe
458-460) non compare mai `omegaE`; (2) `theta` e' etichettata "right ascension" e
la longitudine di ground track si ottiene sottraendo `omegaE*t` (riga 187); (3) la
condizione iniziale ha `v0 = omegaE*RE*cos(lat0)` diversa da zero (riga 67), cosa
che in ECEF sarebbe nulla. I termini `(v^2+w^2)/r`, `-u*v/r`, `v*w*tan(phi)/r`
**non sono** Coriolis: sono termini di trasporto che nascono dal derivare i
versori della terna UEN, che ruota perche' il veicolo si sposta sulla sfera. La
rotazione terrestre entra solo in tre punti espliciti: condizione iniziale,
velocita' relativa (per drag e direzione di spinta), e conversione a longitudine.

**D: Da dove viene `T = T_vac - p_atm*A_ex`? E perche' la portata di massa non
dipende dalla quota?**
R: Dall'equazione della spinta di un ugello, T = mdot*u_e + (p_e - p_amb)*A_e.
Scrivendola in vuoto (T_vac = mdot*u_e + p_e*A_e) e sottraendo, tutti i termini
incogniti (u_e, p_e) si cancellano e resta T(h) = T_vac - p_amb(h)*A_e. Il
termine di quantita' di moto e' indipendente dalla quota, quindi tutta la
variazione e' nel termine di pressione: al livello del mare la perdita e'
p0*Aex = 1.12 MN, cioe' il 13.6% (7.108 MN contro 8.227 MN). La portata invece e'
costante (riga 463) perche' l'ugello e' in condizioni critiche: mdot e' fissata
dalla gola e dalla pressione in camera, non dall'ambiente esterno. Solo la spinta
"sente" l'atmosfera, il consumo no.

**D: Mach 1 a 61.8 s, max-Q a 74.9 s. Perche' in quest'ordine? Sarebbe potuto
essere il contrario?**
R: Sono due condizioni indipendenti. Mach 1 e' una **soglia** su una velocita'
fissa (340.3 m/s con `a_sound` costante). Max-Q e' lo **stazionario** del prodotto
rho(h)*V^2, e cade dove dV/dt / V = u/(2*H), cioe' dove la crescita relativa di
velocita' non riesce piu' a compensare meta' della caduta relativa di densita'.
Con un lanciatore che accelera rapidamente in un'atmosfera con H = 8 km, la
soglia dei 340 m/s viene raggiunta molto prima che il prodotto giri: l'ordine
Mach 1 -> max-Q e' quindi robusto per qualunque lanciatore di questa classe. In
linea di principio un veicolo che accelerasse molto lentamente in quota potrebbe
invertirlo, ma non e' il caso di un primo stadio.

**D: Il modello usa `a_sound` costante e `CD` costante. Cosa si perde?**
R: Due cose diverse. `a_sound` costante discende coerentemente dall'ipotesi di
atmosfera isoterma (se T e' costante, a lo e'), ma quell'ipotesi e' falsa: nella
troposfera reale T scende e a arriva a ~295 m/s a 11 km. Quindi **il Mach in
quota e' sottostimato di circa il 13-15%** e l'attraversamento di Mach 1 e'
riportato piu' tardi del vero. `CD = 0.329` costante e' invece un errore di
modello puro: il coefficiente di resistenza di un lanciatore ha un picco
transonico marcato. Il codice quindi **sottostima il drag proprio attorno a
max-Q**, che e' dove pesa di piu': la quota e la velocita' a MECO risultano
leggermente ottimistiche.

**D: Perche' non e' stata usata una event function di `ode45` per gli switch di
fase, mentre `MaxStep = 1` c'e'?**
R: Le event function servono a localizzare istanti **incogniti** (impatto,
raggiungimento di una quota). Qui gli switch sono a t = 5 s e t = 15 s, noti a
priori: non c'e' nulla da cercare. `MaxStep = 1` risponde a un problema diverso:
max-Q e Mach 1 sono estratti in post-processing con `max()` e `find()` **sui
campioni** (righe 160-161), quindi la loro accuratezza dipende dalla densita'
della griglia di output, non dalla tolleranza dell'integratore. Detto questo, il
RHS **e' realmente discontinuo a t = 15 s** (la direzione di spinta salta dalla
legge di pitchover alla direzione di V_rel) e `ode45` la attraversa a forza di
passi rifiutati. Una soluzione piu' pulita sarebbe stata spezzare l'integrazione
in due chiamate con restart a t = 15 s.

**D: La Figura 7 mostra gammaT che coincide con gammaA dopo t = 15 s. E' una
verifica che il gravity turn e' implementato bene?**
R: **No, e' una tautologia.** In post-processing, alla riga 153, il codice
*assegna* `gammaT(k) = gammaA(k)` per t > 15 s. Le due curve coincidono per
costruzione, non come risultato dell'integrazione. E' vero che la `eom` implementa
davvero la spinta lungo V_rel (righe 443-449), quindi il grafico e' coerente con
il modello -- ma non lo verifica. Una verifica vera richiederebbe di ricalcolare
l'elevazione del vettore spinta in modo indipendente dalla dinamica integrata e
confrontarla.

**D: Cosa succederebbe se si togliesse il termine `(v^2 + w^2)/r` dall'equazione
di `du`?**
R: Il veicolo non potrebbe **mai** andare in orbita. Quel termine e'
l'accelerazione centrifuga associata al moto orizzontale su una superficie di
raggio r: e' cio' che alleggerisce progressivamente il veicolo man mano che
guadagna velocita' orizzontale. La condizione di orbita circolare e' esattamente
`(v^2+w^2)/r = mu/r^2`, cioe' il termine di trasporto che cancella la gravita'.
Rimuovendolo, `du` resterebbe negativa per sempre a spinta nulla e qualunque
traiettoria ricadrebbe. E' il termine che rende "vera" la meccanica orbitale in
queste coordinate.
