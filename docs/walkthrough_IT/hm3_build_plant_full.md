# HM3/build_plant_full.m

## Ruolo del file nel progetto

E' il plant di **Task 2 e Task 3**: il modello rigido di `build_plant_rigid.m`
esteso con il **primo modo di flessione** del corpo (2 stati in piu', `eta` e
`etadot`) e con l'**osservazione INS contaminata** dalla deformazione elastica
(Eq. 2 della traccia). Sei stati `[z, zdot, theta, thetadot, eta, etadot]`, due
ingressi `[delta, alpha_w]`, le stesse 7 uscite del plant rigido -- interfaccia
identica, cosi' che `assemble_loop`, `design_controller`, `classify_margins` e
`simulate_gust_response` funzionino senza modifiche su entrambi.

E' il file che **crea il problema che Task 2 deve risolvere**. Con il plant rigido
il loop e' pulito; sostituendo questo plant, il guadagno di anello alla frequenza
del modo di bending sale a **+29 dB** e il loop diventa instabile. Da li' nascono
il trade fra quattro filtri, il notch profondo scelto e il ri-tuning del PD.

Chi lo chiama: `main_task2.m` (riga 29), `main_task3.m` (righe 21 e 50, per le
quattro corner), `main_montecarlo.m`, `make_pz_figures.m`, `init_simulink_hm3.m`,
`run_simulink_closed_loop.m` (riga 40), `LTV_FULL_ASCENT/main_flex.m`, e i test. Il ramo `meas = 'true'` **non e' usato da
nessun main**: e' esercitato solo da `tests/hm3PlantTest.m` (riga 92) e serve come
ablazione didattica -- ed e' proprio quell'ablazione a dimostrare la tesi centrale
di Task 2 (vedi sotto).

---

## Firma, `arguments`, contratto (righe 1-16)

```matlab
function G = build_plant_full(p, meas)
arguments
    p (1,1) struct
    meas {mustBeTextScalar} = 'ins'   % 'ins' | 'true', case-insensitive
end
```

- Riga 1: due argomenti posizionali. `meas` di default `'ins'`, cioe' **il caso
  realistico** (la piattaforma inerziale sente la deformazione). Chi vuole il
  confronto pulito chiede esplicitamente `'true'`.
- Righe 10-11 (commento): l'header dichiara la tesi del file -- *"INS bending
  contamination (sigma_ins, phi_ins) is what destabilises the loop at wBM and
  motivates the Task-2 notch"*. E' vero, e il codice lo dimostra in modo
  falsificabile: vedi la sezione finale.
- Riga 18: `w = p.wBM; z = p.zBM;` -- attenzione, `z` qui e' lo **smorzamento
  modale**, non lo stato di drift. Shadowing di nome locale, innocuo ma
  potenzialmente confondente in lettura.

---

## La matrice A: le due righe nuove (righe 20-25)

```matlab
A = [0   1            0              0   0      0;
     0   p.a1         p.a1*p.V+p.a4  0   0      0;
     0   0            0              1   0      0;
     0   p.A6/p.V     p.A6           0   0      0;
     0   0            0              0   0      1;
     0   0            0              0  -w^2   -2*z*w];
```

- Righe 20-23: il blocco 4x4 in alto a sinistra e' **identico** a quello di
  `build_plant_rigid.m` (stesse equazioni dei momenti e delle forze laterali,
  stessa convenzione `alpha = theta + zdot/V - alpha_w`).
- Righe 24-25: l'**oscillatore modale**. In forma di equazione differenziale:

      eta_ddot + 2*zeta_BM*omega_BM*etadot + omega_BM^2*eta = f(t)

  cioe' la classica equazione modale disaccoppiata che si ottiene proiettando le
  equazioni della trave elastica sul primo modo. `eta` e' la **coordinata modale
  generalizzata**; qui ha unita' di metri (lo si deduce dalle unita' di
  `sigma_ins` [rad/m] e `phi_tvc` [1/kg], vedi sotto).
  Con `omega_BM = 18.9 rad/s` e `zeta_BM = 0.005` i due autovalori sono

      -zeta*omega +/- j*omega*sqrt(1-zeta^2) = -0.0945 +/- 18.8998j

  (verificato con `pole()`), pinnati dal test `testFullPlantBendingMode`
  (`hm3PlantTest.m` righe 72-80) con tolleranza `1e-9`.

### Il punto chiave: il coupling nella matrice A e' a senso unico

Le colonne 5-6 delle righe rigide (2 e 4) sono **nulle**, e le righe 5-6 hanno
colonne 1-4 **nulle**. Significa:

- la deformazione elastica **non retroagisce** sulla dinamica di corpo rigido (si
  trascurano gli effetti inerziali/aerodinamici del modo sul moto rigido);
- il moto rigido **non eccita** il modo di bending (nessun forcing da `theta`).

L'unico accoppiamento nella dinamica e' l'ingresso `delta` (riga 27), e l'unico
accoppiamento nell'uscita e' la matrice `Cm` (righe 33-36). Quindi **il modo di
bending si chiude in un anello solo passando per il controllore**:

    delta -> eta  (forcing TVC)  ->  theta_m, zdot_m  (contaminazione INS)
          -> controllore -> delta

Non e' un fenomeno strutturale: e' un **loop di controllo indesiderato** creato
dal sensore. Questo e' esattamente il motivo per cui la soluzione e' un filtro
(notch) e non un irrigidimento del corpo.

---

## Le colonne di ingresso (righe 27-28)

```matlab
Bd = [0; p.a3; 0; p.K1; 0; -p.phi_tvc*p.Tc];   % delta column
Bw = [0; -p.a1*p.V; 0; -p.A6; 0; 0];           % alpha_w column
```

- Riga 27: il **forcing del bending da parte del TVC**. La forza laterale
  `Tc*delta` applicata al gimbal eccita il modo con un guadagno pari alla
  pendenza/ampiezza della forma modale in quel punto normalizzata sulla massa
  generalizzata -- che e' il significato di `phi_tvc` [1/kg]. Numericamente:

      phi_tvc * Tc = 4.3095e-5 * 1.5213e6 = 65.56  m/s^2 per rad

  Il **segno negativo** dice che una deflessione positiva flette il corpo nel
  verso di `eta` negativo (coerente con il fatto che la forza al gimbal e'
  opposta al momento generato). Coerenza dimensionale: `[1/kg]*[N]*[rad] =
  [m/s^2]`, quindi `eta` e' in metri e `eta_ddot` in m/s^2. Torna.
- Riga 28: la colonna del vento ha **zeri sulle righe 5-6**: la raffica **non
  eccita direttamente il modo di flessione**. E' un'approssimazione del modello
  della traccia (in un lanciatore vero il carico distribuito di raffica eccita
  eccome i modi elastici). Conseguenza pratica: nelle simulazioni di raffica il
  bending si vede solo *indirettamente*, attraverso il `delta` che il controllore
  comanda. Va detto all'orale se si viene interrogati sui carichi.

> **Possibile domanda d'esame** -- Se `alpha_w` non eccita il bending e il bending
> non retroagisce sul corpo rigido, perche' e' un problema?
> *Risposta:* perche' il problema non e' dinamico ma di **retroazione**. Il TVC
> eccita il modo (`Bd(6)`), l'INS lo misura come se fosse assetto (`Cm(1,5)`), il
> controllore ci reagisce e comanda altro TVC. Se il guadagno complessivo di questo
> giro supera 0 dB con la fase sbagliata alla risonanza, il modo va in
> autoeccitazione: il loop e' instabile anche se la struttura, da sola, e' stabile
> (smorzamento positivo, seppur 0.5 %).

---

## L'osservazione INS: Eq. (2) (righe 30-44)

```matlab
case 'ins'   % Eq. (2): bending leaks into the measurements
    Cm = [0 0 1 0  p.sigma_ins 0;   % theta_m    = theta + sigma*eta
          0 0 0 1  0 p.sigma_ins;   % thetadot_m = thetadot + sigma*etadot
          1 0 0 0 -p.phi_ins 0;     % z_m        = z - phi*eta
          0 1 0 0  0 -p.phi_ins];   % zdot_m     = zdot - phi*etadot
```

- Righe 33-36: e' la **Eq. (2)** della traccia. L'INS e' montato in una stazione
  precisa del corpo e non puo' distinguere fra rotazione rigida e rotazione locale
  dovuta alla flessione:
  - il **giroscopio** misura l'inclinazione locale della struttura, cioe'
    `theta + sigma_ins*eta`, dove `sigma_ins` [rad/m] e' la **pendenza della forma
    modale** alla stazione dell'INS. Con `sigma_ins = 0.178 rad/m`: un metro di
    ampiezza modale generalizzata si presenta al sensore come 0.178 rad = 10.2 deg
    di assetto **falso**;
  - l'**accelerometro/integratore di posizione** misura lo spostamento locale,
    cioe' `z - phi_ins*eta`, con `phi_ins = 0.8` [-] l'**ampiezza** (non la
    pendenza) della forma modale nello stesso punto.
  I segni opposti (`+sigma`, `-phi`) sono quelli della traccia e riflettono la
  geometria della forma modale in quella stazione; il test
  `testInsMeasurementsContaminatedByBending` (righe 82-88) li pinna.
- Righe 37-41: il ramo `'true'` azzera le colonne 5-6 delle misure: retroazione
  sugli stati **veri**, come se esistesse un filtro/osservatore ideale che separa
  rigido ed elastico. Non e' realizzabile, ma e' l'esperimento di controllo che
  isola la causa.
- Righe 42-43: `error('build_plant_full:meas', ...)` con identificatore esplicito,
  testato (`testFullPlantRejectsUnknownMeas`).
- Righe 45-48: le 3 righe di plotting (`theta`, `z`, `zdot` **veri**) sono
  aggiunte sotto in entrambi i rami. Sono i segnali che i main mandano nei grafici:
  quello che il veicolo fa davvero, non quello che il sensore crede.

---

## Da dove vengono i +29 dB (la derivazione che serve all'orale)

Il percorso destabilizzante e' `delta -> eta -> theta_m`. La sua funzione di
trasferimento si legge direttamente dalle matrici:

    theta_m/delta |_bending = sigma_ins * (-phi_tvc*Tc) / (s^2 + 2*zBM*wBM*s + wBM^2)

Alla risonanza (`s = j*wBM`) il denominatore vale `2*j*zBM*wBM^2`, quindi il
modulo del **residuo di bending** e':

    |theta_m/delta|(wBM) = sigma_ins*phi_tvc*Tc / (2*zBM*wBM^2)
                         = 0.178 * 65.56 / (2*0.005*18.9^2)
                         = 11.67 / 3.572 = 3.267   ->  +10.3 dB

A questo si moltiplica il **guadagno del controllore alla stessa frequenza**. Il
PD di Task 1 (`Kp = 1.78`, `Kd = 0.44`) a 18.9 rad/s vale

    |Kp + j*wBM*Kd| = |1.78 + j*8.32| = 8.50   ->  +18.6 dB

perche' a quella frequenza il termine derivativo domina di un fattore ~5. Somma:

    +10.3 dB + 18.6 dB = +28.9 dB

e infatti, valutando il vero guadagno di anello con `assemble_loop` sul plant
completo con attuatore ideale, si ottiene **|L(j*wBM)| ~ 29.0 dB** -- i +29 dB
del README. La derivazione chiude entro ~0.1 dB, ed e' il modo di rispondere alla
domanda "da dove viene quel numero" senza dire "me lo dice MATLAB".

Da qui segue tutto il resto di Task 2: per **gain-stabilizzare** il modo servono
almeno ~41 dB di attenuazione (per portare +29 dB sotto i -12 dB tipicamente
richiesti a un modo gain-stabilizzato); il notch profondo scelto (-51 dB) porta il
lobo a **-18 dB**, cioe' 18 dB di bending gain margin.

### L'ablazione che dimostra la tesi

Valutando lo stesso guadagno di anello con `meas = 'true'`:

    |L(j*wBM)| con 'ins'  = ~+29.0 dB     (instabile senza filtro)
    |L(j*wBM)| con 'true' = ~-19.4 dB
    |L(j*wBM)| plant rigido = ~-19.4 dB   (identico!)

Con la retroazione sugli stati veri il modo di bending, pur essendo **eccitato**
dal TVC ed **esistente** nel plant, e' **completamente invisibile al loop**: `L(s)`
coincide con quella del corpo rigido. Il che prova che il colpevole non e' la
flessibilita' della struttura ne' l'eccitazione da parte dell'attuatore, ma il
**sensor coupling**: e' l'INS che chiude l'anello elastico. Da questa osservazione
seguono due famiglie di rimedi -- filtrare la misura (il notch, scelta di HM3) o
spostare/combinare i sensori in modo che `sigma_ins` si annulli (blending di
giroscopi in stazioni diverse, soluzione industriale).

> **Possibile domanda d'esame** -- Perche' non basta abbassare `Kd` per gain-stabilizzare
> il bending, visto che i +18.6 dB vengono quasi tutti dal termine derivativo?
> *Risposta:* perche' `Kd` e' anche cio' che da' il margine di fase al crossover
> rigido (~2.5 rad/s): ridurlo abbastanza da tagliare 29 dB a 18.9 rad/s (servirebbe
> un fattore ~30) distruggerebbe lo smorzamento del modo rigido e il PM. Serve un
> filtro **selettivo in frequenza** che attenui a `wBM` senza toccare la banda di
> controllo: e' esattamente la definizione di un notch. Il prezzo lo si paga
> comunque in fase (il notch + TVC + ritardo fanno crollare il PM da 30 a 14.6 deg,
> e infatti il PD va ri-tunato sul loop completo).

---

## Possibili domande d'esame

**D: Quali sono i sei stati e cosa rappresentano?**
R: `z` (drift laterale, m), `zdot` (velocita' laterale, m/s), `theta` (assetto di
beccheggio, rad), `thetadot` (velocita' angolare, rad/s), `eta` (coordinata modale
generalizzata del primo modo di flessione, m), `etadot` (sua derivata). I primi
quattro sono il corpo rigido, gli ultimi due un oscillatore del secondo ordine con
`omega_BM = 18.9 rad/s` e `zeta_BM = 0.005`.

**D: Perche' il modo di bending sta a 18.9 rad/s e non a 2 rad/s, e perche' e'
comunque un problema se la banda di controllo e' ~2.5 rad/s?**
R: 18.9 rad/s (3 Hz) e' la prima frequenza propria della struttura, circa un ordine
di grandezza sopra il crossover del controllo -- in teoria "fuori banda". Il
problema e' che con `zeta = 0.005` il picco di risonanza vale `1/(2*zeta) = 100`
(+40 dB) e che il controllore PD ha guadagno **crescente** con la frequenza
(termine derivativo): il roll-off naturale non c'e'. Risultato: nonostante la
separazione in frequenza, il guadagno di anello alla risonanza e' +29 dB invece
che ben sotto 0 dB.

**D: Differenza fra gain stabilization e phase stabilization del modo elastico, e
cosa hai scelto?**
R: *Gain stabilization* = attenuare il modo sotto 0 dB (tipicamente <= -12 dB) cosi'
che la fase alla risonanza sia irrilevante; *phase stabilization* = lasciare
passare il guadagno ma garantire che la fase alla risonanza mantenga il margine
(si usa quando il modo e' troppo vicino in frequenza al crossover per poterlo
attenuare). Qui il modo e' ben separato in frequenza, quindi ho gain-stabilizzato
con un notch profondo, arrivando a `|L(wBM)| = -18 dB`. Il prezzo e' la sensibilita'
alla conoscenza esatta di `omega_BM` (il notch tollera -10 % di detuning ma diventa
instabile a +5 %), che e' il tipico trade-off del notch profondo.

**D: Il modello dice che il bending non retroagisce sul corpo rigido e che il vento
non lo eccita. Sono ipotesi accettabili?**
R: Sono le ipotesi della traccia e vanno dichiarate. Non retroazionare l'elastico
sul rigido e' ragionevole per il **primo** modo di un lanciatore, dove la massa
generalizzata coinvolta e' una frazione piccola della massa totale e l'interesse e'
la stabilita' del loop, non il calcolo dei carichi. Trascurare l'eccitazione del
bending da parte della raffica e' piu' discutibile: significa che le simulazioni di
raffica sottostimano l'ampiezza modale reale e quindi i carichi strutturali. Per il
**progetto del controllore** non cambia nulla (`alpha_w` non entra in `L(s)`), per
un'analisi di carichi cambierebbe.

**D: Perche' `build_plant_full` espone 7 uscite quando il controllore ne usa 4?**
R: Perche' le prime 4 sono le **misure** (contaminate dal bending in modalita'
`'ins'`) e le ultime 3 sono i **valori veri** di `theta`, `z`, `zdot`, che servono
per i grafici e per le metriche di prestazione. La distinzione e' la ragione
d'essere del file: nel plant rigido le due famiglie coincidono, qui no, e mostrare
in figura la differenza fra `theta` e `theta_m` e' il modo piu' diretto di
visualizzare la contaminazione INS. Mantenere le stesse 7 uscite nei due plant
permette poi di passarli indifferentemente a `assemble_loop`.
