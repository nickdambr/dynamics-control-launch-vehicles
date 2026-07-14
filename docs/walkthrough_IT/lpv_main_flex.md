# HM3/LTV_FULL_ASCENT/main_flex.m

## Ruolo del file nel progetto

E' il gemello flessibile di `main_full_ascent.m` (ticket T008, Goal 1). Il Task 2 di HM3
stabilizza il primo modo di bending **in guadagno** (gain stabilisation) mettendo un notch
profondo centrato sulla frequenza modale a max-q, `omega_BM(72) = 18.9 rad/s`. Ma quella
frequenza **non e' una costante**: il primo modo flessionale di un lanciatore e'
essenzialmente `omega ~ sqrt(k/m)`, e mentre il propellente brucia la massa crolla, quindi
`omega` **sale** lungo l'ascesa. Sul dataset del corso `omega(t)` va da 16.6 a 30.1 rad/s
sull'orizzonte simulato (5-140 s): **quasi raddoppia**.

Un notch fisso, quindi, si **detuna**: dopo un po' non copre piu' la risonanza. Questo
script lo dimostra in due modi complementari e li confronta:

- **frequenza** (righe 21-33): congela il plant flessibile completo ad ogni istante e
  legge il guadagno d'anello **alla frequenza di bending istantanea**, `|L(j*omega(t))|`,
  con il notch fisso e con il notch che insegue `omega(t)`; registra anche la stabilita'
  dell'anello chiuso congelato.
- **tempo** (righe 35-41): integra il vero sistema LTV flessibile a 13 stati
  (`ode_lpv_flex`) con il vento del generatore nell'anello, una volta col notch fisso e
  una col notch variabile.

La risposta alla domanda del brief -- *il notch va anch'esso schedulato?* -- e' **si', e il
codice lo fa**: il centro del notch e' un handle `fwn(t)` passato al RHS (riga 40:
`make_flex(S, S.fomega)` per il caso variabile, riga 41: `make_flex(S, @(t) w72)` per il
fisso), e dentro `ode_lpv_flex` (riga 30) `wn = M.fwn(t)` entra nei coefficienti del filtro
ad ogni passo. **Solo il centro** e' schedulato: gli smorzamenti `zN = 0.002` e `zD = 0.7`
restano costanti, e cosi' pure il TVC e i guadagni PD.

Dipende da `init_simulink_lpv`, `ode_lpv_flex`, e dai moduli HM3 `build_plant_full`,
`build_tvc`, `build_notch_filter`, `assemble_loop`, `load_hw3_params`. E' la sorgente di
verita' di `hm3_full_ascent_flex.slx` (autorato da `build_hm3_full_ascent_flex.m`, con il
notch realizzato da blocchi elementari perche' nessun blocco di libreria ha i coefficienti
tempo-varianti).

**Buona notizia rispetto agli altri due script della cartella:** `main_flex` usa i guadagni
PD **frozen** (`sched = false`, riga 105), quindi **non e' toccato** dal problema della gain
schedule descritto in `lpv_main_full_ascent.md`. I suoi risultati si riproducono.

---

## Intestazione e preambolo (righe 1-19)

- Righe 1-11: header. Dichiara `omega(t)` in 16.5 -> 31.8 rad/s. Attenzione: **31.8 rad/s
  e' il valore a t = 150 s**, cioe' fuori dall'orizzonte simulato. Verificato sui dati: sul
  dataset completo [0, 150] s `omega` va da 16.48 a 31.77 rad/s, ma sull'orizzonte
  effettivo dello script -- **[5, 140]**, perche' `init_simulink_lpv` ha `t0 = 5` di
  default e la griglia di sweep e' `5:5:140` -- va da
  **16.65 a 30.12 rad/s**. Il numero che lo script stampa a runtime (riga 45-46) e'
  quest'ultimo; l'header e il README riportano l'altro.
- Riga 14: si spegne il warning di margine instabile (l'anello e' condizionalmente stabile).
- Riga 18: `S = init_simulink_lpv()` -- di nuovo il setup completo (dati LPV, 28 design PD,
  simulazione del generatore di vento).
- Riga 19: `w72 = S.notch.wn72`, cioe' `omega` interpolata a t = 72 s. Verificato: **18.900
  rad/s**, esattamente il valore di Tabella 1 della traccia (il dataset e la tabella
  coincidono al punto di design).

---

## `%% Frozen-time detuning sweep` (righe 21-33)

```matlab
tm = (t0:5:Tend).';
for i = 1:numel(tm)
    p = load_hw3_params('t_ref', tm(i));
    G = build_plant_full(p, 'ins');  Wt = build_tvc(p, 3);
    [Lf, Tf] = assemble_loop(G, S.K0, ...
                Wt*build_notch_filter(w72,   0.002, 0.7, +1));
    [Lv, Tv] = assemble_loop(G, S.K0, ...
                Wt*build_notch_filter(p.wBM, 0.002, 0.7, +1));
    Lfix(i) = 20*log10(abs(squeeze(freqresp(Lf, p.wBM))));
    Lvar(i) = 20*log10(abs(squeeze(freqresp(Lv, p.wBM))));
    stbF(i) = isstable(Tf);  stbV(i) = isstable(Tv);
end
t_unstable = tm(find(~stbF, 1));
```

- Riga 26: `build_plant_full(p, 'ins')` -- plant a **6 stati** `[z zdot theta thetadot eta
  etadot]` con la **contaminazione INS**: la piattaforma inerziale non misura theta pura ma
  `theta_m = theta + sigma_ins*eta` e `z_m = z - phi_ins*eta` (Eq. 2 della traccia). E'
  proprio questa contaminazione a chiudere l'anello sul modo flessibile e a renderlo
  destabilizzabile: senza di essa il bending non tornerebbe nel controllore.
- Riga 26: `build_tvc(p, 3)` -- servo del 2o ordine (`wTVC = 70 rad/s`, `zTVC = 0.7`) in
  serie con l'approssimante di **Pade di ordine 3** del ritardo puro `tau = 0.020 s`. Nota:
  `wTVC`, `zTVC`, `tau` sono **letterali di Tabella 1**, non presenti in
  `GreensiteLPV_DATA.mat` -- quindi `Wt` e' **identico ad ogni iterazione**: l'attuatore
  non e' tempo-variante e non ha senso schedularlo. (Ricostruirlo dentro il loop e' solo
  lavoro sprecato, non un errore.)
- Righe 27-28: i due anelli. La catena serie e' `Wact = TVC * notch` (in LTI l'ordine non
  conta). Il notch e' `build_notch_filter(wx, 0.002, 0.7, +1)`, cioe' **esattamente** il
  "deep notch" che HM3 sceglie in `main_task2.m` (righe 74-78: `zN = 0.002`, `zD = 0.7`,
  `sgn = +1`, min-phase). L'unica differenza fra `Lf` e `Lv` e' il **centro**: `w72` fisso
  contro `p.wBM = omega(t)` variabile.
- Righe 29-30: la metrica e' `|L(j*omega(t))|` in dB, cioe' il guadagno d'anello **alla
  frequenza di risonanza istantanea**. E' la grandezza giusta da guardare per una
  stabilizzazione in guadagno: il criterio di progetto e' "tieni `|L|` ben sotto 0 dB
  dove il modo risuona".
- Riga 31: `isstable(Tf)` sull'anello **chiuso** congelato -- il verdetto vero, che non
  dipende da come leggi i margini.
- Riga 33: `t_unstable` = primo istante in cui l'anello col notch fisso e' instabile.

**Risultati verificati eseguendo il codice:**

| t [s] | 5 | 40 | 55 | 65 | 70 | 75 | 90 | 140 |
|---|---|---|---|---|---|---|---|---|
| omega(t) [rad/s] | 16.65 | 17.83 | 18.33 | 18.67 | 18.83 | 19.39 | 21.87 | 30.12 |
| `\|L\|` notch fisso [dB] | +7.8 | +4.8 | +0.5 | -6.4 | **-15.8** | +0.3 | +14.9 | +23.2 |
| `\|L\|` notch variabile [dB] | -28.1 | -24.6 | -23.2 | -22.4 | -22.0 | -21.9 | -22.2 | -22.7 |
| anello chiuso, notch fisso | stabile | stabile | stabile | stabile | stabile | **INSTABILE** | INSTABILE | INSTABILE |

- `t_unstable = 75 s`. Il notch fisso perde l'anello **3 secondi dopo il punto di design**.
- Il notch variabile tiene `|L(j*omega(t))|` fra **-28.1 e -21.9 dB** su tutta l'ascesa,
  sempre stabile. (Il README dichiara "-12...-18 dB": **non e' quello che esce oggi**, il
  margine reale e' migliore di quanto il README prometta.)

**Quanto e' stretto quel notch?** Il valore al centro esatto e' `|H(j*wn)| = zN/zD =
0.002/0.7 = 2.9e-3`, cioe' **-51 dB**: una tacca profondissima e sottilissima. Il prezzo si
legge nella tabella: a t = 70 s il notch e' disallineato dello **0.4 %** e attenua ancora
-15.8 dB; a t = 65 s il disallineamento e' dell'**1.2 %** e restano solo -6.4 dB; a
t = 75 s siamo al **2.6 %** e l'attenuazione e' **sparita** (+0.3 dB). *Circa il 2-3 % di
errore in frequenza consuma tutto il notch.* E' proprio la ragione per cui HM3, nel Task 2,
propone anche una **tripletta di notch** a `{0.9, 1.0, 1.1} * omega_BM` (`main_task2.m`,
righe 81-83): una tacca larga e' la difesa contro l'incertezza su `omega_BM` **quando non
la si puo' schedulare**. Qui, avendo `omega(t)` dal dataset, si sceglie invece di
inseguirla.

> **Possibile domanda d'esame** -- nella tabella, a t = 5-55 s il notch fisso ha
> `|L(j*omega)| > 0 dB` eppure l'anello chiuso e' **stabile**. Come e' possibile?
> *Risposta:* perche' `|L| > 0 dB` alla frequenza di risonanza **non e' un criterio di
> instabilita'**. La stabilizzazione di un modo flessibile puo' avvenire in due modi: **in
> guadagno** (si tiene `|L|` sotto 0 dB attorno alla risonanza, ed e' quello che fa il
> notch) oppure **in fase** (si lascia passare il guadagno ma si garantisce che il lobo di
> bending sul Nichols passi dalla parte giusta del punto critico). Nell'ascesa iniziale il
> lobo, pur sopra 0 dB, ha una fase che non produce accerchiamento: il modo e' *phase
> stabilised per caso*. Il criterio vero e' quello della riga 31, `isstable` sull'anello
> chiuso. La `yline` alla riga 65, etichettata "0 dB (resonance uncovered)", **semplifica
> troppo** ed e' l'unico punto della figura da difendere con cautela.

---

## `%% Time-domain flexible response` (righe 35-41)

```matlab
tt = (t0:0.02:Tend).';
nt = size(S.tvc.At, 1);
x0 = zeros(6 + 2 + nt, 1);
odeo = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
[~, xv] = ode45(@(t,x) ode_lpv_flex(t,x,make_flex(S,S.fomega)), tt, x0, odeo);
[~, xf] = ode45(@(t,x) ode_lpv_flex(t,x,make_flex(S,@(t) w72)),  tt, x0, odeo);
```

- Riga 37: `nt = 5` (verificato: servo di 2o ordine + Pade 3 = 5 stati), quindi lo stato
  totale e' **6 + 2 + 5 = 13**: 6 di plant flessibile, **2 di notch**, 5 di TVC.
- Righe 40-41: **unica differenza fra le due corse: l'handle `fwn`**. `S.fomega` e' il
  `griddedInterpolant` di `omega(t)`; `@(t) w72` e' la costante. Tutto il resto (plant,
  vento, guadagni, TVC, `zN`, `zD`, `zBM`) e' identico. Confronto controllato pulito.

**Cosa fa `ode_lpv_flex` (il file chiamato, per capire cosa si sta integrando):**

- Righe 24-27: **misura INS contaminata** e legge di controllo.

      theta_m = theta + sigma*eta      z_m    = z - phi*eta
      u_pd    = -(Kp*theta_m + Kd*thetadot_m + Kp_z*z_m + Kd_z*zdot_m)

  Il PD **non vede theta**, vede `theta_m`: e' l'eta' del bending che rientra nell'anello.
- Righe 30-33: **il notch tempo-variante**, realizzato in forma canonica di controllo:

      wn   = fwn(t)
      xn1' = xn2
      xn2' = -wn^2*xn1 - 2*zD*wn*xn2 + u_pd
      v    = 2*(zN - zD)*wn*xn2 + u_pd

  Derivazione (vale la pena saperla): per `H(s) = (s^2 + b1 s + b0)/(s^2 + a1 s + a0)`, la
  forma canonica di controllo da' `y = (b0 - a0)*x1 + (b1 - a1)*x2 + u`. Nel notch
  `b0 = a0 = wn^2` (numeratore e denominatore hanno **lo stesso termine noto**: e' cio' che
  rende la tacca simmetrica attorno a `wn`), quindi il termine in `x1` **si cancella**;
  restano `b1 - a1 = 2*zN*wn - 2*zD*wn = 2*(zN - zD)*wn` e il passante diretto `+u`. Da qui
  la riga 31 esattamente com'e' scritta. Con `zN << zD` il coefficiente e' negativo e
  grande, ed e' quello che scava la tacca.
- Righe 34-35: `v` (uscita notch) entra nel TVC (`At, Bt, Ct, Dt`, LTI) che produce la
  **deflessione fisica** `delta`. Ordine della catena: PD -> notch -> TVC -> plant.
- Righe 38-40: le equazioni del plant flessibile LTV,

      zddot   = a1*zdot + (a1*V + a4)*theta + a3*delta - a1*V*alpha_w
      thddot  = (A6/V)*zdot + A6*theta + K1*delta - A6*alpha_w
      etaddot = -omega^2*eta - 2*zBM*omega*etadot + aqk*delta

  con `aqk = -phi_tvc*Tc` (verificato: `aqk(72) = -65.56`, coerente con
  `-4.31e-5 * 1.52e6` di Tabella 1): e' il **momento che il TVC inietta nel modo
  flessibile** deflettendo l'ugello. Il vento **non** forza direttamente il bending
  (nessun `alpha_w` nella terza equazione), coerentemente con `build_plant_full`
  (`Bw(6) = 0`).

**Limite tecnico da dichiarare (importante).** Realizzare un filtro tempo-variante
congelando la sua forma canonica **non e' esatto**: per un sistema LTV non esiste la
"funzione di trasferimento", e realizzazioni diverse della stessa `H(s, wn)` congelata
danno **sistemi LTV diversi** (i termini `d(wn)/dt` che comparirebbero derivando la
trasformazione di stato sono qui semplicemente **omessi**). L'approssimazione e'
difendibile perche' `omega` varia lentamente: circa 0.1 rad/s al secondo su un `omega ~ 20
rad/s`, cioe' **~0.5 % per secondo**, mentre la dinamica del filtro ha costanti di tempo
dell'ordine di `1/(zD*wn) ~ 0.07 s`. Ma va detto, non nascosto.

---

## `%% Summary` (righe 43-51)

Stampa il range di `omega` sulla griglia di sweep, l'istante di prima instabilita', il
range di `|L|` col notch variabile, e i picchi di `theta`. Verificato eseguendo:

    omega(t) sweeps 16.6 -> 30.1 rad/s ; HM3 notch fixed at 18.9 rad/s
    fixed   notch: loop goes UNSTABLE at t = 75 s
    varying notch: stable all ascent, |L| in [-28.1, -21.9] dB
    peak |theta|: varying = 1.014 deg (bounded) ; fixed = 1.04e+20 deg (diverges)

Sul picco `|eta|`: **variabile 3.5e-3**, **fisso 1.2e20**. Il `1e20` e' ovviamente privo di
senso fisico: il modello e' **lineare, senza saturazione dell'attuatore e senza limiti
strutturali** -- un lanciatore vero si romperebbe molto prima. Va letto come "diverge", non
come un numero.

---

## `%% Figures` (righe 53-81) e `%% Export` (righe 83-90)

- f1 (righe 55-69): lo sweep di detuning, con un `patch` rosso che marca la regione in cui
  il notch fisso rende l'anello instabile (da `t_unstable` in poi). E' la figura che
  racconta la storia.
- f2 (righe 72-81): dominio del tempo. Riquadro sinistro: `theta` con **`ylim([-2 2])`**
  (riga 76) -- la curva del notch fisso **esce dal grafico**, e il titolo lo dichiara
  onestamente ("theta off-scale = divergence"). Riquadro destro: `|eta|` in scala
  **semilogaritmica** (riga 79), l'unico modo per mostrare sullo stesso asse 1e-3 e 1e20;
  il `+1e-12` evita `log(0)`.
- Righe 86-89: export PNG con tema chiaro.

---

## `make_flex` (righe 93-106)

Costruisce la struct per `ode_lpv_flex`. Punti da notare:

- Riga 102: `'fwn', fwn` -- **e' l'unico parametro che distingue le due corse**.
- Riga 105: `'sched', false` -- **i guadagni PD restano quelli frozen di max-q**
  (`S.K0`). I campi `fKp`/`fKd` sono passati ma **mai usati** (codice morto innocuo). Come
  scrive il README, il Goal 1 tiene i guadagni fermi e muove solo il notch; combinare le due
  cose e' un follow-up dichiarato aperto.
- Righe 102-103: `zN`, `zD`, `zBM` sono **costanti**. `zBM = 0.005` viene da Tabella 1 (il
  dataset non fornisce uno smorzamento modale variabile), quindi non c'e' nulla da
  schedulare li'. Il TVC (`At..Dt`) e' l'LTI costruito una volta a t = 72 s.

**Divergenza da HM3 da dichiarare all'orale.** `S.K0` sono i guadagni progettati sul plant
**rigido con attuatore ideale** (Task 1: `design_controller(build_plant_rigid(p0), [])`,
riga 72 di `init_simulink_lpv`). Ma HM3, nel Task 2, dopo aver messo il notch **ri-tara il
PD sull'anello completo** (`main_task2.m`, righe 131-161), perche' -- parole del sorgente --
"the actuator + transport delay + notch lag collapses the rigid phase margin". `main_flex`
**salta quella ritaratura**: usa i guadagni Task-1 con la catena Task-2. Ho verificato la
conseguenza congelando l'anello flessibile a t = 72 s con `S.K0`: aero GM = -6.08 dB (ok),
`|L(omega_BM)| = -21.9 dB` (ok, il notch fa il suo lavoro), anello chiuso **stabile**, ma
**rigid PM = 14.6 deg** contro i 30 deg di target. Cioe' `main_flex` gira nella
configurazione che HM3 chiama "BEFORE re-tuning". Questo **non invalida** il confronto
fisso-vs-variabile (entrambe le corse usano lo stesso PD, quindi la differenza e'
attribuibile al solo notch), ma i **margini assoluti** di questo studio non sono quelli del
design finale di HM3.

---

## La fallacia del frozen-time, vista qui

Questo script e' il posto giusto per capire perche' "ogni LTI congelato e' stabile" **non
e'** un teorema di stabilita' per il sistema tempo-variante -- e vale anche il viceversa.

Il fatto: l'anello col notch fisso e' **congelato-instabile da t = 75 s**. Se la stabilita'
congelata implicasse quella LTV punto per punto, ci si aspetterebbe divergenza immediata.
Invece, verificato sulla simulazione:

| t [s] | 85 | 100 | 140 |
|---|---|---|---|
| `\|eta\|` notch fisso | 3.9e-4 | 6.6e-2 | 1.9e18 |
| `\|eta\|` notch variabile | 3.9e-4 | 5.2e-5 | 1.3e-5 |

A **t = 85 s**, dieci secondi *dopo* che l'anello congelato e' diventato instabile, la
coordinata di bending del caso fisso e' **indistinguibile** da quella del caso stabile
(3.9e-4 contro 3.9e-4). La divergenza diventa visibile solo intorno a **90-100 s** e poi
esplode. Motivo: l'instabilita' congelata dice che esiste un **autovalore a parte reale
positiva** in quel punto, ma la crescita e' **esponenziale con una costante di tempo
finita**; su un orizzonte finito il modo instabile ha bisogno di tempo per emergere dal
rumore. Un'analisi puramente congelata avrebbe dichiarato il veicolo perso a 75 s; il
veicolo vero (LTV) sopravvive ancora un po'. E, simmetricamente, esistono sistemi in cui
**tutti** i congelati sono stabili e l'LTV diverge lo stesso (il controesempio classico:
matrici i cui autovalori stanno sempre a sinistra ma i cui autovettori ruotano
rapidamente).

La regola operativa, ed e' quella che questi tre script applicano: **il design a punti
congelati e' un'euristica di sintesi, non una prova**. La prova e' l'integrazione del vero
LTV (o, se si vuole una garanzia, una funzione di Lyapunov parameter-dependent con un
vincolo sul tasso di variazione dei parametri). Qui la separazione di scale e' generosa --
i coefficienti variano su decine di secondi, la dinamica di bending ha periodo `2*pi/20 ~
0.3 s` -- quindi il congelamento e' ben giustificato; ma la giustificazione va **esibita**,
non presupposta.

---

## Possibili domande d'esame

**D: Perche' la frequenza del primo modo di bending cresce lungo l'ascesa?**
R: Perche' e' essenzialmente `omega ~ sqrt(k/m)` (rigidezza modale su massa modale). La
rigidezza strutturale del corpo cambia poco, ma la **massa crolla** con il consumo di
propellente (il lanciatore brucia gran parte della sua massa nel primo stadio), quindi
`omega` **sale**. Sul dataset del corso passa da 16.5 a ~30 rad/s fra il decollo e t = 140 s:
quasi il doppio. Un notch progettato a max-q e lasciato li' e' fuori bersaglio dopo pochi
secondi.

**D: Il notch va schedulato o basta allargarlo?**
R: Entrambe le strade sono legittime e HM3 le esplora tutte e due. Allargarlo (o metterne
tre, a `0.9/1.0/1.1 * omega_BM`, come in `main_task2.m`) e' la difesa quando `omega_BM` e'
**incerta ma non nota**: si copre un intervallo, pagando piu' ritardo di fase alla
frequenza di crossover rigido e quindi margine di fase. Schedularlo e' la scelta quando
`omega(t)` e' **nota o stimabile** (qui viene dal dataset): si mantiene la tacca stretta e
profondissima -- -51 dB al centro, con `zN/zD = 0.002/0.7` -- senza pagarne la larghezza.
Il codice dimostra il costo del non farlo: **2-3 % di errore in frequenza consuma tutta
l'attenuazione**.

**D: Cosa succede esattamente quando il notch si detuna? Perche' l'anello si perde?**
R: Il modo flessibile rientra nel controllore attraverso la **contaminazione INS**
(`theta_m = theta + sigma*eta`, `z_m = z - phi*eta`): il sensore misura anche la vibrazione
della struttura. Se il guadagno d'anello a quella frequenza non e' attenuato, il PD reagisce
alla vibrazione comandando il TVC, il TVC **eccita** il modo (termine `aqk*delta`
nell'equazione di `eta`, con `aqk = -phi_tvc*Tc`), e si chiude un anello positivo
strutturale: e' il classico **accoppiamento aeroservoelastico** (pogo/tail-wags-dog). Il
notch spezza l'anello attenuando `|L|` proprio li'.

**D: Con il notch fisso l'anello congelato e' instabile da t = 75 s, ma la simulazione
diverge solo dopo ~90-100 s. Contraddizione?**
R: No, e' proprio la **fallacia del frozen-time** vista dal verso benigno. Instabilita'
congelata significa che a quell'istante esiste un autovalore a parte reale positiva:
la crescita e' esponenziale, ma con una costante di tempo finita, e parte da un'ampiezza
piccolissima (a t = 85 s `|eta|` vale ancora 3.9e-4, come nel caso stabile). Su un orizzonte
finito il modo ha bisogno di secondi per emergere. Il corollario e' che **la stabilita' dei
congelati non e' ne' necessaria ne' sufficiente** per la stabilita' dell'LTV: e' un'euristica
valida sotto ipotesi di variazione lenta, che qui e' soddisfatta e che comunque va
verificata integrando il sistema vero -- che e' esattamente cio' che fa questo script.

**D: Il tuo notch tempo-variante e' esatto?**
R: No, ed e' bene dirlo. Realizzo il notch in forma canonica di controllo e poi rendo
tempo-varianti i coefficienti (`wn = wn(t)` nelle righe 30-33 di `ode_lpv_flex`). Per un
sistema LTV il concetto di funzione di trasferimento non esiste, e realizzazioni diverse
dello stesso filtro congelato producono **dinamiche LTV diverse**: i termini in `d(wn)/dt`
sono omessi. L'approssimazione e' difendibile perche' `omega` varia dello ~0.5 % al secondo
mentre il filtro ha costanti di tempo di ~0.07 s -- separazione di scale di oltre due ordini
di grandezza -- ma e' un'approssimazione, non un'identita'.

**D: Perche' `main_flex` non schedula anche i guadagni PD?**
R: Per isolare l'effetto. Se schedulassi contemporaneamente notch e guadagni non saprei a
quale dei due attribuire la differenza. Il Goal 1 (questo file) tiene i guadagni frozen e
muove solo il notch; il Goal 2 (`main_q_scheduling.m`) muove i guadagni sul plant rigido. Il
README dichiara apertamente che la combinazione dei due e' un follow-up ancora aperto. Va
pero' segnalato che i guadagni usati qui sono quelli **Task-1 (attuatore ideale)** e non
quelli ritarati sull'anello completo come fa HM3 nel Task 2: con la catena TVC + ritardo +
notch, il margine di fase rigido scende a **14.6 deg** contro i 30 deg di target.
