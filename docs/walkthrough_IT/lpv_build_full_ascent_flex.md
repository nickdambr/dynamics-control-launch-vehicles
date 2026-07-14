# HM3/LTV_FULL_ASCENT/build_hm3_full_ascent_flex.m

## Ruolo del file nel progetto

Come il suo gemello rigido, questo file **e'** il modello Simulink: un authoring
script che con `new_system` / `add_block` / `add_line` / `set_param` /
`save_system` costruisce da zero `HM3/LTV_FULL_ASCENT/hm3_full_ascent_flex.slx`.
Non produce array di `ss`, non tocca alcun modello preesistente: cancella e
ricrea. Il `.slx` e' un artefatto derivato.

Rispetto a `build_hm3_full_ascent.m`, aggiunge tutta l'**architettura
flessibile** che in HM3 vive solo a max-q congelato (`build_plant_full`):

- primo modo di **bending** con frequenza `omega(t)` tempo-variante;
- accoppiamento **INS**: i sensori misurano `theta_m = theta + sigma*eta`,
  `z_m = z - phi*eta` (e le rispettive derivate), cioe' la piattaforma inerziale
  legge anche la deformazione strutturale;
- **notch** realizzato da blocchi elementari, il cui centro **insegue `omega(t)`**;
- **TVC** (servoattuatore + ritardo di trasporto approssimato con Pade) come
  funzione di trasferimento LTI.

E toglie una cosa: **il gain scheduling**. I guadagni PD sono congelati a max-q
(header, riga 5: "PD gains FROZEN at max-q; the showcase is notch tracking, not
scheduling"). Il messaggio del track e' un altro: la frequenza di bending spazia
su `omega(t) in [16.48, 30.12] rad/s` sull'orizzonte effettivamente simulato
(0-140 s; il dataset grezzo arriva a 31.77 rad/s, ma a t = 150 s, fuori
orizzonte -- `init_simulink_lpv` riga 51 tronca la griglia a `Tstop`), quindi un
notch fisso centrato su `omega(72)` si **detuna** e l'anello perde la
gain-stabilization del modo strutturale.

La controparte MATLAB, dichiarata sorgente di verita', e' `ode_lpv_flex.m`
(13 stati). Chi lo chiama: `run_flex_simulink.m`. Le variabili base necessarie
(`lpv_omega`, `lpv_omega2`, `lpv_2zBMw`, `lpv_aqk`, `lpv_sig`, `lpv_phi`,
`notch_zN`, `notch_zD`, `tvc_num`, `tvc_den`, oltre a tutte quelle del modello
rigido) sono pushate da `init_simulink_lpv.m`.

**Nota fondamentale, dichiarata subito**: a differenza del track rigido, qui il
`.slx` e la RHS `ode45` **non interpolano le stesse tabelle**. Vedi la sezione
"La disciplina dei coefficienti effettivi, rotta".

---

## Firma, header e setup (righe 1-41)

```matlab
lk = @(name, tbl, bp, pos) A(name, LK, pos, ...
        'NumberOfTableDimensions','1', ...
        'BreakpointsForDimension1', bp, 'Table', tbl, ...
        'ExtrapMethod','Clip');
```

- Righe 18-20: unica opzione `o.open`, come nel rigido.
- Righe 28-29: carica `strong_wind.slx` e installa l'`onCleanup` che lo chiude
  con flag `0` (discard). Il file del professore non viene mai salvato.
- Righe 31-35: `new_system` da zero; solver `ode45` variable-step,
  `StopTime = 'Tstop'` (nome di variabile base), `MaxStep = 0.02` per risolvere
  le innovazioni a 0.1 s del generatore di rumore. Come nel rigido, **`RelTol` e
  `AbsTol` non vengono mai impostati**: restano ai default Simulink (1e-3),
  mentre `run_flex_simulink.m` (riga 26) fa la replica con
  `odeset('RelTol',1e-9,'AbsTol',1e-11)`.
- Righe 40-41: qui c'e' una differenza sostanziale rispetto al rigido. L'helper
  `lk` mette **`ExtrapMethod = 'Clip'` su tutte** le lookup, non solo su quelle
  dei guadagni. `Clip` (tieni il valore d'estremo) e' cio' che corrisponde
  all'estrapolazione `'nearest'` usata dalle `griddedInterpolant` della baseline.
  Nel modello rigido questo non era stato fatto per `c1..c7`. Buona correzione,
  anche se su `[0, 140]` con `lpv_t = 0:1:140` non si estrapola comunque mai.

---

## Clock e le 14 lookup (righe 43-51)

```matlab
names = {'c1','c2','c3','c4','c5','c6','c7','invV', ...
         'omega2','w2zBM','aqk','sig','phi','omega'};
tbls  = {'lpv_c1', ... ,'lpv_invV', ...
         'lpv_omega2','lpv_2zBMw','lpv_aqk','lpv_sig', ...
         'lpv_phi','lpv_omega'};
for k = 1:numel(names)
    lk(names{k}, tbls{k}, 'lpv_t', [...]);
    W('Clock/1', [names{k} '/1']);
end
```

- Riga 44: un unico `Clock`: come nel rigido, il parametro di scheduling e' il
  tempo di volo e nient'altro.
- Righe 45-51: **quattordici** lookup 1-D, tutte sugli stessi breakpoint `lpv_t`
  (griglia verificata sul dataset: `0:1:140`, **141 punti, passo 1 s**),
  interpolazione lineare (default `InterpMethod`). Otto sono quelle del plant
  rigido (`c1..c7`, `invV`); le sei nuove sono:

  | blocco   | tabella       | contenuto (da `init_simulink_lpv`)      |
  |----------|---------------|------------------------------------------|
  | `omega2` | `lpv_omega2`  | `omega.^2` -- **gia' elevato al quadrato** |
  | `w2zBM`  | `lpv_2zBMw`   | `2*p0.zBM*omega` -- **gia' moltiplicato**  |
  | `aqk`    | `lpv_aqk`     | forzante di bending del TVC (dal dataset) |
  | `sig`    | `lpv_sig`     | `sigma_ins`, leakage INS su theta         |
  | `phi`    | `lpv_phi`     | `phi_ins`, leakage INS su z               |
  | `omega`  | `lpv_omega`   | `omega(t)` grezzo (serve al notch)        |

- `lpv_2zBMw = 2*p0.zBM*omega`: nota che `zBM`, lo smorzamento del modo, e'
  preso **congelato** dai parametri a `t_ref = 72 s` (`p0 = load_hw3_params()`).
  Solo `omega` e' tempo-variante. E' un'approssimazione mai dichiarata nel codice:
  la frequenza di bending varia dell'83% (16.48 -> 30.12 rad/s sull'orizzonte
  0-140 s) mentre lo smorzamento e' tenuto costante.
- `lpv_aqk` viene dal campo `aqk` del dataset; il commento in `init_simulink_lpv`
  (riga 57) lo descrive come `-phi_tvc*Tc`, cioe' la spinta di controllo
  proiettata sulla forma modale nel punto di applicazione del TVC.

---

## Il generatore di vento e alpha_w (righe 53-58)

Identico al modello rigido: il sottosistema del professore viene **copiato**
(riga 54), il Clock lo alimenta, un `Sum '++'` somma profilo medio e turbolenza,
un `Product` moltiplica per `invV`:

    alpha_w(t) = ( v_wp(t) + turbolenza(t) ) * (1/V(t))

Vento e veicolo condividono lo stesso orologio: e' la differenza chiave rispetto
a HM3, che ritagliava 12 s di vento attorno a max-q.

---

## Plant z/theta (righe 60-83)

Topologia **identica** al modello rigido: `P1..P4 -> zdd('+++-') -> int_zd ->
int_z` e `P5..P8 -> thdd('+++-') -> int_thd -> int_th`, con `c6` che alimenta sia
`P6` (`c6*theta`) sia `P8` (`c6*alpha_w`, con segno meno nel `Sum`).

    zddot     = c1*zdot + c2*theta + c3*delta - c4*alpha_w
    thetaddot = c5*zdot + c6*theta + c7*delta - c6*alpha_w

- Righe 81-83: la retroazione degli stati arriva dalle uscite degli integratori,
  quindi nessun loop algebrico.
- Il commento di riga 60 ("identical to the rigid model") e' vero, e l'header
  (riga 11) aggiunge la cosa importante: **"bending does not feed it"**. Il modo
  flessibile `eta` **non** entra nell'equazione del corpo rigido: entra solo
  attraverso i **sensori** (accoppiamento INS) e viene eccitato dal TVC. E'
  l'ipotesi standard di modo strutturale disaccoppiato dalla dinamica di corpo
  rigido, e chiude l'anello attraverso il controllore, non attraverso la fisica.

---

## Bending (righe 85-96)

```matlab
A('etadd','built-in/Sum',[...],'Inputs','--+', ...);
A('int_etad','built-in/Integrator',[...],'InitialCondition','0');
A('int_eta', 'built-in/Integrator',[...],'InitialCondition','0');
W('omega2/1','Pe1/1'); W('int_eta/1','Pe1/2');
W('w2zBM/1','Pe2/1');  W('int_etad/1','Pe2/2');
```

- Righe 86-88: tre prodotti: `Pe1 = omega^2 * eta`, `Pe2 = (2*zBM*omega) *
  etadot`, `Pe3 = aqk * delta`.
- Riga 89: `Sum` con `Inputs = '--+'`, quindi

      etaddot = -omega(t)^2 * eta - 2*zBM*omega(t) * etadot + aqk(t) * delta

  E' l'oscillatore del secondo ordine del primo modo di bending, forzato dalla
  deflessione dell'ugello. E' esattamente la riga 40 di `ode_lpv_flex.m`.
- Righe 90-91: catena `int_etad -> int_eta`, condizioni iniziali nulle, coerenti
  con `zeros(6+2+nt,1)` della replica (`run_flex_simulink.m`, riga 41).
- Riga 94: `aqk` entra su `Pe3/1`; l'altro fattore, `delta`, arriva molto piu'
  tardi (riga 147) dall'uscita del TVC. E' il **solo** canale con cui il
  controllo eccita la struttura: se `aqk = 0` il modo resta a zero per sempre.

> **Possibile domanda d'esame** -- perche' il bending e' forzato da `delta` e
> non dal vento?
> *Risposta:* Nel modello del corso il modo flessibile e' eccitato dalla spinta
> laterale dell'ugello applicata in coda (`aqk = -phi_tvc*Tc`, cioe' la forma
> modale valutata nel punto di applicazione della spinta). Il vento entra come
> `alpha_w` nelle equazioni aerodinamiche di corpo rigido. E' una scelta di
> modellazione: si trascura l'eccitazione aeroelastica distribuita e si tiene il
> solo percorso che chiude l'anello di controllo, che e' quello che puo' rendere
> instabile il sistema (delta -> eta -> INS -> controllore -> delta).

---

## Accoppiamento INS (righe 98-114)

```matlab
A('theta_m','built-in/Sum',[...],'Inputs','++', ...);  % theta + sig*eta
A('z_m',    'built-in/Sum',[...],'Inputs','+-', ...);  % z - phi*eta
```

- Righe 99-106: quattro prodotti, `sig*eta`, `sig*etadot`, `phi*eta`,
  `phi*etadot`. Sia `sig` sia `phi` sono lookup su `t`: il leakage varia lungo
  l'ascesa.
- Righe 107-114: i quattro sommatori delle **misure**:

      theta_m    = theta    + sigma(t)*eta
      thetadot_m = thetadot + sigma(t)*etadot
      z_m        = z        - phi(t)*eta
      zdot_m     = zdot     - phi(t)*etadot

  I segni (`'++'` per l'assetto, `'+-'` per la deriva) riproducono le righe 24-25
  di `ode_lpv_flex.m`. Fisicamente: il giroscopio, montato in un punto della
  struttura, misura la rotazione locale, che e' quella di corpo rigido **piu'**
  la pendenza locale della deformata; l'accelerometro/posizione vede lo
  spostamento di corpo rigido **meno** lo spostamento modale nel suo punto.
  Segni opposti perche' `sigma` e `phi` sono, rispettivamente, la pendenza e lo
  spostamento della forma modale nel punto di misura.

Questo e' il meccanismo che rende il problema pericoloso: il controllore **non
puo' distinguere** l'assetto vero dalla vibrazione. Se il guadagno d'anello alla
frequenza di bending non e' attenuato, il PD reagisce a `sigma*eta`, comanda
`delta`, che tramite `aqk` rieccita `eta`: e' il flutter di controllo. Da qui il
notch.

---

## Controllore a guadagni congelati (righe 116-123)

```matlab
A('Gth', 'built-in/Gain',[...],'Gain','Kp_th0');
A('Gthd','built-in/Gain',[...],'Gain','Kd_th0');
A('u_pd','built-in/Sum',[...],'Inputs','----', ...);
```

- Righe 117-120: **quattro `Gain` costanti**. Non ci sono `Switch`, non ci sono
  lookup sui guadagni, non c'e' `sched`. Il PD e' il design HM3 Task-1 a max-q.
- Righe 121-123: `u_pd` con `Inputs = '----'`:

      u_pd = -( Kp*theta_m + Kd*thetadot_m + Kp_z*z_m + Kd_z*zdot_m )

  Corrisponde alla riga 27 di `ode_lpv_flex.m`. Attenzione: il controllore chiude
  sulle **misure** `theta_m, z_m, ...`, non sugli stati veri. E' proprio questo
  che porta `eta` dentro l'anello.

**Asimmetria fra i due track, da dichiarare**: `ode_lpv_flex.m` (riga 26)
conserva il ramo `if M.sched, Kp = M.fKp(t); ...` e riceve `fKp`/`fKd` da
`run_flex_simulink.m` (righe 39), che poi imposta `'sched', false`. Il `.slx`
flessibile **non ha proprio l'hardware** per la schedula: non e' selezionabile,
va ri-autorata la build. Quindi il ramo `sched` della RHS flessibile e' capacita'
morta rispetto al modello Simulink. Coerentemente, il README elenca fra i
follow-up rimandati proprio "co-scheduled gains **and** notch in one flexible
run".

---

## Il notch tempo-variante (righe 125-141)

```matlab
A('G2zDw','built-in/Gain',[...],'Gain','2*notch_zD');
A('Gcout','built-in/Gain',[...],'Gain','2*(notch_zN-notch_zD)');
W('omega/1','G2zDw/1'); W('omega/1','Gcout/1');
...
A('xn2d','built-in/Sum',[...],'Inputs','--+', ...);
A('v_sum','built-in/Sum',[...],'Inputs','++', ...);
```

E' il pezzo piu' interessante del file. Il notch **non** e' un blocco di libreria:
e' realizzato a mano in forma canonica di controllabilita', con i coefficienti
presi dalle lookup su `omega(t)`.

- Righe 126-128: due `Gain` costanti applicati al segnale `omega(t)`, che
  producono i coefficienti tempo-varianti `2*zD*omega(t)` e
  `2*(zN - zD)*omega(t)`.
- Righe 129-133: `Pn1 = omega^2 * xn1` (usa la lookup `omega2`),
  `Pn2 = 2*zD*omega * xn2`, `Sum '--+'` con `u_pd`, poi la catena
  `int_xn2 -> int_xn1`. Quindi

      xn1_dot = xn2
      xn2_dot = -omega^2 * xn1 - 2*zD*omega * xn2 + u_pd

- Righe 134-141: `Pn3 = 2*(zN - zD)*omega * xn2`, e `v_sum '++'`:

      v = u_pd + 2*(zN - zD)*omega * xn2

Corrisponde riga per riga alle righe 31-33 di `ode_lpv_flex.m`.

### Derivazione: perche' questa realizzazione E' il notch

Congelando `omega = wn`, la realizzazione e'

    A = [ 0        1          ]   B = [0]
        [ -wn^2   -2*zD*wn    ]       [1]
    C = [ 0   2*(zN - zD)*wn ]   D = 1

Il denominatore e' `den(s) = s^2 + 2*zD*wn*s + wn^2` e si ha
`(sI - A)^-1 * B = [1; s] / den(s)`. Quindi

    v/u = C*(sI-A)^-1*B + D
        = 2*(zN - zD)*wn*s / den(s)  +  1
        = ( s^2 + 2*zD*wn*s + wn^2 + 2*zN*wn*s - 2*zD*wn*s ) / den(s)
        = ( s^2 + 2*zN*wn*s + wn^2 ) / ( s^2 + 2*zD*wn*s + wn^2 )

che e' esattamente il notch classico: stessa `wn` a numeratore e denominatore
(guadagno unitario in DC e ad alta frequenza), `zN` piccolo per scavare, `zD`
grande per allargare. La profondita' e'

    |N(j*wn)| = (2*j*zN*wn^2) / (2*j*zD*wn^2) = zN/zD

Con `notch_zN = 0.002` e `notch_zD = 0.7` (`init_simulink_lpv`, riga 116):
`zN/zD = 2.857e-3`, cioe' **-50.9 dB** al centro. E' il "deep notch" di HM3.

### Il punto delicato: notch LTV, non "notch con parametro che cambia"

Va detto con onesta': quando `wn = wn(t)`, la funzione di trasferimento **non
esiste** piu'. Cio' che il modello realizza e' una **realizzazione di stato
tempo-variante** che, congelata a ogni istante, ha come TF il notch centrato su
`omega(t)`. Due conseguenze:

1. **Il risultato dipende dalla realizzazione.** Forme canoniche diverse danno la
   stessa TF congelata ma sistemi LTV **diversi**. Il `.slx` e `ode_lpv_flex.m`
   usano la stessa realizzazione, quindi coincidono fra loro -- ma "il notch
   variabile" e' definito *da questa scelta*, non da un oggetto matematico unico.
2. **Mancano i termini in `omega_dot`.** Una conversione rigorosa da forma
   congelata a LTV genererebbe termini proporzionali a `d(wn)/dt`. Il codice li
   ignora silenziosamente. La giustificazione (mai scritta nel codice) e'
   quasi-stazionaria: `omega` passa da 16.5 a 30.1 rad/s in 140 s, cioe'
   `omega_dot ~ 0.1 rad/s^2`, del tutto trascurabile rispetto a
   `omega^2 ~ 270-910 rad^2/s^2`.

### Notch a conoscenza perfetta

Il modello alimenta la lookup `omega2` **sia** al bending (`Pe1`, riga 92) **sia**
al notch (`Pn1`, riga 136): centro del notch e frequenza vera del modo sono
letteralmente **la stessa tabella**. E' un notch a conoscenza perfetta: nessun
errore di stima di `omega`, nessun disallineamento fra modello e realta'. E' il
caso migliore possibile e va presentato come tale.

Inoltre: il caso **notch fisso** (quello che nel README diverge dopo `t ~ 85 s`)
**non e' realizzato in questo `.slx`**. Il modello costruisce solo la versione
variabile. Il confronto fisso-vs-variabile esiste unicamente nel track script
(`main_flex.m` / `ode_lpv_flex.m`, che accetta un `M.fwn` qualsiasi). Il
`.slx` non ha uno switch per riprodurlo.

---

## TVC (righe 143-147)

```matlab
A('TVC','built-in/TransferFcn',[...], ...
  'Numerator','tvc_num','Denominator','tvc_den');
W('v_sum/1','TVC/1');
W('TVC/1','P3/2'); W('TVC/1','P7/2'); W('TVC/1','Pe3/2');
```

- Riga 144: il TVC e' un blocco `Transfer Fcn` **LTI**, i cui coefficienti
  arrivano da `init_simulink_lpv` (righe 112-114): `Wtvc = build_tvc(p0, 3)`,
  cioe' servo + **approssimante di Pade di ordine 3** del ritardo di trasporto
  (20 ms), poi `tfdata(tf(Wtvc),'v')`.
- Riga 147: l'uscita `delta` (deflessione fisica) va nel plant rigido (`P3`,
  `P7`) e nel bending (`Pe3`). Il terzo consumatore, il log (`log_delta`), e'
  cablato solo piu' avanti, alla riga 158 -- il commento di riga 146 li elenca
  tutti e tre, il cablaggio della riga 147 ne fa due su tre.

Onesta' obbligatoria su questo blocco:

- **Il TVC e' congelato a max-q**: `p0 = load_hw3_params()` e' valutato a
  `t_ref = 72 s`. In un modello dove *tutto il resto* e' tempo-variante,
  l'attuatore e il ritardo restano LTI. E' una limitazione reale, non dichiarata
  nell'header.
- **Realizzazioni diverse fra i due track.** Il `.slx` usa un `Transfer Fcn`
  (che internamente costruisce la forma canonica di controllabilita' dai
  coefficienti `tvc_num`/`tvc_den`), mentre `ode_lpv_flex.m` (righe 34-35) usa la
  quadrupla `At, Bt, Ct, Dt` presa da `ssdata(ss(Wtvc))`. Sono equivalenti in
  ingresso-uscita ma hanno **stati diversi**: si possono confrontare `delta` e le
  uscite, non gli stati interni del TVC. Il passaggio `ss -> tf -> canonica` su un
  sistema di ordine 5 (servo 2 + Pade 3) e' anche numericamente il punto piu'
  fragile della catena.
- `run_flex_simulink.m` (riga 42) ricostruisce `delta_ode = Ct*x_tvc + Dt*0` con
  il commento "(Dt=0)": e' un **no-op aritmetico** che documenta l'assunzione che
  il TVC sia strettamente proprio. L'assunzione e' plausibile (servo di ordine 2
  a grado relativo 2 in cascata a un Pade a grado relativo 0 -> grado relativo 2),
  ma **non e' verificata da nessuna asserzione**: se `Dt` fosse diverso da zero,
  quel confronto sarebbe semplicemente sbagliato.
- Nessuna **saturazione** e nessun **rate limit** sulla deflessione, ne' qui ne'
  nel modello rigido.

---

## La disciplina dei coefficienti effettivi, rotta (discrepanza reale)

Questa e' la cosa piu' importante della pagina.

Nel **track rigido** il `.slx` legge `lpv_c1..lpv_c7` e la RHS
`ode_lpv_ascent.m` (righe 26-28) legge `M.fc1(t)..M.fc7(t)`, cioe' le
`griddedInterpolant` costruite da `init_simulink_lpv` (righe 102-103) **sulle
stesse identiche array**. Le due catene interpolano gli stessi numeri: integrano
la stessa funzione.

Nel **track flessibile questo non accade**:

| grandezza | cosa interpola il `.slx` | cosa interpola `ode_lpv_flex.m` |
|-----------|--------------------------|----------------------------------|
| coeff. di `theta` in `zddot` | `L{a1*V + a4}` (tabella `lpv_c2`, riga 46) | `L{a1}*L{V} + L{a4}` (righe 19, 38) |
| coeff. di `alpha_w` in `zddot` | `L{a1*V}` (`lpv_c4`) | `L{a1}*L{V}` (riga 38) |
| coeff. di `zdot` in `thetaddot` | `L{A6/V}` (`lpv_c5`) | `L{A6}/L{V}` (riga 39) |
| rigidezza di bending / notch | `L{omega^2}` (`lpv_omega2`) | `L{omega}^2` (righe 20, 33, 40) |

`ode_lpv_flex.m` riceve `fa1, fa3, fa4, fA6, fK1, fV, fomega` -- gli interpolanti
dei **fattori grezzi** -- e li ricombina *dentro* la RHS. Il `.slx` interpola
invece le tabelle **gia' combinate**.

Poiche' l'interpolazione lineare non commuta con il prodotto, i due lati
integrano **funzioni diverse fra i breakpoint**. Su una cella di ampiezza `h`, con
`u = s/h`, `Df = f1-f0`, `Dg = g1-g0`:

    L{f*g}(u) - L{f}(u)*L{g}(u) = Df*Dg * u*(1-u)

nullo solo sui nodi, massimo `Df*Dg/4` a meta' cella. E' un errore di **modello**:
non si riduce stringendo la tolleranza del solver.

Quantificato sul dataset vero (griglia `lpv_t = 0:1:140`, passo 1 s):

| coefficiente | max |scarto| assoluto | relativo |
|--------------|-------------------------|----------|
| `c2 = a1*V + a4` | 2.24e-3 | 0.0034 % |
| `c4 = a1*V`      | 2.24e-3 | 0.016 %  |
| `c5 = A6/V`      | 4.89e-7 | 0.014 %  |
| `omega^2`        | 6.81e-3 | 0.00075 %|

I coefficienti `c1 = a1`, `c3 = a3`, `c6 = A6`, `c7 = K1` sono **grezzi** e quindi
commutano esattamente: per quelli non c'e' discrepanza. Anche
`2*zBM*omega` e' lineare in `omega`, quindi `L{2*zBM*omega} = 2*zBM*L{omega}`:
commuta anch'esso. La discrepanza riguarda solo le quattro righe della tabella
sopra.

**Cosa si puo' e non si puo' dire.** Si puo' dire, con certezza dal codice, che i
due lati del track flessibile **non integrano la stessa funzione** e che quindi il
residuo del confronto ha un **pavimento di modello** che nessuna tolleranza
elimina. Si puo' osservare che il residuo dichiarato nel README per il caso
flessibile (**5e-7 rad su theta**) e' 2-5 volte quello del caso rigido (1.1e-7
frozen, 2.2e-7 scheduled), dove le tabelle invece coincidono -- il che e'
**coerente** con questa spiegazione. Non si puo' affermare che sia *provato* che
la causa sia questa: non ho eseguito la simulazione, e il modello flessibile ha
anche piu' stati, il notch e il TVC.

La correzione sarebbe banale: far consumare a `ode_lpv_flex.m` gli stessi
`S.fc1..S.fc7` e un `S.fomega2 = gi(omega.^2)` invece di ricombinare `a1,V,a4`
e `omega^2` a mano.

> **Possibile domanda d'esame** -- il modello flessibile Simulink e il suo
> `ode45` di riferimento sono lo stesso sistema?
> *Risposta:* No, non esattamente. Coincidono sui breakpoint della griglia
> (passo 1 s) ma differiscono in mezzo, perche' il `.slx` interpola i
> coefficienti **gia' combinati** (`a1*V+a4`, `A6/V`, `omega^2`) mentre la RHS
> interpola i fattori grezzi e li ricombina. L'interpolazione lineare non commuta
> con il prodotto, quindi le due integrano funzioni diverse. Lo scarto e' piccolo
> (fino a ~0.02 % sui coefficienti) ma e' un errore di modello, non di solver: e'
> il motivo strutturale per cui il residuo flessibile non puo' scendere al livello
> del residuo rigido, dove invece la disciplina delle tabelle combinate e'
> rispettata.

---

## Logging (righe 149-158)

```matlab
tw = {'-1','MaxDataPoints','inf'};
A('log_theta', ... ,'VariableName','theta_sl', ...);
A('log_eta',   ... ,'VariableName','eta_sl', ...);
W('int_th/1','log_theta/1');  W('int_eta/1','log_eta/1');
W('TVC/1','log_delta/1');
```

- Riga 150: `SampleTime = '-1'` (ogni passo maggiore del solver) e
  `MaxDataPoints = 'inf'` (via il tetto di 1000 punti).
- Righe 151-158: sei `To Workspace` in formato `Timeseries`: `theta_sl`, `z_sl`,
  `zdot_sl`, `eta_sl`, `delta_sl`, `alpha_w_sl`.
- **Cosa viene loggato come `theta_sl`**: `int_th/1`, cioe' l'**assetto vero**,
  non la misura `theta_m`. E' la scelta corretta, perche' la replica confronta
  con `x(:,3)` di `ode_lpv_flex`, che e' anch'esso lo stato fisico.
- `delta_sl` viene preso dall'**uscita del TVC** (riga 158), quindi e' la
  deflessione fisica, non il comando `v` a valle del notch.

Osservazioni oneste:

- `z_sl` e `zdot_sl` sono loggati ma `run_flex_simulink.m` non li legge mai
  (riga 30 prende solo `theta_sl`, `eta_sl`, `delta_sl`, `alpha_w_sl`).
- Riga 35: `SignalLogging` e' abilitato con `SignalLoggingName = 'logsout'`, ma
  **nessun segnale e' marcato per il logging** in tutta la build. `logsout` resta
  vuoto e non viene letto. Configurazione morta, identica al modello rigido.
- Come nel rigido, `alpha_w_sl` serve a una cosa sola: `run_flex_simulink.m`
  (riga 36) costruisce la `windfun` della replica `ode45` come
  `griddedInterpolant` sull'`alpha_w` **generato dal modello stesso**. Il vento
  non viene rigenerato: e' quindi un **controllo di autoconsistenza** fra due
  implementazioni della stessa matematica, non una validazione indipendente. Il
  ramo del vento e' comune ai due lati e non e' sotto test.

---

## Salvataggio (righe 160-165)

`save_system` scrive `hm3_full_ascent_flex.slx` accanto allo script. Il modello e'
un artefatto derivato: si rigenera lanciando il builder, e ogni modifica fatta a
mano nel canvas viene distrutta dal `new_system` della build successiva (riga 32).
Vale qui la stessa argomentazione del modello rigido: il `.m` e' l'unica fonte di
verita', il modello e' testo diffabile, e il design non puo' divergere
silenziosamente dal codice. Prezzo pagato: il `.slx` **non e' autosufficiente**
(ogni parametro dei blocchi e' un nome di variabile del base workspace, persino
`StopTime = 'Tstop'`), quindi senza `init_simulink_lpv` non compila nemmeno.

---

## Rigido contro flessibile: il riassunto

| | `build_hm3_full_ascent` | `build_hm3_full_ascent_flex` |
|---|---|---|
| stati | 4 (`z, zdot, theta, thetadot`) | 6 + 2 (notch) + `nt` (TVC) = 13 |
| lookup | 8 | 14 |
| `ExtrapMethod` | default `Linear` su `c1..c7`, `Clip` sui guadagni | `Clip` **ovunque** |
| attuatore | **nessuno** (delta ideale) | `Transfer Fcn` LTI: servo + Pade(3), **congelato a max-q** |
| bending | assente | oscillatore 2o ordine, `omega(t)` da lookup, `zBM` **costante** |
| sensori | stati veri | misure INS contaminate: `theta + sigma*eta`, `z - phi*eta` |
| notch | assente | canonica di controllabilita', centro = `omega(t)` (stessa tabella del plant) |
| guadagni | `Switch` frozen / scheduled (28 punti su `tsched = 5:5:140`) | **solo frozen** a max-q |
| tabelle condivise con `ode45` | **si'** (`c1..c7`) | **no** (la RHS ricombina i fattori grezzi) |
| residuo dichiarato vs `ode45` | 1.1e-7 / 2.2e-7 rad | 5e-7 rad |

---

## Limiti, hack e codice stale (riepilogo onesto)

- **Discrepanza tabelle combinate vs fattori grezzi** fra `.slx` e
  `ode_lpv_flex.m` (`c2`, `c4`, `c5`, `omega^2`): i due lati non integrano la
  stessa funzione fra i breakpoint. E' il difetto piu' sostanziale.
- **TVC LTI congelato a max-q** in un modello per il resto tutto tempo-variante.
- **`zBM` costante** (preso a 72 s) mentre `omega` varia dell'83 %.
- **Notch a conoscenza perfetta**: centro del notch e frequenza vera del modo
  sono la stessa lookup. Nessun errore di stima, nessun disallineamento.
- **Il caso "notch fisso" non e' realizzato nel `.slx`**: esiste solo nel track
  script. Il confronto fisso-vs-variabile del README non e' riproducibile con
  questo modello.
- **Nessuno scheduling dei guadagni** nel `.slx` flessibile, mentre
  `ode_lpv_flex.m` conserva il ramo `sched` (capacita' morta lato Simulink).
- **Termini in `omega_dot` ignorati** nella conversione del notch a LTV;
  giustificabile come quasi-stazionario, ma mai dichiarato.
- **Tolleranze del solver mai impostate**: Simulink a `RelTol` di default (1e-3)
  contro `1e-9` della replica; l'accordo si regge su `MaxStep = 0.02`.
- **`Dt = 0` assunto e non verificato** (`run_flex_simulink.m`, riga 42:
  `+ M.Dt*0`, aritmetica morta usata come commento).
- **Nessuna saturazione, nessun rate limit** su `delta`.
- **`SignalLogging`/`logsout` attivi ma vuoti**; `z_sl`, `zdot_sl` loggati e mai
  letti.
- **Validazione autoconsistente**, non indipendente: il vento e' preso dal modello
  stesso.
- **Copia (non link) del sottosistema del professore**: se `strong_wind.slx`
  cambia, serve un rebuild manuale, senza avvisi.

---

## Possibili domande d'esame

**D: Come si realizza un notch il cui centro varia nel tempo, senza blocchi di
libreria?**
R: Si scrive il notch in forma canonica di controllabilita' e si alimentano i
coefficienti con lookup sul parametro di scheduling. Con
`A = [0 1; -wn^2, -2*zD*wn]`, `B = [0;1]`, `C = [0, 2*(zN-zD)*wn]`, `D = 1` si
ottiene `C*(sI-A)^-1*B + D = (s^2 + 2*zN*wn*s + wn^2)/(s^2 + 2*zD*wn*s + wn^2)`,
cioe' il notch classico. In Simulink diventano due integratori, due `Product`
(con `omega^2` e `2*zD*omega` dalle lookup), un `Sum '--+'` per `xn2_dot` e un
`Sum '++'` per l'uscita. Rendendo `wn` un segnale invece che una costante, si ha
un notch che insegue `omega(t)`.

**D: Che profondita' ha il notch e da cosa dipende?**
R: `|N(j*wn)| = zN/zD`. Con `notch_zN = 0.002` e `notch_zD = 0.7` si ottiene
`2.857e-3`, cioe' **-50.9 dB** al centro. `zN` controlla la profondita', `zD` la
larghezza di banda dell'attenuazione. Numeratore e denominatore condividono la
stessa `wn`, quindi il guadagno e' unitario sia in DC sia ad alta frequenza: il
notch non tocca la banda di controllo ne' l'alta frequenza, scava solo attorno a
`omega`.

**D: Perche' serve un notch variabile e non basta quello di HM3?**
R: Perche' HM3 congela tutto a max-q e centra il notch su `omega(72) = 18.9
rad/s`, ma sull'orizzonte simulato (0-140 s) la frequenza del primo modo spazia
in `[16.48, 30.12] rad/s` (fino a 31.77 rad/s a t = 150 s, sul dataset completo,
che il modello non simula). Salendo, `omega(t)` esce dalla campana del notch e il
guadagno d'anello alla risonanza torna sopra 0 dB: il modo non e' piu'
gain-stabilized e l'anello
(che passa da `delta -> eta -> INS -> PD -> delta`) diventa instabile. Il README
riporta divergenza da `t ~ 75-85 s` con notch fisso.

**D: Ma un notch che conosce esattamente `omega(t)` e' realistico?**
R: No, ed e' un limite da dichiarare. Nel modello la lookup `omega2` alimenta
**sia** il bending **sia** il notch: sono letteralmente la stessa tabella. E' il
caso ideale di conoscenza perfetta. In volo `omega(t)` andrebbe stimata (da
modello strutturale, da massa residua, da identificazione online) e l'errore di
stima ridurrebbe l'attenuazione effettiva. Il codice non modella nessun errore di
stima.

**D: Il modello flessibile Simulink e il suo `ode45` sono davvero lo stesso
sistema?**
R: No, non fra i breakpoint. Il `.slx` interpola i coefficienti **gia' combinati**
(`lpv_c2 = a1*V+a4`, `lpv_c5 = A6/V`, `lpv_omega2 = omega^2`), mentre
`ode_lpv_flex.m` interpola i fattori grezzi e li ricombina nella RHS. Poiche'
l'interpolazione lineare non commuta con il prodotto
(`L{f*g} - L{f}*L{g} = Df*Dg*u*(1-u)`), i due integrano funzioni diverse. Lo
scarto sui coefficienti arriva a ~0.02 % su questa griglia a 1 s; e' un errore di
modello, non di solver, quindi mette un pavimento al residuo indipendente dalle
tolleranze. Nel track **rigido** questa disciplina e' invece rispettata (entrambi
leggono `c1..c7`), ed e' proprio questo che rende possibile l'accordo a 1e-7 rad.

**D: Perche' il bending non compare nelle equazioni di `zddot` e `thetaddot`?**
R: Per ipotesi di modellazione: il modo strutturale e' considerato disaccoppiato
dalla dinamica di corpo rigido, e chiude l'anello **attraverso il controllore**,
non attraverso la fisica. Il percorso e': il TVC eccita il modo (`aqk*delta`), il
modo contamina le misure INS (`theta_m = theta + sigma*eta`), il PD reagisce alla
contaminazione e comanda `delta`. E' esattamente questo anello che il notch deve
tagliare. Se il bending retroagisse anche sul corpo rigido, servirebbero i termini
di accoppiamento inerziale nelle prime due equazioni, che qui sono trascurati.
