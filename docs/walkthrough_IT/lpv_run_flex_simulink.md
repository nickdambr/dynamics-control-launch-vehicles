# HM3/LTV_FULL_ASCENT/run_flex_simulink.m

## Ruolo del file nel progetto

Gemello di `run_full_ascent_simulink.m`, ma per il **veicolo flessibile**:
simula `hm3_full_ascent_flex.slx` e sovrappone il risultato alla baseline
`ode45` sul RHS `ode_lpv_flex`. Produce `figures/flex_simulink_vs_script.png` e
stampa i residui su `theta`, `eta` (coordinata di bending) e `delta`.

Il modello flessibile e' l'architettura completa di HM3 Task 2 sollevata al
caso LPV: plant rigido + **primo modo di bending** + **accoppiamento INS**
(l'inerziale legge un `theta` contaminato dalla flessione) + **TVC con ritardo
di trasporto** (Pade del 3o ordine) + **notch che insegue `omega(t)`**. Sono
13 stati:

    x = [ z  zdot  theta  thetadot  eta  etadot | xn1  xn2 | x_tvc(1:5) ]
          ------ plant rigido ------  -bending-   -notch-   ---- TVC ----

Lo showcase di T008 Goal 1 e': la frequenza del primo modo spazza
`omega(t) in [16.5, 31.8] rad/s` durante l'ascesa, mentre il notch di HM3 e'
centrato una volta per tutte su `omega(72) = 18.9 rad/s`; il notch fisso quindi
si **detuna** e il loop diverge (secondo `main_flex`, instabile da `t ~ 75 s`),
mentre il notch che insegue `omega(t)` regge tutta la salita.

Questo file **non** dimostra quel risultato -- lo dimostra `main_flex.m`. Questo
file dimostra una cosa piu' modesta e piu' ingegneristica: che il modello
Simulink e lo script integrano *lo stesso sistema*. Il README dichiara un
accordo di **5e-7 rad su theta** (contro 1.1e-7 / 2.2e-7 del modello rigido);
sotto spiego perche' e' cinque volte peggiore, e la ragione **non e' solo
numerica**.

Chiamato a mano. Dipende da: `init_simulink_lpv`, `build_hm3_full_ascent_flex`,
`ode_lpv_flex`.

---

## Firma, rebuild, setup (righe 1-27)

```matlab
function out = run_flex_simulink(o)
arguments
    o.rebuild (1,1) logical = false
end
...
S    = init_simulink_lpv();
nt   = size(S.tvc.At, 1);              % ordine del TVC = 5
odeo = odeset('RelTol', 1e-9, 'AbsTol', 1e-11);
```

- Righe 20-22: stessa politica del modello rigido -- il `.slx` e' un artefatto
  derivato, rigenerabile con `build_hm3_full_ascent_flex()`. Vale la stessa
  cautela: il `.slx` e' comunque committato, quindi se si modifica il build
  senza `'rebuild', true` si simula un modello vecchio.
- Riga 24: `init_simulink_lpv()` con `push = true` -> spinge nel base workspace
  anche i dati flessibili (`lpv_omega`, `lpv_omega2`, `lpv_2zBMw`, `lpv_aqk`,
  `lpv_sig`, `lpv_phi`, `notch_zN`, `notch_zD`, `tvc_num`, `tvc_den`), che i
  blocchi del modello risolvono per nome.
- Riga 25: `nt = size(S.tvc.At, 1)`. Il TVC e' `build_tvc(p0, 3)` = attuatore del
  2o ordine (`wTVC = 70 rad/s`, `zTVC = 0.7`) moltiplicato per un **Pade del 3o
  ordine** del ritardo `tau = 20 ms` -> **5 stati**. Il codice non lo scrive: lo
  deduce a runtime dalla dimensione di `At`. Buona pratica -- cambiare
  `padeOrder` in `init_simulink_lpv` non rompe questo file.
- Riga 26: tolleranze strette per la baseline. Come nel caso rigido, il **modello
  Simulink NON le usa**: la `set_param` di `build_hm3_full_ascent_flex`
  (righe 33-35) fissa `ode45`, `Variable-step`, `MaxStep = 0.02`, piu'
  `StartTime`/`StopTime` e il logging, ma **non** tocca `RelTol` ne' `AbsTol` -> il
  modello gira con le tolleranze di default (`RelTol = 1e-3`). L'accuratezza del
  lato Simulink e' fissata dal `MaxStep`, non dalla tolleranza.

---

## Simulazione e replay ode45 (righe 29-43)

```matlab
so = sim(mdl, 'StopTime', num2str(S.Tstop));
th = so.theta_sl;  et = so.eta_sl;  de = so.delta_sl;  aw = so.alpha_w_sl;
tt = th.Time;

M = struct('fa1', S.fa1, 'fa3', S.fa3, 'fa4', S.fa4, 'fA6', S.fA6, ...
           'fK1', S.fK1, 'fV', S.fV, 'fomega', S.fomega, ...
           'faqk', S.faqk, 'fsig', S.fsig, 'fphi', S.fphi, ...
           'windfun', griddedInterpolant(aw.Time, squeeze(aw.Data), ...), ...
           'fwn', S.fomega, 'zN', S.notch.zN, 'zD', S.notch.zD, ...
           'zBM', S.notch.zBM, 'At', S.tvc.At, ..., 'sched', false);
[~, x] = ode45(@(t,x) ode_lpv_flex(t,x,M), tt, zeros(6+2+nt,1), odeo);
delta_ode = M.Ct*x(:, 9:end).' + M.Dt*0;      % delta = Ct*x_tvc  (Dt = 0)
```

- Riga 29: nessun `assignin('base','sched', ...)` -- il modello flessibile **non
  ha lo Switch**: i guadagni sono blocchi Gain fissi (`Kp_th0`, `Kd_th0`,
  `Kp_z0`, `Kd_z0`). Coerente con `'sched', false` alla riga 40. Il messaggio del
  build e' esplicito: *"PD gains FROZEN at max-q; the showcase is notch tracking,
  not scheduling"*. **Limite:** il modello flessibile non puo' esercitare lo
  schedule dei guadagni, e il rigido non puo' esercitare il notch: i due
  contributi (Goal 1 e Goal 2 del ticket T008) non sono mai combinati -- il
  README lo ammette nei "follow-ups deferred".
- Riga 36: **stessa mossa chiave del modello rigido**: la baseline riceve
  l'`alpha_w` *loggato dal modello*, non `S.windfun`. Senza questo, si
  confronterebbero due realizzazioni diverse del vento.
- Riga 37: `'fwn', S.fomega` -> il centro del notch e' `omega(t)`, cioe' il
  **notch variabile**. `main_flex` usa lo stesso RHS passando `@(t) w72` per il
  caso fisso; qui il caso fisso **non e' validabile**, perche' il `.slx`
  costruito da `build_hm3_full_ascent_flex` cabla il notch sui lookup di
  `omega(t)` (righe 125-141 di quel file) e non ha un interruttore. Onestamente:
  la validazione copre **solo** la configurazione variabile.
- Riga 41: `zeros(6+2+nt, 1)` = 13 stati a zero, coerenti con gli
  `InitialCondition = 0` di tutti gli integratori del modello (e con la Transfer
  Fcn del TVC, che parte da stato nullo).
- Riga 42: **`M.Dt*0`** -- `delta` esce dal TVC come `Ct*x_tvc + Dt*v`, ma qui
  il termine di feedthrough viene moltiplicato per **zero** invece che per
  l'ingresso `v`. E' corretto *solo perche'* la funzione di trasferimento del TVC
  e' **strettamente propria**: attuatore del 2o ordine per Pade(3) -> numeratore
  di grado 3, denominatore di grado 5, quindi `D = 0`. Il commento lo dichiara
  (`Dt=0`). Ma e' una scorciatoia: se un giorno il TVC diventasse bi-proprio --
  per esempio aggiungendo in serie una rete di anticipo, o un feedthrough
  nell'attuatore -- questa riga produrrebbe silenziosamente un `delta` sbagliato,
  mentre `ode_lpv_flex` (riga 34, che scrive correttamente `M.Ct*xt + M.Dt*v`)
  resterebbe giusta. Meglio sarebbe ricostruire `v` e usarlo. Attenzione pero' a
  non citare `padeOrder` come esempio: `build_tvc` lo vincola a essere positivo
  (riga 12, `mustBePositive`) e comunque nessun ordine di Pade puo' rendere
  bi-proprio il TVC -- per un Pade di ordine n il numeratore ha grado n e il
  denominatore grado n+2 (i due poli dell'attuatore), quindi `D = 0` sempre.

---

## Il notch variabile: perche' e' realizzato "a mano" (contesto, `ode_lpv_flex` righe 30-33 e `build_hm3_full_ascent_flex` righe 125-141)

Entrambi i track realizzano il notch in **forma canonica di controllabilita'**:

    xn1_dot = xn2
    xn2_dot = -omega^2 * xn1 - 2*zD*omega * xn2 + u_pd
    v       = u_pd + 2*(zN - zD)*omega * xn2

Verifica che sia davvero un notch: da `xn1_ddot + 2*zD*omega*xn1_dot +
omega^2*xn1 = u_pd` segue `Xn1(s) = U(s)/(s^2 + 2*zD*omega*s + omega^2)` e
`xn2 = s*Xn1`, quindi

    V(s)/U(s) = 1 + 2*(zN - zD)*omega*s / (s^2 + 2*zD*omega*s + omega^2)
              = (s^2 + 2*zN*omega*s + omega^2)
                / (s^2 + 2*zD*omega*s + omega^2)

che e' esattamente il notch standard: con `zN = 0.002` (numeratore quasi non
smorzato) si scava un buco profondo a `omega`, con `zD = 0.7` (denominatore ben
smorzato) si evita di introdurre risonanza. Questa e' la **stessa** forma usata
nei due track -- e questo e' importante:

> **Possibile domanda d'esame** -- Perche' il notch e' costruito con integratori
> e prodotti invece di usare un blocco Transfer Fcn con coefficienti variabili?
> *Risposta:* Perche' il blocco Transfer Fcn di Simulink accetta coefficienti
> **costanti**: non si puo' far variare `omega` a runtime. Serve quindi una
> realizzazione a stati, con i coefficienti presi da lookup su `omega(t)` e
> moltiplicati per gli stati con blocchi Product. Ma c'e' un punto teorico piu'
> profondo: quando i coefficienti variano nel tempo, **la funzione di
> trasferimento non definisce piu' univocamente il sistema** -- realizzazioni
> diverse della stessa TF "congelata" danno sistemi LTV diversi (compaiono i
> termini in `d omega/dt`). Sia `ode_lpv_flex` sia il modello usano la *stessa*
> forma canonica di controllabilita', quindi coincidono fra loro; ma nessuno dei
> due modella i termini derivativi di `omega(t)`: e' una **sostituzione a
> coefficienti congelati**, un'approssimazione valida perche' `omega(t)` varia
> lentamente rispetto alla banda del notch. Va detto: e' una approssimazione, non
> un'identita'.

Il TVC, invece, e' **LTI**, e per questo puo' essere un blocco Transfer Fcn nel
modello (`tvc_num`/`tvc_den`) e uno `ss` nella baseline (`At, Bt, Ct, Dt`): due
realizzazioni diverse dello stesso sistema tempo-invariante hanno lo stesso
comportamento ingresso-uscita, quindi lo stato interno puo' essere diverso senza
conseguenze. La liberta' che si ha con il TVC e' esattamente quella che **non**
si ha con il notch.

---

## Residui e figura (righe 45-67)

```matlab
err = struct('theta', max(abs(x(:,3) - squeeze(th.Data))), ...
             'eta',   max(abs(x(:,5) - squeeze(et.Data))), ...
             'delta', max(abs(delta_ode - squeeze(de.Data))));
```

- Righe 45-47: si confrontano `theta` (stato 3), `eta` (stato 5, la coordinata
  modale del bending) e `delta`. **`z` non viene confrontata**, pur essendo
  loggata dal modello (`z_sl`): una lacuna gratuita, il dato c'e' gia'.
  Confrontare `eta` e' la scelta giusta: e' la variabile che diverge quando il
  notch e' scentrato, quindi e' il test piu' sensibile.
- Righe 56-61: overlay `ode45` continuo / Simulink tratteggiato, tre tile.

### Perche' 5e-7 e non 1e-7: la discrepanza vera fra i due track

Il modello rigido raggiunge ~1e-7 rad, il flessibile si ferma a ~5e-7. Il codice
non lo spiega; leggendolo, la ragione principale e' **strutturale, non
numerica**:

- `hm3_full_ascent_flex.slx` usa gli **stessi lookup di coefficienti efficaci**
  del modello rigido: `c1..c7`, `invV` (`build_hm3_full_ascent_flex`, righe
  45-51 e 75-83). Cioe' tabelle di `a1*V + a4`, `a1*V`, `A6/V`, `1/V`.
- `ode_lpv_flex` invece prende i coefficienti **grezzi** (righe 19-20:
  `fa1, fa3, fa4, fA6, fK1, fV`) e ricompone i prodotti **a runtime**
  (righe 38-39):

      zdd  = a1*zdot + (a1*V + a4)*theta + a3*delta - a1*V*aw
      thdd = (A6/V)*zdot + A6*theta + K1*delta - A6*aw

Ma **l'interpolazione lineare non commuta con prodotto e divisione**: fra due
breakpoint (qui distanti 1 s) `interp(a1*V) != interp(a1)*interp(V)`. Misurando
sul dataset in repo (massimo scarto relativo, valutato ai punti medi fra i
breakpoint sull'orizzonte 0-140 s):

| termine | modello (`.slx`) | baseline (`ode_lpv_flex`) | scarto relativo |
|---|---|---|---|
| coeff. di `theta` in `zddot` | `interp(a1*V + a4)` | `interp(a1)*interp(V) + interp(a4)` | ~3.9e-5 |
| coeff. di `zdot` in `thetaddot` | `interp(A6/Vsafe)` | `interp(A6)/interp(V)` | ~6.8e-4 |
| rigidezza bending / notch | `interp(omega^2)` | `interp(omega)^2` | ~1.9e-5 |

(Nella riga di mezzo c'e' anche una seconda differenza, di second'ordine: le
tabelle del `.slx` proteggono la divisione con `Vsafe = max(V, 1)` --
`init_simulink_lpv`, righe 59 e 64-68, perche' `V(0) = 0` -- mentre
`ode_lpv_flex`, righe 38-39, usa la `V` grezza. La discrepanza vive solo nel
primo secondo di volo, dove `V < 1 m/s` e il veicolo non e' ancora eccitato.)

I due lati stanno quindi integrando **sistemi LTV leggermente diversi**, e il
residuo non puo' scendere sotto quell'ordine di grandezza. Con
`theta` di picco ~1.7e-2 rad, lo scarto relativo di ~4e-5 sul coefficiente di
`theta` vale ~7e-7 rad -- proprio la scala del residuo osservato. Il modello
rigido non ha questo problema
perche' `ode_lpv_ascent` (righe 26-28) usa `fc1..fc7`, cioe' **gli interpolanti
delle stesse tabelle** che finiscono nei lookup.

Non e' un bug che invalida i risultati (gli scarti relativi stanno fra ~4e-5 e
~7e-4, molto sotto qualunque incertezza fisica del modello), ma **e' una
incoerenza fra i due track**: il residuo del caso flessibile **non** e'
spiegabile come "solver tolerance". Attenzione a citare bene il README: quella
formula la usa (a ragione) solo per l'overlay del modello **rigido** (riga 64,
"frozen 1.1e-7, scheduled 2.2e-7 -- solver tolerance"), mentre per il flessibile
(riga 139) si limita a riportare i 5e-7 rad senza attribuirli. La spiegazione
strutturale qui sopra e' quella che manca. La correzione sarebbe di una riga: far leggere a
`ode_lpv_flex` gli stessi `fc1..fc7` (e `fomega2`) invece dei coefficienti
grezzi. Nota che `2*zBM*omega` e' **lineare** in `omega`, quindi la tabella
`lpv_2zBMw` e il prodotto a runtime coincidono esattamente: la discrepanza vive
solo nei termini non lineari nei coefficienti.

- Termine di smorzamento del bending: `-2*zBM*omega*etadot` con `zBM = 0.005`
  (dalla tabella di HM3): il modo e' **quasi non smorzato** -- ecco perche' un
  notch scentrato lo lascia crescere senza limite invece di vederlo estinguersi.

---

## Cosa deve essere identico perche' l'overlay funzioni (riepilogo)

1. **Stesso vento** -> garantito per costruzione (l'`alpha_w` loggato dal modello
   entra come `windfun` nella baseline, riga 36). Residuo: interpolazione lineare
   fra i campioni loggati.
2. **Stessa realizzazione del notch** -> si', entrambi in forma canonica di
   controllabilita' con centro `omega(t)`. Se uno dei due usasse una
   realizzazione diversa (es. forma modale, o `tf` congelata ricostruita a ogni
   passo), le due traiettorie **divergerebbero** anche a solver perfetto.
3. **Stessa realizzazione del TVC** -> non serve: e' LTI, basta la stessa
   funzione di trasferimento.
4. **Stesse condizioni iniziali** -> 13 zeri contro tutti gli integratori a zero.
5. **Stessi coefficienti fra i breakpoint** -> **NO**: e' l'unico punto in cui
   i due track divergono (vedi tabella sopra). E' il limite principale di questa
   validazione.
6. **Stesso solver** -> stessa famiglia (ode45), **tolleranze diverse**;
   l'accuratezza del lato Simulink viene da `MaxStep = 0.02`.

---

## Possibili domande d'esame

**D: Come si realizza un filtro (il notch) con parametri variabili in Simulink?**
R: Non con un blocco Transfer Fcn (accetta solo coefficienti costanti), ma con
una realizzazione a stati: due integratori in cascata, i coefficienti
`omega(t)^2` e `2*zD*omega(t)` presi da lookup table 1-D sul tempo di volo e
moltiplicati per gli stati con blocchi Product, e l'uscita
`v = u + 2*(zN - zD)*omega(t)*xn2`. Congelando `omega`, la TF ingresso-uscita e'
`(s^2 + 2*zN*omega*s + omega^2)/(s^2 + 2*zD*omega*s + omega^2)`. E' una
sostituzione a coefficienti congelati: i termini in `d omega/dt` non sono
modellati.

**D: Perche' il notch fisso fa divergere il veicolo, mentre il modello rigido
con i guadagni congelati resta stabile per tutta l'ascesa?**
R: Perche' sono due meccanismi diversi. Il PD congelato di max-q resta adeguato
perche' max-q e' l'istante **peggiore** (l'instabilita' aerodinamica `A6` e' li'
massima), quindi il design dimensionato li' e' conservativo altrove. Il notch,
invece, non e' un margine: e' una **cancellazione in frequenza**. Se
`omega(t)` esce dalla banda stretta del notch (`zN = 0.002` -> notch molto
stretto), il modo di bending non e' piu' attenuato e, essendo quasi non
smorzato (`zBM = 0.005`), va instabile per feedback attraverso l'INS (che legge
`theta_m = theta + sigma*eta`) e il TVC (che eccita il modo con `aqk*delta`).
Un design "robusto" in ampiezza non protegge da un errore di sintonia in
frequenza.

**D: L'accordo del caso flessibile e' 5e-7 rad, cinque volte peggiore del rigido.
Colpa del solver?**
R: Solo in parte. La causa dominante e' che i due track **non usano le stesse
tabelle di coefficienti**: il `.slx` usa i coefficienti efficaci gia' combinati
(`interp(a1*V + a4)`), mentre `ode_lpv_flex` interpola i coefficienti grezzi e
li moltiplica a runtime (`interp(a1)*interp(V) + interp(a4)`). Fra un breakpoint
e l'altro sono funzioni diverse (scarto relativo ~4e-5 sul coefficiente di
`theta`, fino a ~7e-4 su quello di `zdot` in `thetaddot`), quindi i due lati
integrano sistemi LTV leggermente diversi. Il modello rigido non ha il problema
perche' li' anche il RHS usa `fc1..fc7`.

**D: Cosa NON e' validato da questo file?**
R: (1) La configurazione a **notch fisso** -- il `.slx` implementa solo il notch
variabile, quindi il risultato piu' spettacolare (la divergenza) e' dimostrato
solo dallo script, non replicato in Simulink. (2) I **guadagni schedulati** sul
plant flessibile: nel modello flessibile il PD e' congelato per costruzione.
(3) La coordinata `z`, loggata ma non confrontata. (4) Il generatore di vento
(usato come ingresso comune ai due lati). (5) Un solo scenario di vento, nessun
Monte Carlo.

**D: A cosa serve `nt = size(S.tvc.At,1)`?**
R: A dimensionare il vettore di stato della baseline (`6 + 2 + nt = 13`) senza
cablare a mano l'ordine del TVC. L'ordine dipende da `padeOrder` (3 in
`init_simulink_lpv`, riga 112) e dall'attuatore del 2o ordine: se si cambia
l'ordine dell'approssimazione di Pade, questo file continua a funzionare. E'
anche il motivo per cui `delta_ode` estrae `x(:, 9:end)`: gli stati del TVC sono
in coda.
