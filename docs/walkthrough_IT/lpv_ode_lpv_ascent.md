# HM3/LTV_FULL_ASCENT/ode_lpv_ascent.m

## Ruolo del file nel progetto

Questo file e' il **cuore numerico** dell'estensione LPV di HM3: e' la
right-hand side (RHS) che `ode45` valuta per propagare il modello rigido di
beccheggio **lungo tutta l'ascesa (0-140 s)**, non piu' congelato a max-q.
La traccia dell'homework chiede solo il punto di progetto a t = 72 s
(pressione dinamica massima), dove il sistema e' LTI e si puo' fare tutta
l'analisi classica in frequenza (Nichols, margini, notch). Qui invece i
coefficienti del plant diventano funzioni del tempo e il sistema e' un LTV
(o LPV, se si legge il tempo come parametro di scheduling): non esiste piu'
una funzione di trasferimento, e l'unico modo onesto di sapere che cosa fa il
veicolo e' **integrare davvero le equazioni**.

Il file e' il gemello dinamico di `HM3/build_plant_rigid.m`. Quello costruisce
la quadrupla (A, B, C, D) con i coefficienti valutati a un istante fissato;
questo prende le stesse identiche equazioni e sostituisce ogni coefficiente
con un interpolante nel tempo, richiamato **a ogni valutazione della RHS**.
La struttura `M` che riceve in ingresso e' preparata da
`HM3/LTV_FULL_ASCENT/init_simulink_lpv.m` (griglie + `griddedInterpolant`) e
impacchettata dalla helper locale `make_model` in
`HM3/LTV_FULL_ASCENT/main_full_ascent.m` (righe 123-134).

Dentro la RHS e' chiuso **anche l'anello di controllo**: non c'e' un blocco
"controllore" separato, la legge PD e' scritta in linea (riga 22). Questo
rende il file l'unico punto in cui coesistono plant, vento e controllore --
per questo `main_full_ascent.m` lo dichiara "source of truth" e il modello
Simulink `hm3_full_ascent.slx` (costruito da `build_hm3_full_ascent.m`) viene
validato **contro** di esso, non viceversa.

Chi lo chiama -- sono tre i punti di ingresso, non uno:
`main_full_ascent.m` (righe 31-32, due integrazioni: guadagni congelati e
guadagni schedulati; piu' la verifica di consistenza a t = 72 s, righe 54-61),
`run_full_ascent_simulink.m` (riga 49: il replay `ode45` sullo stesso vento
prodotto dal modello -- e' esattamente la validazione del `.slx` di cui sopra) e
`main_q_scheduling.m` (riga 48: lo studio T008 sullo scheduling in q(t) invece
che in t, che lavora sul plant rigido e quindi riusa questa stessa RHS). Il file
**non** e' usato dallo studio flessibile: quello ha la sua RHS,
`ode_lpv_flex.m`.

---

## Firma e contratto (righe 1-13)

```matlab
function dx = ode_lpv_ascent(t, x, M)
% LTV rigid pitch-plane RHS, full-ascent LPV baseline.
% ...
% No arguments validation by design: ode45 inner loop.
```

- Riga 1: firma `(t, x, M)`. Lo stato e' `x = [z; zdot; theta; thetadot]`
  (4 stati): `z` [m] deriva laterale rispetto alla traiettoria di
  riferimento, `theta` [rad] assetto di beccheggio. L'uscita `dx` e' 4x1.
- Righe 2-9: la docstring dichiara esattamente cosa e' dentro `M`:
  `fc1..fc7` (interpolanti dei coefficienti), `windfun`, `fKp`/`fKd`
  (schedule dei guadagni), i guadagni congelati e il flag `sched`.
- Righe 11-12: **nessun blocco `arguments`, per scelta**. Questa funzione sta
  nel loop interno di `ode45` e viene chiamata dell'ordine di 1e5 volte per
  simulazione; un blocco di validazione la' dentro costerebbe piu' della
  matematica che esegue. La validazione vive al confine del run, in
  `init_simulink_lpv` (che ha il suo `arguments`, righe 32-37). E' la stessa
  convenzione dichiarata nel `CLAUDE.md` del repo per ODE RHS, residui di
  shooting e callback di `fmincon`.

---

## Perche' il sistema e' LPV (e non solo "non lineare")

Le equazioni sono **lineari nello stato**, ma i coefficienti dipendono dal
tempo. Il modello rigido di HM3 (Eq. 1 della traccia, senza bending) e':

    zddot     = a1*zdot + (a1*V + a4)*theta + a3*delta - a1*V*alpha_w
    thetaddot = (A6/V)*zdot + A6*theta      + K1*delta - A6*alpha_w

Chi varia lungo l'ascesa, e perche':

- **A6 = mu_alpha** [1/s^2 per rad] -- momento aerodinamico. E' il termine
  **destabilizzante**: il sottosistema di solo assetto e' `thetaddot = A6*theta`,
  che ha poli in `+/- sqrt(A6)`. Con A6 > 0 (centro di pressione davanti al
  baricentro) c'e' sempre un polo reale positivo: il razzo e' staticamente
  instabile in aria e **non puo' volare senza controllo attivo**. A6 e'
  proporzionale alla forza normale aerodinamica, quindi cresce con la
  pressione dinamica: il suo massimo cade nella regione di max-q, ed e'
  esattamente per questo che HM3 progetta li'.
- **K1 = mu_c** [1/s^2 per rad] -- efficacia del TVC (spinta per braccio
  diviso inerzia). Varia perche' variano spinta, massa e inerzia mentre il
  primo stadio consuma propellente.
- **V(t)** [m/s] -- velocita' relativa. Entra tre volte: in `a1*V`, in `A6/V`
  e nel termine di vento. Parte da 0 al lift-off (vedi la guardia `Vsafe`,
  piu' sotto) e cresce fino a valori supersonici.
- **a1, a3, a4** -- coefficienti della dinamica traslazionale
  (a1 ~ -N_alpha/(m*V), a1*V + a4 ~ -(N_alpha + T - D)/m). Anche loro sono
  tabulati nel tempo nel dataset.
- **q_bar(t)** (campo `Q` del dataset) non entra nella RHS, ma serve a valle
  come indicatore di carico (`main_full_ascent.m`, riga 154).

Conseguenza concreta: **il plant congelato e' una comodita' di progetto, non
il sistema vero**. Congelare significa scrivere A(t0) e fare finta che
A(t) = A(t0) per sempre. Propagare la dinamica vera vuol dire chiamare
`A(t)` a ogni step del solver -- ed e' quello che fanno le righe 25-28.

> **Possibile domanda d'esame** -- perche' il segno di A6 e' cosi' importante?
> *Risposta:* A6 > 0 significa che il momento aerodinamico amplifica la
> perturbazione di assetto (cp davanti al cg): il polo `+sqrt(A6)` e' in
> semipiano destro. La costante di tempo di divergenza e' `1/sqrt(A6)`, che a
> max-q e' dell'ordine del mezzo secondo: il controllo deve avere banda
> nettamente superiore. Se A6 fosse negativo il veicolo sarebbe
> aerodinamicamente stabile (come una freccia) e il problema di controllo
> sarebbe un altro problema.

---

## Vento: `aw = M.windfun(t)` (riga 14)

```matlab
aw = M.windfun(t);   % wind angle of attack alpha_w(t)
```

- Riga 14: il vento **non e' un forzamento tabulato a parte, e' un ingresso
  del plant valutato dentro la RHS**. `M.windfun` e' un `griddedInterpolant`
  costruito in `init_simulink_lpv.m` alla riga 107 su `wg.t`, la griglia
  temporale prodotta dal generatore del professore
  (`General/hw3-v3/strong_wind.slx`, eseguito una volta sola alle righe 87-90,
  con la sua media `v_wp` piu' la turbolenza di Dryden).
- Il segnale interpolato **non e' la velocita' di vento**, e'
  `alpha_w = v_w / V` gia' adimensionalizzato (righe 89-90 di
  `init_simulink_lpv`: la riga 89 costruisce `Vwg = max(interp1(V), 1)` a
  proteggere la divisione, la riga 90 forma `alphaw = wg.vw ./ Vwg`). Quindi
  `alpha_w` e' un **angolo [rad]**, e la RHS lo tratta come tale.
- Convenzione di segno: nella RHS il vento entra con il **meno** in entrambi i
  canali (`- M.fc4(t)*aw`, `- M.fc6(t)*aw`). Non e' arbitrario: raccogliendo
  il canale di beccheggio si legge

      thetaddot = A6*(theta + zdot/V - alpha_w) + K1*delta

  cioe' l'angolo d'attacco aerodinamico **effettivo** e'
  `alpha_aero = theta + zdot/V - alpha_w`. E' fisicamente corretto: conta la
  velocita' **relativa** all'aria, e la componente laterale relativa e'
  `(zdot - v_w)`; dividendo per V si ottiene proprio `zdot/V - alpha_w`.
  Lo stesso raccoglimento funziona sul canale traslazionale, perche'
  `c2 = a1*V + a4` e `c4 = a1*V` condividono il fattore `a1*V`.

**Incoerenza da segnalare.** L'indicatore di carico in
`main_full_ascent.m`, riga 153, calcola
`alpha = x(:,3) + x(:,2)./V + S.windfun(tt)`, cioe' con il **piu'** davanti a
`alpha_w`, e il README della cartella scrive la stessa formula
(`alpha_total = theta + z_dot/V + alpha_w`). Questa non e' la stessa alpha che
guida il plant: le due convenzioni non possono essere entrambe giuste, e la
RHS (che e' la dinamica) usa il meno. Il grafico `q_bar*alpha` va quindi letto
con cautela -- il modulo del picco cambia a seconda del segno scelto.

---

## Il controllore chiuso dentro la RHS (righe 17-22)

```matlab
if M.sched
    Kp = M.fKp(t);  Kd = M.fKd(t);
else
    Kp = M.Kp_th0;  Kd = M.Kd_th0;
end
delta = -(Kp*x(3) + Kd*x(4) + M.Kp_z*x(1) + M.Kd_z*x(2));
```

- Righe 17-21: l'unico "if" della funzione. Il flag `M.sched` sceglie fra i
  **guadagni congelati** (il PD progettato una volta a max-q,
  `init_simulink_lpv.m` righe 71-72) e la **schedule** `Kp(t), Kd(t)`
  (righe 74-85: un PD ridisegnato per ogni punto della griglia `tsched`).
- **Attenzione: il warm-start della schedule non esiste davvero.** I commenti di
  `init_simulink_lpv.m` promettono una continuation (riga 74, "warm-started
  (continuation)"; riga 84, `Kprev = [Kk.Kp_th Kk.Kd_th]; % continuation`) e la
  riga 82 passa effettivamente `'K0', Kprev` a `design_controller`. Ma
  `design_controller.m` dichiara `o.K0 (1,2) ... = [0 0]   % accepted, unused`
  (riga 32) e la docstring lo dice esplicitamente: "K0 ignored (kept for call
  compatibility)" (riga 19). Il tuner riparte **sempre** dallo stesso punto
  iniziale in forma chiusa di D'Antuono, `x0 = log([2*A6/K1, sqrt(A6)/K1])`
  (riga 55): ogni nodo della griglia e' un `fminsearch` indipendente. Il
  commento nel sorgente e' stale -- da sapere, se all'orale il professore chiede
  "e la continuation dove la vedo?".
- Riga 22: legge di controllo. `theta_ref = 0` non compare esplicitamente:
  e' assorbito nel segno meno (regolazione a zero). La struttura coincide
  esattamente con quella di `HM3/assemble_loop.m` (riga 25), quindi la RHS e
  l'analisi in frequenza chiudono **lo stesso anello**: e' questo che rende
  legittimo confrontare margini congelati e simulazione LTV.
- **Attuatore ideale**: `delta` calcolato alla riga 22 va *direttamente* nel
  plant. Niente TVC, niente ritardo di 20 ms, niente notch, nessuna
  saturazione ne' limite di rate. E' la stessa ipotesi del Task 1 di HM3
  (`assemble_loop(G, K, [])`). Tutta la parte "vera" dell'attuazione compare
  solo in `ode_lpv_flex.m`.
- **I guadagni di drift non sono mai schedulati**: `M.Kp_z` e `M.Kd_z` sono
  costanti anche quando `sched = 1` (vedi `make_model`, riga 133 di
  `main_full_ascent.m`, che passa sempre `S.K0.Kp_z`). Inoltre non sono
  nemmeno *tuned*: sono i default di `design_controller.m` (righe 28-29,
  `Kp_z = Kd_z = -1e-3`). Il progetto ottimizza solo la coppia
  `(Kp_th, Kd_th)`.
- Segno: essendo `Kp_z` negativo, il termine `-(M.Kp_z*x(1))` e' **positivo**
  in `z`. E' una scelta del progetto originale, non un errore di questo file.

---

## Il plant tempo-variante (righe 25-28)

```matlab
dx = [ x(2);
       M.fc1(t)*x(2) + M.fc2(t)*x(3) + M.fc3(t)*delta ...
           - M.fc4(t)*aw;
       x(4);
       M.fc5(t)*x(2) + M.fc6(t)*x(3) + M.fc7(t)*delta ...
           - M.fc6(t)*aw ];
```

- Righe 25-28: e' riga per riga la matrice A di `build_plant_rigid.m`
  (righe 11-14) con `p.a1 -> fc1(t)`, `p.a1*p.V + p.a4 -> fc2(t)`,
  `p.A6/p.V -> fc5(t)`, `p.A6 -> fc6(t)`, e le colonne di ingresso
  `Bd = [0; a3; 0; K1] -> [fc3; fc7]`,
  `Bw = [0; -a1*V; 0; -A6] -> [-fc4; -fc6]` (righe 16-17 di
  `build_plant_rigid.m`). **La corrispondenza e' esatta**: e' questo che
  garantisce che, congelando gli interpolanti a 72 s, si ritrovi il modello
  di HM3 -- ed e' proprio il test di consistenza di `main_full_ascent.m`
  (righe 49-61, `err_consistency` dell'ordine di 1e-10 rad secondo il README).
- Nota di economia: `fc6` e' usato **due volte** (coefficiente di `theta` e
  del vento), coerentemente con il fatto che entrambi i termini valgono A6.
  Una sola tabella, due prodotti -- la stessa scelta e' replicata in Simulink
  (`build_hm3_full_ascent.m`, riga 86: `W('c6/1','P6/1'); W('c6/1','P8/1')`).
- Ogni valutazione della RHS costa **8 valutazioni di interpolante** -- le 7
  tabelle `fc1..fc7`, ma `fc6(t)` compare due volte nella riga 28 --
  **piu' 1 chiamata al vento** (e altre 2 se `sched = 1`). Con `ode45` a
  `RelTol = 1e-8`
  (`main_full_ascent.m`, riga 29) su 135 s di volo sono parecchie centinaia di
  migliaia di lookup: da qui la scelta di `griddedInterpolant` (oggetto
  precompilato) invece di `interp1` dentro il loop.

---

## Che cosa viene interpolato davvero (init_simulink_lpv.m, righe 49-107)

Questo e' il punto sottile che vale la pena saper difendere.

```matlab
c1   = a1;                 % * zdot
c2   = a1.*V + a4;         % * theta
c5   = A6./Vsafe;          % * zdot
gi = @(y) griddedInterpolant(tg, y, 'linear', 'nearest');
```

- Righe 50-52: la griglia `tg` e' quella **del dataset** (`L.V.Time`, troncata
  a `Tstop`); tutti i segnali vengono riportati su `tg` con `interp1` (riga
  52), interpolazione lineare di default. Il codice non specifica il passo
  della griglia: e' quello del `.mat` del professore.
- Righe 61-68: si costruiscono i **coefficienti effettivi gia' combinati**.
  `c2` non e' `a1`, `V`, `a4` separati: e' il **prodotto calcolato sulla
  griglia** e poi tabulato come una singola funzione del tempo. Idem `c5 =
  A6/Vsafe` e `c4 = a1*Vsafe`.
- Righe 101-107: gli interpolanti sono `griddedInterpolant(tg, y, 'linear',
  'nearest')` -- **lineare dentro** la griglia, **nearest fuori**. La
  "nearest" come metodo di *extrapolation* significa: fuori dai breakpoint il
  valore viene **agganciato all'estremo** (clamp). Non c'e' errore, non c'e'
  NaN, non c'e' warning: se per sbaglio si integrasse a t = 500 s si
  otterrebbero i coefficienti di t = 140 s, silenziosamente. E' una guardia
  implicita, non una guardia esplicita -- va detto.
- Riga 59: `Vsafe = max(V, 1)` protegge `A6/V` e `1/V` al lift-off, dove
  V(0) = 0. Il commento e' esplicito. Nota pero' che la protezione e'
  **asimmetrica**: `c4 = a1.*Vsafe` (riga 64) usa Vsafe, mentre nella
  formulazione teorica il termine e' `a1*V`. La giustificazione nel commento
  e' che il fattore si semplifica con `alpha_w = v_w/Vsafe` -- ma la
  semplificazione e' solo **approssimata**, perche' il Vsafe del coefficiente
  vive sulla griglia `tg` mentre quello dentro `alpha_w` e' stato valutato
  sulla griglia del vento `wg.t` (riga 89). Ai fini pratici (si parte da
  t0 = 5 s) e' irrilevante; concettualmente e' un piccolo hack.

### Interpolare i coefficienti grezzi o quelli combinati non e' la stessa cosa

L'interpolazione lineare **non commuta con la moltiplicazione**. Dati due
breakpoint consecutivi t_k, t_{k+1} e s in [0,1]:

    interp(a1*V)(s) = (1-s)*a1_k*V_k + s*a1_{k+1}*V_{k+1}
    interp(a1)(s) * interp(V)(s)
        = [(1-s)*a1_k + s*a1_{k+1}] * [(1-s)*V_k + s*V_{k+1}]

I due risultati coincidono solo nei nodi (s = 0 e s = 1); in mezzo
differiscono di un termine `s*(1-s)*(a1_{k+1}-a1_k)*(V_{k+1}-V_k)`, cioe' di
un O(dt^2) proporzionale alle variazioni dei due fattori. Conseguenza
operativa, ed e' un fatto verificabile nel repo:

- `ode_lpv_ascent.m` (questo file) integra le **combinazioni tabulate**
  (`fc2` = interpolante di `a1*V + a4`);
- `ode_lpv_flex.m` interpola i **coefficienti grezzi** (`fa1`, `fV`, `fa4`,
  righe 19-20) e forma `a1*V + a4` **dentro** la RHS (riga 38).

Quindi i due file **non integrano esattamente lo stesso sistema LTV rigido**,
pur partendo dallo stesso dataset. La differenza e' piccola (i coefficienti
sono lisci sulla griglia del dataset) e nessuno dei due e' "sbagliato" --
sono due discretizzazioni diverse dello stesso LPV continuo. Ma va saputo:
se domani si confrontassero numero per numero le uscite dei due file, un
disaccordo dell'ordine di O(dt^2) e' atteso, non e' un bug.

Una terza variante esiste in Simulink: `build_hm3_full_ascent.m` (righe 56-64)
costruisce `alpha_w = (v_wp + turb) * invV(t)` moltiplicando il vento
*istantaneo* del generatore per la **lookup** `invV`, mentre la RHS interpola
`alpha_w` gia' diviso. Anche qui: stessa fisica, discretizzazioni diverse.
Il fatto che i due si sovrappongano a ~1e-7 rad (README) e' un risultato
numerico, non un'identita' algebrica.

> **Possibile domanda d'esame** -- e' meglio tabulare `a1*V + a4` o
> interpolare a1, V, a4 separatamente?
> *Risposta:* Dipende da che cosa e' liscio. Tabulare il coefficiente
> effettivo garantisce che il plant congelato ricostruito a un nodo sia
> esattamente quello del dataset e toglie dal loop caldo l'aritmetica di
> ricombinazione (un prodotto per termine, niente `a1*V + a4` da riformare a
> ogni step). Non fa pero' risparmiare lookup: qui servono 7 tabelle
> (`fc1..fc7`), mentre `ode_lpv_flex` se la cava con 6 grezze (`fa1, fa3, fa4,
> fA6, fK1, fV`, riga 19). Interpolare i grezzi e' piu' flessibile (posso ricostruire
> qualunque combinazione, per esempio A6/V e a1*V con lo stesso V) e mantiene
> il significato fisico dei singoli termini. Le due scelte differiscono solo
> a O(dt^2) fra i nodi; l'importante e' non mescolarle inconsapevolmente.

---

## La fallacia del frozen-time (il punto centrale)

Tutto l'impianto di HM3 -- Nichols, margine di guadagno 6 dB, margine di fase
30 gradi, notch -- vive su plant **congelati**. La tentazione e' concludere:
"ho verificato che ogni LTI congelato A(t_i) e' stabile per ogni t_i, dunque
il sistema tempo-variante e' stabile". **E' falso**, ed e' un classico
controesempio della teoria dei sistemi LTV.

Il fatto matematico: per un sistema `xdot = A(t)*x`, gli autovalori di A(t)
**non governano** la stabilita'. Esistono matrici A(t) i cui autovalori
congelati hanno parte reale negativa **costante** per ogni t, e la cui
soluzione **diverge** esponenzialmente. L'esempio standard nei testi (Khalil,
*Nonlinear Systems*) e' una A(t) periodica di rotazione+scala i cui autovalori
sono costanti in `(-1 +/- i*sqrt(7))/4` (parte reale -0.25 per ogni t), ma la
cui matrice di transizione contiene un termine `exp(+0.5*t)`. Intuizione: la
stabilita' LTV dipende dalla **matrice di transizione** Phi(t, t0), non dallo
spettro istantaneo; una A(t) che ruota velocemente i propri autovettori puo'
"pompare" energia da una direzione contraente a una espansiva prima che la
contrazione abbia il tempo di agire. Vale anche il viceversa: un sistema i cui
LTI congelati sono **tutti instabili** puo' essere LTV-stabile.

Cosa **salva** l'ingegneria: il teorema di variazione lenta. Se A(t) e'
stabile per ogni t con margini uniformi e `||Adot(t)||` e' sufficientemente
piccola, allora l'LTV e' esponenzialmente stabile. E' esattamente l'ipotesi
implicita del gain scheduling: "il parametro varia lentamente rispetto alla
dinamica dell'anello". Qui la dinamica di anello ha banda dell'ordine di
qualche rad/s (costante di tempo ~1 s) mentre A6(t), K1(t), V(t) cambiano su
scale di decine di secondi: la separazione c'e', ma e' **assunta, non
dimostrata**. Il codice non calcola da nessuna parte una stima di `||Adot||`
ne' una funzione di Lyapunov parameter-dependent, e non risolve nessuna LMI.

**Come il codice affronta la cosa.** In modo pragmatico e onesto:

1. `main_full_ascent.m` righe 40-47 fa lo **sweep dei margini congelati**
   lungo la traiettoria -- che e' *necessario* ma non *sufficiente*;
2. righe 31-32 **propaga davvero l'LTV** con `ode45` sullo stesso anello, con
   il vento vero del generatore, e guarda a posteriori se `theta` e `z`
   restano limitati;
3. le figure mostrano entrambe le cose, e il README dichiara che "la
   integrazione LTV ode45 e' la source of truth".

Il limite residuo, da dire all'orale prima che lo dica l'esaminatore:
**anche il punto 2 non e' una dimostrazione di stabilita'**. E' *una*
traiettoria, con *una* realizzazione di vento, da *una* condizione iniziale
(`x0 = zeros(4,1)`, riga 28 di `main_full_ascent.m`). Un sistema LTV lineare
non ha nemmeno bisogno di ingressi per divergere: se si volesse una prova, si
dovrebbe integrare la **matrice di transizione** (4 run con le colonne
dell'identita' come condizione iniziale e ingressi nulli) e verificare che
`||Phi(t,t0)||` decada, oppure cercare una `P(t) > 0` con
`Pdot + A'P + PA < 0`. Nessuna delle due cose e' nel repo.

> **Possibile domanda d'esame** -- se i margini congelati non bastano, perche'
> li calcoliamo?
> *Risposta:* Perche' sono una condizione necessaria (un LTI congelato
> instabile e' quasi sempre una bandiera rossa) e perche' sono la sola cosa su
> cui si sa *progettare* con metodi classici. Sono uno strumento di sintesi;
> la verifica va fatta sul sistema vero. Il codice rispetta questa divisione
> del lavoro: progetta congelato, verifica integrando l'LTV.

---

## Limiti e punti fragili (riepilogo onesto)

- **Attuatore ideale.** Niente TVC, ritardo, notch, saturazione: `delta` e'
  applicato istantaneamente. Confrontare i margini di questo anello con quelli
  del Task 2 di HM3 e' quindi un confronto fra cose diverse.
- **Nessuna guardia esplicita fuori griglia.** Gli interpolanti fanno clamp
  silenzioso ('nearest'). Idem `fKp`/`fKd`, definiti da `t0 = 5 s`: sotto i
  5 s si usano i guadagni del primo nodo (la lookup Simulink usa
  `ExtrapMethod = 'Clip'`, righe 106-109 di `build_hm3_full_ascent.m`, proprio
  per replicare questo comportamento).
- **Coerenza del segno di alpha_w**: la RHS usa `-alpha_w`, l'indicatore di
  carico `q_bar*alpha` in `main_full_ascent.m` (riga 153) usa `+alpha_w`.
- **I margini vengono letti da `assemble_loop`** (`HM3/assemble_loop.m`,
  riga 38: `getLoopTransfer(T, 'delta', -1)`), cioe' su un anello SISO rotto
  su `delta` che include **anche** la retroazione di drift `z, zdot`. E' la
  stessa convenzione di tutto HM3: qualunque critica a quella scelta si eredita
  qui identica.
- **La linearizzazione su 140 s.** `z, zdot` sono perturbazioni normali alla
  traiettoria di riferimento: su tutto l'arco dell'ascesa l'ipotesi si
  indebolisce (lo dice il README stesso). `z` va letto come indicatore
  relativo, non come cross-range fisico.
- Il vento non forza la deriva in modo indipendente: entra solo attraverso
  `alpha_w`, cioe' come angolo, non come forza laterale addizionale.

---

## Possibili domande d'esame

**D: Che cosa rende questo sistema LPV, e perche' non basta il modello
congelato di `build_plant_rigid`?**
R: I coefficienti A6 (momento aerodinamico destabilizzante), K1 (efficacia
TVC), V, a1, a3, a4 cambiano lungo l'ascesa perche' cambiano pressione
dinamica, massa, inerzia e spinta. Il modello congelato prende una loro
istantanea a t = 72 s e produce un LTI su cui si puo' progettare con metodi in
frequenza. Ma la dinamica vera e' `xdot = A(t)*x + B(t)*u`: per propagarla
bisogna interpolare i coefficienti **dentro** la RHS a ogni valutazione, ed e'
quello che fanno le righe 25-28. Il progetto rimane congelato, la verifica no.

**D: La stabilita' di tutti i plant congelati garantisce la stabilita' del
sistema tempo-variante?**
R: No, ed e' l'errore concettuale piu' comune. Per `xdot = A(t)*x` la
stabilita' dipende dalla matrice di transizione, non dallo spettro istantaneo:
esistono controesempi classici (Khalil) con autovalori congelati costanti a
parte reale -0.25 e soluzioni che divergono come exp(+0.5t). Vale un teorema
di variazione lenta -- se `||Adot||` e' abbastanza piccola la conclusione
diventa valida -- ma "abbastanza piccola" e' un'ipotesi che qui viene assunta,
non verificata. Il codice mitiga il rischio propagando davvero l'LTV con
`ode45` (righe 31-32 di `main_full_ascent.m`) invece di fidarsi dei soli
margini congelati.

**D: E allora la simulazione ode45 dimostra la stabilita'?**
R: Neanche. Dimostra che *quella* traiettoria, con *quel* vento e *quella*
condizione iniziale, resta limitata su 140 s. Una prova richiederebbe di
propagare la matrice di transizione (o esibire una funzione di Lyapunov
tempo-variante / risolvere una LMI parameter-dependent). Il repo non lo fa e
non lo pretende: e' evidenza numerica forte, non un certificato.

**D: Perche' l'interpolazione dei coefficienti "effettivi" (c2 = a1*V + a4) non
e' equivalente a interpolare a1, V, a4 e moltiplicare dentro la RHS?**
R: Perche' l'interpolazione lineare non commuta con il prodotto: fra due
breakpoint la differenza vale `s*(1-s)*(delta_a1)*(delta_V)`, che si annulla
solo nei nodi. `ode_lpv_ascent` interpola i coefficienti combinati (fc1..fc7),
`ode_lpv_flex` interpola i grezzi e li combina in linea: i due file integrano
quindi due sistemi LTV leggermente diversi, entrambi consistenti con il
dataset ai nodi. E' una differenza O(dt^2), non un bug, ma va dichiarata se si
confrontano i risultati.

**D: Perche' non c'e' un blocco `arguments` in questa funzione?**
R: E' un hot loop dentro `ode45`, chiamato ~1e5 volte per simulazione. La
validazione degli input vive nel confine del run (`init_simulink_lpv`, che ha
il suo `arguments`) e nei main. E' una scelta esplicita, dichiarata nelle
righe 11-12 e nella convenzione del repo, non una dimenticanza.

**D: Come entra il vento e come e' chiuso l'anello dentro la RHS?**
R: Il vento e' un `griddedInterpolant` (`M.windfun`) costruito una volta sola
sull'uscita del generatore Simulink del professore, gia' convertito in angolo
`alpha_w = v_w/V`; viene valutato all'istante corrente alla riga 14 ed entra
con segno meno in entrambi i canali, perche' l'incidenza aerodinamica vera e'
`theta + zdot/V - alpha_w` (conta la velocita' relativa all'aria). L'anello e'
chiuso in linea alla riga 22: `delta = -(Kp*theta + Kd*thetadot + Kp_z*z +
Kd_z*zdot)`, attuatore ideale, `theta_ref = 0`. Non c'e' un blocco controllore
separato: plant, vento e controllore stanno tutti dentro la stessa RHS, ed e'
per questo che il file e' la reference contro cui si valida il modello
Simulink.
