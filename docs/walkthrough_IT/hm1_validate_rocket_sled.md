# HM1/validate_rocket_sled.m

## Ruolo del file nel progetto

Questo e' lo **script di validazione della macchina di shooting** di HM1. Non
risolve nessun task dell'homework: risolve un problema-giocattolo diverso -- il
*rocket sled* (slitta a razzo) a energia minima -- che ha una **soluzione in forma
chiusa nota**, e verifica che la catena `ode45 + fsolve` recuperi esattamente
quella soluzione. L'intestazione lo dichiara: "Checks the ode45 + fsolve
single-shooting machinery (same tolerances as the ascent BVP) against the
analytic optimum before using it on the ascent problem" (righe 1-4).

La logica e' quella dello **strumento tarato prima dell'uso**. Il problema di
ascesa indiretta di HM1 non ha soluzione analitica: quando `fsolve` non converge,
o converge a qualcosa di strano, non si sa se la colpa e' (a) della derivazione
delle condizioni di Pontryagin, (b) di un errore di segno nelle equazioni dei
costati, (c) delle tolleranze di `ode45`, (d) del guess iniziale, oppure (e) del
fatto che il problema stesso e' mal posto. Il sled toglie di mezzo tutti gli
imputati tranne quelli implementativi: se il solver **non** recupera i costati
analitici `lam_r0 = lam_v0 = 3/2` su un problema di cui si conosce la risposta,
allora il bug e' nel codice di shooting, non nella fisica dell'ascesa.

Il file e' **autonomo**: non chiama `ode_burn.m` ne' nessun altro file di HM1.
Contiene due funzioni locali (`sled_ode`, `sled_residual`) che sono repliche in
miniatura della struttura usata nei `main_task*.m` -- stato + costati integrati
insieme, residuo terminale valutato in `tf`, `fsolve` sulle condizioni iniziali
incognite dei costati. Riusa **le stesse tolleranze di `ode45`** del problema vero
(`RelTol = 1e-10`, `AbsTol = 1e-12`) e su `fsolve` ne adotta di **piu' strette**
(`1e-12` su funzione e passo, contro `1e-10` dei `main_task*.m`), quindi valida
anche la scelta delle tolleranze, non solo la logica.

Il commento in riga 1 lo etichetta come "Appendix A": e' il materiale di supporto
del report, non parte del flusso di soluzione.

---

## Header e dati del problema (righe 1-22)

```matlab
%  Problem:  rdot = v,  vdot = u,   min J = int_0^2 u^2 dt
%            r(0)=0, v(0)=0, r(2)=1/2, v(2)=0,  tf = 2 fixed.
%  Hamiltonian:  H = -u^2 + lam_r*v + lam_v*u
%  Costates:  lam_r_dot = 0,  lam_v_dot = -lam_r ;  u = lam_v/2.
%  Closed form: lam_r0 = lam_v0 = 3/2  =>  u*(t) = (3/4)(1 - t).
```

- **Righe 6-7 -- il problema.** Doppio integratore, controllo `u` (accelerazione)
  libero e **non vincolato**, costo quadratico di energia. Si parte da fermo
  nell'origine e si deve arrivare a `r = 1/2` con velocita' nulla in `tf = 2`,
  che e' **fissato**. E' il problema LQ piu' semplice che esista con una struttura
  di shooting non banale: due incognite (i due costati iniziali), due condizioni
  terminali.

- **Riga 8 -- l'Hamiltoniana.** Scritta come

      H = -u^2 + lam_r*v + lam_v*u

  Attenzione al **segno del termine di costo**: e' `-u^2`, non `+u^2`. Il codice
  adotta la convenzione di **massimizzazione** di H (PMP nella forma "maximum
  principle"): massimizzare l'integrale di `-u^2` equivale a minimizzare
  l'integrale di `u^2`. Questa e' la stessa convenzione usata in `ode_burn.m`,
  dove l'angolo di spinta viene scelto come `phi = atan2(lam_vy, lam_vx)`, cioe'
  allineando la spinta al primer vector per **massimizzare** `lam_v . u_hat`. La
  coerenza fra i due file non e' un dettaglio: e' proprio cio' che il sled sta
  validando.

- **Riga 9 -- legge di controllo e costati.** Da `dH/du = -2u + lam_v = 0` segue

      u = lam_v/2

  che e' un controllo **interno** (nessuna saturazione, perche' `u` non e'
  vincolato). Le equazioni dei costati vengono da `lam_dot = -dH/dx`:

      lam_r_dot = -dH/dr = 0          (r non compare in H)
      lam_v_dot = -dH/dv = -lam_r

  Quindi `lam_r` e' costante e `lam_v` e' **lineare nel tempo**. Notare
  l'analogia strutturale con HM1: anche li' `lam_vx` e' costante e `lam_vy` e'
  lineare (`lam_vy = lam_vy0 - lam_y*t`), il che produce la *linear tangent law*.

- **Riga 10 -- la soluzione chiusa.** Derivazione. Poniamo `lam_r = a` costante e
  `lam_v = b - a*t`. Allora `u = (b - a*t)/2` e integrando da fermo:

      v(t) = (b*t - a*t^2/2)/2
      r(t) = (b*t^2/2 - a*t^3/6)/2

  Imponendo le due condizioni terminali con `tf = 2`:

      v(2) = (2b - 2a)/2 = b - a = 0        ->  b = a
      r(2) = (2b - (8/6)a)/2 = b - (2/3)a = 1/2

  Sostituendo `b = a`: `a/3 = 1/2`, cioe' `a = 3/2`. Da cui
  `lam_r0 = lam_v0 = 3/2` e

      u*(t) = (3/2 - (3/2)*t)/2 = (3/4)*(1 - t)

  Il controllo ottimo e' una **rampa lineare** che parte da `0.75`, si annulla a
  `t = 1` e diventa negativa (frenata) nella seconda meta': accelera, poi frena
  per arrivare a velocita' nulla. Fisicamente sensato, e questo e' il punto:
  la risposta e' verificabile *a mano*.

- **Riga 12 -- `clear; close all; clc;`.** Questo e' uno **script**, non una
  funzione e non una classe di test. Conseguenza pratica: non puo' essere
  raccolto da `runtests`, e non lancia un errore in caso di fallimento (vedi
  righe 44-48). Va eseguito a mano dalla cartella `HM1/`.

- **Righe 15-17 -- dati.** `tf = 2`, `rf = 1/2`, `vf = 0`. Nessuna
  nondimensionalizzazione: il sled e' gia' adimensionale per costruzione.

- **Righe 19-22 -- tolleranze.** Sono deliberatamente **le stesse** (anzi, su
  `fsolve` piu' strette) di quelle usate nei `main_task*.m`:
  `ode45` con `RelTol = 1e-10, AbsTol = 1e-12`; `fsolve` con
  `FunctionTolerance = StepTolerance = 1e-12`, `Display` spento,
  `MaxIterations = 500`, `MaxFunctionEvaluations = 5000`. Se il sled passasse a
  tolleranze lasche ma non a quelle vere, la validazione non direbbe nulla di
  utile. Il fatto che siano esplicitate (non lasciate ai default) segue la
  convenzione di riproducibilita' del repo.

> **Possibile domanda d'esame** -- Perche' l'Hamiltoniana ha il segno `-u^2` e non
> `+u^2`?
> *Risposta:* Perche' il codice usa la forma del Principio del Massimo: il
> controllo ottimo **massimizza** H. Con `H = -u^2 + lam.f`, massimizzare
> l'integrale di `-u^2` e' identico a minimizzare l'energia. Se scrivessi
> `H = +u^2 + lam.f` dovrei minimizzare H, e la condizione `dH/du = 0` darebbe
> `u = -lam_v/2`, con un segno invertito nei costati. Le due formulazioni sono
> equivalenti purche' si resti coerenti; l'importante e' che il sled usi **la
> stessa convenzione** di `ode_burn.m`, altrimenti non validerebbe nulla.

---

## Single shooting: `fsolve` sui costati iniziali (righe 24-27)

```matlab
lam_guess = [0; 0];   % away from the analytic (3/2, 3/2)
[lam_sol, res, ef] = fsolve(@(L) sled_residual(L, tf, rf, vf, ...
                            opts_ode), lam_guess, opts_fs);
```

- **Riga 25 -- il guess iniziale.** `[0; 0]`, con il commento esplicito "away from
  the analytic (3/2, 3/2)". La scelta e' voluta: partire dal valore corretto
  renderebbe il test vuoto (`fsolve` restituirebbe subito il punto di partenza).

  **Onesta' richiesta**: questo test e' meno severo di quanto sembri. Le
  dinamiche del sled -- stato *e* costati -- sono **lineari**, quindi la mappa di
  flusso e' lineare e il residuo e' una funzione **affine** delle incognite `L`.
  Un metodo di Newton (e `fsolve` in modalita' `trust-region-dogleg` lo e')
  risolve un sistema affine **in un solo passo**, esattamente, qualunque sia il
  guess. Quindi il sled **non** mette alla prova la robustezza di `fsolve` sulla
  non linearita' (che e' esattamente cio' che rende difficile il problema di
  ascesa e che obbliga alla *continuation*). Cio' che valida davvero e':
  l'accuratezza di `ode45`, l'assemblaggio del residuo, le convenzioni di segno
  dei costati, il cablaggio della function handle. Sono i bug piu' comuni, ma non
  sono tutti i bug possibili.

- **Righe 26-27 -- la chiamata.** Le uscite sono `[x, fval, exitflag]` nella
  firma di `fsolve`, quindi `res` e' il **vettore residuo valutato nella
  soluzione** (non il "residuo" nel senso di norma) e `ef` e' l'exitflag
  (`> 0` = convergenza). Le unica incognite sono i due costati iniziali: lo stato
  iniziale e' noto e fissato dentro `sled_residual`.

---

## Report e cross-check sul controllo (righe 29-48)

```matlab
[tt, S] = ode45(@(t,s) sled_ode(t,s), [0 tf], ...
                [0;0;lam_r0;lam_v0], opts_ode);
u_num = S(:,4)/2;                 % u = lam_v/2
u_ana = 0.75*(1 - tt);
```

- **Righe 31-36 -- cosa viene stampato.** Quattro diagnostiche: l'exitflag, i due
  costati recuperati confrontati con `1.5`, la norma del residuo terminale, e
  l'errore sui costati `||lam - [1.5; 1.5]||`. Notare la riga 36:
  `norm(lam_sol - 1.5)` sfrutta l'espansione scalare di MATLAB (sottrae `1.5` a
  entrambe le componenti) -- corretto, ma implicito.

- **Righe 38-42 -- il secondo livello di verifica.** Non basta che i costati
  iniziali siano giusti: si **reintegra** la traiettoria con i costati trovati e
  si confronta il **profilo di controllo intero** `u_num(t) = lam_v(t)/2` contro
  la rampa analitica `u*(t) = 0.75*(1 - t)`, stampando `max |u_num - u*|`.

  Questo e' un controllo piu' forte del solo confronto sui costati iniziali:
  verifica che la **propagazione temporale** dei costati sia corretta lungo tutto
  l'arco, non solo il punto di partenza. Un errore di segno in `lam_v_dot`, per
  esempio, potrebbe (in linea di principio) essere compensato da un valore
  iniziale sbagliato in modo da soddisfare comunque le condizioni terminali; il
  confronto sul profilo `u(t)` lo scoprirebbe.

- **Righe 44-48 -- il criterio di PASS.** Due condizioni **in AND**:
  `ef > 0` (fsolve dichiara convergenza) **e** `norm(lam_sol - 1.5) < 1e-6`.

  La tolleranza di accettazione e' quindi `1e-6` **sui costati**, molto piu' larga
  delle tolleranze numeriche usate (`1e-10 / 1e-12`): e' una soglia di
  *correttezza strutturale*, non di precisione numerica. In pratica il codice, con
  quelle tolleranze, recupera i costati a molte cifre in piu'.

  **Limite da dichiarare**: il fallimento produce solo una `fprintf` con la
  scritta `FAIL`. Lo script **non lancia un errore**, non ritorna un exit code,
  non e' una `matlab.unittest.TestCase`. Quindi non entra in
  `runtests('HM1/tests')` e non puo' rompere una pipeline di CI: se qualcuno
  rompesse il solver, questo file lo direbbe soltanto a chi lo guarda. Trasformarlo
  in una classe di test con una `verifyEqual(lam_sol, [1.5;1.5], 'AbsTol', 1e-6)`
  sarebbe il naturale passo successivo -- il resto della repo (`HM1/tests/`) segue
  gia' quella convenzione.

  **Nota**: il `max |u_num - u*|` stampato alla riga 42 **non entra** nel criterio
  di PASS/FAIL. E' solo informativo.

> **Possibile domanda d'esame** -- Se il solver recupera i costati iniziali
> corretti, perche' serve anche il confronto sul profilo `u(t)`?
> *Risposta:* Perche' le condizioni terminali sono solo due numeri: piu' errori
> possono in linea di principio cancellarsi e produrre comunque un residuo nullo
> in `tf`. Confrontare l'intera storia del controllo verifica anche la dinamica
> dei costati **lungo** l'arco (`lam_v_dot = -lam_r`), cioe' la propagazione, non
> solo il punto iniziale. E' un test sulla funzione, non su due scalari.

---

## `sled_ode` (righe 51-61)

```matlab
function ds = sled_ode(~, s)
    ds = [ s(2);       % rdot = v
           s(4)/2;     % vdot = u = lam_v/2
           0;          % lam_r_dot = 0
          -s(3) ];     % lam_v_dot = -lam_r
end
```

- **Riga 51 -- firma.** Il tempo e' scartato con `~`: il sistema aumentato
  (stato + costati) e' **autonomo**, il tempo non compare esplicitamente. E'
  esattamente la stessa struttura di `ode_burn.m`, con una differenza importante:
  li' i costati `lam_vx, lam_vy` **non** sono integrati, sono ricostruiti in forma
  chiusa da `p` e `t` (per questo `ode_burn` usa `t` esplicitamente e integra solo
  6 stati). Qui invece i due costati sono **stati veri** dell'integrazione, quindi
  il vettore e' `[r; v; lam_r; lam_v]`, 4 componenti.

- **Righe 57-60 -- il sistema aumentato.** La sostituzione della legge di controllo
  `u = lam_v/2` **dentro** la dinamica (riga 58) e' il cuore del metodo indiretto:
  il controllo sparisce come variabile indipendente e sopravvive solo come
  funzione dei costati. Quel che resta e' un sistema di ODE puro, integrabile con
  `ode45`, la cui unica indeterminazione sono le **condizioni iniziali dei
  costati** -- che e' esattamente cio' su cui si fa shooting.

- Notare che il sistema e' **lineare** in `s` (matrice costante 4x4). Da qui la
  conseguenza detta sopra: il residuo di shooting e' affine.

---

## `sled_residual` (righe 63-74)

```matlab
function res = sled_residual(L, tf, rf, vf, opts_ode)
    [~, S] = ode45(@(t,s) sled_ode(t,s), [0 tf], ...
                   [0; 0; L(1); L(2)], opts_ode);
    res = [S(end,1) - rf; S(end,2) - vf];
end
```

- **Riga 63 -- firma.** Prende le due incognite `L = [lam_r0; lam_v0]` e tutti i
  dati del problema passati per valore (nessuna variabile globale, nessuna
  dipendenza dal workspace dello script): la funzione e' pura, e quindi il
  residuo e' riproducibile.

- **Riga 72 -- integrazione forward.** Lo stato iniziale e' `[0; 0; L(1); L(2)]`:
  posizione e velocita' **note** (partenza da fermo nell'origine), costati
  **incogniti**. Questa e' la definizione stessa di *single shooting*: si
  trasforma un problema ai limiti (BVP) in una sequenza di problemi ai valori
  iniziali (IVP), e si cerca il valore iniziale che azzecca le condizioni al
  bordo opposto.

- **Riga 73 -- il residuo.** Due componenti, `r(tf) - rf` e `v(tf) - vf`: tante
  quante le incognite. Sistema quadrato 2x2, come deve essere perche' `fsolve`
  abbia senso.

  **Il conto delle incognite e' il punto pedagogico**. Nel sled: 2 incognite
  (costati iniziali), 2 residui (condizioni terminali), `tf` **fissato**. Nel
  problema vero di HM1 (Task 1): 4 incognite
  (`lam_vx0, lam_vy0, lam_y, tf`) e 4 residui (`y(tf) - yf`, `vx(tf) - 1`,
  `vy(tf)`, e `H(0) = 0`), perche' il tempo finale e' **libero** e questo aggiunge
  una condizione di trasversalita' (`H = 0`) e un'incognita.

  **Limite di copertura da dichiarare apertamente**: il sled **non** valida ne'
  la condizione `H = 0` a tempo libero, ne' la normalizzazione `lam_m0 = 1`, ne'
  la gestione dell'arco di coast di Task 3 (funzione di switching), ne' lo
  staging di Task 4. Valida il **nucleo** -- integrazione + residuo + Newton -- non
  gli strati sopra. Va detto: e' una validazione della *plumbing*, non della
  formulazione completa.

> **Possibile domanda d'esame** -- Perche' il sled ha 2 incognite e HM1 ne ha 4,
> se in entrambi lo stato e' completamente noto all'istante iniziale?
> *Risposta:* Perche' in HM1 il tempo finale e' **libero**. Nel sled `tf = 2` e'
> dato, quindi le incognite sono solo i costati iniziali che non hanno condizione
> iniziale (`lam_r0, lam_v0`). In HM1 `tf` e' esso stesso un'incognita, e la
> corrispondente condizione di trasversalita' e' l'annullamento
> dell'Hamiltoniana; inoltre HM1 ha una condizione terminale in piu' (`y(tf)`
> oltre alle due velocita'). Il conteggio torna: 3 costati incogniti + `tf` = 4
> incognite, 3 condizioni terminali + `H = 0` = 4 residui.

---

## Possibili domande d'esame

**D: Perche' validare un solver indiretto contro un problema con soluzione in
forma chiusa, invece che contro un metodo diretto o un altro codice?**
R: Perche' una soluzione chiusa e' un **oracolo esatto**: non ha errore
numerico, non ha ipotesi implicite, e non puo' essere sbagliata "allo stesso
modo" del codice che sto testando. Confrontare due solver numerici mi dice solo
che sono d'accordo, non che hanno ragione (potrebbero condividere lo stesso
errore di segno nei costati). Il sled invece mi da' il valore vero
`lam_r0 = lam_v0 = 3/2` derivabile a mano in cinque righe: se il codice non lo
riproduce, il codice ha torto, punto. In piu', un oracolo esatto consente di
misurare l'**errore assoluto** (`max |u_num - u*|`), non solo la coerenza interna.

**D: Che cos'e' esattamente il "single shooting" e cosa lo distingue dal
collocation usato in HM2?**
R: Lo shooting converte il BVP in un problema di ricerca di radici sulle sole
condizioni iniziali incognite: si integra **in avanti** l'intero arco con `ode45`
e si annulla il residuo al bordo finale con `fsolve`. Le incognite sono poche
(2 nel sled, 4 in HM1) ma il residuo e' **estremamente non lineare** nelle
incognite, perche' passa attraverso l'intera propagazione -- piccole variazioni
dei costati iniziali possono far esplodere lo stato finale. Il collocation
(HM2) discretizza invece l'intera traiettoria e tratta *tutti* i nodi come
incognite: molte piu' variabili, ma ciascun vincolo e' locale e il problema e'
molto meglio condizionato. Lo shooting compra accuratezza (integratore adattivo,
`RelTol = 1e-10`) al prezzo del condizionamento; il collocation fa il contrario.

**D: Il guess iniziale e' `[0; 0]`, lontano dalla soluzione. Vuol dire che il
solver e' robusto?**
R: **No**, e va detto onestamente. Il sistema aumentato del sled e' lineare,
quindi il residuo e' una funzione affine delle incognite e Newton lo risolve
esattamente in un passo qualunque sia il punto di partenza. La convergenza da
`[0;0]` non dimostra robustezza rispetto alla non linearita'. E' proprio perche'
il problema di ascesa **non** ha questa proprieta' che HM1 deve ricorrere alla
*continuation* (partire da `Q` con T/W ~ 1.8 e fare warm start sui vicini). Il
sled valida la correttezza, non la robustezza.

**D: Quali parti della macchina di HM1 il sled NON valida?**
R: Diverse, e sono le piu' delicate: (1) la condizione di trasversalita'
`H(tf) = 0` per tempo finale libero (nel sled `tf` e' fissato); (2) la
normalizzazione `lam_m0 = 1` che sfrutta l'omogeneita' di grado 1 di `H` nei
costati; (3) la valutazione di `H = 0` all'istante **iniziale** invece che
finale, dove diventa algebrica (`H(0) = -lam_vy0 + T*(|lam_v0| - 1/c)`); (4) la
funzione di switching e la giunzione dell'arco di coast di Task 3; (5) lo staging
di Task 4. Il sled valida il nucleo `ode45 + residuo + fsolve` e le convenzioni
di segno di Pontryagin. Gli strati sopra restano coperti solo dai test in
`HM1/tests/` e dalla plausibilita' fisica dei risultati.

**D: La tolleranza di PASS e' `1e-6` sui costati mentre l'ODE gira a `1e-10`.
Non e' incoerente?**
R: No, sono due cose diverse. `1e-10 / 1e-12` sono le tolleranze **numeriche**
dell'integratore e del root finder: definiscono quanto accuratamente il solver
insegue la propria soluzione. `1e-6` e' invece la soglia di **accettazione del
test**: distingue "il codice implementa il problema giusto" da "il codice
implementa un problema diverso". Un errore di segno o una formula sbagliata
produrrebbe un errore di ordine unitario, non `1e-7`. La soglia larga rende il
test robusto rispetto a differenze di piattaforma/versione MATLAB senza perdere
potere diagnostico.

**D: Come miglioreresti questo file?**
R: Tre cose concrete. (1) Convertirlo in una `matlab.unittest.TestCase` in
`HM1/tests/`, cosi' che `runtests` lo esegua e un fallimento diventi un errore
vero invece che una `fprintf`. (2) Includere `max |u_num - u*|` nel criterio di
PASS, oggi solo stampato. (3) Aggiungere un secondo benchmark **a tempo libero**
(per esempio il classico min-time double integrator con controllo limitato, che
ha soluzione bang-bang nota) per validare anche la condizione `H = 0` e la
struttura di switching, che il sled a energia minima e tempo fissato non tocca.
