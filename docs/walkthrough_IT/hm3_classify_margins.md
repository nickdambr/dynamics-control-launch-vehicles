# HM3/classify_margins.m

## Ruolo del file nel progetto

Questo file e' la **risposta di HM3 alla domanda "quanto margine ho?"** su un
lanciatore che e' open-loop instabile. E' chiamato da `design_controller`
(riga 66 e riga 84, dentro la funzione di costo), da `main_task2.m` (righe 125 e
145) e da `main_task3.m` (riga 56). Prende il loop aperto pieno
L = K_FCS * F * C * G (PD + attuatore TVC + ritardo + notch + impianto con deriva
e bending), chiama `allmargin` e **smista ogni attraversamento nella banda di
frequenza a cui appartiene fisicamente**, restituendo una struct con margini
etichettati: Aero GM, Rigid GM, Rigid PM, Flex GM, Flex PM, |L(omega_BM)| e delay
margin.

Il file esiste perche' `margin()` -- il comando che uno userebbe d'istinto -- **su
questo loop non ha significato**, per due motivi distinti che vanno tenuti
separati.

**Primo motivo: la stabilita' condizionale.** L'aeroshell a max-qbar e'
aerodinamicamente instabile: il polo del corpo rigido sta a
+sqrt(A_6) = +1.84 rad/s (nel loop assemblato, con l'accoppiamento di deriva, si
misura +1.8165 rad/s, piu' un secondo polo instabile lento a +0.0291 rad/s e
l'integratore di posizione laterale in 0). Con poli instabili in catena aperta,
il criterio di Nyquist non e' "non circondare -1" ma "circondare -1 esattamente
P volte in senso antiorario". La conseguenza pratica e' che il loop e'
**condizionalmente stabile**: esiste una **banda** di guadagni ammissibili, e si
esce dalla stabilita' sia **abbassando** sia **alzando** il guadagno. Verificato
eseguendo il codice sul loop completo del Task 2 (TVC + ritardo + notch): il
ciclo chiuso e' stabile per k in **[0.50, 2.39]**, instabile fuori. I due estremi
sono esattamente i due margini che il classificatore etichetta come **Aero GM**
(-6.0 dB, cioe' k = 0.50) e **Rigid GM** (+7.56 dB, cioe' k = 2.39). Un solo
numero -- quello che restituisce `margin()` -- non puo' descrivere una banda a
due estremi.

**Secondo motivo: il feedback di deriva laterale.** Il canale z/delta contiene un
integratore libero (posizione laterale), quindi |L| -> Inf a DC e la curva di
Nichols "viene dall'alto": eseguendo il codice, a omega = 0.01 rad/s il loop e' a
+48 dB. Il lobo di bassa frequenza cosi' generato **attraversa gli 0 dB due
volte** (nel Task 1 a 0.161 e 0.222 rad/s) con fasi che nulla hanno a che vedere
col corpo rigido -- e a queste frequenze `allmargin` diligentemente riporta phase
margin di **-133 deg e -40 deg**. Sono artefatti del canale di deriva, non
margini. La prova che questo non e' un timore teorico: nel corner **V3** del
Task 3 (mu_alpha 1.3, mu_c 0.7), `margin(L)` restituisce **PM = -8.72 deg @ 0.416
rad/s** su un sistema che `isstable` dichiara **stabile**. Il classificatore, sullo
stesso loop, restituisce il phase margin vero: **18.0 deg @ 1.50 rad/s**.

Onesta' doverosa: nei casi **nominali** (Task 1 e Task 2 dopo il re-tune)
`margin()` per pura coincidenza numerica restituisce gli stessi valori del
classificatore (-6.00 dB e 30.0 deg). L'errore di `margin()` non e' sistematico:
morde quando i parametri si spostano, cioe' esattamente nei corner e nel Monte
Carlo, cioe' esattamente dove il numero serve. Questo e' il motivo per cui il
classificatore non e' una raffinatezza estetica ma un requisito.

La struttura per bande segue D'Antuono (Fig. 3.2, Tab. 3.1) e Trotta (Tab. 4.1),
come dichiarato nella docstring righe 4-12.

---

## Docstring e contratto (righe 1-25)

- Righe 4-12: la mappa delle bande. Vale la pena rileggerla perche' e' il
  glossario di tutto l'homework:
  - **Aero GM**: gain-**reduction** margin di bassa frequenza. Il commento a
    riga 7 lo dice esplicitamente: `gmdb < 0`. E' il **bordo inferiore** della
    banda di guadagni condizionalmente stabili.
  - **Rigid PM**: phase margin al crossover del corpo rigido (~ sqrt(A_6)).
  - **Rigid GM**: gain-**increase** margin di media frequenza, prodotta dal
    ritardo di fase di attuatore + delay. **Assente con attuatore ideale**
    (Task 1) -- il codice lo dichiara a riga 11 e il calcolo lo conferma:
    `rigidGM_dB = NaN` nel Task 1.
  - **Flex GM/PM**: margini nella banda del bending.
  - Righe 13-15: gli attraversamenti sotto `w_drift` sono **artefatti di deriva** e
    **non** vengono riportati come margini rigidi; il commento dice testualmente
    che *"taking margin()'s default instead would pick one of these"*.

## Blocco `arguments` e setup (righe 27-37)

- Riga 29: `w_drift` default **0.5** rad/s. Nota: e' un default fisso, ma **tutti
  i chiamanti lo sovrascrivono** con `0.3*sqrt(A6)` (0.5517 al nominale), che ha
  il pregio di **scalare con l'instabilita'** dell'aeroshell. Nel corner V3
  (mu_alpha = 1.3) diventa 0.629 rad/s: senza questo scaling, la crossing di
  deriva a 0.4164 rad/s resterebbe comunque sotto 0.5, ma il margine e' sottile.
- Riga 30: `w_flex` default **Inf** = "nessun modo flessibile". E' il caso Task 1,
  dove tutte le maschere `gf < w_flex` diventano "tutte le frequenze" e quelle
  `gf >= w_flex` diventano vuote -> Flex GM/PM = NaN. Elegante: la stessa
  funzione serve il modello rigido e quello flessibile senza rami `if`.
- Riga 32: `w_bending` default NaN, risolto a riga 34 in `w_flex`.
- Righe 36-37: stesso `onCleanup` di `design_controller` per silenziare
  `Control:analysis:MarginUnstable` e ripristinarlo all'uscita.

## `allmargin` e l'estrazione degli attraversamenti (righe 39-41)

```matlab
am   = allmargin(L);
gmdb = 20*log10(am.GainMargin(:));  gf = am.GMFrequency(:);
pm   = am.PhaseMargin(:);           pf = am.PMFrequency(:);
```

- Riga 39: **`allmargin` invece di `margin`**. E' la scelta che rende possibile
  tutto il resto: `allmargin` restituisce **tutti** gli attraversamenti (vettori),
  `margin` ne sceglie **uno** con una regola che ignora la fisica del problema.
  Sul loop pieno del Task 2, `allmargin` trova **5** gain margins e **3** phase
  margins.
- Riga 40: la conversione in dB. `am.GainMargin` e' un **fattore moltiplicativo**,
  non un dB: per l'Aero GM vale **0.5012**, cioe' "moltiplica il guadagno per 0.50
  e vai instabile". In dB e' -6.0. **Il segno negativo e' il modo in cui il codice
  riconosce una gain-reduction margin** (riga 44: `gmdb < 0`). E' la convenzione di
  Barrows & Orr (cap. 9).
- I margini di fase (`pm`) sono in gradi, le frequenze in rad/s.

> **Possibile domanda d'esame** -- perche' il gain margin aerodinamico e'
> *negativo* in dB? Un margine negativo non vuol dire instabile?
> *Risposta:* no. Il gain margin di MATLAB e' il fattore k per cui il loop
> diventa marginalmente stabile. Se il loop e' condizionalmente stabile con banda
> k in [0.5, 2.39], allora **due** k critici esistono: k = 0.5 (< 1, cioe' < 0 dB
> -> gain **reduction** margin) e k = 2.39 (> 1, > 0 dB -> gain **increase**
> margin). Il segno codifica la **direzione** in cui il guadagno destabilizza,
> non lo stato di stabilita' presente. Il verdetto di stabilita' e' un'altra cosa
> (`isstable`).

## Gain margins: aero, rigid, flex (righe 43-46)

```matlab
[mm.aeroGM_dB,  mm.aeroGM_w]  = pick(gmdb, gf, ...
    gf > 0 & gf < opts.w_flex & gmdb < 0, 'minf');
[mm.rigidGM_dB, mm.rigidGM_w] = pick(gmdb, gf, ...
    gf > 0 & gf < opts.w_flex & gmdb > 0, 'minf');
[mm.flexGM_dB,  mm.flexGM_w]  = pick(gmdb, gf, ...
    gf >= opts.w_flex & gf <= opts.w_flex_hi, 'near', opts.w_bending);
```

- Riga 43 (commento) e la maschera **`gf > 0`**: esclude l'attraversamento a
  **frequenza zero**. Non e' teorico: sul loop pieno del Task 2, `allmargin`
  restituisce un'entrata a `gf = 0` con `gmdb = -283 dB`. Da dove viene? Dal
  fatto che il canale di deriva ha un integratore, quindi |L(0)| = Inf e il
  "fattore per andare instabile" e' numericamente zero -> -283 dB (limite di
  precisione macchina). Senza la maschera, `aeroGM_dB` varrebbe -283 dB e tutta
  la sintonizzazione sarebbe assurda. **L'integratore di deriva e' il colpevole,
  e questa maschera e' la sua neutralizzazione.**
- Riga 44 -- **Aero GM**: fra tutti gli attraversamenti con `gmdb < 0` sotto
  `w_flex`, prende quello a **frequenza minima** (`'minf'`). Traduzione fisica:
  il bordo inferiore della banda di stabilita' condizionale e' quello piu' a
  bassa frequenza, dove la curva incrocia la fase critica scendendo dal lobo di
  deriva. Task 1: **-6.000 dB @ 0.593 rad/s**. Task 2: **-6.000 dB @ 0.542 rad/s**.
  Corner V3: **-0.91 dB @ 0.600 rad/s** -- il margine aerodinamico si e' quasi
  annullato, ed e' il risultato quantitativo piu' importante del Task 3.
- Riga 45 -- **Rigid GM**: stessa banda, ma `gmdb > 0`: il bordo **superiore**
  della banda. Task 1: **NaN** (attuatore ideale, la fase non riattraversa mai il
  punto critico -- verificato: il loop resta stabile fino a k = 1e6). Task 2:
  **+7.561 dB @ 11.11 rad/s**, che corrisponde a k = 2.39, esattamente il bordo
  superiore misurato per forza bruta.
- Riga 46 -- **Flex GM**: nella finestra [`w_flex`, `w_flex_hi`], si prende
  l'attraversamento **piu' vicino a `w_bending`** (modo `'near'`). Nel Task 2
  con il notch profondo il risultato e' **NaN**: il modo e' *gain-stabilizzato*,
  la curva non arriva mai a 0 dB nella banda del bending, quindi **non esiste**
  un flex crossover. Ecco perche' serve `LwBM_dB` (righe 52-58). Senza notch,
  invece, il classificatore trova flexGM = -24.17 dB @ 18.76 rad/s e
  flexPM = -148 deg @ 21.5 rad/s -- numeri di un loop palesemente instabile
  (`stable_am = 0`).

> **Possibile domanda d'esame** -- nel Task 2 `w_flex = 0.6*omega_BM = 11.34
> rad/s` e il Rigid GM cade a 11.11 rad/s. Non e' pericolosamente vicino al
> confine?
> *Risposta:* si', ed e' un punto debole onesto del setup. L'attraversamento di
> guadagno rigido sta appena **2 % sotto** il confine rigido/flex. Se si fosse
> scelto `w_flex = 0.5*omega_BM = 9.45 rad/s` quello stesso attraversamento
> verrebbe **classificato come Flex GM** (+7.56 dB) e il Rigid GM diventerebbe
> NaN. Attenzione pero' a cosa *non* succederebbe: il tuner **non se ne
> accorgerebbe**. La guardia della funzione di costo (`design_controller` riga 85)
> testa solo `isnan(mt.aeroGM_dB) || isnan(mt.rigidPM_deg)` -- **non** il Rigid GM
> -- e con `w_flex = 0.5*omega_BM` l'Aero GM resta -6.00 dB e il Rigid PM resta
> 30.0 deg: la penalita' 1e6 non scatta e l'ottimizzazione converge come se nulla
> fosse. (Prova indipendente: nel Task 1 `rigidGM_dB` **e'** NaN e il tuner gira
> senza problemi.) Il danno sarebbe quindi peggiore di un fallimento rumoroso: una
> **mis-classificazione silenziosa**, un gain margin del corpo rigido riportato in
> tabella come margine del bending. La scelta 0.6 funziona, ma il confine non e'
> generoso: e' un parametro che andrebbe verificato ogni volta che l'attuatore o il
> notch cambiano.

## Phase margins (righe 48-50)

```matlab
[mm.rigidPM_deg, mm.rigidPM_w] = pick(pm, pf, ...
    pf > opts.w_drift & pf < opts.w_flex, 'maxv');
```

- Riga 49 -- **Rigid PM**: la maschera esclude tutto cio' che sta **sotto
  `w_drift`**. E' la riga che salva il progetto. Nel Task 1 `allmargin` riporta
  tre phase margins: **-133.1 deg @ 0.161**, **-40.5 deg @ 0.222** e
  **+30.0 deg @ 2.455**. I primi due sono i due attraversamenti a 0 dB del lobo
  di deriva; solo il terzo e' il crossover del corpo rigido. Con
  `w_drift = 0.5517` la maschera li elimina e resta il 30.0 deg vero.
- Il modo di selezione e' **`'maxv'`**, cioe' il **valore massimo**, non il
  minimo. Questa e' una scelta di robustezza contro gli artefatti: se un
  attraversamento di deriva "sfuggisse" appena sopra `w_drift`, porterebbe un PM
  molto negativo e verrebbe scartato dal `max`. **Il prezzo e' che la regola e'
  ottimistica**: se esistessero *due* crossover genuinamente rigidi,
  `classify_margins` riporterebbe quello col margine **migliore** invece del
  peggiore. Nei casi effettivamente eseguiti (Task 1, Task 2, tutti i corner del
  Task 3) in banda rigida c'e' **un solo** attraversamento, quindi l'ambiguita'
  non si presenta -- ma la regola in se' non e' conservativa e va saputo.
- Riga 50 -- **Flex PM**: analogo, in banda bending, con selezione `'near'` a
  `w_bending`. NaN nel design ritenuto (notch profondo).

## Attenuazione del bending: `LwBM_dB` (righe 52-58)

```matlab
if ~isnan(opts.w_bending) && isfinite(opts.w_bending)
    mm.LwBM_dB = 20*log10(abs(freqresp(L, opts.w_bending)));
else
    mm.LwBM_dB = NaN;
end
```

- Righe 52-53 (commento) spiegano il perche': un modo **gain-stabilizzato** non ha
  crossover, quindi non ha ne' GM ne' PM in senso classico. Il suo "margine" e'
  semplicemente **quanto sta sotto gli 0 dB** alla frequenza del modo. Questo e' il
  numero che nel report vale come *bending gain margin*.
- Riga 55: `freqresp(L, w_bending)` valuta L(j*omega_BM) direttamente -- niente
  griglia, niente interpolazione, valore esatto.
- Valori misurati eseguendo il codice: **+29.0 dB senza filtro** (il modo domina
  il loop e lo destabilizza: e' il "+29 dB bending resonance" del README),
  **-18.2 dB con il notch profondo dopo il re-tune** (il modo e' attenuato ben
  oltre i 12 dB tipicamente richiesti per una gain-stabilisation). Il salto di
  ~47 dB e' tutto il lavoro del notch.
- Righe 56-57 (`else`): NaN nel Task 1 (nessun bending). La guardia `isfinite` e'
  necessaria perche' quando `w_bending` non e' dato, riga 34 lo pone uguale a
  `w_flex` = **Inf**, e `freqresp(L, Inf)` non avrebbe senso.

## Attraversamenti di deriva (righe 60-62)

```matlab
mm.drift_w = pf(pf > 0 & pf <= opts.w_drift);
```

- Riga 62: **non li butta via, li conserva**. Vengono passati a `plot_nichols_lv`
  che li disegna con una croce nera e la label esplicita *"drift 0 dB crossing
  (not a margin)"*. E' una scelta didattica corretta: nascondere gli
  attraversamenti spuri farebbe sembrare la curva di Nichols piu' semplice di
  quello che e'; marcarli come "non margini" **documenta** la scelta di
  classificazione invece di occultarla. Task 1: [0.161, 0.222] rad/s. Task 2:
  [0.158, 0.227] rad/s.

## Delay margin e flag di stabilita' (righe 64-65)

```matlab
mm.DM_s      = min(am.DelayMargin);
mm.stable_am = am.Stable;
```

- Riga 64: **il delay margin e' il minimo su TUTTI gli attraversamenti**,
  artefatti di deriva compresi. In pratica non fa danno, perche' DM = PM/omega e
  gli attraversamenti di deriva stanno a omega ~ 0.2 rad/s, quindi danno delay
  margin **enormi** (24.6 s e 25.2 s nel Task 1) e non vincono mai il `min`:
  vince il crossover rigido (0.2133 s nel Task 1, 0.1652 s nel Task 2, 0.209 s in
  V3, verificati come PM_rad/omega). **Ma la regola non e' filtrata per banda come
  le altre**, quindi in un caso patologico (crossover di deriva a frequenza
  relativamente alta con PM piccolo) potrebbe riportare un delay margin che non e'
  quello rigido. E' l'unica incoerenza metodologica del file e va dichiarata.
- Riga 65: `am.Stable` e' il flag di stabilita' che `allmargin` deriva **dal loop
  aperto L**. Attenzione: nel resto della repo il verdetto usato e' `isstable(T)`
  sul **ciclo chiuso** (`design_controller` riga 67). I due coincidono nei casi
  esaminati, ma quello autorevole e' `isstable(T)`.

## Funzione annidata `pick` (righe 67-79)

```matlab
function [v, w] = pick(vals, freqs, mask, mode, target)
    idx = find(mask);
    if isempty(idx), v = NaN; w = NaN; return; end
    switch mode
        case 'minf'    % lowest crossover frequency
            [w, j] = min(freqs(idx));
        case 'near'    % nearest a target frequency
            [~, j] = min(abs(freqs(idx) - target));  w = freqs(idx(j));
        case 'maxv'    % largest value (rigid PM)
            [~, j] = max(vals(idx));  w = freqs(idx(j));
    end
    v = vals(idx(j));
end
```

- Riga 67: firma. E' il **motore di selezione** condiviso da tutti e cinque i
  margini: una maschera logica (la banda) + una regola di scelta (quale
  attraversamento nella banda).
- Riga 69: **banda vuota -> NaN**, non errore. E' il contratto su cui si regge la
  penalita' 1e6 della funzione di costo di `design_controller` (riga 85-87): un
  NaN significa "questa banda non ha piu' un attraversamento", ed e' un'informazione
  fisica, non un fallimento numerico.
- Righe 71-77: i tre modi. `'minf'` per i bordi della banda condizionale (il primo
  attraversamento incontrato salendo in frequenza), `'near'` per i modi flessibili
  (dove si sa gia' *dove* aspettarsi il fenomeno), `'maxv'` per il PM rigido
  (difesa dagli artefatti di deriva, vedi la nota sopra).
- Nota su `'minf'`: `[w, j] = min(freqs(idx))` restituisce direttamente la
  frequenza in `w`; negli altri due rami `w` va ricostruita come `freqs(idx(j))`.
  Asimmetria stilistica, non un bug.

---

## Possibili domande d'esame

**D: Perche' `margin()` non basta su questo loop? Fai l'esempio concreto.**
R: Per due ragioni. (1) Il loop e' **condizionalmente stabile** (l'aeroshell ha
poli a destra: +1.82 e +0.029 rad/s), quindi la stabilita' vive in una **banda** di
guadagni -- misurata sul loop pieno del Task 2: k in [0.50, 2.39] -- e un singolo
gain margin non puo' descriverla. (2) L'integratore di deriva laterale fa "venire
la curva dall'alto" e genera due attraversamenti a 0 dB a bassa frequenza
(0.16 e 0.22 rad/s) con phase margin di -133 deg e -40 deg, che `allmargin`
riporta assieme al margine vero. Esempio concreto: nel corner V3 del Task 3,
`margin(L)` restituisce **PM = -8.7 deg @ 0.42 rad/s** su un sistema che
`isstable(T)` dichiara **stabile**; il classificatore restituisce il valore vero,
**18.0 deg @ 1.50 rad/s**. Nei casi nominali `margin()` per coincidenza azzecca il
numero, ma non e' affidabile dove serve.

**D: Cosa sono, fisicamente, i cinque margini classificati?**
R: **Aero GM** = bordo *inferiore* della banda di guadagni stabili: quanto puoi
*ridurre* il guadagno prima che l'instabilita' aerodinamica prenda il sopravvento.
E' negativo in dB. **Rigid GM** = bordo *superiore*: quanto puoi *aumentare* il
guadagno prima che il ritardo di fase di attuatore + delay porti il crossover di
media frequenza sul punto critico (assente con attuatore ideale). **Rigid PM** =
margine di fase al crossover del corpo rigido (~ sqrt(A_6)). **Flex GM/PM** =
margini nella banda del bending, se il modo e' *phase*-stabilizzato. **|L(omega_BM)|**
= attenuazione del modo, se invece e' *gain*-stabilizzato (il nostro caso: -18.2 dB
col notch, +29.0 dB senza).

**D: Perche' la maschera `gf > 0` sui gain margins?**
R: Perche' `allmargin` sul loop pieno restituisce un'entrata spuria a frequenza
**zero** con -283 dB, prodotta dall'integratore di posizione laterale (|L(0)| =
Inf). Senza quella maschera il "margine aerodinamico" verrebbe letto come -283 dB
e il tuner inseguirebbe un fantasma. E' la traccia numerica piu' diretta del fatto
che il canale di deriva contamina la lettura dei margini.

**D: `'maxv'` per il Rigid PM: perche' il massimo e non il minimo? Non e'
ottimistico?**
R: Si', e' ottimistico, ed e' una scelta deliberata di difesa contro gli
artefatti: se un attraversamento del lobo di deriva sfuggisse appena sopra
`w_drift`, avrebbe un PM molto negativo e il `max` lo scarterebbe automaticamente.
Il rovescio della medaglia e' che, se in banda rigida esistessero *due*
attraversamenti genuini, verrebbe riportato il piu' favorevole invece del piu'
critico. Nei casi effettivamente calcolati (Task 1, Task 2, tutti i corner del
Task 3) in banda rigida c'e' sempre **un solo** attraversamento, quindi il punto
resta teorico -- ma la regola in se' non e' conservativa.

**D: Come e' definito il delay margin qui, ed e' classificato per banda?**
R: `DM_s = min(am.DelayMargin)`, cioe' il **minimo su tutti** gli attraversamenti,
**non** filtrato per banda -- unica eccezione alla filosofia del file, e vale
segnalarla. In pratica il minimo cade sempre sul crossover rigido perche'
DM = PM_rad/omega e gli attraversamenti di deriva stanno a frequenze ~10 volte
piu' basse (delay margin di 24-25 s contro 0.17-0.21 s). Valori: 213 ms nel
Task 1, 165 ms nel Task 2 dopo il re-tune, 209 ms nel corner V3 -- tutti sopra la
soglia tipica di 100 ms per un lanciatore.

**D: Perche' il Rigid GM e' NaN nel Task 1?**
R: Perche' l'attuatore e' ideale e non aggiunge ritardo di fase: la fase del loop
tende a -90 deg alle alte frequenze e **non riattraversa mai** la fase critica.
Non esiste quindi un bordo superiore alla banda di guadagni: verificato
numericamente, il ciclo chiuso resta stabile fino a k = 1e6. La stabilita'
condizionale nel Task 1 e' quindi **a un solo lato** (solo la riduzione di
guadagno destabilizza). Il bordo superiore compare nel Task 2, quando i 20 ms di
ritardo di trasporto e il TVC del 2o ordine introducono il lag che crea
l'attraversamento a 11.1 rad/s (+7.56 dB).
