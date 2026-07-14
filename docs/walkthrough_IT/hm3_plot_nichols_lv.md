# HM3/plot_nichols_lv.m

## Ruolo del file nel progetto

Questa funzione disegna la **carta di Nichols del loop pieno con il punto
critico a (-180 deg, 0 dB)** -- la convenzione degli appunti del corso, da
1 + L = 0 <=> L = -1 = e^{-j*180 deg} -- curva che "viene dall'alto", e i
margini **classificati** da `classify_margins` marcati ciascuno alla sua
frequenza di attraversamento. E' chiamata da `main_task1.m` riga 63 e da
`main_task2.m` riga 197; il Task 3 (`main_task3.m` riga 81) e la figura di
trade dei filtri (`main_task2.m` righe 222-243) rifanno lo shift a mano per
poter sovrapporre piu' loop con **un unico** riferimento di fase comune.

Perche' non basta `nichols(L)` e via? Perche' la figura deve **raccontare la
storia** della stabilita' condizionale, e la figura di default la nasconde:

1. **La convenzione di fase.** Il punto critico di Nyquist (-1, 0) si traduce
   sulla carta di Nichols in **(-180 deg, 0 dB)**: e' dove gli appunti del corso
   leggono il GM (single-sheet chart con fase in [-360, 0]). La fisica del
   lanciatore rende la lettura immediata: la dinamica rotazionale instabile ha
   denominatore s^2 - A_6, che in continua vale -A_6 < 0 -- e quel segno meno
   **e'** la fase critica (-180 deg, ossia +180: stesso punto, fase mod 360) --
   quindi la curva **parte** esattamente sulla fase critica e il margine
   aerodinamico e' l'altezza della curva sopra il punto critico
   all'attraversamento di bassa frequenza. Il `PhaseMatching` serve a garantire
   che quell'attraversamento cada davvero a -180 e non su una copia a
   -180 + 360k. (Nota storica: D'Antuono Fig. 3.2 e Trotta mostrano la
   **stessa** carta rietichettata di +360, punto critico a +180; il codice
   usava quella rietichettatura fino a poco fa ed e' stato riallineato alla
   convenzione del corso. E' una pura rietichettatura dell'asse: **nessun
   margine cambia**.)
2. **La curva viene dall'alto.** Il canale di deriva laterale ha un integratore
   libero, quindi |L| -> Inf a DC: eseguendo il codice, a omega = 0.01 rad/s il
   loop del Task 2 sta a **+47.9 dB**. Una Nichols di libro parte da destra e
   scende; questa parte in cima.
3. **I margini vanno etichettati.** Su questa curva ci sono *cinque* possibili
   attraversamenti e due di essi (gli 0 dB del lobo di deriva) **non sono
   margini**. La funzione li marca esplicitamente come artefatti.

---

## Docstring e contratto (righe 1-22)

- Riga 1: firma `ax = plot_nichols_lv(L, mm, opts)`. Prende il loop `L` e la
  struct dei **margini gia' classificati** `mm`. La funzione **non ricalcola
  niente**: e' un puro renderer, e questo e' corretto -- la classificazione e'
  responsabilita' di `classify_margins`, il disegno di questa.
- Righe 4-6: la scelta di usare **`nicholsplot` nativo** e non un `plot` a mano.
  Il motivo dichiarato e' la **griglia M/N** standard con le etichette in dB
  (0, +/-0.25, +/-0.5, +/-1, +/-3, +/-6, +/-12, +/-20 dB). Quelle non sono linee
  di guadagno del loop aperto: sono i **contorni a modulo costante del ciclo
  chiuso** |T| = |L/(1+L)|, cioe' le M-circles di Hall trasportate sul piano di
  Nichols. Il contorno a 6 dB e' quello che nel report si dice "sfiorato" o
  "evitato": significa picco di risonanza del ciclo chiuso <= 6 dB.
- Righe 6-8: il punto critico rigido a **-180 deg** -- "the course convention,
  from 1 + L = 0 <=> L = -1" -- con le copie "wrappate" dei punti critici flex a
  -540 / -900. Con l'`xlim` di default [-720, 0] ne restano visibili **due**:
  -180 e -540.
- Righe 9-10: la nota di raccordo con la letteratura: D'Antuono Fig. 3.2 mostra
  la **stessa carta rietichettata di +360 deg** (la fase e' definita mod 360,
  quindi i due display sono la stessa curva). E' la risposta pronta se all'orale
  chiedono perche' nella tesi il punto critico sta a +180.

## Blocco `arguments` (righe 24-30)

- Riga 27: `wrange` default `[1e-3 1e3]`; entrambi i main passano `[1e-2 1e2]`.
- Riga 28: `xlim` default `[-720 0]` -- la finestra che contiene il punto
  critico a -180 e la sua copia wrappata a -540, e regge le escursioni di fase
  misurate (la fase visualizzata del loop Task 1 copre [-448, -90] deg, quella
  del Task 2 [-914, -149] deg).
- Righe 32-33: al solito, `Control:analysis:MarginUnstable` silenziato con
  ripristino via `onCleanup` (il loop e' open-loop instabile, `nicholsplot` e
  `bode` lo urlerebbero a ogni chiamata).

## Scelta della frequenza di riferimento (riga 35)

```matlab
wref = mm.rigidPM_w;  if isnan(wref), wref = mm.aeroGM_w; end
```

- Riga 35: **tutto lo shift di fase e' ancorato al crossover del corpo rigido.**
  E' la scelta giusta: e' il punto della curva di cui si vuole leggere il phase
  margin, quindi e' quello che deve cadere dove ci si aspetta. Il fallback su
  `aeroGM_w` copre il caso patologico in cui il PM rigido non esiste (banda
  vuota -> NaN). **Non c'e' un terzo fallback**: se anche `aeroGM_w` fosse NaN,
  `wref` resterebbe NaN e lo shift a riga 49 produrrebbe NaN su tutta la curva
  (la Nichols nativa verrebbe comunque disegnata, ma i marker sparirebbero). Non
  succede nei casi dell'homework, ma il codice non lo protegge.

## Nichols nativa e phase matching a -180 (righe 37-44)

```matlab
L.InputName = '';  L.OutputName = '';
h = nicholsplot(L, {opts.wrange(1), opts.wrange(2)});
setoptions(h, 'PhaseMatching','on', 'PhaseMatchingFreq', wref, ...
             'PhaseMatchingValue', -180, ...
             'Grid','on', 'XLimMode','manual','YLimMode','manual', ...
             'XLim', {opts.xlim}, 'YLim', {[-40 40]});
```

- Riga 38: azzerare `InputName`/`OutputName` serve solo a **togliere il sottotitolo
  "From: delta To: delta"** che MATLAB stampa sui modelli con IO nominati (il loop
  esce da `assemble_loop`, che nomina tutto). MATLAB passa gli argomenti per
  valore, quindi **la `L` del chiamante non viene toccata**.
- Riga 39: `nicholsplot` (non `nichols`) perche' restituisce un **handle di plot**
  su cui si puo' chiamare `setoptions`.
- Riga 40: **il cuore della convenzione** (nel sorgente le tre opzioni
  `PhaseMatching` / `PhaseMatchingFreq` / `PhaseMatchingValue` stanno tutte sulla
  stessa riga; il blocco qui sopra le manda a capo per leggibilita').
  `PhaseMatching` istruisce MATLAB a
  sommare a tutta la curva un multiplo di 360 deg tale che, alla frequenza
  `PhaseMatchingFreq = wref`, la fase risulti il piu' vicino possibile a
  `PhaseMatchingValue = -180`. Non e' una deformazione: **sommare 360 deg alla
  fase e' un'identita'** (la risposta in frequenza e' la stessa), e' solo una
  scelta di ramo. Nel Task 1 il ramo naturale di MATLAB mette gia' il crossover
  rigido a **-150 deg**, quindi lo shift risulta **0**: la convenzione del corso
  coincide con il ramo naturale. Nel Task 2, invece, la fase grezza (unwrapped)
  al crossover rigido esce a **+570 deg** e lo shift e' **-720**, che la riporta
  a **-150**. In entrambi i casi il crossover cade **30 deg a destra del punto
  critico -180** -- e il phase margin di 30 deg si **legge** dal grafico invece
  di doverlo dedurre.
- Riga 42: `YLim [-40 40]` dB. Nota: la curva a bassissima frequenza arriva a
  +48 dB, quindi **il ramo iniziale esce dall'inquadratura**. E' voluto: quel tratto
  non porta informazione (e' l'asintoto dell'integratore) e comprimerlo
  schiaccerebbe la regione interessante.

**Che aspetto ha la curva risultante** (numeri misurati eseguendo il codice sul
loop pieno del Task 2, gia' shiftati nella convenzione del corso, punto critico
a -180):

| omega [rad/s] | fase [deg] | |L| [dB] | cosa e' |
|---|---|---|---|
| 0.01  | -430.0 | +47.9  | la curva entra dall'alto (integratore di deriva) |
| 0.158 | -320.5 | 0.00   | 1o attraversamento 0 dB del lobo di deriva -- **non un margine** |
| 0.227 | -214.4 | 0.00   | 2o attraversamento 0 dB del lobo di deriva -- **non un margine** |
| 0.542 | **-180.00** | **+6.00** | **Aero GM**: attraversa la fase critica **sopra** il punto critico |
| 3.170 | -150.0 | 0.00   | **Rigid PM**: 0 dB a 30 deg dal punto critico |
| 11.11 | **-179.99** | **-7.56** | **Rigid GM**: rientra sulla fase critica **sotto** il punto critico |
| 18.90 | -230.0 | -18.2  | omega_BM: il notch tiene il modo ben sotto 0 dB (gain-stabilizzato) |
| 100   | -488.0 | -33.0  | coda ad alta frequenza |

Ed **ecco disegnata la stabilita' condizionale**: la curva attraversa la retta
fase = -180 deg **due volte**, una a +6.00 dB e una a -7.56 dB. Il punto critico
(-180 deg, 0 dB) resta **incastrato fra i due attraversamenti**. Abbassare il
guadagno di 6 dB fa scendere il primo attraversamento sul punto critico; alzarlo
di 7.56 dB fa salire il secondo. Entrambe le direzioni destabilizzano, ed e'
esattamente la banda k in [0.50, 2.39] verificata per forza bruta con
`isstable(feedback(k*L,1))`. **La curva deve tenere il punto critico incastrato
fra i due passaggi sulla fase critica** -- e su questa figura lo si vede.

> **Possibile domanda d'esame** -- negli appunti del corso il punto critico sta
> a (-180 deg, 0 dB), ma nella tesi di D'Antuono (Fig. 3.2) sta a +180. Chi ha
> ragione? E il codice quale convenzione usa?
> *Risposta:* hanno ragione entrambi, perche' +180 e -180 sono **lo stesso
> punto** (la fase e' definita modulo 360: e^{+j180} = e^{-j180} = -1). Il -180
> discende direttamente da 1 + L = 0 <=> L = -1 ed e' il ramo su cui il corso
> legge il GM (single-sheet chart con fase in [-360, 0]); D'Antuono mostra la
> stessa carta **rietichettata di +360**, una scelta diffusa nella letteratura
> lanciatori perche' il denominatore s^2 - A_6 a DC vale -A_6 < 0, che si puo'
> leggere come "+180 deg di fase". La fisica non cambia con l'etichetta: in
> entrambe le convenzioni la curva **parte** dalla fase critica e il gain margin
> aerodinamico e' "quanti dB sopra il punto critico sono a bassa frequenza".
> Nota onesta: il codice usava la rietichettatura a +180 fino a poco fa ed e'
> stato allineato alla convenzione del corso; e' una pura rietichettatura
> dell'asse e **nessun margine cambia** (`allmargin` non guarda il display --
> verificato: 6.00 dB / 30.0 deg identici prima e dopo, suite di test 35/35).
> Bonus pratico: la griglia M/N nativa di `ngrid` vive in [-360, 0], quindi le
> figure disegnate a mano (trade del Task 2, corner del Task 3) ora cascano
> dentro la griglia.

## Replica manuale dello shift per gli overlay (righe 46-49)

```matlab
wv = logspace(log10(opts.wrange(1)), log10(opts.wrange(2)), 4000);
[mag, ph] = bode(L, wv);  mag = squeeze(mag);  ph = squeeze(ph);
gdb = 20*log10(mag);
sh = ph + 360*round((-180 - interp1(wv, ph, wref))/360);
```

- Righe 47-48: si ricampiona `L` su 4000 punti log-spaziati **nella stessa
  wrange** della Nichols nativa. Serve perche' i marker vanno messi sopra la curva
  gia' disegnata, e per farlo bisogna sapere dove si trova -- ma il PhaseMatching
  di MATLAB e' interno all'oggetto plot e non e' interrogabile.
- Riga 49: **si replica a mano lo stesso shift.** `interp1(wv, ph, wref)` da' la
  fase grezza al riferimento; `(-180 - quella)/360` e' di quanti giri va
  spostata; `round` sceglie il multiplo intero **piu' vicino**, che e'
  esattamente la regola di `PhaseMatching`. Verificato numericamente: Task 1 ->
  shift **0** (fase grezza -150, gia' sul ramo giusto), Task 2 -> shift **-720**
  (fase grezza +570 -> -150). Se questa riga divergesse dalla regola interna di
  MATLAB, i marker cadrebbero **fuori** dalla curva -- ed e' l'unico modo in cui
  questa funzione puo' mentire.
- Nota sull'unico difetto reale: `interp1` fuori dal range restituisce **NaN**, e
  `plot` con NaN non disegna nulla **senza errore**. Se un margine classificato
  cadesse fuori da `wrange` (per esempio gli attraversamenti ad alta frequenza a
  129 e 847 rad/s che `allmargin` trova sul loop del Task 2, ben oltre il
  `wrange = [1e-2 1e2]` usato dai main) il suo marker **sparirebbe in silenzio**.
  Nei casi dell'homework nessun margine *classificato* cade fuori, ma e' una
  fragilita' da conoscere.

## Marker dei margini classificati (righe 51-56)

```matlab
hleg = addmark(hleg, mm.aeroGM_w,  'Aero |GM|', 'rs');
hleg = addmark(hleg, mm.rigidPM_w, 'Rigid PM',  'rd');
hleg = addmark(hleg, mm.rigidGM_w, 'Rigid GM',  'r^');
hleg = addmark(hleg, mm.flexGM_w,  'Flex GM',   'ro');
hleg = addmark(hleg, mm.flexPM_w,  'Flex PM',   'rv');
```

- Righe 52-56: un marker per **banda**, con forma diversa: quadrato (Aero GM),
  rombo (Rigid PM), triangolo su (Rigid GM), cerchio (Flex GM), triangolo giu'
  (Flex PM). Tutti rossi. I margini **assenti** (NaN) sono saltati da `addmark`
  (riga 77) e non compaiono in legenda: nel Task 1 restano solo Aero GM e Rigid PM
  (niente Rigid GM perche' l'attuatore e' ideale, niente Flex perche' non c'e'
  bending); nel Task 2 compaiono Aero GM, Rigid PM e Rigid GM, ma **non** Flex
  GM/PM, perche' il notch **gain-stabilizza** il modo e non esiste alcun
  attraversamento in banda flessibile. **L'assenza dei marker Flex e' quindi
  informazione, non una dimenticanza**: e' la firma grafica della gain
  stabilisation.

## Attraversamenti di deriva marcati come non-margini (righe 58-67)

```matlab
for wd = mm.drift_w(:)'
    hd = plot(ax, interp1(wv, sh, wd), interp1(wv, gdb, wd), 'kx', ...
              'MarkerSize', 10, 'LineWidth', 1.8, 'HandleVisibility', 'off');
end
set(hd, 'HandleVisibility', 'on', ...
        'DisplayName', 'drift 0 dB crossing (not a margin)');
```

- Righe 59-67: gli attraversamenti a 0 dB del lobo di deriva vengono disegnati con
  una **croce nera** e una label che dice a chiare lettere *"not a margin"*. E' la
  scelta editoriale piu' importante della figura: quei punti **ci sono** sulla
  curva, `allmargin` **li vede**, `margin()` potrebbe **sceglierli**; segnalarli e
  dichiararli non-margini e' piu' onesto che nasconderli.
- **Sottigliezza da conoscere** (righe 60-64): dentro il `for`, `hd` viene
  **sovrascritto** a ogni iterazione (`hd = plot(...)`, non un append). Tutti i
  marker sono creati con `HandleVisibility = 'off'`, e alla fine **solo l'ultimo**
  viene riacceso e messo in legenda. Il risultato e' quello voluto (tutte le croci
  disegnate, **una sola** voce di legenda), ma per effetto della sovrascrittura,
  non per costruzione esplicita: `hleg(end+1) = hd` funziona proprio perche' `hd`
  e' scalare. Se qualcuno "correggesse" il loop accumulando gli handle, quella riga
  darebbe errore. Il `hd = gobjects(0)` di riga 60 e' di fatto inutile (il ciclo e'
  gia' protetto dall'`isempty` di riga 59).

## Legenda (righe 68-72)

```matlab
% Call the built-in legend on the marker handles (avoid the Nichols chart's
% own overloaded legend method, which would treat handles as labels).
legend(hleg, 'Location', 'southwest', 'FontSize', 8, 'Box', 'on');
```

- Righe 69-71: **`legend` viene chiamata sugli handle dei marker**, non sugli assi.
  Il commento spiega perche': gli assi di un chart Nichols hanno un metodo
  `legend` **sovraccaricato** (quello dei `resppack` del Control System Toolbox)
  che interpreta il primo argomento come **elenco di etichette**, non di handle.
  Passando gli handle si forza il `legend` built-in di MATLAB. E' un bug scoperto
  sul campo, ed e' documentato nel codice -- bene cosi'.

## Funzione annidata `addmark` (righe 74-82)

```matlab
function hl = addmark(hl, w, name, style)
    if isnan(w), return; end
    h2 = plot(ax, interp1(wv, sh, w), interp1(wv, gdb, w), style, ...
              'MarkerSize', 10, 'LineWidth', 1.8, ...
              'DisplayName', sprintf('%s (%.2g rad/s)', name, w));
    hl(end+1) = h2;
end
```

- Riga 77: **`if isnan(w), return; end`** -- la guardia che rende inutili tutti i
  rami `if` nel corpo principale. Un margine NaN (banda senza attraversamento)
  semplicemente **non viene disegnato**, e per contratto non finisce in legenda.
- Righe 78-80: la posizione del marker viene **interpolata dalla curva
  ricampionata**, non ricalcolata: cosi' il marker cade **sulla curva disegnata**
  per costruzione, anche se la griglia `wv` non contiene esattamente la frequenza
  di attraversamento. La label include la frequenza (`%.2g rad/s`), che e' l'altra
  meta' dell'informazione: dire "PM = 30 deg" senza dire "a 2.45 rad/s" non
  identifica il crossover.
- La funzione e' **annidata**, quindi vede `ax`, `wv`, `sh`, `gdb` per chiusura --
  ecco perche' la firma e' cosi' corta.

---

## Possibili domande d'esame

**D: Cosa fa esattamente il `PhaseMatching` e perche' non e' un imbroglio?**
R: Somma a tutta la curva un multiplo intero di 360 deg scelto in modo che, alla
frequenza di riferimento (il crossover rigido), la fase risulti il piu' vicino
possibile a -180 deg. Non e' un imbroglio perche' la fase di una risposta in
frequenza e' definita **modulo 360**: sommare un giro intero non cambia il
sistema, cambia solo il ramo su cui lo si legge. Il codice replica lo stesso shift
a mano a riga 49 (`sh = ph + 360*round((-180 - ph(wref))/360)`) per poter piazzare
i marker sulla curva gia' disegnata. Task 1: shift 0 (fase grezza -150 deg, gia'
sul ramo del corso). Task 2: shift -720 (fase grezza +570 -> -150). In entrambi i
casi il PM di 30 deg si legge come distanza dal punto critico a -180.

**D: Come si legge la stabilita' condizionale su questa figura?**
R: Si guarda **quante volte** la curva attraversa la retta fase = -180 deg e **a
che quota**. Sul loop pieno del Task 2 la attraversa **due volte**: a 0.542 rad/s
con +6.00 dB (sopra il punto critico) e a 11.11 rad/s con -7.56 dB (sotto). Il
punto critico (-180 deg, 0 dB) e' quindi **bracketed** fra i due attraversamenti:
ridurre il guadagno di 6 dB fa collassare il primo sul punto critico, aumentarlo
di 7.56 dB fa salire il secondo. La banda di guadagni stabili e' k in [0.50, 2.39],
verificata per forza bruta. Nel Task 1 (attuatore ideale) il secondo
attraversamento **non esiste** -- la fase tende a -90 deg e non torna mai sulla
critica -- quindi la conditional stability e' a un solo lato.

**D: Perche' la curva "viene dall'alto"?**
R: Perche' il canale di deriva laterale z/delta ha un **integratore libero**
(posizione = doppio integrale dell'accelerazione, meno un polo cancellato):
|L(j*omega)| -> Inf per omega -> 0. Misurato: +47.9 dB a 0.01 rad/s, oltre il
limite superiore dell'asse (`YLim [-40 40]`). E' anche il motivo per cui esistono
i due attraversamenti a 0 dB a bassa frequenza (0.158 e 0.227 rad/s), che sono
**artefatti** del lobo di deriva e vengono marcati con una croce nera e la label
"not a margin".

**D: Nel Task 2 non compaiono i marker Flex GM e Flex PM. E' un bug?**
R: No, e' il risultato. Il notch profondo (zeta_N = 0.002, zeta_D = 0.7)
**gain-stabilizza** il primo modo di bending: la curva a omega_BM = 18.9 rad/s sta
a **-18.2 dB** dopo il re-tune (e' il `LwBM_dB` che `classify_margins` calcola con
`freqresp(L, w_bending)` alla riga 55, cioe' esattamente il punto della curva), quindi
**non c'e' alcun attraversamento** in banda flessibile e `classify_margins`
restituisce NaN per entrambi. `addmark` (riga 77) salta i NaN. L'assenza dei
marker Flex e' la firma grafica della gain stabilisation: il "margine" del modo non
e' un GM/PM ma l'attenuazione, riportata a parte come `LwBM_dB` e messa nel titolo
della figura da `main_task2.m` riga 199.

**D: Che cos'e' la griglia di sfondo e cosa vuol dire "sfiorare il contorno a
6 dB"?**
R: La griglia M/N di Nichols non riguarda il loop aperto: sono i **contorni a
modulo costante del ciclo chiuso** |T| = |L/(1+L)| (le M-circles di Hall
trasformate nel piano guadagno-fase) e i contorni a fase costante di T. Il
contorno etichettato "6 dB" e' il luogo dei punti in cui il ciclo chiuso ha un
picco di 6 dB. Dire che la curva "resta fuori dal contorno a 6 dB" significa che
il picco di risonanza del ciclo chiuso e' sotto 6 dB, cioe' un requisito di
smorzamento -- **un'informazione diversa** dal gain margin, che invece e' una
distanza orizzontale/verticale dal punto critico. Confonderli e' l'errore piu'
comune sulla carta di Nichols.
