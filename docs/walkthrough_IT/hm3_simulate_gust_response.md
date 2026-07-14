# HM3/simulate_gust_response.m

## Ruolo del file nel progetto

`simulate_gust_response.m` e' il **verificatore nel dominio del tempo** di HM3.
Il progetto e' fatto tutto in frequenza (carta di Nichols, `classify_margins`),
ma un margine di 6 dB / 30 gradi non dice quanto beccheggia il razzo o quanto
ugello serve. Questa funzione prende il **loop chiuso `T`** costruito da
`assemble_loop`, gli passa il profilo di raffica costruito da
`load_wind_profile`, e restituisce le storie temporali piu' i **picchi** che
diventano i numeri del report.

E' chiamata da `main_task1.m` (riga 45), `main_task2.m` (righe 183-184, sia sul
modello completo sia su quello rigido per il confronto), `main_task3.m` (per i
quattro vertici della scatola di incertezza) e `main_montecarlo.m`.

Il file e' brevissimo -- 35 righe -- ma condensa tre cose che vanno capite:

1. **Come il disturbo entra nell'anello.** `alpha_w` NON e' nel loop aperto `L`:
   e' un ingresso di disturbo che entra nel plant (seconda colonna di B) in un
   punto diverso da dove si rompe l'anello. E' per questo che serve `T` e non
   basta `L` (vedi `hm3_assemble_loop.md`).
2. **Cosa si monitora**: assetto `theta`, deflessione `delta` (il budget di
   attuatore), deriva `z`, e soprattutto il **bilancio di angolo d'attacco**
   `alpha`, che a max-q e' la grandezza dimensionante.
3. **L'indicatore di carico** `q_bar*alpha`, calcolato dai chiamanti a partire
   da `r.peak_alpha`.

> **Il punto della pagina** e' la **riga 29**, dove `alpha` viene assemblata con
> un **segno MENO** su `alpha_w`:
>
>     r.alpha = r.theta + r.zdot/w.V - r.alphaw;
>
> Il meno non e' una scelta di gusto: e' l'unico segno compatibile con l'Eq. (1)
> della traccia e con il plant che la funzione sta simulando. La sezione
> "Il segno di alpha_w" lo dimostra per algebra e ne tira le conseguenze
> fisiche, che sono controintuitive e valgono l'orale.

---

## `simulate_gust_response` (righe 1-35)

```matlab
function r = simulate_gust_response(T, w)
t = w.t(:);
u = [w.alphaw(:), zeros(numel(t),1)];   % [alpha_w, theta_ref]
y = lsim(T, u, t);
```

- **Riga 1**: firma `r = simulate_gust_response(T, w)`. Un solo output, una
  struct con tutto dentro.

- **Righe 9-12 (docstring)**: la funzione **documenta esplicitamente il segno**.
  Dice che il meno e' "the plant's own convention" e cita la colonna di disturbo
  `Bw = [0; -a1*V; 0; -A6]`. Non e' un dettaglio cosmetico: e' l'unico posto
  della repo dove la convenzione di `alpha` e' scritta a parole, ed e' li' apposta
  perche' e' il punto in cui e' facile sbagliare (vedi la nota storica piu' sotto).

- **Righe 14-17**: `arguments`. `T` deve essere `lti` (il `genss` restituito da
  `assemble_loop` lo e'), `w` una struct scalare. **Nessun controllo sui campi
  di `w`**: la funzione assume `w.t`, `w.alphaw`, `w.V`, cioe' assume che `w`
  venga da `load_wind_profile`.

### Righe 19-21 -- costruzione dell'ingresso e simulazione

- **Riga 19**: `t = w.t(:)` -- forza colonna. Il vettore tempi arriva da
  `load_wind_profile` (riga 51: `t = 0:o.dt:o.Tend`), quindi e' **uniforme** con
  `dt = 0.005 s` di default. I chiamanti passano `Tend = 80` (main_task1 riga 44,
  main_task2 riga 182) perche' il modo lento dell'anello chiuso ha
  `tau ~ 18-20 s` e serve un orizzonte lungo per vedere l'assetto tornare a zero.

- **Riga 20**: la matrice degli ingressi.

      u = [ w.alphaw(:) ,  zeros(N,1) ]
            |               |
            alpha_w         theta_ref = 0

  **L'ordine delle colonne e' posizionale e deve combaciare con l'ordine degli
  `InputName` di `T`**, che `assemble_loop` fissa alla riga 35:
  `{'alpha_w','theta_ref'}`. Se qualcuno invertisse quell'ordine in
  `assemble_loop`, questa funzione **non darebbe errore**: darebbe risultati
  sbagliati in silenzio. E' una fragilita' reale (accoppiamento implicito fra due
  file), da segnalare.

  `theta_ref = 0` per tutta la simulazione: la prova e' un **test di reiezione
  del disturbo** (regolatore), non di inseguimento. Il ramo `Kp_th*theta_ref`
  del controllore (`assemble_loop.m` riga 25) non viene mai eccitato qui.

- **Riga 21**: `y = lsim(T, u, t)`. `T` e' un `genss` con l'analysis point su
  `delta`: `lsim` lo valuta **con l'AP chiuso al valore nominale**, cioe' simula
  il vero anello chiuso. Le colonne di `y` seguono l'ordine degli `OutputName` di
  `T`, di nuovo fissato in `assemble_loop` riga 35:
  `{'theta','z','zdot','delta'}` -- stesso accoppiamento posizionale implicito.

  Nota sulla discretizzazione: `lsim` su un sistema continuo con griglia uniforme
  discretizza internamente a `dt = 5 ms` (frequenza di Nyquist ~628 rad/s). I poli
  piu' veloci del modello sono gli zeri/poli di Pade-3 a ~230 rad/s e il servo TVC
  a 70 rad/s: la griglia e' **adeguata ma non abbondante**. Il notch, con un null
  largo solo ~0.08 rad/s a 18.9 rad/s, e' invece risolto benissimo.

### Righe 23-29 -- lo spacchettamento e il bilancio di alpha

```matlab
r.theta  = y(:,1);
r.z      = y(:,2);
r.zdot   = y(:,3);
r.delta  = y(:,4);
r.alpha  = r.theta + r.zdot/w.V - r.alphaw;   % total angle of attack
```

- **Righe 25-28**: unita' di misura. `theta` e `delta` sono in **radianti** (il
  plant e' in rad; `alpha_w = v_w/V` e' un rapporto adimensionale, quindi rad),
  `z` in metri, `zdot` in m/s. Tutti i chiamanti convertono con `*180/pi` in
  stampa. `r.delta` e' la deflessione **post-attuatore** (uscita di `Wa` in
  `assemble_loop`), non il comando `u_pd`: in Task 1, con attuatore ideale,
  coincidono; in Task 2 no.

- **Riga 29 -- il bilancio di angolo d'attacco.** Questa e' la riga concettuale
  del file:

      alpha = theta + zdot/V - alpha_w

  Tre contributi, ognuno con un significato fisico distinto:

  | termine | cosa rappresenta | chi lo controlla |
  |---|---|---|
  | `+theta` | l'assetto del corpo rispetto all'inerziale | il **controllore** |
  | `+zdot/V` | l'incidenza generata dalla **deriva laterale propria** del veicolo (angolo di traiettoria relativa: velocita' laterale su velocita' assiale) | l'anello, solo indirettamente e lentamente |
  | `-alpha_w` | l'incidenza indotta dal **vento** (`= v_w/V`, `load_wind_profile.m` riga 67), **col segno meno** | nessuno: e' il disturbo |

  **Perche' il vento entra col meno.** `alpha` non e' un angolo rispetto al suolo:
  e' l'incidenza rispetto alla velocita' **relativa all'aria**. Se l'aria stessa
  si muove lateralmente con `v_w` nello stesso verso in cui il veicolo deriva,
  la velocita' laterale *relativa* e' `zdot - v_w`, e dividendo per `V`:

      alpha = theta + (zdot - v_w)/V = theta + zdot/V - alpha_w

  Il vento **sottrae** perche' un veicolo che vola gia' alla velocita' dell'aria
  non vede nessuna incidenza. Il caso limite lo conferma: con un gradino di vento
  sostenuto e il feedback di deriva spento, la simulazione converge a
  `zdot_ss = +V*alpha_w` (rapporto verificato 1.0000) con `theta_ss = 0`, e quindi
  `alpha_ss = 0`. Il veicolo viene spinto sottovento finche' non **vola insieme
  all'aria**, e a quel punto l'incidenza relativa e' esattamente zero. Con il
  segno opposto (`+alpha_w`) lo stato stazionario darebbe `alpha = 2*alpha_w`, che
  e' fisicamente assurdo.

  A max-q questo `alpha` e' la grandezza **dimensionante**: il momento
  aerodinamico e il carico di flessione sul corpo scalano con `q_bar*alpha`.

### Righe 31-34 -- le metriche di picco

```matlab
r.peak_theta = max(abs(r.theta));
r.peak_z     = max(abs(r.z));
r.peak_delta = max(abs(r.delta));
r.peak_alpha = max(abs(r.alpha));
```

Sono i quattro numeri che finiscono nel report e nelle tabelle di Task 3.
`max(abs(...))` su **tutto** l'orizzonte, quindi il **segno e' perso**: sapere
che `peak_theta = 0.261 gradi` non dice se il razzo si e' beccheggiato **verso** o
**contro** il vento -- distinzione che e' proprio quella fra load relief e load
aggravation. I chiamanti recuperano il segno solo guardando i grafici.

`r.peak_alpha` e' poi usato dai chiamanti per l'indicatore di carico
(`main_task1.m` riga 51):

    q_bar*alpha = p.qbar/1000 * peak_alpha*180/pi     [kPa * gradi]

Nota: e' un'unita' **mista** (kPa moltiplicati per gradi), non SI. Va bene come
indicatore relativo e per il confronto fra i casi di Task 3, ma per confrontarlo
con un limite di letteratura bisogna verificare che il limite sia espresso nelle
stesse unita' (spesso e' in psf-deg o Pa-rad).

---

## Il segno di alpha_w (riga 29) -- perche' e' un MENO

E' il punto di questa pagina, ed e' una domanda d'orale quasi garantita perche'
si risolve con **pura algebra sul codice**, senza interpretazioni.

### Cosa fa il plant

Dalle matrici di `build_plant_rigid.m` (righe 11-17), l'equazione del momento e':

    thetaddot = A6*theta + (A6/V)*zdot + K1*delta + (-A6)*alpha_w

              = A6 * ( theta + zdot/V - alpha_w )  +  K1*delta
                       \_______________________/
                        alpha "vista" dal plant

La colonna del disturbo e' `Bw = [0; -a1*V; 0; -A6]` (riga 17 di
`build_plant_rigid.m`, riga 28 di `build_plant_full.m`): il **segno e' negativo**
su `alpha_w`, ed e' esattamente cio' che scrive l'Eq. (1) della traccia. Lo stesso
vale nell'equazione della forza normale (seconda riga di A e di Bw):
`zddot = a1*(zdot + V*theta - V*alpha_w) + a4*theta + a3*delta`, che di nuovo si
raggruppa attorno a `theta + zdot/V - alpha_w`.

**Entrambe le righe aerodinamiche del plant -- laterale e beccheggio -- portano i
termini `a1*V*alpha` e `A6*alpha` con QUESTO alpha.** Il plant e' quindi
internamente coerente, e la riga 29 riporta **la stessa quantita' che il plant
integra davvero**. Un `+` contraddirebbe il plant che la funzione sta simulando.

### Nota storica (utile all'orale)

Fino a poco fa la riga 29 aveva un **`+ r.alphaw`**. Era un **bug di
post-processing**, non di modello: il plant e' sempre stato giusto (la colonna
`Bw` non e' mai cambiata), sbagliava solo la formula con cui si *ri-costruiva*
`alpha` a valle di `lsim` per graficarla e per calcolare `q_bar*alpha`. Il bug e'
stato corretto e la docstring ora fissa la convenzione a parole.

Vale la pena saperlo raccontare, per due ragioni:

- **Cosa NON toccava**: margini, Nichols, stabilita', tutte le time history di
  `theta`, `z`, `zdot`, `delta`, e l'intero Task 3. `alpha_w` non entra in `L(s)`,
  e la riga 29 e' pura post-elaborazione: non retroagisce sulla simulazione.
- **Cosa toccava**: i due numeri di carico (`peak_alpha` e `q_bar*alpha`), la
  figura del budget di alpha... e la **narrazione fisica**, che con il segno
  giusto si ribalta (sotto).

### Quanto pesa il segno (numeri della run reale)

Modello rigido di Task 1, raffica `severe`, `Vg` dalla dispersione `drywind`,
`Tend = 80`:

| grandezza | valore |
|---|---|
| picco `alpha_w` (vento da solo) | 0.390 gradi |
| picco `\|theta\|` | 0.261 gradi (negativo al picco di raffica: -0.178 gradi) |
| picco `\|z\|` / `\|delta\|` | 2.27 m / 0.528 gradi |
| **picco `\|alpha\|` -- codice attuale (`-alpha_w`)** | **0.577 gradi -> q_bar*alpha = 46.8 kPa deg** |
| picco `\|alpha\|` -- vecchio codice buggato (`+alpha_w`) | 0.255 gradi -> q_bar*alpha = 20.7 kPa deg |

Sul modello completo di Task 2 i valori corretti sono **0.565 gradi** e
**45.8 kPa deg**: praticamente identici, come deve essere (attuatore, ritardo e
notch cambiano l'anello, non il bilancio cinematico di `alpha`).

### La fisica: un puro attitude hold e' LOAD-AGGRAVATING

Ecco il punto controintuitivo, ed e' il vero contenuto d'orale.

Il picco di incidenza totale (**0.577 gradi**) **supera** il contributo del solo
vento (**0.390 gradi**). Non lo attenua: lo **amplifica** di circa il 50%.

Il perche' e' semplice una volta scritto il segno giusto. Il vento genera un
momento `-A6*alpha_w` che tende a spingere il muso in negativo; l'anello di
assetto, che ha come unico obiettivo `theta -> 0`, reagisce riportando l'assetto
verso zero -- e cosi' facendo **beccheggia il muso DENTRO il vento relativo**. Il
`theta` che ne risulta (-0.18 gradi al picco) entra in `alpha = theta + zdot/V -
alpha_w` con lo **stesso segno** del termine `-alpha_w`: si **somma** al vento
invece di cancellarlo.

Una legge di **puro attitude-hold e' quindi load-aggravating**, non
load-relieving. Ed e' esattamente il motivo per cui i lanciatori reali aggiungono
un **termine esplicito di load relief** (retroazione di accelerometro laterale o
di `alpha` stimata) invece di limitarsi al PD di assetto: senza quel termine, il
controllo di assetto lavora *contro* il carico strutturale.

Il velivolo, da parte sua, non offre nessun aiuto: con `A6 > 0` il centro di
pressione sta **davanti** al baricentro, il momento aerodinamico e' **divergente**,
e non esiste nessuna stabilita' a banderuola che riallinei spontaneamente il muso
col vento relativo.

### E allora a cosa servono i guadagni di deriva `Kp_z = Kd_z = -1e-3`?

Non sono un dispositivo di load relief -- o meglio, non in questa misura. Sono
anzitutto un **requisito di stabilita'**: chiudono l'**integratore libero della
posizione laterale** (il plant ha un polo esatto in `s = 0`, perche' `z` e'
l'integrale di `zdot`).

Verifica diretta, spegnendo i due guadagni (`Kp_z = Kd_z = 0`):

| | poli di anello chiuso | picco `\|alpha\|` | picco `\|z\|` |
|---|---|---|---|
| con drift feedback | -0.056 +/- 0.233i, -0.953 +/- 1.905i | 0.577 gradi | 2.27 m |
| senza (`Kp_z = Kd_z = 0`) | **0**, -0.076, -0.98 +/- 1.93i | 0.584 gradi | **9.54 m** |

Senza i guadagni di deriva resta un polo **esattamente nell'origine**: il sistema
e' solo **marginalmente stabile** e `z` non torna mai a zero. Sull'incidenza quei
guadagni incidono per l'**1%** (0.584 -> 0.577 gradi: irrilevante), ma tagliano la
deriva **da 9.5 m a 2.3 m**. La loro funzione e' contenere la deriva e chiudere
l'integratore, non alleggerire il carico.

> **Possibile domanda d'esame** -- Mi ricavi il bilancio di angolo d'attacco e mi
> dici quale dei tre termini il controllore puo' effettivamente usare?
> *Risposta:* `alpha = theta + zdot/V - alpha_w`. Il meno sul vento perche' `alpha`
> e' misurata rispetto alla velocita' **relativa all'aria**: se il veicolo deriva
> insieme al vento, l'incidenza relativa si annulla. `alpha_w = v_w/V` e' il
> disturbo puro, non controllabile. `zdot/V` e' l'incidenza dovuta alla deriva
> propria: il controllore la influenza solo indirettamente e lentamente (la deriva
> e' l'integrale dell'accelerazione laterale). L'unico termine con autorita'
> diretta e' **theta**: il TVC produce un momento (`K1*delta`) che ruota il corpo in
> ~1 secondo. **Ma qui sta il punto**: un PD che porta `theta -> 0` mette il muso
> dentro il vento relativo, e il suo contributo si **somma** a `-alpha_w` invece di
> cancellarlo. Il picco di `alpha` (0.577 gradi) risulta **maggiore** del solo
> vento (0.390 gradi): l'attitude hold **aggrava** il carico. Un vero load relief
> richiede un termine dedicato (accelerometro laterale / alpha stimata), che questo
> controllore **non ha**. I guadagni `Kp_z = Kd_z = -1e-3` non lo sostituiscono:
> servono a chiudere l'integratore di posizione laterale (senza, un polo resta in
> `s = 0`) e tagliano la deriva da 9.5 m a 2.3 m, ma sull'incidenza pesano l'1%.

---

## Possibili domande d'esame

**D: Perche' la simulazione della raffica ha bisogno di `T` e non si puo' fare
con `L`?**
R: Perche' `alpha_w` **non compare in `L`**. Il loop aperto e' la catena
controllore -> attuatore -> plant -> sensore rotta in `delta`: e' un oggetto SISO
che descrive la stabilita' e i margini. La raffica invece entra nel plant sulla
**seconda colonna di B**, in un punto diverso dalla rottura. `T`, costruito con
`connect` in `assemble_loop`, mantiene `alpha_w` come ingresso esterno proprio per
questo. La regola generale: `L` risponde a "e' stabile e con che margine",
`T` risponde a "quanto beccheggia e quanto ugello serve". Corollario utile: siccome
`alpha_w` non entra in `L`, **nessun margine dipende dal segno della riga 29**.

**D: Che raffica e', e da dove viene l'ampiezza?**
R: Di default e' una **raffica 1-cosine** (`load_wind_profile.m` riga 56):
`v_w(t) = 0.5*Vg*(1 - cos(2*pi*(t-t0)/Tg))`, durata `Tg = 3 s`, onset `t0 = 1 s`.
E' il profilo di raffica discreta standard (nessun contenuto ad alta frequenza,
transizione dolce). L'ampiezza `Vg` non e' inventata: viene interpolata dalla
dispersione `sigma` del file `drywind.mat` all'altitudine di max-q (15143 m) per
la severita' richiesta (`severe` di default) -- righe 39-48 di
`load_wind_profile.m`. Ne risulta un picco di `alpha_w` di 0.390 gradi. Esiste
anche il profilo `strongwind`, che gira il modello Simulink del professore
(vento medio + Dryden schedulato in quota) e finestra il risultato intorno a
max-q.

**D: Che unita' hanno le uscite, e dove si converte?**
R: Tutto **SI e radianti** dentro la funzione: `theta`, `delta`, `alpha`,
`alpha_w` in rad; `z` in m; `zdot` in m/s. La conversione in gradi la fanno i
chiamanti (`*180/pi` nelle `fprintf` e nei plot). L'indicatore di carico
`q_bar*alpha` viene calcolato **fuori** da questa funzione, in `main_task*.m`,
come `p.qbar/1000 * peak_alpha*180/pi`, quindi in kPa*gradi -- unita' mista, non
SI, ottima per confrontare i casi di Task 3 fra loro ma da maneggiare con cura se
la si confronta con un limite di letteratura.

**D: `max(abs(...))` -- che informazione perdi?**
R: Il **segno**, e con esso la lettura fisica del transitorio. `peak_theta =
0.261 gradi` non dice che al picco di raffica `theta = -0.178 gradi`, cioe' che il
veicolo ruota **nel verso opposto** al segno del disturbo. E siccome `alpha =
theta + zdot/V - alpha_w`, quel `theta` negativo si **somma** a `-alpha_w`: e'
proprio il meccanismo di load aggravation. Il picco scalare non lo mostra --
bisogna guardare `task1_alpha_load.png`, dove i tre contributi sono disegnati
separatamente e si vede che `theta` e `-alpha_w` vanno nella stessa direzione.

**D: `theta_ref` e' sempre zero. Non dovresti provare anche l'inseguimento?**
R: Nella finestra di max-q il compito e' **regolazione**, non tracking: il
programma di beccheggio nominale e' gia' stato eseguito e in quei 10-20 secondi il
riferimento e' congelato. La prova che conta e' la **reiezione del disturbo**
(raffica), che e' quello che questa funzione fa. Il ramo `Kp_th*theta_ref` esiste
in `assemble_loop` (riga 25) e `T` espone `theta_ref` come secondo ingresso, quindi
un test di step di assetto **sarebbe possibile senza toccare nulla** -- basta
mettere una colonna non nulla alla riga 20. Il codice, semplicemente, non lo fa.

**D: Il modello di Task 2 ha un notch largo 0.08 rad/s a 18.9 rad/s. La griglia
di `lsim` a 5 ms basta?**
R: Si', ampiamente. `dt = 5 ms` da' una frequenza di Nyquist di ~628 rad/s, contro
i poli piu' veloci del modello (Pade-3 a ~230 rad/s, servo TVC a 70 rad/s). Il
null del notch, essendo **stretto in frequenza ma centrato in basso** (18.9 rad/s),
e' risolto con enorme margine. Il vincolo vero su `dt` viene dai poli di Pade, non
dal notch. Nota comunque che l'orizzonte lungo (`Tend = 80 s`, imposto dal modo
lento di deriva con `tau ~ 20 s`) rende il vettore lungo 16001 punti: la
simulazione e' dominata dalla dinamica lenta, non da quella veloce.
