# Reverb and effect algorithms and the factory programs — simulation reference

This reference describes all ten algorithms of the reverb unit and all 200 factory programs at
the level needed to re-implement them sample by sample.

* **Part I** (sections 1 – 10): the reverb algorithms Chamber, Plate, Concert Hall, Inverse and Chorus+Rvb, the
  modulation sources, notes on 34 of the programs (section 7), accuracy, presentation and the analog input path.
* **Part II** (sections 11 – 19): the effect algorithms M-Band+Rvb, Glide>Hall, Res1>Plate, Res2>Plate and
  Infinite, tempo and further modulation details, a description of every one of the 200 programs (17.1) and the
  program data of all 200 programs (17.2: `effect-programs.json`, read with `effect_programs.py`).

# Part I — Reverb algorithms

Part I describes the Chamber, Plate, Concert Hall and Inverse reverb algorithms and the Chorus+Rvb algorithm
of the reverb unit, its 34 documented factory programs (section 7), the modulation sources that drive them, and the
unit's analog input path (section 10), at the level needed to re-implement them sample by sample. It describes only
the signal processing and its control, not how the original unit implements them. Rules
marked **exact** reproduce the original for every value of the parameters involved.

Sample rate: **fs = 48 kHz**. All delays are in samples at 48 kHz unless stated otherwise. The unit can also run at
44.1 kHz; section 4.15 lists everything that changes then.

## 1. Conventions

* Signals are real numbers; full scale = ±1. The original uses fixed point (24-bit coefficients/data in the
  input/output section, 18-bit delay memory and 4- or 8-bit coefficients in the tank). A double-precision float
  implementation is indistinguishable from the original.
* **Saturation.** Every signal value that the program stores saturates at ±1: every write into a delay memory (tank
  and input/output section), the program's input and output samples, the tank returns retL/retR, and the stored
  states of the input/output section (output low-pass, voice feedback, chorus line damping). A value never wraps
  around. Sums and products inside one expression are not limited; only the stored result is.
* `x[n − d]` = the signal d samples ago. Delay memories are circular buffers of 65 536 samples. Every offset is
  taken modulo 65 536, so an offset that reaches past the memory reads correspondingly old data.
* Coefficients in the tank have limited resolution: "4-bit" values are multiples of 1/16 (|k| ≤ 15/16, optionally
  plus an implicit +1), "8-bit" values multiples of 1/256. The tables give the exact effective values.
* Channels: L = left, R = right.

## 2. Building blocks

**Allpass (Schroeder), gain g, delay d.** With v = the line output (the value written d samples ago):

    w = x − g·v ;  write w into the line ;  y = v + g·w

All allpasses in the tanks have this form, or the equivalent "±g" form (w = x + g·v, y = v − g·w). Several are
**nested** (Concert Hall): the delay line of an outer allpass contains an inner allpass.

**One-pole low-pass (damping, crossover).** `y[n] = a·x[n] + b·y[n−1]` with b = 1 − a. The Rt HC and Crossover
parameters select a from a fixed table (section 4.4).

**Interpolated tap.** `y = (1 − f)·x[n − a] + f·x[n − (a+1)]`, f = fractional coefficient. Used for the modulated
Concert Hall taps (section 6.2).

**Fractional delay of the input/output section.** Lengths in units of 1/128 sample: `L = 128·D` gives
`y = f·x[n − ⌊D⌋ − 1] + (1 − f)·x[n − ⌊D⌋]` with f = D − ⌊D⌋, i.e. linear interpolation for a delay of D samples.
A length change glides or cross-fades by the rule of 4.13.

## 3. Signal flow and the input/output section

```
            ┌─────────────────────────── voices (4 delay taps, optional diffusion, feedback) ────────┐
 in L,R ──► send matrix ──► reverb send ──► [1 sample] ──► TANK (Chamber / Plate / Hall) ──► [1 sample]
                                                                                                │ retL, retR
                  return gain ──► post-delays + post mix ──► reverb width ───────────────────────┤
                                                                                                ▼
                          FX mix (voices fv, reverb fr) ──► output low-pass ──► output width ──► ×4·go ──► wet
 out = wet·(wet gain) + in·(dry gain)
```

Per sample n (all coefficients in the program tables of section 7; names in brackets):

```
# send
mL = aLL·inL + aLR·inR ;   mR = aRL·inL + aRR·inR                    [send matrix]
sendL[n] = gs·mL ;          sendR[n] = gs·mR                          [send gain gs]
# tank: receives send[n−1], produces retL/retR (section 5); the output section uses ret[n−1]
rL = gr·retL[n−1] ;  rR = gr·retR[n−1]                               [return gain gr]
# post-delays (interpolated, lengths PL, PR samples) and post mix
tL = delay(rL, PL) ;  tR = delay(rR, PR)
pL = pd·tL + pr·rL ;  pR = pd·tR + pr·rR                              [post mix]
qL = c1·pL + c2·pR ;  qR = c1·pR + c2·pL                              [reverb width]
# delay voices: two voice lines (left-fed: V1, V2; right-fed: V3, V4)
aL = AP(AP(mL; g, 741); g, 1995) ;   aR = AP(AP(mR; g, 1083); g, 1653)   [voice diffusion g]
vL = vd·aL + vp·mL + fbL ;           vR = vd·aR + vp·mR + fbR            [voice input]
t1 = delay(vL-line, V1) ;  t2 = delay(vL-line, V2) ;  t3 = delay(vR-line, V3) ;  t4 = delay(vR-line, V4)
busL = Σ t_j·(V_j → L) ;   busR = Σ t_j·(V_j → R)                     [voice levels/pans]
fbL(next sample) = fb1·t1 + fb2·t2 ;  fbR(next sample) = fb3·t3 + fb4·t4  [voice feedback]
# mix and output
oL = fv·busL + fr·qL ;  oR = fv·busR + fr·qR                          [FX mix]
sL += a·(oL − sL) ;     sR += a·(oR − sR)                             [output low-pass]
uL = w1·sL + w2·sR ;    uR = w1·sR + w2·sL                            [output width]
wetL = 4·go·uL ;        wetR = 4·go·uR
outL = wet·wetL + dry·inL ;  outR = wet·wetR + dry·inR               [Mix]
```

The voice lines and post-delays write their input before they are read (a delay of D samples returns the value
written D samples ago). The two one-sample delays around the tank are part of the algorithm: they set the timing of
the reverb relative to the voices. An overall output latency common to all paths is
irrelevant and omitted. The stored values inL/inR, retL/retR, fbL/fbR, sL/sR and outL/outR saturate at ±1 (section 1).

### 3.1 Chorus+Rvb input/output section

Chorus+Rvb runs the **Plate tank** (5.2, its own coefficient values in section 7) and replaces the voice and
post-delay section of 3 by six modulated chorus voices on two delay lines. Names in brackets are the table rows of
the Chorus+Rvb program in section 7.

```
# send
sendR[n] = c94·(c90·inR + c92·inL) ;   sendL[n] = c94·(c91·inR + c93·inL)
# tank (5.2): receives send[n−1]
# chorus lines (A = left, B = right): 2 × this sample's send (13.1); v_j = voice outputs of the previous sample
fA = 2·sendL[n] + c96·v1 + c97·v2 + c98·v3 ;    fB = 2·sendR[n] + c99·v4 + c9A·v5 + c9B·v6   [voice feedback]
lA += c9C·(fA − lA) ;                     lB += c9C·(fB − lB)                           [line damping]
write lA into line A, lB into line B (written before the reads below)
v_j = line A (j = 1..3) / line B (j = 4..6) read at delay D_j = (round(b_j) + round(m_j))/128 − 1 samples, linear interp.
chL = Σ_j cL_j·v_j ;   chR = Σ_j cR_j·v_j                                                [voice levels]
oL = cAC·chL + cAD·retL[n−1] ;   oR = cAC·chR + cAD·retR[n−1]                         [FX mix]
wL = cAE·oL + cAF·oR ;   wR = cAE·oR + cAF·oL                                           [output width]
outL = wet·4·cB0·wL + dry·inL ;   outR = wet·4·cB0·wR + dry·inR                         [output gain, Mix]
```

**Voice modulation** (all per sample, after the voices are read; 24-bit fixed point as stated):

```
b_j  += g·(B_j − b_j)                 delay glide to the target B_j (units 1/128 sample) [glide g]
m_j   = |a_j · φ_j|                   triangle LFO in 1/128 sample, used in the NEXT sample (a_j = A·2^−23 a
                                      fraction, φ_j a signed 24-bit integer)
φ_j   = (φ_j + r_j) mod 2^24          LFO phase, increment r_j (signed 24-bit integer)
a_j  += c38·(A_j − a_j)               depth smoothing, only for one voice per sample (j cycles 1..6)
rate updates: a 24-bit Weyl sequence u += U (mod 2^24, every sample); a counter N counts down every sample; when it
reaches 0, the next voice k (cycling 1..6) gets r_k = s·R_k + (|u|/2^26)·R_k with s = +1 if bit 1 of |u| is set,
else −1 (r_k rounded to 24 bits), and N = 4096 + (|u| >> 11).   (|u| = magnitude of u as a signed 24-bit number)
Order per sample: glides and LFOs (all voices), then the counter / rate update, then one depth smoothing, then u += U.
```

The targets B_j, A_j and base rates R_j follow from the voice parameters (4.12); the glide g is G[GldResp]/2^23 of
the program's DelayTime glide (4.13), and the line damping c9C follows High Cut (4.6). The Weyl step U, its start u,
the counter start N and the depth smoothing c38 are constants of the algorithm (values in the Chorus+Rvb table of
section 7). At a settled state
b_j = B_j and a_j = A_j; the LFO phases, the current rates and the Weyl state depend on the time since the program
was loaded (start values after loading: φ_j = 0, r_j = +R_j, N, u). The stored values lA, lB and the outputs saturate
at ±1 (section 1).

## 4. Parameter laws

The tables in section 7 contain the resulting values for the documented programs; the laws below produce them and
every other setting. Parameter ranges are listed in 4.14.

### 4.1 Size

**Plate and Chamber** (Chorus+Rvb uses the Plate). Every Size-dependent tank offset follows

    A(v) = A0 + ⌊(A0 − O)·v / 16⌋          (Size v = 0 … 144; O = 20 Plate, 33 Chamber)

so A(144) − O = 10·(A0 − O). Size-0 values A0:

Plate Size-0 values A0 (O = 20; 63 tank offsets):

```
A0: 1150         A1: 430          A3: 1379         A5: 2001         A6: 1151         A8: 1379
A11: 2286        A13: 2729        A14: 2001        A16: 2286        A29: 2998        A32: 2729
A33: 3853        A34: 2998        A36: 4079        A38: 4703        A39: 3853        A41: 4079
A44: 5131        A45: 937         A46: 1056        A47: 4703        A48: 453         A49: 5131
A50: 1316        A51: 1646        A52: 1999        A53: 2383        A54: 2787        A55: 3202
A56: 3576        A57: 4000        A58: 4366        A59: 4812        A60: 5277        A64: 623
A66: 1454        A67: 1885        A68: 2061        A69: 2530        A70: 2705        A71: 3018
A72: 3392        A73: 3790        A74: 4101        A75: 4532        A76: 5028        A77: 5421
A84: 33          A92: 78          A95: 33          A96: 152         A99: 78          A100: 161
A105: 152        A106: 199        A107: 257        A109: 161        A110: 288        A113: 199
A114: 5553       A123: 429        A126: 300
```

Chamber Size-0 values A0 (O = 33; 62 tank offsets):

```
A0: 420          A2: 1141         A3: 291          A5: 420          A6: 1369         A8: 1991
A9: 1142         A11: 1369        A15: 2276        A17: 2719        A18: 1991        A20: 2276
A34: 2988        A37: 2719        A38: 3843        A39: 2988        A40: 4069        A43: 3843
A44: 4693        A45: 4069        A49: 5121        A50: 928         A51: 614         A52: 4693
A53: 496         A54: 5121        A55: 1307        A56: 1637        A57: 1990        A58: 2374
A59: 2611        A60: 2778        A61: 3193        A62: 3567        A63: 3991        A64: 4357
A65: 4803        A66: 5168        A67: 5447        A71: 324         A73: 1047        A74: 1445
A75: 1876        A76: 2052        A77: 2521        A78: 2696        A79: 3009        A80: 3383
A81: 3781        A82: 4092        A83: 4523        A84: 4899        A85: 5179        A89: 78
A95: 152         A98: 78          A99: 190         A104: 152        A105: 279        A108: 190
A109: 5543       A122: 513
```

**Input taps** (Plate and Chamber): the taps of the input lines that feed the tank start at the pre-delay tap
PD = 1 + 96·(Pre Delay) (or the tempo-locked length + 1, 4.11) and scale with Size:

    Plate:    q = ⌊(16 + v)/2⌋,            A = PD + ⌊S·q/100⌋,     S = 16444 (A2, A35), 878 (A10, A43), 2390 (A28, A122)
    Chamber:  q = ⌊Spread·(16 + v)/160⌋,   A_j = PD + ⌊M_j·q/255⌋  (seven Spread taps, M_j in 4.8)

The input taps scale with Size only when **Link is on**. With Link off they stay at their Size-144 spacing: use
v = 144 in the formulas above (for the Chamber q = Spread). The tank offsets follow Size either way. The other
input-line taps (Plate A61, A78; Chamber A47, A68, A86) sit at PD; the early-reflection and echo taps (RefDly,
EkoDly) follow their own parameters (4.5). Exact for every offset and every Size.

**Concert Hall**: Size 0 … 304. The 49 unmarked Size-dependent offsets follow A(v) = A0 + ⌊(A1 − A0)·v/304⌋; the four
offsets marked * follow A(v) = A0 + ⌊(2·(A1 − A0) + 1)·v/608⌋ (their end point lies half a sample further out).
Exact for every Size. The Concert Hall input taps follow Pre Delay only (A17, A22, A31, A32, A33, A75, A76, A77 =
26 626 + PD; A24 / A26 = 26 626 + the RefDly R / L offset).

Concert Hall Size end points (A0 → A1, samples):

```
A5: 3366 → 27420         A9: 3774 → 35580         A10: 3499 → 30080        A12: 2819 → 16480
A13: 2515 → 10400        A15: 3105 → 22200        A16: 3466 → 29420        A18: 3774 → 35580
A19: 3496 → 30020        A20: 3223 → 24560        A21: 2849 → 17080        A38: 2184 → 2982
A42*: 2323 → 4441        A44: 2702 → 14140        A45*: 2185 → 2992        A46: 2817 → 16440
A47: 2514 → 10380        A48: 2741 → 14920        A50: 2850 → 17100        A51: 2741 → 14920
A53: 2742 → 14940        A54: 3095 → 22000        A55: 2818 → 16460        A56: 2966 → 19420
A58: 2879 → 17680        A61: 2879 → 17680        A63: 2880 → 17700        A64: 3221 → 24520
A65: 3095 → 22000        A66: 2967 → 19440        A67: 3096 → 22020        A79: 2396 → 5208
A86: 2324 → 4452         A87*: 2513 → 6436        A89: 3364 → 27380        A90*: 2397 → 5218
A91: 3464 → 29380        A92: 3222 → 24540        A93: 3415 → 28400        A95: 3497 → 30040
A96: 3415 → 28400        A98: 3416 → 28420        A99: 3772 → 35540        A100: 3465 → 29400
A101: 3637 → 32840       A103: 3559 → 31280       A106: 3559 → 31280       A108: 3560 → 31300
A110: 3772 → 35540       A111: 3638 → 32860       A113: 3773 → 35560       A119: 2703 → 14160
A125: 2851 → 17120
```

Size changes no coefficient directly; with Link off it changes the decay gains through the loop time K (4.2).

### 4.2 Decay (Mid Rt, Link)

Each tank segment gain is `g = ⌊256·10^(−3·K / RT[MidRt])⌋ / 256` (RT in ms; rounded down to 1/256, exact for every
value):

RT[v] in ms for Mid Rt v = 0 … 63:

```
  0-15:  242   291   329   363   395   425   454   483   512   541   570   600   629   660   691   723
 16-31:  755   789   824   860   897   935   975  1017  1061  1106  1154  1203  1256  1311  1369  1430
 32-47: 1495  1564  1637  1715  1799  1888  1984  2088  2200  2321  2453  2598  2757  2932  3126  3344
 48-63: 3588  3864  4179  4543  4967  5468  6069  6802  7719  8897 10468 12666 15963 21456 32441 65393
```

The loop time K is a whole number of milliseconds, computed from an algorithm constant L:

```
b     = ⌊L / f⌋                     f = 48 (at 44.1 kHz: f = 44, 4.15)
K_raw = b·M                         Link on
K_raw = ⌊b·(16 + Size) / 16⌋        Link off
K     = ⌊K_raw / u⌋  ms
```

| tank | L | M | u | b | K Link on | K Link off, Size 0 |
|---|---|---|---|---|---|---|
| Plate (also Chorus+Rvb) | 5253 | 10 | 6 | 109 | 181 ms | 18 ms |
| Chamber | 5252 | 10 | 6 | 109 | 181 ms | 18 ms |
| Concert Hall | 2186 | 20 | 2 | 45 | 450 ms | 22 ms |

Example: Chamber, Link off, Size 24: ⌊⌊109·40/16⌋/6⌋ = ⌊272/6⌋ = 45 ms. At Size max (144, Hall 304) the Link-off K
equals the Link-on K.

* **Link on**: K is the same for every Size, so the per-pass gain is fixed while the loop grows with Size: the decay
  time follows Size. Plate and Chamber decay with RT60 ≈ RT[Mid Rt]·(16 + Size)/160; e.g. Size 90 gives
  ≈ 0.66 × RT[Mid Rt].
* **Link off**: K grows with the loop length, so the decay time stays ≈ RT[Mid Rt] at every Size.
* The Concert Hall decays faster than its displayed time (≈ 0.7×), because its K does not match its gain-to-gain loop
  length (≈ 295 ms at Size max). That is how the original behaves.

The Inverse has no Mid Rt decay (its slopes set the gains directly, 4.9). The ADJUST knob of most programs is
patched to Mid Rt; the tables of section 7 give the effective value.

### 4.3 Low Rt (low-frequency decay multiplier)

Low Rt v selects m = 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0, 4.0. Each crossover adds `k_low = g·(g^e − 1)`
to its low band, with g the **unrounded** decay gain of 4.2, so the low band decays with g^(1+e) where the mid band
decays with g. Plate and Chamber: e = 3·(1/m − 1); Concert Hall: e = 1/m − 1. The coefficient is an 8-bit value:
truncate 256·k_low toward zero, then clamp it to ±255 (in units of 1/256). The clamp matters for long low-band decays
combined with short mid decays (Locker Room: Low Rt 8, g = 0.601, g·(g^e − 1) = 1.064, stored as 255/256). The law is
accurate to within 1.2/256 and exact for 96 % of the values (8).

### 4.4 Damping (Rt HC) and Crossover

One-pole filter `y += a·(x − y)` with b = round(256·exp(−2π·f/48 000))/256 and a = 1 − b (the two coefficients of
the filter steps are a and b; exact for every value); entry "Full": a = 1, b = 0. Table (entry: f in Hz):

`0:30 1:60 2:90 3:120 4:151 5:181 6:212 7:243 8:273 9:336 10:398 11:461 12:525 13:589 14:654 15:818 16:986 17:1158
18:1333 19:1513 20:1697 21:1886 22:2079 23:2278 24:2481 25:2691 26:2906 27:3127 28:3355 29:3591 30:3833 31:4084
32:4343 33:4611 34:4888 35:5177 36:5476 37:5788 38:6113 39:6453 40:6808 41:7181 42:7573 43:7986 44:8423 45:8886
46:9379 47:9906 48:10472 49:11084 50:11748 51:12476 52:13281 53:14181 54:15201 55:16379 56:17772 57:19476 58:21674
59:24772 60:Full (a = 1, b = 0)`

Crossover v uses entry v; Rt HC v uses entry v + 12 (v = 48 → Full = no damping).

### 4.5 Levels, mixes, widths, tank parameters

| quantity | law |
|---|---|
| input level InLvl (0–160) | v > 80: 10^((v − 160)/20); 80: 0; v < 80: −10^(−v/20) |
| input pan InPan (0–100) | own channel 1 − v/100, opposite v/100 (InPan L 0 and InPan R 100 = each input on its own channel) |
| Rvb In (send gain gs, 0–80) | 0.4375·10^((v − 80)/20), quantized to 16 bits (24-bit word = n·257); within 0.2 dB of the law for v ≥ 20, below that n = 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4, 5 for v = 0 … 11. Chorus+Rvb: no effect (c94 = 0.4375) |
| FX Adjust (output gain, 0–80) | output gain 4·go with go = (256·T[v] + (T[v] mod 256))/2^23, T from the table below (the 16-bit word rule of 11; ≈ T[v]/32768, 1/256 higher; v = 80 gives the word 0x7FFFFF). **Exact** |
| Mix (program controller, 0–100) | wet = sin(π·v/200), dry = cos(π·v/200) (equal power, like FX Mix) |
| mixes (FX Mix, PstMix, 0–100) | equal power: first path cos(π·v/200), second sin(π·v/200) (FX Mix: voices fv = cos, reverb fr = sin; PstMix: direct pr = cos, delayed pd = sin); exactly the 16-bit table ⌊32767·sin(π·v/200)⌋, second path at v, first at 100 − v (14.1; Mix the same) |
| width (Rvb Width, FX Width, 0–720, 405 = unchanged) | same = sin θ, cross = cos θ with θ = 45° + v° − 0.1°·(v mod 9) (a rotation of the two outputs; 405 gives θ = 450° ≡ 90°: same = 1, cross = 0) |
| tank level Rvb Out, RefLvl (0–15) | v/16 |
| tank times Pre Delay, EkoDly, RefDly (2-ms units) | offset = 1 + 96·v samples; tempo-locked: 4.11 |
| post-delays, delay voices, chorus voices | 4.12 |
| Voice Dif (0–100) | voice diffusion g = v/200; diffused path on for v > 0 (vd = 1, vp = 0), else vd = 0, vp = 1 |
| High Cut (output low-pass, 0–121) | a = H[v]/2^23 from the table in 4.6 (v = 121: no filter, a = 1) |
| EkoFbk (Plate, Chamber echo feedback, 0–30) | k = sign(v − 15)·2·round(127·|v − 15|/15)/256 (15 = no feedback) |
| Shape (Chamber, Inverse), Attack (Plate) | curve law of 4.7; the Chamber's "a" coefficients also follow Rt HC |
| Spread (Chamber) | tap law of 4.8 |
| Spread (Plate, Concert Hall, Inverse), Attack (Concert Hall), AP Sweep (Inverse) | no effect on the signal |
| Duration, Low Slope, Mid Slope (Inverse) | 4.9 |
| Diffusion (tank input allpasses; Plate, Chamber, Concert Hall) | two allpass gains per input, g1 = ⌊192·v/100⌋/256 (→ 0.75) and g2 = ⌊168·v/100⌋/256 (→ 0.656); the negative coefficient of each allpass is −g, except where it is ×2-coded (Hall: all four; Plate/Chamber: the g2 allpasses), there its magnitude is 2·⌊256·g/2⌋/256 (g rounded down to an even multiple of 1/256). Exact for every v |
| Diffusion (Inverse) | all six allpasses of the program (k88/k90 … k108/k110) use g1, as −g1 / +g1 |
| Def (Concert Hall, 0–100) | rest gains of the AP Sweep pairs (6.1): pairs with |g0| = 0.5 at Def 0 get ⌊128·(100 − v)/100⌋/256, those with 0.25 get ⌊64·(100 − v)/100⌋/256 (pattern over the 10 pairs: large, small, large, small, small, large, small, large, small, small); Def 100 = no allpass gain. A change takes effect through the AP Sweep steps (6.1). Exact for every v; nothing else depends on Def |
| Depth (Concert Hall, 0–15) | two 2×2 mixers in the tank (k13/k14/k15/k16 and k20/k19/k18/k21), table below; Depth also shifts the left/right balance of the tank output |
| Chorus (Concert Hall, 0–10) | modulated tap step divider D (6.2): Chorus 1 … 10 → D = 53, 42, 32, 27, 22, 17, 12, 7, 4, 1; Chorus 0 stops the taps |

T[v] of the FX Adjust gain, v = 0 … 80:

```
  0-  9:     0     1     2     3     4     5     6     7     8     9
 10- 19:    10    11    13    14    16    18    20    23    26    29
 20- 29:    32    36    41    46    51    58    65    73    82    92
 30- 39:   103   116   130   146   164   184   206   231   260   292
 40- 49:   327   367   412   462   519   582   653   733   823   923
 50- 59:  1036  1162  1304  1463  1642  1842  2067  2319  2602  2920
 60- 69:  3276  3676  4125  4628  5193  5826  6537  7335  8230  9234
 70- 79: 10361 11626 13044 14636 16422 18426 20674 23197 26027 29203
 80- 80: 32767
```

Concert Hall Depth table (coefficients in sixteenths; k13 = k20, k15 = k18, k16 = k21 = −k14 = −k19):

```
Depth v   0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
k13·16   15 15 15 14 14 14 13 12 11 10  9  7  5  3  1  0
k15·16    0  0  1  1  2  2  3  4  5  6  7  9 10 12 13 15
k16·16    5  5  6  6  7  7  8  8  9  9 10  8  6  4  2  0
```

### 4.6 High Cut tables

H[v] = a·2^23 of the output low-pass for High Cut v = 0 … 121 (−3 dB from 20 Hz at v = 0 to 19.9 kHz at v = 120):

```
  0-  9:   21845   22616   23387   24415   26214   26985   28013   28784   31611   32382
 10- 19:   35209   37008   39578   41377   44204   46003   48573   50372   53970   55769
 20- 29:   59367   62965   66307   69905   75302   78900   84297   87895   93292   97918
 30- 39:  102287  106913  115651  120277  129272  138011  147006  156001  164996  173991
 40- 49:  182986  191981  205089  214084  227448  241069  259059  272167  290157  308147
 50- 59:  326137  343871  361861  384220  406323  428939  451298  473658  496017  540736
 60- 69:  563095  607814  630173  674892  719867  764329  809048  853767  898742  943461
 70- 79:  987923 1077617 1122336 1211773 1256492 1345929 1435623 1524804 1614498 1703935
 80- 89: 1793373 1927529 2017223 2151379 2285535 2419948 2554104 2688004 2867135 2990497
 90- 99: 3107948 3225399 3372405 3519411 3666417 3813167 3930874 4077624 4371636 4518642
100-109: 4665392 4897209 4993329 5186339 5379092 5571845 5765111 5957864 6105385 6232088
110-119: 6485236 6611682 6738128 6937306 7020319 7186857 7353139 7467505 7576218 7684931
120-121: 7767686 8388607
```

At 44.1 kHz (4.15) the output low-pass uses its own table, H44[v] = a·2^23. It gives the same −3 dB frequency as H
for v ≤ 104 (up to ≈ 9 kHz); above that its corner frequencies lie higher (v = 110: 15.6 kHz instead of 14.5 kHz).
The Chorus+Rvb line damping c9C uses the same table as the output low-pass.

H44[v] = a·2^23 at 44.1 kHz, v = 0 … 121:

```
  0-  9:    23901    24672    25700    26728    28527    29555    30583    31354    34438    35209
 10- 19:    38293    40092    42919    44975    47802    49601    52685    54484    58339    60395
 20- 29:    64250    67592    71447    75302    81213    85068    90722    94577   100488   105114
 30- 39:   109997   114880   124389   129272   138525   148291   157800   167309   177075   186584
 40- 49:   196350   205603   219995   229504   243896   258545   277307   291699   310974   329993
 50- 59:   349268   368543   387561   411463   435364   459266   483167   507325   530970   579029
 60- 69:   602931   650990   674635   722438   770497   818300   866103   914162   961965  1009768
 70- 79:  1057571  1153176  1200979  1296841  1344644  1440506  1536112  1631974  1727580  1823442
 80- 89:  1919048  2062713  2158318  2301727  2445392  2589057  2732722  2876387  3067598  3195586
 90- 99:  3316635  3437684  3589059  3740435  3891810  4043185  4164234  4315866  4618617  4769992
100-109:  4921111  5155498  5251104  5442572  5634296  5825764  6017232  6208700  6352365  6473414
110-119:  6715255  6836304  6957353  7143938  7220782  7373956  7527130  7628904  7726051  7822942
120-121:  7894903  8388607
```

### 4.7 Shape and Attack curves

Each coefficient controlled by Shape (Chamber, Inverse) or Attack (Plate) follows one of five curves c0 … c4 (the
coefficient tables of section 5 list which coefficients follow Shape / Attack, section 7 gives their values), given
at x = 0, 32, …, 256 (×2^−15):

```
c0: 23000 23500 15000  3500  2000  2000  2000  2000  2000
c1:     0  5000 17000 22000 20000 20000 14500  7500  6000
c2:     0  2500  5000 10000 16350 15000 14500 15500 13000
c3:     0  1500  3000  4500  6360 10000 14500 19000 21000
c4:  9000     0     0     0     0     0     0     0     0
```

With x = v (x = 256 for v ≥ 254), i = min(⌊x/32⌋, 7) and C = c[i] + ⌊(c[i+1] − c[i])·(x − 32i)/32⌋ (floor division),
the coefficient is k = ⌊C/2048⌋/16. Curves per coefficient:

| tank | c0 | c1 | c2 | c3 | c4 |
|---|---|---|---|---|---|
| Chamber (Shape) | k13, k47 | k24, k113 ("a") | k14, k48 | k25, k114 ("a") | k68, k86 |
| Plate (Attack) | k10, k43 | k28, k122 | — | k2, k35 | k61, k78 |
| Inverse (Shape) | k11, k45 | k12, k46 | k13, k47 | k14, k48 | k70, k87 |

The Chamber coefficients marked "a" are scaled by the Rt HC damping coefficient a of the same tank (4.4; a = k20 =
k109):

    k = ⌈C·A / 2^19⌉ / 16,     A = 256·a      (A = 256 at Rt HC 48, no damping)

so they change when either Shape or Rt HC changes. Exact for every Shape and Rt HC value.

### 4.8 Spread (Chamber)

Spread moves seven taps of the input lines (Chamber A13, A14, A24, A25, A48, A113, A114) away from the pre-delay
tap PD = 1 + 96·(Pre Delay):

    q = ⌊v·(16 + Size)/160⌋ ,   A_j = PD + ⌊M_j·q/255⌋,
    M_j = 878, 13488, 4497, 20736, 8154, 2390, 16444   (A13, A14, A24, A25, A48, A113, A114)

With Link off Size counts as 144 (q = v, 4.1). Exact for every Spread and Size.

### 4.9 Inverse: Duration, Low Slope, Mid Slope

**Duration** (0–112) scales the delay offsets like Size: A(v) = A0 + ⌊(A1 − A0)·v/112⌋ for the offsets below. The
seven Spread/Shape taps A7, A8, A13, A14, A42, A46, A47 follow their own rounding and are listed in full for
Pre Delay 3; at another Pre Delay they and the taps A39, A70, A87 move by exactly 96·(Pre Delay − 3) samples
(A39 = A87 = 32 768 + PD, A70 = PD, PD = 1 + 96·(Pre Delay)). The two reflection taps follow Duration, Pre Delay and
RefDly (all offsets modulo 65 536):

    A69 = 2701 + ⌊18914·Duration/112⌋ + PD + (RefDly R offset)        = 2703 + ⌊18914·Duration/112⌋ + 96·(Pre Delay + RefDly R)
    A86 = A69 with the RefDly L offset in place of the RefDly R offset, + 32 768

(RefDly offset = 1 + 96·RefDly, or tempo-locked, 4.11.) Exact for every Duration, Pre Delay, RefDly L and RefDly R.

**Slopes** (0–32, coefficients in 1/256, clamped to ±255): k4 = k30 = k38 = k121 = 16·(Mid − 16);
k32 = k123 = 16·(Low − Mid); k17 = k51 = J[Mid] with

```
J[Mid], Mid = 0 ... 32: 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 128 124 120 116 112 102 92 84 74 66 58 48 42 40 34 32 28
```

Duration end points (A0 at 0 → A1 at 112):

```
A0: 429 → 3292         A1: 1150 → 9060        A3: 300 → 2260         A5: 429 → 3292         A6: 1378 → 10884
A9: 1151 → 9068        A11: 1378 → 10884      A15: 2285 → 18140      A19: 2728 → 21684      A20: 2000 → 15860
A22: 2285 → 18140      A25: 257 → 1916        A34: 2997 → 23836      A37: 2728 → 21684      A38: 3852 → 30676
A40: 4078 → 32484      A41: 2997 → 23836      A43: 3852 → 30676      A45: 4078 → 32484      A49: 5130 → 40900
A52: 1056 → 8308       A53: 937 → 7356        A54: 4702 → 37476      A55: 453 → 3484        A56: 5130 → 40900
A57: 1316 → 10388      A58: 1646 → 13028      A59: 1999 → 15852      A60: 2383 → 18924      A61: 2787 → 22156
A62: 3202 → 25476      A63: 3576 → 28468      A64: 4000 → 31860      A65: 4366 → 34788      A66: 4812 → 38356
A67: 5277 → 42076      A68: 603 → 4684        A71: 623 → 4844        A73: 1454 → 11492      A74: 1885 → 14940
A75: 2061 → 16348      A76: 2530 → 20100      A77: 2705 → 21500      A78: 3018 → 24004      A79: 3392 → 26996
A80: 3790 → 30180      A81: 4101 → 32668      A82: 4532 → 36116      A83: 5028 → 40084      A84: 5421 → 43228
A85: 4703 → 37484      A88: 33 → 124          A92: 78 → 484          A95: 33 → 124          A96: 152 → 1076
A99: 78 → 484          A100: 161 → 1148       A103: 152 → 1076       A104: 199 → 1452       A107: 161 → 1148
A108: 288 → 2164       A110: 5552 → 44276     A111: 199 → 1452
```

Listed taps (value for Duration 0, 1, …, 112):

```
A7:
    316   316   319   319   323   323   326   326   330   330   333   333   337   337   340   340
    344   344   347   347   350   350   354   354   357   357   361   361   364   364   368   368
    371   371   375   375   378   378   381   381   385   385   388   388   392   392   395   395
    399   399   402   402   406   406   409   409   412   412   416   416   419   419   423   423
    426   426   430   430   433   433   437   437   440   440   443   443   447   447   450   450
    454   454   457   457   461   461   464   464   468   468   471   471   474   474   478   478
    481   481   485   485   488   488   492   492   495   495   499   499   502   502   505   505
    509
A8:
    430   430   447   447   465   465   482   482   500   500   518   518   535   535   553   553
    571   571   588   588   606   606   624   624   641   641   659   659   676   676   694   694
    712   712   729   729   747   747   765   765   782   782   800   800   818   818   835   835
    853   853   870   870   888   888   906   906   923   923   941   941   959   959   976   976
    994   994  1012  1012  1029  1029  1047  1047  1064  1064  1082  1082  1100  1100  1117  1117
   1135  1135  1153  1153  1170  1170  1188  1188  1206  1206  1223  1223  1241  1241  1258  1258
   1276  1276  1294  1294  1311  1311  1329  1329  1347  1347  1364  1364  1382  1382  1400  1400
   1417
A13:
    712   712   765   765   817   817   870   870   923   923   976   976  1029  1029  1082  1082
   1135  1135  1188  1188  1241  1241  1293  1293  1346  1346  1399  1399  1452  1452  1505  1505
   1558  1558  1611  1611  1664  1664  1717  1717  1770  1770  1822  1822  1875  1875  1928  1928
   1981  1981  2034  2034  2087  2087  2140  2140  2193  2193  2246  2246  2298  2298  2351  2351
   2404  2404  2457  2457  2510  2510  2563  2563  2616  2616  2669  2669  2722  2722  2775  2775
   2827  2827  2880  2880  2933  2933  2986  2986  3039  3039  3092  3092  3145  3145  3198  3198
   3251  3251  3303  3303  3356  3356  3409  3409  3462  3462  3515  3515  3568  3568  3621  3621
   3674
A14:
    939   939  1020  1020  1102  1102  1183  1183  1264  1264  1346  1346  1427  1427  1508  1508
   1590  1590  1671  1671  1752  1752  1834  1834  1915  1915  1996  1996  2077  2077  2159  2159
   2240  2240  2321  2321  2403  2403  2484  2484  2565  2565  2647  2647  2728  2728  2809  2809
   2891  2891  2972  2972  3053  3053  3135  3135  3216  3216  3297  3297  3379  3379  3460  3460
   3541  3541  3623  3623  3704  3704  3785  3785  3866  3866  3948  3948  4029  4029  4110  4110
   4192  4192  4273  4273  4354  4354  4436  4436  4517  4517  4598  4598  4680  4680  4761  4761
   4842  4842  4924  4924  5005  5005  5086  5086  5168  5168  5249  5249  5330  5330  5412  5412
   5493
A42:
    803   803   868   868   932   932   997   997  1061  1061  1126  1126  1190  1190  1255  1255
   1319  1319  1384  1384  1448  1448  1513  1513  1577  1577  1642  1642  1706  1706  1771  1771
   1835  1835  1900  1900  1964  1964  2029  2029  2093  2093  2158  2158  2222  2222  2287  2287
   2351  2351  2416  2416  2480  2480  2545  2545  2609  2609  2673  2673  2738  2738  2802  2802
   2867  2867  2931  2931  2996  2996  3060  3060  3125  3125  3189  3189  3254  3254  3318  3318
   3383  3383  3447  3447  3512  3512  3576  3576  3641  3641  3705  3705  3770  3770  3834  3834
   3899  3899  3963  3963  4028  4028  4092  4092  4157  4157  4221  4221  4286  4286  4350  4350
   4415
A46:
  33131 33131 33141 33141 33150 33150 33160 33160 33169 33169 33178 33178 33188 33188 33197 33197
  33206 33206 33216 33216 33225 33225 33235 33235 33244 33244 33253 33253 33263 33263 33272 33272
  33281 33281 33291 33291 33300 33300 33310 33310 33319 33319 33328 33328 33338 33338 33347 33347
  33356 33356 33366 33366 33375 33375 33385 33385 33394 33394 33403 33403 33413 33413 33422 33422
  33431 33431 33441 33441 33450 33450 33460 33460 33469 33469 33478 33478 33488 33488 33497 33497
  33506 33506 33516 33516 33525 33525 33535 33535 33544 33544 33553 33553 33563 33563 33572 33572
  33581 33581 33591 33591 33600 33600 33609 33609 33619 33619 33628 33628 33638 33638 33647 33647
  33656
A47:
  33312 33312 33344 33344 33376 33376 33408 33408 33440 33440 33472 33472 33504 33504 33536 33536
  33568 33568 33600 33600 33632 33632 33664 33664 33696 33696 33728 33728 33760 33760 33792 33792
  33824 33824 33856 33856 33888 33888 33920 33920 33952 33952 33984 33984 34016 34016 34048 34048
  34080 34080 34112 34112 34144 34144 34176 34176 34208 34208 34240 34240 34272 34272 34304 34304
  34336 34336 34368 34368 34400 34400 34431 34431 34463 34463 34495 34495 34527 34527 34559 34559
  34591 34591 34623 34623 34655 34655 34687 34687 34719 34719 34751 34751 34783 34783 34815 34815
  34847 34847 34879 34879 34911 34911 34943 34943 34975 34975 35007 35007 35039 35039 35071 35071
  35103
```

### 4.10 Patches: the ADJUST knob and the modulation sources

A program can patch up to ten **sources** to parameters. A patch maps the source value x (0 … 127) through a
breakpoint curve of 2 … 8 points (x_i, y_i), y in the destination's own units, and **replaces** the stored value of
the destination; several patches on one destination add up (clamped to the parameter range, 4.14). The curve is
evaluated in integer arithmetic (exact for every position):

    below the first point: y_first;  from the last point on: y_last;
    between x_i and x_(i+1), rising:   S = ⌊(y_(i+1) − y_i)·65536 / (x_(i+1) − x_i)⌋,   y = y_i + ⌊S·(x − x_i)/65536⌋
                             falling:  S = ⌊(y_i − y_(i+1))·65536 / (x_(i+1) − x_i)⌋,   y = y_i − ⌈S·(x − x_i)/65536⌉

(a floor of the negative slope would be 1 too low at some positions, e.g. 100 → 66 over 0 … 50 at x = 25: 83,
not 82). **Destination resolution — exact:**
the curve is evaluated in the destination's own internal scale M. Each point is normalised to
u_i = ⌊(65535·y_i + max/2)/max⌋ (max = the parameter's maximum), rescaled to T_i = ⌊(u_i·M + 32767)/65535⌋, the law
above runs on T, and the result is T/(M/max) parameter units. For most parameters M = max and T_i = y_i (the law
above in parameter units). Finer scales: the delay voices, post-delays and GldRange of the reverb algorithms and
Chorus+Rvb and the chorus voice depths keep 1/32 of a unit (M = 32·max; the voices' 1/32-ms count of 4.12), and the
Glide>Hall tap delays A/B Dly 1/128 of a unit (14.3). Example: Dyna Vibrato, ADJUST 0 … 420 over the knob 0 … 127 on
B Dly, stored knob 63: B Dly = 208 + 44/128, a glide target of 1 001.046875 samples (an integer 208 would give
999.398). Several patches add in T. With the standard delay memory a patch into a post-delay is limited to 682 ms
(11.1). The front-panel **ADJUST** knob is the source of most patches; its
position is stored with the program, and the tables of section 7 give the values at that position. The other
sources are the modulators of 6.3 (LFO waveforms, envelope followers, AR envelope, time switches, input levels) and
MIDI controllers. Destinations are page-qualified parameters (e.g. "Rvb Time:Mid Rt", "Panning:Voice2"); a patch
can also drive a modulator's own settings ("Mod:" destinations, e.g. the LFO rate). A source that has no value
(6.3: a time switch in Mode 0) contributes nothing; its destination keeps the value of the other patches or its
stored value. Patches are evaluated in the 20 ms control step and reach the signal as described in 6.4.

### 4.11 Tempo locks and beat length

Each program stores a tempo as a 16-bit word (BPM = word + 40, 40 … 400; section 7.1 lists it in BPM; four
programs of Part II are stored above 250 BPM, 16.3). A time parameter is tempo-locked when its flag is set (marked
* in 7.1); its value v = 32·hi + lo then means lo/hi beats (a patch into a locked parameter uses another code,
16.6). With the tempo rate R = ⌊4 096 000/⌊960 000/BPM⌋⌋ (16.3) the beat fraction in samples is

    S = ⌊256·fs·lo / (R·hi)⌋        (= ⌊12 288 000·lo/(R·hi)⌋ at 48 kHz; one rounding, **exact**)

and B[BPM] below is S for lo/hi = 1. At 44.1 kHz S is computed with fs = 44 100, so a locked length keeps its
time (unlike the plain tank offsets and post-delays of 4.15, which keep their sample counts). Then:

* tank offsets (Pre Delay, RefDly, EkoDly): S + 1 samples;
* post-delays: S samples, not halved (a length of 131 072 samples or more wraps modulo 131 072, the range of the
  post-delay's target word; 3 beats at 40 BPM give 85 775 samples; no factory program comes near it);
* delay voices and chorus voices: the beat fraction is first converted to whole milliseconds,
  t = ⌊1000·S/fs⌋ ms (= ⌊S/48⌋ at 48 kHz), and then goes through the voice law of 4.12.

A form with two roundings, ⌊B·lo/hi⌋, is 1 sample short in some cases: 5/4 beat at 94 BPM is S = 38 304
samples = 798 ms, not 38 303 = 797 ms (16.3).

A length longer than the parameter's range (Pre Delay 465, RefDly 600 / Inverse 400 and EkoDly 600 in 2-ms units,
i.e. 96 samples each; voices 1365 ms, i.e. 48 samples each) is halved until it fits. Examples:
Deep Space Voice 3 = 9/8 beat at 120 BPM = 27 000 samples → 562 ms; BowedChords Voice 2 = 1/4 beat at 44 BPM →
342 ms (16 416 samples).

B[BPM] in samples for the program tempo 40 … 250 BPM. **Exact rule** (16.3): d = ⌊960 000/BPM⌋ (beat period in
1/16 ms), R = ⌊4 096 000/d⌋, B = ⌊12 288 000/R⌋; it reproduces every entry (tempos above 250 BPM follow the same
rule). A MIDI-clock tempo sets R directly and BeatValue does not change B (16.3):

```
 40- 49:  72282  70620  68648  67146  65710  64000  62692  61440  60234  58794
 50- 59:  57690  56626  55601  54371  53426  52512  51630  50567  49748  48956
 60- 69:  48000  47261  46545  45850  45010  44361  43729  43115  42372  41795
 70- 79:  41234  40688  40026  39511  39009  38400  37925  37463  37012  36462
 80- 89:  36035  35617  35209  34711  34324  33944  33573  33121  32768  32422
 90- 99:  32000  31670  31346  31030  30643  30340  30044  29753  29397  29118
100-109:  28845  28576  28248  27990  27738  27428  27185  26947  26713  26425
110-119:  26200  25978  25761  25493  25283  25077  24874  24625  24429  24236
120-129:  24000  23813  23630  23450  23228  23054  22882  22713  22505  22341
130-139:  22180  22021  21825  21671  21520  21333  21186  21041  20897  20721
140-149:  20582  20445  20310  20144  20013  19883  19755  19598  19473  19351
150-159:  19200  19080  18962  18846  18703  18590  18478  18367  18231  18123
160-169:  18017  17886  17782  17680  17579  17454  17355  17258  17162  17042
170-179:  16948  16855  16763  16650  16560  16471  16362  16275  16189  16104
180-189:  16000  15917  15835  15753  15653  15574  15495  15417  15321  15245
190-199:  15170  15095  15003  14930  14858  14769  14698  14628  14559  14473
200-209:  14405  14338  14271  14189  14124  14059  13995  13916  13853  13791
210-219:  13714  13653  13592  13533  13458  13400  13342  13269  13212  13156
220-229:  13100  13030  12975  12921  12867  12800  12746  12694  12641  12577
230-239:  12525  12475  12412  12362  12312  12263  12202  12154  12106  12058
240-249:  12000  11953  11906  11861  11804  11758  11714  11658  11614  11570
250-250:  11527
```

### 4.12 Post-delays, delay voices and chorus voices

**Post-delays** (PstDly L/R, ms): PL = 48·t samples (no master; tempo-locked: 4.11). The same sample count at
44.1 kHz.

**Delay voices** (Chamber, Plate, Concert Hall, Inverse: voices 1–4) and **chorus voices** (Chorus+Rvb: voices 1–6),
with t the voice time in ms (tempo-locked: 4.11) and M the Delay Master (%):

    n = ⌊32·t·M/100⌋                       the delay time in 1/32 ms
    V_j = 1.5·n samples                    (section 3 voices; at 44.1 kHz ⌊44.1·n⌋/32 samples, 4.15)
    B_j = 128 + 192·n                      (chorus delay target in 1/128 sample, D = B/128 − 1 = 1.5·n samples;
                                            at 44.1 kHz B_j = 128 + 4·⌊44.1·n⌋)

**Voice levels, pans and feedback** (both kinds of voices), with ℓ, p, f the voice's Level (0–160), Pan (0–100,
0 = left) and Feedback (0–200), and M_L, M_P, M_F the Level, Pan and Feedback Masters:

    L     = input-level law (4.5) of clamp(ℓ + M_L − 40, 0, 160), times 32767/32768
    p′    = clamp(p + M_P − 50, 0, 100)
    V_j → L = L·(1 − p′/100),   V_j → R = L·p′/100          (chorus: cL_j, cR_j)
    fb_j  = (f − 100)/100 · M_F/100                        (chorus: c96 … c9B for voices 1 … 6)

Level 80 is off (L = 0), Feedback 100 is no feedback. In section 3, fb1/fb2 feed the left line, fb3/fb4 the right.

**Chorus depth and rate** (Chorus+Rvb), with d, r the voice's Depth (0–500) and Rate (0–100), M_D, M_R the masters
MstDepth and MstRate (0–200):

```
A_j  = 128 + ⌊96·⌊16·d·M_D/50⌋ / 5⌋        depth target: peak excursion = 1 sample + 2.4·d·M_D/50 samples
                                            (the product d·M_D/50 is first floored to 1/16)
R_j  = ⌊⌊H(r)·M_R/100⌋ · 2^24 / 4 800 000⌋  base rate; H(r) = LFO frequency in 0.01 Hz:
       H = r (r ≤ 20), 20 + 2(r − 20) (≤ 50), 80 + 4(r − 50) (≤ 55), 100 + 5(r − 55) (≤ 95), 300 + 10(r − 95)
       (0.01 … 3.5 Hz); M_R scales the frequency in whole 0.01 Hz; R = f·2^24/fs
```

At 44.1 kHz: A_j = 128 + ⌊0.4·44.1·⌊16·d·M_D/50⌋⌋ (the excursion stays the same in ms) and R_j uses fs = 44 100.
MstDepth and MstRate act only through these products. Exact for every value of every voice parameter and master.
Example (Deep Space, MstDepth 50, MstRate 110): d = 52, 48, 47, 48, 41, 37 and r = 65, 50, 52, 55, 62, 68 give
A_j = 16102, 14873, 14566, 14873, 12723, 11494 and R_j = 576, 307, 335, 384, 517, 632; at 44.1 kHz A_j = 14804,
13675, 13393, 13675, 11699, 10570 and R_j = 627, 334, 365, 418, 563, 688.

A level, feedback or rate whose result is exactly 0 is handled by 6.4 (it is not written).

### 4.13 Glide and Clear

There are two GldResp / GldRange pairs: the Rvb Time page's pair acts on the post-delays, the DelayTime page's pair
on the delay voices (Chorus+Rvb: one pair, on the six chorus voices). They matter only when a delay length changes:
an edit, a patch, or a tempo change of a tempo-locked time.

* A target length that changes by at most the glide range (GldRange in ms, i.e. a range word of 48·v + 1 samples)
  is approached per sample as `len += r·(target − len)`, with r = G[GldResp]/2^23 and G = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 22, 24, 26, 28, 30, 32, 33, 36, 39, 42, 45, 48, 51, 54, 57, 60, 64, 68, 72, 77, 83, 90, 98, 107, 117, 128, 140, 153, 167, 182, 198, 215, 233, 252, 272, 293, 315, 338, 362, 387, 413, 440, 468, 497, 527, 558, 590, 623, 657, 693, 731, 769, 809, 851, 895, 941, 989, 1039, 1092, 1148, 1207, 1269, 1334, 1402, 1473, 1548, 1627, 1710, 1797, 1888, 1983, 2083, 2191, 2301, 2416, 2536, 2661, 2791, 2926, 3066, 3211, 3356, 3506) for GldResp 0 … 100.
  GldResp 50 gives 233 (time constant 36 000 samples = 0.75 s), GldResp 100 gives 3 506 (2 393 samples = 50 ms).
* A larger change jumps to the target with a 171-sample cross-fade from the old tap to the new one. GldRange 0 means
  1 sample, so practically every change is cross-faded.
* The chorus voices glide their targets B_j with g = G[GldResp]/2^23 of the DelayTime pair (3.1).

**Clear** (DelayTime page, 0/1) acts on the delay-voice levels and feedback gains (and the chorus voice levels and
feedback), which reach the signal through the input/output section's smoothing (6.4):

* **Clear = 1** removes all of them from the smoothing: each keeps the value it has at that moment, and later edits
  of these parameters do not reach the signal.
* **Clear = 0** puts them back: each glides (6.4) from its frozen value to the value the current settings give.
* Voices at 0 ms are not affected by Clear. The same hold follows every program load: a voice's level and feedback
  words stay at their load value until the time since the load reaches the voice's delay (16.4).

### 4.14 Parameter ranges

Range 0 … max of every parameter, per algorithm and edit page (values in the parameter's own units), and of the Mod
page. Parameters on no edit page are fixed by the program. Tank offsets are not clamped: an offset that reaches past
the 65 536-sample memory wraps (section 1; Concert Hall: Pre Delay 465 = 44 641 samples plus the input-tap distances).

**Chamber**

| page | parameters (max) |
|---|---|
| Controls | InLvl R 160, InLvl L 160, InPan R 100, InPan L 100, FX Mix 100, FX Adjust 80, FX Width 720, Voice Dif 100, High Cut 121 |
| Rvb Time | Rt HC 48, Low Rt 9, Crossover 60, Mid Rt 63, Pre Delay 465, EkoDly L 600, EkoDly R 600, EkoFbk L 30, EkoFbk R 30, RefDly L 600, RefDly R 600, RefLvl L 15, RefLvl R 15, PstDly L 1365, PstDly R 1365, PstMix 100, GldResp 100, GldRange 1365 |
| RvbDesign | Size 144, AP Sweep 50, Diffusion 100, Spread 255, Shape 255, Link 1, Rvb Out 15, Rvb In 80, Rvb Width 720 |
| Levels | Master 80, Voice1 160, Voice2 160, Voice3 160, Voice4 160 |
| DelayTime | Master 200, Voice1 1365, Voice2 1365, Voice3 1365, Voice4 1365, Clear 1, GldResp 100, GldRange 1365 |
| Feedback | Master 100, Voice1 200, Voice2 200, Voice3 200, Voice4 200 |
| Panning | Master 100, Voice1 100, Voice2 100, Voice3 100, Voice4 100 |

**Plate**

| page | parameters (max) |
|---|---|
| Controls | InLvl R 160, InLvl L 160, InPan R 100, InPan L 100, FX Mix 100, FX Adjust 80, FX Width 720, Voice Dif 100, High Cut 121 |
| Rvb Time | Rt HC 48, Low Rt 9, Crossover 60, Mid Rt 63, Pre Delay 465, EkoDly L 600, EkoDly R 600, EkoFbk L 30, EkoFbk R 30, RefDly L 600, RefDly R 600, RefLvl L 15, RefLvl R 15, PstDly L 1365, PstDly R 1365, PstMix 100, GldResp 100, GldRange 1365 |
| RvbDesign | Size 144, AP Sweep 50, Diffusion 100, Attack 100, Link 1, Rvb Out 15, Rvb In 80, Rvb Width 720 |
| Levels | Master 80, Voice1 160, Voice2 160, Voice3 160, Voice4 160 |
| DelayTime | Master 200, Voice1 1365, Voice2 1365, Voice3 1365, Voice4 1365, Clear 1, GldResp 100, GldRange 1365 |
| Feedback | Master 100, Voice1 200, Voice2 200, Voice3 200, Voice4 200 |
| Panning | Master 100, Voice1 100, Voice2 100, Voice3 100, Voice4 100 |
| (no page) | Spread 100 |

**Concert Hall**

| page | parameters (max) |
|---|---|
| Controls | InLvl R 160, InLvl L 160, InPan R 100, InPan L 100, FX Mix 100, FX Adjust 80, FX Width 720, Voice Dif 100, High Cut 121 |
| Rvb Time | RefDly L 600, RefDly R 600, RefLvl L 15, RefLvl R 15, Pre Delay 465, Rt HC 48, Low Rt 9, Crossover 60, Mid Rt 63, PstDly L 1365, PstDly R 1365, PstMix 100, GldResp 100, GldRange 1365 |
| RvbDesign | Diffusion 100, Def 100, Depth 15, Size 304, Rvb Out 15, Link 1, AP Sweep 50, Chorus 10, Rvb In 80, Rvb Width 720 |
| Levels | Master 80, Voice1 160, Voice2 160, Voice3 160, Voice4 160 |
| DelayTime | Master 200, Voice1 1365, Voice2 1365, Voice3 1365, Voice4 1365, Clear 1, GldResp 100, GldRange 1365 |
| Feedback | Master 100, Voice1 200, Voice2 200, Voice3 200, Voice4 200 |
| Panning | Master 100, Voice1 100, Voice2 100, Voice3 100, Voice4 100 |
| (no page) | Spread 255, Attack 255 |

**Inverse**

| page | parameters (max) |
|---|---|
| Controls | InLvl R 160, InLvl L 160, InPan R 100, InPan L 100, FX Mix 100, FX Adjust 80, FX Width 720, Voice Dif 100, High Cut 121 |
| Rvb Time | Rt HC 48, Crossover 60, Low Slope 32, Mid Slope 32, Pre Delay 465, RefDly L 400, RefDly R 400, RefLvl L 15, RefLvl R 15, PstDly L 1365, PstDly R 1365, PstMix 100, GldResp 100, GldRange 1365 |
| RvbDesign | Duration 112, Diffusion 100, Shape 255, Rvb In 80, Rvb Width 720 |
| Levels | Master 80, Voice1 160, Voice2 160, Voice3 160, Voice4 160 |
| DelayTime | Master 200, Voice1 1365, Voice2 1365, Voice3 1365, Voice4 1365, Clear 1, GldResp 100, GldRange 1365 |
| Feedback | Master 100, Voice1 200, Voice2 200, Voice3 200, Voice4 200 |
| Panning | Master 100, Voice1 100, Voice2 100, Voice3 100, Voice4 100 |
| (no page) | AP Sweep 50, Spread 255 |

**Chorus+Rvb**

| page | parameters (max) |
|---|---|
| Controls | InLvl R 160, InLvl L 160, InPan R 100, InPan L 100, High Cut 121, FX Width 720, FX Mix 100, FX Adjust 80 |
| Rvb Time | Rt HC 48, Low Rt 9, Crossover 60, Mid Rt 63, Pre Delay 465, EkoDly L 600, EkoDly R 600, EkoFbk L 30, EkoFbk R 30, RefDly L 600, RefDly R 600, RefLvl L 15, RefLvl R 15 |
| RvbDesign | Size 144, AP Sweep 50, Diffusion 100, Attack 100, Link 1, Rvb Out 15 |
| Levels | Master 80, Voice1 160, Voice2 160, Voice3 160, Voice4 160, Voice5 160, Voice6 160 |
| DelayTime | Master 200, Voice1 1365, Voice2 1365, Voice3 1365, Voice4 1365, Voice5 1365, Voice6 1365, Clear 1, GldResp 100, GldRange 1365 |
| Chorus | MstDepth 200, V1 Depth 500, V2 Depth 500, V3 Depth 500, V4 Depth 500, V5 Depth 500, V6 Depth 500, MstRate 200, V1 Rate 100, V2 Rate 100, V3 Rate 100, V4 Rate 100, V5 Rate 100, V6 Rate 100 |
| Feedback | Master 100, Voice1 200, Voice2 200, Voice3 200, Voice4 200, Voice5 200, Voice6 200 |
| Panning | Master 100, Voice1 100, Voice2 100, Voice3 100, Voice4 100, Voice5 100, Voice6 100 |
| (no page) | Rvb In 80, Spread 100 |

**Mod page (all programs)**

| page | parameters (max) |
|---|---|
| Mod | Mix 100, LFO Rate 2500, LFO Shape 5, LFO PWidth 98, Env L Release 500, Env R Release 500, AR Attack 500, AR Release 500, AR Mode 3, Sw 1 Mode 2, Sw 1 Rate 2500, Sw 1 PWidth 98, Sw 2 Mode 2, Sw 2 Rate 2500, Sw 2 PWidth 98, LFO Depth 100 |
| trigger block | trigger levels 0 … 127, compared with the source values (6.3) |

### 4.15 Operation at 44.1 kHz

When the unit runs at 44.1 kHz, exactly the following changes. Everything is exact for every parameter value, and
every quantity not listed keeps the value it has at 48 kHz.

* **Decay gains** (4.2): the loop time K uses f = 44 (not 44.1): b = 119 (Plate, Chamber) and 49 (Concert Hall), K
  Link on = 198 / 490 ms, Link off at Size 0 = 19 / 24 ms. The low-band terms of Low Rt (4.3) follow from the new g
  with the unchanged law. The Inverse has no Mid Rt decay, so its coefficients do not change.
* **High Cut** (4.6): the output low-pass and the Chorus+Rvb line damping use the H44 table.
* **Delay voices and chorus voices** (4.12): their times are kept in milliseconds: V_j = ⌊44.1·n⌋/32 samples,
  B_j = 128 + 4·⌊44.1·n⌋, A_j = 128 + ⌊0.4·44.1·⌊16·d·M_D/50⌋⌋, R_j with fs = 44 100 (tempo-locked voices take t from
  4.11 as at 48 kHz).
* **Kept in samples** (8.8 % longer in time): every tank offset (Size, Pre Delay, RefDly, EkoDly, Spread, Attack,
  the modulated-tap windows), tempo-locked offsets (B[BPM] is a sample count), the post-delays (48·t samples), the
  glide range and the glide rates, the per-sample smoothing of the input/output section (6.4). Every other tank
  coefficient is unchanged too: damping and crossover (their corner frequencies scale by 44.1/48), diffusion,
  Shape / Attack, AP Sweep rest gains, Depth and Def. The AP Sweep and modulated-tap motion run in sample time (their
  speed in Hz scales by 44.1/48). The modulators keep their 50 Hz control rate (6.3).

Because the tank delays keep their sample counts, one pass through the loop takes 48/44.1 = 1.088 times longer,
while K grows by 198/181 = 1.094 (Plate, Chamber) or 490/450 = 1.089 (Concert Hall). The mid-band decay times change
by the ratio of the two: 0.995 and 1.000. With the damping corners lowered by 44.1/48 as well, the broadband T20 of a
whole program changes a little more (Good olPlate 0.981, ConcertHall2 0.985). Apart from the decay gains, the High
Cut filter and the millisecond-based voices, the 44.1 kHz response is the 48 kHz response stretched by 48/44.1
(section 8).

## 5. Tank programs

Each tank is a fixed per-sample program of delay-memory reads and writes. The program below is exact (it reproduces the reference implementation to rounding error); section 7 supplies the program values of its symbols, the laws of section 4 compute them for any setting.

## 5.1 Chamber tank — exact program

### Structure

The Chamber is built like the Plate (input diffusers into the input lines M1/M2, one figure-of-eight tank loop in
M0, two damping one-poles, two crossovers, multi-tap outputs), with these differences:

* the six tank allpasses have gain **0.5625** (AP Sweep pairs k0/k2, k6/k8, k15/k17, k34/k36, k40/k42, k49/k51);
* the input-line taps inject into the tank at the start of two loop segments, so the early density builds up faster;
* the output taps are spread differently (15 taps per side), and the tank is shorter (Size max 144, origin 33).

### Program

Delay memory: three circular buffers of 65 536 samples, **M0** (tank), **M1** (right input line), **M2** (left input line). Per sample n, execute the operations below in order. `R s: X_s = Mb[n − A]` reads (offset A samples back), `W s: Mb[n − A] = w_s` writes the value given in the equations. Offsets written as `A<s>` vary with the program (laws of section 4, values in section 7); all others are fixed. X′ denotes the value the same read produced in the **previous** sample (partial sums that the program carries across the sample boundary). `sendL/sendR` = reverb send into the tank, `retL/retR` = tank output.

```
R0:M0[n-A0]  R2:M0[n-A2]  W3:M0[n-A3]  W5:M0[n-A5]  R6:M0[n-A6]  R8:M0[n-A8]
W9:M0[n-A9]  W11:M0[n-A11]  R13:M1[n-A13]  R14:M1[n-A14]  R15:M0[n-A15]  R17:M0[n-A17]
W18:M0[n-A18]  W20:M0[n-A20]  R22:M0[n-6]  R23:M0[n-6]  R24:M1[n-A24]  R25:M1[n-A25]
R26:M0[n-2]  W27:M0[n-5]  R28:M0[n-6]  W31:M0[n-1]  R34:M0[n-A34]  W37:M0[n-A37]
R38:M0[n-A38]  W39:M0[n-A39]  R40:M0[n-A40]  W43:M0[n-A43]  R44:M0[n-A44]  W45:M0[n-A45]
R47:M2[n-A47]  R48:M2[n-A48]  R49:M0[n-A49]  R50:M0[n-A50]  R51:M0[n-A51]  W52:M0[n-A52]
R53:M0[n-A53]  W54:M0[n-A54]  R55:M0[n-A55]  R56:M0[n-A56]  R57:M0[n-A57]  R58:M0[n-A58]
R59:M0[n-A59]  R60:M0[n-A60]  R61:M0[n-A61]  R62:M0[n-A62]  R63:M0[n-A63]  R64:M0[n-A64]
R65:M0[n-A65]  R66:M0[n-A66]  R67:M0[n-A67]  R68:M1[n-A68]  R69:M1[n-A69]  R71:M0[n-A71]
R73:M0[n-A73]  R74:M0[n-A74]  R75:M0[n-A75]  R76:M0[n-A76]  R77:M0[n-A77]  R78:M0[n-A78]
R79:M0[n-A79]  R80:M0[n-A80]  R81:M0[n-A81]  R82:M0[n-A82]  R83:M0[n-A83]  R84:M0[n-A84]
R85:M0[n-A85]  R86:M2[n-A86]  R87:M2[n-A87]  R89:M0[n-A89]  R91:M1[n-A91]  W94:M0[n-33]
R95:M0[n-A95]  W98:M0[n-A98]  R99:M0[n-A99]  W100:M1[n-0]  R101:M2[n-A101]  W104:M0[n-A104]
R105:M0[n-A105]  R106:M0[n-4]  W108:M0[n-A108]  R109:M0[n-A109]  W110:M2[n-0]  R111:M0[n-8]
R113:M2[n-A113]  R114:M2[n-A114]  W116:M0[n-7]  R118:M0[n-8]  W121:M0[n-3]  R122:M0[n-A122]
```

Values written (w_s) and tank outputs:

```
w3   (M0) = X0*k0 + X106′*k116*k122 + X111′*k120 + X118′*k118*k122
w5   (M0) = X0*k0*k2 + X0 + X106′*k116*k122*k2 + X111′*k120*k2 + X118′*k118*k122*k2
w9   (M0) = X2*k4 + X6*k6
w11  (M0) = X2*k4*k8 + X6*k6*k8 + X6
w18  (M0) = X13*k13 + X14*k14 + X15*k15 + X8*k11
w20  (M0) = X13*k13*k17 + X14*k14*k17 + X15*k15*k17 + X15 + X8*k11*k17
w27  (M0) = X17*k20 + X22*k22 + X24*k24 + X25*k25
w31  (M0) = X26*k26 + X28*k28
w37  (M0) = X23*k30 + X26*k26*k32 + X28*k28*k32 + X34*k34
w39  (M0) = X23*k30*k36 + X26*k26*k32*k36 + X28*k28*k32*k36 + X34*k34*k36 + X34
w43  (M0) = X38*k38 + X40*k40
w45  (M0) = X38*k38*k42 + X40*k40*k42 + X40
w52  (M0) = X44*k45 + X47*k47 + X48*k48 + X49*k49
w54  (M0) = X44*k45*k51 + X47*k47*k51 + X48*k48*k51 + X49*k49*k51 + X49
w94  (M0) = X89*k89 + X91*k91 + sendR
w98  (M0) = X89*k89*k93 + X89 + X91*k91*k93 + X95*k95 + k93*sendR
w100 (M1) = X89*k89*k93*k97 + X89*k97 + X91*k91*k93*k97 + X95*k95*k97 + X95 + k93*k97*sendR
w104 (M0) = X101*k101 + X99*k99 + sendL
w108 (M0) = X101*k101*k103 + X105*k105 + X99*k103*k99 + X99 + k103*sendL
w110 (M2) = X101*k101*k103*k107 + X105*k105*k107 + X105 + X99*k103*k107*k99 + X99*k107 + k103*k107*sendL
w116 (M0) = X109*k109 + X111*k111 + X113*k113 + X114*k114
w121 (M0) = X106*k116 + X118*k118
retL      = 7*X51*k88/16 + X71*k88/2 + 3*X73*k88/8 + X74*k88/2 - 3*X75*k88/8 + X76*k88/2 - 7*X77*k88/16 - 3*X78*k88/8 + X79*k88/2 - 7*X80*k88/16 + 3*X81*k88/8 - X82*k88/2 + 3*X83*k88/8 + X84*k88/2 - 3*X85*k88/8 + X86*k86*k88 + X87*k87
retR      = -3*X50*k70/8 + X53*k70/2 - X55*k70/2 - 7*X56*k70/16 - 3*X57*k70/8 + X58*k70/2 + 3*X59*k70/8 - X60*k70/2 - 7*X61*k70/16 + 3*X62*k70/8 + X63*k70/2 + 3*X64*k70/8 + X65*k70/2 - 7*X66*k70/16 + 3*X67*k70/8 + X68*k68*k70 + X69*k69
```

Coefficients `k<s>` set by the parameters (laws of section 4; values per program in section 7; every other coefficient is a number in the program):

| k | set by (law) |
|---|---|
| k0 | AP Sweep pair 1, first coefficient (6.1) |
| k2 | AP Sweep pair 1, second coefficient (6.1) |
| k4 | Size with Link off, Mid Rt, Link (4.2) |
| k6 | AP Sweep pair 2, first coefficient (6.1) |
| k8 | AP Sweep pair 2, second coefficient (6.1) |
| k11 | Size with Link off, Mid Rt, Link (4.2) |
| k13 | Shape (4.7) |
| k14 | Shape (4.7) |
| k15 | AP Sweep pair 3, first coefficient (6.1) |
| k17 | AP Sweep pair 3, second coefficient (6.1) |
| k20 | Rt HC (4.4) |
| k22 | Rt HC (4.4) |
| k24 | Rt HC (4.4); Shape (4.7) |
| k25 | Rt HC (4.4); Shape (4.7) |
| k26 | Crossover (4.4) |
| k28 | Crossover (4.4) |
| k30 | Size with Link off, Mid Rt, Link (4.2) |
| k32 | Size with Link off, Mid Rt, Link (4.2); Low Rt (4.3) |
| k34 | AP Sweep pair 4, first coefficient (6.1) |
| k36 | AP Sweep pair 4, second coefficient (6.1) |
| k38 | Size with Link off, Mid Rt, Link (4.2) |
| k40 | AP Sweep pair 5, first coefficient (6.1) |
| k42 | AP Sweep pair 5, second coefficient (6.1) |
| k45 | Size with Link off, Mid Rt, Link (4.2) |
| k47 | Shape (4.7) |
| k48 | Shape (4.7) |
| k49 | AP Sweep pair 6, first coefficient (6.1) |
| k51 | AP Sweep pair 6, second coefficient (6.1) |
| k68 | Shape (4.7) |
| k69 | RefLvl R (4.5) |
| k70 | Rvb Out (4.5) |
| k86 | Shape (4.7) |
| k87 | RefLvl L (4.5) |
| k88 | Rvb Out (4.5) |
| k89 | Diffusion (4.5) |
| k91 | EkoFbk R (4.5) |
| k93 | Diffusion (4.5) |
| k95 | Diffusion (4.5) |
| k97 | Diffusion (4.5) |
| k99 | Diffusion (4.5) |
| k101 | EkoFbk L (4.5) |
| k103 | Diffusion (4.5) |
| k105 | Diffusion (4.5) |
| k107 | Diffusion (4.5) |
| k109 | Rt HC (4.4) |
| k111 | Rt HC (4.4) |
| k113 | Rt HC (4.4); Shape (4.7) |
| k114 | Rt HC (4.4); Shape (4.7) |
| k116 | Crossover (4.4) |
| k118 | Crossover (4.4) |
| k120 | Size with Link off, Mid Rt, Link (4.2) |
| k122 | Size with Link off, Mid Rt, Link (4.2); Low Rt (4.3) |

## 5.2 Plate tank — exact program

### Structure

* **Input diffusion** (per channel): the send passes three series allpasses in M0 (right: lengths 13, 45, 74 at the
  stock settings, gains −k86 (0.75·Diffusion) and 0.656·Diffusion; left: 9, 38, 89). The result is written into the
  input line M1 (right) / M2 (left). Each chain also receives `2·k88 ×` (right) or `2·k102 ×` (left) a tap of its
  own input line: an echo with feedback (EkoDly sets the tap, EkoFbk the gain).
* **Tank**: one figure-of-eight loop in M0 made of delay segments and **six allpasses**, all with gain 0.5
  (the AP Sweep pairs, section 6.1: k123/k125, k3/k5, k11/k13, k29/k31, k36/k38, k44/k46). Four segment gains carry the
  decay (k0, k8, k33, k41 = g of 4.2). Two damping one-poles (Rt HC; k16/k18, k114/k116) and two crossovers (Low Rt;
  k20/k22 with k24/k26/k28 and k81/k83 with k118/k120/k122) sit in the loop. Taps of the input lines inject the
  signal at four points (k2, k10, k35, k43; tap positions set by Pre Delay, Size and Link, 4.1; gains by Attack, 4.7).
* **Outputs**: retR = k63·(Σ of 14 fixed-gain tank taps + k61·M1 tap) + k62·M1[RefDly R tap]; retL likewise (k78,
  k79, k80). k63/k80 = Rvb Out (v/16), k62/k79 = RefLvl R / L.

### Program

Delay memory: three circular buffers of 65 536 samples, **M0** (tank), **M1** (right input line), **M2** (left input line). Per sample n, execute the operations below in order. `R s: X_s = Mb[n − A]` reads (offset A samples back), `W s: Mb[n − A] = w_s` writes the value given in the equations. Offsets written as `A<s>` vary with the program (laws of section 4, values in section 7); all others are fixed. X′ denotes the value the same read produced in the **previous** sample (partial sums that the program carries across the sample boundary). `sendL/sendR` = reverb send into the tank, `retL/retR` = tank output.

```
R0:M0[n-A0]  W1:M0[n-A1]  R2:M1[n-A2]  R3:M0[n-A3]  R5:M0[n-A5]  W6:M0[n-A6]
W8:M0[n-A8]  R10:M1[n-A10]  R11:M0[n-A11]  R13:M0[n-A13]  W14:M0[n-A14]  W16:M0[n-A16]
R18:M0[n-6]  R19:M0[n-6]  R20:M0[n-2]  W21:M0[n-5]  R22:M0[n-6]  W25:M0[n-1]
R28:M2[n-A28]  R29:M0[n-A29]  W32:M0[n-A32]  R33:M0[n-A33]  W34:M0[n-A34]  R35:M2[n-A35]
R36:M0[n-A36]  R38:M0[n-A38]  W39:M0[n-A39]  W41:M0[n-A41]  R43:M2[n-A43]  R44:M0[n-A44]
R45:M0[n-A45]  R46:M0[n-A46]  W47:M0[n-A47]  R48:M0[n-A48]  W49:M0[n-A49]  R50:M0[n-A50]
R51:M0[n-A51]  R52:M0[n-A52]  R53:M0[n-A53]  R54:M0[n-A54]  R55:M0[n-A55]  R56:M0[n-A56]
R57:M0[n-A57]  R58:M0[n-A58]  R59:M0[n-A59]  R60:M0[n-A60]  R61:M1[n-A61]  R62:M1[n-A62]
R64:M0[n-A64]  R66:M0[n-A66]  R67:M0[n-A67]  R68:M0[n-A68]  R69:M0[n-A69]  R70:M0[n-A70]
R71:M0[n-A71]  R72:M0[n-A72]  R73:M0[n-A73]  R74:M0[n-A74]  R75:M0[n-A75]  R76:M0[n-A76]
R77:M0[n-A77]  R78:M2[n-A78]  R79:M2[n-A79]  R81:M0[n-4]  R83:M0[n-8]  R84:M0[n-A84]
W86:M0[n-3]  R88:M1[n-A88]  W91:M0[n-20]  R92:M0[n-A92]  W95:M0[n-A95]  R96:M0[n-A96]
W99:M0[n-A99]  R100:M0[n-A100]  W101:M1[n-0]  R102:M2[n-A102]  W105:M0[n-A105]  R106:M0[n-A106]
R107:M0[n-A107]  W109:M0[n-A109]  R110:M0[n-A110]  R112:M0[n-8]  W113:M0[n-A113]  R114:M0[n-A114]
W115:M2[n-0]  R116:M0[n-8]  W119:M0[n-7]  R122:M1[n-A122]  R123:M0[n-A123]  W126:M0[n-A126]
```

Values written (w_s) and tank outputs:

```
w1   (M0) = X112′*k118*k125 + X122′*k122*k125 + X123′*k123*k125 + X123′ + X81′*k120*k125*k81 + X83′*k120*k125*k83
w6   (M0) = X0*k0 + X2*k2 + X3*k3
w8   (M0) = X0*k0*k5 + X2*k2*k5 + X3*k3*k5 + X3
w14  (M0) = X10*k10 + X11*k11 + X5*k8
w16  (M0) = X10*k10*k13 + X11*k11*k13 + X11 + X5*k13*k8
w21  (M0) = X13*k16 + X18*k18
w25  (M0) = X20*k20 + X22*k22
w32  (M0) = X19*k24 + X20*k20*k26 + X22*k22*k26 + X28*k28 + X29*k29
w34  (M0) = X19*k24*k31 + X20*k20*k26*k31 + X22*k22*k26*k31 + X28*k28*k31 + X29*k29*k31 + X29
w39  (M0) = X33*k33 + X35*k35 + X36*k36
w41  (M0) = X33*k33*k38 + X35*k35*k38 + X36*k36*k38 + X36
w47  (M0) = X38*k41 + X43*k43 + X44*k44
w49  (M0) = X38*k41*k46 + X43*k43*k46 + X44*k44*k46 + X44
w86  (M0) = X81*k81 + X83*k83
w91  (M0) = X84*k86 + X88*k88 + sendR
w95  (M0) = X84*k86*k90 + X84 + X88*k88*k90 + X92*k92 + k90*sendR
w99  (M0) = X84*k86*k90*k94 + X84*k94 + X88*k88*k90*k94 + X92*k92*k94 + X92 + X96*k96 + k90*k94*sendR
w101 (M1) = X84*k86*k90*k94*k98 + X84*k94*k98 + X88*k88*k90*k94*k98 + X92*k92*k94*k98 + X92*k98 + X96*k96*k98 + X96 + k90*k94*k98*sendR
w105 (M0) = X100*k100 + X102*k102 + sendL
w109 (M0) = X100*k100*k104 + X100 + X102*k102*k104 + X106*k106 + k104*sendL
w113 (M0) = X100*k100*k104*k108 + X100*k108 + X102*k102*k104*k108 + X106*k106*k108 + X106 + X110*k110 + k104*k108*sendL
w115 (M2) = X100*k100*k104*k108*k112 + X100*k108*k112 + X102*k102*k104*k108*k112 + X106*k106*k108*k112 + X106*k112 + X110*k110*k112 + X110 + k104*k108*k112*sendL
w119 (M0) = X114*k114 + X116*k116
w126 (M0) = X112*k118 + X122*k122 + X123*k123 + X81*k120*k81 + X83*k120*k83
retL      = 3*X46*k80/8 + X64*k80/2 + X66*k80/2 - 3*X67*k80/8 + X68*k80/2 - 7*X69*k80/16 - 3*X70*k80/8 + X71*k80/2 - 7*X72*k80/16 + 3*X73*k80/8 - X74*k80/2 + 3*X75*k80/8 + X76*k80/2 - 3*X77*k80/8 + X78*k78*k80 + X79*k79
retR      = -3*X45*k63/8 + X48*k63/2 - X50*k63/2 - 7*X51*k63/16 - 3*X52*k63/8 + X53*k63/2 - X54*k63/2 - 7*X55*k63/16 + 3*X56*k63/8 + X57*k63/2 + 3*X58*k63/8 + X59*k63/2 - 3*X60*k63/8 + X61*k61*k63 + X62*k62
```

Coefficients `k<s>` set by the parameters (laws of section 4; values per program in section 7; every other coefficient is a number in the program):

| k | set by (law) |
|---|---|
| k0 | Size with Link off, Mid Rt, Link (4.2) |
| k2 | Attack (4.7) |
| k3 | AP Sweep pair 2, first coefficient (6.1) |
| k5 | AP Sweep pair 2, second coefficient (6.1) |
| k8 | Size with Link off, Mid Rt, Link (4.2) |
| k10 | Attack (4.7) |
| k11 | AP Sweep pair 3, first coefficient (6.1) |
| k13 | AP Sweep pair 3, second coefficient (6.1) |
| k16 | Rt HC (4.4) |
| k18 | Rt HC (4.4) |
| k20 | Crossover (4.4) |
| k22 | Crossover (4.4) |
| k24 | Size with Link off, Mid Rt, Link (4.2) |
| k26 | Size with Link off, Mid Rt, Link (4.2); Low Rt (4.3) |
| k28 | Attack (4.7) |
| k29 | AP Sweep pair 4, first coefficient (6.1) |
| k31 | AP Sweep pair 4, second coefficient (6.1) |
| k33 | Size with Link off, Mid Rt, Link (4.2) |
| k35 | Attack (4.7) |
| k36 | AP Sweep pair 5, first coefficient (6.1) |
| k38 | AP Sweep pair 5, second coefficient (6.1) |
| k41 | Size with Link off, Mid Rt, Link (4.2) |
| k43 | Attack (4.7) |
| k44 | AP Sweep pair 6, first coefficient (6.1) |
| k46 | AP Sweep pair 6, second coefficient (6.1) |
| k61 | Attack (4.7) |
| k62 | RefLvl R (4.5) |
| k63 | Rvb Out (4.5) |
| k78 | Attack (4.7) |
| k79 | RefLvl L (4.5) |
| k80 | Rvb Out (4.5) |
| k81 | Crossover (4.4) |
| k83 | Crossover (4.4) |
| k86 | Diffusion (4.5) |
| k88 | EkoFbk R (4.5) |
| k90 | Diffusion (4.5) |
| k92 | Diffusion (4.5) |
| k94 | Diffusion (4.5) |
| k96 | Diffusion (4.5) |
| k98 | Diffusion (4.5) |
| k100 | Diffusion (4.5) |
| k102 | EkoFbk L (4.5) |
| k104 | Diffusion (4.5) |
| k106 | Diffusion (4.5) |
| k108 | Diffusion (4.5) |
| k110 | Diffusion (4.5) |
| k112 | Diffusion (4.5) |
| k114 | Rt HC (4.4) |
| k116 | Rt HC (4.4) |
| k118 | Size with Link off, Mid Rt, Link (4.2) |
| k120 | Size with Link off, Mid Rt, Link (4.2); Low Rt (4.3) |
| k122 | Attack (4.7) |
| k123 | AP Sweep pair 1, first coefficient (6.1) |
| k125 | AP Sweep pair 1, second coefficient (6.1) |

## 5.3 Concert Hall tank — exact program

### Structure

* **Input**: both inputs are written into the input lines M1/M2. Taps of those lines feed the tank (Pre Delay), and two of them form the early reflections (RefDly / RefLvl).
* **Tank**: one figure-of-eight loop in M0, two halves. Each half starts with the decay gain
  (k27, k71; low-band terms k29/k73 with the crossovers), a damping one-pole (Rt HC) and two
  diffusion allpasses. It continues with **two nested allpass blocks**: an outer allpass (±g form) whose delay line
  contains an inner allpass, whose own delay line runs through a short **modulation line** (written at offsets 255, 767,
  1279, 1791). The modulation line is read back by one of four **interpolated, slowly moving taps** (section 6.2) and
  re-written into the main loop. The ten allpass pairs of the nested blocks are the AP Sweep pairs (k46/k52, k48/k50,
  k54/k64, k56/k62, k58/k60, k91/k97, k93/k95, k99/k109, k101/k107, k103/k105). Their rest gains are set by
  Def (4.5), e.g. 0.43 / 0.21 for Concert Hall.
* **Outputs**: retR = k23·(Σ of 4 tank taps + k17·M1 tap) + k24·M1[RefDly R tap]; retL likewise (k22, k25, k26).
* Loop length at Size 304 ≈ 27 300–29 600 samples, with two decay gains per circulation.

### Program

Delay memory: three circular buffers of 65 536 samples, **M0** (tank), **M1** (right input line), **M2** (left input line). Per sample n, execute the operations below in order. `R s: X_s = Mb[n − A]` reads (offset A samples back), `W s: Mb[n − A] = w_s` writes the value given in the equations. Offsets written as `A<s>` vary with the program (laws of section 4, values in section 7); all others are fixed. X′ denotes the value the same read produced in the **previous** sample (partial sums that the program carries across the sample boundary). `sendL/sendR` = reverb send into the tank, `retL/retR` = tank output.

```
R3:M0[n-A3]  R4:M0[n-A4]  W5:M0[n-A5]  R9:M0[n-A9]  W10:M0[n-A10]  R11:M0[n-2]
R12:M0[n-A12]  R13:M0[n-A13]  W14:M0[n-1]  R15:M0[n-A15]  R16:M0[n-A16]  R17:M1[n-A17]
R18:M0[n-A18]  R19:M0[n-A19]  R20:M0[n-A20]  R21:M0[n-A21]  R22:M2[n-A22]  R24:M1[n-A24]
R26:M2[n-A26]  R31:M1[n-A31]  R32:M1[n-A32]  R33:M1[n-A33]  R34:M0[n-6]  R38:M0[n-A38]
W39:M0[n-5]  W41:M0[n-2100]  R42:M0[n-A42]  R44:M0[n-A44]  W45:M0[n-A45]  R46:M0[n-A46]
W47:M0[n-A47]  R48:M0[n-A48]  R50:M0[n-A50]  R51:M0[n-A51]  W52:M0[n-255]  W53:M0[n-A53]
R54:M0[n-A54]  W55:M0[n-A55]  R56:M0[n-A56]  R58:M0[n-A58]  R61:M0[n-A61]  W62:M0[n-767]
W63:M0[n-A63]  R64:M0[n-A64]  R65:M0[n-A65]  W66:M0[n-A66]  W67:M0[n-A67]  R68:M0[n-4]
W71:M0[n-3]  R74:M0[n-257]  R75:M2[n-A75]  R76:M2[n-A76]  R77:M2[n-A77]  R78:M0[n-8]
R79:M0[n-A79]  W83:M0[n-7]  W86:M0[n-A86]  R87:M0[n-A87]  R89:M0[n-A89]  W90:M0[n-A90]
R91:M0[n-A91]  W92:M0[n-A92]  R93:M0[n-A93]  R95:M0[n-A95]  R96:M0[n-A96]  W97:M0[n-1279]
W98:M0[n-A98]  R99:M0[n-A99]  W100:M0[n-A100]  R101:M0[n-A101]  R103:M0[n-A103]  R106:M0[n-A106]
W107:M0[n-1791]  W108:M0[n-A108]  R109:M0[n-A109]  R110:M0[n-A110]  W111:M0[n-A111]  R112:M0[n-A112]
W113:M0[n-A113]  W116:M1[n-26626]  R117:M0[n-A117]  R118:M0[n-A118]  W119:M0[n-A119]  W122:M2[n-26626]
R123:M0[n-A123]  R124:M0[n-A124]  W125:M0[n-A125]
```

Values written (w_s) and tank outputs:

```
w5   (M0) = -X123′*k1 + X123′ + X124′*k1
w10  (M0) = -X3*k7 + X3 + X4*k7
w14  (M0) = X11*k11 + X9*k9
w39  (M0) = X11*k11*k29*k36 + X31*k31*k36 + X32*k32*k36 + X33*k33*k36 + X34*k34 + X9*k27*k36 + X9*k29*k36*k9
w41  (M0) = X11*k11*k29*k36 + X31*k31*k36 + X32*k32*k36 + X33*k33*k36 + X34*k34 + X38*k38 + X9*k27*k36 + X9*k29*k36*k9
w45  (M0) = X11*k11*k29*k36*k40 + X31*k31*k36*k40 + X32*k32*k36*k40 + X33*k33*k36*k40 + X34*k34*k40 + X38*k38*k40 + X38 + X42*k42 + X9*k27*k36*k40 + X9*k29*k36*k40*k9
w47  (M0) = X11*k11*k29*k36*k40*k44 + X31*k31*k36*k40*k44 + X32*k32*k36*k40*k44 + X33*k33*k36*k40*k44 + X34*k34*k40*k44 + X38*k38*k40*k44 + X38*k44 + X42*k42*k44 + X42 + X9*k27*k36*k40*k44 + X9*k29*k36*k40*k44*k9
w52  (M0) = X44 + X46*k46 + X48*k48
w53  (M0) = X44*k50 + X46*k46*k50 + X48*k48*k50 + X51
w55  (M0) = X44*k52 + X46*k46*k52 + X46
w62  (M0) = X50 + X54*k54 + X56*k56 + X58*k58
w63  (M0) = X50*k60 + X54*k54*k60 + X56*k56*k60 + X58*k58*k60 + X61
w66  (M0) = X50*k62 + X54*k54*k62 + X56*k56*k62 + X56
w67  (M0) = X50*k64 + X54*k54*k64 + X65
w71  (M0) = X64*k66 + X68*k68
w83  (M0) = X64*k66*k73*k80 + X64*k71*k80 + X68*k68*k73*k80 + X75*k75*k80 + X76*k76*k80 + X77*k77*k80 + X78*k78
w86  (M0) = X64*k66*k73*k80 + X64*k71*k80 + X68*k68*k73*k80 + X75*k75*k80 + X76*k76*k80 + X77*k77*k80 + X78*k78 + X79*k83
w90  (M0) = X64*k66*k73*k80*k85 + X64*k71*k80*k85 + X68*k68*k73*k80*k85 + X75*k75*k80*k85 + X76*k76*k80*k85 + X77*k77*k80*k85 + X78*k78*k85 + X79*k83*k85 + X79 + X87*k87
w92  (M0) = X64*k66*k73*k80*k85*k89 + X64*k71*k80*k85*k89 + X68*k68*k73*k80*k85*k89 + X75*k75*k80*k85*k89 + X76*k76*k80*k85*k89 + X77*k77*k80*k85*k89 + X78*k78*k85*k89 + X79*k83*k85*k89 + X79*k89 + X87*k87*k89 + X87
w97  (M0) = X89 + X91*k91 + X93*k93
w98  (M0) = X89*k95 + X91*k91*k95 + X93*k93*k95 + X96
w100 (M0) = X89*k97 + X91*k91*k97 + X91
w107 (M0) = X101*k101 + X103*k103 + X95 + X99*k99
w108 (M0) = X101*k101*k105 + X103*k103*k105 + X106 + X95*k105 + X99*k105*k99
w111 (M0) = X101*k101*k107 + X101 + X95*k107 + X99*k107*k99
w113 (M0) = X110 + X95*k109 + X99*k109*k99
w116 (M1) = 3*sendR/4
w119 (M0) = -X109*k115 + X109 + X112*k115
w122 (M2) = 3*sendL/4
w125 (M0) = -X117*k121 + X117 + X118*k121
retL      = X18*k18*k25 + X19*k19*k25 + X20*k20*k25 + X21*k21*k25 + X22*k22*k25 + X26*k26
retR      = X12*k14*k23 + X13*k13*k23 + X15*k15*k23 + X16*k16*k23 + X17*k17*k23 + X24*k24
```

Coefficients `k<s>` set by the parameters (laws of section 4; values per program in section 7; every other coefficient is a number in the program):

| k | set by (law) |
|---|---|
| k1 | modulated tap fraction (6.2) |
| k7 | modulated tap fraction (6.2) |
| k9 | Crossover (4.4) |
| k11 | Crossover (4.4) |
| k13 | Depth (4.5) |
| k14 | Depth (4.5) |
| k15 | Depth (4.5) |
| k16 | Depth (4.5) |
| k18 | Depth (4.5) |
| k19 | Depth (4.5) |
| k20 | Depth (4.5) |
| k21 | Depth (4.5) |
| k23 | Rvb Out (4.5) |
| k24 | RefLvl R (4.5) |
| k25 | Rvb Out (4.5) |
| k26 | RefLvl L (4.5) |
| k27 | Size with Link off, Mid Rt, Link (4.2) |
| k29 | Size with Link off, Mid Rt, Link (4.2); Low Rt (4.3) |
| k31 | fixed (0.6875) |
| k32 | fixed (0.125) |
| k33 | fixed (0.0625) |
| k34 | Rt HC (4.4) |
| k36 | Rt HC (4.4) |
| k38 | Diffusion (4.5) |
| k40 | Diffusion (4.5) |
| k42 | Diffusion (4.5) |
| k44 | Diffusion (4.5) |
| k46 | AP Sweep pair 1, first coefficient (6.1); rest gain: Def (4.5) |
| k48 | AP Sweep pair 2, first coefficient (6.1); rest gain: Def (4.5) |
| k50 | AP Sweep pair 2, second coefficient (6.1); rest gain: Def (4.5) |
| k52 | AP Sweep pair 1, second coefficient (6.1); rest gain: Def (4.5) |
| k54 | AP Sweep pair 3, first coefficient (6.1); rest gain: Def (4.5) |
| k56 | AP Sweep pair 4, first coefficient (6.1); rest gain: Def (4.5) |
| k58 | AP Sweep pair 5, first coefficient (6.1); rest gain: Def (4.5) |
| k60 | AP Sweep pair 5, second coefficient (6.1); rest gain: Def (4.5) |
| k62 | AP Sweep pair 4, second coefficient (6.1); rest gain: Def (4.5) |
| k64 | AP Sweep pair 3, second coefficient (6.1); rest gain: Def (4.5) |
| k66 | Crossover (4.4) |
| k68 | Crossover (4.4) |
| k71 | Size with Link off, Mid Rt, Link (4.2) |
| k73 | Size with Link off, Mid Rt, Link (4.2); Low Rt (4.3) |
| k75 | fixed (0.6875) |
| k76 | fixed (0.125) |
| k77 | fixed (0.0625) |
| k78 | Rt HC (4.4) |
| k80 | Rt HC (4.4) |
| k83 | Diffusion (4.5) |
| k85 | Diffusion (4.5) |
| k87 | Diffusion (4.5) |
| k89 | Diffusion (4.5) |
| k91 | AP Sweep pair 6, first coefficient (6.1); rest gain: Def (4.5) |
| k93 | AP Sweep pair 7, first coefficient (6.1); rest gain: Def (4.5) |
| k95 | AP Sweep pair 7, second coefficient (6.1); rest gain: Def (4.5) |
| k97 | AP Sweep pair 6, second coefficient (6.1); rest gain: Def (4.5) |
| k99 | AP Sweep pair 8, first coefficient (6.1); rest gain: Def (4.5) |
| k101 | AP Sweep pair 9, first coefficient (6.1); rest gain: Def (4.5) |
| k103 | AP Sweep pair 10, first coefficient (6.1); rest gain: Def (4.5) |
| k105 | AP Sweep pair 10, second coefficient (6.1); rest gain: Def (4.5) |
| k107 | AP Sweep pair 9, second coefficient (6.1); rest gain: Def (4.5) |
| k109 | AP Sweep pair 8, second coefficient (6.1); rest gain: Def (4.5) |
| k115 | modulated tap fraction (6.2) |
| k121 | modulated tap fraction (6.2) |

## 5.4 Inverse tank — exact program

### Structure

The Inverse program has no recirculating tank. It builds a reverse / gated reverb from delay taps:

* **Input**: as in the Plate and Chamber, each input passes input diffusers (allpasses) into its input line
  (M1 right, M2 left); taps of these lines (Pre Delay, RefDly) feed the rest of the program.
* **Envelope path**: four taps of each input line are weighted by the Shape coefficients (four per side) and summed;
  the sum is diffused by two allpasses and written into a long delay region of M0.
* **Output**: each output is a fixed sum of 16 taps of that region with fixed gains (±3/8, ±7/16, ±1/2) plus two
  direct Shape-weighted terms. The tap positions scale with Duration, so the response is a block of dense
  reflections whose envelope (rising, falling or flat) is set by Shape and the two Slope parameters.
* There is no AP Sweep sequencer and no modulation; the program is static.

### Program

Delay memory: three circular buffers of 65 536 samples, **M0** (tank), **M1** (right input line), **M2** (left input line). Per sample n, execute the operations below in order. `R s: X_s = Mb[n − A]` reads (offset A samples back), `W s: Mb[n − A] = w_s` writes the value given in the equations. Offsets written as `A<s>` vary with the program (laws of section 4, values in section 7); all others are fixed. X′ denotes the value the same read produced in the **previous** sample (partial sums that the program carries across the sample boundary). `sendL/sendR` = reverb send into the tank, `retL/retR` = tank output.

```
R0:M0[n-A0]  R1:M0[n-A1]  W3:M0[n-A3]  W5:M0[n-A5]  R6:M0[n-A6]  R7:M1[n-A7]
R8:M1[n-A8]  W9:M0[n-A9]  W11:M0[n-A11]  R13:M1[n-A13]  R14:M1[n-A14]  R15:M0[n-A15]
R19:M0[n-A19]  W20:M0[n-A20]  W22:M0[n-A22]  R24:M0[n-6]  R25:M0[n-A25]  R26:M0[n-2]
W27:M0[n-5]  R28:M0[n-6]  W31:M0[n-1]  R34:M0[n-A34]  R35:M2[n-32769]  W37:M0[n-A37]
R38:M0[n-A38]  R39:M2[n-A39]  R40:M0[n-A40]  W41:M0[n-A41]  R42:M0[n-A42]  W43:M0[n-A43]
W45:M0[n-A45]  R46:M2[n-A46]  R47:M2[n-A47]  R49:M0[n-A49]  R52:M0[n-A52]  R53:M0[n-A53]
W54:M0[n-A54]  R55:M0[n-A55]  W56:M0[n-A56]  R57:M0[n-A57]  R58:M0[n-A58]  R59:M0[n-A59]
R60:M0[n-A60]  R61:M0[n-A61]  R62:M0[n-A62]  R63:M0[n-A63]  R64:M0[n-A64]  R65:M0[n-A65]
R66:M0[n-A66]  R67:M0[n-A67]  R68:M0[n-A68]  R69:M1[n-A69]  R70:M1[n-A70]  R71:M0[n-A71]
R73:M0[n-A73]  R74:M0[n-A74]  R75:M0[n-A75]  R76:M0[n-A76]  R77:M0[n-A77]  R78:M0[n-A78]
R79:M0[n-A79]  R80:M0[n-A80]  R81:M0[n-A81]  R82:M0[n-A82]  R83:M0[n-A83]  R84:M0[n-A84]
R85:M0[n-A85]  R86:M2[n-A86]  R87:M2[n-A87]  R88:M0[n-A88]  W91:M0[n-20]  R92:M0[n-A92]
W95:M0[n-A95]  R96:M0[n-A96]  W99:M0[n-A99]  R100:M0[n-A100]  W101:M1[n-0]  W103:M0[n-A103]
R104:M0[n-A104]  W107:M0[n-A107]  R108:M0[n-A108]  R109:M0[n-4]  R110:M0[n-A110]  W111:M0[n-A111]
W113:M2[n-32768]  R114:M0[n-8]  W117:M0[n-7]  R119:M0[n-8]  W122:M0[n-3]
```

Values written (w_s) and tank outputs:

```
w3   (M0) = -X0/2 + X109′*k117*k123 + X114′*k121 + X114′ + X119′*k119*k123
w5   (M0) = 3*X0/4 + X109′*k117*k123/2 + X114′*k121/2 + X114′/2 + X119′*k119*k123/2
w9   (M0) = X1*k4 + X1 - X6/2
w11  (M0) = X1*k4/2 + X1/2 + 3*X6/4
w20  (M0) = X13*k13*k17 + X14*k14*k17 - X15/2 + X7*k11*k17 + X8*k12*k17
w22  (M0) = X13*k13*k17/2 + X14*k14*k17/2 + 3*X15/4 + X7*k11*k17/2 + X8*k12*k17/2
w27  (M0) = X19*k22 + X24*k24
w31  (M0) = X26*k26 + X28*k28
w37  (M0) = X24*k30 + X24 + X26*k26*k32 + X28*k28*k32 - X34/2
w41  (M0) = X24*k30/2 + X24/2 + X26*k26*k32/2 + X28*k28*k32/2 + 3*X34/4
w43  (M0) = X38*k38 + X38 - X40/2
w45  (M0) = X38*k38/2 + X38/2 + 3*X40/4
w54  (M0) = X35*k48*k51 + X39*k45*k51 + X46*k46*k51 + X47*k47*k51 - X49/2
w56  (M0) = X35*k48*k51/2 + X39*k45*k51/2 + X46*k46*k51/2 + X47*k47*k51/2 + 3*X49/4
w91  (M0) = X88*k88 + sendR
w95  (M0) = X88*k88*k90 + X88 + X92*k92 + k90*sendR
w99  (M0) = X88*k88*k90*k94 + X88*k94 + X92*k92*k94 + X92 + X96*k96 + k90*k94*sendR
w101 (M1) = X88*k88*k90*k94*k98 + X88*k94*k98 + X92*k92*k94*k98 + X92*k98 + X96*k96*k98 + X96 + k90*k94*k98*sendR
w103 (M0) = X100*k100 + sendL
w107 (M0) = X100*k100*k102 + X100 + X104*k104 + k102*sendL
w111 (M0) = X100*k100*k102*k106 + X100*k106 + X104*k104*k106 + X104 + X108*k108 + k102*k106*sendL
w113 (M2) = X100*k100*k102*k106*k110 + X100*k106*k110 + X104*k104*k106*k110 + X104*k110 + X108*k108*k110 + X108 + k102*k106*k110*sendL
w117 (M0) = X110*k112 + X114*k114
w122 (M0) = X109*k117 + X119*k119
retL      = 3*X52/8 + X71/2 + X73/2 - 3*X74/8 + X75/2 - 7*X76/16 - 3*X77/8 + X78/2 - 7*X79/16 + 3*X80/8 - X81/2 + 3*X82/8 + X83/2 - 3*X84/8 - 3*X85/8 + X86*k86 + X87*k87
retR      = -3*X53/8 + X55/2 - X57/2 - 7*X58/16 - 3*X59/8 + X60/2 - X61/2 - 7*X62/16 + 3*X63/8 + X64/2 + 3*X65/8 + X66/2 - 3*X67/8 - 3*X68/8 + X69*k69 + X70*k70
```

Coefficients `k<s>` set by the parameters (laws of section 4; values per program in section 7; every other coefficient is a number in the program):

| k | set by (law) |
|---|---|
| k4 | Mid Slope (4.9) |
| k11 | Shape (4.7) |
| k12 | Shape (4.7) |
| k13 | Shape (4.7) |
| k14 | Shape (4.7) |
| k17 | Mid Slope (4.9) |
| k22 | Rt HC (4.4) |
| k24 | Rt HC (4.4) |
| k26 | Crossover (4.4) |
| k28 | Crossover (4.4) |
| k30 | Mid Slope (4.9) |
| k32 | Low Slope, Mid Slope (4.9) |
| k38 | Mid Slope (4.9) |
| k45 | Shape (4.7) |
| k46 | Shape (4.7) |
| k47 | Shape (4.7) |
| k48 | Shape (4.7) |
| k51 | Mid Slope (4.9) |
| k69 | RefLvl R (4.5) |
| k70 | Shape (4.7) |
| k86 | RefLvl L (4.5) |
| k87 | Shape (4.7) |
| k88 | Diffusion (4.5) |
| k90 | Diffusion (4.5) |
| k92 | Diffusion (4.5) |
| k94 | Diffusion (4.5) |
| k96 | Diffusion (4.5) |
| k98 | Diffusion (4.5) |
| k100 | Diffusion (4.5) |
| k102 | Diffusion (4.5) |
| k104 | Diffusion (4.5) |
| k106 | Diffusion (4.5) |
| k108 | Diffusion (4.5) |
| k110 | Diffusion (4.5) |
| k112 | Rt HC (4.4) |
| k114 | Rt HC (4.4) |
| k117 | Crossover (4.4) |
| k119 | Crossover (4.4) |
| k121 | Mid Slope (4.9) |
| k123 | Low Slope, Mid Slope (4.9) |

## 6. Modulation

Two slow processes change tank coefficients and offsets while the reverb runs (6.1, 6.2); both are described as
algorithms in sample time. In addition, patched parameters follow their modulation sources (6.3), the Chorus+Rvb
voices are modulated (3.1), and every parameter change reaches the signal as described in 6.4. Programs with
AP Sweep = 0, no modulated taps and no active modulator (marked in section 7) are completely static.

### 6.1 AP Sweep: sequential allpass-gain sweep (Plate, Chamber, Concert Hall)

The tank allpasses are grouped in pairs (a, b): the two coefficients of one allpass. In sweep order:

| tank | pairs (first / second coefficient) |
|---|---|
| Chamber | k0/k2, k6/k8, k15/k17, k34/k36, k40/k42, k49/k51 |
| Plate (also Chorus+Rvb) | k123/k125, k3/k5, k11/k13, k29/k31, k36/k38, k44/k46 |
| Concert Hall | k46/k52, k48/k50, k54/k64, k56/k62, k58/k60, k91/k97, k93/k95, k99/k109, k101/k107, k103/k105 |

Each pair has a rest gain g0 (Plate 0.5, Chamber 0.5625, Concert Hall set by Def, 4.5) and a stock sign (the signs
of section 7). The gain is handled as an 8-bit code `c = 255 − 256·|g|` (g = (255 − c)/256), the second coefficient
of the pair b directly, the first a derived from it: for an allpass whose first coefficient is −g/2 doubled (Plate,
Chamber, and Concert Hall pairs 2, 5, 7 and 10) `c_a = (c_b − [bit1(c_b) = 0]) >> 1` (bit1 = the bit of value 2),
which after doubling modulo 1 equals g within 1/256. Otherwise `c_a = c_b`. So every pair always stays a valid
allpass.

**Timing.** The AP Sweep parameter v (0 … 50, displayed as 2v %) sets a base period byte

    x = ⌊v·32768/50⌋,   B = max(1, 255 − ⌊178·x/32768⌋)

While AP Sweep runs, the period P is not constant: every 1.10 s (52 800 samples; 1 090 ticks of an ≈ 1 ms timer) it
moves to the next entry of an 11-step multiplier cycle

    T = (128, 192, 179, 177, 94, 162, 147, 64, 124, 109, 132)      (k = 0, 1, …, 10, 0, …)
    P_k = 255 if ⌊B·T_k/128⌋ ≥ 255, else ⌊B·T_k/128⌋ + 1

so the sweep speed changes stepwise (by up to a factor 3 at AP Sweep Full) within each 12.1 s cycle. One step (one call
of the state routine below, including the LOAD and END calls, i.e. 2·256·g0 + 2 calls per visit) takes
**0.130·P_k samples** (one constant for every AP Sweep setting, within ±0.7 %). Over a whole cycle the mean step time is
0.130·P̄ with P̄ = 11 / Σ(1/P_k) (harmonic mean; table in section 7). The phase of the multiplier cycle relative to
the program load is not defined (free-running timer); a simulation starts it at k = 0.

```
for each pair in list order, cyclically:
    repeat: c_b += 1 (gain moves 1/256 toward 0)      until c_b = 255  (g = 0)
    flip the sign of both coefficients of the pair
    repeat: c_b −= 1                                   until c_b = c_rest (|g| = g0, sign now inverted)
    (AP Sweep = 0: the process stops here after the pair it was visiting)
```

One visit takes 2·256·g0 + 2 steps ≈ (2·256·g0 + 2)·0.130·P̄ samples on average; e.g. Vox Plate (g0 = 0.5, AP Sweep 10 %,
P̄ = 215.0): 258 steps ≈ 7 210 samples = 0.15 s per pair, 0.90 s for all six. Each pair's sign alternates from one
visit to the next.

**New rest gains** (a Def change, Concert Hall) take effect through the sequencer's 1/256 steps: every pair outside
the current visit moves its code by one step per sequencer step toward its new rest code, also while AP Sweep is 0
(the step timer runs until every pair has arrived); the pair being visited ends its visit at its new rest gain. The
coefficients never jump.

**At rest (AP Sweep = 0)** all pairs hold their rest gain. After a program change all pairs fade to zero
and back, and the signs of the first one or two pairs of the list are then not defined (they may return inverted).
This changes the fine structure of the impulse response but not its decay or spectrum. A program load sets all
pairs to their rest gain with the stock signs, as tabulated.

### 6.2 Concert Hall modulated taps

Four interpolated taps read the four modulation lines of the Hall tank (section 5.3). Tap j reads two adjacent samples
at offsets A_a, A_b (low and high byte; the high byte, i.e. the 256-sample window, never changes) and interpolates
with the fraction coefficient k (`y = (1 − k)·M0[n − A_a] + k·M0[n − A_b]`, k = (255 − c)/256 from an 8-bit code c):

| tap | fraction | reads | window (offsets) |
|---|---|---|---|
| 1 | k115 | A109, A112 | 256 … 511 |
| 2 | k121 | A117, A118 | 768 … 1023 |
| 3 | k1   | A123, A124 | 1280 … 1535 |
| 4 | k7   | A3, A4     | 1792 … 2047 |

**Motion** (exact rule): one step every T = 0.2575·(D + 4.1) samples, D = divider set by Chorus (4.5;
Chorus 0: no steps), within 1 % over Chorus 1 … 10 (tap speed 728/(D + 4.1) samples/s);
equivalently, a step every D + 4.1 periods of a clock running at ≈ 4·fs.

```
state: c_j (8-bit fraction code), dir_j = ±1, s_j = ±2 (address step), countdown N, next
on every step:
    N −= 1
    if N >= 0:  c_j = (c_j + dir_j) mod 256  for all four taps
    else:
        N = 64;  j = next;  next = (next + 1) mod 4  (when next wraps to 0: N = 63)
        dir_j = −dir_j
        if c_j == 255:                       # fraction at its end: move the second read
            x = lowbyte(A_b,j)
            if s_j > 0: lowbyte(A_b,j) = x + s_j
            elif (x + s_j) mod 256 == 255: s_j = −s_j     # window edge: reverse the address stepping
            else: lowbyte(A_b,j) = x + s_j
        else:                                # otherwise move the first read
            x = lowbyte(A_a,j)
            if s_j < 0: lowbyte(A_a,j) = x + s_j
            elif (x + s_j) mod 256 == 0: s_j = −s_j
            else: lowbyte(A_a,j) = x + s_j
start state:
    c    = (0xC0, 0x7F, 0x40, 0xFF)      (= 0x3F, 0x7F, 0xBF, 0xFF, complemented where dir = +1)
    dir  = (+1, −1, +1, −1),   s = (+2, −2, +2, −2),   N = 63,   next = 0
    read pair of each tap: lowbyte(A_b) odd, A_a = A_b + 1 (inside its 256-sample window)
```

With this start state each tap is visited exactly when its fraction code has reached an end (c = 255 or 0). The read
that moves therefore carries a weight of at most 1/256, so the motion is click-free, and the two reads of every tap
stay adjacent (|A_a − A_b| = 1). All four taps sweep their delay as triangle waves across the whole 256-sample window
at a constant speed of 1/(256·T) samples per sample, i.e. 728/(D + 4.1) samples per second (Concert Hall, D = 1:
143 samples/s; ConcertHall2, D = 7: 66 samples/s), with staggered phases. The phase of the motion depends on the time since the program
was loaded; the tables of section 7 give one snapshot, and a simulation places each tap's start read pair at the
tabulated window position.

### 6.3 Modulation sources (patch sources)

Patches (4.10) take their input from one of these sources; every source delivers 0 … 127:

* **MIDI** (16.1): controllers 1–31 and 33–119, P Bend, A Touch, Velocity, Last Note, Low Note, High Note,
  Clk Cmnds
* **LFO** (the waveform chosen by LFO Shape) and the six waveforms **Sine**, **Cosine**, **Square**, **Sawtooth**,
  **Pulse**, **Triangle**
* **Env L**, **Env R** (left / right input envelope), **AR Env** (attack/release envelope), **Latch**
* **Sw 1**, **Sw 2** (time switches), **Sw 1&2**
* **Mono Lvl** (input level, L+R), **Left Lvl**, **Right Lvl** (left / right input level)
* **Footpedal**, **Foot Sw 1**, **Foot Sw 2**, **ADJUST** knob, **Tempo**
* **On**, **Off**

The six waveforms and the LFO source come from one LFO (one phase, one rate). Each program stores the modulator
settings ("Mod" page) and a trigger block (7.1). MIDI controllers, Footpedal, Foot switches and Tempo take the value
supplied from outside; until one is supplied they have no value, and a patch from them leaves its destination at
the stored value (16.1). On = 127, Off = 0.

**Control rate.** All modulators, the input-level sources and the patch evaluation run once per **20 ms control step
(50 Hz)**, at both sample rates. In each step the generators are updated first, then every patch is re-evaluated
(curve law in 4.10) and the results go to the destinations (6.4). When a step is late, the next one updates every
generator twice (one catch-up step: the second update of each generator is noted in its rule below), so the average
rate stays 50 Hz. A simulation that runs every step on time never takes the catch-up path.

**Trigger block.** The AR envelope and the two time switches each have a trigger source s and a trigger level T
(7.1, "trigger"; "none" = no source). In a control step the generator is **triggered** when s has a value and
s ≥ T, or when there is no source and T = 0.

**LFO** (settings Rate 0–2500, Shape 0–5, PWidth 0–98, Depth 0–100):

```
inc   = floor((6710·Rate + 128) / 256)          # phase increment per step; f = inc·50 / 2^17 Hz ≈ Rate / 100 Hz
pw    = floor(512·(PWidth + 1) / 100)           # pulse threshold, 5 … 506
dm    = round(2.56·Depth)                       # depth multiplier, 0 … 256
per step:
  i  = ((φ + 128) >> 8) mod 512                 # φ: unsigned phase, one cycle = 2^17
  Sine     = W[i]            Cosine   = W[(i + 128) mod 512]
  Square   = 127 if i < 256 else 0              Sawtooth = i >> 2
  Pulse    = 127 if i < pw  else 0              Triangle = i >> 1 if i < 256 else 255 − (i >> 1)
  LFO      = floor(w·dm / 256),  w = the waveform selected by Shape
              (0 Sine, 1 Cosine, 2 Square, 3 Sawtooth, 4 Pulse, 5 Triangle)
  φ += inc                                      # after the outputs: the outputs use the phase before the step
```

W is a 512-entry sine table on 0 … 127: W[k] for k = 0 … 127 is

    64, 65, 66, 66, 67, 68, 69, 69, 70, 71, 72, 72, 73, 74, 75, 76, 76, 77, 78, 79, 79, 80, 81, 82, 82, 83, 84, 84, 85, 86, 87, 87, 88, 89, 90, 90, 91, 92, 92, 93, 94, 94, 95, 96, 96, 97, 98, 98, 99, 100, 100, 101, 102, 102, 103, 103, 104, 105, 105, 106, 106, 107, 107, 108, 109, 109, 110, 110, 111, 111, 112, 112, 113, 113, 114, 114, 115, 115, 116, 116, 116, 117, 117, 118, 118, 118, 119, 119, 120, 120, 120, 121, 121, 121, 122, 122, 122, 122, 123, 123, 123, 124, 124, 124, 124, 125, 125, 125, 125, 125, 125, 126, 126, 126, 126, 126, 126, 126, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127

and W[255 − k] = W[k], W[k + 256] = 127 − W[k] (W ≈ round(63.5 + 63.35·sin(2π(k + ½)/512)), 16 entries differ by 1).
Example: JetChamber (Rate 8 → inc 210 → 0.0801 Hz).

**Input levels** (Mono Lvl, Left Lvl, Right Lvl): per step, p_L and p_R are the input peak magnitudes since the previous step, in
1/128 of full scale (0 … 127); Left Lvl = V[p_L], Right Lvl = V[p_R], Mono Lvl = floor((Left Lvl + Right Lvl) / 2), with
the 128-entry level map V (piecewise-linear log curve, ≈ 4 units per dB above p = 16):

    0, 3, 6, 10, 14, 18, 22, 25, 28, 31, 34, 37, 40, 43, 46, 50, 52, 53, 55, 56, 58, 59, 61, 62, 64, 65, 67, 68, 70, 71, 73, 75, 76, 77, 78, 78, 79, 80, 81, 81, 82, 83, 84, 84, 85, 86, 87, 87, 88, 89, 90, 90, 91, 92, 93, 93, 94, 95, 96, 96, 97, 98, 99, 100, 100, 101, 101, 102, 102, 102, 103, 103, 104, 104, 105, 105, 105, 106, 106, 107, 107, 108, 108, 108, 109, 109, 110, 110, 111, 111, 111, 112, 112, 113, 113, 114, 114, 114, 115, 115, 116, 116, 117, 117, 117, 118, 118, 119, 119, 120, 120, 120, 121, 121, 122, 122, 123, 123, 123, 124, 124, 125, 125, 126, 126, 126, 127, 127

**AR Env** (settings Attack N_A, Release N_R, Mode): a 16-bit state e (0 … 65535), output e >> 9 (0 … 127). The step
sizes are a_A = floor(65535 / N_A), a_R = floor(65535 / N_R) (65535 when the setting is 0; settings above 500 count as
500), so attack and release take N or N + 1 steps of 20 ms (N = 500 → 10.02 s). Per step, with "triggered" from the
trigger block:

```
Mode 2 (Retrigger): triggered and not at the top → e = min(e + a_A, 65535)        (attack, also from a partial release)
                    not triggered and e > 0      → e = max(e − a_R, 0)            (release starts at once)
Mode 1 (One Shot):  at rest and triggered → attack; the attack always completes; then hold at 65535 while triggered;
                    when no longer triggered, a full release (not interrupted by a new trigger)
Mode 3 (Repeat):    at rest and triggered → attack; at the top → full release at once; at 0, attack again if triggered
                    (triangle wave when N_A = N_R)
Mode 0: off
```

Example: BowedChords (AR Attack 63, Mode 2, trigger Mono Lvl ≥ 5, AR → input level −40 … 0 dB) rises linearly in
dB over 63 steps = 1.26 s from the onset of a signal. Modes 1 and 3: 16.2.

**Env L / Env R** (input envelope followers, setting Release 0–500). Each follower has a 16-bit state e (0 … 65535).
Its input is **Left Lvl** (Env L) or **Right Lvl** (Env R), the 0 … 127 level-map value above, not the raw peak. The
release step is

    r = ⌊65535 / N⌋     N = min(Release, 500);  r = 65535 when Release = 0

Per control step, with x = the input level of this step and t = 516·x (x = 127 → 65532):

```
if e ≤ t:            e = t                        # attack: instant, to the input level
elif e − t < r:      e = t                        # release would pass the input level: stop there
else:                e = e − r                    # linear release
                     (catch-up step: a second  e = e − r, or e = 0 if r > e; this one is not floored at t)
output = e >> 9                                   # 0 … 127; equals x during the attack
```

So the output follows the input level instantly upward and falls linearly from 127 to 0 in N control steps
(N·20 ms; Release 500 = 10 s), never below the current input level. **Exact.** Used by BigBoomRoom (Env L Release 38,
Env R Release 28, Env L/R → post-delays) and Drum Gate (Env L → Rvb Width, ADJUST → Env L Release).

**Sw 1 / Sw 2** (time switches; settings Mode, Rate 0–2500, PWidth 0–98). Each switch has a 32-bit phase φ (one
cycle = 2^17), a ramp value q (16-bit, 0 … 0x7F00), a running flag, a restart flag (set after a program load) and these
settings:

```
inc  = ⌊(6710·Rate + 128) / 256⌋                    # the LFO rate law: f = inc·50/2^17 Hz ≈ Rate/100 Hz
pw   = ⌊512·(PWidth + 1) / 100⌋                     # the LFO pulse law
up   = ⌊127·inc / pw⌋ ,  down = ⌊127·inc / (512 − pw)⌋   # ramp steps, computed when Mode = 2 (Ramp)
Mode 0 = off, 1 = Switch (pulse), 2 = Ramp
```

Per control step (Mode 1 or 2):

```
if triggered: running = 1
if running:
    if restart:  φ = 0, q = 0, restart = 0, rising = 1, i = 0
    else:        i = (φ >> 8) & 0x1FF                         # 0 … 511
    if Mode == 1:  out = 127 if i < pw else 0
    else:                                                     # Ramp
        out = q >> 8                                          # value before this step's update
        if i < pw:  (if not rising: rising = 1, q = 0);  q = min(q + up, 0x7F00)
        else:       (if rising: rising = 0, q = 0x7F00);  q = q − down, and q = 0 if that goes below 0
    φ += inc   (twice in a catch-up step)
    if φ ≥ 2^17:  φ −= 2^17;  running = 0;  if not triggered: out = 0
    write out to the source
else:
    restart = 1                                               # output not written: it keeps its last value
```

* A trigger starts a full cycle at phase 0. The cycle always runs to its end, and further cycles follow while the
  trigger is still present at the end of a cycle. When the trigger has gone, the cycle that ends writes 0 and the
  switch stops.
* **Idle output:** 0 after a program load (Mode 1/2); Mode 0 marks the source as having no value (4.10). After a
  completed cycle the output stays at the last written value, which is 0.
* **Sw 1&2** outputs the value of one of the two switches, starting with Sw 1. When the selected switch completes a
  cycle (its running flag clears in that step), the selection passes to the other switch, so Sw 1&2 plays one cycle
  of Sw 1, then one cycle of Sw 2, and so on (at most one hand-over per step; the step in which a cycle ends still
  outputs the switch that ended it).

**Exact.** Piano Hall: Sw 1 Mode 2 (Ramp), Rate 26, PWidth 49, trigger Mono Lvl ≥ 75 → inc = 681, pw = 256,
up = down = 337: a 192.5-step cycle (3.85 s) that rises 0 → 127 in 97 steps (1.93 s) and falls back in 97 steps,
driving Rvb Width through its patch curve.

**Latch** (no program of section 7 uses it): the exact state machine of the original is given in 16.5.

### 6.4 Parameter changes while running

Every patch destination is re-evaluated in the 20 ms control step in which its source changes and applied by the
parameter's own handler in the **next** control step (one step of latency). An edit
of a parameter, or a tempo change, goes to the parameter's handler directly. From the handler, the new values reach
the signal as follows.

* **Tank coefficients** move by a linear ramp. When a ramp starts, Δ = |target − current| in 16-bit units
  (1.0 = 32 768) sets the rate 2·(⌊Δ/128⌋ | 1) per ramp step (`|` = bitwise or); each ramp step moves the coefficient
  by that rate, clamped at the target, so a ramp takes about 64 ramp steps whatever its size. At each step the
  coefficient takes the ramp value at its own resolution (4-bit coefficients — those set by Shape / Attack, Rvb Out,
  RefLvl and Depth — in 1/16, all others in 1/256, truncated toward zero). One ramp step takes
  about 2.5 samples in Plate, Chamber, Inverse and Chorus+Rvb, about 54 samples in Concert Hall
  (whose modulated-tap processing lengthens every step), so a ramp takes ≈ 3 ms and ≈ 70 ms. These durations are
  accurate to ±15 %. The AP Sweep pairs (6.1) and the modulated taps (6.2) are written by their
  own processes, not ramped.
* **Tank offsets** are not ramped: the coefficient of the step that uses the offset fades out in steps of one
  coefficient LSB, the offset changes, and the coefficient fades back in. Changing the offset at once instead is
  heard as a jump of the tap (9.5).
* **Input/output-section coefficients** (send matrix, gains, mixes, widths, output low-pass, voice levels and
  feedback, chorus voice levels and feedback) glide toward their new value through one shared smoothing with the
  same rate for every coefficient: per sample `c += κ·(target − c)` with κ = 1 − exp(−1/960) (time constant ≈ 20 ms
  at 48 kHz).
* **Delay lengths** (post-delays, delay voices, chorus voice targets) glide or cross-fade by the rule of 4.13. The
  chorus depth targets A_j are smoothed by c38 (3.1); a new base rate R_j takes effect at the voice's next rate update.
* **Results of exactly zero.** A voice level of exactly 80 (off), a feedback of exactly 100 (none), a master at 0, a
  pan that gives one side exactly 0, or a chorus rate that rounds to 0 is not written: that coefficient is taken out
  of the smoothing instead and keeps the value it has. After a program load all coefficients start at 0, so a stored
  "off" is silent, while a live edit that lands exactly on "off" leaves the coefficient where it was.
* **Clear** freezes and releases the voice levels and feedback gains (4.13).

## 7. Programs

The stored data and the resulting values of the 34 programs of Part I are in the program data file (17.2), in the same
form as for the other 166 programs: the stored values (`stored`), the parameters in effect (`parameters`: stored →
effective value after the patch at the stored knob position; Mid Rt also as RT[v] in `rt_s`), the modulation state
(`modulation`), the tank coefficients k and program-dependent offsets A (`tank`; symbols of section 5) and the
input/output-section values (`io`; section 3 or 3.1, keys as in 17.2). The key `part1_section` names the group
below. Everything not given there is identical for all programs of a tank and appears as a number in the tank program.

Notes on the resulting values:
* The AP Sweep pairs are given at their rest gain with the stock signs (6.1). Programs marked "sweeping" keep
  modulating them; their values are the rest values.
* Concert Hall modulated tap coefficients (k1, k7, k115, k121) and their read offsets (A3, A4, A109, A112, A117, A118,
  A123, A124) are a snapshot of a continuously moving state (6.2).
* Patches from sources other than the ADJUST knob (LFO, envelopes, switches, input levels; 6.3) are evaluated at
  source value 0 (input silent, envelopes at rest); the affected values are marked "(patch at rest)". Exception:
  BowedChords is given with its AR envelope at the top (input level 0 dB), see 7.5. The stored patches
  (`stored.patches`) list every patch with its curve.
* Chamber, Plate, Concert Hall and Inverse programs share the input/output section of 3; the Chorus+Rvb program
  (Deep Space) has its own (3.1) and runs the Plate tank.
* The voice and post-delay lengths are given in samples at 48 kHz.

### 7.1 Stored program data

Key `stored` of each program record (17.2): algorithm, tempo (BPM), ADJUST knob position and label, the stored value of every parameter by edit
page (`*` = tempo-locked, value = 32·hi + lo, 4.11), the parameters on no edit page ("(no page)": stored, not
editable), the Mod page, the trigger block (6.3) and the patches with page-qualified destinations and their curves
[x:y …] (4.10). An ADJUST patch replaces its destination at the stored knob position (`parameters` gives the effective
values). The delay voices V1 … V4 of section 3 are Voice1 … Voice4 of the Levels, DelayTime, Feedback and Panning
pages; the chorus voices of 3.1 are Voice1 … Voice6 there and V1 … V6 of the Chorus page.

### 7.2 Chamber programs: Small Room, Living Room, Brick Kick, Large Room, SnareChamber, Tiled Room, Rich Chamber

Values: program data (17.2), records with `part1_section` = "7.2".

### 7.3 Chamber programs: Vox Chamber, Locker Room, Wide Chamber, BigBoomRoom, JetChamber, EnvChamber

Values: program data (17.2), records with `part1_section` = "7.3".

### 7.4 Chamber program BowedChords

Values: program data (17.2), records with `part1_section` = "7.4".

### 7.5 BowedChords: complete description

BowedChords is a Chamber program whose input level follows an attack/release envelope: the reverb and the delay
voices **swell in** over 1.26 s while the input is above a trigger level and fall back to −40 dB within 0.18 s when it
drops below it. Everything needed to reproduce it is in 7.1, 7.4 and this section.

**Reproduction recipe.**

1. Tank: the Chamber program of 5.1 with the coefficients k and offsets A of 7.4. AP Sweep is stopped: all
   six allpass pairs at their rest gain with the stock signs, as tabulated.
2. Input/output section: section 3 with the values of 7.4, except the send matrix, which is multiplied
   by the time-varying input gain g of step 3: aLL = aRR = 0.999969·g, aLR = aRL = 0. Because the voice lines take
   their input from mL and mR (section 3), g scales the reverb send and the delay voices alike.
3. Input gain g from the envelope rule below, run at the 50 Hz control rate, with a one-step latency and a one-pole
   smoothing per sample.
4. Start (program loaded, silence): e = 0, g = 0.01 (−40 dB), all delay memories zero.

All other values are static. The tables give the state with the envelope **at the top** (g = 1), unlike the
"(patch at rest)" convention of the other programs; at rest the send matrix is 0.01 times the tabulated values.

**Input gain (exact rule; the latency and smoothing are approximate).** Per sample n and per 20 ms
control step (960 samples at 48 kHz; the step phase relative to the input is arbitrary):

```
every sample:
    mL_peak = max(mL_peak, |inL|) ;  mR_peak = max(mR_peak, |inR|)
    g += c·(g_target − g)                          # c = κ = 1 − exp(−1/960) (6.4; time constant ≈ 20 ms)
    aLL = aRR = 0.999969·g                         # send matrix of section 3 in this sample
at the end of every control step:
    g_target = g_next                              # the value computed at the end of the previous step
    pL = min(127, floor(128·mL_peak)) ;  pR = min(127, floor(128·mR_peak)) ;  mL_peak = mR_peak = 0
    Mono = floor((V[pL] + V[pR]) / 2)              # level map V of 6.3 (input-level sources)
    if Mono >= 5:  e = min(e + 1040, 65535)        # attack:  floor(65535/63), AR Attack 63
    else:          e = max(e − 7281, 0)            # release: floor(65535/9),  AR Release 9 (Mode 2, Retrigger)
    x       = e >> 9                               # envelope output 0 … 127
    InLvl   = 160 if x >= 127 else 120 + floor(20641·x / 65536)    # patch curve [0:120 127:160] (4.10)
    g_next  = 10^((InLvl − 160)/20)                # input-level law (4.5): −40 … 0 dB in 1 dB steps
```

The rule follows the general AR envelope (6.3) and patch law (4.10). The two patches InLvl L ← AR Env and
InLvl R ← AR Env are identical, so one gain serves both inputs. The latency (one control step) and the smoothing
time constant (≈ 20 ms) are approximate; combinations with the same total delay of ≈ 40 ms (e.g. 30 ms
latency, 10 ms smoothing) are equivalent.

Consequences:
* A sustained input reaches full level after 63 steps = **1.26 s**, rising almost linearly in dB (≈ 0.64 dB per step,
  31.7 dB/s).
* After the input stops, the gain falls in **9 steps (0.18 s)** to −40 dB: −5 dB after the first step, then ≈ 4.5 dB
  per step.
* A new note re-attacks from wherever the release has got to.
* In silence, and for any input below the trigger, the program still passes the input at **−40 dB**. It is not muted.
* Trigger level: Mono ≥ 5 ⇔ V[p] ≥ 5 ⇔ p ≥ 2 when both inputs carry the signal, i.e. an input peak of 2/128 of full
  scale (−36.1 dB). A signal on one input only needs p ≥ 3 (V = 10, −32.6 dB).
* The gain starts to move 40–60 ms after the input starts or stops.

**Static part.**
* Tank: a large chamber (Size 90 of 144, Link 1, so the decay follows Size: ≈ 0.66 × RT[41] ≈ 1.54 s; Mid Rt 41,
  Low Rt 6, Crossover 8, Rt HC 25), dense input diffusion (Diffusion 86), and Shape 90 with Spread 255 (maximum) for a late, wide build-up.
  No pre-delay; early reflections RefDly 7 / 8 (14 / 16 ms) at RefLvl 10/16; no output low-pass (High Cut 121).
* Direct voice taps: V1 (left line, 0 ms: no delay) to the left and V3 (right line, 5 ms) to the right, both at
  0 dB and without feedback.
* Echo voices: V2 (left line) to the left and V4 (right line) to the right, at ¼ beat of the program tempo (44 BPM:
  16 416 samples = 342 ms, tempo-locked, 4.11), at −6 dB with feedback −0.456, so each echo is 6.8 dB lower and
  inverted.
* The ADJUST knob drives **FX Mix** (patch [0:0 100:100], so FX Mix = knob position; the stored FX Mix 43 has no
  effect). At the stored knob 92 the voices are mixed at fv = cos(92π/200) = 0.125 (−18 dB) against the reverb at
  fr = 0.992. Turning the knob down brings the direct voice taps and the echoes forward (knob 0: voices only, no
  reverb; 50: equal power). Like every input/output coefficient, the mix follows a knob change with the ≈ 20 ms
  smoothing of 6.4.

**Glide and Clear** (4.13) are silent in this program: they matter only when a delay length changes (an
edit, or a tempo change of the tempo-locked V2/V4), and BowedChords modulates no delay. Its post-delay
pair is GldResp 50 / GldRange 0, its voice pair GldResp 50 / GldRange 55; Clear is 0.

**Every parameter of the program.** Columns: the stored value; the value in effect; what the parameter changes when
it alone is moved (to 0, half, maximum and the stored value ± 1, patches removed), as symbols of sections 3 and 5;
the law; and whether its present value shapes the sound.

| page | parameter | stored (range) | in effect | writes (symbols of sections 3, 5) | law | audible at the stored setting |
|---|---|---|---|---|---|---|
| Controls | InLvl R | 160 (0–160) | AR Env → [0:120 127:160] | aRR | input level law 4.5; driven by the AR envelope (above) | yes (AR: 120 … 160 = −40 … 0 dB) |
| Controls | InLvl L | 160 (0–160) | AR Env → [0:120 127:160] | aLL | input level law 4.5; driven by the AR envelope (above) | yes (AR: 120 … 160 = −40 … 0 dB) |
| Controls | InPan R | 100 (0–100) | = | aRR, aLR | 4.5 | yes (100 = own channel only) |
| Controls | InPan L | 0 (0–100) | = | aRL, aLL | 4.5 | yes (0 = own channel only) |
| Controls | FX Mix | 43 (0–100) | 92 (ADJUST knob, curve [0:0 100:100]) | fv, fr | equal-power mix 4.5 | replaced by the ADJUST patch (stored 43 has no effect) |
| Controls | FX Adjust | 71 (0–80) | = | go | output gain go, 4.5 | yes |
| Controls | FX Width | 405 (0–720) | = | w1, w2 | width law 4.5 | no: 405 = unchanged |
| Controls | Voice Dif | 0 (0–100) | = | g (voice diffusion), vd, vp | 4.5 | no: 0 = diffusion off |
| Controls | High Cut | 121 (0–121) | = | a (output low-pass) | 4.6 | no: 121 = no filter (a = 1) |
| Rvb Time | Rt HC | 25 (0–48) | = | k20, k22, k24, k25, k109, k111, k113, k114 | 4.4 | yes |
| Rvb Time | Low Rt | 6 (0–9) | = | k32, k122 | 4.3 | yes |
| Rvb Time | Crossover | 8 (0–60) | = | k26, k28, k116, k118 | 4.4 | yes |
| Rvb Time | Mid Rt | 41 (0–63) | = | k4, k11, k30, k32, k38, k45, k120, k122 | 4.2 | yes |
| Rvb Time | Pre Delay | 0 (0–465) | = | A13, A14, A24, A25, A47, A48, A68, A86, A113, A114 | 4.5 | yes (0 → offset 1) |
| Rvb Time | EkoDly L | 0 (0–600) | = | A101 | 4.5 | no while EkoFbk L = 15 |
| Rvb Time | EkoDly R | 0 (0–600) | = | A91 | 4.5 | no while EkoFbk R = 15 |
| Rvb Time | EkoFbk L | 15 (0–30) | = | k101 | 4.5 | no: 15 = no feedback |
| Rvb Time | EkoFbk R | 15 (0–30) | = | k91 | 4.5 | no: 15 = no feedback |
| Rvb Time | RefDly L | 7 (0–600) | = | A87 | 4.5 | yes |
| Rvb Time | RefDly R | 8 (0–600) | = | A69 | 4.5 | yes |
| Rvb Time | RefLvl L | 10 (0–15) | = | k87 | 4.5 | yes |
| Rvb Time | RefLvl R | 10 (0–15) | = | k69 | 4.5 | yes |
| Rvb Time | PstDly L | 0 (0–1365) | = | PL | 4.12 | no while PstMix = 0 |
| Rvb Time | PstDly R | 0 (0–1365) | = | PR | 4.12 | no while PstMix = 0 |
| Rvb Time | PstMix | 0 (0–100) | = | pd, pr | 4.5 | no: 0 = post-delays bypassed |
| Rvb Time | GldResp | 50 (0–100) | = | glide rate PR, glide rate PL | glide rule 4.13 | no (only when a post-delay length changes) |
| Rvb Time | GldRange | 0 (0–1365) | = | glide range PR, glide range PL | glide rule 4.13 | no (only when a post-delay length changes) |
| RvbDesign | Size | 90 (0–144) | = | 69 tank offsets | 4.1 | yes |
| RvbDesign | AP Sweep | 0 (0–50) | = | AP Sweep sequencer (6.1) | 6.1 | no: 0 = stopped, pairs at rest |
| RvbDesign | Diffusion | 86 (0–100) | = | k89, k93, k95, k97, k99, k103, k105, k107 | 4.5 | yes |
| RvbDesign | Spread | 255 (0–255) | = | A13, A14, A24, A25, A48, A113, A114 | 4.8 | yes (maximum; 254 gives the same offsets) |
| RvbDesign | Shape | 90 (0–255) | = | k13, k14, k24, k25, k47, k48, k68, k86, k113, k114 | 4.7 | yes |
| RvbDesign | Link | 1 (0–1) | = | k4, k11, k30, k32, k38, k45, k120, k122; A13, A14, A24, A25, A48, A113, A114 | 4.2 | yes (1: K fixed, so the decay follows Size: Size 90 gives ≈ 0.66 × RT[Mid Rt]) |
| RvbDesign | Rvb Out | 15 (0–15) | = | k70, k88 | 4.5 | yes |
| RvbDesign | Rvb In | 80 (0–80) | = | gs | 4.5 | yes (80 = 0 dB) |
| RvbDesign | Rvb Width | 405 (0–720) | = | c1, c2 | 4.5 | no: 405 = unchanged |
| Levels | Master | 40 (0–80) | = | V4→R, V3→R, V2→L, V1→L | 4.12 | yes |
| Levels | Voice1 | 160 (0–160) | = | V1→R, V1→L; level smoothing | 4.12 | yes |
| Levels | Voice2 | 154 (0–160) | = | V2→R, V2→L; level smoothing | 4.12 | yes |
| Levels | Voice3 | 160 (0–160) | = | V3→R, V3→L; level smoothing | 4.12 | yes |
| Levels | Voice4 | 154 (0–160) | = | V4→R, V4→L; level smoothing | 4.12 | yes |
| DelayTime | Master | 100 (0–200) | = | V2, V3, V4 | 4.12 | no at 100 % |
| DelayTime | Voice1 | 0 (0–1365) | = | V1 | 4.12 | yes (0 ms: no delay) |
| DelayTime | Voice2 | 129 [tempo] (0–1365) | = | V2 | tempo-locked 4.11 / 4.12 | yes (¼ beat) |
| DelayTime | Voice3 | 5 (0–1365) | = | V3 | 4.12 | yes |
| DelayTime | Voice4 | 129 [tempo] (0–1365) | = | V4 | tempo-locked 4.11 / 4.12 | yes (¼ beat) |
| DelayTime | Clear | 0 (0–1) | = | fb4, V4→R, V4→L, V3→R, V3→L, fb2, V2→R, V2→L; level smoothing | glide rule 4.13 | no: 0 |
| DelayTime | GldResp | 50 (0–100) | = | glide rate V1, glide rate V2, glide rate V3, glide rate V4 | glide rule 4.13 | no (only when a voice delay changes) |
| DelayTime | GldRange | 55 (0–1365) | = | glide range V1, glide range V2, glide range V3, glide range V4 | glide rule 4.13 | no (only when a voice delay changes, e.g. a tempo change) |
| Feedback | Master | 48 (0–100) | = | fb4, fb2; level smoothing | 4.12 | yes |
| Feedback | Voice1 | 100 (0–200) | = | fb1; level smoothing | 4.12 | no: 100 = no feedback |
| Feedback | Voice2 | 5 (0–200) | = | fb2; level smoothing | 4.12 | yes |
| Feedback | Voice3 | 100 (0–200) | = | fb3; level smoothing | 4.12 | no: 100 = no feedback |
| Feedback | Voice4 | 5 (0–200) | = | fb4; level smoothing | 4.12 | yes |
| Panning | Master | 50 (0–100) | = | V4→R, V4→L, V3→R, V3→L, V2→R, V2→L, V1→R, V1→L | 4.12 | yes |
| Panning | Voice1 | 0 (0–100) | = | V1→R, V1→L | 4.12 | yes |
| Panning | Voice2 | 0 (0–100) | = | V2→R, V2→L | 4.12 | yes |
| Panning | Voice3 | 100 (0–100) | = | V3→R, V3→L | 4.12 | yes |
| Panning | Voice4 | 100 (0–100) | = | V4→R, V4→L | 4.12 | yes |

The parameters that are silent at the stored setting (High Cut, Voice Dif, FX Width, Rvb Width, EkoFbk / EkoDly,
PstMix / PstDly, AP Sweep, both Glide pairs, Clear, Delay Master, Feedback V1 / V3, and Spread above 254) all act as soon
as they are edited. FX Mix acts only if its ADJUST patch is removed.

### 7.6 Plate programs: Vox Plate, Good olPlate, Slap Plate, Brass Plate, Drum Plate, Rich Plate, Plate4 Horns

Values: program data (17.2), records with `part1_section` = "7.6".

### 7.7 Plate program Alley Slap

Values: program data (17.2), records with `part1_section` = "7.7".

### 7.8 Chorus+Rvb program Deep Space (Plate tank 5.2, input/output section 3.1)

Values: program data (17.2), records with `part1_section` = "7.8".

### 7.9 Concert Hall programs: Small+Stage, Concert Wave, Concert Hall, ConcertHall2, Piano Hall, Medium Hall, Vocal Hall

Values: program data (17.2), records with `part1_section` = "7.9".

### 7.10 Concert Hall program Whammy Hall

Values: program data (17.2), records with `part1_section` = "7.10".

### 7.11 Inverse programs: Gate Chamber, Drum Gate, Slope Down

Values: program data (17.2), records with `part1_section` = "7.11".

## 8. Accuracy

The tank programs (5), the input/output sections (3, 3.1) and every law marked **exact** reproduce the original for
every value of the parameters involved. The remaining limits:

* **Low-band decay terms** (4.3): within 1.2/256, exact for 96 % of the values.
* **AP Sweep step time** (6.1): 0.130·P_k samples per step, within ±0.7 %.
* **Modulated tap step time** (6.2): 0.2575·(D + 4.1) samples, within 1 % over Chorus 1 … 10.
* **AP Sweep-pair signs after a program change** are not defined (6.1); a program load uses the stock signs.
* **Parameter changes** (6.4): the ramp durations are accurate to ±15 %. A live edit that lands exactly on "off"
  follows the parameter path.
* **BowedChords input gain** (7.5): the latency and smoothing are approximate; combinations with the same total
  delay of ≈ 40 ms are equivalent.
* The Concert Hall decays faster than its displayed time (4.2); this is part of the algorithm, not an error.
* The delay memory (standard or extended) affects only the effect algorithms of Part II (11.1).

## 9. Presenting the parameters in a simulation

A program has 40–80 parameters, but a listener controls a reverb through a handful of perceptual qualities: how long
it lasts, how big and how far away it sounds, how bright it is, how dense or grainy it is, whether it moves, how wide
it is, and how loud it is against the dry signal. This section groups the parameters by those qualities and ranks them
in four layers:

* **Layer 0, program controls**: the same for every program.
* **Layer 1, primary**: the few controls that define the sound of a given algorithm.
* **Layer 2, shaping**: character and fine balance.
* **Layer 3, fine adjustment**: the delay-voice network, levels, and the modulation patches.

A fifth group should not be shown at all.

The ranking follows how the original's own factory programs use the parameters. Each program stores a 10-slot
quick-access selection of its most relevant parameters and one ADJUST-knob macro; how often a parameter appears there
(table 9.6), and how much its value varies between the factory programs, decides its layer. Where that data is thin,
the parameter's audible effect decides.

### 9.1 Layer 0: program controls (all algorithms)

| control | function | notes |
|---|---|---|
| Program | selects algorithm + all values | the factory programs of section 7 are the natural starting points |
| **Mix** | wet/dry, equal power (4.5) | in the quick-access selection of every factory program |
| **ADJUST** | one macro per program: drives one or more parameters through the program's patch curves (4.10) | show it with the program's own label (e.g. "Decay", "Liveness", "FX Mix") and the patched parameters' resulting values; most programs map it to Mid Rt, many Chamber programs to Shape / Rt HC / Crossover, Chorus+Rvb programs to FX Mix or the chorus depth/rate |
| Tempo (BPM) | beat length of tempo-locked times (4.11) | show only for programs with tempo-locked delays (Vox Plate, Vocal Hall, Alley Slap, Deep Space, BowedChords) |
| Input / output level | gain before the program, trim after it | not part of the programs; the program's own input level (InLvl) and output gain (FX Adjust) belong to layer 3 / 2 |

### 9.2 Layers 1–3 per algorithm

`1`, `2`, `3` = layer; `–` = the parameter does not exist or has no effect in that algorithm (9.4).

| quality | parameter | Plate | Chamber | Inverse | Concert Hall | Chorus+Rvb |
|---|---|---|---|---|---|---|
| length | **Mid Rt** (decay time) | 1 | 1 | – | 1 | 1 |
| length | Duration (Inverse: length of the gated/reverse envelope) | – | – | 1 | – | – |
| length | Low Slope / Mid Slope (Inverse: envelope shape per band) | – | – | 1 | – | – |
| length | Low Rt (bass decay multiplier) | 2 | 1 | – | 1 | 2 |
| length | Crossover (bass/mid split frequency) | 2 | 2 | 2 | 2 | 3 |
| length | Link (decay follows Size: off / on) | 3 | 3 | – | 3 | 3 |
| size / distance | **Size** | 1 | 1 | – | 1 | 1 |
| size / distance | **Pre Delay** | 1 | 1 | 1 | 1 | 2 |
| size / distance | RefDly + RefLvl (early reflection, one pair per side) | 2 | 2 | 2 | 2 | 3 |
| size / distance | PstMix + PstDly L/R (post-delay of the reverb output) | 3 | 2 | 2 | 2 | – |
| brightness | **High Cut** (output low-pass) | 1 | 1 | 1 | 1 | 1 |
| brightness | **Rt HC** (high-frequency damping of the decay) | 1 | 1 | 2 | 1 | 2 |
| density / texture | Diffusion | 2 | 1 | 1 | 1 | 2 |
| density / texture | Shape (Chamber, Inverse: build-up envelope) | – | 1 | 1 | – | – |
| density / texture | Spread (Chamber: spacing of the injection taps) | – | 2 | – | – | – |
| density / texture | Attack (Plate: injection balance, slow build-up) | 2 | – | – | – | 3 |
| motion | AP Sweep (allpass sweep) | 2 | 2 | – | 3 | 3 |
| motion | Depth (Concert Hall: tank mixers, also balance) | – | – | – | 1 | – |
| motion | Def (Concert Hall: AP Sweep allpass gain, "definition") | – | – | – | 2 | – |
| motion | Chorus (Concert Hall: modulated tap speed) | – | – | – | 2 | – |
| motion | MstDepth / MstRate (chorus voices) | – | – | – | – | 1 |
| width | Rvb Width (rotation of the reverb outputs) | 2 | 2 | 2 | 2 | – |
| width | FX Width (rotation of the final output) | 3 | 3 | 3 | 3 | 2 |
| balance | FX Mix (delay voices / chorus vs reverb) | 2 | 2 | 3 | 2 | 1 |
| balance | FX Adjust (output gain) | 3 | 2 | 3 | 3 | 2 |
| voices | Delay Master, Feedback Master | 2 | 3 | 3 | 3 | 1 |
| voices | Level Master | 3 | 3 | 3 | 3 | 2 |
| voices | Voice Dif (diffusion of the voice inputs) | 2 | 3 | 3 | 3 | – |
| voices | per voice: delay, level, pan, feedback; Pan Master | 3 | 3 | 3 | 3 | 3 |
| voices | per voice: chorus depth, rate | – | – | – | – | 3 |
| echo | EkoDly + EkoFbk (tank input echo, per side) | 3 | 3 | – | – | 3 |
| levels | Rvb In, Rvb Out (tank send / return level) | 3 | 3 | 3 | 3 | 3 (Rvb Out only) |
| levels | InLvl, InPan (program input level and pan) | 3 | 3 | 3 | 3 | 3 |
| delay changes | GldResp, GldRange (two pairs: post-delays, voices) | 3 | 3 | 3 | 3 | 3 (voices) |
| modulation | modulator settings, patches, trigger levels (6.3, 4.10) | 3 | 3 | 3 | 3 | 3 |

Practical layouts:
* **Layer 1** fits one screen of five to eight controls per algorithm:
  * Plate: Mid Rt, Size, Pre Delay, High Cut, Rt HC.
  * Chamber: Mid Rt, Size, Pre Delay, High Cut, Rt HC, Low Rt, Diffusion, Shape.
  * Inverse: Duration, Low Slope, Mid Slope, Shape, Diffusion, Pre Delay, High Cut.
  * Concert Hall: Mid Rt, Size, Pre Delay, High Cut, Rt HC, Low Rt, Diffusion, Depth.
  * Chorus+Rvb: Mid Rt, Size, High Cut, MstDepth, MstRate, Delay Master, Feedback Master, FX Mix.
* **Layer 2** is a second page, "character".
* **Layer 3** is an editor page per section: delay voices, levels, modulation.
* The ADJUST knob (layer 0) stays visible on every page. It is the program designer's choice of the one control that
  matters most for that program.

### 9.3 Couplings to make visible

* **Size and decay**: with Link on, the decay time grows with Size (RT ≈ RT[Mid Rt]·(16 + Size)/160 for Plate and
  Chamber, 4.2). Display the resulting time in seconds next to Mid Rt, or offer "Decay follows Size" as a
  switch.
* **Link** also decides whether the input taps (Plate early build-up, Chamber Spread taps) follow Size (4.1, 4.8).
* **Rt HC and Shape (Chamber):** four Shape coefficients are scaled by the damping coefficient (4.7), so changing
  Rt HC also changes the early build-up.
* **Concert Hall** decays at ≈ 0.7 × its displayed time (4.2). Show the original's display value for compatibility
  and, optionally, the actual decay time.
* **Controls that act only together**:
  * EkoDly needs EkoFbk ≠ 15 (15 = no feedback).
  * PstDly needs PstMix > 0.
  * The voice delays need a non-zero voice level and FX Mix < 100.
  * Voice Dif 0 switches the voice diffusion off.
  * RefDly needs RefLvl > 0.

  Grey out the dependent control while its partner disables it.
* **ADJUST overrides**: a parameter driven by an ADJUST patch ignores its stored value (e.g. FX Mix in BowedChords,
  Mid Rt in most programs). Show it as driven by the knob, not as an independent control.
* **Modulated parameters**: parameters driven by an LFO, envelope or input level (6.3) change on their own. Show the
  live value and mark the control as modulated.
* **AP Sweep = 0** stops the sweep. The pair signs left after a program change are not defined (6.1), so a "reset" of the
  AP Sweep state has no audible meaning beyond the fine structure.

### 9.4 Not to present

* **No effect on the sound**: Spread in Plate and Concert Hall, Attack in Concert Hall, AP Sweep and Spread in Inverse,
  Spread and Rvb In in Chorus+Rvb (they sit on no edit page).
* **Clear** is an action rather than a setting: it freezes the voice levels and feedback (4.13).
  Offer it as a momentary command at most.
* **Fixed in practice** and safe to hide behind defaults:
  * InLvl 160 (0 dB) / InPan own channel in nearly all programs (a few Plate and Chorus+Rvb programs use −6 dB or
    centred pans; some programs drive them from a modulator).
  * Rvb In 80 (0 dB).
  * Pan Master 50.
  * Delay Master 100 % in Chamber and Inverse.

### 9.5 Display units and behaviour when changed

| parameters | display | when changed |
|---|---|---|
| Mid Rt | seconds: RT[v] (4.2), × size factor with Link on | tank coefficient ramp (6.4) |
| Low Rt | multiplier ×0.2 … ×4.0 (4.3) | tank coefficient ramp (6.4) |
| Rt HC, Crossover | Hz from the table in 4.4 (Rt HC 48 = "Full") | tank coefficient ramp (6.4) |
| High Cut | −3 dB frequency (4.6), 121 = off | input/output coefficient: one-pole glide ≈ 20 ms (6.4) |
| Size, Duration | value (no physical unit is defined) | **tank offsets change without a ramp** (6.4): heard as a jump of the taps; offer as non-automatable, or cross-fade two instances |
| Pre Delay, RefDly, EkoDly | milliseconds, 2·v (4.5) | tank offsets change (as Size) |
| Spread | value 0 … 255 | tank offsets change |
| Diffusion, AP Sweep, Def | percent of maximum, ⌊100·v/max⌋ %, 0 = "Off", max = "Full" | Diffusion: tank coefficient ramp; AP Sweep: next sweep step; Def: through the sweep steps (6.1) |
| Shape, Attack | value 0 … 255 / 0 … 100 | tank coefficient ramps (6.4) |
| Depth, Chorus | value | tank coefficient ramps / tap speed at once |
| Rvb Width, FX Width | rotation θ = 45° + v° (405 = unchanged stereo, 4.5) | input/output coefficients, glide ≈ 20 ms |
| Mix, FX Mix, PstMix | percent (equal-power law) | glide ≈ 20 ms (6.4) |
| Rvb In, FX Adjust, InLvl, voice levels | dB (4.5) | glide ≈ 20 ms (6.4) |
| voice / post-delay times | milliseconds or beats (tempo-locked) | glide per GldResp up to GldRange, larger jumps cross-fade over 171 samples (4.13) |
| voice feedback, pans | percent | glide ≈ 20 ms (6.4) |
| Rvb Out, RefLvl, EkoFbk | value 0 … 15 / 0 … 30 | tank coefficient ramps (6.4) |

### 9.6 Quick-access selections of the factory programs

How many factory programs of each algorithm put the parameter into their 10-slot quick-access selection (Plate 25
programs, Chamber 38, Inverse 7, Concert Hall 16, Chorus+Rvb 35; a slot holding a parameter pair counts for both, e.g.
RefDly/RefLvl), and the number of ADJUST-knob patches that target it:

| parameter | Plate | Chamber | Inverse | Concert Hall | Chorus+Rvb | ADJUST targets |
|---|---|---|---|---|---|---|
| Mix | 25 | 38 | 7 | 16 | 35 | — |
| Mid Rt | 11 | 26 | – | 12 | 28 | 31 |
| Size | 13 | 37 | – | 14 | 22 | — |
| Duration / Low Slope / Mid Slope | – | – | 2 / 5 / 5 | – | – | Duration 3 |
| High Cut | 22 | 30 | 6 | 13 | 21 | 4 |
| Rt HC | 13 | 23 | 2 | 13 | 13 | 7 |
| Low Rt | 12 | 30 | – | 15 | 7 | 1 |
| Pre Delay | 12 | 22 | 4 | 12 | 9 | 2 |
| Diffusion | 7 | 21 | 4 | 12 | 9 | 4 |
| Shape | – | 16 | 4 | – | – | 10 |
| Depth / Def / Chorus (Hall) | – | – | – | 11 / 1 / 3 | – | Depth 3, Def 3 |
| MstDepth / MstRate | – | – | – | – | 26 / 0 | 7 / 7 |
| Feedback Master / Delay Master | 9 / 12 | 3 / 4 | 0 / 0 | 2 / 2 | 29 / 27 | 1 / 5 |
| FX Mix | 11 | 8 | 1 | 3 | 19 | 18 |
| FX Width | 13 | 9 | 4 | 1 | 20 | 2 |
| Rvb Width | 6 | 14 | 3 | 8 | – | 3 |
| Crossover | 5 | 9 | 1 | 3 | 2 | 5 |
| RefDly / RefLvl | 6 | 8 | 2 | 0 | 2 | 6 (RefLvl) |
| FX Adjust | 2 | 11 | 1 | 0 | 10 | — |
| Voice Dif | 17 | 9 | 0 | 1 | – | — |
| Attack | 3 | – | – | – | 2 | 1 |
| Spread (Chamber) | – | 3 | – | – | – | 1 |
| AP Sweep | 0 | 1 | – | 0 | 1 | — |
| PstMix | 2 | 6 | 2 | 1 | – | 3 |
| Level Master | 4 | 2 | 1 | 0 | 11 | 3 |
| EkoDly / EkoFbk | 4 | 0 | – | – | 2 | 2 (EkoFbk) |
| per-voice parameters | 0–4 | 0 | 0 | 0–1 | 0–4 | 32 |
| Link, Rvb In, Rvb Out, InLvl/InPan, Glide, Clear | 0 | 0–10 (Rvb In) | 0 | 0 | 0–3 (InLvl) | Glide 8, Rvb Out 2, InPan 2, Rvb In 1 |

Some placements in 9.2 deviate from the counts. AP Sweep rarely appears in the quick-access selections, but its value
varies widely between programs (up to 15 different settings) and it is the only motion control of Plate and Chamber,
so it is layer 2 there. FX Width and Voice Dif appear often in Plate selections, but almost every program leaves them
at their defaults, so FX Width is layer 3 (layer 2 only in Chorus+Rvb, where its value varies) and Voice Dif is
layer 2 only in Plate. Rvb In appears in 10 Chamber selections, but every program keeps 0 dB, so it is layer 3.

## 10. Analog input path

How the unit turns an analog input into the 48 kHz samples its program receives, including all overload behaviour:
two analog hard clips, the INPUT knob between them, and a hard clip on the ADC's output samples, which aliases.
Scope: the dry path at 48 kHz, from the input jack to the
ADC's output samples, before *Analog In Lvl* and the program. The output stage and D/A are not covered.

### 10.1 Signal chain

```
input jack
  → input stage (rear bal/unbal switch: Out 0 dB / In +20 dB)   ── analog hard clip, asymmetric, no aliasing
  → INPUT knob (front panel)                                    ── gain only
  → ADC input buffer                                            ── analog hard clip at 1.16 × ADC full scale
  → ADC anti-alias / decimation filter (48 kHz)
  → ADC full scale                                              ── hard clip on the 48 kHz samples: aliases
  → Analog In Lvl (digital gain, 0–100 %)
  → program (dry/wet mix, effects)
```

Two consequences shape all overload behaviour:

* With the knob above about **−23 dB from its maximum**, the ADC clips first. The clip is digital (on the samples),
  so every harmonic above 24 kHz folds back into the audio band. The buffer limit 1.3 dB higher stops the harmonics
  growing under heavy drive.
* With the knob further down, the **input stage** clips first. That clip is analog: its harmonics above 24 kHz are
  removed by the ADC filter, so nothing aliases.

### 10.2 Analog drive and digital drive

Only two gains matter: the one in front of the input-stage clip and the one between it and the ADC. A simulation
can expose exactly those two, each as a level relative to the clip it drives:

| parameter | definition | on the unit |
|---|---|---|
| **Analog drive A** (dB) | input level relative to the input-stage clip (0 dB = sine peak just reaching it) | source level and rear switch (In adds 20 dB) |
| **Digital drive K** (dB) | input-stage clip level relative to ADC full scale | the INPUT knob: at most +23.2 dB (switch Out) / +24.25 dB (In) at full knob |

The small-signal gain into the ADC is **A + K**. How the two interact:

* **Raising A** adds input-stage (analog) distortion once A > 0. Until then it drives the ADC like any other gain,
  and the ADC clips (digitally, with aliasing) once A + K > 0. Above A = 0 the input stage caps the ADC drive: its
  clipped output, band-limited by the ADC filter, reaches up to 4/π of the clip level (+2.1 dB; a clipped tone above
  8 kHz keeps only its fundamental). So the ADC never clips when K is below about −2.5 dB, and all overload is then
  analog; at 7 kHz this holds up to K = 0 dB.
* **Raising K** is exactly turning the INPUT knob up. It scales whatever the input stage delivers (clean or clipped)
  into the ADC, so it controls how much digital clipping there is.
* **Make-up gain:** to keep the loudness of clean signals constant while turning the two drives, scale the output by
  −(A + K) dB.
* **The rear switch** becomes a character choice. Its 20 dB of gain is absorbed into A, and it selects the input-stage
  variant: asymmetry 1.24 % (Out) or 3.0 % (In), and the maximum K. The In variant's clip is also 1.0 dB higher
  relative to the switched signal, which the A definition absorbs too.

In these terms the model is (per channel):

```
s = x * 10^(A/20)                        # x: program signal; A = 0 dB puts its peak at the input-stage clip
s = clip(s, -(1-a), +(1+a))              # input stage, oversampled; a = 0.0124 (Out) or 0.030 (In)
s = s * 10^(K/20)                        # INPUT knob
s = clip(s, -1.16, +1.16)                # ADC input buffer, oversampled
y = decimate(s)                          # ADC filter -> 48 kHz
y = clip(y, -1, +1)                      # ADC full scale on the samples (aliases)
y = y * 10^(-(A+K)/20)                   # optional make-up gain
```

It is identical, sample for sample, to the hardware-units version of 10.3.

**Hardware constants.** All levels are sine **peak** values on one scale: **1.0 = ADC full scale, with the knob at
maximum and the switch at Out.** The absolute input level in dBu is not specified.

| parameter | value | notes |
|---|---|---|
| Switch gain, In | ×10.0 (+20.0 dB) | |
| Input-stage clip, switch Out | ±14.45 (+23.2 dB), asymmetry a = +1.24 % | positive limit 14.45·(1+a), negative 14.45·(1−a) |
| Input-stage clip, switch In | ±16.31 (+24.25 dB), asymmetry a = +3.0 % | level after the +20 dB switch gain |
| INPUT knob gain | 1.0 at maximum (0 < knob_gain ≤ 1) | flat with frequency; the dial law is not part of the model |
| ADC input buffer | hard clip ±1.16 (+1.3 dB) | |
| ADC filter | flat to at least 21 kHz; ≥ 90 dB stopband above about 35 kHz | modelled as flat to 21.8 kHz, ≥ 100 dB from 24 kHz |
| ADC full scale | hard clip ±1.0 on the 48 kHz samples | |
| Clip onset vs frequency | flat ±0.1 dB, 200 Hz – 16 kHz | |

The asymmetry signs are in end-to-end polarity; for an inverting dry path, flip both
signs. Both channels use the same asymmetry.

### 10.3 Algorithm in hardware units

The same model, written with the unit's own controls. Per channel, independently (both channels use the same
parameters):

```
# x: input samples at 48 kHz, in the units of 10.2
xo = upsample(x, OS)                    # analog-domain stand-in; OS = 16
if switch_in: xo = 10.0 * xo
L, a = (16.31, 0.030) if switch_in else (14.45, 0.0124)
xo = clip(xo, -L*(1-a), +L*(1+a))       # input stage
xo = knob_gain * xo                     # INPUT knob, 0 < knob_gain <= 1
xo = clip(xo, -1.16, +1.16)             # ADC input buffer
y  = decimate(xo, OS)                   # ADC filter -> 48 kHz
y  = clip(y, -1.0, +1.0)                # ADC full scale on the samples
y  = analog_in_lvl * y                  # then the program
```

* **Input stage.** Hard clip, asymmetric, with no knee. In the In position the stage clips 1.0 dB higher (after
  the extra 20 dB) and about 2.5 times as asymmetric. This clip is analog, so it must run in the oversampled domain:
  anything it creates above 24 kHz must be removed by the ADC filter, not folded.
* **INPUT knob.** A plain gain: the digital drive K is `20·log10(knob_gain · 14.45)` (Out) or
  `20·log10(knob_gain · 16.31)` (In). Below `knob_gain ≈ 1/14.45` (−23.2 dB) for Out, or −24.3 dB for In, K is
  negative and the input stage clips before the ADC can. The dial law (gain against rotation) is not part of the
  model: a simulation sets knob_gain directly.
* **ADC input buffer.** Hard clip at ±1.16 of ADC full scale, before the anti-alias filter. It only matters above
  +1.3 dB of drive, where it limits what reaches the sample clip. That is why the folded harmonics level off under
  heavy drive instead of growing like a plain 48 kHz clip.
* **ADC filter.** Decimation from the oversampled domain to 48 kHz. It passes 21 kHz at full level and suppresses
  components above about 35 kHz by at least 90 dB.
  The 24–26.4 kHz region, where a real half-band ADC filter may let some aliasing through, is not specified.
* **ADC full scale.** A hard clip of the 48 kHz samples at ±1.0, after the filter. **This is the only aliasing stage,
  and it must run at 48 kHz, not oversampled.** Its aliasing is part of the sound. Example: a 7001 Hz tone's
  5th, 7th and 9th harmonics land at 13.0, 1.0 and 15.0 kHz.
* **After the ADC.** *Analog In Lvl* (0–100 %) is a digital gain after the clip, so it cannot prevent or cause input
  clipping. 100 % is unity; the law for other values is not specified. The program then saturates its own
  stores at ±1 (section 1). The overload indicator lights about 1 dB below full scale.

### 10.4 Implementation notes

* **Oversampling.** 16× keeps the aliasing of the analog clips themselves negligible. A cascade of four halfband
  stages (flat to 21.8 kHz, ≥ 93 dB stopband from 26.2 kHz at 48 kHz) matches a linear-phase 768 kHz filter (flat to
  21.8 kHz, ≥ 100 dB from 24.0 kHz) within 0.1 dB in harmonics and folded products. The test is that no folded
  harmonics appear when only the input stage clips. The final ±1.0 clip stays at 48 kHz in every case.
* **Filters.** The upsampling filter stands in for a band-limited analog source; for an input that is already a
  48 kHz signal it only has to be flat to about 21 kHz. The decimation filter is the ADC's, and its stopband sets the
  alias floor of the analog clips.
* **Latency.** The two filters add their delay to the path; a real-time version reports it.
* **What is not modelled:** the 2nd harmonic of the input stage at heavy drive, the knob's dial law, and the
  input level in dBu.

# Part II — Effect algorithms and the factory programs

Part II completes Part I. It describes the effect algorithms of the reverb unit — **M-Band+Rvb**, **Glide>Hall**, **Res1>Plate**,
**Res2>Plate** and **Infinite** — and all 166 factory programs that Part I does not cover (programs 1 – 150
except BowedChords (program 139), and programs 184 – 200), at the level needed to re-implement them sample by sample.
Section 17.1 describes all 200 factory programs, including the 34 whose data is in Part I.
Programs of the Chamber, Plate, Concert Hall, Inverse and Chorus+Rvb algorithms use the descriptions of Part I unchanged, with the correction of 13.1 already applied in 3.1. As in Part I, only the signal
processing is described. Rules marked **exact** reproduce the original for every value of the parameters involved.
The values of the individual programs are not printed here: they are in the data file `effect-programs.json`,
which the script `effect_programs.py` prints per program or side by side (section 17.2).

| algorithm | programs | of these in section 7 | tank (section) | effect section |
|---|---|---|---|---|
| M-Band+Rvb | 50 | – | Chamber (5.1) | six band-filtered delay taps with feedback (13.2) |
| Chorus+Rvb | 35 | 1 | Plate (5.2) | six modulated chorus voices (3.1, corrected in 13.1) |
| Glide>Hall | 18 | – | Concert Hall (5.3) | four gliding interpolated taps, six echo voices (13.3) |
| Res1>Plate, Res2>Plate | 6 + 4 | – | Plate (5.2) | six tuned resonators fed by six delay taps (13.4) |
| Infinite | 1 | – | own (15.1) | reverb input/output section of 3 (13.5) |
| Chamber, Plate, Concert Hall, Inverse | 38, 25, 16, 7 | 14, 8, 8, 3 | 5.1 – 5.4 | 3 |

Sample rate **fs = 48 kHz**; delays in samples at 48 kHz (44.1 kHz: 14.7).

## 11. Conventions

Everything in 1 applies. In addition:

* **16-bit coefficient words.** The coefficients of the input/output sections are computed as 16-bit integers n
  (1.0 = 32768) and written as the 24-bit word 256·n + (n mod 256), value word/2^23 ≈ n/32768 (the low byte is
  repeated; n = 32767 gives 0x7FFFFF). The program data (17.2) gives the resulting values. Where this document states
  a law in 16-bit arithmetic, `rmul(a, b)` = ⌊a·b/32768 + ½⌋, `P16(p)` = ⌊32767·p/100⌋, and `T[v]` is the dB table of
  4.5 (T[v] = 32768·10^((v − 80)/20), T[80] = 32767), extended to the signed input-level law
  `L16(v)` = T[v − 80] for v > 80, 0 for v = 80, −T[80 − v] for v < 80.
* **Saturation** as R1: every stored value saturates at ±1, including the stored states of the band filters and
  resonators and every write into a delay line.

### 11.1 Delay memory size and the resolution r of the effect delay lines

The long delay lines of the effect algorithms (M-Band band taps, Glide>Hall echo voices, Res input taps) live in a
separate delay memory whose size depends on the unit's memory fitting. With B banks of 65 536 words, a line holds
65 536·r samples (Glide: 63 488·r, Res: 57 344·r) and **a tap can only be set in steps of r = B/2 samples**:

| memory | r | longest M-Band delay | Glide / Res voice range |
|---|---|---|---|
| standard (4 banks) | 2 samples | 2.730 s | 2.645 s / 2.389 s |
| extended (64 banks) | 32 samples | 43.69 s | 42.325 s / 38.229 s |

A voice time of t ms (whole ms after the delay master, 14.1) gives the tap **n = ⌊48·t/r⌋** and the delay
**D = r·n samples** (exact; e.g. 445, 449, 455 ms → 21 344, 21 536, 21 824 samples, while
even-ms values such as 224, 448, 450, 600, 900, 1300 ms land on 48·t). With the standard memory **D = 48·t exactly**, for every t up to the range. The program data (17.2) gives D for
the extended memory and, with the suffix `_std`, for the standard memory. A tap reads the value written exactly D
samples earlier ("write before read": D = 0 is the current input). The short lines (Glide taps, allpasses, resonator
loops, all reverb-algorithm lines and the tank) do not depend on the memory size, with one exception: with the
standard memory the post-delays of the reverb algorithms (3) are limited to 682 ms (32 736 samples) — a plain value
above it is not sent (the post-delay keeps its previous length), and a patch result is clamped there (4.10). No factory
program uses a post-delay above 682 ms. With the extended memory they reach 1 365 ms.

## 12. Building blocks

In addition to R2:

**Two-integrator state-variable section** (band filters, 13.2). Input x, states lp, bp, damping d, integrator gain f:

    hp = x − lp − d·bp ;   lp ← lp + f·bp ;   bp ← bp + f·hp

(hp uses the old lp and bp; lp is updated with the old bp). f = 2·sin(π·fc/fs) sets the corner, d = 1/Q the damping.

**Coarse delay tap.** `t = line[n − D]`, D a multiple of r (11.1).

**Gliding interpolated tap** (Glide>Hall, 13.3). Position p in units of 1/4096 sample, n = ⌊p⌋, f = p − n:

    t = f·x[n − n_p] − (f + c)·x[n − n_p + 1]

c = −1 is linear interpolation for a delay of n_p − 1 + f samples; the program also allows other c (a first-order
weighting), but no factory program uses one (14.3). Each sample the position glides `p += ρ·(P − p)` toward its target P (the tap uses the position before the
update).

**Resonator** (13.4): a delay loop with a one-pole low-pass in the loop and a linearly interpolated read, see 13.4.

**Tap changes.** When a coarse tap changes, the old and the new tap are
cross-faded linearly: `t = w·line[n − D_old] + (1 − w)·line[n − D_new]`, w falling from 1 to 0 in equal steps, over
**128 samples** in M-Band+Rvb and Glide>Hall (step 1/128) and **171 samples** in Res1/Res2>Plate (step 3/512, the
length of the reverb voices' cross-fade, 4.13). Only one tap of an algorithm fades at a time: the six taps are
checked in turn and no new fade starts while one runs, so a tap that changes meanwhile waits, and a change of all six
(e.g. the Delay Master) takes up to 6 × 128 (16 ms) or 6 × 171 samples. The same sample counts apply at 44.1 kHz.

## 13. Signal flow of the effect algorithms

### 13.1 The band and chorus line inputs (correction to 3.1)

Both M-Band+Rvb and Chorus+Rvb feed their delay lines with **2 × the reverb send of the same sample**
(`2·sendL[n]` into the left line, `2·sendR[n]` into the right line). The send is the value written to the tank input
(3, `sendL = gs·mL`), not the tank's diffused input line. 3.1 therefore changes to

    fA = 2·sendL[n] + c96·v1 + c97·v2 + c98·v3 ;   fB = 2·sendR[n] + c99·v4 + c9A·v5 + c9B·v6

(with the send of 3.1's first line). The M-Band band output and the Chorus+Rvb voices at 0 ms therefore arrive with the dry signal, not after the tank's
input diffusers.

### 13.2 M-Band+Rvb

Tank: Chamber (5.1), fed by send[n−1] as in 3. Per sample n (A = left, B = right; names in brackets = the program-data
quantities, 17.2):

```
# send: input matrix plus band feedback (h_j of the previous sample)
mL = aLL·inL + aLR·inR + fb1·h1 + fb2·h2 + fb3·h3                    [send matrix, band feedback]
mR = aRL·inL + aRR·inR + fb4·h4 + fb5·h5 + fb6·h6
sendL[n] = gs·mL ;   sendR[n] = gs·mR                                [send gain gs]
# band lines (written before the reads)
write 2·sendL[n] into line L, 2·sendR[n] into line R
t_j = line L (j = 1..3) / line R (j = 4..6) delayed by D_j samples     [band delays, 11.1]
# two state-variable sections per band (12); h_j is the band output
for j = 1..6 (H = HiCut set (d1, f1, c), L = LoCut set (e, d2, f2) of band j, 14.2):
    hp = g_j·t_j − lp − d1·bp ;  lp' = lp + f1·bp ;  bp' = bp + f1·hp ;  x = 8·(lp' + c·hp)
    hq = e·x − lq − d2·bq ;      lq' = lq + f2·bq ;  bq' = bq + f2·hq ;  h_j = hq
# mix and output
bandL = Σ pL_j·h_j ;  bandR = Σ pR_j·h_j                              [band pans]
oL = fv·bandL + fr·retL[n−1] ;  oR = fv·bandR + fr·retR[n−1]          [FX mix]
wL = w1·oL + w2·oR ;  wR = w1·oR + w2·oL                              [output width]
outL = wet·4·go·wL + dry·inL ;  outR = wet·4·go·wR + dry·inR          [output gain, Mix]
```

Every stored value (lp, bp, lq, bq, x, h_j, the line inputs) saturates at ±1. The band feedback closes a loop of
D_j + 1 samples. The first section's output mixes low-pass and high-pass (`lp' + c·hp`) and is scaled by 8 to undo
the band input level g_j (= 0.125 × the level law, 14.2); the second section's high-pass output is the band output.

### 13.3 Glide>Hall

Tank: Concert Hall (5.3) with its own decay loop time K in the Mid Rt law (4.2; 225 ms with Link on, 14.3), fed by send[n−1]. Unlike the other
algorithms, the tank is fed by the echo voices, not by the input.

```
# glide section: two 2048-sample lines (L, R), two gliding taps each (A, B)
gL = aLL·inL + aLR·inR + φL ;   gR = aRL·inL + aRR·inR + φR          [input matrix; φ = glide feedback, previous sample]
write gL into glide line L, gR into glide line R
tAL = tap(L, p_AL, cL) ; tBL = tap(L, p_BL, −1) ; tAR = tap(R, p_AR, cR) ; tBR = tap(R, p_BR, −1)   [2: gliding tap]
oL = λAL·tAL + λBL·tBL ;   oR = λAR·tAR + λBR·tBR                    [glide tap levels]
φL' = κLL·oL + κRL·oR ;    φR' = κLR·oL + κRR·oR                      [glide feedback: Fbk L, X-Fbk R, X-Fbk L, Fbk R]
dL = gg·oL ;  dR = gg·oR                                             [glide output level gg]
# echo section: two coarse lines (r-sample resolution, 11.1), three voices each
write dL + ηL into echo line L, dR + ηR into echo line R             [η = echo feedback, previous sample]
t_j = echo line L (j = 1..3) / R (j = 4..6) delayed by D_j
S1 = Σ_{j=1..3} sf_j·t_j ;  S2 = Σ_{j=1..3} sx_j·t_j ;  S3 = Σ_{j=4..6} sf_j·t_j ;  S4 = Σ_{j=4..6} sx_j·t_j
                                                                    [echo feedback: own line sf, other line sx]
aL = AP(AP(S1; g, 190); g, 445) ;   aR = AP(AP(S3; g, 265); g, 370)                  [voice diffusion g]
ηL' = vd·aL + vp·S1 + S4 ;          ηR' = vd·aR + vp·S3 + S2
# voices to the output and to the tank
eL = Σ eL_j·t_j ;  eR = Σ eR_j·t_j                                   [echo voice levels/pans]
sendL[n] = gs·eL ;  sendR[n] = gs·eR                                 [send gain gs (Rvb In)]
oL2 = fv·eL + fr·retL[n−1] ;  oR2 = fv·eR + fr·retR[n−1]             [FX mix]
width w1/w2, output 4·go, Mix as in 13.2
```

AP is the allpass of 2 with the delays shown (fixed). The glide taps take the positions of the program data (17.2; static) and
glide with ρ toward new targets (14.3); several programs move them with modulators (flanger, Doppler "car-by"
effects). The stored values φ, η, the line inputs and the allpass states saturate at ±1.

### 13.4 Res1>Plate and Res2>Plate

Tank: Plate (5.2), fed by send[n−1]. Both algorithms share this section; they differ only in how the six resonator
pitches are assigned (14.4): **Res2** harmonises the Assign note (Key, Scale, Root, Rule and one interval per voice)
and retunes whenever one of these changes; **Res1** hands each new Assign value to the next resonator in turn
(round robin). In a Res1 program at rest no resonator is assigned, because the first Assign after a program load is
ignored: every voice has the word 6 011 780 (C0) and feedback 0, the rest state in the program data (17.2).

```
# input lines (r-sample resolution, 11.1), three taps each; feedback of the taps of the previous sample
xL = aLL·inL + aLR·inR + ψ1·t1 + ψ2·t2 + ψ3·t3 ;   xR = aRL·inL + aRR·inR + ψ4·t4 + ψ5·t5 + ψ6·t6
write xL into input line L, xR into input line R
t_k = input line L (k = 1..3) / R (k = 4..6) delayed by D_k            [voice delays]
# resonator k = 1..6: loop delay N_k = ⌊P_k/2048⌋ − 2, fraction f_k = (P_k mod 2048)/2048  [tuning word P_k]
v_k[n] = v_k[n−1] + a_k·(u_k[n−1] − v_k[n−1] + t_k[n]/2)               [loop low-pass a_k]
d_k[n] = v_k[n − 1 − N_k]
u_k[n] = −b_k·e_k[n]·(f_k·d_k[n] + (1 − f_k)·d_k[n−1])                 [loop feedback b_k; envelope e_k]
# output
rL = Σ rL_k·v_k[n] ;  rR = Σ rR_k·v_k[n]                              [resonator levels/pans]
sendL[n] = gs·rL ;  sendR[n] = gs·rR                                  [the tank is fed by the resonators]
oL = fv·rL + fr·retL[n−1] ;  oR = fv·rR + fr·retR[n−1]                [FX mix]
width w1/w2, output 4·go, Mix as in 13.2
```

**Pitch.** The nominal resonator period is P_k/2048 samples (P_k: 14.4, program data 17.2), but the loop also contains
the low-pass, whose delay is not compensated: the loop resonates where the phase of
`b·(f + (1 − f)·z⁻¹)·z^−(N+1)·a·z⁻¹/(1 − (1 − a)·z⁻¹)` (z = e^jω) is a multiple of 2π, i.e. below 48 000·2048/P_k Hz by
about the low-pass group delay (≈ (1 − a)/a samples; with a = 0.23, C3 sounds 12 cents and C4 33 cents flat). With
b > 0 the resonator rings at that frequency and its harmonics, with b < 0 one octave lower with odd harmonics only.

e_k is −1 in steady state
(full feedback). After a program load or a retune it ramps linearly from 0 to −1 in 256 samples; on a retune at
full feedback it first ramps to 0 (256 samples), the new delay is installed, and it ramps back. The stored values
v_k, u_k, the line inputs saturate at ±1.

### 13.5 Infinite

The Infinite algorithm uses the reverb input/output section of 3 unchanged (send, post-delays, four delay voices,
output section) with its own tank program (15.1), a Chamber-like tank with an **Infinite** switch (14.5). Its only
factory program, AutoInfinite, switches Infinite with an AR envelope: while the input is above the trigger level the
tank runs normally; when the input stops, the tank freezes and sustains indefinitely.

## 14. Parameter laws

The laws of 4 apply to every parameter that has the same name and function (tank parameters, FX Mix, FX Width,
FX Adjust, High Cut, Rvb In of the reverb algorithms, glide, patches 4.10, tempo locks 4.11 with the additions of 6).
The laws below are those of the effect sections (exact where marked).

### 14.1 Common input laws (all effect algorithms) — exact

* **Input level and pan** (InLvl 0–160, InPan 0–100, 0 = left): an input reaches the left send with
  `rmul(L16(InLvl), P16(100 − InPan))` and the right send with `rmul(L16(InLvl), P16(InPan))` (exact for every
  value; word rule of 1).
* **Voice levels with a level master** (M-Band bands, Glide echo voices, resonators; Level ℓ 0–160, master M 0–80,
  40 = neutral): the level moves away from "off" by M − 40: ℓ′ = clamp(ℓ + M − 40, 80, 160) for ℓ > 80,
  clamp(ℓ − (M − 40), 0, 80) for ℓ < 80; gain = L16(ℓ′) (M-Band: × 0.125, see 14.2). A result of exactly 0 is not
  written (6.4). **Exact** (M-Band; Glide and Res within 1/32768).
* **Equal-power mixes** (FX Mix, PstMix, and the Mix program controller; 0–100): the second path (FX Mix: reverb fr;
  PstMix: delayed pd; Mix: wet) gets n = S[v] and the first path (voices fv, direct pr, dry) n = S[100 − v], with the
  101-entry table S[v] = ⌊32767·sin(π·v/200)⌋ (rounded down), word rule of 1. **Exact**.
  The plain
  cos/sin of 4.5 differ from it by up to 4·10⁻⁵.
* **Pans with a pan master** (0–100, 0 = left; master 50 = neutral): p′ = clamp(p + M − 50, 0, 100), left
  P16(100 − p′), right P16(p′). **Exact.**
* **Voice delays with a delay master** (ms; master M 0–200 %; range T = the longest voice time of the memory, 11.1:
  M-Band 43 690 / 2 730 ms, Glide 42 325 / 2 645 ms, Res 38 229 / 2 389 ms for the extended / standard memory):

  ```
  t0 = min(v, T)                                  plain value v (ms)
  t0 = ⌊1000·S/fs⌋, then t0 ← ⌊t0/2⌋ while t0 > T    tempo-locked value (num/den beats; S of 4.11 / 16.3)
  t  = min(⌊t0·M/100⌋, T)                         delay master, clamped once more
  n  = ⌊fs·t/(1000·r)⌋ = ⌊48·t/r⌋,  D = r·n          11.1 (44.1 kHz: 14.7)
  ```

  **Exact** (M-Band, Glide and Res, both memories, 48 and 44.1 kHz). So the delay master also scales tempo-locked voices, and the
  "halve until it fits" step comes before the master and works on milliseconds against the range of the memory in
  use (RezClimber: 4 beats at 81 BPM = 2 968 ms → 1 484 ms with the standard memory, then Master 7 % → 103 ms; with
  the extended memory 2 968 → 207 ms). Six Trips' 5/4 beat at 94 BPM is S = 38 304 samples → 798 ms. A patch into a
  locked voice: 16.6.

### 14.2 M-Band+Rvb

| quantity | law |
|---|---|
| band input level g_j | 0.125 · level law of 14.1: n = rmul(L16(ℓ′), 4096). **Exact** |
| band delay D_j | 14.1, 11.1 (range 0 … 43 690 ms; the standard memory clamps at 2 730 ms). **Exact** |
| HiCut / LoCut (0–121) with Mstr HC / Mstr LC (0–240, 121 = neutral) | table index i = clamp(V + M − 121, 0, 121); HiCut i selects the triple (d1 − 1, f1 − 1, c) and LoCut i the triple (e, d2 − 1, f2 − 1) of the HiCut / LoCut tables (24-bit words, value = word/2^23, signed; 122 triples each, in the program data: `python effect_programs.py table hicut`, `… table locut`). **Exact** |
| band feedback fb_j (Feedback 0–200, master 0–100) | n = rmul(rmul(F16(f), P16(M)), 27852) with F16(f) = ⌊32767·(f − 100)/100⌋, i.e. fb ≈ 0.85·(f − 100)/100·M/100 (fb into the same side's send: bands 1–3 left, 4–6 right). **Exact** for the master; voice side within 1/32768 |
| band pans pL_j, pR_j | 14.1. **Exact** |
| send gain gs | fixed 0.4375 (Rvb In is not stored in M-Band programs) |
| FX Mix, FX Width, FX Adjust, tank parameters | 4.5, 4.2 … (Chamber tank) |
| Clear | as 4.13 for the band levels and feedback |

HiCut table (`table hicut`): word1 = d1 − 1, word2 = f1 − 1, word3 = c; index 121 = open (f1 ≈ 1.87, about 18 kHz),
0: f1 ≈ 0.003, about 24 Hz. LoCut table (`table locut`): word1 = e, word2 = d2 − 1, word3 = f2 − 1; index 0 = open
(corner about 21 Hz), 121: about 16.6 kHz.

### 14.3 Glide>Hall

| quantity | law |
|---|---|
| glide tap position (A Dly L, B Dly L, A Dly R, B Dly R; 0–420 in 0.1 ms) | target P = 1 + ⌊256·fs·v/10 000⌋/256 samples = 1 + ⌊1228.8·v⌋/256 at 48 kHz (≈ 1 + 4.8·v); v = 0 → 1 sample. **Exact** (every v, 14.7). ρ = 0x400/2^23 (time constant 8 192 samples) |
| tap form c (12) and glide rate ρ | no parameter sets them: all 18 factory programs have cL = cR = −1 (linear interpolation) and ρ = 0x400/2^23 for all four taps (program data, 17.2) |
| glide tap levels (A Lvl, B Lvl; 0–160) | L16 (input-level law, 1.0 at 160). **Exact** |
| glide feedback (Fbk L/R, X-Fbk L/R; 0–200) | n = ±rmul(P16(\|v − 100\|), 16383), sign of v − 100, word rule of 1; ≈ (v − 100)/200 (± 0.5). **Exact** (every v; all 72 values of the program data) |
| Gld Lvl (0–80) | T[v] (0 = off, 80 = 1.0). **Exact** |
| echo voice delays (Voice1–6, ms, Delay Master %) | 14.1, 11.1 |
| echo voice levels / pans (Levels, Panning pages with masters) | 14.1 |
| echo feedback (V Fbk, V X-Fbk 0–200; Mstr Fbk, Mstr XFbk 0–100) | n = ±rmul(rmul(P16(\|f − 100\|), P16(M)), 32767), sign of f − 100; sf_j (own line) from V Fbk and Mstr Fbk, sx_j (other line) from V X-Fbk and Mstr XFbk; word 256·n + (n mod 256) (11), i.e. ≈ (f − 100)/100·M/100; f = 200, M = 100 gives 32 765/32 768 = 0.99994 (−0.0005 dB per pass). **Exact** (all values) |
| Voice Dif (0–100) | g = 0.328125·v/100 (21/64 at 100); diffused path on (vd = 1, vp = 0) for v > 0 |
| Rvb In (send gain gs) | rmul(T[v], 14336) = 0.4375·10^((v − 80)/20). **Exact** |
| tank: Mid Rt, Size, Link, Def, Depth, Chorus, AP Sweep … | Concert Hall laws 4, except the loop time K of 4.2: Link on K = ⌊b·M/u⌋/2 = 225 ms for every Size (half the Concert Hall's 450 ms); Link off K = ⌊⌊b·(16 + Size)/16⌋/u⌋, the Concert Hall's own value (22 … 450 ms). **Exact** |

### 14.4 Res1>Plate and Res2>Plate

| quantity | law |
|---|---|
| resonator levels / pans (V Lvl with Mstr Lvl; Panning) | 14.1 (within 1/32768) |
| voice feedback ψ_k into the input lines (V Fbk 0–200, Mstr Fbk 0–100) | ≈ (f − 100)/100 · M/100, within 1/32768 |
| loop low-pass a_k (V HiCut 0–121, Mstr HC 0–240) | a = H[clamp(V + M − 121, 0, 121)], H = the High Cut table of 4.6 (index 121 = 8 388 607, the loop low-pass open). **Exact** |
| loop feedback b_k (V Res 0–200, Mstr Res 0–100) | n = rmul(rmul(RES[M], ±RES[abs(V − 100)]), 32766), b = n/32768; the sign is that of V − 100; RES = the resonance table (101 entries, 0 … 32767; `table resonance`). A muted voice (below) gets no feedback word. A result of 0 (Mstr Res 0, V Res 100) is not written: the voice keeps its previous feedback (6.4). **Exact** |
| voice delays D_k | 14.1, 11.1 (range 0 … 38 229 ms extended, 2 389 ms standard) |
| tuning word P of note n (0 = C0, 57 = A4; Tuning 0–200 = 430.0 … 450.0 Hz) | q = 10·⌊fs·2^15/(4300 + Tuning)⌋, w = ⌊(q·T[(n − 57) mod 12] + 2^14)/2^15⌋, P = w shifted right by 4 + ⌊(n − 57)/12⌋ octaves (left for negative counts); T = 2^(−s/12) in Q15 (s = semitones above A: 32768, 30928, 29192, 27554, 26007, 24548, 23170, 21870, 20642, 19483, 18390, 17358). A at Tuning 100 is exact (440 Hz: P = 223 418); other notes lie within 40 ppm of equal temperament. **Exact** (every Tuning, every note) |

**Pitch assignment** (control step, 20 ms; **exact**):

* **Notes.** A voice's note is limited to 0 … 115: below 0 it is 0, above 115 it drops by octaves. A voice whose note
  comes out as 0 (C0) is not retuned (it keeps its previous pitch). A Res2 voice that has never been tuned since the
  power-up has the tuning word 0 (a loop of −2 samples, not implementable): a simulation
  treats it as silent. No factory program reaches it at rest: the patch sources without value (MIDI, 16.1) leave Assign
  at its stored value (MIDIchords: Assign 48, the chord on C4).
* **Res2.** Voice k with interval index iv_k (0 … 70; 35 = unison, 7 steps per octave) gets the note
  `n_k = Assign + 12·(⌊iv_k/7⌋ − 5) + I[Scale][Root][Rule][iv_k mod 7][(Assign − Key) mod 12]` (I = 0 for
  iv_k mod 7 = 0), where I is the interval table (2 scales × 7 roots × 4 rules × 6 scale degrees × 12 pitch classes;
  `table res2_intervals`): the semitone distance of the scale degree above the Assign note, with out-of-key notes
  handled by Rule (Round Down, Round Up, Shift Down, Shift Up) as the parameter descriptions state. Every change of
  Assign, Key, Scale, Root, Rule or an interval recomputes all voices. Voices k ≤ Active sound; voices above Active are
  muted.
* **Muted voices:** the tuning word is still updated, but the
  loop-feedback word b_k of a muted voice is never written; nothing else changes (input tap, levels, pans). After a program
  load that word is 0, so a muted voice does not resonate, but its delayed
  input still reaches the output through the loop low-pass (it is quiet, not silent). The program data shows b_k = 0
  for every muted voice (MajorMinor and MIDIsustain at rest, all Res1 programs at rest); a simulation uses b_k = 0.
* **Res1.** Each new value of Assign (from an edit or a patch) is given to one voice: a pointer advances by one and
  wraps at Active (after a load it points at the last voice, so the first note goes to voice 1), that voice takes the
  note and is unmuted. With Unison on, the note goes to every voice up to Active. With Active 0, Assign is ignored. The
  first Assign after a program load is ignored (the stored value does not sound). Lowering Active mutes the voices
  above it; they stay muted, also after Active rises again, until they receive a new note.
* **Transposition.** None.
* **Sample rate.** fs in the tuning-word law is the sample rate: 48 000, or 44 100 at 44.1 kHz (exact, 14.7), so the
  resonators keep their pitch at both rates.

### 14.5 Infinite

The Infinite tank follows the Chamber laws of 4 for Size, Mid Rt, Low Rt, Crossover, Rt HC, Diffusion, Pre Delay,
Spread, Shape, Link, Rvb Out, RefDly/RefLvl, EkoDly/EkoFbk (the steps of 15.1 named in its coefficient table).

**Size — exact:** the 61 Size-dependent tank offsets follow the
Chamber rule of 4.1, `A(v) = A0 + ⌊(A0 − 33)·v/16⌋`, with these Size-0 values:

```
A0: 420          A2: 1141         A3: 291          A5: 420          A7: 1369         A9: 1991
A10: 1142        A12: 1369        A15: 2276        A17: 2719        A18: 1991        A20: 2276
A34: 2988        A37: 2719        A38: 3843        A39: 2988        A41: 4069        A44: 3843
A45: 4693        A46: 4069        A50: 5121        A51: 928         A52: 614         A53: 4693
A54: 496         A55: 5121        A56: 1307        A57: 1637        A58: 1990        A59: 2374
A60: 2778        A61: 3193        A62: 3567        A63: 3991        A64: 4357        A65: 4803
A66: 5168        A67: 5447        A73: 324         A75: 1047        A76: 1445        A77: 1876
A78: 2052        A79: 2696        A80: 3009        A81: 3383        A82: 3781        A83: 4092
A84: 4523        A85: 4899        A86: 5179        A87: 324         A92: 78          A98: 152
A100: 190        A101: 78         A112: 152        A113: 279        A114: 257        A116: 190
A117: 5543
```

and the seven input-line taps follow the Chamber Spread rule of 4.8 (`q = ⌊Spread·(16 + Size)/160⌋`,
`A = PD + ⌊M·q/255⌋`; Size = 144 with Link off) with M = 4497 (A6), 13488 (A14), 16444 (A25), 2390 (A40), 8154 (A49),
878 (A121), 20736 (A122). A68 = 3563 does not depend on Size. AP Sweep changes k0, k2 and the step-87 tap (6.1
not studied for this tank; the factory program has AP Sweep 0).

**Decay — exact:** the five tank-loop gains k4, k12,
k30, k38, k123 are the decay gain g of 4.2 with the Chamber constants (L = 5252, M = 10, u = 6; e.g. Size 106, Link
off, Mid Rt 60: K = 138 ms, g = 241/256); k47 = 0.5 − g for g ≥ 1/2 and 0.5 − g − 1/256 for
g < 1/2; k46 = −0.5 at every setting; k32 and k125 are the Low Rt terms of 4.3 with this g.

**Infinite = 1** (the freeze) changes 20 tank steps; 15.1 lists every one of them as a symbol, with its value in both
states. In words: the AP Sweep allpass pair k0/k2 and the decay gains become 0 / 1 (the loop neither diffuses at the
first allpass nor loses energy), the damping (Rt HC: k20, k117 → 1, k22, k119 → 0) and the low-band terms (k32, k125
→ 0) are switched off, k24 and k121 become 0.75, k47 becomes −0.5, three output taps move (k54 → 0; k68 → 0.171875 with
A68 = 3563; k73 → 0.5 and k87 → 0 with A73 = A87 = 3002), and the send gain gs of the input section becomes 0 (no new
input). Infinite = 0 restores the values the laws give. The coefficients change with the shared coefficient smoothing
of 6.4. FlyTapRng and FlyTapDly are not stored in programs (no
effect in the factory program).

**Moving tap — no parameter law.** The coefficient k87 and the offsets A73 and
A87 of the running tank (Infinite = 0) are not set by a parameter: they belong to a modulated tap that moves while
the program runs (a tap list containing steps 73 and 87 with a period of 121, like the Concert Hall's modulated
taps of 6.2). Their values change from one recall to the
next with every edit (k87 0.01 … 0.20, A73 / A87 among a few values such as 2 102, 2 251, 2 614 at Size 106) without
following Mid Rt, Low Rt or any other parameter: each recall captures the motion at another moment. 15.1 and the program
data give the factory snapshot (k87 = 0.08203125, A73 = A87 = 2 251); a simulation keeps these values, also after edits
of the tank parameters. The motion itself is not described. Infinite = 1 replaces the tap by the fixed
values of 15.1.

**Notation of 15.1.** Each `k<s>` is the whole coefficient of step s, including its integer part on
some steps, so one set of equations covers both switch states. "Set by" names the parameter whose law gives the value;
"Infinite only" = a constant of the running tank that only the
switch changes.

### 14.6 Parameter ranges

Every parameter runs from 0 to a maximum; the program data (17.2) gives it as
`max` for every stored value, Mod-page value and parameter in effect (4.14 lists the reverb algorithms and the Mod
page). A patch result is clamped to 0 … max (4.10). For the effect algorithms:

| algorithm | parameters (max) |
|---|---|
| all three | InLvl L/R 160, InPan L/R 100, FX Mix 100, FX Width 720, FX Adjust 80, Rvb In 80; tank parameters as 4.14 (Chamber / Concert Hall / Plate) |
| M-Band+Rvb | DelayTime Master 200, Voice1–6 43 690 ms; Levels Master 80, Voice 160; Feedback Master 100, Voice 200; Filters Mstr HC / Mstr LC 240, V HiCut / V LoCut 121; Panning Master, Voice 100 |
| Glide>Hall | A/B Lvl 160, Fbk / X-Fbk 200, A/B Dly 420 (42 ms), Gld Lvl 80; DelayTime Master 200, Voice1–6 42 325 ms; Levels Master 80, Voice 160; Mstr Fbk / Mstr XFbk 100, V Fbk / V X-Fbk 200; Voice Dif 100; Panning 100 |
| Res1/Res2>Plate | DelayTime Master 200, Voice1–6 38 229 ms; Mstr Lvl 80, V Lvl 160; Mstr Fbk 100, V Fbk 200; Mstr HC 240, V HiCut 121; Mstr Res 100, V Res 200; Tuning 200, Active 6, Assign 127; Res1: Unison 1; Res2: Scale 1, Key 11, Root 6, Rule 3, Voice1–6 70; Panning 100; Plate tank with Spread 100, Attack 100 |
| Infinite | the Chamber list of 4.14 (tank, input/output section, voices, PstDly, GldResp / GldRange) plus FlyTapRng 31, FlyTapDly 5, Infinite 1 |

The voice-delay maxima are those of the extended memory; with the standard memory the delay law clamps them at 2 730 /
2 645 / 2 389 ms (14.1). A user interface should offer the range of the memory it simulates.

### 14.7 Operation at 44.1 kHz — exact

At 44.1 kHz the effect algorithms change exactly the following. The tanks, the reverb input/output section, the delay voices, the
output High Cut and Chorus+Rvb follow 4.15. Every quantity not named here or in 4.15 keeps its 48 kHz value (levels,
pans, feedback words, the resonance table, FX Mix / Width, the glide feedback and levels).

| quantity | at fs = 44 100 |
|---|---|
| effect voice delays (M-Band, Glide, Res; 14.1, 11.1) | times kept in ms, same ranges T in ms: n = ⌊fs·t/(1000·r)⌋, D = r·n (t = 999 ms → 44 054 samples, standard memory). Tempo-locked: S with fs = 44 100 (4.11), so they keep their time. **Exact** (both memories: tempo 40 … 250 BPM, every lo/hi of 16.3, plain values across both ranges) |
| M-Band band filters (14.2) | own tables, `table hicut44` / `table locut44` (same layout and index law as hicut / locut). The corners stay within about 1.5 % of the 48 kHz ones up to index 110; index 121 (open) lies at about 19.3 kHz instead of 18.4 kHz. **Exact** (all 122 indices of each) |
| glide tap positions (14.3) | P = 1 + ⌊256·fs·v/10 000⌋/256 samples (v in 0.1 ms; at 48 kHz this is the 1 + 4.8·v of 14.3, now exact). ρ is unchanged (in samples: the glide takes 8.8 % longer). **Exact** (v = 0 … 420 at both rates) |
| Glide>Hall loop time K (14.3) | the law of 14.3 with b = 49 (4.15): Link on 245 ms (48 kHz: 225 ms), Link off 24 … 490 ms. **Exact** (Size 0 and 304) |
| resonator tuning word (14.4) | the law of 14.4 with fs = 44 100: the resonators keep their pitch. **Exact** (Tuning 0 … 200, Assign 0 … 127) |
| resonator loop low-pass (14.4) | a = H44[i] (the 44.1 kHz High Cut table of 4.15) with the same index i. **Exact** |
| Infinite | 4.15 (tank decay with f = 44, output High Cut H44, delay voices in ms) |

The modulators keep their 50 Hz control rate (4.15).

## 15. Tank programs

M-Band+Rvb runs the Chamber tank (5.1), Glide>Hall the Concert Hall tank (5.3), Chorus+Rvb and Res1/Res2>Plate the
Plate tank (5.2), unchanged: in all 165 programs every tank coefficient and offset not listed as a symbol in 5 has
exactly the value 5 shows. Section 17 lists the symbols' values per program.

### 15.1 Infinite tank — exact program

### Program

Delay memory and notation as in section 5. Three circular buffers of 65 536 samples, **M0** (tank), **M1** (right input line), **M2** (left input line); `R s` reads, `W s` writes, X′ = the value the same read produced one sample earlier, `A<s>` program-dependent offsets (values in the program data, 17.2, and below). Every coefficient that a tank parameter or the Infinite switch changes is kept as a symbol `k<s>`; all other numbers are fixed. Per sample n, execute in order:

```
R0:M0[n-A0]  R2:M0[n-A2]  W3:M0[n-A3]  W5:M0[n-A5]  R6:M1[n-A6]  R7:M0[n-A7]
R9:M0[n-A9]  W10:M0[n-A10]  W12:M0[n-A12]  R14:M1[n-A14]  R15:M0[n-A15]  R17:M0[n-A17]
W18:M0[n-A18]  W20:M0[n-A20]  R22:M0[n-6]  R23:M0[n-6]  R24:M2[n-A24]  R25:M2[n-A25]
R26:M0[n-2]  W27:M0[n-5]  R28:M0[n-6]  W31:M0[n-1]  R34:M0[n-A34]  W37:M0[n-A37]
R38:M0[n-A38]  W39:M0[n-A39]  R40:M2[n-A40]  R41:M0[n-A41]  W44:M0[n-A44]  R45:M0[n-A45]
W46:M0[n-A46]  R49:M2[n-A49]  R50:M0[n-A50]  R51:M0[n-A51]  R52:M0[n-A52]  W53:M0[n-A53]
R54:M0[n-A54]  W55:M0[n-A55]  R56:M0[n-A56]  R57:M0[n-A57]  R58:M0[n-A58]  R59:M0[n-A59]
R60:M0[n-A60]  R61:M0[n-A61]  R62:M0[n-A62]  R63:M0[n-A63]  R64:M0[n-A64]  R65:M0[n-A65]
R66:M0[n-A66]  R67:M0[n-A67]  R68:M0[n-A68]  R70:M1[n-A70]  R71:M1[n-A71]  R73:M0[n-A73]
R75:M0[n-A75]  R76:M0[n-A76]  R77:M0[n-A77]  R78:M0[n-A78]  R79:M0[n-A79]  R80:M0[n-A80]
R81:M0[n-A81]  R82:M0[n-A82]  R83:M0[n-A83]  R84:M0[n-A84]  R85:M0[n-A85]  R86:M0[n-A86]
R87:M0[n-A87]  R89:M2[n-A89]  R90:M2[n-A90]  R92:M0[n-A92]  R94:M1[n-A94]  W97:M0[n-33]
R98:M0[n-A98]  R100:M0[n-A100]  W101:M0[n-A101]  R102:M0[n-4]  W103:M1[n-0]  R104:M0[n-8]
W107:M0[n-3]  R109:M2[n-A109]  W112:M0[n-A112]  R113:M0[n-A113]  R114:M0[n-A114]  R115:M0[n-8]
W116:M0[n-A116]  R117:M0[n-A117]  W118:M2[n-0]  R119:M0[n-8]  R121:M1[n-A121]  R122:M1[n-A122]
W124:M0[n-7]
```

Values written (w_s) and tank outputs:

```
w3   (M0) = X0*k0 + X102′*k102*k125 + X104′*k104*k125 + X115′*k123
w5   (M0) = X0*k0*k2 + X0 + X102′*k102*k125*k2 + X104′*k104*k125*k2 + X115′*k123*k2
w10  (M0) = X2*k4 + X6*k6
w12  (M0) = X2*k4*k9 + X6*k6*k9 + X7
w18  (M0) = X14*k14 + X9*k12
w20  (M0) = X15
w27  (M0) = X17*k20 + X22*k22 + X24*k24 + X25*k25
w31  (M0) = X26*k26 + X28*k28
w37  (M0) = X23*k30 + X26*k26*k32 + X28*k28*k32
w39  (M0) = X34
w44  (M0) = X38*k38 + X40*k40
w46  (M0) = X41
w53  (M0) = X45*k46 + X45*k47 + X49*k49
w55  (M0) = X50
w97  (M0) = X92*k92 + X94*k94 + sendR
w101 (M0) = X92*k92*k96 + X92 + X94*k94*k96 + X98*k98 + k96*sendR
w103 (M1) = X92*k100*k92*k96 + X92*k100 + X94*k100*k94*k96 + X98*k100*k98 + X98 + k100*k96*sendR
w107 (M0) = X102*k102 + X104*k104
w112 (M0) = X100*k107 + X109*k109 + sendL
w116 (M0) = X100*k107*k111 + X100 + X109*k109*k111 + X113*k113 + k111*sendL
w118 (M2) = X100*k107*k111*k115 + X100*k115 + X109*k109*k111*k115 + X113*k113*k115 + X113 + k111*k115*sendL
w124 (M0) = X117*k117 + X119*k119 + X121*k121 + X122*k122
retL      = 7*X52*k91/16 + X73*k73*k91 + 3*X75*k91/8 + X76*k91/2 - 3*X77*k91/8 + X78*k91/2 - 3*X79*k91/8 + X80*k91/2 - 7*X81*k91/16 + 3*X82*k91/8 - X83*k91/2 + 3*X84*k91/8 + X85*k91/2 - 3*X86*k91/8 + X87*k87*k91 + X89*k89*k91 + X90*k90
retR      = -3*X51*k72/8 + X54*k54*k72 - X56*k72/2 - 7*X57*k72/16 - 3*X58*k72/8 + X59*k72/2 - X60*k72/2 - 7*X61*k72/16 + 3*X62*k72/8 + X63*k72/2 + 3*X64*k72/8 + X65*k72/2 - 7*X66*k72/16 + 3*X67*k72/8 + X68*k68*k72 + X70*k70*k72 + X71*k71
```

Coefficients `k<s>`: what sets them, and their values in AutoInfinite (stored values, patches off) with Infinite = 0 (running) and Infinite = 1 (frozen). Values are the whole coefficients, integer part included, so the same equations hold in both states. Bold = changed by the switch.

| k | set by | Infinite 0 | Infinite 1 |
|---|---|---|---|
| k0 | AP Sweep pair (6.1) | −0.5625 | **0** |
| k2 | AP Sweep pair (6.1) | 0.5625 | **0** |
| k4 | decay gain g: Mid Rt, Link, Size with Link off (4.2 with the Chamber constants) | 0.921875 | **1** |
| k6 | Shape (4.7) | 0.0625 | 0.0625 |
| k9 | Diffusion (4.5); Shape (4.7) | 0 | 0 |
| k12 | decay gain g: Mid Rt, Link, Size with Link off (4.2 with the Chamber constants) | 0.921875 | **1** |
| k14 | Shape (4.7) | 0 | 0 |
| k20 | Rt HC (4.4) | 0.765625 | **1** |
| k22 | Rt HC (4.4) | 0.234375 | **0** |
| k24 | Rt HC (4.4); Shape (4.7) | 0.5625 | **0.75** |
| k25 | Shape (4.7) | 0.0625 | 0.0625 |
| k26 | Crossover (4.4) | 0.97265625 | 0.97265625 |
| k28 | Crossover (4.4) | 0.02734375 | 0.02734375 |
| k30 | decay gain g: Mid Rt, Link, Size with Link off (4.2 with the Chamber constants) | 0.921875 | **1** |
| k32 | Low Rt with the decay gain g (4.3) | 0.03515625 | **0** |
| k38 | decay gain g: Mid Rt, Link, Size with Link off (4.2 with the Chamber constants) | 0.921875 | **1** |
| k40 | Shape (4.7) | 0.0625 | 0.0625 |
| k46 | fixed −0.5 (written with the decay gains) | −0.5 | −0.5 |
| k47 | 0.5 − g, or 0.5 − g − 1/256 for g < 1/2 (Mid Rt, Link, Size with Link off; 14.5) | −0.421875 | **−0.5** |
| k49 | Shape (4.7) | 0 | 0 |
| k54 | Infinite only | 0.5 | **0** |
| k68 | Infinite only | 0 | **0.171875** |
| k70 | Shape (4.7) | 0 | 0 |
| k71 | RefLvl R (4.5) | 0 | 0 |
| k72 | Rvb Out (4.5) | 0.9375 | 0.9375 |
| k73 | Infinite only | 0 | **0.5** |
| k87 | moving tap of the algorithm's own modulation, with A73 / A87 (14.5): factory snapshot | 0.08203125 | **0** |
| k89 | Shape (4.7) | 0 | 0 |
| k90 | RefLvl L (4.5) | 0 | 0 |
| k91 | Rvb Out (4.5) | 0.9375 | 0.9375 |
| k92 | Diffusion (4.5) | −0.56640625 | −0.56640625 |
| k94 | EkoFbk R (4.5) | 0 | 0 |
| k96 | Diffusion (4.5) | 0.56640625 | 0.56640625 |
| k98 | Diffusion (4.5) | −0.4921875 | −0.4921875 |
| k100 | Diffusion (4.5) | 0.49609375 | 0.49609375 |
| k102 | Crossover (4.4) | 0.97265625 | 0.97265625 |
| k104 | Crossover (4.4) | 0.02734375 | 0.02734375 |
| k107 | Diffusion (4.5) | −0.56640625 | −0.56640625 |
| k109 | EkoFbk L (4.5) | 0 | 0 |
| k111 | Diffusion (4.5) | 0.56640625 | 0.56640625 |
| k113 | Diffusion (4.5) | −0.4921875 | −0.4921875 |
| k115 | Diffusion (4.5) | 0.49609375 | 0.49609375 |
| k117 | Rt HC (4.4) | 0.765625 | **1** |
| k119 | Rt HC (4.4) | 0.234375 | **0** |
| k121 | Rt HC (4.4); Shape (4.7) | 0.5625 | **0.75** |
| k122 | Shape (4.7) | 0.0625 | 0.0625 |
| k123 | decay gain g: Mid Rt, Link, Size with Link off (4.2 with the Chamber constants) | 0.921875 | **1** |
| k125 | Low Rt with the decay gain g (4.3) | 0.03515625 | **0** |

Offsets that the switch changes (the other `A<s>` follow 14.5 and the laws of 4):

| A | Infinite 0 | Infinite 1 |
|---|---|---|
| A68 | 0 | 3563 |
| A73 | 2251 | 3002 |
| A87 | 2251 | 3002 |


## 16. Modulation and tempo

Everything in 6 applies. Additions and corrections:

### 16.1 Patch sources (completes 6.3)

The MIDI sources of 6.3:

| source | value |
|---|---|
| controllers 1–31 and 33–119 | the controller value (controller 32 is not a source) |
| P Bend | pitch-bend MSB |
| A Touch | aftertouch |
| Velocity | velocity of the last note-on |
| Last Note | note number of the last note-on |
| Low Note / High Note | lowest / highest note held |
| Clk Cmnds | MIDI Start/Continue → 127, Stop → 0 |

A patch slot can also have no source ("none").

The MIDI sources start **without a value** after power-up: a patch from them contributes nothing (its destination
keeps the stored value) until a message arrives. Factory programs use Velocity and High Note (MIDIchords), Last Note
(MIDIsustain). The footpedal input behaves as "no value" or full value when no pedal is connected (Pedal Swell,
patched Footpedal → input level, passes the input at full level without a pedal).

### 16.2 AR envelope, all modes (6.3)

Modes 1 (One Shot) and 3 (Repeat) behave as 6.3 states (examples: OneShotCarBy, AmbientCarBy; AR Attack 500, Release 0): One Shot attacks once (10.0 s), holds while triggered,
releases when the input stops and attacks again on the next trigger; Repeat cycles attack → full release while
triggered, an attack in progress continues through a 1.5 s input gap. The new envelope value reaches the signal one control step later, as in Mode 2.

### 16.3 Tempo (completes 4.11)

* **Program tempo.** A program stores its tempo as a 16-bit word, BPM = word + 40 (40 … 400). Four factory programs are stored above 250 BPM: Dly>EQ>Pans 333 BPM, SweptNotches, Env Notches and
  BandReject 4 300 BPM. Only Dly>EQ>Pans has
  tempo-locked values: its "2 and 3 beats" are 360 and 540 ms at 333 BPM (equally 1/2 and 3/4 beat of
  83.25 BPM — the unit counts the beats of the stored tempo).
* **Tempo representation.** The tempo is held as a beat period d in 1/16 ms. A program tempo of BPM (integer)
  gives d = ⌊960 000/BPM⌋; a tapped tempo sets d from the measured interval (fractional tempos are possible).
  From d: R = ⌊4 096 000/d⌋ (a 16-bit rate), displayed BPM = ⌊(60·R + 128)/256⌋. A MIDI-clock tempo sets R
  directly (below); everything tempo-locked depends on R only.
* **Beat length in samples.** B = ⌊12 288 000/R⌋. **Exact**: it reproduces the beat table of 4.11 for every tempo
  40 … 250 BPM (the table is exact only where 4 096 000/d is an integer, hence its "up to 0.25 % longer" values).
* **Locked lengths — exact.** A locked value of num/den beats is S = ⌊256·fs·num/(R·den)⌋ samples with one rounding
  (shared by every tempo-locked length: voices, post-delays, tank offsets, effect voices); tank offsets
  S + 1, post-delays S, voices ⌊1000·S/fs⌋ ms (4.11), effect voices then 14.1 and 11.1. The stored code is
  v = 32·den + num (the "lo/hi" of 4.11: lo = num, hi = den), num and den 1 … 24.
* **Tempo-locked modulator rates — exact.** LFO Rate, Sw 1 Rate and Sw 2 Rate can be tempo-locked (43 factory
  programs, e.g. the Tap and Tremolo programs, Concert Wave). The stored value v = 32·hi + lo means one cycle per lo/hi
  beats (Tremolo Tap1: hi = 2, lo = 1 = two cycles per beat; lo and hi are limited to 1 … 24). The phase increment per
  control step (20 ms; phase cycle 2^17, as for the free-running rate of 6.3) is
  `inc = ⌊K·R/256⌋` with `K = ⌊131 072·hi/(50·lo)⌋` and the tempo rate R = ⌊4 096 000/d⌋ (below), i.e.
  f = inc·50/131 072 Hz ≈ (hi/lo)·BPM/60. It is recomputed whenever the tempo changes. Example: Concert Wave
  (LFO Rate 34 = 2 beats per cycle, 86 BPM): R = 366, K = 1 310, inc = 1 872, f = 0.7141 Hz.
* **MIDI-clock tempo (Tempo Source = MIDI Clock) — exact.** Each clock
  (F8) is time-stamped with the 1 ms tick counter (20 ticks = one control step). About every 500 ticks
  (the first control-loop pass ≥ 500 ticks after the previous evaluation, with n ≥ 1 clocks received since), the
  mean clock interval in 1/256 tick is x = ⌊(256·Δ + ⌊n/2⌋)/n⌋, Δ = ticks between the last clock of this window and
  the last clock of the previous one. Smoothing: y = x at the first window after clocks start (and whenever
  |x − y| ≥ 21), otherwise y = ⌊(31·y + x + 16)/32⌋ (so a change of less than about 16 is ignored). Then
  R = ⌊(2 730 667 + ⌊y/2⌋)/y⌋ (= 4 096 000·2/3/y rounded; 24 clocks = 1 beat), accepted when 171 ≤ R ≤ 1 707 and
  changed; tempo-locked values are then recomputed. When the clocks stop, R holds.
* **BeatValue** (system setting, stored with a program) only sets the rate of the MIDI clock the unit **sends**
  (Int Clock On, Source Internal): one clock every P/256 ticks, P = ⌊(256·q + 12)/24⌋, q = ⌊4·m/t⌋, m = ⌊d/16⌋
  (beat in whole ms), t = 2, 3, 4, 6, 8, 12, 16 for 1/8, dotted 1/8, 1/4, dotted 1/4, 1/2, dotted 1/2, whole. It
  does not change R, B or any rate.
* The other tempo settings (Tempo Mode Pgm/Global, Rate 40–400 BPM, tap with averaging) are system settings; a
  simulation that follows a host tempo uses d = 960 000/BPM (rounded down) and the rules above.

### 16.4 Effect-section changes while running

Patches on effect parameters take the path of 6.4: coefficient words glide with the shared smoothing; coarse taps
cross-fade when their target changes (128 / 171 samples, one tap at a time, section 12); glide taps glide with ρ (14.3);
resonators retune with the envelope fade of 13.4. A new tap target is computed in the control step (20 ms, 6.3), so a
patch that sweeps a delay (Super Ball!, LevlStutter, StereoEqEkos with the knob) moves it in 20 ms steps, each a
short cross-fade; the effect taps have no GldResp/GldRange glide (4.13 is for the reverb voices only).

**Voice gate after a program load, and Clear — exact.**
Each delay voice has a gate on its level and feedback words: while the gate is closed the words are held (their glide
rate is 0, 6.4: each keeps the value it has); when it opens they glide to the values the settings give.

* After a program load the gate of a voice opens when the time since the load reaches the voice's delay time (after
  the delay master). A voice therefore stays silent until the first input has travelled through its line — the same
  as a cleared line. Its held words are 0 after a load.
* While **Clear** = 1 every gate of a voice whose delay is not 0 stays closed (voices at 0 ms are not affected); when
  Clear returns to 0 they open again.
* The gated words are, per algorithm: the reverb algorithms' delay-voice levels and feedback gains (4.13) and the
  chorus voices'; M-Band+Rvb the band input levels g_j (the band feedback fb_j is not gated); Glide>Hall the echo voice
  levels eL_j / eR_j and the echo feedback sf_j / sx_j; Res1/Res2>Plate the resonator levels rL_k / rR_k and the input
  feedback ψ_k.

The program data (17.2) gives every program with all gates open (settled beyond its longest voice), except where Clear
is 1 at rest (TapEkoClear, MIDIsustain): there the gated words are 0. The earlier data had settled only 2.3 s and showed
0 for the voices longer than that (FinishLine voices 2 and 5, MajorMinor voices 4 and 5; corrected).

### 16.5 Latch (completes 6.3) — exact

The Latch source passes another source x (0 … 127; trigger block: source, Low L, High H) through a small
state machine, evaluated in the control step **only when x changes** (a source without value gives the Latch no
value). Output y:

* **No latching** when L ≥ H, or L = 0 and H = 127: y = x.
* **High latch** (L = 0, H < 127): following: y = x until x ≥ H, then y = H (latched). Latched: y stays H while x
  falls below H, rises to ≥ H again and falls back to ≤ H; at that last step y = x and the latch follows again.
* **Low latch** (L > 0, H = 127): mirrored: following y = x until x ≤ L, then y = L; released after x rises to ≥ L,
  falls to ≤ L and rises to ≥ L again (y = x at that step).
* **Both** (L > 0, H < 127): following y = x while L < x < H; x ≤ L latches low and x ≥ H latches high, each released
  as above, back to following.

A momentary footswitch (0 / 127) with L = 0 and H < 127 thus becomes push-on / push-off with the values H / 0
(FSw2 Freeze, LongestLoop: L = 0, H = 100, source Foot Sw 2; FSw2 Rotary and the other rotary programs: H = 33, the
Latch then triggers the AR envelope). The Latch can also be the trigger source of AR, Sw 1 and Sw 2 (Prime Blue:
Sw 1 and Sw 2 triggered by a high latch at 60 on Mono Lvl).

### 16.6 Patches into tempo-locked parameters — exact

A patch whose destination is tempo-locked (delay voices, post-delays, tank offsets, LFO / Sw rates with the lock
flag) does not use its curve of 4.10. Each patch keeps two 128-entry tables, one for a plain destination and one for a
locked destination (the destination's lock flag selects the table at run time):

* **Point code.** A patch stored while its destination was locked (patch record mode byte 1; all such factory patches)
  has points y = 32·num + den, meaning num/den beats (num = y >> 5, den = y & 31; a point with num or den 0 or above
  24 becomes 1/1). This is the reverse of the stored parameter code v = 32·den + num (16.3): the stored value 35* is
  3 beats, a patch point 35 is 1/3 beat. A patch stored for a plain destination has num = den = 1 in its locked table
  (it gives 1/1 beat if the destination is locked later), and a locked-code patch has 0 in its plain table.
* **Interpolation.** num and den are interpolated separately between the points, each in 8.8 fixed point: start
  256·num_i, add ⌊256·(num_{i+1} − num_i)/Δx⌋ (truncated toward zero) per source step, take the integer part (the same
  for den); the end values hold outside the points.
* **Several patches.** The patch words 256·num + den of all patches with a value are added (16-bit, no clamp) and
  replace the stored value, as for plain patches (6.4); when no patch has a value (6.3, 16.1) the
  destination keeps its value (the stored one after a program load).
* The result num/den then follows 16.3 and 14.1 (voices: S, halving, delay master, 11.1; LFO / Sw: one cycle per
  num/den beats).

Examples: StereoEqEkos (ADJUST "Ekos/Beat" k, points (1, 32·j + 1) … (24, 32·j + 24)
on voice j) gives voice j = j/k beat, so the stored knob 3 gives 1/3, 2/3, … 2 beats (247 ms steps at 81 BPM);
Detuned Ekos and [n]Ekos/Beat the same for 1/k beat; RotorWaa (ADJUST → LFO Rate, points (0, 257), (10, 40)) gives an
LFO period of 8 beats at knob 0, 5/3 beat at the stored knob 4 and 1/8 beat at knob 10. At rest (source value 0) the
AR-driven patches of LevlStutter put voice j at j beats (53 BPM), those of FSw2 Elevate voices 1 and 3 at 1 beat
(75 BPM); Chord Walk's Sw 1 / Sw 2 patches 4 beats (8/2), Freez 2over3's Sw 1 patch 1/2 beat. The program data (17.2)
follows this rule.

## 17. Programs

The per-program data is not printed in this document. It is in the data file **`effect-programs.json`**, read with
the query script **`effect_programs.py`** (both next to this document; Python 3.6 or later, standard library only).
17.1 describes every one of the 200 programs, the 34 programs of Part I included, so that the
programs can be reproduced one by one; 17.2 says what the data holds and how to get it.


### 17.1 Program descriptions

What each of the 200 factory programs does and how it uses its algorithm, in program order. The exact values (every
stored parameter, the resulting coefficients and delays) are in the program data (17.2) for all 200 programs; the 34
programs of Part I are marked *Part I* (section 7 keeps their group notes and the full BowedChords
description). A description names the parts of the
algorithm that are audible, the role of the reverb, what the ADJUST knob changes, and what every modulator and patch
does, including the controls a simulation has to supply from outside (MIDI, footswitches, pedal).

Conventions used below (all from 4, 6.3 and sections 14 and 16 of this document):

* **Levels** 0 … 160 with 80 = off; above 80 in phase, below 80 phase-inverted (0 = full level inverted).
  **Feedback** 0 … 200 with 100 = none, below 100 negative. **Pan** 0 = left, 100 = right.
* **FX Mix** is an equal-power crossfade: 0 = effect voices (bands, chorus voices, resonators, glide/echo voices) only,
  100 = reverb only. The program **Mix** (Mod page) is the wet/dry crossfade; 100 = wet only.
* "Beat" is a beat at the program tempo; tempo-locked values follow the tempo rules of 16.3 (a host tempo can replace
  the stored one). Delays and rates are given at the stored tempo.
* Chorus voices: the delay sweeps from the voice time up by the peak excursion (in samples) with a triangle LFO whose
  rate and direction change pseudo-randomly about every 0.75 s (3.1, 4.12).
* Modulators run at the 20 ms control rate (6.3). **AR** = attack/release envelope (Mode 1 One Shot, 2 Retrigger,
  3 Repeat), **Sw 1 / Sw 2** = time switches (Mode 1 pulse, 2 ramp), **Env L / Env R** = input envelope followers,
  **Mono Lvl** = input level (L+R), **Latch** = 16.5. "At rest" = the state of a modulator before any input arrives.
* External sources (Footpedal, Foot Sw 1/2, MIDI controllers, aftertouch, notes) have no value until something is
  supplied; a patch from them leaves its destination at the stored value until then (16.1).

#### Programs 1 – 50

**1. Prime Blue** (Chorus+Rvb, 81 BPM). A wide chorus with tempo echoes into a 2.6 s reverb. Four short voices
(9–19 ms, peak sweep 55–90 samples at 1.1–2.0 Hz, slight negative feedback) make the chorus; voices 2 and 4 are long
tempo echoes (20/19 and 1 beat, 781 / 742 ms) with +0.29 feedback. The LFO (0.33 Hz) swaps the two echoes between the
sides (V2 pan 0→100, V4 100→0). ADJUST ("Efx/Rvb X") crossfades effect and reverb (FX Mix 0 → 100) and below 32
pulls the echoes down (level 152 → 113, feedback to none). Sw 1 (ramp, 6 beats per cycle) is triggered through a
high Latch (60) on the input level and pans voices 3 and 6 across while the input stays loud; Sw 2 is set up but
drives nothing. The Footpedal narrows the output (FX Width 405 → 360). Pre-delay 1/12 beat, early reflections at 1/12
and 1/8 beat.

**2. EkoChorus** (Chorus+Rvb, 120 BPM). Chorus and two slow echoes, little reverb (1.98 s, FX mix 0.95/0.32).
Voices 1–2 (25 / 33 ms, panned 30 / 70) are a deep chorus (213-sample sweep, 0.42 / 0.34 Hz); voices 4 and 5 are
echoes at 356 ms (left, +0.21 feedback, swept at 0.76 Hz, so its repeats waver in pitch) and 408 ms (right, −0.23
feedback, unmodulated). Voices 3 and 6 are off. ADJUST sets FX Mix (0 … 100 over the knob).

**3. Wet Chorus** (Chorus+Rvb, 120 BPM). Six chorus voices at 0 ms (126-sample sweep at 0.4–1.1 Hz each, voices
1, 2, 5, 6 with +0.30 feedback, short flanged resonance), alternately panned hard left/right, into a long 6.8 s
reverb. ADJUST = FX Mix (stored 38: voices 0.83, reverb 0.56).

**4. Wet Chorus 2** (Chorus+Rvb, 120 BPM). Wet Chorus with voices 1 and 2 moved to 400 ms (+0.15 feedback, hard
left / right): a chorus plus a soft doubled slap-back echo, with the same 6.8 s reverb. ADJUST = FX Mix (stored 20).

**5. Wet Chorus 3** (Chorus+Rvb, 120 BPM). As Wet Chorus 2 with the echoes panned 25 / 75 and the 0-ms voices
3 / 4 hard left / right; ADJUST = FX Mix.

**6. Chorus Plate** (Plate, 120 BPM). A plate (3.59 s) whose four delay voices (3–9 ms, diffused, Voice Dif 51)
and post-delays are swept by the triangle LFO: the LFO moves PstDly L / R in opposite directions (16 ↔ 32), scales
all voice times through the Delay Master (20 … 100 %) and modulates its own rate (0.18 … 0.50 Hz), so the chorus
speeds up and slows down. ADJUST ("Korus Mix") sets PstMix, the share of the post-delayed (chorused) output.

**7. Chorus Room** (Chamber, 120 BPM). A chamber (Mid Rt via ADJUST "Decay": 0 … 60, stored 43 = 2.60 s) with no
delay voices (FX Mix 100, reverb only). The chorus is the left post-delay swept by the triangle LFO (9 ↔ 27); the
LFO rate follows the input level (Mono Lvl → 0.20 … 1.00 Hz), so louder playing makes a faster chorus. The AR
envelope is set up (trigger Mono Lvl ≥ 49) but drives nothing.

**8. ChorusToRvb** (Chorus+Rvb, 120 BPM). A chorus that gives way to the reverb as the input gets louder. Voices
at 0 / 43 / 74 ms each side (Delay Master 51 %), the 74-ms voices with strong negative feedback (−0.35, −0.31), wide
mixed sweeps. The input level drives FX Mix (soft 70 → loud 12, i.e. from mostly reverb to mostly chorus) and
narrows FX Width at higher levels. ADJUST ("Decay") sets Mid Rt 0 … 63 (stored 25 = 1.11 s).

**9. Funkus Room** (Chamber, 120 BPM). A 4.2 s chamber with four short delay voices (0–18 ms, panned across) that
the triangle LFO sweeps in opposite directions (V1 1 ↔ 24 ms, V3 18 ↔ 2 ms, …) together with both post-delays
(cosine on the left one): a strong flanging, warbling room. ADJUST sets the LFO rate (0.02 … 3.5 Hz, stored 0.29 Hz).

**10. Detune&Room** (Plate, 125 BPM). Pitch detune by gliding delays: the four plate voices are moved by the LFO
(0.20 Hz; triangle on V1 50 → 0 ms and V4 0 → 40 ms, sine / cosine on V2 / V3 between 10 and 40 ms) and glide toward
their targets, which shifts their pitch slightly up and down. ADJUST sets the glide response GldResp (0 … 100;
slower glide = smoother, smaller detune). Reverb 2.76 s, FX mix 0.59 / 0.81.

**11. Detune&Dbl** (Plate, 125 BPM). Detune as Detune&Room (program 10) on voices 1 and 4 only, plus voices 2 / 3 as a fixed 62-ms
doubling hard left / right; short 0.60 s plate. ADJUST = GldResp.

**12. Tight Space** (Plate, 120 BPM). A small, short plate (1.31 s, Size 7) with a doubling: voices 2 / 3 at 34 ms
hard left / right and voices 1 / 4 moved by the LFO (triangle, 0.33 Hz) in opposite directions (0 … 11 ms). ADJUST
("SplitWide") stretches all voice times (Delay Master 0 … 200 %) and pans voices 1 / 4 apart (50 → 0 / 100).

**13. Flange>Rvb** (Glide>Hall, 117 BPM). A through-zero style stereo flanger into the hall: the glide taps of each
side (A, B) are moved by the triangle LFO (0.10 Hz) in opposite directions (A: 0 ↔ 190, B: 0 ↔ 120 in 0.1 ms, i.e.
up to 19 / 12 ms), mixed at 0.5 each with glide feedback −0.25 / −0.30 (left) and +0.25 / +0.30 (right). The echo
voices pass the glide output at 0 ms (V1 left, V4 right). Hall 1.64 s. ADJUST = FX Mix.

**14. Flange+Rvb** (Chorus+Rvb, 120 BPM). Six chorus voices at 0–2 ms (Delay Master 1 %) with deep sweeps
(≈ 200 samples) and alternating ±0.27 feedback: a stereo flanger, voices 1–3 left, 4–6 right, into a 1.43 s reverb.
ADJUST = FX Mix.

**15. X Eko Flange** (Glide>Hall, 81 BPM). A cross-feedback flanger with tempo echoes. The glide taps sweep up to
≈ 20 ms (sine on A, cosine on B, 0.11 Hz, left and right in opposite directions; glide output −5 dB). The echo
voices 1 / 4 pass the line at 0 ms; voices 3 (3/2 beat) and 6 (5/2 beat) are not heard directly (level off) but feed
each line into the other (cross feedback ≈ −0.27 / +0.27 at the stored Mstr XFbk 30), so tempo echoes bounce between
the sides. ADJUST sets Mstr XFbk (0 … 100), i.e. how long the echoes last.

**16. Glide > Verb** (Glide>Hall, 75 BPM). Gliding taps (sine / cosine, 0.21 Hz; A up to 20 ms, B up to
14–15 ms) with glide feedback +0.19 / −0.12 give a slow chorus/flange that feeds the hall through diffused echo
voices at 0 ms (Voice Dif 100). The tempo voices 3 / 6 are off. Hall 0.79 s. ADJUST = FX Mix.

**17. Glide X-Ekos** (Glide>Hall, 130 BPM). As Glide > Verb with the echo voices 1 / 4 at 1000 ms × Delay Master
(stored 74 %: 140 ms after ADJUST) and cross feedback 0.6 between the lines: a modulated ping-pong echo. ADJUST
("Eko Time") sets the Delay Master 1 … 200 %, i.e. the echo time.

**18. Detuned Ekos** (Plate, 80 BPM). Two echoes, left voice 1 and right voice 3, both at a tempo fraction chosen
with ADJUST ("Ekos/Beat", 1 … 24 echoes per beat; stored 3 = 1/3 beat = 250 ms) and with negative feedback
(−0.12). The LFO (0.10 Hz) stretches all voice times by up to 7 % (Delay Master 100 → 107) and swaps the
post-delays (46 ↔ 0), detuning the repeats. Short plate (0.54 s), diffused voices (Voice Dif 75).

**19. StereoEqEkos** (M-Band+Rvb, 81 BPM). A six-tap tempo echo through the band filters: band delays at 1 … 6
units of a tempo fraction set by ADJUST ("Ekos/Beat", 1 … 24 per beat; stored 3 → 247 ms steps) with levels falling
5 dB per tap and the taps panned alternately left (1, 2, 3) and right (4, 5, 6) toward the centre. All bands share
the master filters (HiCut 83, LoCut 29); bands 3 and 6 feed back at −0.14. Reverb off (FX mix 1.00/0.00).

**20. ADJUpMyEchos** (Concert Hall, 61 BPM). A 3.1 s hall with four tempo echoes (1 beat left, 2 beats right,
1/2 and 1/4 beat inside) with slight positive feedback (+0.07 … +0.14). ADJUST ("EchosLvl", 0 … 10) brings the
echoes in: Levels Master 0 → 40; at the stored 0 the echoes are silent and only the hall (with tempo reflections at
1/4 and 1/2 beat) sounds.

**21. FSw2 Elevate** (Plate, 75 BPM). An echo that speeds up from a footswitch. The AR patches (tempo-coded,
16.6) replace the stored tempo values of voices 1 (left) and 3 (right): at rest both are 1 beat (800 ms; with −0.4
feedback after the master, diffused by Voice Dif 100, reverb out of the mix). Foot Sw 2 triggers the AR envelope
(attack 64 steps ≈ 1.3 s, release 61 steps), which divides both times step by step down to 1/12 beat (66 ms; the
denominator runs 1 → 12 with the envelope), raises the feedback master to 80, narrows the image (FX Width 405 → 315)
and lowers the input slightly; releasing the switch reverses it. ADJUST sets the feedback master (0 … 99, stored 45).
The LFO (0.10 Hz) stretches the voices by up to 7 % and swaps the post-delays as in Detuned Ekos (program 18). Sw 1 (ramp, triggered by
Foot Sw 1) is set up but drives nothing.

**22. Chorus&Pan** (Chorus+Rvb, 120 BPM). A chorus (voices at 0 / 3 / 5 ms each side, small depth, fast rates
≈ 1.9 Hz, slight negative feedback) auto-panned by the triangle LFO across the full width (Panning Master 0 ↔ 100).
ADJUST ("LFO Rate") sets both the chorus rate (MstRate 6 … 107 %) and the pan rate (0.08 … 3.74 Hz). Reverb 1.31 s,
not in the mix (FX mix 1.00/0.00).

**23. Chorus&Amb** (Chorus+Rvb, 120 BPM). A chorus with ambience: six voices at 12–22 ms with strong negative
feedback (−0.1 … −0.2 after the master), slowly auto-panned (triangle, one cycle per 3 beats, 25 ↔ 75) into a 2.2 s
reverb (FX mix 0.81/0.59). ADJUST ("MstChorus") sets the chorus rate (MstRate 6 … 150 %) and, above 50, the depth
(MstDepth 86 → 125).

**24. 6 Vox Chorus** (Chorus+Rvb, 120 BPM). Six chorus voices at 12, 18, … 42 ms spread from left to right, positive
feedback (≈ +0.1), voices only (reverb 1.72 s but FX mix 1.00/0.00). ADJUST ("Less/More", 0 … 50) raises depth
(MstDepth 15 → 85) and rate (MstRate 80 → 150 %) together.

**25. Split C&E** (Chorus+Rvb, 112 BPM). Chorus on one input, echoes on the other: line A (left input) has three
chorus voices (1 / 3 / 36 ms, increasing depth, slow rates 0.1–0.2 Hz), line B (right input) has three tempo echoes
(1, 1/2 and 1/4 beat, panned right, little modulation, feedback −0.19 / +0.16 / −0.16). ADJUST ("Input Pan") swaps
the inputs between the two lines (InPan L 0 → 100, InPan R 100 → 0). Reverb not in the mix.

**26. Env:PanKorus** (Chorus+Rvb, 120 BPM). A chorus (voices 1–33 ms, small depth, ≈ 1.5 Hz, alternating ±0.2
feedback) whose output is panned by the sine LFO (15 ↔ 85); the LFO rate follows an AR envelope triggered by the
input (Mono Lvl ≥ 36; attack 5 steps, release 56 steps): the pan speeds up from 0.21 Hz to several Hz while playing
loud and slows down over ≈ 1.1 s after. ADJUST ("Chorus") scales chorus rate (0 … 200 %) and depth (33 … 166).

**27. 6 Vox Flange** (Chorus+Rvb, 120 BPM). Six flanger voices at 0–2 ms with deep sweeps (≈ 200 samples) and
alternating ±0.3 feedback (masters 74 %), three left and three right in alternation, into a 1.3 s reverb. ADJUST sets
the depth (MstDepth 0 … 200, stored 70).

**28. Mod Max** (Chamber, 120 BPM). Everything moves with one LFO (0.25 Hz, rate set by ADJUST "Mod Knob" 0.10 …
5.0 Hz): the delay voices (V1 27 ↔ 42 ms, V2 17 ↔ 5 ms, V3 0 ↔ 7 ms; V3 / V4 nearly inverted at full level), the
input panning (triangle swaps InPan L / R), the reverb width (sawtooth 315 ↔ 405) and the output width (cosine
429 ↔ 291). The sawtooth also flips the LFO shape between sine and cosine once per cycle. Chamber 2.45 s, FX mix
0.71/0.71, voice feedback ±0.2.

**29. PreciseGlide** (Glide>Hall, 97 BPM). A slow, tempo-locked stereo glide: the LFO (triangle, one cycle per
8 beats) moves the four glide taps in opposite directions over 0 … 6 / 8 ms, mixed at about −3 dB. ADJUST
("Resonance") turns the glide feedback from none (stored) to Fbk L +0.24, Fbk R −0.24 with the cross feedback
removed, giving a resonant, pitched flange. The echo voices pass the glide output (0 ms, delay master 0). Sw 1 /
Sw 2 and Env L are set up but drive nothing. Hall 1.72 s.

**30. Round 147** (Chorus+Rvb, 120 BPM). A rotary-speaker imitation: six short voices (1–5 ms) with fast chorus
(≈ 2 Hz × MstRate) and mixed ±0.2 feedback, plus the LFO (cosine) panning the output 25 ↔ 75 and the sine narrowing
the width 405 ↔ 315 at the same rate (0.9 Hz stored). ADJUST ("RotorRate", 0 … 20) sets both the rotor rate (LFO
0.26 … 6.97 Hz) and the chorus rate (MstRate 25 … 200 %).

**31. FSw2 Rotary** (Chorus+Rvb, 120 BPM). A rotary speaker with slow/fast switching from Foot Sw 2: the switch
goes through a Latch (high 33), so each press toggles. The latched value triggers the AR envelope (attack 83 steps
≈ 1.7 s, release 56 steps ≈ 1.1 s), which ramps the rotor (LFO 0.40 → 6.0 Hz), the chorus rate (30 → 175 %) and
depth (60 → 125) and reduces the LFO depth, like a motor speeding up and down. The LFO pans the output 0 ↔ 100 and
sweeps the High Cut (91 ↔ 113). Voices at 2 / 5 / 11 ms (Delay Master 1 %), V1 / V4 with ∓0.5 feedback. ADJUST sets
FX Width (315 … 405). Without a footswitch the rotor stays slow.

**32. RotorCabinet** (Chorus+Rvb, 120 BPM). As FSw2 Rotary, switched by the ADJUST knob instead (0 = slow, 1 = fast:
the AR is triggered while ADJUST ≥ 1). A Latch on Foot Sw 1 is defined but nothing uses it.

**33. MIDI Rotary** (Chorus+Rvb, 111 BPM). As FSw2 Rotary, switched by MIDI aftertouch through the Latch (high 33),
with a shorter reverb (0.45 s) and tempo echoes in the reverb (EkoDly 1/3 and 1/4 beat). ADJUST lowers the echo
feedback (EkoFbk 15 → 0).

**34. Tiled Rotary** (Chorus+Rvb, 120 BPM). RotorCabinet (ADJUST 0 / 1 = slow / fast) with longer voice times
(28 / 70 / 154 ms, Delay Master 14 %) and a 0.94 s room in the mix (0.89/0.45).

**35. RotoWood** (Chorus+Rvb, 120 BPM). A small rotary sound without reverb: voices at 2–5 ms (Delay Master 10 %),
the LFO (2.22 Hz stored) pans 20 ↔ 80 and widens 380 ↔ 410. ADJUST ("Speed", 0 … 10) sets rotor rate (0.49 … 6.7 Hz),
chorus depth (45 → 119) and rate (83 → 174 %). The AR (triggered through a Latch on Foot Sw 1) drives nothing.

**36. RandomImages** (Chorus+Rvb, 120 BPM). Six voices with deep, slow chorus (MstDepth 200) whose pans are moved by
the LFO; the LFO rate follows the input level through a jagged curve (silence → 2.17 Hz, any signal → almost stopped,
loud → 0.03–0.05 Hz), so the image jumps when the input starts. Voices 3 and 6 are tempo echoes (2 and 4 beats ×
Delay Master) with feedback +0.76 / −0.9. ADJUST ("ImageKnob") walks through a set of combinations: delay master,
chorus depth, voice-3/6 feedback and output level change along multi-point curves.

**37. Under Water** (M-Band+Rvb, 120 BPM). Band echoes (200 / 70 ms each side, +0.1 … +0.4 feedback) whose master
filters are swept by the sine LFO at audio-like rates (11 Hz stored; ADJUST "DrownKnob" 10 … 25 Hz): High Cut
106 ↔ 161 and Low Cut 157 ↔ 206 move together, a bubbling, gargling sound. Long reverb (7.7 s), FX mix 0.89/0.45.

**38. Thunder FX** (M-Band+Rvb, 125 BPM). Six tempo-locked band echoes (from 1/7 to 24/23 beat, alternating
feedback −0.3 … +0.5) behind master filters driven by the AR envelope in Repeat mode: every time the input passes
Mono Lvl 31 the envelope jumps open (attack 0) and closes over AR Release (stored 440 steps ≈ 8.8 s), so the filters
open fully (High Cut 0 → 121, Low Cut 223 → 122) and then darken slowly, a flash followed by rumble. ADJUST
("FlashTime") sets the release (0 … 500 steps). Reverb 2.1 s, FX mix 0.99/0.13.

**39. Thunder FX 2** (M-Band+Rvb, 125 BPM). Thunder FX with all band times doubled (Delay Master 200 %), a longer
third band (8/5 beat), the output narrowed by the envelope (FX Width 360 → 405 on the flash), reverb out of the mix
and a short stored release (ADJUST 4 → 15 steps).

**40. ChaosImpuls** (Plate, 120 BPM). Four delay voices swept over very long ranges by the triangle LFO (0.11 Hz):
V1 1365 → 0 ms, V4 0 → 1365 ms, V2 0 → 715 ms, V3 523 → 0 ms, gliding (GldResp, GldRange 341) so their repeats
(feedback +0.16, V1 / V4 −0.5) sweep in pitch, diffused (Voice Dif 100). ADJUST sets the glide response. Plate 4.5 s.

**41. DemonDescent** (Chorus+Rvb, 120 BPM). Very slow, very deep chorus voices (MstDepth 200, MstRate 6 %) with
their delay times swept by the LFO (sawtooth 0.13 Hz: V1 0 → 1365 ms, V4 0 → 493 ms; triangle: V2 517 ↔ 274 ms,
V6 133 ↔ 843 ms) and gliding: continuously falling and rising pitch. Feedback ±0.3. ADJUST = GldResp. Reverb 3.1 s.

**42. RoughIdle FX** (Chamber, 50 BPM). A stutter / engine-idle effect on the reverb input: Sw 1 (ramp, 7.3 Hz
stored; ADJUST "Speed" 1.5 … 10 Hz), triggered while the input is above Mono Lvl 24, chops the input level
(96 ↔ 160), pans the input across and moves both post-delays and PstMix in each cycle. Only the chamber is heard
(FX Mix 100, 2.2 s). The AR (Repeat) is triggered too but drives nothing.

**43. RoxStutter** (Chorus+Rvb, 120 BPM). A chorus (voices 3–33 ms) plus two echo voices (770 / 884 ms, +0.31
feedback) whose delay follows the input envelope: Env L / Env R (release 60 steps = 1.2 s) pull the echo time down to
97 / 83 ms while the input is loud and let it glide back up as it decays, so echoes stutter and swoop. ADJUST = FX Mix.

**44. Split Pitch?** (Plate, 120 BPM). Two gliding delay voices, hard left (V1 1365 → 0 ms) and hard right
(V4 0 → 1365 ms), swept in opposite directions by the LFO (0.21 Hz) with feedback −1 × master 0 (none): the two sides
glide in pitch in opposite directions. ADJUST ("SplitWide") raises LFO depth (0 … 100) and rate (0.20 … 0.50 Hz) and
the Delay Master (0 … 200 %), so at 0 there is no split. No reverb (Mid Rt 0, FX Mix 0).

**45. TryTalk'n** (Plate, 120 BPM). ChaosImpuls with the triangle LFO at 0.21 Hz on V1 / V4 and the sawtooth on
V2 / V3, no voice feedback: gliding, talking pitch sweeps. ADJUST = GldResp.

**46. AutoInfinite** (Infinite, 120 BPM). An automatic freeze: the AR envelope, triggered while the input is above
Mono Lvl 64 (attack and release 0 = instantly), sets the tank's Infinite switch off while playing and on as soon as the
input drops, so the last sound is held indefinitely (16 s decay otherwise); new loud input releases it.
Delay voices 19 / 18 ms (LFO-modulated slightly, +0.32 feedback) and 1-beat voices 2 / 3 (−0.37 feedback) and LFO-swept
post-delays add movement. ADJUST = FX Mix (stored 70).

**47. Remove Cntr** (M-Band+Rvb, 120 BPM). Centre-channel removal: band 1 passes the left input and band 6 the right
input inverted (level 0 = −1), both at 0 ms and panned to the middle, so the output is L − R and anything panned to
the centre cancels. ADJUST ("Cntr Fc") raises the low cut of the inverted band (V6 LoCut 0 … 121) so that the centre
is removed only above that frequency (bass and kick stay). No reverb, no feedback.

**48. V-Eliminate** (M-Band+Rvb, 120 BPM). Vocal elimination: both inputs pass straight (bands 1 and 4) and an
inverted band-passed copy of each (bands 3 and 6, HiCut 117, LoCut 10, levels −1) is added on the other side,
cancelling the centre only inside the pass band. ADJUST ("Cut Zone") narrows the band (LoCut 10 → 49, HiCut 117 →
92).

**49. NoCenter Eko** (M-Band+Rvb, 143 BPM). Remove Cntr applied to a tempo echo: all four bands at 1/2 beat, the
straight ones with −0.08 feedback and the inverted ones with +0.08 (after the 40 % master), master high cut 92: the
echoes repeat only the side (L − R) content. ADJUST = V3 / V6 low cut as in Remove Cntr (program 47).

**50. Rvb On L-R** (Res1>Plate, 214 BPM). A reverb on the side signal only, using the resonator algorithm without
resonance: the six resonators are unassigned (feedback 0) and act only as their loop low-pass. Input voices 1 and 4
pass the left and right input at the centre, voices 3 and 6 add the left input inverted on the right and the right
input inverted on the left, so the plate (which in this algorithm is fed by the resonator section) receives
±(L − R)/2: centre-panned sources get no reverb. FX Mix 100 = reverb only. ADJUST ("Decay") sets Mid Rt 0 … 63
(stored 46 = 3.1 s). Reverb echoes at 1/6 and 1/4 beat.

#### Programs 51 – 100

**51. [n]Ekos/Beat** (Plate, 81 BPM). A stereo echo whose rate is chosen in echoes per beat: ADJUST ("[n]", 1 … 24)
sets voice 1 (left), voice 3 (right) and the pre-delay to 1/n beat together (stored 2: 1/2 beat = 371 ms), with
slight negative feedback (−0.12). The plate (1.31 s) is mostly out of the mix (FX mix 0.99/0.13). The LFO has depth 0.

**52. StereoTapDly** (M-Band+Rvb, 81 BPM). Two 1-beat echoes (742 ms), band 1 hard left and band 4 hard right,
with a dark master high cut (95) and feedback −1 × the master: ADJUST ("Mstr Fbk", 0 … 100, stored 20) sets the number
of repeats. A triangle LFO (8 beats per cycle) is patched to swap the pans but its depth is 0, so the echoes stay put.
No reverb.

**53. OffBeat Eko** (M-Band+Rvb, 81 BPM). Off-beat echoes: band 1 left at 3/4 beat × 50 % (278 ms) and band 4 right
at 3/2 beat × 50 % (556 ms), feedback −0.1, dark tape-like filtering (LoCut 17, master LoCut 138). ADJUST
("DarkKnob") sets the master high cut 0 … 240 (stored 85). Band 5 (1/4 beat) is set but off. No reverb in the mix.

**54. OffBeatEko2** (M-Band+Rvb, 81 BPM). OffBeat Eko with two extra quieter taps (−6 dB): band 2 at 1/2 beat
× 50 % (185 ms, left) and band 5 at 1/4 beat × 50 % (92 ms, right), all with tiny feedback; a denser off-beat
pattern. ADJUST = master high cut.

**55. Dotted 8ths** (M-Band+Rvb, 80 BPM). A six-tap rhythm in dotted eighths: taps at 3/4, 4/4, 6/4, 7/4, 9/4 and
10/4 beat (563 ms … 1.88 s) spread across the stereo field, each 4–5 dB quieter than the one before, from a mono
input (InPan 50 / 50); the feedback master is 0 (no repeats). ADJUST ("EchoSlope") tilts the level slope: at 0 the
first tap is loudest, at 127 the last (rising echoes). No reverb.

**56. Shuffle Ekos** (M-Band+Rvb, 81 BPM). A shuffle (triplet) echo: taps at 2, 3, 5, 6, 8 and 9 thirds of a beat ×
50 %, the on-beat taps (1, 3, 5) at −3 dB and the shuffle taps (2, 4, 6) at −22 dB at the stored knob. ADJUST
("Shufl Lvl") crossfades between the two groups (0: shuffle taps loud, main taps down; 127: the reverse).

**57. Shuffles** (M-Band+Rvb, 120 BPM). The shuffle tap pattern at full length (2/3 … 3 beats, 333 ms … 1.5 s) with
a falling slope (0 … −24 dB). ADJUST ("EchoSlope") sets how steeply taps 2–6 fall (0: steep, down to off; 100: all
equal).

**58. Sliding Eko** (Plate, 80 BPM). Two 1-beat echoes (voices 1 left, 3 right; feedback master 0, so single
repeats) with a 1-beat reflection pair in a very short plate. ADJUST ("Slide It!") moves the Delay Master 50 … 150 %
(stored 75 % = 562 ms); turning it makes the echoes glide in pitch while their delay slides. FX mix 0.71/0.71.

**59. Six Trips** (M-Band+Rvb, 94 BPM). Six taps at 1/4, 2/4 … 6/4 beat (159 … 957 ms), panned across, rising in
level from −10 to −3 dB at the stored knob; ADJUST ("EchoSlope") tilts the slope as in Dotted 8ths. No reverb.

**60. 6 StrokeRoll** (M-Band+Rvb, 86 BPM). A drum-roll build-up: five taps at 1/4 … 5/4 beat × 33 % (57 … 288 ms;
tap 6 is off) rising from −19 to −5 dB, a crescendo of strokes after each hit. ADJUST ("Attack") tilts the slope
(stored 85: rising). No reverb.

**61. TapGated Rvb** (Chamber, 81 BPM). A tempo-gated chamber: the LFO (pulse, 2 beats per cycle, short pulse
width 12) switches the reverb input (Rvb In) between off and the level set by ADJUST, so the chamber (2.3 s) is fed
only during a short window at the start of every 2 beats; reverb only (FX Mix 100), tempo echoes and reflections in
the tank (1/3 and 2/3 beat). ADJUST ("Gate dB") sets the open level (Rvb In 80 → 1, i.e. 0 to −79 dB).

**62. TapEkoGate 1** (M-Band+Rvb, 81 BPM). The OffBeatEko2 pattern (3/4, 1/2, 3/2, 1/4 beat; feedback +0.2 /
−0.14) behind a tempo gate on the input: the LFO (pulse, 4 beats per cycle, width 86) pulls the input level from
160 down toward 81 (off) for most of each 4-beat cycle. ADJUST ("GateWidth") sets the LFO depth (0 = no gating,
stored 47 = partial). Sw 1 / Sw 2 (pulse, 4 beats) are set up but drive nothing. Reverb out of the mix.

**63. TapEkoGate 2** (M-Band+Rvb, 81 BPM). As TapEkoGate 1 with other tap times (1/2, 3/2, 1/2, 3/4 beat) and the
gate patched from Sw 1 (pulse, 4 beats per cycle, width by ADJUST 0 … 98). Sw 1's trigger is "none ≥ 1", which by
6.3 never triggers, so the switch stays at its idle 0 (input open) unless it is restarted from outside (the Tap key
restarts the time switches). FX mix 0.99/0.16.

**64. TapEkoClear** (M-Band+Rvb, 81 BPM). The TapEkoGate 1 echoes with more feedback (master 71: +0.30 / −0.22),
cleared rhythmically: the LFO pulse (4 beats per cycle) drives the Clear switch, which holds the band input levels
(16.4) while the pulse is low, so the lines empty out, and lets the echoes build up again when it is high. ADJUST ("GateWidth")
sets the pulse width (stored 78 %). At rest (LFO output 0) the lines are cleared. FX mix 0.89/0.45.

**65. Tap Chamber1** (Chamber, 81 BPM). A short chamber (0.42 s) with two early voices (17 ms left, 29 ms right)
whose colour moves with the tempo: Sw 1 (ramp, 4 beats per cycle, always running) sweeps Rt HC (10 → 34), Crossover
and the Shape (10 → 149) once every 4 beats. ADJUST sets the reverb width (360 … 450).

**66. Tap Chamber2** (Chamber, 81 BPM). A long chamber (6.8 s, Size 127) with four early voices (5–15 ms) and
tempo echoes in the tank (3/4 and 12/15 beat). Sw 1 (ramp, 4 beats) breathes the reverb width (360 ↔ 460); the input
level triggers the AR (Mono Lvl ≥ 46, attack 3, release 50 steps), which skews the echo feedback (left down to 9,
right up to 21) while playing. ADJUST ("Decay") sets the Shape 0 … 254.

**67. LatchedEkos** (M-Band+Rvb, 81 BPM). Two echo pairs that take turns in tempo: 1/2-beat echoes (bands 1 / 4,
feedback −0.28) and 1/3-beat echoes (bands 3 / 6, +0.52). The LFO pulse (4 beats per cycle, width by ADJUST
"EchoWidth") gates the input (81 ↔ 160) and also feeds a Latch (low 11, high 110) that triggers the AR envelope
(attack / release 5 steps); the AR crossfades the levels between the pairs (V1 / V4 up, V3 / V6 down). No reverb.

**68. X-PanEQ BPM** (M-Band+Rvb, 97 BPM). Tempo auto-pan of the inputs through two band filters: bands 1 and 4 pass
the left and right lines at 0 ms with different band-pass settings (HiCut 80 / 87, LoCut 69 / 62); the sine and cosine
of the LFO (2/3 beat per cycle) pan each input across, so the two filtered colours circle the stereo field in time.
ADJUST ("Low Tone") sets the master low cut (61 … 121). Short reverb in the background (0.76 s, 0.98/0.22).

**69. Pan>Eko BPM** (M-Band+Rvb, 120 BPM). The inputs are panned by the LFO (sine, 2/3 beat per cycle) before
dark echoes: band 2 / 5 dry at 0 ms and band 1 / 4 at 200 ms (Delay Master 20 %, −6 dB, HiCut 112, feedback −0.53).
ADJUST ("DarkKnob") = master high cut 0 … 121. Reverb 1.31 s at 0.41.

**70. Tempo Verb** (Chamber, 97 BPM). A chamber whose decay follows the Tempo source (16.1): source 0 → Mid Rt 56
(7.7 s), 127 → 1; with no tempo value supplied the decay stays at 56. ADJUST ("Liveness") raises the reflections
(RefLvl 0 → 15), Spread (43 → 255) and Rt HC (14 → 45) together. Reverb only.

**71. Tempo Gate** (Inverse, 110 BPM). An inverse/gated reverb whose length follows the Tempo source (Duration
112 → 32 → 0 over source 0 … 30 … 100; stored 44). The four short voices (20–36 ms) are at −40 dB (Levels Master 0)
and FX Mix is 100, so only the reverb is heard. ADJUST sets the High Cut (68 … 118).

**72. Tape Echo** (M-Band+Rvb, 80 BPM). A dark tape echo: band 1 left and band 4 right at 1/2 beat (375 ms) with
+0.42 feedback through LoCut 17 and the master filters (HiCut 85, LoCut 138), so each repeat gets darker and thinner.
ADJUST ("DarkKnob") sets the master high cut 0 … 240. No reverb in the mix.

**73. NonLinear 1** (Chorus+Rvb, 80 BPM). A shaped echo burst from the chorus voices used as plain delays: six
voices at 1/2 … 6/2 beat × 50 % (188 … 563 ms), no modulation (MstDepth 0), panned around, with an ADJUST level curve
("EchoSlope"): at 64 all equal, toward 127 a rising slope (last voice loudest, like a reversed decay), toward 0 a
falling slope. No reverb.

**74. NonLinear 2** (Chorus+Rvb, 176 BPM). NonLinear 1 at 176 BPM (85 … 511 ms) with some chorus movement (MstDepth
50, slow MstRate 17 %), pre-delay 1/8 beat and a little reverb (0.9 s, FX mix 0.95/0.31).

**75. NonLinear 3** (Chorus+Rvb, 111 BPM). The shuffle pattern (2/3 … 3 beats × 25 %: 90 … 270 ms) on the chorus
voices with slight movement (MstDepth 50, MstRate 38 %), from a mono input; ADJUST level slope as NonLinear 1.

**76. Six Across** (M-Band+Rvb, 120 BPM). Six band echoes spaced 120 ms apart (120 … 720 ms) and panned from right
(longest) to left (shortest), each in its own frequency band (bass at 120 ms … treble at 720 ms); the feedback master
is 0 (single repeats). ADJUST ("Spacing") scales all times (Delay Master 0 … 200 %). No reverb.

**77. BandEko Rvb** (M-Band+Rvb, 143 BPM). A band-pass echo: bands 1 / 4 at 1 beat (419 ms) with −0.5 feedback,
through master filters that ADJUST ("Center Fc") moves together (HiCut 0 → 120, LoCut 122 → 240): a sweepable
band-pass centre. Reverb 1.2 s with tempo echoes at 3/2 beat, pre-delay 1/2 beat.

**78. BandEkoSweep** (M-Band+Rvb, 120 BPM). Three echo times (1/2, 1 and 3/2 beat) per side; the 1/2- and
3/2-beat bands have their filters swept by the time switches (Sw 1 ramp 0.15 Hz on bands 1 / 4, Sw 2 ramp 0.27 Hz
on bands 3 / 6), left and right in opposite directions, so the repeats sweep through the spectrum; the LFO breathes
the output width (315 ↔ 405). Feedback +0.8 / +0.4 / −0.4 × master (ADJUST "Mstr Fbk", stored 46). No reverb.

**79. LFO EQ Echo** (M-Band+Rvb, 80 BPM). 1/2-beat echoes (375 ms, −0.5 feedback) through master filters swept by
the LFO in tempo (4/3 beat per cycle): the low cut follows the sine (121 ↔ 205), the high cut the cosine (100 ↔ 146).
ADJUST ("LC Depth") = LFO depth. Reverb 1.9 s at 0.31.

**80. ChaseEcho** (M-Band+Rvb, 80 BPM). LFO EQ Echo with 1-beat echoes whose pans are patched from Sw 1 (ramp, rate
by ADJUST "ChaseRate" 0 … 5 Hz). Sw 1's trigger is "none ≥ 30", which never triggers (6.3), so at rest both echoes
sit at the left; the chase runs only when the switch is restarted from outside.

**81. Panned Dlys** (M-Band+Rvb, 125 BPM). Two echoes, 1 beat (480 ms, feedback −0.1) and 1/2 beat (240 ms,
feedback −0.5), auto-panned by the triangle LFO (0.25 Hz) in opposite directions. ADJUST = feedback master. No reverb.

**82. X-Pan Delays** (M-Band+Rvb, 176 BPM). Two 1-beat echoes (340 ms) crossing in the stereo field (triangle LFO,
0.25 Hz, opposite directions), feedback −0.25 at the stored ADJUST ("Mstr Fbk"). No reverb.

**83. Dly>EQ>Pans** (M-Band+Rvb, 333 BPM). Two filtered echoes, a dark low band at 2 beats (left, HiCut 65, LoCut
36, feedback −0.6) and a brighter band at 3 beats (right, HiCut 105, LoCut 55, feedback −0.36), i.e. 360 and 540 ms at
the stored 333 BPM (16.3), auto-panned in opposite directions (triangle 0.25 Hz). ADJUST = feedback master. No reverb.

**84. Ekoz 4 Drums** (Plate, 115 BPM). A drum ambience with moving echoes: four diffused voices (Voice Dif 90) at
25–92 ms (× 144 %) with feedback (+0.1, +0.05), slowly modulated (0.09 Hz) by sine, triangle and cosine together with
both post-delays; a 1.5 s plate with tempo echoes and reflections at 1/4 and 3/13 beat. ADJUST = FX Mix (stored 70).

**85. Haas PanKnob** (Plate, 125 BPM). Panning by time difference (Haas effect) of a mono input: voice 1 hard left
and voice 3 hard right at full level, voices only (FX Mix 0). ADJUST ("L=0,R=127") delays one side by up to 2 ms
(below 64 the right voice, above 64 the left one; 64 = centre), with the fastest glide (GldResp 100) so the knob
moves smoothly. The plate is not heard.

**86. Dial a Delay** (M-Band+Rvb, 80 BPM). A simple stereo echo whose time is "dialled": bands 1 / 4 at 1250 ms ×
Delay Master, ADJUST ("Delay") sets the master 0 … 200 % (stored 30 % = 375 ms); feedback −0.25. Reverb out of the
mix.

**87. PrecisionDly** (Glide>Hall, 117 BPM). A precise sub-millisecond delay between the channels: glide tap A of
each side at full level, the left one set by ADJUST ("L ms/100") from 0 to 1.0 ms in 0.1-ms steps (tap positions in
1/4096 sample, gliding). No feedback, no reverb. Several modulators are set (LFO, Env L, Sw 1, Sw 2) but nothing is
patched from them.

**88. PhoneOrRoom?** (M-Band+Rvb, 120 BPM). A two-way switch on ADJUST ("Pick One", 0 / 1). At 0 the six bands
pass the inputs at 0 ms through a narrow telephone band (master HiCut 74, master LoCut 200) with a narrowed image
(FX Width 360) and no reverb; at 1 FX Mix goes to 100 and only a small, short room is heard (Size 5, 1.1 s) at full
width.

**89. CheapTV Room** (M-Band+Rvb, 120 BPM). The telephone band (HiCut 80, master LoCut 202) collapsed to the centre
(all pans 50) plus a tiny boxy room (0.29 s, Size 2) at 0.41: a small speaker in a small room. ADJUST ("The Walls")
sets the low-frequency decay (Low Rt 0 … 9).

**90. Empty Stage** (Chamber, 120 BPM). A stage ambience: two early voices (49 ms left, 40 ms right) and a chamber
whose liveness is set by ADJUST: Mid Rt 0 … 39 (stored 0.86 s), Crossover, Shape (10 … 206) and High Cut
(97 … 114) together. FX mix 0.23/0.97.

**91. Tomb Room** (M-Band+Rvb, 120 BPM). A stone room: reverb only (FX Mix 100, the bands at −12 dB are not
heard), 1.89 s. ADJUST ("The Walls") makes the walls harder/softer: Low Rt 0 … 9, Rvb Out 15 → 6 and Shape 0 → 100.

**92. Comb Room** (Chorus+Rvb, 120 BPM). Static comb filtering: the chorus voices without modulation (MstDepth 0)
at 0–16 ms × Delay Master, with feedback ±0.3, colour the sound like a resonant tube; ADJUST ("Tone") sets the Delay
Master 0 … 200 % (stored 68 %), which tunes the comb. A 1.3 s reverb at 0.26.

**93. Zoom Over** (Chamber, 71 BPM). Something flying over, repeating: one sawtooth LFO cycle (0.07 Hz; ADJUST
"Speed" 0.01 … 0.25 Hz) raises the input level to the maximum at mid-cycle and back down (96 → 160 → 96), widens the
image (360 → 450), moves both post-delays through 0 at the centre (Doppler-like pitch shift), opens the High Cut at the
centre (59 → 120) and shortens the decay as it passes (Mid Rt 39 → 21 → 39). Reverb only. The AR (Repeat, trigger
Mono Lvl ≥ 24) drives nothing.

**94. OneShotCarBy** (Glide>Hall, 47 BPM). A car driving past, once per note: the input (above Mono Lvl 12) starts
the AR envelope in One Shot mode; during the attack (ADJUST "Speed": 500 → 25 steps = 10 … 0.5 s, stored 10 s) the
input is panned from left to right, its level rises to full at the middle and falls again (81 → 160 → 81), and glide
tap A moves from 42 ms to 0 at the middle and back to 42 ms: the approaching car is pitched up, the receding one down
(Doppler). Release 0: back to the start at once when the input stops. Hall 0.6 s. One Shot envelope: 16.2.

**95. AmbientCarBy** (Chamber, 120 BPM). The car-by on a chamber (reverb only): the AR in Repeat mode (attack set by
ADJUST, stored 10 s) pans the input across, swells its level and swings the post-delays (46 → 0 → 107 left,
98 → 0 → 46 right) for the Doppler shift; while the input stays above Mono Lvl 24, passes repeat. Repeat envelope: 16.2.

**96. HeadOn CarBy** (Chamber, 71 BPM). A continuous head-on pass driven by the sawtooth LFO (ADJUST "Speed"
0.10 … 2.0 Hz, stored 0.10 Hz): input level up to the middle of the cycle and down, width 360 → 450, post-delays
through 0 at the centre, PstMix, FX Adjust and Mid Rt (0 → 58 → 62) all follow the cycle. The input is not panned
(head-on). Reverb only.

**97. 2WayStreet** (Chamber, 50 BPM). Traffic both ways: Sw 1 in ramp mode (up then down; rate by ADJUST
0.10 … 2.0 Hz, stored 0.21 Hz), triggered while the input is above Mono Lvl 24, swells the input level, pans it
across and back and swings the post-delays; one ramp cycle = one pass out and one back. 2.2 s chamber, reverb only.

**98. IntoTunnel** (Chamber, 71 BPM). Driving into a tunnel and out, repeating with the sawtooth LFO (ADJUST
0.05 … 0.50 Hz): the input swells and pans across, the post-delays swing, the High Cut opens at the middle, and the
decay grows from 0 to 58 (5 s) toward the middle of the cycle. Reverb only.

**99. 2WayTunnel** (Chamber, 50 BPM). 2WayStreet in a tunnel: the triangle LFO (ADJUST 0.05 … 2.0 Hz, stored
0.08 Hz) drives the pass (level, pan, post-delays, PstMix) and its cosine opens the High Cut (56 → 121) and shortens
the very long decay (Mid Rt 62 → 33) in the middle of each pass.

**100. FinishLine** (Glide>Hall, 47 BPM). The OneShotCarBy pass in Repeat mode (cars keep passing while the input
is above Mono Lvl 25), plus long echo voices (V3 1.59 s and V2 2.88 s left, V6 1.74 s and V5 2.88 s right after the Delay Master) that
echo each pass like a grandstand. ADJUST ("Speed") shortens the echoes (Delay Master 200 → 37 %) and the pass (attack
500 → 25 steps) together. Hall 1.9 s at 0.26.

#### Programs 101 – 150

**101. StereoLoPass** (M-Band+Rvb, 120 BPM). A stereo low-pass filter: all six bands pass the input at 0 ms (−6 dB,
three per side; feedback master 0) and ADJUST ("CutOff") sweeps the master high cut 0 … 120 (stored 56). The reverb
(2.9 s) is not in the mix.

**102. Low Pass HiQ** (M-Band+Rvb, 94 BPM). A resonant low-pass: the band delays are 0 (Delay Master 0), so each band
feeds back around a one-sample loop with −0.82 feedback, which peaks the response at the cutoff. ADJUST ("CutOff")
sweeps the master high cut 0 … 120 and, in its top range, takes the resonance out (feedback master 100 → 90 at 90,
→ 0 at 127).

**103. Low Pass LFO** (M-Band+Rvb, 80 BPM). The resonant low-pass swept by the LFO in tempo (triangle, 8 beats per
cycle): master high cut 120 → 23 and feedback master 0 → 100 along the cycle. ADJUST ("Depth") = LFO depth (stored
89 %).

**104. StereoHiPass** (M-Band+Rvb, 176 BPM). A stereo high-pass: six bands at 0 ms, ADJUST ("CutOff") sweeps the
master low cut 122 … 240 (stored 190); feedback master 0.

**105. StereoB-Pass** (M-Band+Rvb, 120 BPM). A resonant stereo band-pass: ADJUST ("Center") moves the master high and
low cut together (HiCut 0 → 120, LoCut 122 → 240), with one-sample-loop feedback (−0.5, −0.2) for a peaked response.

**106. Stereo Notch** (M-Band+Rvb, 120 BPM). A notch filter: the six bands are phase-inverted (levels 0) band-passes
centred by ADJUST ("Center": HiCut 0 → 120, LoCut 122 → 240 together) and are added to the reverb path, which at
Size 0 / Mid Rt 0 passes the signal almost directly; the inverted band cancels the centre frequency (FX Mix 60 sets
the balance).

**107. SweptNotches** (M-Band+Rvb, 300 BPM). Four inverted resonant bands at 1 ms (feedback +0.3 … +0.55) whose
filters are swept continuously by the two time switches (Sw 1 ramp on bands 1 / 4, Sw 2 ramp on bands 3 / 6, left
and right in opposite directions), combined with the near-direct reverb path: moving notches like a phaser. ADJUST
("Rate") sets the Sw 1 rate (0 … 1.7 Hz); Sw 2 modulates its own rate.

**108. Env Notches** (M-Band+Rvb, 300 BPM). As SweptNotches with the filter sweeps driven by the input envelope
followers (Env L on bands 1 / 4, Env R on bands 3 / 6) from a mono input: louder playing moves the notches (auto-wah
style). ADJUST ("Release") sets both follower releases (0 … 50 steps).

**109. BandReject 4** (M-Band+Rvb, 300 BPM). Four 2-ms bands with ±0.4 feedback whose filters are swept by Sw 1 / Sw 2
(ramps, 0.16 / 0.24 Hz, opposite directions per side), the output width breathing with the LFO (315 ↔ 405), into a
4.2 s reverb. ADJUST = FX Mix.

**110. WaaPedalEko** (M-Band+Rvb, 120 BPM). A pedal wah: three band-passes at 0 / 4 / 0 ms (mono input) with
feedback near +0.9 resonate like a wah filter; the Footpedal moves their centres (band 1 HiCut/LoCut 46 → 66, band 2
84 → 71, band 3 89 → 86) and lowers the resonance in the middle of its travel. Only the bands are heard (FX mix
1.00/0.06); the reverb holds tempo echoes (2 and 3 beats). ADJUST ("FeedBack") sets the feedback master 90 … 100.
Without a pedal the stored settings apply (or those of a full pedal, 16.1).

**111. RotorWaa** (M-Band+Rvb, 120 BPM). The wah of WaaPedalEko (program 110) moved by the sine LFO in tempo instead of the pedal, with the
output panned 10 ↔ 90 by the same LFO; reverb 3.6 s at 0.79. The stored LFO Rate is 2 beats per cycle, but ADJUST
("Rate", tempo-coded patch, 16.6) replaces it: 8 beats per cycle at 0 (0.25 Hz), 5/3 beat at the stored 4 (1.20 Hz),
1/8 beat at 10 (16 Hz).

**112. Movable Echo** (M-Band+Rvb, 67 BPM). A left echo at 1/2 beat (449 ms, +0.21 feedback) and a right echo at
1/4 beat (224 ms, +0.5), dark filtering (master HiCut 109, LoCut 138). ADJUST ("Position") sets the output width
360 … 450, moving the echoes inward or outward. Tap delays are in 32-sample steps (11.1).

**113. Movable Hall** (Concert Hall, 120 BPM). A 2.2 s hall, nearly all reverb (FX mix 0.06/1.00), with two 80-ms
voices (left / right, small negative feedback, slightly diffused). ADJUST ("Position") sets the output width 360 … 450.

**114. Go Away** (Concert Hall, 120 BPM). Distance on one knob: ADJUST ("How Far?") darkens the sound (High Cut
111 → 28), raises the reverb share (FX Mix 10 → 100), the Def and Depth of the hall, the reflection levels and the
reverb width, all reaching their end at 96. At 0 the dry-like 0-ms voices dominate (FX mix 0.99/0.16). Hall 3.6 s.

**115. Circular Pan** (Chamber, 158 BPM). A sound circling the listener: the sine of the LFO pans the input (InPan L
and R together), the cosine moves the "front/back" cues at the same rate (High Cut 121 → 74, reverb share FX Mix 30 →
71, width 405 → 450, level −8 at the back). ADJUST ("Rate") sets the LFO rate (0.01 … 4.0 Hz, stored 0.08 Hz).
The AR, Sw 1 and Latch are set up but drive nothing.

**116. Spatial Hall** (Concert Hall, 136 BPM). A 2.1 s hall, reverb only, with tempo-locked reflections (1/6 and 2/9
beat) and post-delays (1/4 and 5/19 beat). ADJUST ("Position", 0 … 35) raises Depth (0 → 15) and widens the reverb
(390 → 425).

**117. Rear of Hall** (Concert Hall, 120 BPM). A large hall heard from the back: reverb only (Size 295), High Cut 109;
ADJUST ("Decay") sets Mid Rt 7 … 57 (stored 40 = 2.2 s).

**118. Backstage** (Concert Hall, 120 BPM). Go Away without the reverb crossfade (reverb only): ADJUST ("How Far?")
darkens the hall (High Cut 111 → 28) and raises Def, Depth and the reflections. Hall 2.45 s.

**119. Steered Rear** (Plate, 120 BPM). A long plate (10.5 s; ADJUST "Decay" Mid Rt 51 … 61) whose position follows
the playing: the AR (trigger Mono Lvl ≥ 44, attack 3, release 63 steps) brings the dry-like 0-ms voices forward while
playing (FX Mix 100 → 30, PstMix 100 → 50) and narrows the reverb (450 → 375), and lets the plate come back as the
input stops. The LFO (0.34 Hz) moves the post-delays and the voice times (0 … 13 / 15 ms) for slow chorusing.

**120. Too Deep!** (Concert Hall, 120 BPM). A deep hall (Size 296; ADJUST "Decay" Mid Rt 23 … 60, stored 3.6 s) with
two 0-ms voices and two echoes (409 / 433 ms, feedback ±0.2), diffused (Voice Dif 67). The input envelopes push the
post-delays apart (Env R → PstDly L, Env L → PstDly R, 0 → 47) while playing, and the AR narrows the reverb.

**121. Dyna-Hall** (Concert Hall, 120 BPM). A hall (Size 174; ADJUST "Decay" Mid Rt 28 … 58, stored 3.9 s) that gets
wider when played loud: the AR (Mono Lvl ≥ 16, attack 10, release 64 steps) moves the reverb width 372 → 416.

**122. RotoRox** (Chorus+Rvb, 120 BPM). RoxStutter's envelope-pulled echoes (770 / 884 ms → 97 / 82 ms while loud,
release 2 s) with a chorus (7–43 ms), a 12.7 s reverb and a moving width (triangle LFO 405 ↔ 315) whose rate jumps
from 0.30 to 3.65 Hz while the input is loud (AR, trigger Mono Lvl ≥ 37). ADJUST = FX Mix.

**123. RotoRoomVox** (Chamber, 79 BPM). A chamber (reverb only; ADJUST "Decay" Mid Rt 5 … 56) whose output width is
swept by a sawtooth LFO over the whole range (0 … 720, i.e. mono through to wide and swapped), slowly (0.07 Hz) while
the input is loud and fast (0.60 Hz) when it is quiet (AR, trigger Mono Lvl ≥ 43, attack 4, release 7 steps).

**124. RotoRoom** (Chamber, 120 BPM). A chamber (ADJUST "Decay" Mid Rt 0 … 62) with the output width swept by the
triangle LFO (360 ↔ 425, 0.35 Hz). The AR is set up but drives nothing.

**125. RotoRoom#2** (Chamber, 120 BPM). A chamber with a rotating image: sawtooth LFO on the output width (0 … 720),
sine / cosine on the two voice times (26 ↔ 41 ms; voice 1 almost fully inverted, left; voice 3 right) and the triangle
moving the post-delays in opposite directions. ADJUST ("Rt & Rate") sets the decay (Mid Rt 6 … 61) and, inversely,
the LFO rate (1.17 … 0.05 Hz).

**126. Tremolo Tap1** (M-Band+Rvb, 81 BPM). A tempo tremolo: the LFO (sine, 1/2 beat per cycle) moves both input
levels through 160 → 120 → 160, a level dip at mid-cycle, so the tremolo runs at twice the LFO rate (four dips per
beat); bands 1 / 4 pass the inputs at 0 ms (FX mix 0.99/0.16). ADJUST ("Depth") is also patched to the input levels
(40 → 0).

**127. X-TremoloTap** (M-Band+Rvb, 143 BPM). A stereo (alternating) tremolo: the LFO (1 beat per cycle) moves the
left input 120 → 160 and the right 160 → 120, so the level swings between the sides once per beat. ADJUST ("Depth")
is patched to the input levels as in Tremolo Tap1 (program 126).

**128. Panner BPM** (M-Band+Rvb, 61 BPM). A tempo auto-pan of the inputs (2 beats per cycle): the left input is panned
by the sine, the right by the LFO output, whose shape ADJUST ("Pan Phase") selects: 1 = cosine (the two inputs circle
90° apart, stored), 0 = sine (both move together).

**129. Nice Pan!** (Plate, 106 BPM). A diffused plate (Voice Dif 99, 3.9 s) whose four short voices (18–38 ms, ±0.2 …
0.4 feedback) are panned as a whole by the triangle LFO (Panning Master 0 ↔ 100, 0.24 Hz) while the post-delays
move with it. ADJUST = FX Mix.

**130. Spin & Duck** (Chorus+Rvb, 81 BPM). The rotary voices of FSw2 Rotary with the rotor speed on ADJUST ("Spin
Rate": LFO 0.40 … 6.0 Hz, chorus depth 156 → 48, rate 30 → 175 %), pan and High Cut moved by the sine, and ducked
echoes: the input level sets the reverb echo feedback (EkoFbk 11 when quiet → 15 = none when loud), so 1/2-beat
echoes appear when the playing stops. The AR (through a Latch on Foot Sw 1) drives nothing.

**131. MultiFx Fade** (Chorus+Rvb, 80 BPM). Chorus voices and short echoes (30–129 ms after the Delay Master 10 %,
mixed depths, ±0.4 feedback) with an automatic swell: the AR (Mono Lvl ≥ 16, attack 7, release 34 steps) raises the
input level (136 → 160), opens the High Cut (47 → 111) and swaps the pans of voices 3 and 6 while playing. ADJUST =
FX Mix.

**132. Ghost** (Inverse, 120 BPM). An inverse reverb (Duration 64), reverb only, whose two post-delays are moved in
opposite directions by the LFO (0.14 Hz, 32 ↔ 0): a slow, ghostly flanging of the reverse swell. ADJUST ("FadeShape")
sets Diffusion (0 … 100) and Shape (0 … 255) together.

**133. Ghost Flange** (Inverse, 120 BPM). Ghost with a longer Duration (89).

**134. GhostVibrato** (Inverse, 120 BPM). The inverse reverb (Duration 89) with a vibrato: the LFO (triangle,
3.72 Hz) moves both post-delays (0 → 32), and its depth is brought in by the AR while playing (Mono Lvl ≥ 12; attack
203 steps ≈ 4 s, release 15: LFO Depth 0 → 45). ADJUST ("FadeShape") as in Ghost.

**135. AutoFadeIn 1** (Chorus+Rvb, 80 BPM). An automatic fade-in: the AR (Mono Lvl ≥ 31, attack 76 steps ≈ 1.5 s at
first, release 39) raises the input level 124 → 160, widens the image (315 → 405) and opens the High Cut (47 → 111);
the AR also lengthens its own attack as it rises (14 → 54 steps). Chorus voices (0 ms, deep) and tempo echoes (1/2 beat
left, 3/2 beat right, ∓0.6 / +0.4 feedback). ADJUST = FX Mix.

**136. AutoFadeIn 2** (Chorus+Rvb, 80 BPM). The fade-in with its attack on ADJUST ("Fade Rate": AR attack 0 … 150
steps, stored 45 ≈ 0.9 s): input level 136 → 160, High Cut 47 → 111, voices 3 / 6 (1/2-beat echoes) swapping sides.
Chorus voices at 0–17 ms. Reverb 4.5 s at 0.71.

**137. AutoFadeIn 3** (Chorus+Rvb, 70 BPM). The fade-in (attack 7 → 19 steps as it rises, release 4) on tempo echoes
(2, 3, 4, 5 beats × 11 % and two 1/2-beat × 11 % voices, 47–142 ms) with a short reverb (0.5 s). ADJUST ("Rvb Lvl")
sets Rvb Out 0 … 15.

**138. ChordSwells** (Chorus+Rvb, 80 BPM). Swelling chords: slow chorus voices at 200 ms … 1.2 s with feedback
(+0.25 / −0.2 / −0.4) build up sustained layers; the AR (attack 14 → 100 steps as it rises, release 74 steps) fades the
input in (136 → 160), opens the High Cut and swaps voices 3 / 6. Reverb 4.5 s. ADJUST ("Rvb Lvl") = Rvb Out 0 … 15.

**139. BowedChords** (Chamber, 44 BPM; *Part I*). Bowed attacks: the AR (trigger Mono Lvl ≥ 5, attack 63 steps =
1.26 s, release 9) fades the input in from −40 dB to 0 dB at every new note, removing the attack. Delay voices: 0 ms
left, 5 ms right and two 1/4-beat echoes (341 ms) with −0.46 feedback into a 2.3 s chamber. ADJUST = FX Mix (stored
92).

**140. BowedEchoes** (Chamber, 81 BPM). Bowed echoes: voices 2 / 4 at 1/2 beat (371 ms, left / right) whose feedback
falls while playing and rises after (AR, Mono Lvl ≥ 5, attack 27, release 9 steps: feedback master 85 → 22, i.e.
about −0.8 → −0.2), while the input level swells 130 → 160 with the same envelope; voices 1 / 3 at 0 ms. ADJUST =
FX Mix (stored 65).

**141. Pedal Swell** (Chorus+Rvb, 120 BPM). A volume-pedal swell into chorus and echoes: the Footpedal sets the input
level (80 = off → 160). Voices 1, 2, 5, 6 are 400-ms echoes (+0.2 feedback, chorused) and voices 3 / 4 a 0-ms chorus,
all with 145-sample sweeps at 0.4 Hz; reverb 3.9 s at 0.16. ADJUST = FX Mix. Without a pedal the input
passes at full level (16.1).

**142. RezClimber** (Res1>Plate, 81 BPM). A climbing arpeggio of resonances: each note above Mono Lvl 7 starts the AR
(One Shot, attack 50 steps ≈ 1 s), which steps Assign through 0, 7, 12, 17 (plus the ADJUST offset "Tuning" 0 … 60,
stored 19), and each new Assign value tunes the next resonator in turn (Res1 round robin, 14.4): root, fifth, octave,
eleventh build up and ring. The input level sets the resonance (Mono Lvl → Mstr Res 0 … 100). Input voices at tempo
fractions × 7 %, unison on. At rest no resonator is assigned (13.4).

**143. Chord Walk** (Res2>Plate, 90 BPM). A walking chord: the Res2 harmoniser builds a six-note chord on Assign;
Sw 1 (ramp, 8 beats) steps Assign through 36, 38, 40, 41 (C, D, E, F), Sw 2 (pulse, 8 beats) shifts the Root by 4,
and the LFO bends voice 6's interval (44 → 48), so the chord changes in tempo; strong negative loop feedback (−0.975)
makes each resonator ring an octave lower with odd harmonics (13.4). The input voices 1 and 4 are set by tempo-coded
patches from Sw 1 / Sw 2 (16.6): 4 beats × Delay Master 50 % = 1 333 ms, 8/3 beat (888 ms) while the switch is at its
top. ADJUST ("High Cut") = resonator high cut.

**144. Mars Bars** (Res1>Plate, 81 BPM). Evolving clusters: the LFO (cosine, 6 beats per cycle) steps Assign through
36, 44, 59–66, and each new value is handed to the next resonator (round robin), so the six resonators collect a
changing set of notes; resonators alternate positive (bright) and negative (hollow, octave-down) resonance. Input
voices at 1, 1/2, 3/2 and 2/3 beat with feedback. ADJUST ("High Cut") = resonator high cut. Reverb 3.3 s at 0.89.

**145. ModalImpulse** (Res2>Plate, 80 BPM). An arpeggiated chord of near-infinite resonators (loop feedback ±0.999):
E3, G3, B3, C3, C4, G3, excited one after the other by input voices at 0, 2/3, 1, 4/3, 5/3 and 2 beats, so an impulse
rings out as a slowly building chord. ADJUST ("Tone") opens the resonator high cut (71 → 121) and lowers the
resonance (100 → 57).

**146. MajorMinor** (Res2>Plate, 81 BPM). A chord that grows and changes mode: the sawtooth LFO (24 beats per cycle)
raises Active from 0 to 6 resonators, the pulse of the same LFO switches the scale between major and minor each
half-cycle, and input voices at 3/2 … 7/2 beats feed them. Sw 1 (patched to Root) has trigger "none ≥ 1" and does not
run on its own. ADJUST ("Tone") = resonator high cut and resonance. At rest Active is 0: no resonator rings, but the
input voices still reach the output through the loop low-passes (14.4).

**147. MIDIchords** (Res2>Plate, 120 BPM). Resonant chords played from MIDI: the highest held note (High Note) sets
Assign, the harmoniser tunes six resonators to the chord (full resonance), note velocity sets their level, the sine LFO
widens and narrows the image (135 ↔ 585), and the AR on the input moves FX Adjust. Until MIDI notes arrive the patches
from High Note and Velocity have no value (16.1): the stored Assign 48 sets the chord (on C4) and ADJUST alone sets the
level. ADJUST ("Sparkle", 0 … 5) opens the high cuts and raises the resonator level.

**148. LvlSweeper** (Res1>Plate, 120 BPM). Level-driven pitch: loud input (Mono Lvl ≥ 35) starts the AR, whose
release (62 steps) walks Assign down through a stepped curve (60 … 24), assigning a falling series of notes to the
resonators in turn as each note decays; the AR also sets the resonator level and resonance. Sw 1&2 (2 Hz pulses,
triggered by the input) pan the output. Input voices at 1/3 … 1/8 beat. ADJUST ("Color") offsets Assign 0 … 24.

**149. Sweeper** (Res1>Plate, 120 BPM). A rising sweep of resonances: the sawtooth LFO (1 beat per cycle) sweeps
Assign 36 → 60, so each beat the round robin hands a rising run of notes to the six resonators; Sw 1 (ramp, 4 beats)
pans the output. Input voices at 1/2 … 3/5 beat, reverb with 2-beat echoes. ADJUST ("Pitch") offsets Assign 0 … 36.
The AR (Foot Sw 1) drives nothing.

**150. MIDIsustain** (Res1>Plate, 120 BPM). A sustain-pedal resonator: the last MIDI note sets Assign (round robin);
Foot Sw 2 (sustain) triggers the AR (attack 2, release 33), which raises Active 0 → 6, stops clearing the input lines and
switches unison off, and the switch itself sets the resonance (98 → 100). Without a footswitch Active stays 0 and no
resonator sounds. The LFO (2.1 Hz) moves the reverb echo times. ADJUST ("Bright") = Rt HC, Low Rt and
resonator level.

#### Programs 151 – 200

Programs 151 – 183 are Part I's programs (*Part I*); their values are in the program data (17.2) like all
others, and the group notes and the tank programs are in Part I. Unless stated, they use no modulators, and ADJUST is a decay or liveness control.

**151. Small+Stage** (Concert Hall, 120 BPM; *Part I*). A small hall (2.3 s, Size 136) with a stage: four early voices
(9–37 ms × 70 %, diffused) around the image. ADJUST ("Stage Lvl", 0 … 50) raises the voices (Levels Master 20 → 50)
and the voice share (FX Mix 100 → 66). The input levels pan voices 2 and 4 toward the centre as the left / right
signal gets louder.

**152. Small Room** (Chamber, 120 BPM; *Part I*). A small room (Size 27, 1.56 s at the knob), reverb only, with early
reflections at 2 / 3 ms; ADJUST ("Liveness") sets Mid Rt and Shape together along a two-segment curve.

**153. Living Room** (Chamber, 120 BPM; *Part I*). A short, soft room (Size 42, 0.63 s), reverb only; ADJUST ("Decay")
= Mid Rt 0 … 60.

**154. Brick Kick** (Chamber, 120 BPM; *Part I*). A tight, hard room for kick drum (0.40 s, Link off) with two early voices
(17 ms left, 29 ms right) and reflections; ADJUST ("Liveness") brightens and reshapes the decay (Rt HC 10 → 34,
Crossover, Shape 10 → 186).

**155. Large Room** (Chamber, 120 BPM; *Part I*). A large room (Size 64, 1.89 s at the knob), reverb only; ADJUST
("Decay") = Mid Rt 0 … 59.

**156. SnareChamber** (Chamber, 120 BPM; *Part I*). A short bright chamber for snare (0.42 s) with the Brick Kick early
voices and reflections; ADJUST ("Liveness") as Brick Kick over 0 … 60.

**157. Tiled Room** (Chamber, 120 BPM; *Part I*). A small tiled room (Size 23, 2.32 s) with two very early voices
(2 ms left, 4 ms right); ADJUST ("Decay") = Mid Rt 0 … 62.

**158. Rich Chamber** (Chamber, 120 BPM; *Part I*). A dense, dark chamber (Diffusion 86, 2.45 s), reverb only; ADJUST
("Decay") = Mid Rt 0 … 61.

**159. Vox Chamber** (Chamber, 100 BPM; *Part I*). A long vocal chamber (6.8 s, Size 127) with four early voices
(5–15 ms) and tank echoes at 354 / 350 ms whose feedback rises while singing (AR on Mono Lvl ≥ 46, attack 3, release
50 steps: EkoFbk 15 → 20). ADJUST ("Liveness") sets Rt HC and Shape (0 … 254).

**160. Locker Room** (Chamber, 120 BPM; *Part I*). A short hard room (Size 21, 2.45 s) with two early voices, reflections
and a 14-ms post-delay mixed in (PstMix 100); ADJUST ("Depth") = Rt HC, Crossover and Shape.

**161. Wide Chamber** (Chamber, 94 BPM; *Part I*). A wide chamber (Rvb Width 413, AP Sweep 19, 1.64 s), reverb only;
ADJUST ("Decay") = Mid Rt 0 … 60.

**162. Gate Chamber** (Inverse, 95 BPM; *Part I*). A gated reverb (Inverse algorithm, flat-ish Shape 16): ADJUST
("Duration", 0 … 64) sets the gate length (Duration 0 … 112, stored 56).

**163. Vox Plate** (Plate, 97 BPM; *Part I*). A vocal plate (3.34 s at the knob, Attack 23) with two 3-ms voices and
tempo-locked tank echoes (1 beat and 7/8 beat); ADJUST ("Decay", 0 … 30) = Mid Rt 30 … 60.

**164. Good olPlate** (Plate, 100 BPM; *Part I*). A classic plate (3.86 s), reverb only, post-delay mix 100 with a
small post-delay offset on the right; ADJUST = Mid Rt 0 … 60.

**165. Slap Plate** (Plate, 106 BPM; *Part I*). A plate (2.93 s) with two slap-back voices (72 ms left, 40 ms right);
ADJUST = Mid Rt 26 … 59.

**166. Brass Plate** (Plate, 100 BPM; *Part I*). A bright, dense plate (Diffusion 93, 2.93 s) with two 5-ms voices;
ADJUST = Mid Rt 30 … 61.

**167. Drum Plate** (Plate, 100 BPM; *Part I*). A short drum plate (2.2 s, Attack 25) with voices at 23 / 7 ms and
reflections; ADJUST = Mid Rt 34 … 59.

**168. Rich Plate** (Plate, 100 BPM; *Part I*). A very dense plate (Diffusion 98, Attack 58, 3.86 s), reverb only; ADJUST
= Mid Rt 0 … 62.

**169. Concert Wave** (Concert Hall, 86 BPM; *Part I*). A large hall (Size 304) that "waves": the tempo-locked LFO (one
cycle per 2 beats, 0.714 Hz; 16.3) moves both post-delays (0 → 21 / 23 → 0, PstMix 88) and, through its square output,
the reverb width (360 ↔ 405); its triangle also selects the LFO shape. ADJUST ("Wave Knob", 0 … 40) sets Mid Rt
(18 → 53), Diffusion (100 → 2) and Def (1 → 100).

**170. Concert Hall** (Concert Hall, 120 BPM; *Part I*). A large concert hall (Size 271, 3.59 s) with early voices at 0 /
135 / 139 ms and post-delays of 139 ms; ADJUST ("Decay", 0 … 100) = Mid Rt 0 … 60. It decays faster than
displayed (4.2).

**171. ConcertHall2** (Concert Hall, 120 BPM; *Part I*). A concert hall with reflections and 147-ms post-delays; ADJUST
("Empty/Full", 0 … 15) shortens and darkens it as the hall "fills" (Mid Rt 50 → 35, Rt HC 29 → 16).

**172. Piano Hall** (Concert Hall, 120 BPM; *Part I*). A hall for piano (3.59 s) with four early voices (0 / 31 / 0 /
26 ms, slight feedback, diffused 67). The AR (Mono Lvl ≥ 31, attack 36, release 67 steps) opens the post-delays
(0 → 26 / 27) while playing, and Sw 1 (ramp, 3.85-s cycle, triggered at Mono Lvl ≥ 75) narrows the reverb on loud
passages. ADJUST = Mid Rt 0 … 61.

**173. Medium Hall** (Concert Hall, 120 BPM; *Part I*). A medium hall (Size 192, 2.32 s), reverb only; ADJUST = Mid Rt
0 … 60.

**174. Vocal Hall** (Concert Hall, 120 BPM; *Part I*). A vocal hall (4.54 s at the knob) with tempo-locked reflections
(1/6 beat) and post-delays (1/2 and 3/7 beat); ADJUST = Mid Rt 0 … 58.

**175. Deep Space** (Chorus+Rvb, 120 BPM; *Part I*). A huge, dark ambience (12.7 s, pre-delay 137) behind a chorus
(voices 7–43 ms) and two half-second echoes (+0.3 feedback); the AR (Mono Lvl ≥ 37, attack 10, release 83 steps)
brightens the reverb while playing (Rt HC 11 → 27) and speeds the width LFO (0.02 → 0.47 Hz, triangle on FX Width
388 ↔ 410). ADJUST = FX Mix (stored 69).

**176. Plate4 Horns** (Plate, 120 BPM; *Part I*). A horn plate with two early voices (28 / 22 ms) and reflections;
ADJUST ("Impact", 0 … 100) raises the decay (Mid Rt 0 → 59) and the attack (0 → 100) together.

**177. Alley Slap** (Plate, 120 BPM; *Part I*). A slap-back alley: a very short plate (Mid Rt 0) with tempo-locked tank
echoes (1/8 and 1/9 beat) and feedback; ADJUST sets the pre-delay (6 … 100).

**178. Drum Gate** (Inverse, 120 BPM; *Part I*). A gated drum reverb (Shape 255) whose width follows the left input
envelope (Env L → Rvb Width 416 → 361 while loud); ADJUST ("Duration", 0 … 70) sets both the gate length (Duration
22 … 112) and the envelope release (17 … 35).

**179. Slope Down** (Inverse, 120 BPM; *Part I*). An inverse reverb (Duration 94) whose post-delays are swept by the AR
in Repeat mode while the input is loud (Mono Lvl ≥ 64, attack 75 steps): PstDly 10 → 682 with PstMix 0 → 70, a
long sliding echo. ADJUST = Duration 0 … 112.

**180. BigBoomRoom** (Chamber, 120 BPM; *Part I*). A boomy room (1.31 s) whose post-delays open with the input envelopes
(Env L → PstDly R, Env R → PstDly L); the input level is also patched to the LFO rate, but nothing is
patched from the LFO; ADJUST ("BoomKnob", 0 … 36)
raises the low-frequency decay (Low Rt 4 → 9) and crossover (8 → 20).

**181. Whammy Hall** (Concert Hall, 120 BPM; *Part I*). A hall (4.54 s at the knob) with a fast vibrato on the
post-delays (LFO cosine, 4.98 Hz) that loud playing switches off: the AR (Mono Lvl ≥ 73, release 114 steps) takes the
LFO depth to 0 and lets it return over 2.3 s. ADJUST = Mid Rt 30 … 60.

**182. JetChamber** (Chamber, 120 BPM; *Part I*). A jet-flanged chamber (7.7 s): the sine LFO (0.08 Hz)
moves the two post-delays in opposite directions (23 ↔ 0); ADJUST ("Jet Knob") sets PstMix 0 … 50, i.e. the
flange depth. AP Sweep 31, tank echoes at 23 / 23.

**183. EnvChamber** (Chamber, 120 BPM; *Part I*). JetChamber with the post-delays moved by the AR (Mono Lvl ≥ 28,
attack 20, release 13 steps) instead of the LFO: the flange sweeps with each phrase.

**184. Wizz & Wazz** (Plate, 81 BPM). Everything on one tempo LFO (4 beats per cycle): voice 1 (12 ↔ 23 ms, sine) and
voice 3 (34 ↔ 12 ms, triangle) sweep and flange, the post-delays rise and fall, and the sawtooth pans the whole
output from left to right and crossfades voices → reverb (FX Mix 0 → 100) each cycle; voices 2 / 4 are long tempo
echoes (17/18 and 9/10 beat) with +0.3 feedback, all diffused (Voice Dif 100). ADJUST ("Decay") = Mid Rt 0 … 60.

**185. Sci Fi** (Chamber, 120 BPM). A tiny, short chamber (Size 0, 0.6 s, reverb only) with a vibrato on the
post-delays (triangle LFO, opposite directions, 0 ↔ 57) whose rate falls while playing (AR: 3.2 → 1.26 Hz). ADJUST
("Wiggle") = LFO depth.

**186. Wobble Plate** (Plate, 120 BPM). A 4.5 s plate (reverb only) with its post-delays moved by triangle and cosine
(0.26 Hz); ADJUST ("Wobble") sets the delay glide response (GldResp 30 … 60), i.e. how strongly the moving delays
bend the pitch.

**187. Dyna Vibrato** (Glide>Hall, 117 BPM). A dynamic vibrato on gliding taps: the triangle LFO moves the four
glide taps (up to 12–15 ms), but its rate and depth follow the input level (quiet: 5.3 Hz at depth 38; loud: 2.0 Hz at
depth 0), so the vibrato appears as notes decay. ADJUST ("Glide") sets the B taps' base delay (0 … 42 ms). Glide
feedback +0.19 / −0.12, diffused echo voices at 0 ms, short hall (0.76 s).

**188. VibroVerb** (Chamber, 120 BPM). A 1.4 s chamber (reverb only) with a vibrato on the post-delays (triangle,
0 ↔ 57); the AR in Repeat mode (while Mono Lvl ≥ 24) cycles the LFO rate between 4.1 and 1.6 Hz. ADJUST ("V-Depth")
= LFO depth.

**189. SweepVerb** (Chamber, 120 BPM). A 3.1 s chamber (reverb only) whose post-delays sweep apart on each attack
(AR, Mono Lvl ≥ 32, attack 13, release 20 steps: PstDly L 0 → 114, R 114 → 0); ADJUST sets the glide response
(GldResp 0 … 100), i.e. how fast and how pitched the sweep is.

**190. EnveloVerb** (Chamber, 120 BPM). A ducking, blooming reverb: at rest a 21 s chamber; while the input is loud
(AR, Mono Lvl ≥ 64, attack 1, release 75 steps) the decay drops (Mid Rt 61 → 51) and the reverb is pulled down
(FX Mix 100 → 38), and after the input stops it returns over 1.5 s. ADJUST ("Rt HC") = 5 … 45.

**191. Super Ball!** (Glide>Hall, 111 BPM). A bouncing ball: turning ADJUST to 1 triggers the AR in One Shot mode
(attack 390 steps ≈ 7.8 s): the input to the echoes is cut (Gld Lvl → 0), the echo voices (1/2 beat left, right
voice) get strong negative feedback (−1 → −0.9), and the Delay Master falls from 200 % to 4 %, so the captured sound
bounces faster and faster; at the top the lines are cleared. Back to 0 resets it. A Latch on Foot Sw 1 is defined but
not used.

**192. LevlStutter** (M-Band+Rvb, 53 BPM). A six-tap stutter: the AR patches (tempo-coded, 16.6) replace the tempo
values, so at rest tap j sits at j beats (1 132, 2 265, … 6 796 ms at 53 BPM; with the standard memory the taps beyond
2 730 ms are halved, 14.1), panned left to right; a loud input (Mono Lvl ≥ 40) starts the AR (One Shot, attack 5,
release 50 steps), which divides every tap time by a denominator rising from 1 to 24 (down to j/24 beat, 47 ms steps)
and back. The cosine LFO (16 beats) pans the output slowly. ADJUST = feedback master (stored 0). No reverb.

**193. Freez 2 of 8** (Glide>Hall, 81 BPM). A freeze looper with 1/2-beat loops (371 ms, left / right). Feedback on
the loops is set by ADJUST (V1 / V4 Fbk 100 → 200, i.e. none → maximum, 0.99994 per pass: practically lossless, 14.3).
Foot Sw 2 clears the loops. Sw 1 (pulse, 8 beats, width 12) is patched to freeze for part of every 8 beats (input cut
via Gld Lvl, feedback to maximum), but its trigger is "none ≥ 1", which does not run on its own (6.3); a Latch on
Foot Sw 1 is defined but not used.

**194. Freez 2over3** (Glide>Hall, 81 BPM). As Freez 2 of 8 with the right loop set by a tempo-coded Sw 1 patch
(16.6): 1/2 beat, 1/3 beat while Sw 1 is at its top (two against three), feedback set by ADJUST (left 100 → 200,
right 100 → 178).

**195. Freeze&Speed** (Glide>Hall, 125 BPM). ADJUST is the freeze switch (0 / 1, through the AR: attack 0,
release 56 steps): at 1 the 1/4-beat loops (left / right) take the maximum feedback, the input is cut from them,
the dry-through voices 2 / 5 go off, the reverb input opens (Rvb In 0 → 80) and the Delay Master steps 80 → 125 %;
back at 0 the release moves the Delay Master through 100 → 80 % (the loop repeats speed up) before normal operation
resumes. Sw 1 (Foot Sw 1) drives nothing.

**196. ADJToFreeze** (Glide>Hall, 81 BPM). ADJUST as a direct freeze switch: 0 = the input passes (voices 2 / 5 at
0 ms, reverb input closed); 1 = 1-beat loops (voices 1 / 4) at maximum feedback with the input cut and the reverb
fed by the loops.

**197. ADJToFreeze2** (Glide>Hall, 81 BPM). ADJToFreeze with loops of 1/3 beat (left) and 1/2 beat (right), a
cross-rhythm when frozen, and the reverb out of the mix.

**198. FSw2 Freeze** (Glide>Hall, 81 BPM). A footswitch freeze: Foot Sw 2 through a Latch (low 0 / high 100) gives
push-on / push-off; on, the 1/2-beat loops go to maximum feedback and the input is cut from them. Foot Sw 1 clears.
ADJUST sets the loop feedback when not frozen (100 → 200, stored 131).

**199. FSw2 Speedup** (Glide>Hall, 79 BPM). Freeze&Speed on Foot Sw 2: the AR (attack 0) freezes the 1/4-beat loops
at once with the Delay Master at 200 %, and on release (length set by ADJUST "How Long?": AR release 0 … 500 steps,
stored 155) the Delay Master runs down 200 → 100 → 16 %, so the frozen loop repeats faster and faster before
normal operation resumes.

**200. LongestLoop** (Glide>Hall, 120 BPM). The longest loop: the input (both channels panned left) is written into
2.5-s echo lines; voice 1 feeds the left line into the right line at maximum cross feedback. At rest Mix is 0 (dry
only) while the lines keep recording. Foot Sw 2 through a Latch (push-on / push-off) closes the loop: voice 4's cross
feedback goes to maximum, so the two 2.5-s lines form one 5-s loop, the input to the loop is cut, voice 4 is heard
and Mix goes to 100 (wet). Foot Sw 1 clears. ADJUST = FX Mix (stored 0: no reverb). With the standard delay memory
the 2.5-s lines also fit (11.1).


### 17.2 Program data (`effect-programs.json`, `effect_programs.py`)

The data file holds one record for **each of the 200 factory programs** (numbered 1 – 200), in the same
form for all of them. It is plain JSON (UTF-8); the script `effect_programs.py` (Python 3.6+, standard library only)
prints records as Markdown tables or JSON and can be imported. Section numbers in the file follow the merged
reference (1 – 10 as 1 – 10, this part's 1 – 9 as 11 – 19).

**Top level.** `programs` (the 200 records, in program order), `io_layouts` (for each input/output-section type the
rows of its quantities: label, keys, print format), `tables` (the law tables `hicut`, `locut` of 14.2, their
44.1 kHz versions `hicut44`, `locut44` (14.7), and `pitch_ratio`, `resonance`, `res2_intervals` of 14.4), `notes` (conventions, below), `fs` = 48000.

**A program record:**

* `number` (1 … 200), `name`, `algorithm`, `tempo_bpm` (the 16-bit program tempo,
  16.3), `adjust` (label, stored knob position, knob maximum);
* `description`: the text of 17.1 for this program;
* `stored`: the stored value of every parameter by edit page (`pages`: page, name, value, `locked` = tempo-locked,
  value = 32·hi + lo = num/den beats with num = lo, den = hi, 4.11 and 16.3; `max` = the parameter's maximum, the
  minimum is 0, 14.6), the Mod page (`mod`, with `max`), the trigger block of AR, Sw 1 and Sw 2 (`trigger`: source
  and level, 6.3), the Latch (`latch`: source, low and high threshold, 16.5) and the patches (`patches`: source,
  page-qualified destination, curve points [x, y], 4.10; sources 16.1; `tempo_coded`: the points are beat codes
  32·num + den, also given as `points_beats`, 16.6). The effect voices V1 … V6 of 13.2 – 13.4 are
  Voice1 … Voice6 (Res: V1 … V6) of the Levels, DelayTime, Feedback, Filters/Resonance and Panning pages;
* `parameters`: every parameter with its stored value, its `max` and, where a patch changes it at rest, `effective`
  and `via` (`ADJUST` at the stored knob position, `patch at rest`, or `patches (sum)` when several patches add on
  one destination); a value set by a tempo-coded patch also as `effective_beats` (num/den, 16.6); Mid Rt also as RT[v]
  in seconds (`rt_s`);
* `modulation`: the state of the moving parts (AP Sweep, Concert Hall / Glide>Hall modulated taps);
* `tank`: `program` (the tank listing to use: 5.1 Chamber, 5.2 Plate, 5.3 Concert Hall, 5.4 Inverse, or 15.1
  Infinite), `k` (coefficient k<s> of listing step s) and `A` (program-dependent offset A<s> in samples). Everything
  else in a tank is the same for all its programs and appears as a number in the listing;
* `io`: `section` (the input/output section: 3, 3.1, 13.2, 13.3 or 13.4) and `values`, exact: coefficient words as
  word/2^23 (11), delays in samples;
* AutoInfinite (program 46): `tank` uses the notation of 15.1 (k = whole coefficients, integer part included). It gives the
  rest state, Infinite = 1 (the AR patch gives 1 at source value 0: frozen), and the running state, Infinite = 0, as
  `k_infinite0` / `A_infinite0`.
* `part1_section`: for the 34 programs tabulated in 7 (151 – 183 and BowedChords), the 7 subsection of their
  group; null for the others. Their coefficient values were taken from 7 and set to the exact coefficient words of
  the program.

Conventions (as 7):

* AP Sweep pairs at their rest gain with the stock signs (6.1); programs marked "sweeping" keep modulating them.
* Concert Hall / Glide>Hall modulated tap coefficients and offsets are a snapshot (6.2).
* Patches from sources other than ADJUST are evaluated at source value 0 (input silent, envelopes at rest, LFO at the
  start of its curve); such values are marked "patch at rest". The MIDI sources (and "none") have no value after
  power-up (16.1): their patches leave the destination at its stored value. Several patches on one destination add
  ("patches (sum)"). Every value follows the exact patch law of 4.10, including the destination's resolution, so an
  `effective` value can be fractional (Glide A/B Dly: 1/128 step).
* Every effect voice is given with its gate open (16.4); with Clear = 1 at rest its gated words are 0. Exception: BowedChords is
  given with its AR envelope at the top (input level 0 dB), 7.5.
* Delays of the effect lines are given for the extended memory (r = 32) and, with keys ending in `_std`, for the
  standard memory (r = 2); 11.1.
* Res2 tuning words P_k are those in effect at rest (at the stored ADJUST position). Res1 programs assign their
  resonators only when Assign changes after the load (13.4, 14.4): at rest every Res1 resonator shows the unassigned
  word 6 011 780 with feedback 0.

**Configuring a simulation from a record.**

1. Tank: run the listing named in `tank.program` with the record's `k` and `A` (5, or 15.1 for Infinite).
2. Input/output section: run the listing of `io.section` (3, 3.1 with the correction of 13.1, or 13.2 – 13.4) with
   `io.values`; the key names are the symbols of the listing (table below).
3. Control layer (20 ms steps, 6.3): the modulators from `stored.mod`, `stored.trigger` and `stored.latch`, the
   patches from `stored.patches` evaluated with the patch law (4.10); a patch result goes through the parameter law of
   its destination (the laws of 4 and of section 14) and replaces the value in `io` / `tank` as described in 6.4 and 16.4. Tempo-locked values
   follow 4.11 and 16.3; a patch into a tempo-locked parameter follows 16.6 (key `tempo_coded`).
4. The ADJUST knob, external sources (MIDI, pedal, footswitches) and the tempo are the simulation's inputs.

**Commands.** PROG is a program number (`72`) or a program name or a unique part of
it (case-insensitive: `"tape echo"`, `modal`). ALG is part of an algorithm name (`m-band`, `glide`, `res2`, `chamber`).

| command | prints |
|---|---|
| `python effect_programs.py list [ALG] [--long]` | the programs (number, name, algorithm, ADJUST label); `--long`: their descriptions |
| `python effect_programs.py show PROG …` | everything about a program, as Markdown tables |
| `python effect_programs.py compare PROG …` | several programs of one algorithm family side by side (parameters in effect, modulation, k, A, patches, input/output section) |
| `python effect_programs.py compare --alg ALG` | all programs of an algorithm, in tables of seven |
| `python effect_programs.py value KEY [ALG]` | one quantity for every program: an input/output key (`D1`, `P3`, `gs`), `k12`, `A40`, or a parameter (`Mid Rt`, `Levels:Voice1`) |
| `python effect_programs.py keys PROG` | the keys available for a program |
| `python effect_programs.py json PROG` | the program's complete record as JSON |
| `python effect_programs.py table NAME` | the tables of the laws: `hicut`, `locut` (14.2), `hicut44`, `locut44` (14.7); `pitch_ratio`, `resonance`, `res2_intervals` (14.4) |

From Python (e.g. to configure a simulation):

```
import effect_programs as E
p = E.program('Tape Echo')                    # one record (dict)
p['io']['values']['D1'], p['io']['values']['fb1']    # input/output-section quantities of 13.2
p['tank']['k']['13'], p['tank']['A']['40']   # tank coefficient k13, offset A40
E.value(p, 'Mid Rt')                          # effective parameter value
[q['name'] for q in E.programs('res2')]       # all Res2>Plate programs
E.res_word(57, 100)                           # tuning word of A4 at Tuning 100 (14.4): 223418
E.res2_note(36, 37)                           # Res2 voice note: Assign C3, interval +3rd -> 40 (E3)
E.res_feedback(91, 190)                       # loop feedback n (b = n/32768) for Mstr Res 91, V Res 190
E.voice_delay(133, 'M-Band+Rvb', locked=True, bpm=94)   # effect voice delay (14.1): 5/4 beat -> 38304 samples
E.beat_samples(94, 5, 4)                      # tempo-locked length S (16.3): 38304
```

Without the script: `json.load(open('effect-programs.json', encoding='utf-8'))['programs'][i]`.

**Keys of the input/output section.** The keys are the symbols of the signal-flow listings, written in ASCII:

| section | keys |
|---|---|
| 3 (Chamber, Plate, Concert Hall, Inverse, Infinite) | `aLL aLR aRL aRR` send matrix, `gs` send gain, `gr` return gain, `postL postR` post-delays (samples), `pd pr` post mix, `c1 c2` reverb width, `g` voice diffusion, `vd vp` voice input, `DV1 … DV4` voice delays (samples), `V1L V1R V1fb` … `V4L V4R V4fb` voice levels and feedback, `fv fr` FX mix, `a` output low-pass, `w1 w2` output width, `out_gain` = 4·go, `wet dry` |
| 3.1 with 13.1 (Chorus+Rvb) | `c90 … c94` send, `c96 … c9B` voice feedback, `c9C` line damping, `B1 … B6` delay targets (1/128 sample), `g` delay glide, `A1 … A6` LFO depth targets (×2^−23), `R1` … `R6` LFO base rates, `c38` depth smoothing, `U u N` Weyl step / start / counter start, `cL1 … cL6`, `cR1 … cR6` voice levels, `cAC cAD` FX mix, `cAE cAF` width, `out_gain` = 4·cB0, `wet dry` |
| 13.2 (M-Band+Rvb) | `aLL aLR aRL aRR`, `gs`, `fb1 … fb6`, `D1 … D6` (and `D1_std …`), `g1 … g6`, `hicut1 … hicut6` / `locut1 … locut6` (indices into the HiCut / LoCut tables: H and L sets of 13.2), `pL1 … pL6`, `pR1 … pR6`, `fv fr`, `w1 w2`, `out_gain` = 4·go, `wet dry` |
| 13.3 (Glide>Hall) | `aLL aLR aRL aRR`, `pAL pBL pAR pBR` glide targets P (samples), `rhoAL …` ρ, `cL cR` tap form, `lamAL lamBL lamAR lamBR` λ, `kapLL kapLR kapRR kapRL` κ (L→L, L→R, R→R, R→L), `gg`, `D1 … D6` (and `_std`), `sf1 … sf6` / `sx1 … sx6` echo feedback, `g vd vp` voice diffusion, `eL1 … eL6`, `eR1 … eR6`, `gs`, `fv fr`, `w1 w2`, `out_gain`, `wet dry` |
| 13.4 (Res1>Plate, Res2>Plate) | `aLL aLR aRL aRR`, `psi1 … psi6` ψ, `D1 … D6` (and `_std`), `P1 … P6` tuning words, `a1 … a6` loop low-pass, `b1 … b6` loop feedback, `rL1 … rL6`, `rR1 … rR6`, `gs`, `fv fr`, `w1 w2`, `out_gain`, `wet dry` |


## 18. Accuracy

The effect sections (13.2 – 13.5), the Infinite tank (15.1) and every law marked **exact** reproduce the original for
every value of the parameters involved. The remaining limits:

* Some 16-bit laws are exact only to within 1/32768 (one least significant bit), as marked in 4.
* The Glide>Hall Voice Dif law (14.3) is approximate; the echo feedback law (14.3) is exact.
* **Freeze and loop programs** (Super Ball!, Freez 2 of 8, Freez 2over3, Freeze&Speed, ADJToFreeze, ADJToFreeze2,
  FSw2 Freeze, FSw2 Speedup, LongestLoop): the frozen loop runs at the maximum echo feedback, 0.99994 per pass
  (14.3), i.e. it loses 0.0005 dB per pass (a 1/2-beat loop at 81 BPM: 60 dB in about 11 hours); saturation (1)
  bounds it.
* Infinite: the tank listing (15.1) is exact for the factory program; FlyTapRng / FlyTapDly are not described.
* **Modulator time base:** one control step is 20 ms (16.3).
* The effect-line resolution depends on the delay memory (11.1).

## 19. Presenting the parameters in a simulation

9 applies to the reverb part of every algorithm. For the effect sections:

* **Layer 1 (primary):** M-Band band delay times (or beats), band HiCut/LoCut (as corner frequencies) and feedback;
  Glide tap delays (A/B Dly) with Gld Lvl and feedback; Res tuning (Assign / key) and Mstr Res; Chorus depth and rate
  masters; FX Mix.
* **Layer 2:** per-voice levels and pans, the masters (Delay, Level, Pan, Feedback, HC, LC), Voice Dif, glide
  response.
* **Layer 3:** the six individual voices, cross-feedback, the trigger and modulator settings.
* **Couplings to show:** delay times locked to tempo (16.3); the 32-sample delay resolution of the extended memory
  (11.1); M-Band and Chorus lines take the send, so Rvb In-like gains and the input matrix affect both the effect and
  the reverb (13.1); in Glide>Hall and Res the reverb is fed by the effect, not by the input (13.3, 13.4).
* **Not to present:** FlyTapRng / FlyTapDly (not stored), Rvb In of M-Band and Chorus+Rvb (fixed).
