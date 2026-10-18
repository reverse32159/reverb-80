# reverb-80

A complete, simulation-grade description of a classic stereo digital reverb/effects unit: its ten algorithms and all
200 factory programs, written so that the unit can be re-implemented sample by sample. The files do not implement
the reverb unit: they contain no signal processing code, and it is not possible to run a simulation of any of the
programs, or to process audio, with them directly.

The repository contains one reference document, one data file with every factory program, and a small Python script
that reads the data.

## Files

| file | what it is |
|---|---|
| [`algorithms-and-programs.md`](algorithms-and-programs.md) | the reference: number formats, building blocks, signal flow, exact step-by-step listings of every reverb tank, all parameter laws, the modulation and control layer, tempo rules, a description of every program and notes on accuracy |
| [`effect-programs.json`](effect-programs.json) | the program data: one record per factory program (200) with its stored parameters, patches and modulator settings, and the resulting tank coefficients, delays and input/output-section values |
| [`effect_programs.py`](effect_programs.py) | command-line tool and Python module for the data: lists programs, prints a program as Markdown tables, compares programs, lists one value across all programs, and computes the main laws (delays, tempo-locked lengths, resonator tuning) |

## The reference document

`algorithms-and-programs.md` has two parts.

**Part I — reverb algorithms (sections 1 – 10)**

| section | content |
|---|---|
| 1 | conventions: number formats, coefficient words, notation of the listings |
| 2 | building blocks: allpass, interpolated taps, saturation |
| 3 | signal flow of the reverb algorithms and the shared input/output section; 3.1 the Chorus+Rvb section |
| 4 | parameter laws (stored value → coefficient or delay): Size, decay, damping, levels, High Cut tables, Shape/Attack, Spread, Inverse, patches and the ADJUST knob, tempo locks, post-delays and voices, Glide/Clear, ranges, operation at 44.1 kHz |
| 5 | the tank programs as exact listings: Chamber, Plate, Concert Hall, Inverse |
| 6 | modulation: AP Sweep, Concert Hall modulated taps, the modulation sources (LFO, envelopes, AR, switches), parameter changes while running |
| 7 | the reverb programs: stored data and notes, with a full worked description of BowedChords |
| 8 | accuracy: what is exact and what is approximate |
| 9 | presenting the parameters in a simulation (edit layers, couplings, display units) |
| 10 | the analog input path (drive and clipping), optional |

**Part II — effect algorithms and the factory programs (sections 11 – 19)**

| section | content |
|---|---|
| 11 | conventions; 11.1 delay memory size (standard vs extended) and the delay-line resolution |
| 12 | building blocks of the effect sections, including tap cross-fades |
| 13 | signal flow of M-Band+Rvb, Glide>Hall, Res1>Plate / Res2>Plate and Infinite |
| 14 | parameter laws of the effect algorithms, ranges, operation at 44.1 kHz |
| 15 | the Infinite tank as an exact listing |
| 16 | patch sources, AR envelope modes, tempo (including MIDI clock), effect-section changes while running, Latch, patches into tempo-locked parameters |
| 17 | the programs: 17.1 a description of every one of the 200 programs, 17.2 the format of the program data and how to use it |
| 18 | accuracy |
| 19 | presenting the effect parameters in a simulation |

### Suggested reading order for building a simulator

1. Sections 1 – 3: number formats, building blocks, signal flow.
2. Section 5 (and 15.1 for Infinite): the tank listings.
3. Sections 4 and 14: the parameter laws, needed for live edits and modulation.
4. Sections 6 and 16: the 20 ms control layer: LFO, envelopes, switches, Latch, patches, tempo.
5. Sections 11 – 13: the effect algorithms and the delay memory.
6. Section 17: the programs.
7. Sections 8 and 18 (accuracy and limits), 9 and 19 (user interface), 10 (analog input, optional).

## The program data

`effect-programs.json` is plain UTF-8 JSON. Its top level holds `programs` (the 200 records, numbered 1 … 200),
`io_layouts` (the quantities of each input/output-section type), `tables` (law tables such as `hicut`, `locut`,
`pitch_ratio`), `notes` (the conventions of the data) and `fs` = 48000.

A program record contains:

* `number` (1 … 200), `name`, `algorithm`, `tempo_bpm`, `adjust` (the ADJUST knob's label and
  position);
* `description`: the program's text from section 17.1;
* `stored`: every stored parameter by edit page, the modulator settings, triggers, Latch and patches;
* `parameters`: every parameter with its stored value, its maximum and, where a patch changes it at rest, the
  effective value;
* `tank`: which tank listing to run (`program`) and its program-dependent coefficients `k` and offsets `A`;
* `io`: which input/output-section listing to run (`section`) and its exact `values` (coefficients, delays in
  samples).

The programs are given at rest: the ADJUST knob at its stored position, every other modulation source at 0. Effect
delays are given for the extended delay memory; keys ending in `_std` give the standard-memory values. Section 17.2
documents every field and convention.

| algorithm | programs |
|---|---|
| M-Band+Rvb | 50 |
| Chamber | 38 |
| Chorus+Rvb | 35 |
| Plate | 25 |
| Glide>Hall | 18 |
| Concert Hall | 16 |
| Inverse | 7 |
| Res1>Plate | 6 |
| Res2>Plate | 4 |
| Infinite | 1 |

## Using the script

`effect_programs.py` needs Python 3.6 or later and only the standard library. Keep it in the same folder as
`effect-programs.json`.

A program (PROG) can be given by its number (`72`) or by name or a unique part of it, case-insensitive
(`"tape echo"`, `modal`). An algorithm (ALG) is part of an algorithm name (`m-band`, `glide`, `res2`, `chamber`).

```
python effect_programs.py list                     # all 200 programs
python effect_programs.py list res2 --long         # the Res2>Plate programs with their descriptions
python effect_programs.py show "Tape Echo"         # one program: stored data, values in effect, tank, io
python effect_programs.py compare --alg chamber    # all Chamber programs side by side, in tables of seven
python effect_programs.py value "Mid Rt" plate     # one parameter (or io key, k13, A40) across programs
python effect_programs.py keys 72                  # the value keys available for a program
python effect_programs.py json 157                 # the complete record as JSON
python effect_programs.py table hicut              # a law table
```

As a module, for example to configure a simulator:

```python
import effect_programs as E

p = E.program('Tape Echo')                  # one record (dict)
p['io']['values']['D1']                     # an input/output-section value
p['tank']['k']['13'], p['tank']['A']['40']  # tank coefficient k13 and offset A40
E.value(p, 'Mid Rt')                        # a parameter's effective value
E.programs('res2')                          # all Res2>Plate records

E.voice_delay(133, 'M-Band+Rvb', locked=True, bpm=94)   # effect voice delay in samples (section 14.1)
E.beat_samples(94, 5, 4)                                # tempo-locked length (section 16.3)
E.res_word(57, 100)                                     # resonator tuning word (section 14.4)
```

Without the script: `json.load(open('effect-programs.json', encoding='utf-8'))['programs'][i]`.

## Building a simulation from a record

1. **Tank:** run the listing named in `tank.program` (section 5, or 15.1 for Infinite) with the record's `k` and `A`.
2. **Input/output section:** run the listing named in `io.section` (3, 3.1, or 13.2 – 13.4) with `io.values`; the
   keys are the symbols of the listing.
3. **Control layer** (20 ms steps): the modulators, Latch and patches from `stored`, evaluated with the patch law
   (4.10) and passed through each destination's parameter law (sections 4 and 14). Tempo-locked values follow 4.11
   and 16.3.
4. **Inputs of the simulation:** the ADJUST knob, external controllers (MIDI, pedal, footswitches) and the tempo.

Two choices are left to the implementer:

* **Delay memory.** The extended memory (voice delays up to about 43.7 s) is the reference and the default in the
  data; the standard memory (about 2.4 – 2.7 s) suits a small target (section 11.1).
* **Sample rate.** 48 kHz is the reference; 44.1 kHz is fully covered (sections 4.15 and 14.7).

## Accuracy

The tank listings, the input/output sections, the control layer and almost all parameter laws are exact: they
reproduce the original for every parameter value. Sections 8 and 18 list the few approximations.
