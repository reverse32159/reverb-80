"""Program data of *Reverb and effect algorithms and the factory programs — simulation reference*
(algorithms-and-programs.md, section 17.2).

Reads effect-programs.json (same folder): all 200 factory programs. Prints the per-program tables: stored data,
Latch, patches, parameters in effect, modulation state, tank coefficients k and offsets A and the
input/output-section values, for one program, several side by side, or one value across programs. Section numbers
refer to algorithms-and-programs.md. Python 3.6+, standard library only.

usage:
  python effect_programs.py list [ALG] [--long]       programs (ALG: part of an algorithm name, e.g. res2, m-band)
  python effect_programs.py show PROG [PROG ...]      everything about each program (Markdown)
  python effect_programs.py compare PROG [PROG ...]   side-by-side tables (as the document's 7.x tables)
  python effect_programs.py compare --alg ALG         all programs of an algorithm, in groups of 7
  python effect_programs.py value KEY [ALG]           one value for every program (io key, k12, A40, [page:]parameter)
  python effect_programs.py keys PROG                 the value keys of a program
  python effect_programs.py json PROG                 the program's record as JSON (for a simulator)
  python effect_programs.py table NAME                tables: hicut, locut (14.2), hicut44, locut44 (14.7);
                                                      pitch_ratio, resonance, res2_intervals (14.4)

PROG is a program number 1..200 ("72") or a program name (case-insensitive; a unique part of a name is enough,
e.g. "tape echo", "modal").

In Python:  import effect_programs as E;  p = E.program('Tape Echo');  p['io']['values']['D1'], p['tank']['k']['13']
            E.res_word(note, tuning, fs), E.res2_note(assign, iv, scale, key, root, rule), E.res_feedback(mstr, vres),
            E.voice_delay(value, algorithm, locked, bpm, master, extended, fs), E.beat_samples(bpm, num, den, fs)
"""
import json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'effect-programs.json')
PARAMS = ['Size', 'Duration', 'Mid Rt', 'Rt HC', 'Low Rt', 'Crossover', 'Diffusion', 'Pre Delay', 'Spread', 'Shape', 'Attack',
          'Link', 'AP Sweep', 'Rvb Out', 'Def', 'Depth', 'Chorus', 'Infinite', 'RefDly L', 'RefDly R', 'RefLvl L', 'RefLvl R',
          'EkoDly L', 'EkoDly R', 'EkoFbk L', 'EkoFbk R', 'Rvb In', 'FX Mix', 'Rvb Width', 'FX Width', 'High Cut', 'FX Adjust']
_D = None


def data():
    global _D
    if _D is None:
        with open(DATA, encoding='utf-8') as f: _D = json.load(f)
    return _D


def programs(alg=None):
    """all program records, optionally only those whose algorithm name contains alg (case-insensitive)"""
    ps = data()['programs']
    return [p for p in ps if alg is None or alg.lower() in p['algorithm'].lower()]


def program(q):
    """one program record by number (1..200) or (part of) its name"""
    ps = data()['programs']; s = str(q).strip()
    if s.isdigit():
        r = [p for p in ps if p['number'] == int(s)]
    else:
        r = [p for p in ps if p['name'].lower() == s.lower()] or [p for p in ps if s.lower() in p['name'].lower()]
    if len(r) == 1: return r[0]
    if not r: raise KeyError('no program %r in %s' % (q, os.path.basename(DATA)))
    raise KeyError('%r is ambiguous: %s' % (q, ', '.join('%s (%d)' % (p['name'], p['number']) for p in r)))


def value(p, key):
    """io value, tank k<s> / A<s>, or parameter in effect (by name)"""
    v = p['io']['values'] if p['io'] else {}
    if key in v: return v[key]
    if key[:1] in 'kA' and key[1:].isdigit() and p['tank']:
        return p['tank']['k' if key[0] == 'k' else 'A'].get(key[1:])
    for e in p['parameters']:
        if key.lower() in (e['name'].lower(), ('%s:%s' % (e['page'], e['name'])).lower()):
            return e.get('effective', e['stored'])
    return None


# ------------------------------------------------------------------------------- resonator laws (document 14.4)
def res_word(n, tuning=100, fs=48000):
    """tuning word P (period in 1/2048 sample) of note n (0 = C0, 57 = A4); Tuning 0..200 = 430.0..450.0 Hz"""
    T = data()['tables']['pitch_ratio']['values']
    q = 10 * ((fs << 15) // (tuning + 4300))
    w = (q * T[(n - 57) % 12] + 0x4000) >> 15
    oc = 4 + (n - 57) // 12
    return w << -oc if oc < 0 else w >> oc


RANGE_MS = {'M-Band+Rvb': (43690, 2730), 'Glide>Hall': (42325, 2645), 'Res1>Plate': (38229, 2389),
            'Res2>Plate': (38229, 2389)}


def beat_samples(bpm, num=1, den=1, fs=48000):
    """tempo-locked length S of num/den beats at an integer program tempo (16.3): one rounding"""
    R = 4096000 // (960000 // bpm)
    return (256 * fs * num) // (R * den)


def voice_delay(value, algorithm, locked=False, bpm=120, master=100, extended=True, fs=48000):
    """effect voice delay in samples (14.1, 11.1, 14.7): value = stored ms, or the code 32*den + num when locked;
    master = Delay Master (%); extended = memory with r = 32 (False: standard memory, r = 2)"""
    T = RANGE_MS[algorithm][0 if extended else 1]
    r = 32 if extended else 2
    if locked:
        t = 1000 * beat_samples(bpm, value & 31, value >> 5, fs) // fs
        while t > T: t //= 2
    else:
        t = min(value, T)
    t = min(t * master // 100, T)
    return r * (fs * t // 1000 // r)


def res2_note(assign, iv, scale=0, key=0, root=0, rule=0):
    """note of a Res2 voice: Assign note, interval index iv (35 = unison, 7 steps per octave), Scale 0/1, Key 0..11,
    Root 0..6, Rule 0..3. A result of 0 leaves the voice unchanged."""
    d = iv % 7
    off = data()['tables']['res2_intervals']['values'][scale][root][rule][d - 1][(assign + 12 - key) % 12] if d else 0
    n = max(0, assign + (iv // 7 - 5) * 12 + off)
    while n > 115: n -= 12
    return n


def _rm(a, b): return ((a * b) >> 15) + (((a * b) >> 14) & 1)


def res_feedback(mstr, vres):
    """loop feedback b of a voice as a 16-bit value n (b = n/32768): Mstr Res 0..100, V Res 0..200 (100 = 0)"""
    t = data()['tables']['resonance']; R = t['values']
    v = R[vres - 100] if vres >= 100 else -R[100 - vres]
    return _rm(_rm(R[mstr], v), t['scale'])


# ------------------------------------------------------------------------------------------------- formatting
def _param_cell(e):
    st = e['stored']; ef = e.get('effective')
    t = ('%d → **%s**' % (st, repr(ef) if isinstance(ef, float) else ef)) if ef is not None else str(st)
    if 'rt_s' in e: t += ' (%.2f s)' % e['rt_s']
    if 'effective_beats' in e: t += ' (%s beat)' % e['effective_beats']
    if e.get('via') == 'patch at rest': t += ' (patch at rest)'
    return t


def _curve(pt):
    if pt.get('tempo_coded'):          # points 32*num + den = num/den beats (16.6)
        return ' '.join('%d:%s beat' % tuple(x) for x in pt['points_beats'])
    return ' '.join('%d:%d' % tuple(x) for x in pt['points'])


def _patch(pt): return '%s → %s [%s]' % (pt['source'], pt['destination'], _curve(pt))


def _fmt(v, fmt):
    if isinstance(v, str): return v
    return fmt % v


def lpf_hz(a):
    """-3 dB frequency of the one-pole low-pass y += a·(x − y) at 48 kHz (approximation used in the tables)"""
    return 'flat' if a >= 0.999 else '%.0f Hz' % (-math.log(1 - a) * 48000 / (2 * math.pi))


def _io_cell(p, r):
    v = p['io']['values']; t = ' '.join(_fmt(v[k], r['fmt']) for k in r['keys'])
    if r['label'].startswith('output low-pass a'): t += ' (%s)' % lpf_hz(v['a'])
    return t


def _io_rows(sec):
    return data()['io_layouts'][sec]


def stored_table(p):
    s = p['stored']; L = ['| page | stored values |', '|---|---|']
    for g in s['pages']:
        L.append('| %s | %s |' % (g['page'], ', '.join('%s %d%s' % (q['name'], q['value'], '*' if q['locked'] else '') for q in g['params'])))
    L.append('| Mod | %s |' % ', '.join('%s %d%s' % (q['name'], q['value'], '*' if q['locked'] else '') for q in s['mod']))
    L.append('| trigger | %s |' % '; '.join('%s: %s ≥ %d' % (n, t['source'], t['level']) for n, t in s['trigger'].items()))
    if 'latch' in s:
        la = s['latch']
        L.append('| latch | %s |' % ('none' if la['source'] == 'none' else
                                     '%s, low %d, high %d' % (la['source'], la['low'], la['high'])))
    L.append('| patches | %s |' % ('; '.join(_patch(x) for x in s['patches']) or '—'))
    return '\n'.join(L)


def show(p):
    a = p['adjust']
    L = ['## %s (%s; program %d; tempo %d BPM; ADJUST %d, label "%s")' % (p['name'], p['algorithm'], p['number'], p['tempo_bpm'],
                                                                  a['position'], a['label']), '', p['description'], '',
         '### Stored data (`*` = tempo-locked, value = 32·hi + lo)', '', stored_table(p), '',
         '### Parameters in effect (stored → effective)', '', '| parameter | value |', '|---|---|']
    L += ['| %s | %s |' % ('%s:%s' % (e['page'], e['name']) if e['page'] else e['name'], _param_cell(e)) for e in p['parameters']]
    if p['io'] is None:
        L += ['', 'Parameters in effect, modulation state, tank coefficients and input/output-section values: Part I '
                  'section %s.' % p['part1_section']]
        return '\n'.join(L) + '\n'
    L += ['', '### Modulation state', '', '| | |', '|---|---|'] + ['| %s | %s |' % kv for kv in p['modulation'].items()]
    t = p['tank']
    L += ['', '### Tank coefficients k (%s, effective values)' % t['program'], '', '| k | value |', '|---|---|']
    L += ['| k%s | %.6f |' % kv for kv in t['k'].items()]
    L += ['', '### Program-dependent offsets A (samples)', '', '| A | value |', '|---|---|']
    L += ['| A%s | %d |' % kv for kv in t['A'].items()]
    L += ['', '### Input/output section values (%s)' % p['io']['section'], '', '| quantity | values |', '|---|---|']
    for r in _io_rows(p['io']['section']):
        L.append('| %s | %s |' % (r['label'], _io_cell(p, r)))
    return '\n'.join(L) + '\n'


def compare(ps):
    t1 = [p for p in ps if p['io'] is None]
    if t1:
        raise ValueError('tabulated in Part I: ' + ', '.join('%s (section %s)' % (p['name'], p['part1_section']) for p in t1))
    secs = {p['io']['section'] for p in ps}
    if len(secs) > 1:
        raise ValueError('compare needs programs of one input/output section (algorithm family); got ' + ', '.join(sorted(secs)))
    hdr = '| | ' + ' | '.join(p['name'] for p in ps) + ' |'; sep = '|---|' + '---|' * len(ps)
    L = ['#### Parameters in effect', '', hdr, sep]
    for n in PARAMS:
        cells = []
        for p in ps:
            e = [x for x in p['parameters'] if x['name'] == n]
            cells.append(_param_cell(e[0]) if e else '')
        if any(cells): L.append('| %s | %s |' % (n, ' | '.join(cells)))
    keys = []
    for p in ps:
        for k in p['modulation']:
            if k not in keys: keys.append(k)
    if keys:
        L += ['', '#### Modulation state', '', hdr, sep]
        L += ['| %s | %s |' % (k, ' | '.join(p['modulation'].get(k, '') for p in ps)) for k in keys]
    tk = {p['tank']['program'] for p in ps}
    L += ['', '#### Tank coefficients k (%s, effective values)' % ' / '.join(sorted(tk)), '', hdr, sep]
    ks = sorted({int(s) for p in ps for s in p['tank']['k']})
    L += ['| k%d | %s |' % (s, ' | '.join(('%.6f' % p['tank']['k'][str(s)]) if str(s) in p['tank']['k'] else '' for p in ps)) for s in ks]
    L += ['', '#### Program-dependent offsets A (samples at 48 kHz)', '', hdr, sep]
    As = sorted({int(s) for p in ps for s in p['tank']['A']})
    L += ['| A%d | %s |' % (s, ' | '.join(str(p['tank']['A'].get(str(s), '')) for p in ps)) for s in As]
    L += ['', '#### Patches', '', hdr, sep]
    L.append('| ADJUST knob (label = stored position) | %s |' % ' | '.join('%s = %d' % (p['adjust']['label'], p['adjust']['position']) for p in ps))
    L.append('| patches (source → destination [x:y …]) | %s |' % ' | '.join('<br>'.join(_patch(x) for x in p['stored']['patches']) or '—' for p in ps))
    sec = secs.pop()
    L += ['', '#### Input/output section values (%s)' % sec, '', hdr, sep]
    for r in _io_rows(sec):
        L.append('| %s | %s |' % (r['label'], ' | '.join(_io_cell(p, r) for p in ps)))
    return '\n'.join(L) + '\n'


def table(name):
    t = data()['tables'][name]
    if name == 'res2_intervals':
        L = ['# ' + t['note'], '# scale root rule degree: offsets for pitch class 0..11 (re Key)']
        for sc, a in enumerate(t['values']):
            for ro, b in enumerate(a):
                for ru, c in enumerate(b):
                    for d, row in enumerate(c):
                        L.append('%d %d %d %d: %s' % (sc, ro, ru, d + 1, ' '.join('%2d' % x for x in row)))
        return '\n'.join(L)
    if 'values' in t:
        return '# %s\n%s' % (t['note'], '\n'.join('%5d %6d' % (i, v) for i, v in enumerate(t['values'])))
    L = ['index ' + ' '.join('%10s' % c for c in t['columns'])]
    L += ['%5d %s' % (i, ' '.join('%10d' % w for w in row)) for i, row in enumerate(t['rows'])]
    return '\n'.join(L)


def main(av):
    if not av or av[0] in ('-h', '--help', 'help'):
        print(__doc__); return 0
    cmd, args = av[0], av[1:]
    if cmd == 'list':
        long_ = '--long' in args; args = [a for a in args if a != '--long']
        for p in programs(args[0] if args else None):
            if long_: print(p['description'] + '\n')
            else: print('%3d  %-13s %-12s %s' % (p['number'], p['name'], p['algorithm'], p['adjust']['label']))
    elif cmd == 'show':
        print('\n'.join(show(program(q)) for q in args))
    elif cmd == 'compare':
        if args[:1] == ['--alg']:
            ps = [p for p in programs(' '.join(args[1:])) if p['io']]
            for i in range(0, len(ps), 7): print(compare(ps[i:i + 7]))
        else:
            print(compare([program(q) for q in args]))
    elif cmd == 'value':
        for p in programs(args[1] if len(args) > 1 else None):
            v = value(p, args[0])
            if v is not None: print('%3d  %-13s %s' % (p['number'], p['name'], v))
    elif cmd == 'keys':
        p = program(args[0])
        if p['io']:
            print('io (%s): %s' % (p['io']['section'], ' '.join(p['io']['values'])))
            print('tank: %s' % ' '.join(['k' + s for s in p['tank']['k']] + ['A' + s for s in p['tank']['A']]))
        else:
            print('io / tank: Part I section %s' % p['part1_section'])
        print('parameters: %s' % ', '.join('%s:%s' % (e['page'], e['name']) if e['page'] else e['name'] for e in p['parameters']))
    elif cmd == 'json':
        print(json.dumps(program(args[0]), ensure_ascii=False, indent=1))
    elif cmd == 'table':
        print(table(args[0].lower()))
    else:
        print(__doc__); return 2
    return 0


if __name__ == '__main__':
    try:
        sys.stdout.reconfigure(encoding='utf-8')
    except AttributeError:
        pass
    try:
        sys.exit(main(sys.argv[1:]))
    except (KeyError, ValueError) as e:
        print('error:', e.args[0] if e.args else e, file=sys.stderr); sys.exit(1)
