"""Benchmark summary over all cases complete on both arms (<project_root>/results/<S>/, copied from RESULTS on the cluster).
Per case: read pairs, per-arm step times and speedups, flagstat identity, read-level BAM comparison
(<S>_bam_comparison.txt, if present), and CPU-vs-GPU VCF concordance (downsampling off) genome-wide and on target
(IDT V5 exons +/-10 bp), raw and PASS. Writes to analysis/tables/:
  bench_summary.tsv         one row per case
  bench_timing_long.tsv     sample, read_pairs, arm, step, seconds
  bench_target_calls.tsv    on-target union of calls: CPU and GPU filter, DP, AF (for VAF agreement)
  bench_discordant.tsv      genome-wide calls present in one arm only: arm, filter, DP, alt reads, AF, on_target
  bench_dp_all.tsv          depth distribution of all CPU calls (DP bin, count), the background for discordant calls
Usage: python -I analysis/bench_summary.py <project_root>"""
import gzip, bisect, sys, os, re, csv
from collections import defaultdict
ROOT = sys.argv[1]
SNAP = f'{ROOT}/results'

bed = defaultdict(list)
for l in open(f'{ROOT}/analysis/tables/panel_targets.bed'):
    c, s, e = l.split()[:3]; bed[c].append((int(s), int(e)))
for c in bed: bed[c].sort()
starts = {c: [s for s, _ in v] for c, v in bed.items()}
def on_target(c, p):
    i = bisect.bisect_right(starts.get(c, []), p - 1) - 1
    return i >= 0 and p <= bed[c][i][1]

def load(f):
    d = {}
    for l in gzip.open(f, 'rt'):
        x = l.rstrip('\n').split('\t'); d[(x[1], int(x[2]), x[3], x[4])] = tuple(x[5:11])
    return d
def timing(f):
    return {r['step']: int(r['seconds']) for r in csv.DictReader(open(f), delimiter='\t')}
J = lambda a, b: len(a & b) / len(a | b) if a | b else 1.0

rows, tlong, tcalls, disc = [], [], [], []
dp_all = defaultdict(int)
DP_BINS = [1, 2, 3, 5, 10, 20, 50, 100, 200, 500, 1000, 10**9]
def dp_bin(d):
    return next(b for b in DP_BINS if d < DP_BINS[DP_BINS.index(b) + 1]) if d < DP_BINS[-2] else DP_BINS[-2]
for S in sorted(os.listdir(SNAP)):
    p = f'{SNAP}/{S}/{S}'
    if not all(os.path.exists(f'{p}_{a}_nods.calls.tsv.gz') and os.path.exists(f'{p}_{a}_timing.tsv') for a in ('cpu', 'gpu')):
        continue
    tc, tg = timing(f'{p}_cpu_timing.tsv'), timing(f'{p}_gpu_timing.tsv')
    cpu_bam = tc['align_sort'] + tc['markdup'] + tc['index']
    fs = [open(f'{p}_{a}_flagstat.txt').read() for a in ('cpu', 'gpu')]
    r = {'sample': S, 'read_pairs': int(re.search(r'(\d+) \+ \d+ read1', fs[0]).group(1)),
         'cpu_bam_s': cpu_bam, 'gpu_bam_s': tg['fq2bam'], 'bam_speedup': round(cpu_bam / tg['fq2bam'], 1),
         'cpu_m2_s': tc['mutect2'], 'gpu_m2_s': tg['mutect2'], 'm2_speedup': round(tc['mutect2'] / tg['mutect2'], 1),
         'cpu_total_s': cpu_bam + tc['mutect2'], 'gpu_total_s': tg['fq2bam'] + tg['mutect2']}
    r['total_speedup'] = round(r['cpu_total_s'] / r['gpu_total_s'], 1)
    r['flagstat_identical'] = fs[0] == fs[1]
    bc = f'{p}_bam_comparison.txt'
    if os.path.exists(bc):
        m = dict(re.findall(r'(\w+)=(\S+)', open(bc).read()))
        r['bam_records'] = int(m['both']); r['bam_nonidentical'] = int(m['both']) - int(m['identical'])
        r['bam_only_one_arm'] = int(m['cpu_only']) + int(m['gpu_only'])
    for a, t in (('cpu', tc), ('gpu', tg)):
        tlong += [{'sample': S, 'read_pairs': r['read_pairs'], 'arm': a, 'step': k, 'seconds': v} for k, v in t.items()]
    X0, Y0 = load(f'{p}_cpu_nods.calls.tsv.gz'), load(f'{p}_gpu_nods.calls.tsv.gz')
    for v in X0.values(): dp_all[dp_bin(int(v[3]))] += 1
    for a, A, B in (('cpu', X0, Y0), ('gpu', Y0, X0)):
        for k in set(A) - set(B):
            v = A[k]
            disc.append({'sample': S, 'arm': a, 'chrom': k[0], 'pos': k[1], 'ref': k[2], 'alt': k[3], 'filter': v[0],
                         'DP': v[3], 'alt_reads': v[4], 'AF': v[5], 'on_target': on_target(k[0], k[1])})
    for k in sorted(set(X0) | set(Y0)):
        if not on_target(k[0], k[1]): continue
        x, y = X0.get(k, ('',) * 6), Y0.get(k, ('',) * 6)
        tcalls.append({'sample': S, 'chrom': k[0], 'pos': k[1], 'ref': k[2], 'alt': k[3],
                       'cpu_filter': x[0], 'gpu_filter': y[0], 'cpu_DP': x[3], 'gpu_DP': y[3], 'cpu_AF': x[5], 'gpu_AF': y[5]})
    for scope in ('genome', 'target'):
        X = X0 if scope == 'genome' else {k: v for k, v in X0.items() if on_target(k[0], k[1])}
        Y = Y0 if scope == 'genome' else {k: v for k, v in Y0.items() if on_target(k[0], k[1])}
        px = {k for k, v in X.items() if v[0] == 'PASS'}; py = {k for k, v in Y.items() if v[0] == 'PASS'}
        r[f'{scope}_raw_n'] = len(set(X) | set(Y)); r[f'{scope}_raw_J'] = round(J(set(X), set(Y)), 5)
        r[f'{scope}_raw_disc'] = len(set(X) ^ set(Y))
        r[f'{scope}_pass_n'] = len(px | py); r[f'{scope}_pass_J'] = round(J(px, py), 5); r[f'{scope}_pass_disc'] = len(px ^ py)
    rows.append(r)

def write(name, recs):
    cols = list(dict.fromkeys(k for r in recs for k in r))
    with open(f'{ROOT}/analysis/tables/{name}', 'w', newline='') as f:
        w = csv.DictWriter(f, cols, delimiter='\t', restval=''); w.writeheader(); w.writerows(recs)
write('bench_summary.tsv', rows); write('bench_timing_long.tsv', tlong)
write('bench_target_calls.tsv', tcalls); write('bench_discordant.tsv', disc)
write('bench_dp_all.tsv', [{'DP_bin': b, 'n': dp_all[b]} for b in DP_BINS[:-1]])
show = ['sample', 'read_pairs', 'bam_speedup', 'm2_speedup', 'total_speedup', 'flagstat_identical', 'bam_records', 'bam_nonidentical', 'bam_only_one_arm',
        'genome_raw_J', 'genome_raw_disc', 'genome_pass_disc', 'target_raw_J', 'target_pass_J', 'target_pass_disc']
print('\t'.join(show))
for r in rows: print('\t'.join(str(r.get(k, '')) for k in show))
print(f'{len(rows)} cases')
