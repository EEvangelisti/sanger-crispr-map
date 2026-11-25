#!/usr/bin/env python3
"""Reference-guided visualisation of CRISPR edits from Sanger reads.

Each Sanger read (an ``.ab1`` chromatogram or a FASTA-formatted ``.seq`` file)
is aligned locally to a reference gene. Deletions and insertions are called
from the alignment and drawn, one row per sample, under an annotated map of
the reference (CDS, primers, CRISPR targets with their PAM and crRNA).

Inputs
------
* A template FASTA file whose record IDs end with a type suffix:

  - ``<id>|gene``    exactly one record: the reference sequence;
  - ``<id>|cds``     optional; the coding sequence (exact match or, failing
                     that, inferred exon blocks by alignment);
  - ``<name>|primer``  optional; located by exact match on either strand;
  - ``<name>|crispr``  optional; a CRISPR target written as 5' PAM + protospacer
                     (or the protospacer alone when ``--pam-motif`` is given).

* Read files, given as files or folders (searched recursively):
  ``.seq`` (FASTA) by default, ``.ab1`` with ``--ab1``.

Read file names
---------------
Reads are grouped into samples from their file stem:

* ``1_17_1_623`` -> sample ``1_17`` (``--group-regex``, last two fields are the
  replicate / primer);
* ``WT_618``, ``WT_2_623``, ``wild-type_618`` -> sample ``WT``.

With ``--split-bands`` (gel bands sequenced separately), the third field is the
band: ``2_75_1_623`` -> sample ``2_75``, band 1; ``2_75_2_623`` -> band 2.
Numbered names without a band field (``2_75_623``) go to band 1, WT reads go to
both bands, and any other name is excluded.

Outputs (in ``--out``; ``_1`` / ``_2`` suffixes with ``--split-bands``)
---------------------------------------------------------------------
* ``edit_map*.pdf|png``     the figure (plus ``edit_map_combined`` for bands);
* ``summary*.tsv``          one row per sample;
* ``events*.tsv``           one row per distinct indel per sample;
* ``report*.txt``           human-readable report, per read and per event;
* ``qc_reads.tsv|pdf|png``  per-read quality control (``--ab1`` only).

All coordinates are 0-based, end-exclusive internally and 1-based, inclusive
in every output file.
"""
import argparse
import csv
import re
import struct
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib
matplotlib.use('Agg')  # headless backend: figures are only written to files
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrow, Patch, Rectangle


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Quality-control thresholds for the 'status' column of qc_reads.tsv.
QC_MIN_USABLE_BP = 100      # shorter usable window -> FAILED
QC_MIN_MEAN_Q = 30          # lower mean Phred over the usable window -> LOW QUALITY
QC_MAX_MIXED_PCT = 10.0     # more double peaks (% of assessable bases) -> MIXED TRACE

# Rescue of short aligned blocks at the border of a large deletion
# (see summarize_alignment).
EDGE_GAP = 50       # reference gap (bp) that counts as a large-deletion border
EDGE_MIN_INFO = 8   # informative columns a short block needs to be rescued
EDGE_NEAR = 10      # max gap (bp, reference and read) to its reliable neighbour
EDGE_MIN_IDENTITY = 0.9

# Read-to-reference alignment scores (gap scores are command-line options).
ALIGN_MATCH = 2
ALIGN_MISMATCH = -3

# An event within this distance (bp) of a CRISPR target is reported as nearby.
NEARBY_TARGET_BP = 50

# Figure colours.
GREY = '#92979f'          # aligned region / reference backbone
GREY_LIGHT = '#ced0d4'    # legend swatch for the white veil over grey
VIOLET = '#8652b4'        # CRISPR target region (--highlight)
PURPLE_DARK = '#653a91'   # crRNA arrows and target names
CDS_GREY = '#545961'
PAM_GOLD = '#d4a34f'
INK = '#303030'

IUPAC = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'R': 'AG', 'Y': 'CT', 'S': 'CG', 'W': 'AT', 'K': 'GT', 'M': 'AC',
    'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT',
}


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def records(path):
    """Return [(id, sequence)] for a FASTA file, upper-cased, with U read as T."""
    return [(r.id, str(r.seq).upper().replace('U', 'T')) for r in SeqIO.parse(str(path), 'fasta')]


def revcomp(seq):
    """Reverse complement of a DNA string."""
    return str(Seq(seq).reverse_complement())


def merge(intervals):
    """Merge overlapping or touching half-open intervals; returns them sorted."""
    out = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(e, out[-1][1]))
        else:
            out.append((s, e))
    return out


def subtract_intervals(spans, remove):
    """Parts of the merged, sorted intervals ``spans`` not covered by ``remove``
    (also merged and sorted)."""
    out = []
    for a, b in spans:
        cur = a
        for ra, rb in remove:
            if rb <= cur or ra >= b:
                continue
            if ra > cur:
                out.append((cur, ra))
            cur = max(cur, rb)
        if cur < b:
            out.append((cur, b))
    return out


def informative(base):
    """True for an unambiguous base (A, C, G or T)."""
    return base in 'ACGT'


# ---------------------------------------------------------------------------
# .ab1 (ABIF) reading
#
# A small self-contained parser is used instead of Bio.SeqIO's 'abi' reader so
# that the result does not depend on the Biopython version, and so that peak
# positions and the processed trace channels are available for the
# double-peak estimate.
# ---------------------------------------------------------------------------

def read_abif(path):
    """Parse the directory of an ABIF file.

    Returns ``{(tag_name, tag_number): (element_type, n_elements, raw_bytes)}``.
    Raises ValueError if the file does not start with the ABIF signature.
    """
    data = Path(path).read_bytes()
    if data[:4] != b'ABIF':
        raise ValueError(f'{path}: not an ABIF (.ab1) file')
    # Directory entry: name, number, element type, element size, number of
    # elements, data size, data offset, data handle (28 bytes, big-endian).
    entry = struct.Struct('>4sihhiiii')
    # The header (after the 4-byte signature and 2-byte version) is itself a
    # directory entry pointing to the tag directory.
    _, _, _, _, n_entries, _, dir_offset, _ = entry.unpack_from(data, 6)
    tags = {}
    for k in range(n_entries):
        o = dir_offset + k * entry.size
        name, num, etype, _, n_elem, size, offset, _ = entry.unpack_from(data, o)
        # Data of 4 bytes or less is stored in the offset field itself (byte 20).
        raw = data[o + 20:o + 20 + size] if size <= 4 else data[offset:offset + size]
        tags[(name.decode('ascii', 'replace'), num)] = (etype, n_elem, raw)
    return tags


def ab1_record(path):
    """Base calls, Phred qualities, peak positions and traces of one .ab1 file.

    Uses the edited base calls (PBAS1/PCON1/PLOC1) when present, otherwise the
    basecaller output (PBAS2/PCON2/PLOC2). The traces are the processed
    channels DATA9-12, keyed by base following the FWO_ channel order.

    Returns a dict with keys ``name`` (file stem), ``bases`` (str), ``quals``
    (list of int), ``ploc`` (tuple of trace indices, or None) and ``traces``
    ({base: tuple of int}, possibly empty).
    """
    tags = read_abif(path)

    def raw(name, num):
        return tags[(name, num)][2] if (name, num) in tags else None

    for num in (1, 2):
        if raw('PBAS', num) is not None and raw('PCON', num) is not None:
            break
    else:
        raise ValueError(f'{path}: no base calls (PBAS/PCON) found')
    bases = raw('PBAS', num).decode('ascii', 'replace').upper()
    quals = list(raw('PCON', num))
    n = min(len(bases), len(quals))
    bases, quals = bases[:n], quals[:n]

    ploc = None
    if raw('PLOC', num) is not None:
        ploc_raw = raw('PLOC', num)
        ploc = struct.unpack('>%dh' % (len(ploc_raw) // 2), ploc_raw)[:n]

    traces = {}
    order = raw('FWO_', 1)
    if order:
        for i, base in enumerate(order.decode('ascii', 'replace').upper()):
            d = raw('DATA', 9 + i)
            if d is not None:
                traces[base] = struct.unpack('>%dh' % (len(d) // 2), d)

    return dict(name=Path(path).stem, bases=bases, quals=quals, ploc=ploc, traces=traces)


def mott_window(quals, error_cutoff):
    """Modified Mott trimming, as in Phred and Biopython.

    Finds the contiguous window maximising the sum of
    ``error_cutoff - 10**(-Q/10)``. Returns ``(start, end)``, end exclusive;
    ``(0, 0)`` if no base passes.
    """
    best, cur, start = 0.0, 0.0, 0
    best_start = best_end = 0
    for i, q in enumerate(quals):
        cur += error_cutoff - 10 ** (-q / 10)
        if cur <= 0:
            cur, start = 0.0, i + 1
        elif cur > best:
            best, best_start, best_end = cur, start, i + 1
    return best_start, best_end


def window_trim(quals, min_quality, window=10):
    """Trim the two ends of a read on a sliding-window mean Phred.

    The kept region starts at the first window whose mean reaches
    ``min_quality`` and ends at the last such window; it is then shrunk so that
    its first and last bases themselves reach ``min_quality``. Unlike Mott
    trimming this keeps ONE long region and tolerates isolated low-quality
    bases inside it (they are masked individually afterwards instead of
    cutting the read).

    Returns ``(start, end)``, end exclusive; ``(0, 0)`` if nothing passes.
    """
    n = len(quals)
    if n == 0:
        return 0, 0
    w = max(1, min(window, n))
    csum = [0]
    for q in quals:
        csum.append(csum[-1] + q)
    start = next((i for i in range(n - w + 1) if (csum[i + w] - csum[i]) / w >= min_quality), None)
    if start is None:
        return 0, 0
    # A passing window exists, so this search always succeeds.
    end = next(j for j in range(n, w - 1, -1) if (csum[j] - csum[j - w]) / w >= min_quality)
    while start < end and quals[start] < min_quality:
        start += 1
    while end > start and quals[end - 1] < min_quality:
        end -= 1
    return (start, end) if end > start else (0, 0)


def trim_bounds(quals, mode, min_quality, window=10):
    """Region ``(start, end)`` of a read kept before per-base masking.

    mode ``'window'``: sliding-window trimming of the ends (window_trim);
    mode ``'mott'``  : modified Mott algorithm (mott_window), a single best
                       window, harsh on reads with scattered low-quality bases;
    mode ``'none'``  : no trimming.
    ``min_quality <= 0`` also disables trimming.
    """
    n = len(quals)
    if min_quality <= 0 or mode == 'none':
        return 0, n
    if mode == 'mott':
        return mott_window(quals, 10 ** (-min_quality / 10))
    return window_trim(quals, min_quality, window)


def secondary_peak_ratios(bases, ploc, traces):
    """Per base: height of the strongest competing trace channel divided by
    the height of the called channel, at the called peak position.

    The value is None where it cannot be assessed (no peak positions, called
    base not A/C/G/T or without a trace, position out of range, or called
    channel <= 0). Returns an empty list when ``ploc`` is None.
    """
    ratios = []
    for b, p in zip(bases, ploc or ()):
        if b not in traces or not 0 <= p < len(traces[b]):
            ratios.append(None)
            continue
        main = traces[b][p]
        if main <= 0:
            ratios.append(None)
            continue
        other = max(traces[x][p] for x in traces if x != b)
        ratios.append(max(other, 0) / main)
    return ratios


def load_ab1_read(path, min_quality=20, mixed_ratio=0.3, trim_mode='window', trim_window=10):
    """Read an .ab1 file, trim its low-quality ends and mask low-quality bases.

    1. The ends are trimmed according to ``trim_mode`` (see trim_bounds).
    2. Within the kept region, bases with Phred < ``min_quality`` become N, so
       they align to anything at no cost and never count as covered.
    3. Quality metrics are computed over the kept region, including the
       percentage of double peaks (competing channel >= ``mixed_ratio`` of the
       called one), and a QC status is assigned (OK, FAILED, LOW QUALITY or
       MIXED TRACE; see the QC_* constants).

    ``min_quality = 0`` disables both trimming and masking.

    Returns ``(record name, sequence for alignment, qc dict)``. The qc dict
    holds the QC_FIELDS columns (except ``sample``, ``strand`` and
    ``aligned_bp``, filled in later) plus the private keys ``_quals`` and
    ``_window`` used by the QC figure.
    """
    rec = ab1_record(path)
    bases, quals = rec['bases'], rec['quals']
    n = len(bases)
    s, e = trim_bounds(quals, trim_mode, min_quality, trim_window)
    window_bases, window_quals = bases[s:e], quals[s:e]
    seq = ''.join(b if q >= min_quality else 'N' for b, q in zip(window_bases, window_quals))
    masked = sum(1 for b, q in zip(window_bases, window_quals) if q < min_quality and b != 'N')

    ratios = secondary_peak_ratios(bases, rec['ploc'], rec['traces'])[s:e]
    assessable = [r for b, r in zip(window_bases, ratios) if informative(b) and r is not None]
    mixed = sum(1 for r in assessable if r >= mixed_ratio)
    mixed_pct = 100 * mixed / len(assessable) if assessable else float('nan')

    usable = e - s
    mean_raw = sum(quals) / n if n else 0.0
    mean_usable = sum(window_quals) / usable if usable else 0.0
    if usable < QC_MIN_USABLE_BP:
        status = f'FAILED (usable window < {QC_MIN_USABLE_BP} bp)'
    elif mean_usable < QC_MIN_MEAN_Q:
        status = f'LOW QUALITY (mean Q < {QC_MIN_MEAN_Q})'
    elif mixed_pct > QC_MAX_MIXED_PCT:   # False when mixed_pct is NaN
        status = f'MIXED TRACE (> {QC_MAX_MIXED_PCT:g}% double peaks)'
    else:
        status = 'OK'
    qc = dict(
        file=str(path),
        raw_length=n,
        mean_q_raw=round(mean_raw, 1),
        trim_start_1based=s + 1 if usable else 0,
        trim_end=e,
        trimmed_5prime=s,
        trimmed_3prime=n - e,
        usable_length=usable,
        mean_q_usable=round(mean_usable, 1),
        pct_q30_usable=round(100 * sum(q >= 30 for q in window_quals) / usable, 1) if usable else 0.0,
        n_calls_in_window=window_bases.count('N'),
        masked_low_q=masked,
        double_peaks_pct=round(mixed_pct, 1) if assessable else '',
        status=status,
        _quals=quals,
        _window=(s, e),
    )
    return rec['name'], seq, qc


# ---------------------------------------------------------------------------
# Read quality control output (--ab1)
# ---------------------------------------------------------------------------

QC_FIELDS = ['file', 'sample', 'raw_length', 'mean_q_raw', 'trim_start_1based', 'trim_end',
             'trimmed_5prime', 'trimmed_3prime', 'usable_length', 'mean_q_usable',
             'pct_q30_usable', 'n_calls_in_window', 'masked_low_q', 'double_peaks_pct',
             'strand', 'aligned_bp', 'status']


def plot_read_quality(qc_rows, min_quality, out_prefix, font_size=10):
    """Write ``<out_prefix>.pdf|png``: one horizontal strip per read, each base
    coloured by Phred class (>= 30, between the threshold and 29, below the
    threshold), trimmed ends in light grey."""
    rows = sorted(qc_rows, key=lambda r: Path(r['file']).name)
    if not rows:
        return
    palette = {'trim': '#d9d9d9', 'high': '#4c78a8', 'mid': '#f2a541', 'low': '#c0392b'}
    low_cut = max(min_quality, 1)
    max_len = max(r['raw_length'] for r in rows)
    fig, ax = plt.subplots(figsize=(10, max(2.4, 1.3 + 0.24 * len(rows))))
    for y, r in enumerate(rows):
        s, e = r['_window']
        classes = []
        for i, q in enumerate(r['_quals']):
            if i < s or i >= e:
                classes.append('trim')
            elif q >= 30:
                classes.append('high')
            elif q >= low_cut:
                classes.append('mid')
            else:
                classes.append('low')
        # Run-length encode the classes into (start, length) segments.
        runs = defaultdict(list)
        prev, start = None, 0
        for i, c in enumerate(classes + [None]):
            if c != prev:
                if prev is not None:
                    runs[prev].append((start, i - start))
                prev, start = c, i
        for c, segs in runs.items():
            ax.broken_barh(segs, (y - 0.38, 0.76), facecolors=palette[c], linewidth=0)
    ax.set_yticks(range(len(rows)), [Path(r['file']).stem for r in rows])
    ax.set_ylim(len(rows) - 0.4, -0.6)
    ax.set_xlim(0, max_len)
    ax.set_xlabel('Position in read (bases)', fontsize=font_size)
    ax.set_title('Read quality (Phred) from .ab1 files', fontsize=font_size + 2)
    ax.tick_params(axis='y', labelsize=font_size * 0.8, length=0)
    ax.tick_params(axis='x', labelsize=font_size)
    handles = [
        Patch(facecolor=palette['high'], label='Q ≥ 30'),
        Patch(facecolor=palette['mid'], label=f'Q {low_cut}-29'),
        Patch(facecolor=palette['low'], label=f'Q < {low_cut} (masked as N)'),
        Patch(facecolor=palette['trim'], label='Trimmed ends'),
    ]
    fig.legend(handles=handles, loc='lower center', ncol=4, frameon=False,
               fontsize=font_size * 0.85, bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0, 0.5 / fig.get_figheight(), 1, 1))
    fig.savefig(str(out_prefix) + '.pdf', bbox_inches='tight')
    fig.savefig(str(out_prefix) + '.png', dpi=200, bbox_inches='tight')
    plt.close(fig)


def write_qc(args, qc_rows):
    """Write qc_reads.tsv and the quality-strip figure, and print a summary
    listing every read whose status is not OK."""
    if not qc_rows:
        return
    rows = sorted(qc_rows, key=lambda r: Path(r['file']).name)
    with (args.out / 'qc_reads.tsv').open('w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=QC_FIELDS, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        for r in rows:
            writer.writerow(dict(r, file=Path(r['file']).name))
    plot_read_quality(rows, args.min_quality, args.out / 'qc_reads', args.font_size * 0.8)
    bad = [r for r in rows if r['status'] != 'OK']
    print(f'QC: {len(rows)} .ab1 reads, {len(rows) - len(bad)} OK, {len(bad)} to review (see qc_reads.tsv)')
    for r in bad:
        print(f"  {Path(r['file']).name}: {r['status']}; usable {r['usable_length']} bp, mean Q {r['mean_q_usable']}")


def load_reads(path, args):
    """Return [(record name, sequence, qc dict or None)] for one read file."""
    if args.ab1:
        name, seq, qc = load_ab1_read(path, args.min_quality, args.mixed_ratio,
                                      args.trim_mode, args.trim_window)
        return [(name, seq, qc)]
    return [(name, seq, None) for name, seq in records(path)]


def read_diagnostics(seq, window=50):
    """Print (--verbose) the length and the ambiguous-base profile of a read,
    in windows of ``window`` bases (one '#' per 5 %)."""
    n = len(seq)
    ambiguous = sum(not informative(b) for b in seq)
    print(f'  Raw length: {n} bp')
    print(f'  Ambiguous bases: {ambiguous}')
    print(f'  Ambiguity: {100 * ambiguous / max(n, 1):.2f}%')
    for start in range(0, n, window):
        chunk = seq[start:start + window]
        count = sum(not informative(b) for b in chunk)
        pct = 100 * count / len(chunk)
        bar = '#' * round(pct / 5)
        print(f'  {start + 1:5d}-{start + len(chunk):5d} | {pct:6.1f}% | {bar}')


# ---------------------------------------------------------------------------
# Alignment and indel calling
# ---------------------------------------------------------------------------

def reliable_intervals(covered, span, ref_len, window=20, min_frac=0.7):
    """Parts of the aligned span that are well supported by informative bases.

    A position is reliable when, in a sliding window of ``window`` bp centred
    on it, at least ``min_frac`` of the positions that lie inside the aligned
    span are covered by an informative (not masked, not ambiguous) base.
    Normalising by the in-span part of the window keeps the read edges from
    being penalised for the empty space beyond them.

    The result is only used to lighten weakly supported parts of the figure;
    it never removes coverage or events.
    """
    flags = np.zeros(ref_len)
    for a, b in covered:
        flags[a:b] = 1.0
    inside = np.zeros(ref_len)
    for a, b in span:
        inside[a:b] = 1.0
    kernel = np.ones(max(1, int(window)))
    num = np.convolve(flags, kernel, mode='same')
    den = np.convolve(inside, kernel, mode='same')
    density = np.divide(num, den, out=np.zeros_like(num), where=den > 0)
    ok = density >= min_frac
    out = []
    for a, b in span:
        start = None
        for i in range(a, b):
            if ok[i] and start is None:
                start = i
            elif not ok[i] and start is not None:
                out.append((start, i))
                start = None
        if start is not None:
            out.append((start, b))
    return out


def _aligned_blocks(coords, ref, read):
    """Gap-free blocks of a Biopython alignment.

    Returns [(ref_start, ref_end, read_start, read_end, n_informative)], where
    n_informative counts columns with an unambiguous base on both sequences.
    """
    blocks = []
    for i in range(coords.shape[1] - 1):
        r0, r1 = int(coords[0, i]), int(coords[0, i + 1])
        q0, q1 = int(coords[1, i]), int(coords[1, i + 1])
        if r1 > r0 and q1 > q0:
            n = min(r1 - r0, q1 - q0)
            info = sum(1 for j in range(n) if informative(read[q0 + j]) and informative(ref[r0 + j]))
            blocks.append((r0, r0 + n, q0, q0 + n, info))
    return blocks


def _is_rescued_edge_block(i, blocks, solid, ref, read):
    """True if the short block ``blocks[i]`` is the genuine end of a deletion
    flank and must be kept even though it is below --min-flank.

    A 2-bp slip or a homopolymer error a few bases from a deletion junction can
    cut the last ~12 bp of the flank off from the main block; dropping them
    would leave an uncovered hole beside the deletion. The block is rescued
    only if ALL of these hold, so that isolated random matches in noisy
    stretches never create a deletion:

    * it has at least EDGE_MIN_INFO informative columns;
    * at least EDGE_MIN_IDENTITY of them match the reference;
    * on one side it is joined to a reliable block (gap <= EDGE_NEAR on both
      the reference and the read), and on the other side it borders a large
      reference gap (>= EDGE_GAP).
    """
    r0, r1, q0, _, info = blocks[i]
    if info < EDGE_MIN_INFO:
        return False
    matches = sum(1 for j in range(r1 - r0)
                  if informative(read[q0 + j]) and informative(ref[r0 + j]) and read[q0 + j] == ref[r0 + j])
    if matches / info < EDGE_MIN_IDENTITY:
        return False

    def ordered(n):
        return (blocks[n], blocks[i]) if n < i else (blocks[i], blocks[n])

    def joined(n):
        if not (0 <= n < len(blocks)) or not solid[n]:
            return False
        lo, hi = ordered(n)
        return hi[0] - lo[1] <= EDGE_NEAR and abs(hi[2] - lo[3]) <= EDGE_NEAR

    def big_gap(n):
        if not (0 <= n < len(blocks)):
            return False
        lo, hi = ordered(n)
        return hi[0] - lo[1] >= EDGE_GAP

    return (joined(i - 1) and big_gap(i + 1)) or (joined(i + 1) and big_gap(i - 1))


def summarize_alignment(coords, ref, read, min_flank=0, conf_window=20, conf_min_frac=0.7):
    """Coverage, indels and mismatches from Biopython alignment coordinates.

    Bases that are not A/C/G/T (masked or ambiguous) never count as covered.

    With ``min_flank > 0``, an aligned block with fewer than ``min_flank``
    informative columns is considered unreliable and dropped (unless rescued
    as a deletion edge, see _is_rescued_edge_block). The gaps around it then
    merge into ONE net event between the neighbouring kept blocks; when the
    same number of reference and read bases is skipped there is no indel at
    all. This removes the spurious deletion + insertion pairs that appear
    where a noisy stretch of a read is bridged by two gaps instead of being
    aligned through.

    Between two consecutive kept blocks, ``skipped_ref - skipped_read`` > 0 is
    a deletion starting right after the first block, and < 0 an insertion of
    that many bases placed right after the first block.

    Returns ``(covered, deletions, insertions, mismatches, n_informative,
    span, reliable)``:

    * covered    merged [(start, end)] of informative reference positions;
    * deletions  [(start, end)] reference intervals;
    * insertions [(position, size)], inserted before reference index position;
    * mismatches [position];
    * n_informative  informative columns in the kept blocks;
    * span       merged [(start, end)] of the kept blocks (masked bases included);
    * reliable   sub-intervals of span (see reliable_intervals).
    """
    blocks = _aligned_blocks(coords, ref, read)
    if min_flank > 0:
        solid = [b[4] >= min_flank for b in blocks]
        kept = [b for i, b in enumerate(blocks)
                if solid[i] or _is_rescued_edge_block(i, blocks, solid, ref, read)]
    else:
        kept = blocks

    covered, mismatches, deletions, insertions = [], [], [], []
    for r0, r1, q0, _, _ in kept:
        for j in range(r1 - r0):
            if informative(read[q0 + j]) and informative(ref[r0 + j]):
                covered.append((r0 + j, r0 + j + 1))
                if read[q0 + j] != ref[r0 + j]:
                    mismatches.append(r0 + j)
    for prev, nxt in zip(kept, kept[1:]):
        skipped_ref = nxt[0] - prev[1]
        skipped_read = nxt[2] - prev[3]
        if skipped_ref > skipped_read:
            deletions.append((prev[1], prev[1] + skipped_ref - skipped_read))
        elif skipped_read > skipped_ref:
            insertions.append((prev[1], skipped_read - skipped_ref))
    covered = merge(covered)
    span = merge([(b[0], b[1]) for b in kept])
    reliable = reliable_intervals(covered, span, len(ref), conf_window, conf_min_frac)
    return covered, deletions, insertions, mismatches, sum(b[4] for b in kept), span, reliable


def align(ref, seq, gap_open=-2.0, gap_extend=-0.05, min_flank=0, conf_window=20, conf_min_frac=0.7):
    """Align one read to the reference and call its events.

    Local affine-gap alignment (match +2, mismatch -3) of the read and of its
    reverse complement; the higher-scoring strand is kept. Terminal Ns are
    removed first; internal Ns are kept and score 0 against any base, so a
    masked stretch is aligned through rather than bridged by gaps.

    Returns a dict with keys ``score``, ``strand`` ('+', '-' or '?' when
    nothing aligned), ``coverage``, ``span``, ``reliable``, ``deletions``,
    ``insertions``, ``mismatches``, ``informative``, ``aligned_bases``
    (number of covered reference positions) and ``read_length``.
    """
    seq = seq.upper().strip('N')
    aligner = PairwiseAligner(mode='local')
    aligner.match_score = ALIGN_MATCH
    aligner.mismatch_score = ALIGN_MISMATCH
    aligner.open_gap_score = gap_open
    aligner.extend_gap_score = gap_extend
    try:
        aligner.wildcard = 'N'
    except (AttributeError, TypeError, ValueError):
        pass  # very old Biopython: N scores as a mismatch (--min-flank still helps)
    candidates = []
    for strand, s in [('+', seq), ('-', revcomp(seq))]:
        if not s or not any(informative(b) for b in s):
            continue
        try:
            al = next(iter(aligner.align(ref, s)), None)
        except IndexError:
            al = None
        if al is None:
            continue
        candidates.append((al.score, strand, s, al))
    if not candidates:
        return dict(score=0, strand='?', coverage=[], span=[], reliable=[], deletions=[],
                    insertions=[], mismatches=[], informative=0, aligned_bases=0,
                    read_length=len(seq))
    score, strand, s, al = max(candidates, key=lambda x: x[0])
    covered, deletions, insertions, mismatches, n_info, span, reliable = summarize_alignment(
        al.coordinates, ref, s, min_flank, conf_window, conf_min_frac)
    return dict(score=score, strand=strand, coverage=covered, span=span, reliable=reliable,
                deletions=deletions, insertions=insertions, mismatches=mismatches,
                informative=n_info, aligned_bases=sum(e - b for b, e in covered),
                read_length=len(seq))


def passing_hits(hits, min_aligned):
    """Alignments with at least ``min_aligned`` covered reference positions."""
    return [hit for hit in hits if hit['aligned_bases'] >= min_aligned]


def collect_events(hits):
    """Distinct indels of a sample and the read files supporting each.

    Returns ``{(kind, a, b): {file, ...}}`` in first-seen order, where
    ``(kind, a, b)`` is ``('deletion', start, end)`` or
    ``('insertion', position, size)``.
    """
    events = defaultdict(set)
    for hit in hits:
        for a, b in hit['deletions']:
            events['deletion', a, b].add(hit['file'])
        for pos, size in hit['insertions']:
            events['insertion', pos, size].add(hit['file'])
    return events


def compute_wt_events(groups, min_aligned):
    """Set of indel keys (see collect_events) found in passing WT reads.

    The same event in an edited sample is flagged with '*' and 'also_in_wt',
    as it is more likely an alignment artefact than an edit.
    """
    wt_events = set()
    for sample, hits in groups.items():
        if sample.upper() == 'WT':
            wt_events.update(collect_events(passing_hits(hits, min_aligned)))
    return wt_events


# ---------------------------------------------------------------------------
# Sample naming
# ---------------------------------------------------------------------------

WT_STEM = r'(?i)(WT|WILD[_-]?TYPE)'


def sample_name(path, pattern):
    """Sample ID of a read file, from its stem.

    ``pattern`` (--group-regex) is tried first: its named group ``sample``, or
    else its first group, is the ID. Otherwise WT stems (``WT_618``,
    ``WT_2_623``) give 'WT', ``<a>_<b>_<c>`` gives ``<a>_<b>``, and anything
    else is returned unchanged.
    """
    stem = path.stem
    m = re.fullmatch(pattern, stem)
    if m:
        return m.group('sample') if 'sample' in m.groupdict() else m.group(1)
    if re.fullmatch(WT_STEM + r'(?:_\d+)?_\d+', stem):
        return 'WT'
    m = re.fullmatch(r'(\d+_\d+)_(\d+)', stem)
    if m:
        return m.group(1)
    return stem


def band_name(path, pattern):
    """``(sample ID, band)`` of a read file for --split-bands.

    * WT stems -> ('WT', None): WT reads are used in every band;
    * ``<a>_<b>_<1|2>[_<primer>]`` -> ('<a>_<b>', '1' or '2');
    * ``<a>_<b>[_<n>]`` -> ('<a>_<b>', '1'): no band field means band 1;
    * anything else -> (sample_name(path, pattern), None): excluded by the caller.
    """
    stem = path.stem
    if re.fullmatch(WT_STEM + r'(?:_\d+){0,2}', stem):
        return ('WT', None)
    m = re.fullmatch(r'(\d+_\d+)_([12])(?:_\d+)?', stem)
    if m:
        return (m.group(1), m.group(2))
    m = re.fullmatch(r'(\d+_\d+)(?:_\d+)?', stem)
    if m:
        return (m.group(1), '1')
    return (sample_name(path, pattern), None)


def sorted_samples(samples):
    """Sample IDs in display order: alphabetical, WT last."""
    return sorted(samples, key=lambda s: (s.upper() == 'WT', s))


def select_plotted_samples(args, groups, band):
    """Sample rows to draw for one band (``band`` is '' without --split-bands).

    Without --samples every sample of the band is drawn. With --samples and
    --split-bands, the list is used as given, so samples absent from a band
    still get an empty 'No aligned reads' row and the rows line up between
    bands. --samples-N then restricts (and orders) band N, and
    --hide-samples-N removes samples from it.
    """
    if args.samples is None:
        plotted = sorted_samples(groups)
    elif args.split_bands:
        plotted = list(args.samples)
    else:
        plotted = [s for s in args.samples if s in groups]

    if band:
        selected = getattr(args, f'samples_{band}')
        hidden = set(getattr(args, f'hide_samples_{band}'))
        if selected is not None:
            unknown = [s for s in selected if s not in plotted]
            if unknown:
                print(f'WARNING: band {band}: samples outside global selection: {unknown}')
            plotted = [s for s in selected if s in plotted]
        plotted = [s for s in plotted if s not in hidden]
    return plotted


# ---------------------------------------------------------------------------
# Reference annotation: CDS, primers, CRISPR targets and PAM
# ---------------------------------------------------------------------------

@dataclass
class Reference:
    """The reference gene and its annotations (0-based, end-exclusive)."""
    gene_id: str
    seq: str
    cds_length: int = None                          # None if no |cds record
    cds_note: str = ''
    cds_blocks: list = field(default_factory=list)  # [(start, end)]
    cds_strand: str = None                          # '+', '-' or None if unknown
    primers: list = field(default_factory=list)     # [(name, [(start, end, strand)])]
    guides: list = field(default_factory=list)      # [(name, [(start, end, strand)])]
    pam_sites: list = field(default_factory=list)   # [(start, end, strand)]

    @property
    def name(self):
        return self.gene_id.split('|')[0]


def pam_interval(a, b, strand, pam_length=4):
    """Reference interval of the 5' PAM of a target site (a, b, strand).

    The site spans PAM + protospacer, so on '+' the PAM is its first
    ``pam_length`` bases and on '-' its last. None if the site is too short.
    """
    if b - a <= pam_length:
        return None
    return (a, a + pam_length) if strand == '+' else (b - pam_length, b)


def protospacer_interval(a, b, strand, pam_length=4):
    """Reference interval of the protospacer of a target site (PAM excluded)."""
    if b - a <= pam_length:
        return None
    return (a + pam_length, b) if strand == '+' else (a, b - pam_length)


def crrna_display_interval(a, b, strand, pam_length=4, min_display_bp=30):
    """Interval over which the crRNA arrow is drawn: the protospacer, widened
    symmetrically to ``min_display_bp`` so that it stays visible at gene scale."""
    proto = protospacer_interval(a, b, strand, pam_length)
    if proto is None:
        return None
    ps, pe = proto
    if pe - ps >= min_display_bp:
        return (ps, pe)
    center = (ps + pe) / 2
    half = min_display_bp / 2
    return (center - half, center + half)


def crrna_strand(site_strand):
    """Strand the crRNA hybridises to: the one opposite the PAM-carrying strand."""
    return '-' if site_strand == '+' else '+'


def pam_matches(seq, motif):
    """True if ``seq`` matches the IUPAC ``motif`` (same length, any case).
    Always True when ``motif`` is None."""
    if motif is None:
        return True
    if len(seq) != len(motif):
        return False
    return all(b in IUPAC[m] for b, m in zip(seq.upper(), motif.upper()))


def pam_regex(motif):
    """Regex for an IUPAC motif; a lookahead, so overlapping hits are found."""
    return '(?=(' + ''.join('[' + IUPAC[m] + ']' for m in motif) + '))'


def find_pam_sites(ref, motif):
    """All occurrences of a PAM motif in the reference, on both strands.

    Returns sorted [(start, end, strand)] in forward-reference coordinates;
    on '-' the motif reads 5'->3' on the reverse strand. Only used for the
    count in the report.
    """
    sites = [(m.start(), m.start() + len(motif), '+') for m in re.finditer(pam_regex(motif), ref)]
    sites += [(m.start(), m.start() + len(motif), '-') for m in re.finditer(pam_regex(revcomp(motif)), ref)]
    return sorted(sites)


def pam_sequence(ref, a, b, strand, pam_length=4):
    """PAM sequence of a target site, read 5'->3' on the PAM-carrying strand."""
    iv = pam_interval(a, b, strand, pam_length)
    if iv is None:
        return None
    seq = ref[iv[0]:iv[1]]
    return seq if strand == '+' else revcomp(seq)


def locate_guide(ref, seq, pam_length, pam_motif=None):
    """Locate a CRISPR target with a 5' PAM (Cas12a-type) in the reference.

    Without ``pam_motif``: exact search of the whole template sequence
    (PAM + protospacer) on both strands.

    With ``pam_motif``: the PAM written in the template is not trusted. The
    protospacer is searched on both strands, the ``pam_length`` reference bases
    immediately 5' of it (on the protospacer's strand) are read and tested
    against the motif, and only matching sites are kept. Two readings of the
    template are tried, the second only if the first gives no valid site:

    1. the template is PAM + protospacer (its first ``pam_length`` bases are
       dropped before searching);
    2. the template is the protospacer alone.

    Returns ``(hits, rejected)``: hits = sorted [(start, end, strand)] spanning
    PAM + protospacer; rejected = [(start, end, strand, pam_seen)] for
    protospacer matches whose PAM does not match (empty when hits is not).
    """
    hits, rejected = [], []
    if not pam_motif:
        for strand, target in [('+', seq), ('-', revcomp(seq))]:
            start = ref.find(target)
            while start >= 0:
                hits.append((start, start + len(target), strand))
                start = ref.find(target, start + 1)
        return hits, rejected

    if len(seq) <= pam_length:
        return hits, rejected
    for proto in (seq[pam_length:], seq):
        cur_hits, cur_rejected = [], []
        for strand, target in [('+', proto), ('-', revcomp(proto))]:
            start = ref.find(target)
            while start >= 0:
                end = start + len(target)
                if strand == '+':
                    a, b = start - pam_length, end
                    pam_5to3 = ref[a:start] if a >= 0 else ''
                else:
                    a, b = start, end + pam_length
                    pam_ref = ref[end:b] if b <= len(ref) else ''
                    pam_5to3 = revcomp(pam_ref) if pam_ref else ''
                if len(pam_5to3) == pam_length and pam_matches(pam_5to3, pam_motif):
                    cur_hits.append((a, b, strand))
                else:
                    cur_rejected.append((max(a, 0), min(b, len(ref)), strand, pam_5to3 or 'out of reference'))
                start = ref.find(target, start + 1)
        if cur_hits:
            return sorted(cur_hits), []
        if not rejected:
            rejected = cur_rejected  # report the first reading that found anything
    return [], rejected


def format_guide_hit(ref, start, end, strand, pam_length, pam_motif):
    """One-line description of a target site for the text report (1-based)."""
    pam_iv = pam_interval(start, end, strand, pam_length)
    proto_iv = protospacer_interval(start, end, strand, pam_length)
    pam_txt, pam_seq_txt, crrna_txt = '?-?', '', '?-?'
    if pam_iv is not None:
        pam_txt = '%d-%d' % (pam_iv[0] + 1, pam_iv[1])
        seq = pam_sequence(ref, start, end, strand, pam_length)
        ok = 'matches' if pam_matches(seq, pam_motif) else 'DOES NOT match'
        pam_seq_txt = ' [%s%s]' % (seq, (', ' + ok + ' ' + pam_motif) if pam_motif else '')
    if proto_iv is not None:
        crrna_txt = '%d-%d' % (proto_iv[0] + 1, proto_iv[1])
    return '%d-%d (%s); PAM %s%s; crRNA %s (%s)' % (
        start + 1, end, strand, pam_txt, pam_seq_txt, crrna_txt, crrna_strand(strand))


def locate_cds(ref, cds_seq):
    """Place the CDS on the reference.

    An exact match on either strand gives one block. Otherwise the CDS is
    aligned locally (gaps allowed, so introns become reference gaps) and its
    aligned segments of >= 12 bp are taken as exons; the strand is left
    unknown (None) if both orientations score the same.

    Returns ``(blocks, strand, note)``.
    """
    for strand, target in [('+', cds_seq), ('-', revcomp(cds_seq))]:
        start = ref.find(target)
        if start >= 0:
            return [(start, start + len(target))], strand, ''
    aligner = PairwiseAligner(mode='local')
    aligner.match_score = 2
    aligner.mismatch_score = -4
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.3
    options = [aligner.align(ref, x)[0] for x in (cds_seq, revcomp(cds_seq))]
    best_index = max(range(len(options)), key=lambda i: options[i].score)
    strand = ('+', '-')[best_index] if options[0].score != options[1].score else None
    blocks = [(int(s), int(e)) for s, e in options[best_index].aligned[0] if e - s >= 12]
    return blocks, strand, 'CDS segments inferred by alignment; verify exon boundaries.'


def locate_primers(ref, primers):
    """[(name, [(start, end, strand)])]: first exact match on each strand."""
    out = []
    for name, seq in primers:
        hits = []
        for strand, s in [('+', seq), ('-', revcomp(seq))]:
            at = ref.find(s)
            if at >= 0:
                hits.append((at, at + len(s), strand))
        out.append((name, hits))
    return out


def locate_guides(ref, guides, pam_length, pam_motif):
    """[(name, [(start, end, strand)])] for every |crispr record, printing a
    warning for each target that is too short, absent or has a bad PAM."""
    out = []
    for name, seq in guides:
        if len(seq) <= pam_length:
            print(f'WARNING: CRISPR target {name} too short to contain PAM and protospacer')
        hits, rejected = locate_guide(ref, seq, pam_length, pam_motif)
        for rs, r_end, rstrand, pam_seen in rejected:
            print(f"WARNING: CRISPR target {name}: protospacer found at {rs + 1}-{r_end} ({rstrand}) "
                  f"but PAM '{pam_seen}' does not match {pam_motif}; site ignored")
        if pam_motif and not hits:
            if rejected:
                print(f'WARNING: CRISPR target {name}: no site with a valid {pam_motif} PAM '
                      f'(see the PAM(s) observed above)')
            else:
                print(f'WARNING: CRISPR target {name}: protospacer NOT found in the reference '
                      f'(exact match required, both strands, with or without the {pam_length}-nt PAM '
                      f'prefix); check the template sequence, not the PAM')
        out.append((name, hits))
    return out


def load_reference(path, pam_length, pam_motif):
    """Read the template FASTA and annotate the reference.

    Raises ValueError unless the template has exactly one ``|gene`` record.
    Only the first ``|cds`` record is used.
    """
    template = records(path)
    genes = [(name, seq) for name, seq in template if name.lower().endswith('|gene')]
    if len(genes) != 1:
        raise ValueError('Template must contain exactly one |gene record')
    gene_id, seq = genes[0]
    ref = Reference(gene_id=gene_id, seq=seq)

    cds = [s for n, s in template if n.lower().endswith('|cds')]
    if cds:
        ref.cds_length = len(cds[0])
        ref.cds_blocks, ref.cds_strand, ref.cds_note = locate_cds(seq, cds[0])

    primers = [(n.split('|')[0], s) for n, s in template if n.lower().endswith('|primer')]
    ref.primers = locate_primers(seq, primers)

    guides = [(n.split('|')[0], s) for n, s in template if n.lower().endswith('|crispr')]
    ref.guides = locate_guides(seq, guides, pam_length, pam_motif)
    if pam_motif:
        ref.pam_sites = find_pam_sites(seq, pam_motif)
    return ref


def get_highlight_intervals(mode, guides):
    """Violet overlay for --highlight: 'between' -> from the first to the last
    target (or the single target), 'targets' -> each target, 'none' -> nothing."""
    intervals = sorted((a, b) for _, hits in guides for a, b, _ in hits)
    if mode == 'targets':
        return intervals
    if mode == 'between' and len(intervals) >= 2:
        return [(intervals[0][0], intervals[-1][1])]
    if mode == 'between':
        return intervals
    return []


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def draw_reference_panel(ax, ref, highlight, font_size, pam_length=4):
    """Draw the reference track at y = 0 (y grows downwards).

    Backbone (grey) with the --highlight overlay (violet); CDS blocks above it,
    the last block in reading direction drawn as an arrow; primers below, as
    arrows in the direction of extension; for each CRISPR target, its PAM
    (gold box plus a tick visible at gene scale), its crRNA as an arrow along
    the hybridised strand, and its name.
    """
    n = len(ref.seq)
    ax.add_patch(Rectangle((0, -0.11), n, 0.22, color=GREY, zorder=2))
    for a, b in highlight:
        ax.add_patch(Rectangle((a, -0.11), b - a, 0.22, color=VIOLET, zorder=3))

    blocks = sorted(ref.cds_blocks)
    terminal = None
    if blocks and ref.cds_strand in ('+', '-'):
        terminal = blocks[-1 if ref.cds_strand == '+' else 0]
    for a, b in blocks:
        if (a, b) == terminal:
            direction = 1 if ref.cds_strand == '+' else -1
            ax.add_patch(FancyArrow(
                a if direction == 1 else b, 0.1675, direction * (b - a), 0,
                width=0.055, head_width=0.105,
                head_length=min(n * 0.015, (b - a) * 0.4),
                length_includes_head=True, color=CDS_GREY, linewidth=0, zorder=2))
        else:
            ax.add_patch(Rectangle((a, 0.14), b - a, 0.055, color=CDS_GREY, zorder=2))

    for name, hits in ref.primers:
        for a, b, strand in hits:
            forward = strand == '+'
            ax.annotate(
                '', xy=(b if forward else a, -0.16), xytext=(a if forward else b, -0.16),
                arrowprops=dict(arrowstyle='-|>', color='#333333', lw=1.5, mutation_scale=11),
                zorder=6)
            ax.text((a + b) / 2, -0.25, name, ha='center', va='bottom', fontsize=font_size * 0.72)

    for name, hits in ref.guides:
        for a, b, strand in hits:
            pam = pam_interval(a, b, strand, pam_length)
            if pam is not None:
                ps, pe = pam
                ax.add_patch(Rectangle((ps, -0.11), pe - ps, 0.22,
                                       facecolor=PAM_GOLD, edgecolor='none', zorder=7))
                ax.plot([(ps + pe) / 2] * 2, [-0.13, 0.12], color=PAM_GOLD, lw=2.0, zorder=8)

            crrna = crrna_display_interval(a, b, strand, pam_length)
            if crrna is not None:
                cs, ce = crrna
                direction = 1 if crrna_strand(strand) == '+' else -1
                ax.add_patch(FancyArrow(
                    cs if direction == 1 else ce, 0.29, direction * (ce - cs), 0,
                    width=0.018, head_width=0.06,
                    head_length=min(max((ce - cs) * 0.35, 10), n * 0.014),
                    length_includes_head=True, color=PURPLE_DARK, linewidth=0, zorder=6))
            ax.text((a + b) / 2, 0.41, name, ha='center', va='top',
                    fontsize=font_size * 0.8, fontweight='bold', color=PURPLE_DARK)


def draw_sample_row(ax, y, hits, ref_len, args, highlight, wt_events):
    """Draw one sample at height ``y`` from its read alignments ``hits``.

    The union of the aligned spans of passing reads is a solid band (masked
    bases do not punch holes in it), overlaid in violet by --highlight, and
    veiled in white where read support is weak (see reliable_intervals).
    Deletions >= --min-indel are hatched boxes labelled 'Δn', insertions are
    vertical ticks labelled '+n'; '*' marks an event also called in WT.
    """
    valid = passing_hits(hits, args.min_aligned)
    if not valid:
        ax.text(ref_len / 2, y, 'No aligned reads', ha='center', va='center',
                fontsize=args.font_size * 0.75, color='#777777')
        return

    allspan = merge([iv for hit in valid for iv in hit['span']])
    reliable = merge([iv for hit in valid for iv in hit['reliable']])
    for a, b in allspan:
        ax.add_patch(Rectangle((a, y - 0.105), b - a, 0.21, color=GREY, zorder=2))
        for hs, he in highlight:
            left, right = max(a, hs), min(b, he)
            if right > left:
                ax.add_patch(Rectangle((left, y - 0.105), right - left, 0.21, color=VIOLET, zorder=3))
    for a, b in subtract_intervals(allspan, reliable):
        ax.add_patch(Rectangle((a, y - 0.105), b - a, 0.21, facecolor='white',
                               edgecolor='none', alpha=0.55, zorder=4))

    for key in collect_events(valid):
        kind, a, b = key
        size = b - a if kind == 'deletion' else b
        if size < args.min_indel:
            continue
        mark = '*' if key in wt_events else ''
        if kind == 'deletion':
            ax.add_patch(Rectangle((a, y - 0.13), b - a, 0.26, facecolor='white',
                                   edgecolor='#3c3c3c', linewidth=0, hatch='///', zorder=5))
            ax.text((a + b) / 2, y, f'Δ{size}{mark}', ha='center', va='center',
                    fontsize=args.font_size * 0.75, color=INK,
                    bbox=dict(facecolor='white', edgecolor='none', boxstyle='square,pad=0.15'),
                    zorder=6)
            ax.add_patch(Rectangle((a, y - 0.13), b - a, 0.26, facecolor='none',
                                   edgecolor='#3c3c3c', linewidth=1.1, zorder=7))
        else:
            ax.plot([a, a], [y - 0.14, y + 0.14], color='#333333', lw=2, zorder=5)
            ax.text(a, y - 0.17, f'+{size}{mark}', ha='center', va='bottom',
                    fontsize=args.font_size * 0.75, color=INK)


def legend_handles(pam_motif=None):
    """Legend entries shared by all edit maps."""
    return [
        Patch(facecolor=GREY, label='Aligned region'),
        Patch(facecolor=GREY_LIGHT, label='Low-confidence read quality'),
        Patch(facecolor=VIOLET, label='CRISPR target region'),
        Patch(facecolor='white', edgecolor='#333333', hatch='///', label='Deletion'),
        Patch(facecolor=CDS_GREY, label='CDS'),
        Patch(facecolor=PAM_GOLD, label=f'PAM ({pam_motif})' if pam_motif else 'PAM'),
        Line2D([0, 1], [0, 0], color=PURPLE_DARK, lw=2.2, marker='>', markevery=[1],
               label='crRNA target strand'),
    ]


def finish_edit_map(fig, ax, ticks, labels, y_bottom, ref, args, stem):
    """Axes, title and legend common to all edit maps; write ``<stem>.pdf|png``."""
    n = len(ref.seq)
    ax.set_yticks(ticks, labels)
    ax.set_xlim(-n * 0.012, n * 1.012)
    ax.set_ylim(y_bottom, -0.48)
    ax.tick_params(axis='both', labelsize=args.font_size, length=0)
    ax.tick_params(axis='x', length=4)
    ax.set_xlabel(f'{ref.name} reference position (bp)', fontsize=args.font_size)
    ax.set_title('CRISPR editing — Sanger alignment overview', fontsize=args.font_size + 2, pad=14)
    ax.grid(axis='x', alpha=0.1, zorder=0)
    fig.legend(handles=legend_handles(args.pam_motif), loc='lower center',
               bbox_to_anchor=(0.5, 0.015), ncol=4, frameon=False,
               fontsize=args.font_size * 0.82)
    fig.tight_layout(rect=(0, 0.72 / fig.get_figheight(), 1, 1))
    fig.savefig(args.out / f'{stem}.pdf', bbox_inches='tight')
    fig.savefig(args.out / f'{stem}.png', dpi=220, bbox_inches='tight')
    plt.close(fig)


def sample_label(sample):
    """Y-axis label of a sample: '2_75' -> '2.75'."""
    return sample.replace('_', '.')


# ---------------------------------------------------------------------------
# Tables and report
# ---------------------------------------------------------------------------

SUMMARY_FIELDS = ['sample', 'reads', 'passing_reads', 'covered_bp', 'reference_bp', 'events', 'status']
EVENT_FIELDS = ['sample', 'event', 'start_1based', 'end_1based', 'size', 'supporting_reads',
                'nearby_crispr', 'also_in_wt']


def build_summary_and_details(samples, groups, min_aligned, guides, wt_events, ref_length):
    """Rows of summary.tsv (one per sample) and events.tsv (one per event).

    Event coordinates are 1-based: a deletion covers start..end inclusive; an
    insertion lies immediately before reference position start (start = end).
    Every event is reported, whatever --min-indel. An event is 'nearby' a
    CRISPR target when within NEARBY_TARGET_BP of it.
    """
    summary, details = [], []
    for sample in samples:
        hits = groups[sample]
        valid = passing_hits(hits, min_aligned)
        allcov = merge([iv for hit in valid for iv in hit['coverage']])
        events = collect_events(valid)

        for (kind, a, b), files in events.items():
            nearby = []
            for gname, ghits in guides:
                for gs, ge, _ in ghits:
                    end = b if kind == 'deletion' else a
                    if max(gs - end, a - ge, 0) <= NEARBY_TARGET_BP:
                        nearby.append(gname)
            details.append(dict(
                sample=sample,
                event=kind,
                start_1based=a + 1,
                end_1based=b if kind == 'deletion' else a + 1,
                size=b - a if kind == 'deletion' else b,
                supporting_reads=len(files),
                nearby_crispr=','.join(sorted(set(nearby))),
                also_in_wt='yes' if (kind, a, b) in wt_events else 'no',
            ))

        if not valid:
            status = 'No passing reads'
        elif events:
            status = 'Review alignments'
        else:
            status = 'No indel detected in aligned regions'
        summary.append(dict(
            sample=sample,
            reads=len(hits),
            passing_reads=len(valid),
            covered_bp=sum(e - a for a, e in allcov),
            reference_bp=ref_length,
            events=len(events),
            status=status,
        ))
    return summary, details


def write_tsv(path, rows, fields):
    with path.open('w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def write_report(path, ref, args, groups, summary, details):
    """Write the human-readable report of one band."""
    n = len(ref.seq)
    with path.open('w') as f:
        f.write(
            f"Reference: {ref.gene_id} ({n} bp)\n"
            f"CDS: {ref.cds_length if ref.cds_length is not None else 'not supplied'} bp. {ref.cds_note}\n"
            f"PAM motif: {args.pam_motif or 'not specified'} (length {args.pam_length} bp); "
            + ('searched in the reference 5-prime of each protospacer; non-matching sites are excluded\n'
               if args.pam_motif else 'not used for searching\n')
        )
        if args.pam_motif:
            n_plus = sum(1 for _, _, s in ref.pam_sites if s == '+')
            f.write(f'{args.pam_motif} occurrences in the reference: {len(ref.pam_sites)} '
                    f'({n_plus} on + strand, {len(ref.pam_sites) - n_plus} on - strand)\n')
        for name, hits in ref.guides:
            if hits:
                desc = '; '.join(format_guide_hit(ref.seq, s, e, st, args.pam_length, args.pam_motif)
                                 for s, e, st in hits)
            elif args.pam_motif:
                desc = 'no site with a valid PAM found'
            else:
                desc = 'not found exactly'
            f.write(f'CRISPR target {name}: {desc}\n')
        for name, hits in ref.primers:
            where = ', '.join(f'{s + 1}-{e} ({st})' for s, e, st in hits) if hits else 'not found exactly'
            f.write(f'Primer {name}: {where}\n')
        f.write('An asterisk (*) on the figure denotes an identical indel also called in WT; '
                'such events may be alignment artefacts.\n')

        for row in summary:
            sample = row['sample']
            f.write(f"\n{sample}: {row['reads']} reads; {row['passing_reads']} passing; "
                    f"{row['covered_bp']}/{n} bp covered\n")
            for hit in groups[sample]:
                intervals = ', '.join(f'{a + 1}-{b}' for a, b in hit['coverage'])
                qc = hit.get('qc')
                qc_text = ''
                if qc:
                    qc_text = (f"; quality: mean Q {qc['mean_q_usable']} over {qc['usable_length']} usable bp "
                               f"(ends trimmed {qc['trimmed_5prime']}/{qc['trimmed_3prime']}, "
                               f"{qc['masked_low_q']} low-Q bases masked), status {qc['status']}")
                f.write(f"  {Path(hit['file']).name}: {hit['strand']}, aligned {hit['aligned_bases']} bp, "
                        f"mismatches {len(hit['mismatches'])}; terminal Ns trimmed (read 5prime/3prime): "
                        f"{hit['terminal_n_5prime']}/{hit['terminal_n_3prime']}{qc_text}; "
                        f"reference intervals {intervals}\n")
            for event in details:
                if event['sample'] == sample:
                    f.write(f"  {event['event']}: {event['start_1based']}-{event['end_1based']}; "
                            f"{event['size']} bp; supporting reads {event['supporting_reads']}; "
                            f"nearby CRISPR target(s): {event['nearby_crispr'] or 'none'}; "
                            f"also in WT: {event['also_in_wt']}\n")
            if row['passing_reads'] and not row['events']:
                f.write('  No indel detected within aligned regions (uncovered positions are unassessed).\n')

        f.write('\nCAUTION: Indels are candidate calls, not validated genotypes. '
                'Local alignment can shift breakpoints in repeats; discordant reads, mixed '
                'chromatograms, and heterokaryosis require manual review.'
                + ('' if args.ab1 else ' FASTA alone does not encode Sanger quality scores.')
                + '\n')


def render(args, groups, plotted, band, ref, highlight):
    """Edit map, summary/events tables and report of one band ('' = no bands)."""
    suffix = '_' + band if band else ''
    if not groups:
        print(f'WARNING: no reads for band {band}; no output generated')
        return
    if not plotted:
        print(f'WARNING: no selected samples for band {band}; reference-only figure')

    wt_events = compute_wt_events(groups, args.min_aligned)
    row = args.row_height
    fig, ax = plt.subplots(figsize=(10, max(3.6, 2.05 + row * len(plotted))))
    draw_reference_panel(ax, ref, highlight, args.font_size, args.pam_length)
    ticks = [0]
    for i, sample in enumerate(plotted):
        y = 0.85 + i * row
        draw_sample_row(ax, y, groups.get(sample, []), len(ref.seq), args, highlight, wt_events)
        ticks.append(y)
    y_bottom = 0.85 + (max(1, len(plotted)) - 1) * row + 0.35
    finish_edit_map(fig, ax, ticks, ['Reference'] + [sample_label(s) for s in plotted],
                    y_bottom, ref, args, f'edit_map{suffix}')

    summary, details = build_summary_and_details(
        sorted_samples(groups), groups, args.min_aligned, ref.guides, wt_events, len(ref.seq))
    write_tsv(args.out / f'summary{suffix}.tsv', summary, SUMMARY_FIELDS)
    write_tsv(args.out / f'events{suffix}.tsv', details, EVENT_FIELDS)
    write_report(args.out / f'report{suffix}.txt', ref, args, groups, summary, details)
    print('Results:', args.out.resolve())


def render_combined(args, band_groups, plotted, ref, highlight):
    """edit_map_combined.pdf|png: upper band (1) above lower band (2), each
    sample row computed against its own band's WT."""
    if not band_groups['1'] and not band_groups['2']:
        return
    if not plotted['1'] and not plotted['2']:
        print('WARNING: no selected samples for combined view; no output generated')
        return

    row = args.row_height
    n_rows = max(1, len(plotted['1']) + len(plotted['2']) + 1)
    fig, ax = plt.subplots(figsize=(10, max(4.4, 2.25 + row * n_rows)))
    draw_reference_panel(ax, ref, highlight, args.font_size, args.pam_length)

    ticks, labels = [0], ['Reference']
    label_x = -0.08  # in axes coordinates: just left of the y tick labels

    def draw_band(band, start):
        groups = band_groups[band]
        wt_events = compute_wt_events(groups, args.min_aligned)
        positions = []
        for i, sample in enumerate(plotted[band]):
            y = start + i * row
            draw_sample_row(ax, y, groups.get(sample, []), len(ref.seq), args, highlight, wt_events)
            positions.append(y)
            ticks.append(y)
            labels.append(sample_label(sample))
        return positions

    def band_label(text, positions):
        ax.text(label_x, sum(positions) / len(positions), text, transform=ax.get_yaxis_transform(),
                ha='right', va='center', fontsize=args.font_size * 0.88, fontweight='bold',
                color='#666666', clip_on=False)

    upper_start = 0.85
    upper = draw_band('1', upper_start)
    if upper:
        band_label('Upper\nband', upper)
    separator_y = upper[-1] + row * 0.58 if upper else upper_start + row * 0.25
    ax.hlines(separator_y, 0, len(ref.seq), colors='#bbbbbb', linewidth=1.1, zorder=1)
    lower = draw_band('2', separator_y + row * 0.72)
    if lower:
        band_label('Lower\nband', lower)

    last_y = (lower or upper or [0])[-1]
    finish_edit_map(fig, ax, ticks, labels, last_y + 0.45, ref, args, 'edit_map_combined')


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------

def build_parser():
    p = argparse.ArgumentParser(
        description=__doc__.split('\n\n')[0],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    io = p.add_argument_group('input / output')
    io.add_argument('--template', required=True, type=Path,
                    help='FASTA with one <id>|gene record and optional |cds, |primer and |crispr records')
    io.add_argument('--reads', required=True, type=Path, nargs='+',
                    help='Read files or folders (searched recursively): .seq (FASTA) files, '
                         'or .ab1 chromatograms with --ab1. *.unclipped.seq files are ignored')
    io.add_argument('--out', type=Path, default=Path('results'), help='Output folder')
    io.add_argument('--verbose', action='store_true',
                    help='Print the length and ambiguous-base profile of every read')

    q = p.add_argument_group('.ab1 quality control (--ab1)')
    q.add_argument('--ab1', action='store_true',
                   help='Read .ab1 chromatograms instead of .seq files: low-quality ends are '
                        'trimmed, remaining low-quality bases are masked as N, and qc_reads.tsv/.pdf/.png '
                        'are written')
    q.add_argument('--min-quality', type=int, default=20,
                   help='Phred threshold for trimming and masking; 0 disables both')
    q.add_argument('--trim-mode', choices=['window', 'mott', 'none'], default='window',
                   help='End trimming: window = sliding-window mean >= --min-quality at both ends '
                        '(isolated low-quality bases are masked, not cut); mott = modified Mott '
                        'algorithm (single best window, harsh on scattered low-quality bases); '
                        'none = per-base masking only')
    q.add_argument('--trim-window', type=int, default=10, help='Window size (bases) for --trim-mode window')
    q.add_argument('--mixed-ratio', type=float, default=0.3,
                   help='A base counts as a double peak when the strongest competing channel reaches '
                        'this fraction of the called one (QC table only)')

    g = p.add_argument_group('grouping')
    g.add_argument('--group-regex', default=r'^(?P<sample>\d+_\d+)_\d+_\d+$',
                   help='Regex on the file stem; its "sample" group (or first group) is the sample ID. '
                        'WT_* files are always grouped as WT. Ignored for numbered names with --split-bands')
    g.add_argument('--split-bands', action='store_true',
                   help='Separate gel bands 1 and 2 (third field of the file name); numbered samples '
                        'without a band field go to band 1, WT goes to both, other names are excluded')

    a = p.add_argument_group('alignment and event calling')
    a.add_argument('--gap-open', type=float, default=-2.0, help='Gap opening score (negative)')
    a.add_argument('--gap-extend', type=float, default=-0.05, help='Gap extension score (negative)')
    a.add_argument('--min-flank', type=int, default=20,
                   help='Aligned blocks with fewer unambiguous bases are dropped and the gaps around '
                        'them merged into one net event, which removes spurious small indels in noisy '
                        'regions; 0 disables')
    a.add_argument('--min-aligned', type=int, default=10,
                   help='Minimum unambiguous aligned bases for a read to pass; 1 keeps any alignment')

    t = p.add_argument_group('CRISPR targets')
    t.add_argument('--pam-length', type=int, default=4,
                   help="Length of the 5' PAM at the start of each |crispr sequence")
    t.add_argument('--pam-motif', default=None,
                   help="IUPAC motif of a 5' PAM (e.g. TTTV for Cas12a), checked in the reference "
                        "immediately 5' of each protospacer; sites whose PAM does not match are "
                        "dropped. Overrides --pam-length. 3' PAMs (e.g. NGG for SpCas9) are not supported")

    f = p.add_argument_group('figure')
    f.add_argument('--samples', nargs='+', default=None,
                   help='Samples to plot, in this order (e.g. 3_14 3.16 WT); tables and reports keep all samples')
    for band_id in ('1', '2'):
        f.add_argument(f'--samples-{band_id}', nargs='+', default=None,
                       help=f'Only plot these samples in band {band_id} (needs --split-bands)')
        f.add_argument(f'--hide-samples-{band_id}', nargs='+', default=[],
                       help=f'Hide these samples in band {band_id} (needs --split-bands)')
    f.add_argument('--min-indel', type=int, default=1,
                   help='Smallest indel drawn (e.g. 5 hides 1-4 bp); tables keep all events')
    f.add_argument('--highlight', choices=['between', 'targets', 'none'], default='between',
                   help='Violet overlay: span between the first and last target, targets alone, or none')
    f.add_argument('--conf-window', type=int, default=20,
                   help='Window (bp) used to judge read support along the aligned region')
    f.add_argument('--conf-min-frac', type=float, default=0.7,
                   help='Parts of the aligned region with a lower fraction of unmasked bases in the '
                        'window are drawn lighter (never removed)')
    f.add_argument('--font-size', type=float, default=13, help='Base font size')
    f.add_argument('--row-height', type=float, default=0.44, help='Spacing between sample rows')
    f.add_argument('--hide-snps', action='store_true',
                   help='No effect (SNPs are never drawn); kept for compatibility with older scripts')
    return p


def parse_args(argv=None):
    """Parse and validate the command line; normalise sample IDs ('.' -> '_',
    duplicates removed) and the PAM motif."""
    p = build_parser()
    args = p.parse_args(argv)
    if args.pam_motif:
        args.pam_motif = args.pam_motif.upper().replace('U', 'T')
        if not re.fullmatch(r'[ACGTRYSWKMBDHVN]+', args.pam_motif):
            p.error('--pam-motif must use the IUPAC DNA alphabet (ACGTRYSWKMBDHVN)')
        args.pam_length = len(args.pam_motif)
    if args.pam_length < 1:
        p.error('--pam-length must be positive')
    if not args.split_bands and any((args.samples_1, args.samples_2, args.hide_samples_1, args.hide_samples_2)):
        p.error('Band-specific filters require --split-bands')
    for attr in ('samples', 'samples_1', 'samples_2', 'hide_samples_1', 'hide_samples_2'):
        values = getattr(args, attr)
        if values is not None:
            setattr(args, attr, list(dict.fromkeys(v.replace('.', '_') for v in values)))
    if args.gap_open >= 0 or args.gap_extend >= 0:
        p.error('Gap scores must be negative')
    if args.min_aligned < 1:
        p.error('--min-aligned must be at least 1')
    if args.min_flank < 0:
        p.error('--min-flank must be >= 0')
    if args.conf_window < 1 or not 0 < args.conf_min_frac <= 1:
        p.error('--conf-window must be >= 1 and --conf-min-frac in (0, 1]')
    if args.min_quality < 0:
        p.error('--min-quality must be >= 0')
    if args.mixed_ratio <= 0:
        p.error('--mixed-ratio must be positive')
    if args.trim_window < 1:
        p.error('--trim-window must be at least 1')
    return p, args


def collect_read_paths(inputs, ab1):
    """Read files among ``inputs`` (files, or folders searched recursively).
    Raises FileNotFoundError for an input that does not exist."""
    ext = '.ab1' if ab1 else '.seq'

    def is_read_file(f):
        return f.suffix.lower() == ext and not f.name.lower().endswith('.unclipped.seq')

    paths = []
    for path in inputs:
        if path.is_file():
            if is_read_file(path):
                paths.append(path)
            else:
                print(f'WARNING: non-{ext} read file ignored: {path}')
        elif path.is_dir():
            paths += [f for f in path.rglob('*') if f.is_file() and is_read_file(f)]
        else:
            raise FileNotFoundError(f'Cannot access {path}')
    return sorted(set(paths))


def load_and_align_reads(paths, args, ref_seq):
    """Load, align and group every read.

    Returns ``(band_groups, qc_rows)``: band_groups maps a band ('1', '2', or
    '' without --split-bands) to ``{sample: [alignment dict]}``; qc_rows holds
    one qc dict per .ab1 file (empty in FASTA mode).
    """
    qc_rows = []
    bands = ('1', '2') if args.split_bands else ('',)
    band_groups = {b: defaultdict(list) for b in bands}
    for f in paths:
        if args.split_bands:
            sample, band = band_name(f, args.group_regex)
            if sample == 'WT':
                destinations = ['1', '2']
            elif band is not None:
                destinations = [band]
            else:
                print(f'WARNING: unassigned band, excluded: {f}')
                continue
        else:
            sample = sample_name(f, args.group_regex)
            destinations = ['']
        try:
            read_records = load_reads(f, args)
        except (ValueError, struct.error, OSError) as exc:
            print(f'WARNING: cannot read {f}: {exc}')
            continue
        if not read_records:
            print(f'WARNING: no FASTA records in {f}')
        for name, seq, qc in read_records:
            if args.verbose:
                print(f'\nREAD: {f.name}')
                read_diagnostics(seq)
            if qc is not None:
                qc.update(sample=sample, strand='', aligned_bp=0)
                qc_rows.append(qc)
            if not seq or not any(informative(b) for b in seq):
                print(f'WARNING: empty or uninformative read: {f} ({name})'
                      + (f" - quality status: {qc['status']}" if qc else ''))
                continue
            hit = align(ref_seq, seq.strip('N'), args.gap_open, args.gap_extend, args.min_flank,
                        args.conf_window, args.conf_min_frac)
            if hit['aligned_bases'] < args.min_aligned:
                print(f"WARNING: insufficient alignment: {f} ({name}), {hit['aligned_bases']} bp")
            hit.update(
                file=str(f), record=name, qc=qc,
                terminal_n_5prime=len(seq) - len(seq.lstrip('N')),
                terminal_n_3prime=len(seq) - len(seq.rstrip('N')),
            )
            if qc is not None:
                qc.update(strand=hit['strand'], aligned_bp=hit['aligned_bases'])
            for band in destinations:
                band_groups[band][sample].append(hit)
    return band_groups, qc_rows


def main(argv=None):
    parser, args = parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    try:
        ref = load_reference(args.template, args.pam_length, args.pam_motif)
        paths = collect_read_paths(args.reads, args.ab1)
    except (ValueError, FileNotFoundError) as exc:
        parser.error(str(exc))
    if not paths:
        ext = '.ab1' if args.ab1 else '.seq'
        parser.error(f'No {ext} reads found' + ('' if args.ab1 else ' (use --ab1 for .ab1 chromatograms)'))

    band_groups, qc_rows = load_and_align_reads(paths, args, ref.seq)
    if args.ab1:
        write_qc(args, qc_rows)

    found = set().union(*band_groups.values())
    if args.split_bands and args.samples is None:
        # Same rows, in the same order, in every band.
        args.samples = sorted_samples(found)
    if args.samples:
        missing = [s for s in args.samples if s not in found]
        if missing:
            print('WARNING: samples not found:', ', '.join(missing))

    highlight = get_highlight_intervals(args.highlight, ref.guides)
    plotted = {band: select_plotted_samples(args, groups, band) for band, groups in band_groups.items()}
    for band, groups in band_groups.items():
        render(args, groups, plotted[band], band, ref, highlight)
    if args.split_bands:
        render_combined(args, band_groups, plotted, ref, highlight)


if __name__ == '__main__':
    main()
