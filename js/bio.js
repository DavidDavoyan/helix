/**
 * bio.js — the biology. No rendering, no THREE, no DOM.
 *
 * Every number in this file is a measured property of B-form DNA or a
 * property of the standard genetic code. If a value here is wrong, the
 * simulation is wrong, so each one carries its source in a comment.
 */

/* ------------------------------------------------------------------ B-DNA */

/**
 * B-DNA helical parameters (Watson-Crick, fibre-diffraction/crystal averages).
 *
 * RISE and TWIST are the canonical pair: 3.38 A of climb and 34.3 deg of
 * rotation per base pair, which is 10.5 bp per full turn and 35.5 A of pitch.
 *
 * MINOR_GROOVE_SPAN is the angle, looking down the helix axis, from one
 * strand's glycosidic bond to the other's, going the short way round. 127 deg
 * is not a free parameter — it is pinned twice over:
 *
 *   1. The two grooves then span 127 deg and 233 deg, a ratio of 1.83, and the
 *      measured groove widths are 12 A and 22 A, a ratio of 1.83.
 *   2. With C1' atoms on a 5.9 A radius, the chord between them comes out at
 *      2 * 5.9 * sin(63.5 deg) = 10.56 A, and the measured C1'-C1' distance in
 *      a Watson-Crick pair is 10.5 A.
 *
 * Getting this angle wrong is the single most common error in DNA artwork: put
 * the strands 180 deg apart and you get two identical grooves, which erases the
 * whole reason proteins can read the sequence without opening the helix.
 */
export const BDNA = {
  RISE: 3.38,               // A per base pair, along the axis
  TWIST: 34.3,              // deg per base pair -> 10.5 bp/turn
  MINOR_GROOVE_SPAN: 127,   // deg between glycosidic bonds, short way round
  R_PHOSPHATE: 8.9,         // A, phosphorus from the axis
  R_SUGAR: 7.1,             // A, deoxyribose ring centre from the axis
  R_GLYCOSIDIC: 5.9,        // A, C1' from the axis
  PROPELLER: -11,           // deg, the two bases of a pair are not coplanar
  HELIX_RADIUS: 10.0,       // A, van der Waals half-width -> the famous 20 A
};

export const BP_PER_TURN = 360 / BDNA.TWIST;      // 10.5
export const PITCH = BP_PER_TURN * BDNA.RISE;     // 35.5 A

/* ------------------------------------------------------------------- bases */

/**
 * Purines are bicyclic, pyrimidines monocyclic. A purine always pairs with a
 * pyrimidine, so every rung of the ladder is the same width — which is why the
 * backbone can be a smooth helix at all.
 *
 * `reach` is how far the ring system extends from C1' towards its partner, and
 * it is the length of the drawn plate, measured off the generated geometry:
 * 4.7 A across a fused bicycle built at the 1.39 A aromatic bond length, and
 * 2.65 A across a single six-ring.
 *
 * Those two numbers are what leaves the gap in the middle. With C1' atoms
 * 10.56 A apart, the space between the plates comes out at
 * 10.56 - 4.7 - 2.65 = 3.2 A, against a real Watson-Crick hydrogen bond of
 * 2.8-2.9 A. No bond length is imposed anywhere — the gap falls out of the
 * ring sizes, and it lands 0.3 A long because the plate is drawn from C1'
 * rather than from the glycosidic nitrogen 1.5 A further in.
 */
export const BASES = {
  A: { name: 'Adenine',  kind: 'purine',     pairsWith: 'T', rnaPairsWith: 'U', hbonds: 2, reach: 4.7, colour: 0x53e08a, ring: 'fused' },
  T: { name: 'Thymine',  kind: 'pyrimidine', pairsWith: 'A', rnaPairsWith: 'A', hbonds: 2, reach: 2.65, colour: 0xff6b6b, ring: 'single' },
  G: { name: 'Guanine',  kind: 'purine',     pairsWith: 'C', rnaPairsWith: 'C', hbonds: 3, reach: 4.7, colour: 0xffc04d, ring: 'fused' },
  C: { name: 'Cytosine', kind: 'pyrimidine', pairsWith: 'G', rnaPairsWith: 'G', hbonds: 3, reach: 2.65, colour: 0x5aa9ff, ring: 'single' },
  U: { name: 'Uracil',   kind: 'pyrimidine', pairsWith: 'A', rnaPairsWith: 'A', hbonds: 2, reach: 2.65, colour: 0xc98bff, ring: 'single' },
};

export const DNA_LETTERS = ['A', 'T', 'G', 'C'];

/** Watson-Crick partner on a DNA strand. */
export function complement(b) {
  return BASES[b] ? BASES[b].pairsWith : 'N';
}

/** Watson-Crick partner when the new strand is RNA: A is read as U. */
export function rnaComplement(b) {
  return { A: 'U', T: 'A', G: 'C', C: 'G' }[b] || 'N';
}

/** Reverse complement, i.e. the other strand written 5'->3'. */
export function reverseComplement(seq) {
  let out = '';
  for (let i = seq.length - 1; i >= 0; i--) out += complement(seq[i]);
  return out;
}

/** Strip anything that is not ATGC and upper-case the rest. */
export function cleanSequence(text) {
  return (text || '').toUpperCase().replace(/[^ATGC]/g, '');
}

/**
 * GC content drives melting temperature, because G-C is a three-hydrogen-bond
 * pair and A-T only two. This is the Wallace rule, valid for short oligos.
 */
export function meltingPoint(seq) {
  let gc = 0;
  for (const b of seq) if (b === 'G' || b === 'C') gc++;
  const at = seq.length - gc;
  if (seq.length === 0) return 0;
  if (seq.length < 14) return 2 * at + 4 * gc;                       // Wallace
  return 64.9 + 41 * (gc - 16.4) / seq.length;                       // Marmur
}

export function gcContent(seq) {
  if (!seq.length) return 0;
  let gc = 0;
  for (const b of seq) if (b === 'G' || b === 'C') gc++;
  return gc / seq.length;
}

/* ---------------------------------------------------------- genetic code */

/**
 * The standard genetic code, keyed by RNA codon. Sixty-four codons, twenty
 * amino acids, three stops. The redundancy is nearly all in the third
 * position — that is the wobble, and it is why many point mutations are
 * silent.
 */
export const CODON_TABLE = {
  UUU: 'F', UUC: 'F', UUA: 'L', UUG: 'L',
  CUU: 'L', CUC: 'L', CUA: 'L', CUG: 'L',
  AUU: 'I', AUC: 'I', AUA: 'I', AUG: 'M',
  GUU: 'V', GUC: 'V', GUA: 'V', GUG: 'V',
  UCU: 'S', UCC: 'S', UCA: 'S', UCG: 'S',
  CCU: 'P', CCC: 'P', CCA: 'P', CCG: 'P',
  ACU: 'T', ACC: 'T', ACA: 'T', ACG: 'T',
  GCU: 'A', GCC: 'A', GCA: 'A', GCG: 'A',
  UAU: 'Y', UAC: 'Y', UAA: '*', UAG: '*',
  CAU: 'H', CAC: 'H', CAA: 'Q', CAG: 'Q',
  AAU: 'N', AAC: 'N', AAA: 'K', AAG: 'K',
  GAU: 'D', GAC: 'D', GAA: 'E', GAG: 'E',
  UGU: 'C', UGC: 'C', UGA: '*', UGG: 'W',
  CGU: 'R', CGC: 'R', CGA: 'R', CGG: 'R',
  AGU: 'S', AGC: 'S', AGA: 'R', AGG: 'R',
  GGU: 'G', GGC: 'G', GGA: 'G', GGG: 'G',
};

/**
 * The twenty proteinogenic amino acids. `class` is the property that decides
 * how the residue behaves when the chain folds: hydrophobic residues bury
 * themselves in the core, charged ones stay on the surface. That is the whole
 * basis of the folding heuristic in fold.js.
 */
export const AMINO = {
  A: { name: 'Alanine',       tla: 'Ala', class: 'hydrophobic', colour: 0x9aa7b8 },
  R: { name: 'Arginine',      tla: 'Arg', class: 'basic',       colour: 0x4d7bff },
  N: { name: 'Asparagine',    tla: 'Asn', class: 'polar',       colour: 0x37c9a6 },
  D: { name: 'Aspartate',     tla: 'Asp', class: 'acidic',      colour: 0xff5a6e },
  C: { name: 'Cysteine',      tla: 'Cys', class: 'polar',       colour: 0xf5d76e },
  Q: { name: 'Glutamine',     tla: 'Gln', class: 'polar',       colour: 0x37c9a6 },
  E: { name: 'Glutamate',     tla: 'Glu', class: 'acidic',      colour: 0xff5a6e },
  G: { name: 'Glycine',       tla: 'Gly', class: 'hydrophobic', colour: 0xd7dee8 },
  H: { name: 'Histidine',     tla: 'His', class: 'basic',       colour: 0x6f8cff },
  I: { name: 'Isoleucine',    tla: 'Ile', class: 'hydrophobic', colour: 0x8f9db0 },
  L: { name: 'Leucine',       tla: 'Leu', class: 'hydrophobic', colour: 0x8f9db0 },
  K: { name: 'Lysine',        tla: 'Lys', class: 'basic',       colour: 0x4d7bff },
  M: { name: 'Methionine',    tla: 'Met', class: 'hydrophobic', colour: 0xc9b45e },
  F: { name: 'Phenylalanine', tla: 'Phe', class: 'aromatic',    colour: 0xb07cff },
  P: { name: 'Proline',       tla: 'Pro', class: 'hydrophobic', colour: 0x7f8a99 },
  S: { name: 'Serine',        tla: 'Ser', class: 'polar',       colour: 0x37c9a6 },
  T: { name: 'Threonine',     tla: 'Thr', class: 'polar',       colour: 0x37c9a6 },
  W: { name: 'Tryptophan',    tla: 'Trp', class: 'aromatic',    colour: 0xb07cff },
  Y: { name: 'Tyrosine',      tla: 'Tyr', class: 'aromatic',    colour: 0xa06fe8 },
  V: { name: 'Valine',        tla: 'Val', class: 'hydrophobic', colour: 0x8f9db0 },
  '*': { name: 'Stop',        tla: 'Ter', class: 'stop',        colour: 0xff3b52 },
};

/** Side-chain bulk, in cubic angstroms. Used to size the residue spheres. */
export const RESIDUE_VOLUME = {
  G: 60, A: 89, S: 89, C: 109, D: 111, P: 113, N: 114, T: 116, E: 138,
  V: 140, Q: 144, H: 153, M: 163, I: 167, L: 167, K: 169, R: 174, F: 190,
  Y: 194, W: 228, '*': 100,
};

/**
 * Chou-Fasman helix propensity, P(alpha). Above 1.0 the residue is found in
 * alpha helices more often than chance; below, less.
 *
 * The two at the bottom of the list are the interesting ones. Glycine has no
 * side chain, so the backbone can go anywhere and there is no reason to pick
 * the helix. Proline's side chain loops back and bonds to its own backbone
 * nitrogen, so it has no N-H to donate and it kinks the chain — a proline in
 * the middle of a helix breaks it. Both come out at 0.57, and in the folding
 * here that is what decides where a helix stops.
 */
export const HELIX_PROPENSITY = {
  E: 1.51, M: 1.45, A: 1.42, L: 1.21, K: 1.16, F: 1.13, Q: 1.11, W: 1.08,
  I: 1.08, V: 1.06, D: 1.01, H: 1.00, R: 0.98, T: 0.83, S: 0.77, C: 0.70,
  Y: 0.69, N: 0.67, P: 0.57, G: 0.57, '*': 0.5,
};

/**
 * Alpha-carbon distances that define an alpha helix, in angstroms. A chain
 * holding all three at once has nowhere to go but 3.6 residues per turn with a
 * 1.5 A rise — the helix is not imposed, it is what these distances mean.
 */
export const HELIX_SPACING = {
  i3: 5.0,   // Ca(i) to Ca(i+3)
  i4: 6.2,   // Ca(i) to Ca(i+4) — the one the backbone hydrogen bond makes
};

/** Kyte-Doolittle hydropathy. Positive is water-hating. */
export const HYDROPATHY = {
  I: 4.5, V: 4.2, L: 3.8, F: 2.8, C: 2.5, M: 1.9, A: 1.8, G: -0.4, T: -0.7,
  S: -0.8, W: -0.9, Y: -1.3, P: -1.6, H: -3.2, E: -3.5, Q: -3.5, D: -3.5,
  N: -3.5, K: -3.9, R: -4.5, '*': 0,
};

/** Transcribe a DNA coding strand to mRNA: T becomes U, nothing else moves. */
export function transcribe(dna) {
  return dna.replace(/T/g, 'U');
}

/** Read an mRNA into codons from a given offset. */
export function codonsOf(rna, offset = 0) {
  const out = [];
  for (let i = offset; i + 3 <= rna.length; i += 3) out.push(rna.slice(i, i + 3));
  return out;
}

/**
 * Translate mRNA. Ribosomes begin at the first AUG and run until a stop codon,
 * so that is what this does unless told to read straight through.
 */
export function translate(rna, { requireStart = true } = {}) {
  const start = requireStart ? rna.indexOf('AUG') : 0;
  if (start < 0) return { start: -1, peptide: [], stopped: false, stopCodon: null };

  const peptide = [];
  let i = start;
  for (; i + 3 <= rna.length; i += 3) {
    const codon = rna.slice(i, i + 3);
    const aa = CODON_TABLE[codon];
    if (aa === undefined) break;
    if (aa === '*') return { start, peptide, stopped: true, stopCodon: codon, end: i + 3 };
    peptide.push({ codon, aa, index: peptide.length, rnaIndex: i });
  }
  return { start, peptide, stopped: false, stopCodon: null, end: i };
}

/**
 * Classify a single-base substitution by what it does to the protein. This is
 * the payoff of the whole pipeline: it is why a one-letter typo in a sequence
 * can mean nothing at all, or can mean sickle-cell.
 */
export function classifyMutation(dnaBefore, dnaAfter) {
  const a = translate(transcribe(dnaBefore));
  const b = translate(transcribe(dnaAfter));

  if (a.start < 0 || b.start < 0) return { type: 'no-orf', label: 'no reading frame' };

  const pa = a.peptide.map((p) => p.aa).join('');
  const pb = b.peptide.map((p) => p.aa).join('');

  if (pa === pb) return { type: 'silent', label: 'silent — same protein', from: null, to: null };
  if (pb.length < pa.length && b.stopped) {
    return { type: 'nonsense', label: `nonsense — truncated at residue ${pb.length + 1}`, at: pb.length };
  }
  if (pa.length === pb.length) {
    for (let i = 0; i < pa.length; i++) {
      if (pa[i] !== pb[i]) {
        return {
          type: 'missense',
          label: `missense — ${AMINO[pa[i]].tla}${i + 1}${AMINO[pb[i]].tla}`,
          at: i, from: pa[i], to: pb[i],
        };
      }
    }
  }
  return { type: 'frameshift', label: 'reading frame changed', at: 0 };
}

/* --------------------------------------------------------------- sequences */

/**
 * Real sequences, so the thing on screen is not gibberish. Each is short
 * enough to see end to end and carries an intact open reading frame.
 */
export const PRESETS = [
  {
    id: 'insulin',
    name: 'Human insulin, B-chain start',
    note: 'The first residues of the insulin B chain — Phe-Val-Asn-Gln-His-Leu.',
    seq: 'ATGTTTGTGAACCAACATCTGTGCGGCTCACACCTGGTGGAAGCTCTCTACCTAGTGTGCGGGTAA',
  },
  {
    id: 'sickle',
    name: 'Beta-globin, start of the gene',
    // The real thing: this translates to MVHLTPEEKSAVTALWGK, the N-terminus of
    // human haemoglobin beta. Clicking position 20 from A to T turns GAG into
    // GTG and that single letter is sickle-cell disease.
    //
    // The readout will call it Glu7Val while every textbook calls it Glu6Val.
    // Both are right: haemoglobin is numbered from the mature protein, after
    // the initiator methionine has been cut off, and this counts the chain as
    // the ribosome actually built it.
    note: 'Click position 20 (A→T): GAG becomes GTG, and that one letter is sickle-cell.',
    seq: 'ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGTAA',
  },
  {
    id: 'villin',
    name: 'Villin headpiece HP36',
    // MLSDEDFKAVFGMTRSAFANLPLWKQQNLKKEKGLF — 36 residues, three short helices,
    // and the standard benchmark for folding simulations because it is about
    // the smallest thing that folds into a real protein rather than a blob.
    // Long enough to have an inside, which the other presets are not.
    note: 'Thirty-six residues that genuinely fold. Use the Proteins view on this one.',
    seq: 'ATGCTGAGCGATGAAGATTTTAAAGCGGTGTTTGGCATGACCCGTAGCGCGTTTGCGAACCTGCCGCTGTGGAAACAGCAGAACCTGAAAAAAGAAAAAGGCCTGTTTTAA',
  },
  {
    id: 'tata',
    name: 'TATA box promoter',
    note: 'The AT-rich block a polymerase melts first, because A-T has only two hydrogen bonds.',
    seq: 'GCGCGCTATAAAAGGGCGCGCATGAAACGCGAAGCTTTTCGCGCGCTAA',
  },
  {
    id: 'gcrich',
    name: 'GC-rich island',
    note: 'Three hydrogen bonds per rung. Hard to melt — watch the helicase struggle.',
    seq: 'ATGGCGGCGCCGCCGGGCGCGGGCCCGCGCGGCGCCGGCTGCTAA',
  },
];

export const DEFAULT_SEQUENCE = PRESETS[0].seq;
