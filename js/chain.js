/**
 * chain.js — a nucleic acid strand as data, plus the B-DNA placement maths.
 *
 * A Chain is a flat array of residue slots. Each slot carries a letter, a
 * present flag, and the three vectors a renderer needs: where the sugar's C1'
 * sits, which way the base points from it, and which way the helix axis runs
 * there. Nothing here knows about THREE beyond Vector3 arithmetic on plain
 * typed arrays, and nothing here knows what process is driving it — the
 * processes write these slots, the renderer reads them.
 */

import { BDNA, BASES } from './bio.js';

const DEG = Math.PI / 180;

/* ------------------------------------------------------------------ helices */

/**
 * Place one residue on a B-DNA helix.
 *
 * The helix runs along +z through the point (ax, ay). `phase` is the angle of
 * this strand's glycosidic bond about that axis; the partner strand on the
 * same helix is offset by MINOR_GROOVE_SPAN.
 *
 * Writes into `out`, which is a plain object of arrays so callers can pour the
 * result straight into a Chain slot without allocating.
 */
export function helixResidue(out, i, {
  ax = 0, ay = 0, z0 = 0,
  phaseOffset = 0,
  handed = 1,          // +1 right-handed, which is what B-DNA is
  radiusScale = 1,
  rise = BDNA.RISE,
  twist = BDNA.TWIST,
} = {}) {
  const phase = (phaseOffset + i * twist * handed) * DEG;
  const z = z0 + i * rise;

  const cp = Math.cos(phase);
  const sp = Math.sin(phase);

  const rC1 = BDNA.R_GLYCOSIDIC * radiusScale;
  const rSug = BDNA.R_SUGAR * radiusScale;
  const rP = BDNA.R_PHOSPHATE * radiusScale;

  out.pos = [ax + cp * rC1, ay + sp * rC1, z];
  out.sugar = [ax + cp * rSug, ay + sp * rSug, z];
  out.backbone = [ax + cp * rP, ay + sp * rP, z];
  // The base points inward, at the axis — that is the whole architecture:
  // hydrophobic bases stacked in the core, charged backbone out in the water.
  out.out = [-cp, -sp, 0];
  out.up = [0, 0, 1];
  out.phase = phase;
  return out;
}

/* -------------------------------------------------------------------- Chain */

const SCRATCH = { pos: null, sugar: null, backbone: null, out: null, up: null, phase: 0 };

export class Chain {
  /**
   * @param {number} capacity  maximum residues
   * @param {object} opts
   * @param {'dna'|'rna'} opts.kind
   * @param {string} opts.label  shown in the UI when this strand is highlighted
   */
  constructor(capacity, { kind = 'dna', label = '', role = '' } = {}) {
    this.capacity = capacity;
    this.kind = kind;
    this.label = label;
    this.role = role;      // 'template' | 'leading' | 'lagging' | 'mrna' | ''
    this.length = 0;

    /**
     * Whether successive residues are joined by phosphodiester bonds. A pool
     * of loose nucleotides waiting in the nucleoplasm is a Chain with this
     * turned off: same residues, no backbone drawn, because there isn't one.
     */
    this.covalent = true;

    this.letters = new Array(capacity).fill('');
    this.present = new Uint8Array(capacity);
    this.highlight = new Float32Array(capacity);   // 0..1, glow amount

    this.pos = new Float32Array(capacity * 3);
    this.sugar = new Float32Array(capacity * 3);
    this.backbone = new Float32Array(capacity * 3);
    this.out = new Float32Array(capacity * 3);
    this.up = new Float32Array(capacity * 3);

    /** Extra per-residue scale, used to fade a nucleotide in as it docks. */
    this.scale = new Float32Array(capacity).fill(1);
  }

  clear() {
    this.present.fill(0);
    this.highlight.fill(0);
    this.scale.fill(1);
    this.length = 0;
  }

  setLetters(seq) {
    this.length = Math.min(seq.length, this.capacity);
    for (let i = 0; i < this.length; i++) this.letters[i] = seq[i];
  }

  /** Copy a computed residue frame into slot i. */
  write(i, f) {
    const o = i * 3;
    this.pos[o] = f.pos[0]; this.pos[o + 1] = f.pos[1]; this.pos[o + 2] = f.pos[2];
    this.sugar[o] = f.sugar[0]; this.sugar[o + 1] = f.sugar[1]; this.sugar[o + 2] = f.sugar[2];
    this.backbone[o] = f.backbone[0]; this.backbone[o + 1] = f.backbone[1]; this.backbone[o + 2] = f.backbone[2];
    this.out[o] = f.out[0]; this.out[o + 1] = f.out[1]; this.out[o + 2] = f.out[2];
    this.up[o] = f.up[0]; this.up[o + 1] = f.up[1]; this.up[o + 2] = f.up[2];
  }

  /** Place slot i on a helix described by `opts`. */
  placeHelical(i, opts) {
    this.write(i, helixResidue(SCRATCH, i, opts));
  }

  /** Place slot i explicitly, for strands that are not on a helix at all. */
  placeFree(i, pos, out, up) {
    const o = i * 3;
    this.pos[o] = pos[0]; this.pos[o + 1] = pos[1]; this.pos[o + 2] = pos[2];
    this.out[o] = out[0]; this.out[o + 1] = out[1]; this.out[o + 2] = out[2];
    this.up[o] = up[0]; this.up[o + 1] = up[1]; this.up[o + 2] = up[2];
    // Backbone sits behind the base, sugar between the two.
    const bx = pos[0] - out[0] * BDNA.R_PHOSPHATE * 0.34;
    const by = pos[1] - out[1] * BDNA.R_PHOSPHATE * 0.34;
    const bz = pos[2] - out[2] * BDNA.R_PHOSPHATE * 0.34;
    this.backbone[o] = bx; this.backbone[o + 1] = by; this.backbone[o + 2] = bz;
    this.sugar[o] = (pos[0] + bx) / 2;
    this.sugar[o + 1] = (pos[1] + by) / 2;
    this.sugar[o + 2] = (pos[2] + bz) / 2;
  }

  getPos(i, target) {
    const o = i * 3;
    target.set(this.pos[o], this.pos[o + 1], this.pos[o + 2]);
    return target;
  }

  getBackbone(i, target) {
    const o = i * 3;
    target.set(this.backbone[o], this.backbone[o + 1], this.backbone[o + 2]);
    return target;
  }

  /**
   * Contiguous runs of present residues. Okazaki fragments are exactly this:
   * a lagging strand is one Chain holding several unjoined runs, and ligase
   * finishing the job is those runs merging into one.
   */
  runs() {
    const out = [];
    let start = -1;
    for (let i = 0; i < this.length; i++) {
      if (this.present[i]) {
        if (start < 0) start = i;
      } else if (start >= 0) {
        out.push([start, i - 1]);
        start = -1;
      }
    }
    if (start >= 0) out.push([start, this.length - 1]);
    return out;
  }

  countPresent() {
    let n = 0;
    for (let i = 0; i < this.length; i++) n += this.present[i];
    return n;
  }
}

/* ------------------------------------------------------- duplex convenience */

/**
 * Lay two chains out as an intact double helix.
 *
 * `a` is written 5'->3' with increasing index. `b` is antiparallel, so index i
 * of chain b is the partner of index i of chain a but chain b's own 5' end is
 * at high i. The renderer does not care; the processes do, and this is where
 * that convention is fixed.
 */
export function layDuplex(a, b, n, opts = {}) {
  const {
    ax = 0, ay = 0, z0 = 0, phaseOffset = 0, radiusScale = 1,
    twist = BDNA.TWIST, rise = BDNA.RISE,
  } = opts;

  for (let i = 0; i < n; i++) {
    a.placeHelical(i, { ax, ay, z0, phaseOffset, radiusScale, twist, rise });
    b.placeHelical(i, {
      ax, ay, z0, phaseOffset: phaseOffset + BDNA.MINOR_GROOVE_SPAN, radiusScale, twist, rise,
    });
  }
}

/**
 * Point two paired bases at each other.
 *
 * `placeHelical` aims each base at the helix axis, which is right for a strand
 * with no partner but wrong for a pair: the two bases of a Watson-Crick pair
 * are coplanar and their long axes lie along the line between the two C1'
 * atoms, not along two radii 127 deg apart. Aiming them radially leaves the
 * hydrogen bonds nearly 4 A long and visibly skew.
 *
 * Also applies propeller twist. The two bases of a pair are not quite
 * coplanar — they turn against each other by about 11 deg, like a propeller,
 * which is what lets the stack stay tight while the backbone spirals.
 */
export function facePair(a, i, b, j, propeller = BDNA.PROPELLER) {
  const oa = i * 3, ob = j * 3;
  let dx = b.pos[ob] - a.pos[oa];
  let dy = b.pos[ob + 1] - a.pos[oa + 1];
  let dz = b.pos[ob + 2] - a.pos[oa + 2];
  const l = Math.hypot(dx, dy, dz);
  if (l < 1e-6) return;
  dx /= l; dy /= l; dz /= l;

  a.out[oa] = dx; a.out[oa + 1] = dy; a.out[oa + 2] = dz;
  b.out[ob] = -dx; b.out[ob + 1] = -dy; b.out[ob + 2] = -dz;

  // Turn each base's plane by half the propeller angle, in opposite senses,
  // about its own long axis.
  const half = (propeller / 2) * Math.PI / 180;
  const roll = (chain, o, sign) => {
    const ux = chain.up[o], uy = chain.up[o + 1], uz = chain.up[o + 2];
    const k = sign * half;
    const c = Math.cos(k), s = Math.sin(k);
    // Rodrigues, about the (already unit) long axis d.
    const dot = dx * ux + dy * uy + dz * uz;
    const cx = dy * uz - dz * uy, cy = dz * ux - dx * uz, cz = dx * uy - dy * ux;
    chain.up[o] = ux * c + cx * s + dx * dot * (1 - c);
    chain.up[o + 1] = uy * c + cy * s + dy * dot * (1 - c);
    chain.up[o + 2] = uz * c + cz * s + dz * dot * (1 - c);
  };
  roll(a, oa, +1);
  roll(b, ob, -1);
}

/**
 * Where the hydrogen bonds of a pair should be drawn: between the far tips of
 * the two base plates, two bonds for A-T and three for G-C, spread across the
 * width of the pair.
 */
export function hydrogenBonds(chainA, i, chainB, j, count, target) {
  const oa = i * 3, ob = j * 3;

  // Each bond starts at the tip of its own base plate, so a purine's bonds
  // start further in than a pyrimidine's and the gap between them lands at the
  // real 2.8-3.0 A. Nothing here assumes a bond length; it is what is left over.
  const ra = (BASES[chainA.letters[i]] || { reach: 3.7 }).reach;
  const rb = (BASES[chainB.letters[j]] || { reach: 3.7 }).reach;

  const pax = chainA.pos[oa] + chainA.out[oa] * ra;
  const pay = chainA.pos[oa + 1] + chainA.out[oa + 1] * ra;
  const paz = chainA.pos[oa + 2] + chainA.out[oa + 2] * ra;
  const pbx = chainB.pos[ob] + chainB.out[ob] * rb;
  const pby = chainB.pos[ob + 1] + chainB.out[ob + 1] * rb;
  const pbz = chainB.pos[ob + 2] + chainB.out[ob + 2] * rb;

  // Spread the bonds along the pair's short axis, which is the helix axis
  // crossed with the long axis.
  let lx = pbx - pax, ly = pby - pay, lz = pbz - paz;
  const ll = Math.hypot(lx, ly, lz) || 1;
  lx /= ll; ly /= ll; lz /= ll;
  const ux = chainA.up[oa], uy = chainA.up[oa + 1], uz = chainA.up[oa + 2];
  let sx = ly * uz - lz * uy, sy = lz * ux - lx * uz, sz = lx * uy - ly * ux;
  const sl = Math.hypot(sx, sy, sz) || 1;
  sx /= sl; sy /= sl; sz /= sl;

  target.length = 0;
  const spread = 1.15;
  for (let k = 0; k < count; k++) {
    const t = count === 1 ? 0 : (k / (count - 1) - 0.5) * 2 * spread;
    target.push({
      from: [pax + sx * t, pay + sy * t, paz + sz * t],
      to: [pbx + sx * t, pby + sy * t, pbz + sz * t],
    });
  }
  return target;
}
