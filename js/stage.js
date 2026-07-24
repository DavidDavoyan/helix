/**
 * stage.js — the shared world, and the conventions every scene obeys.
 *
 * Conventions, fixed here once:
 *
 *   - `seq` is the coding strand, written 5'->3', index 0 at the 5' end.
 *   - Base pair i sits at z = (i - (n-1)/2) * RISE, so the molecule is centred
 *     on the origin and grows along +z.
 *   - Chain `a` is the coding strand. Chain `b` is the template, antiparallel:
 *     slot i of b is the partner of slot i of a, but b's own 5' end is at high
 *     i. Every scene that cares about direction reads that from here.
 *   - +i is "downstream". A replication fork travels towards -i; RNA
 *     polymerase travels towards +i.
 */

import * as THREE from './lib/three.module.js';
import { BDNA, complement } from './bio.js';
import { Chain, helixResidue, facePair } from './chain.js';
import { NucleicRenderer } from './render.js';
import { NucleotidePool } from './pool.js';
import { Peptide } from './peptide.js';
import {
  makeHelicase, makePolymerase, makeRnaPolymerase, makeRibosome,
  makeSSB, makePrimase, makeLigase, makeLabel,
} from './machines.js';

export const MAX_BP = 260;

const SCRATCH = {};

/** z of base pair i for a molecule of n pairs. */
export function zOf(i, n) {
  return (i - (n - 1) / 2) * BDNA.RISE;
}

/** 0 below edge0, 1 above edge1, smooth between. */
export function smoothstep(edge0, edge1, x) {
  const t = Math.min(1, Math.max(0, (x - edge0) / (edge1 - edge0 || 1e-6)));
  return t * t * (3 - 2 * t);
}

/**
 * Place one residue of a melted region — position only.
 *
 * When a base pair opens the backbone swings out to a wider radius and bulges
 * sideways. `open` runs 0 (paired) to 1 (fully separated), `side` is +1 or -1
 * for the two strands. The base's *direction* is not set here: pairs have to
 * be aimed at each other first (facePair) and then swung apart (swingOut), or
 * the orientation jumps the moment the bonds break.
 */
export function placeMelted(chain, i, {
  n, phaseOffset, open, side, splay = 11, z0 = null, ax = 0, ay = 0,
}) {
  helixResidue(SCRATCH, i, {
    ax, ay, phaseOffset,
    z0: z0 !== null ? z0 : zOf(0, n),
    radiusScale: 1 + open * 0.42,
  });

  const s = open * splay * side;
  SCRATCH.pos[1] += s;
  SCRATCH.sugar[1] += s;
  SCRATCH.backbone[1] += s * 1.08;

  chain.write(i, SCRATCH);
}

/**
 * Turn a base away from its lost partner and out into the water. At `open` = 1
 * it points straight out from the helix axis, which is where an unpaired base
 * ends up once there is nothing to stack against.
 */
export function swingOut(chain, i, open, { ax = 0, ay = 0 } = {}) {
  if (open <= 0) return;
  const o = i * 3;
  let tx = chain.pos[o] - ax;
  let ty = chain.pos[o + 1] - ay;
  const tl = Math.hypot(tx, ty) || 1;
  tx /= tl; ty /= tl;

  const ox = chain.out[o] * (1 - open) + tx * open;
  const oy = chain.out[o + 1] * (1 - open) + ty * open;
  const oz = chain.out[o + 2] * (1 - open);
  const l = Math.hypot(ox, oy, oz) || 1;
  chain.out[o] = ox / l; chain.out[o + 1] = oy / l; chain.out[o + 2] = oz / l;
}

/* -------------------------------------------------------------------- Stage */

export class Stage {
  constructor(scene) {
    this.root = new THREE.Group();
    scene.add(this.root);

    this.renderer = new NucleicRenderer(this.root, {
      maxResidues: MAX_BP * 4,
      maxBonds: MAX_BP * 8,
      maxRibbons: 48,
    });

    this.pool = new NucleotidePool({ count: 96, bounds: [52, 52, 120] });
    this.root.add(this.pool.sparks.points);

    this.peptide = new Peptide(this.root, { capacity: 240 });

    // The four chains every scene draws from. Not every scene uses all four.
    this.coding = new Chain(MAX_BP, { kind: 'dna', label: 'coding strand', role: 'coding' });
    this.template = new Chain(MAX_BP, { kind: 'dna', label: 'template strand', role: 'template' });
    this.newTop = new Chain(MAX_BP, { kind: 'dna', label: 'new strand', role: 'leading' });
    this.newBottom = new Chain(MAX_BP, { kind: 'dna', label: 'new strand', role: 'lagging' });
    this.rna = new Chain(MAX_BP + 20, { kind: 'rna', label: 'mRNA', role: 'mrna' });

    this.machines = {
      helicase: makeHelicase(),
      polyLeading: makePolymerase(0xffb454, 0.92),
      polyLagging: makePolymerase(0xff8a5c, 0.86),
      rnap: makeRnaPolymerase(),
      ribosome: makeRibosome(),
      primase: makePrimase(),
      ligase: makeLigase(),
      ssb: [makeSSB(), makeSSB(), makeSSB(), makeSSB()],
    };
    for (const k of Object.keys(this.machines)) {
      const m = this.machines[k];
      if (Array.isArray(m)) m.forEach((x) => this.root.add(x));
      else this.root.add(m);
    }

    this.labels = {};
    this.showLabels = true;

    this.seq = '';
    this.n = 0;
  }

  /** Attach a floating tag to a machine. Created lazily, reused after that. */
  label(key, text, target, offset = [0, 26, 0]) {
    let l = this.labels[key];
    if (!l || l.userData.text !== text) {
      if (l) this.root.remove(l);
      l = makeLabel(text);
      l.userData.text = text;
      this.labels[key] = l;
      this.root.add(l);
    }
    l.visible = this.showLabels && target.visible;
    if (target.visible) {
      l.position.set(
        target.position.x + offset[0],
        target.position.y + offset[1],
        target.position.z + offset[2],
      );
    }
    return l;
  }

  hideAllMachines() {
    for (const k of Object.keys(this.machines)) {
      const m = this.machines[k];
      if (Array.isArray(m)) m.forEach((x) => { x.visible = false; });
      else m.visible = false;
    }
    for (const k of Object.keys(this.labels)) this.labels[k].visible = false;
  }

  /** Load a sequence. The template is derived, never stored independently. */
  setSequence(seq) {
    this.seq = seq;
    this.n = Math.min(seq.length, MAX_BP);

    this.coding.clear();
    this.template.clear();
    this.newTop.clear();
    this.newBottom.clear();
    this.rna.clear();

    this.coding.length = this.n;
    this.template.length = this.n;
    this.newTop.length = this.n;
    this.newBottom.length = this.n;

    for (let i = 0; i < this.n; i++) {
      this.coding.letters[i] = seq[i];
      this.template.letters[i] = complement(seq[i]);
      this.coding.present[i] = 1;
      this.template.present[i] = 1;
    }
  }

  /** Lay the whole molecule out as an undisturbed B-DNA duplex. */
  layIntact(openFn = null) {
    const n = this.n;
    const z0 = zOf(0, n);
    for (let i = 0; i < n; i++) {
      const open = openFn ? openFn(i) : 0;
      if (open < 0.001) {
        this.coding.placeHelical(i, { z0, phaseOffset: 0 });
        this.template.placeHelical(i, { z0, phaseOffset: BDNA.MINOR_GROOVE_SPAN });
        facePair(this.coding, i, this.template, i);
      } else {
        placeMelted(this.coding, i, { n, phaseOffset: 0, open, side: +1, z0 });
        placeMelted(this.template, i, { n, phaseOffset: BDNA.MINOR_GROOVE_SPAN, open, side: -1, z0 });
        // Aim them at each other first, then swing them apart by however far
        // the pair has opened — so the turn is continuous as it melts.
        facePair(this.coding, i, this.template, i, BDNA.PROPELLER * (1 - open));
        swingOut(this.coding, i, open);
        swingOut(this.template, i, open);
      }
    }
  }

  dispose() {
    this.root.traverse((o) => {
      if (o.geometry) o.geometry.dispose();
      if (o.material) {
        if (Array.isArray(o.material)) o.material.forEach((m) => m.dispose());
        else o.material.dispose();
      }
    });
  }
}
