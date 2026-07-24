/**
 * scene-expression.js — reading the sequence out, and building the protein.
 */

import * as THREE from './lib/three.module.js';
import {
  BDNA, BASES, CODON_TABLE, AMINO, transcribe, translate, rnaComplement,
} from './bio.js';
import { Chain, hydrogenBonds, facePair } from './chain.js';
import { makeTRNA } from './machines.js';
import { zOf, smoothstep, swingOut } from './stage.js';

const _bonds = [];
const _v = new THREE.Vector3();
const _t = new THREE.Vector3();
const _n = new THREE.Vector3();
const _s = new THREE.Vector3();

const NT_SPACING = 6.3;   // A between consecutive bases in a single strand

/* ------------------------------------------------------------ free strands */

/**
 * A single strand with nothing holding it in shape.
 *
 * A transcript coming out of a polymerase is not a helix and is not straight;
 * it is a floppy chain being extruded into water. So this is a chain of points
 * held at a fixed spacing, pushed apart where they collide, and otherwise left
 * to drift. The newest residue is pinned wherever the machine put it.
 */
class FreeStrand {
  constructor(capacity) {
    this.capacity = capacity;
    this.count = 0;
    this.pos = [];
    this.vel = [];
    for (let i = 0; i < capacity; i++) {
      this.pos.push(new THREE.Vector3());
      this.vel.push(new THREE.Vector3());
    }
  }

  clear() { this.count = 0; }

  append(at) {
    if (this.count >= this.capacity) return;
    const i = this.count++;
    this.pos[i].set(at[0], at[1], at[2]);
    this.vel[i].set(0, 0, 0);
  }

  /** `pinned` residues at the head are held fast by whatever is making them. */
  relax(dt, { spacing = NT_SPACING, pinned = 1, drift = 0 } = {}) {
    const n = this.count;
    if (n < 2) return;
    const step = Math.min(dt, 1 / 30);

    for (let i = 0; i < n - pinned; i++) {
      const p = this.pos[i];
      const v = this.vel[i];
      for (const j of [i - 1, i + 1]) {
        if (j < 0 || j >= n) continue;
        _v.copy(this.pos[j]).sub(p);
        const d = _v.length() || 1e-4;
        v.addScaledVector(_v.multiplyScalar(1 / d), (d - spacing) * 55 * step);
      }
      // Self-avoidance, so a long transcript does not fold through itself.
      for (let j = 0; j < n; j++) {
        if (Math.abs(j - i) < 3) continue;
        _v.copy(this.pos[j]).sub(p);
        const d2 = _v.lengthSq();
        if (d2 > 90 || d2 < 1e-6) continue;
        const d = Math.sqrt(d2);
        v.addScaledVector(_v.multiplyScalar(1 / d), -(9.5 - d) * 16 * step);
      }
      if (drift) {
        v.x += (Math.random() - 0.5) * drift * step;
        v.y += (Math.random() - 0.5) * drift * step;
        v.z += (Math.random() - 0.5) * drift * step;
      }
      v.multiplyScalar(Math.exp(-5.5 * step));
      if (v.lengthSq() > 2500) v.setLength(50);
      p.addScaledVector(v, step);
    }
  }

  /**
   * Pour the strand into a Chain. The base plate's normal follows the chain's
   * tangent — bases stack perpendicular to the backbone — and the direction
   * the base points is carried along by parallel transport so it turns
   * smoothly instead of flipping.
   */
  writeTo(chain, twistPerResidue = 32) {
    const n = this.count;
    chain.length = n;
    const ref = _s.set(0, 1, 0);
    let carry = null;

    for (let i = 0; i < n; i++) {
      // A one-residue transcript has no tangent to take; give it an arbitrary
      // one rather than reading off the end of the array.
      if (n < 2) _t.set(0, 0, 1);
      else if (i === 0) _t.copy(this.pos[1]).sub(this.pos[0]);
      else if (i === n - 1) _t.copy(this.pos[i]).sub(this.pos[i - 1]);
      else _t.copy(this.pos[i + 1]).sub(this.pos[i - 1]);
      if (_t.lengthSq() < 1e-10) _t.set(0, 0, 1);
      _t.normalize();

      if (!carry) {
        carry = new THREE.Vector3().copy(ref);
        if (Math.abs(carry.dot(_t)) > 0.9) carry.set(1, 0, 0);
      }
      _n.copy(carry).projectOnPlane(_t);
      if (_n.lengthSq() < 1e-8) _n.set(1, 0, 0).projectOnPlane(_t);
      _n.normalize();
      // Twist the base direction round the backbone as we go, the way a
      // stacked single strand actually spirals.
      _n.applyAxisAngle(_t, twistPerResidue * Math.PI / 180);
      carry.copy(_n);

      chain.placeFree(i,
        [this.pos[i].x, this.pos[i].y, this.pos[i].z],
        [_n.x, _n.y, _n.z],
        [_t.x, _t.y, _t.z]);
      chain.present[i] = 1;
    }
  }
}

/* ======================================================= TRANSCRIPTION ===== */

const BUBBLE = 7;   // half-width of the melted region, in base pairs

/**
 * Transcription.
 *
 * RNA polymerase does not need a helicase and does not need a primer. It opens
 * about fourteen base pairs itself, copies one strand, and lets the DNA close
 * again behind it, so the bubble travels along the molecule as a moving
 * blister. The transcript leaves through a channel of its own — if it stayed
 * paired to the template the DNA could never re-form.
 *
 * The copy is made from the template strand, which means the transcript comes
 * out reading like the *coding* strand, with U wherever that had T.
 */
export class TranscriptionScene {
  constructor(stage) {
    this.stage = stage;
    this.id = 'transcription';
    this.title = 'Transcription';
    this.strand = new FreeStrand(stage.rna.capacity);
    this.reset();
  }

  reset() {
    this.pos = -BUBBLE;      // polymerase position, in base pairs
    this.job = null;
    this.made = 0;
    this.done = false;
    this.strand.clear();
    this._t = 0;
  }

  enter() {
    const s = this.stage;
    s.hideAllMachines();
    s.pool.reset();
    s.pool.setVisible(true);
    s.peptide.clear();
    for (let i = 0; i < s.n; i++) {
      s.coding.present[i] = 1;
      s.template.present[i] = 1;
    }
    s.newTop.present.fill(0);
    s.newBottom.present.fill(0);
    s.rna.clear();
    this.reset();
  }

  exit() {}

  /** 1 inside the bubble, 0 outside, with soft edges. */
  _open(i) {
    const d = Math.abs(i - this.pos);
    return 1 - smoothstep(BUBBLE - 2.5, BUBBLE + 1.5, d);
  }

  /**
   * Where the transcript leaves the polymerase, in world coordinates. Well
   * clear of the body, or the RNA is born inside the protein and never seen.
   */
  _exitPoint() {
    const s = this.stage;
    const z0 = zOf(0, s.n);
    return [-40, -29, z0 + this.pos * BDNA.RISE - 10];
  }

  update(dt, speed) {
    const s = this.stage;
    const n = s.n;
    this._t += dt;
    const z0 = zOf(0, n);

    /* --- polymerase walks the template ---------------------------------- */
    if (!this.done) {
      this.pos += dt * speed * 3.0;
      if (this.pos > n - 1 + BUBBLE) { this.pos = n - 1 + BUBBLE; this.done = true; }
    }

    /* --- and lays down one nucleotide per base pair it passes ------------ */
    const wantIndex = this.made;
    if (!this.job && wantIndex < n && wantIndex < this.pos) {
      const letter = rnaComplement(s.template.letters[wantIndex]);
      const exit = this._exitPoint();
      this.job = s.pool.incorporate({
        letter: { A: 'A', U: 'T', G: 'G', C: 'C' }[letter] || 'A',
        site: () => ({ pos: this._exitPoint(), out: [0, -1, 0], up: [0, 0, 1] }),
        duration: 0.24 / Math.max(0.25, speed),
        wrongTries: Math.random() < 0.45 ? 1 : 0,
        onLock: () => {
          this.strand.append(exit);
          s.rna.letters[this.made] = letter;
          s.rna.highlight[this.made] = 1;
          this.made++;
        },
      });
    }
    if (this.job && this.job.done) this.job = null;

    /* --- geometry: a travelling bubble ---------------------------------- */
    for (let i = 0; i < n; i++) {
      const o = this._open(i);
      if (o < 0.02) {
        s.coding.placeHelical(i, { z0, phaseOffset: 0 });
        s.template.placeHelical(i, { z0, phaseOffset: BDNA.MINOR_GROOVE_SPAN });
        facePair(s.coding, i, s.template, i);
      } else {
        // The coding strand is pushed aside; the template is held open and read.
        s.coding.placeHelical(i, { z0, phaseOffset: 0, radiusScale: 1 + o * 0.55, ay: o * 8 });
        s.template.placeHelical(i, { z0, phaseOffset: BDNA.MINOR_GROOVE_SPAN, radiusScale: 1 + o * 0.3, ay: -o * 6 });
        facePair(s.coding, i, s.template, i, BDNA.PROPELLER * (1 - o));
        swingOut(s.coding, i, o, { ay: o * 8 });
        swingOut(s.template, i, o, { ay: -o * 6 });
      }
      s.coding.highlight[i] = 0;
      s.template.highlight[i] = o > 0.5 ? 0.35 : 0;
    }

    // The transcript hangs off the exit channel and drifts.
    if (this.strand.count) {
      const exit = this._exitPoint();
      this.strand.pos[this.strand.count - 1].set(exit[0], exit[1], exit[2]);
      this.strand.relax(dt, { pinned: 1, drift: 140 });
      this.strand.writeTo(s.rna);
      for (let i = 0; i < s.rna.length; i++) s.rna.highlight[i] *= Math.exp(-dt * 2.2);
    }

    /* --- machines -------------------------------------------------------- */
    const m = s.machines;
    m.rnap.visible = true;
    m.rnap.place([0, 0, z0 + this.pos * BDNA.RISE], [0, 0, 1]);
    m.rnap.rotation.z = Math.sin(this._t * 1.1) * 0.12;
    s.label('rnap', 'RNA polymerase', m.rnap, [0, 34, 0]);
    s.label('rna', 'mRNA transcript', { visible: this.made > 2, position: new THREE.Vector3(-34, -22, z0 + this.pos * BDNA.RISE) }, [0, 0, 0]);
  }

  draw(r) {
    const s = this.stage;
    r.addChain(s.coding, { backboneColour: 0x93a9c6 });
    r.addChain(s.template, { backboneColour: 0x6f8299 });
    r.addChain(s.rna, { backboneColour: 0xc084fc, radius: 1.2 });

    for (let i = 0; i < s.n; i++) {
      const o = this._open(i);
      if (o > 0.45) continue;
      const info = BASES[s.coding.letters[i]];
      if (!info) continue;
      hydrogenBonds(s.coding, i, s.template, i, info.hbonds, _bonds);
      const fade = 1 - o / 0.45;
      for (const b of _bonds) r.addBond(b.from, b.to, { radius: 0.22 * fade, colour: 0xdfe9ff });
    }
  }

  state() {
    const s = this.stage;
    const rna = s.rna.letters.slice(0, this.made).join('');
    return {
      caption: this.done
        ? 'Transcript complete. The DNA is untouched — this is a copy, and the original goes back in the drawer.'
        : 'The bubble travels: fourteen pairs open ahead of the polymerase and shut again behind it.',
      progress: this.made / Math.max(1, s.n),
      stats: [
        ['polymerase at', `${Math.max(0, Math.min(s.n, Math.round(this.pos) + 1))} / ${s.n}`],
        ['transcript', `${this.made} nt`],
        ['bubble', `${BUBBLE * 2} bp open`],
        ['reading', 'template strand, 3′→5′'],
      ],
      rna,
      track: {
        mark: [Math.max(0, Math.min(s.n - 1, Math.round(this.pos)))],
        doneRange: [0, Math.max(0, Math.round(this.pos) - 1)],
      },
    };
  }
}

/* ========================================================= TRANSLATION ===== */

const CODON_LEN = 3;

/**
 * Translation.
 *
 * The ribosome holds three tRNA sites in a row. A charged tRNA arrives at A;
 * if its anticodon pairs with the codon there, the chain hanging off the tRNA
 * in P is handed to it, and the whole ribosome then steps exactly three
 * nucleotides — so the tRNA that was in A is now in P and the empty one leaves
 * from E. Repeat until a codon arrives that no tRNA reads. That is a stop, and
 * the protein comes off.
 *
 * The three-at-a-time step is why the reading frame matters, and why deleting
 * a single base wrecks everything downstream of it.
 */
export class TranslationScene {
  constructor(stage) {
    this.stage = stage;
    this.id = 'translation';
    this.title = 'Translation';

    // A working set of tRNAs: one per site plus two for the rejects.
    this.trnas = [];
    for (let i = 0; i < 5; i++) {
      const t = makeTRNA('M', 'UAC');
      t.visible = false;
      stage.root.add(t);
      const anti = new Chain(3, { kind: 'rna', label: 'anticodon' });
      anti.covalent = false;
      anti.length = 3;
      this.trnas.push({ mesh: t, anticodon: anti, busy: false });
    }
    this.mrna = new Chain(stage.rna.capacity, { kind: 'rna', label: 'mRNA' });
    this.reset();
  }

  reset() {
    this.rna = '';
    this.reading = null;
    this.codon = 0;          // which codon of the ORF we are on
    this.phase = 'idle';
    this.phaseT = 0;
    this.occupancy = { A: null, P: null, E: null };
    this.rejects = [];
    this.done = false;
    this.released = false;
    this.foldT = 0;
    this.events = [];
    for (const t of this.trnas) { t.mesh.visible = false; t.busy = false; }
  }

  enter() {
    const s = this.stage;
    s.hideAllMachines();
    s.pool.reset();
    s.pool.setVisible(false);
    s.peptide.clear();
    s.coding.present.fill(0);
    s.template.present.fill(0);
    s.newTop.present.fill(0);
    s.newBottom.present.fill(0);
    s.rna.present.fill(0);
    this.reset();

    this.rna = transcribe(s.seq);
    this.reading = translate(this.rna);

    this.mrna.clear();
    this.mrna.length = this.rna.length;
    for (let i = 0; i < this.rna.length; i++) {
      this.mrna.letters[i] = this.rna[i];
      this.mrna.present[i] = 1;
    }

    if (this.reading.start < 0) {
      this.phase = 'noorf';
      this._log('No AUG in this sequence — nothing to translate.');
    } else {
      this.phase = 'scan';
      this.phaseT = 0;
      this._log(`Start codon AUG found at position ${this.reading.start + 1}`);
    }
  }

  exit() {
    for (const t of this.trnas) t.mesh.visible = false;
  }

  _log(text) {
    this.events.push(text);
    if (this.events.length > 5) this.events.shift();
  }

  /** World position of mRNA nucleotide k. */
  _nt(k) {
    const L = this.rna.length;
    return [0, 0, (k - (L - 1) / 2) * NT_SPACING];
  }

  /** World position of the ribosome when reading codon starting at nt `k`. */
  _ribosomeAt(k) {
    const p = this._nt(k + 1);
    return [p[0], p[1] + 18, p[2]];
  }

  _freeTrna() {
    return this.trnas.find((t) => !t.busy) || null;
  }

  _dressTrna(slot, aa, anticodon) {
    const info = AMINO[aa] || AMINO.G;
    slot.mesh.userData.aa = aa;
    slot.mesh.userData.anticodon = anticodon;
    slot.mesh.userData.residue.material.color.setHex(info.colour);
    slot.mesh.userData.residue.material.emissive.setHex(info.colour);
    slot.mesh.userData.residue.material.emissive.multiplyScalar(0.2);
    slot.mesh.userData.residue.visible = aa !== null;
    for (let i = 0; i < 3; i++) slot.anticodon.letters[i] = anticodon[i] || 'A';
    slot.anticodon.present.fill(1);
  }

  /**
   * Put a tRNA at a ribosome site, or somewhere on the way to one.
   *
   * The quarter turn about y is load-bearing. tRNA is modelled with its long
   * arm and its anticodon spread along local +x; the mRNA runs along world z.
   * Without the turn the L points at the camera and — worse — the three
   * anticodon bases lie across the message instead of along it, so they cannot
   * line up with a codon at all. Turning +90° maps local +x to world -z, which
   * both shows the L in profile and lays the anticodon antiparallel to the
   * codon, which is the direction it actually reads in.
   */
  _placeTrna(slot, worldPos, lean = 0, tilt = 0) {
    slot.mesh.visible = true;
    slot.mesh.position.set(worldPos[0], worldPos[1], worldPos[2]);
    slot.mesh.rotation.set(tilt, Math.PI / 2, lean);

    slot.mesh.updateMatrixWorld();
    for (let i = 0; i < 3; i++) {
      _v.set((i - 1) * NT_SPACING, -25.5, 0).applyMatrix4(slot.mesh.matrixWorld);
      slot.anticodon.placeFree(i, [_v.x, _v.y, _v.z], [0, -1, 0], [0, 0, 1]);
    }
  }

  update(dt, speed) {
    const s = this.stage;
    const m = s.machines;
    this.phaseT += dt * speed;

    const orf = this.reading;
    if (!orf || orf.start < 0) {
      m.ribosome.visible = false;
      return;
    }

    const codonStart = orf.start + this.codon * CODON_LEN;
    const rib = this._ribosomeAt(codonStart);

    m.ribosome.visible = true;
    m.ribosome.place(rib, [0, 0, 1]);
    m.ribosome.updateMatrixWorld();
    s.label('ribosome', 'ribosome', m.ribosome, [0, 52, 0]);

    const site = (name) => {
      const p = m.ribosome.userData.sites[name].clone();
      m.ribosome.localToWorld(p);
      return [p.x, p.y, p.z];
    };

    /* --- the elongation cycle -------------------------------------------- */
    const D = { scan: 0.7, testing: 0.45, docking: 0.5, bonding: 0.45, translocating: 0.55, release: 0.9 };

    switch (this.phase) {
      case 'scan': {
        if (this.phaseT >= D.scan) { this.phase = 'testing'; this.phaseT = 0; this._beginTesting(); }
        break;
      }
      case 'testing': {
        // One or two wrong tRNAs try the A site and fall off. Real ribosomes
        // reject far more than they accept; the fidelity is in the rejection.
        const k = Math.min(1, this.phaseT / D.testing);
        for (const r of this.rejects) {
          const from = r.from;
          const to = site('A');
          const e = k < 0.55 ? k / 0.55 : 1 - (k - 0.55) / 0.45;
          const t = e * 0.8;
          this._placeTrna(r.slot, [
            from[0] + (to[0] - from[0]) * t,
            from[1] + (to[1] - from[1]) * t,
            from[2] + (to[2] - from[2]) * t,
          ], (1 - e) * 0.7, (1 - e) * 0.5);
        }
        if (k >= 1) {
          for (const r of this.rejects) { r.slot.mesh.visible = false; r.slot.busy = false; }
          this.rejects = [];
          this.phase = 'docking';
          this.phaseT = 0;
          this._beginDocking(codonStart);
        }
        break;
      }
      case 'docking': {
        const k = Math.min(1, this.phaseT / D.docking);
        const e = k * k * (3 - 2 * k);
        const slot = this.occupancy.A;
        if (slot) {
          const to = site('A');
          this._placeTrna(slot, [
            slot.from[0] + (to[0] - slot.from[0]) * e,
            slot.from[1] + (to[1] - slot.from[1]) * e,
            slot.from[2] + (to[2] - slot.from[2]) * e,
          ], (1 - e) * 0.6, (1 - e) * 0.4);
        }
        if (k >= 1) { this.phase = 'bonding'; this.phaseT = 0; }
        break;
      }
      case 'bonding': {
        const slotA = this.occupancy.A;
        if (slotA) this._placeTrna(slotA, site('A'));
        if (this.occupancy.P) this._placeTrna(this.occupancy.P, site('P'));
        if (this.occupancy.E) this._placeTrna(this.occupancy.E, site('E'));

        if (this.phaseT >= D.bonding) {
          // The peptide bond forms and the chain moves to the A-site tRNA.
          const aa = slotA ? slotA.mesh.userData.aa : null;
          if (aa) {
            const exit = site('exit');
            s.peptide.push(aa, exit, [0.25, 1, -0.35]);
            if (slotA) slotA.mesh.userData.residue.visible = false;   // handed over
          }
          this.phase = 'translocating';
          this.phaseT = 0;
        }
        break;
      }
      case 'translocating': {
        const k = Math.min(1, this.phaseT / D.translocating);
        const e = k * k * (3 - 2 * k);
        // Everything slides back one site while the ribosome moves forward one
        // codon; on screen that is the tRNAs holding still and the ribosome
        // walking past them.
        const lerpSite = (from, to) => {
          const a = site(from), b = site(to);
          return [a[0] + (b[0] - a[0]) * e, a[1] + (b[1] - a[1]) * e, a[2] + (b[2] - a[2]) * e];
        };
        if (this.occupancy.A) this._placeTrna(this.occupancy.A, lerpSite('A', 'P'));
        if (this.occupancy.P) this._placeTrna(this.occupancy.P, lerpSite('P', 'E'));
        if (this.occupancy.E) {
          const from = site('E');
          this._placeTrna(this.occupancy.E, [from[0] - 26 * e, from[1] - 10 * e, from[2] - 20 * e], e * 1.2, e * 0.8);
        }

        if (k >= 1) {
          if (this.occupancy.E) { this.occupancy.E.mesh.visible = false; this.occupancy.E.busy = false; }
          this.occupancy.E = this.occupancy.P;
          this.occupancy.P = this.occupancy.A;
          this.occupancy.A = null;
          this.codon++;

          const nextStart = orf.start + this.codon * CODON_LEN;
          const next = this.rna.slice(nextStart, nextStart + 3);
          if (!next || next.length < 3 || CODON_TABLE[next] === '*' || CODON_TABLE[next] === undefined) {
            this.phase = 'release';
            this.phaseT = 0;
            this._log(next && CODON_TABLE[next] === '*'
              ? `Stop codon ${next} — no tRNA reads it. Release factor moves in.`
              : 'Ran off the end of the message.');
          } else {
            this.phase = 'testing';
            this.phaseT = 0;
            this._beginTesting();
          }
        }
        break;
      }
      case 'release': {
        if (this.occupancy.P) this._placeTrna(this.occupancy.P, site('P'));
        if (this.occupancy.E) this._placeTrna(this.occupancy.E, site('E'));
        if (this.phaseT >= D.release && !this.released) {
          this.released = true;
          this.done = true;
          for (const t of this.trnas) { t.mesh.visible = false; t.busy = false; }
          this.occupancy = { A: null, P: null, E: null };
          this._log('Chain released. Now it folds — and the folding is what makes it a protein.');
        }
        break;
      }
      default: break;
    }

    /* --- the peptide ----------------------------------------------------- */
    if (s.peptide.count) {
      if (!this.released) {
        // Still attached: it dangles from the exit tunnel and only loosely
        // starts to collapse, because most of it does not exist yet.
        s.peptide.relax(dt * speed, 0.35, 0.45);
        const exit = site('exit');
        const last = s.peptide.count - 1;
        s.peptide.pos[last].lerp(_v.set(exit[0], exit[1], exit[2]), Math.min(1, dt * 12));
      } else {
        // Released, and annealing: hot at first so the chain can leave the rod
        // it was extruded as, then cooling so the hydrophobic term can hold on
        // to what it finds.
        this.foldT += dt * speed;
        const heat = Math.max(0.05, 1 - this.foldT / 7);
        s.peptide.relax(dt * speed * 1.6, 1, heat);
        s.peptide.folding = Math.min(1, this.foldT / 12);
      }
    }

    /* --- the mRNA -------------------------------------------------------- */
    for (let i = 0; i < this.mrna.length; i++) {
      const p = this._nt(i);
      this.mrna.placeFree(i, p, [0, 1, 0], [0, 0, 1]);
      const inCodon = i >= codonStart && i < codonStart + 3;
      this.mrna.highlight[i] = inCodon ? 1 : (i >= orf.start && (i - orf.start) % 3 === 0 ? 0.25 : 0);
    }

    m.ribosome.userData.large.rotation.y = Math.sin(this.phaseT * 3) * 0.03;
  }

  _beginTesting() {
    const orf = this.reading;
    const codonStart = orf.start + this.codon * CODON_LEN;
    const codon = this.rna.slice(codonStart, codonStart + 3);
    const correct = CODON_TABLE[codon];

    this.rejects = [];
    const tries = Math.random() < 0.65 ? 1 : 2;
    for (let k = 0; k < tries; k++) {
      const slot = this._freeTrna();
      if (!slot) break;
      slot.busy = true;
      // A wrong tRNA: some other anticodon, carrying some other amino acid.
      let wrongCodon = codon;
      let guard = 0;
      while (wrongCodon === codon && guard++ < 20) {
        const keys = Object.keys(CODON_TABLE);
        wrongCodon = keys[(Math.random() * keys.length) | 0];
      }
      const aa = CODON_TABLE[wrongCodon];
      this._dressTrna(slot, aa === '*' ? 'G' : aa, this._anticodonFor(wrongCodon));
      const from = [
        (Math.random() - 0.5) * 120,
        40 + Math.random() * 60,
        this._nt(codonStart)[2] + (Math.random() - 0.5) * 90,
      ];
      this.rejects.push({ slot, from });
    }
    if (correct && correct !== '*') this._log(`Codon ${codon} — testing arrivals`);
  }

  _beginDocking(codonStart) {
    const codon = this.rna.slice(codonStart, codonStart + 3);
    const aa = CODON_TABLE[codon];
    const slot = this._freeTrna();
    if (!slot) return;
    slot.busy = true;
    this._dressTrna(slot, aa === '*' ? 'G' : aa, this._anticodonFor(codon));
    slot.from = [
      (Math.random() - 0.5) * 110,
      55 + Math.random() * 45,
      this._nt(codonStart)[2] + (Math.random() - 0.5) * 80,
    ];
    this.occupancy.A = slot;
    const info = AMINO[aa];
    if (info) this._log(`${codon} → ${info.tla}. Peptide bond ${this.codon + 1}.`);
  }

  /** The anticodon is the reverse complement of the codon. */
  _anticodonFor(codon) {
    const comp = { A: 'U', U: 'A', G: 'C', C: 'G' };
    return codon.split('').reverse().map((c) => comp[c] || 'A').join('');
  }

  draw(r) {
    const s = this.stage;
    r.addChain(this.mrna, { backboneColour: 0xc084fc, radius: 1.25 });

    // Anticodon-codon pairing, drawn only for a tRNA actually seated.
    for (const key of ['A', 'P']) {
      const slot = this.occupancy[key];
      if (!slot || !slot.mesh.visible) continue;
      r.addChain(slot.anticodon, { baseScale: 0.85 });
      if (key !== 'P' && this.phase !== 'bonding') continue;
      const codonStart = this.reading.start + (key === 'P' ? this.codon - 1 : this.codon) * CODON_LEN;
      for (let k = 0; k < 3; k++) {
        const mi = codonStart + k;
        if (mi < 0 || mi >= this.mrna.length) continue;
        const info = BASES[this.mrna.letters[mi]];
        if (!info) continue;
        hydrogenBonds(this.mrna, mi, slot.anticodon, 2 - k, info.hbonds, _bonds);
        for (const b of _bonds) r.addBond(b.from, b.to, { radius: 0.2, colour: 0xffe0a8 });
      }
    }
    for (const rj of this.rejects) {
      if (rj.slot.mesh.visible) r.addChain(rj.slot.anticodon, { baseScale: 0.8 });
    }
  }

  state() {
    const s = this.stage;
    const orf = this.reading;
    if (!orf || orf.start < 0) {
      return {
        caption: 'No AUG anywhere in this sequence, so there is no reading frame and nothing gets built.',
        progress: 0,
        stats: [['start codon', 'not found']],
        track: {},
      };
    }
    const peptide = s.peptide;
    const seqSoFar = [];
    for (let i = 0; i < peptide.count; i++) seqSoFar.push(peptide.aa[i]);

    let caption;
    if (this.released) {
      const rg = peptide.radiusOfGyration();
      caption = peptide.folding < 0.99
        ? 'Collapsing. Heat shakes the chain through shapes and the greasy residues hold on to each other.'
        : `Collapsed to a globule — radius of gyration ${rg.toFixed(1)} Å, down from about ${(peptide.count * 3.8 / Math.sqrt(12)).toFixed(0)} Å extended. A chain this short has no real core; that needs about fifty residues.`;
    } else if (this.phase === 'testing') {
      caption = 'Wrong tRNAs keep arriving. They cannot pair, so they fall straight back out.';
    } else if (this.phase === 'bonding') {
      caption = 'The chain is handed from the tRNA in P to the one in A. That transfer is the peptide bond.';
    } else if (this.phase === 'translocating') {
      caption = 'The ribosome steps exactly three nucleotides. Three, every time — that is the reading frame.';
    } else {
      caption = 'A charged tRNA reaches the A site and its anticodon is checked against the codon.';
    }

    return {
      caption,
      progress: orf.peptide.length ? Math.min(1, peptide.count / orf.peptide.length) : 0,
      stats: [
        ['reading frame', `starts at nt ${orf.start + 1}`],
        ['codon', `${this.codon + 1} / ${orf.peptide.length}`],
        ['residues', `${peptide.count}`],
        ['stop codon', orf.stopCodon || '—'],
        ['radius of gyration', peptide.count > 2 ? `${peptide.radiusOfGyration().toFixed(1)} Å` : '—'],
      ],
      peptide: seqSoFar,
      events: this.events,
      track: {
        mark: [0, 1, 2].map((k) => orf.start + this.codon * CODON_LEN + k),
        doneRange: [orf.start, orf.start + this.codon * CODON_LEN - 1],
        frame: { start: orf.start, codon: this.codon },
      },
    };
  }
}
