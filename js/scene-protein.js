/**
 * scene-protein.js — the protein: how it folds, how it works, and what
 * happens when you break it.
 *
 * The argument the scene is making, in three stages:
 *
 *   1. A chain of amino acids is not a protein. It comes out of the ribosome
 *      as a floppy string with no shape and no function.
 *   2. It collapses. Water squeezes the greasy residues inward, the backbone
 *      closes into helices, and a pocket appears on the surface that was in no
 *      way present in the sequence — it is made by residues that sit far apart
 *      along the chain and are brought together by the fold.
 *   3. That pocket is the active site, and it is the whole job. Substrate
 *      diffuses in, is held, is cut, and the pieces leave. Turn the heat up and
 *      the fold comes apart, the pocket stops existing, and the protein stops
 *      working — not because the chemistry changed, but because the shape did.
 *
 * Stage 3 is why this is the answer to "how do proteins work". Shape is the
 * mechanism, and the only way to show that is to take the shape away.
 */

import * as THREE from './lib/three.module.js';
import { translate, transcribe } from './bio.js';
import { Ribbon } from './geometry.js';
import { makeLabel } from './machines.js';

// Scratch vectors. `_w` in particular must not be handed out to anything that
// keeps it across a loop iteration — the protein centre lives in a field of its
// own for exactly that reason.
const _v = new THREE.Vector3();
const _w = new THREE.Vector3();
const _u = new THREE.Vector3();
const _m = new THREE.Matrix4();
const _col = new THREE.Color();

function rand(a, b) { return a + Math.random() * (b - a); }

/* -------------------------------------------------------------- substrate */

// A between substrate units. Sized against the enzyme rather than for its own
// sake: this protein is only 36 residues and about 25 A across, so a substrate
// that reads correctly as "small thing being held by a big thing" has to be
// under half that. A hexasaccharide in lysozyme's cleft is the proportion to
// aim for, and lysozyme has 129 residues to play with.
const SUB_LINK = 3.5;

/**
 * A substrate: a short chain of units that the enzyme cuts in half.
 *
 * Drawn as beads on a thread rather than as any particular molecule, because
 * the point being made is about the enzyme, not the chemistry. Lysozyme cutting
 * a cell-wall sugar chain is the picture to have in mind.
 */
class Substrate {
  constructor(parent, units = 6) {
    this.units = units;
    this.pos = [];
    this.vel = [];
    for (let i = 0; i < units; i++) {
      this.pos.push(new THREE.Vector3());
      this.vel.push(new THREE.Vector3());
    }
    this.state = 'free';       // free | binding | bound | cutting | leaving
    this.t = 0;
    this.cutAt = Math.floor(units / 2);
    this.severed = false;

    const geo = new THREE.IcosahedronGeometry(1.45, 2);
    this.mesh = new THREE.InstancedMesh(geo, new THREE.MeshStandardMaterial({
      roughness: 0.3, metalness: 0.1, envMapIntensity: 1.3,
    }), units);
    this.mesh.instanceColor = new THREE.InstancedBufferAttribute(new Float32Array(units * 3), 3);
    this.mesh.frustumCulled = false;
    parent.add(this.mesh);

    this.ribbonA = new Ribbon(units, 0.55, 6);
    this.ribbonB = new Ribbon(units, 0.55, 6);
    const mat = () => new THREE.MeshStandardMaterial({
      color: 0x63d3ff, roughness: 0.35, metalness: 0.2, envMapIntensity: 1.2,
    });
    this.meshA = new THREE.Mesh(this.ribbonA.geometry, mat());
    this.meshB = new THREE.Mesh(this.ribbonB.geometry, mat());
    this.meshA.frustumCulled = false;
    this.meshB.frustumCulled = false;
    parent.add(this.meshA, this.meshB);

    this._pts = [];
    for (let i = 0; i < units; i++) this._pts.push(new THREE.Vector3());
    this.respawn([0, 0, 0], 70);
  }

  respawn(around, radius) {
    const dir = new THREE.Vector3(rand(-1, 1), rand(-1, 1), rand(-1, 1)).normalize();
    const start = new THREE.Vector3(around[0], around[1], around[2]).addScaledVector(dir, radius);
    const along = new THREE.Vector3(rand(-1, 1), rand(-1, 1), rand(-1, 1)).normalize();
    for (let i = 0; i < this.units; i++) {
      this.pos[i].copy(start).addScaledVector(along, (i - (this.units - 1) / 2) * SUB_LINK);
      this.vel[i].set(0, 0, 0);
    }
    this.state = 'free';
    this.t = 0;
    this.severed = false;
  }

  /**
   * Hold the chain together and let it tumble.
   *
   * `bounds` keeps the substrate in the neighbourhood of the enzyme. A real
   * cell does this with concentration — there are a great many substrate
   * molecules and they are everywhere — but three of them diffusing freely in
   * an unbounded box would essentially never find a pocket 10 A across, and
   * the scene would show nothing happening for ever. Confining them to a
   * volume the enzyme sits in the middle of buys the same encounter rate
   * without pretending diffusion is directed.
   */
  _relax(step, drift, bounds = null) {
    for (let i = 0; i < this.units; i++) {
      // A severed substrate is two pieces: the bond across the cut is gone.
      for (const j of [i - 1, i + 1]) {
        if (j < 0 || j >= this.units) continue;
        if (this.severed && Math.min(i, j) === this.cutAt - 1) continue;
        _v.copy(this.pos[j]).sub(this.pos[i]);
        const d = _v.length() || 1e-4;
        this.vel[i].addScaledVector(_v.multiplyScalar(1 / d), (d - SUB_LINK) * 40 * step);
      }
      if (drift) {
        this.vel[i].x += (Math.random() - 0.5) * drift * step;
        this.vel[i].y += (Math.random() - 0.5) * drift * step;
        this.vel[i].z += (Math.random() - 0.5) * drift * step;
      }
      this.vel[i].multiplyScalar(Math.exp(-3.4 * step));
      if (this.vel[i].lengthSq() > 3600) this.vel[i].setLength(60);
      this.pos[i].addScaledVector(this.vel[i], step);

      if (bounds) {
        _v.copy(this.pos[i]).sub(bounds.centre);
        const d = _v.length();
        if (d > bounds.radius) {
          _v.multiplyScalar(1 / (d || 1));
          this.pos[i].copy(bounds.centre).addScaledVector(_v, bounds.radius);
          this.vel[i].addScaledVector(_v, -2 * this.vel[i].dot(_v) * 0.8);
        }
      }
    }
  }

  centre(target) {
    target.set(0, 0, 0);
    for (const p of this.pos) target.add(p);
    return target.multiplyScalar(1 / this.units);
  }

  render() {
    for (let i = 0; i < this.units; i++) {
      const cut = this.severed && i >= this.cutAt;
      _m.makeScale(1, 1, 1);
      _m.setPosition(this.pos[i]);
      this.mesh.setMatrixAt(i, _m);
      // setHex, not hand-unpacked bytes: colour management converts sRGB to
      // the linear working space, and writing the raw bytes skips that and
      // comes out washed out next to everything else in the scene.
      _col.setHex(cut ? 0x8affc4 : 0x63d3ff);
      this.mesh.instanceColor.setXYZ(i, _col.r, _col.g, _col.b);
      this._pts[i].copy(this.pos[i]);
    }
    this.mesh.instanceMatrix.needsUpdate = true;
    this.mesh.instanceColor.needsUpdate = true;

    // Two ribbons, so a cut chain visibly comes apart into two.
    const split = this.severed ? this.cutAt : this.units;
    this.ribbonA.update(this._pts, Math.max(0, split));
    this.meshA.visible = split >= 2;
    if (this.severed && this.units - split >= 2) {
      const tail = this._pts.slice(split);
      this.ribbonB.update(tail, tail.length);
      this.meshB.material.color.setHex(0x8affc4);
      this.meshB.visible = true;
    } else {
      this.meshB.visible = false;
    }
  }

  setVisible(v) {
    this.mesh.visible = v;
    if (!v) { this.meshA.visible = false; this.meshB.visible = false; }
  }
}

/* =========================================================== the scene ==== */

const STAGE = { CHAIN: 'chain', FOLD: 'fold', WORK: 'work' };

export class ProteinScene {
  constructor(stage) {
    this.stage = stage;
    this.id = 'protein';
    this.title = 'Proteins';

    /** Shared with the structure scene's slider: 0 is cool, 1 is boiling. */
    this.temperature = 0;

    this.substrates = [];
    for (let i = 0; i < 3; i++) this.substrates.push(new Substrate(stage.root, 5));

    // A marker for the active site, and its label.
    this.siteMarker = new THREE.Mesh(
      new THREE.IcosahedronGeometry(3.4, 2),
      new THREE.MeshBasicMaterial({
        color: 0xffd27a, transparent: true, opacity: 0.22,
        depthWrite: false, blending: THREE.AdditiveBlending,
      }),
    );
    this.siteMarker.visible = false;
    stage.root.add(this.siteMarker);

    this.siteLabel = makeLabel('active site');
    this.siteLabel.visible = false;
    stage.root.add(this.siteLabel);

    this._centre = new THREE.Vector3();
    this.bounds = null;

    this.reset();
  }

  reset() {
    this.stageName = STAGE.CHAIN;
    this.t = 0;
    this.foldT = 0;
    this.turnovers = 0;
    this.pocket = null;
    this.pocketAge = 0;
    this.events = [];
    this.peptideSeq = '';
    this.denatured = false;
    this._wasFolded = false;
  }

  _log(text) {
    this.events.push(text);
    if (this.events.length > 5) this.events.shift();
  }

  enter() {
    const s = this.stage;
    s.hideAllMachines();
    s.pool.reset();
    s.pool.setVisible(false);
    s.coding.present.fill(0);
    s.template.present.fill(0);
    s.newTop.present.fill(0);
    s.newBottom.present.fill(0);
    s.rna.present.fill(0);

    this.reset();

    // The protein is whatever the current sequence codes for — the through
    // line from the DNA views is the whole point.
    const orf = translate(transcribe(s.seq));
    this.peptideSeq = orf.start >= 0 ? orf.peptide.map((p) => p.aa).join('') : '';

    s.peptide.clear();
    s.peptide.secondary = true;
    // Ball-and-stick rather than space-filling here: the helices and the
    // pocket are the point, and full-size residues bury both.
    s.peptide.residueScale = 0.52;
    if (this.peptideSeq.length >= 3) {
      s.peptide.fromSequence(this.peptideSeq, [0, 0, -this.peptideSeq.length * 1.9]);
      this._log(`${this.peptideSeq.length} residues, extended. No shape, no function yet.`);
    } else {
      this._log('This sequence codes for nothing — no AUG, or nothing after it.');
    }

    for (const sub of this.substrates) { sub.setVisible(false); sub.respawn([0, 0, 0], rand(22, 30)); }
    this.siteMarker.visible = false;
    this.siteLabel.visible = false;
  }

  exit() {
    this.stage.peptide.secondary = false;
    this.stage.peptide.residueScale = 1;
    for (const sub of this.substrates) sub.setVisible(false);
    this.siteMarker.visible = false;
    this.siteLabel.visible = false;
  }

  /**
   * How much fold the protein is currently holding, 0..1.
   *
   * Denaturation is cooperative and it is sharp: proteins do not soften
   * gradually as you warm them, they hold their shape and then let go over a
   * few degrees, because the fold is held by many weak interactions that have
   * to give way together. Hence the narrow window rather than a slope.
   */
  _foldStrength() {
    const t = this.temperature;
    const melt = 0.52;
    const width = 0.10;
    return 1 - Math.min(1, Math.max(0, (t - (melt - width)) / (2 * width)));
  }

  update(dt, speed) {
    const s = this.stage;
    const p = s.peptide;
    this.t += dt;
    if (!p.count) return;

    const hold = this._foldStrength();
    const denatured = hold < 0.35;

    /* --- folding, and unfolding ----------------------------------------- */
    if (this.stageName === STAGE.CHAIN) {
      this.stageName = STAGE.FOLD;
      this.foldT = 0;
    }

    if (this.stageName === STAGE.FOLD || this.stageName === STAGE.WORK) {
      this.foldT += dt * speed;
      // Anneal on the way in; above the melting point, heat wins and stays won.
      const anneal = Math.max(0.06, 1 - this.foldT / 8);
      const heat = denatured ? 0.35 + this.temperature * 0.9 : anneal;
      p.relax(dt * speed * 1.6, hold, heat);
      p.folding = Math.min(1, this.foldT / 12);
    }

    /* --- the pocket ------------------------------------------------------ */
    // Re-scanning every frame is wasteful and jitters; a few times a second is
    // plenty, and it lets the site disappear promptly when the fold goes.
    this.pocketAge += dt;
    if (this.pocketAge > 0.25) {
      this.pocketAge = 0;
      this.pocket = hold > 0.5 ? p.findPocket() : null;
    }

    const folded = !!this.pocket && p.folding > 0.35 && !denatured;

    if (folded && !this._wasFolded) {
      this._wasFolded = true;
      this.stageName = STAGE.WORK;
      const lining = this.pocket.lining || [];
      const spread = lining.length
        ? Math.max(...lining) - Math.min(...lining)
        : 0;
      this._log(`Folded. A pocket has appeared, lined by residues up to ${spread} apart in the chain.`);
    }
    if (!folded && this._wasFolded && denatured) {
      this._wasFolded = false;
      this._log('Denatured. The pocket is gone, so the enzyme is gone — same atoms, no shape.');
      for (const sub of this.substrates) {
        if (sub.state !== 'free') { sub.state = 'leaving'; sub.t = 0; }
      }
    }

    if (this.pocket) {
      this.siteMarker.visible = true;
      this.siteMarker.position.copy(this.pocket.point);
      this.siteMarker.scale.setScalar(1 + Math.sin(this.t * 2.6) * 0.07);
      this.siteLabel.visible = s.showLabels;
      this.siteLabel.position.copy(this.pocket.point).add(_v.set(0, 9, 0));
    } else {
      this.siteMarker.visible = false;
      this.siteLabel.visible = false;
    }

    /* --- the catalytic cycle --------------------------------------------- */
    const step = Math.min(dt * speed, 1 / 20);
    // Held in a field, not a scratch vector: the substrate loop below reuses
    // the scratches freely and would otherwise trample it partway through.
    const centre = this._centre.set(0, 0, 0);
    for (let i = 0; i < p.count; i++) centre.add(p.pos[i]);
    centre.multiplyScalar(1 / p.count);

    // Tight enough that substrate stays in shot — if it wanders off-screen to
    // diffuse, the viewer sees an enzyme doing nothing.
    if (!this.bounds) this.bounds = { centre: new THREE.Vector3(), radius: 30 };
    this.bounds.centre.copy(centre);

    for (const sub of this.substrates) {
      sub.setVisible(true);
      sub.t += step;

      // Every state after 'free' is holding on to a pocket that may have just
      // stopped existing — denature the enzyme mid-reaction and the complex has
      // nothing left to sit in. Let go before anything reads through it.
      if (!this.pocket && sub.state !== 'free' && sub.state !== 'leaving') {
        sub.state = 'leaving';
        sub.t = 0;
      }

      switch (sub.state) {
        case 'free': {
          sub._relax(step, 190, this.bounds);

          if (folded) {
            const c = sub.centre(_v);
            const dist = c.distanceTo(this.pocket.point);

            // Electrostatic steering. Free diffusion alone would take a very
            // long time to put a 30 A substrate into a 10 A pocket — most of
            // the volume available to it is out near the edge. Real enzymes
            // have the same problem and solve it the same way: charged
            // residues around the mouth of the site set up a field that pulls
            // substrate in, which is how the fastest enzymes reach the
            // diffusion limit. Weak enough here that arrival still looks like
            // wandering rather than being reeled in.
            if (dist < 44) {
              for (let i = 0; i < sub.units; i++) {
                sub.vel[i].addScaledVector(
                  _u.copy(this.pocket.point).sub(sub.pos[i]).normalize(), 30 * step,
                );
              }
            }
            // Most passes still miss.
            if (dist < 20 && Math.random() < 0.05) {
              sub.state = 'binding';
              sub.t = 0;
              sub.from = sub.pos.map((q) => q.clone());
            }
          }
          break;
        }
        case 'binding': {
          const k = Math.min(1, sub.t / 0.9);
          const e = k * k * (3 - 2 * k);
          // Slide into the cleft, laid along it.
          const axis = _v.copy(this.pocket.point).sub(centre).normalize();
          const along = _u.set(-axis.y, axis.x, axis.z * 0.3).normalize();
          for (let i = 0; i < sub.units; i++) {
            const target = this.pocket.point.clone()
              .addScaledVector(along, (i - (sub.units - 1) / 2) * SUB_LINK);
            sub.pos[i].lerpVectors(sub.from[i], target, e);
            sub.vel[i].set(0, 0, 0);
          }
          if (k >= 1) { sub.state = 'bound'; sub.t = 0; }
          break;
        }
        case 'bound': {
          if (!folded) { sub.state = 'leaving'; sub.t = 0; break; }
          // Induced fit: the enzyme closes a little around what it is holding.
          for (let i = 0; i < sub.units; i++) {
            sub.pos[i].addScaledVector(
              _v.copy(this.pocket.point).sub(sub.pos[i]),
              Math.min(1, step * 2.2) * 0.25,
            );
          }
          if (sub.t > 0.75) { sub.state = 'cutting'; sub.t = 0; }
          break;
        }
        case 'cutting': {
          if (!sub.severed) {
            sub.severed = true;
            this.turnovers++;
            this.stage.pool.sparks.burst(
              [this.pocket.point.x, this.pocket.point.y, this.pocket.point.z], 10,
            );
            this._log(`Cut. Turnover ${this.turnovers} — the enzyme itself is unchanged.`);
          }
          // The two halves separate along the chain's own axis, so the break
          // reads as a break rather than as an explosion.
          _v.copy(sub.pos[sub.units - 1]).sub(sub.pos[0]);
          if (_v.lengthSq() < 1e-6) _v.set(1, 0, 0);
          _v.normalize();
          for (let i = 0; i < sub.units; i++) {
            const side = i < sub.cutAt ? -1 : 1;
            sub.vel[i].addScaledVector(_v, 26 * side);
          }
          sub._relax(step, 40);
          if (sub.t > 0.5) { sub.state = 'leaving'; sub.t = 0; }
          break;
        }
        case 'leaving': {
          sub._relax(step, 120, this.bounds);
          for (let i = 0; i < sub.units; i++) {
            sub.vel[i].addScaledVector(_v.copy(sub.pos[i]).sub(centre).normalize(), 40 * step);
          }
          if (sub.t > 2.0) sub.respawn([centre.x, centre.y, centre.z], rand(24, 30));
          break;
        }
        default: break;
      }
      sub.render();
    }

    this.denatured = denatured;
  }

  draw() { /* the peptide and substrates draw themselves */ }

  state() {
    const p = this.stage.peptide;
    if (!p.count) {
      return {
        caption: 'Nothing to fold — this sequence has no reading frame. Try the villin headpiece preset.',
        progress: 0,
        stats: [['residues', '0']],
        track: {},
      };
    }

    const hold = this._foldStrength();
    const helix = p.helixContent();
    const rg = p.radiusOfGyration();
    const extended = p.count * 3.8 / Math.sqrt(12);

    let caption;
    if (this.denatured) {
      caption = 'Denatured. Every atom and every bond is still there — only the shape has gone, and with it the pocket and the whole function.';
    } else if (!this._wasFolded) {
      caption = 'A chain with no shape. Heat shakes it through conformations while the greasy residues hold on to whatever brings them together.';
    } else if (this.turnovers > 0) {
      caption = 'Working. Substrate diffuses in, is held in the cleft, is cut, and the pieces leave. The enzyme comes out of it unchanged, which is what makes it a catalyst.';
    } else {
      caption = 'Folded. The pocket on the surface is made by residues far apart in the sequence — it exists only because the chain folded.';
    }

    return {
      caption,
      progress: Math.min(1, p.folding),
      stats: [
        ['residues', `${p.count}`],
        ['radius of gyration', `${rg.toFixed(1)} Å (extended ${extended.toFixed(0)})`],
        ['helix content', `${(helix * 100).toFixed(0)}%`],
        ['fold held', `${(hold * 100).toFixed(0)}%`],
        ['active site', this.pocket ? 'open' : '—'],
        ['turnovers', `${this.turnovers}`],
      ],
      events: this.events,
      peptide: p.aa.slice(0, p.count),
      track: {},
    };
  }
}
