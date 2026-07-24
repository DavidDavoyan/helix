/**
 * peptide.js — the protein as it comes off the ribosome, and what happens next.
 *
 * Residues are drawn at their real spacing: 3.8 A between consecutive alpha
 * carbons, which is fixed by the geometry of a trans peptide bond and does not
 * vary with the amino acid. Sphere sizes come from measured side-chain volumes.
 *
 * The folding is a toy — a hydrophobic collapse, run as a few hundred steps of
 * relaxation with three forces:
 *
 *   1. the chain holds itself at 3.8 A between neighbours,
 *   2. residues cannot overlap,
 *   3. hydrophobic residues attract each other and charged ones do not.
 *
 * That is not how you predict a structure. But it is the actual reason
 * proteins fold at all — water squeezes the greasy residues together — and the
 * globule it produces has a hydrophobic core and a charged surface, which is
 * the thing worth seeing.
 */

import * as THREE from './lib/three.module.js';
import { AMINO, RESIDUE_VOLUME, HYDROPATHY } from './bio.js';
import { Ribbon } from './geometry.js';

const CA_SPACING = 3.8;   // A, alpha-carbon to alpha-carbon across a peptide bond

const _m = new THREE.Matrix4();
const _c = new THREE.Color();
const _d = new THREE.Vector3();

export class Peptide {
  constructor(parent, { capacity = 220 } = {}) {
    this.capacity = capacity;
    this.count = 0;
    this.folding = 0;         // 0 = extended nascent chain, 1 = fully relaxed

    this.aa = new Array(capacity).fill('G');
    this.pos = [];
    this.vel = [];
    for (let i = 0; i < capacity; i++) {
      this.pos.push(new THREE.Vector3());
      this.vel.push(new THREE.Vector3());
    }

    const geo = new THREE.IcosahedronGeometry(1, 2);
    this.mesh = new THREE.InstancedMesh(geo, new THREE.MeshStandardMaterial({
      roughness: 0.36, metalness: 0.08, envMapIntensity: 1.0,
    }), capacity);
    this.mesh.instanceColor = new THREE.InstancedBufferAttribute(new Float32Array(capacity * 3), 3);
    this.mesh.frustumCulled = false;
    this.mesh.count = 0;
    parent.add(this.mesh);

    this.ribbon = new Ribbon(capacity, 1.05, 6);
    this.ribbonMesh = new THREE.Mesh(this.ribbon.geometry, new THREE.MeshStandardMaterial({
      color: 0xd8e3f2, roughness: 0.42, metalness: 0.12, envMapIntensity: 1.0,
    }));
    this.ribbonMesh.frustumCulled = false;
    parent.add(this.ribbonMesh);

    this._pts = [];
    for (let i = 0; i < capacity; i++) this._pts.push(new THREE.Vector3());
  }

  clear() {
    this.count = 0;
    this.folding = 0;
    this.mesh.count = 0;
    this.ribbon.update([new THREE.Vector3(), new THREE.Vector3()], 0);
  }

  /**
   * Append a residue. It appears at the ribosome's exit tunnel and pushes the
   * chain out ahead of it — a nascent chain is extended, not folded, because
   * it has not all been made yet.
   */
  push(aa, at, dir) {
    if (this.count >= this.capacity) return;
    const i = this.count++;
    this.aa[i] = aa;
    if (i === 0) {
      this.pos[i].set(at[0], at[1], at[2]);
    } else {
      _d.set(dir[0], dir[1], dir[2]).normalize();
      // A little wander, so a fresh chain looks like a coil rather than a rod.
      _d.x += (Math.random() - 0.5) * 0.75;
      _d.y += (Math.random() - 0.5) * 0.75;
      _d.z += (Math.random() - 0.5) * 0.75;
      _d.normalize().multiplyScalar(CA_SPACING);
      this.pos[i].copy(this.pos[i - 1]).add(_d);
    }
    this.vel[i].set(0, 0, 0);
  }

  /** Radius drawn for residue i, from its side-chain volume. */
  radiusOf(i) {
    const v = RESIDUE_VOLUME[this.aa[i]] || 120;
    return 0.62 * Math.cbrt(v * 3 / (4 * Math.PI));   // ~1.9-2.5 A
  }

  /**
   * One relaxation step. `strength` scales the folding forces, so a chain can
   * be left loose while it is still being made and collapsed once released.
   *
   * The thermal term is not decoration. A chain extruded straight out of a
   * ribosome sits in an unstable equilibrium — every pairwise force along it is
   * collinear, so it will stay a rod forever no matter how strong the
   * hydrophobic attraction is. Something has to break the symmetry, and in a
   * cell that something is heat. Folding is thermal motion exploring shapes and
   * the hydrophobic effect keeping the ones that bury the greasy residues;
   * without the first half, the second half has nothing to select from.
   */
  relax(dt, strength = 1, heat = 1) {
    const n = this.count;
    if (n < 3) return;
    const step = Math.min(dt, 1 / 30);
    // Heat has to start above the hydrophobic force and end below it. Too cold
    // and the chain never leaves the rod it was extruded as; too hot and it is
    // simply denatured, with the greasy residues no more buried than anything
    // else. The caller anneals this down as the fold settles.
    const kT = 42 * heat;

    // Centroid, for the weak centring that stands in for solvent pressure.
    const cx = { x: 0, y: 0, z: 0 };
    for (let i = 0; i < n; i++) { cx.x += this.pos[i].x; cx.y += this.pos[i].y; cx.z += this.pos[i].z; }
    cx.x /= n; cx.y /= n; cx.z /= n;

    for (let i = 0; i < n; i++) {
      const p = this.pos[i];
      const v = this.vel[i];
      const hi = HYDROPATHY[this.aa[i]] || 0;

      // 1. Chain bonds — a stiff spring at 3.8 A either side.
      for (const j of [i - 1, i + 1]) {
        if (j < 0 || j >= n) continue;
        _d.copy(this.pos[j]).sub(p);
        const d = _d.length() || 1e-4;
        _d.multiplyScalar(1 / d);
        v.addScaledVector(_d, (d - CA_SPACING) * 60 * step);
      }

      // 2. Excluded volume, and 3. hydrophobic attraction. Both are pairwise,
      //    and at these chain lengths the O(n^2) sweep is free.
      for (let j = 0; j < n; j++) {
        if (j === i || Math.abs(j - i) === 1) continue;
        _d.copy(this.pos[j]).sub(p);
        const d2 = _d.lengthSq();
        if (d2 < 1e-6) continue;
        const d = Math.sqrt(d2);
        _d.multiplyScalar(1 / d);

        const rsum = this.radiusOf(i) + this.radiusOf(j) + 0.9;
        if (d < rsum) {
          v.addScaledVector(_d, -(rsum - d) * 42 * step);   // push apart
        } else if (d < 18 && strength > 0) {
          const hj = HYDROPATHY[this.aa[j]] || 0;
          // Two greasy residues seek each other; a greasy and a charged one do
          // not. Scaled so the strongest attraction is Ile-Ile and the weakest
          // is Arg-Asp.
          const affinity = (hi * hj) / 20.25;
          if (affinity > 0) {
            v.addScaledVector(_d, affinity * 60 * strength * step * (1 - d / 18));
          }
        }
      }

      // Heat. Isotropic, so it adds no drift — it only lets the chain find
      // conformations the hydrophobic term can then hold on to.
      if (kT > 0) {
        v.x += (Math.random() - 0.5) * kT * step;
        v.y += (Math.random() - 0.5) * kT * step;
        v.z += (Math.random() - 0.5) * kT * step;
      }

      // Charged residues drift outward; the core is no place for a charge.
      if (hi < -2.5 && strength > 0) {
        _d.set(p.x - cx.x, p.y - cx.y, p.z - cx.z);
        const d = _d.length() || 1e-4;
        v.addScaledVector(_d.multiplyScalar(1 / d), 7 * strength * step);
      }
    }

    const damp = Math.exp(-6.5 * step);
    for (let i = 0; i < n; i++) {
      this.vel[i].multiplyScalar(damp);
      // Clamp, because a stiff bond spring plus a big dt is a rocket.
      if (this.vel[i].lengthSq() > 900) this.vel[i].setLength(30);
      this.pos[i].addScaledVector(this.vel[i], step);
    }
  }

  /** Move the whole chain so residue 0 sits at `at`. */
  anchor(at) {
    if (!this.count) return;
    _d.set(at[0], at[1], at[2]).sub(this.pos[0]);
    for (let i = 0; i < this.count; i++) this.pos[i].add(_d);
  }

  /** Longest dimension, in angstroms — a crude radius of gyration readout. */
  radiusOfGyration() {
    const n = this.count;
    if (n < 2) return 0;
    let cx = 0, cy = 0, cz = 0;
    for (let i = 0; i < n; i++) { cx += this.pos[i].x; cy += this.pos[i].y; cz += this.pos[i].z; }
    cx /= n; cy /= n; cz /= n;
    let s = 0;
    for (let i = 0; i < n; i++) {
      s += (this.pos[i].x - cx) ** 2 + (this.pos[i].y - cy) ** 2 + (this.pos[i].z - cz) ** 2;
    }
    return Math.sqrt(s / n);
  }

  /** Fraction of hydrophobic residues that ended up buried. */
  buriedFraction() {
    const n = this.count;
    if (n < 4) return 0;
    let cx = 0, cy = 0, cz = 0;
    for (let i = 0; i < n; i++) { cx += this.pos[i].x; cy += this.pos[i].y; cz += this.pos[i].z; }
    cx /= n; cy /= n; cz /= n;
    const rg = this.radiusOfGyration() || 1;
    let hydro = 0, buried = 0;
    for (let i = 0; i < n; i++) {
      if ((HYDROPATHY[this.aa[i]] || 0) <= 0) continue;
      hydro++;
      const d = Math.hypot(this.pos[i].x - cx, this.pos[i].y - cy, this.pos[i].z - cz);
      if (d < rg) buried++;
    }
    return hydro ? buried / hydro : 0;
  }

  render() {
    const n = this.count;
    this.mesh.count = n;
    if (!n) { this.ribbonMesh.visible = false; return; }

    for (let i = 0; i < n; i++) {
      const r = this.radiusOf(i);
      _m.makeScale(r, r, r);
      _m.setPosition(this.pos[i]);
      this.mesh.setMatrixAt(i, _m);
      const info = AMINO[this.aa[i]] || AMINO.G;
      _c.setHex(info.colour);
      this.mesh.instanceColor.setXYZ(i, _c.r, _c.g, _c.b);
      this._pts[i].copy(this.pos[i]);
    }
    this.mesh.instanceMatrix.needsUpdate = true;
    this.mesh.instanceColor.needsUpdate = true;

    if (n >= 2) {
      this.ribbon.update(this._pts, n);
      this.ribbonMesh.visible = true;
    } else {
      this.ribbonMesh.visible = false;
    }
  }

  setEnvironment(env) {
    this.mesh.material.envMap = env;
    this.mesh.material.needsUpdate = true;
    this.ribbonMesh.material.envMap = env;
    this.ribbonMesh.material.needsUpdate = true;
  }
}
