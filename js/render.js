/**
 * render.js — draws Chains.
 *
 * Everything on screen that is a nucleic acid comes through here. The renderer
 * is stateless between frames: processes rewrite chain slots, then each frame
 * we walk the chains and refill instance buffers. Base plates, sugars,
 * phosphates and hydrogen bonds are four instanced meshes; backbones are a
 * pool of Ribbons whose vertices are rewritten rather than reallocated.
 *
 * A run of present residues gets one ribbon. That is how Okazaki fragments
 * draw themselves for free: three unjoined runs are three ribbons, and when
 * ligase seals them the runs merge and the pool hands two ribbons back.
 */

import * as THREE from './lib/three.module.js';
import { BASES } from './bio.js';
import {
  makeBaseGeometries, makeSugarGeometry, makePhosphateGeometry,
  makeUnitBond, Ribbon,
} from './geometry.js';

const _m = new THREE.Matrix4();
const _q = new THREE.Quaternion();
const _p = new THREE.Vector3();
const _s = new THREE.Vector3();
const _c = new THREE.Color();
const _x = new THREE.Vector3();
const _y = new THREE.Vector3();
const _z = new THREE.Vector3();
const _a = new THREE.Vector3();
const _b = new THREE.Vector3();

/** Build a matrix whose +x is `out`, +z is `up`, positioned at `pos`. */
function frameMatrix(target, px, py, pz, ox, oy, oz, ux, uy, uz, scale) {
  _x.set(ox, oy, oz).normalize();
  _z.set(ux, uy, uz);
  // Re-orthogonalise: the up vector is advisory, the out vector is exact.
  _z.sub(_x.clone().multiplyScalar(_z.dot(_x)));
  if (_z.lengthSq() < 1e-8) _z.set(0, 0, 1).sub(_x.clone().multiplyScalar(_x.z));
  _z.normalize();
  _y.crossVectors(_z, _x);

  target.set(
    _x.x * scale, _y.x * scale, _z.x * scale, px,
    _x.y * scale, _y.y * scale, _z.y * scale, py,
    _x.z * scale, _y.z * scale, _z.z * scale, pz,
    0, 0, 0, 1,
  );
  return target;
}

export class NucleicRenderer {
  constructor(parent, { maxResidues = 900, maxBonds = 1400, maxRibbons = 40 } = {}) {
    this.parent = parent;
    this.maxResidues = maxResidues;

    // Aromatic rings are planar. Drawn thin, so a stack of them reads as the
    // stack of flat plates it is, with the 3.38 A rise showing between them.
    const baseGeo = makeBaseGeometries(0.62);
    const mat = () => new THREE.MeshStandardMaterial({
      roughness: 0.42, metalness: 0.08, envMapIntensity: 1.1,
    });

    this.purines = new THREE.InstancedMesh(baseGeo.fused, mat(), maxResidues);
    this.pyrimidines = new THREE.InstancedMesh(baseGeo.single, mat(), maxResidues);

    const sugarMat = new THREE.MeshStandardMaterial({
      color: 0x93a6c0, roughness: 0.5, metalness: 0.12, envMapIntensity: 1.1,
    });
    this.sugars = new THREE.InstancedMesh(makeSugarGeometry(0.8), sugarMat, maxResidues);

    const phosMat = new THREE.MeshStandardMaterial({
      color: 0xffa53d, roughness: 0.3, metalness: 0.4,
      emissive: 0x341800, envMapIntensity: 1.4,
    });
    this.phosphates = new THREE.InstancedMesh(makePhosphateGeometry(1.5), phosMat, maxResidues);

    // Hydrogen bonds are drawn as thin bright rods. They are not covalent and
    // they are what melts, so they get their own look: emissive, unlit-ish.
    const bondMat = new THREE.MeshStandardMaterial({
      color: 0xdfe9ff, emissive: 0x4a6ea8, roughness: 0.9, metalness: 0,
      transparent: true, opacity: 0.85,
    });
    this.bonds = new THREE.InstancedMesh(makeUnitBond(6), bondMat, maxBonds);

    for (const m of [this.purines, this.pyrimidines, this.sugars, this.phosphates, this.bonds]) {
      m.frustumCulled = false;
      m.castShadow = false;
      m.receiveShadow = false;
      parent.add(m);
    }
    this.purines.instanceColor = new THREE.InstancedBufferAttribute(new Float32Array(maxResidues * 3), 3);
    this.pyrimidines.instanceColor = new THREE.InstancedBufferAttribute(new Float32Array(maxResidues * 3), 3);
    this.bonds.instanceColor = new THREE.InstancedBufferAttribute(new Float32Array(maxBonds * 3), 3);

    // Ribbon pool.
    this.ribbons = [];
    for (let i = 0; i < maxRibbons; i++) {
      const r = new Ribbon(320, 1.0, 8);
      const m = new THREE.Mesh(r.geometry, new THREE.MeshStandardMaterial({
        color: 0x8ea3bd, roughness: 0.38, metalness: 0.2, envMapIntensity: 1.2,
        side: THREE.DoubleSide,
      }));
      m.frustumCulled = false;
      m.visible = false;
      parent.add(m);
      this.ribbons.push({ ribbon: r, mesh: m, points: [] });
    }

    this._counts = { pur: 0, pyr: 0, sug: 0, pho: 0, bond: 0, rib: 0 };
    this.showBases = true;
    this.showBackbone = true;
    this.showSugars = true;
    this.showHBonds = true;
  }

  begin() {
    const c = this._counts;
    c.pur = 0; c.pyr = 0; c.sug = 0; c.pho = 0; c.bond = 0; c.rib = 0;
  }

  /** Draw one chain: its backbone ribbon runs and its per-residue pieces. */
  addChain(chain, { backboneColour = 0x8ea3bd, radius = 1.0, baseScale = 1 } = {}) {
    const c = this._counts;

    if (this.showBackbone && chain.covalent) {
      for (const [s, e] of chain.runs()) {
        if (e - s < 1) continue;                       // a lone residue has no tube
        if (c.rib >= this.ribbons.length) break;
        const slot = this.ribbons[c.rib++];
        const pts = slot.points;
        let n = 0;
        for (let i = s; i <= e; i++) {
          if (!pts[n]) pts[n] = new THREE.Vector3();
          chain.getBackbone(i, pts[n]);
          n++;
        }
        slot.ribbon.radius = radius;
        slot.ribbon.update(pts, n);
        slot.mesh.material.color.setHex(backboneColour);
        slot.mesh.visible = true;
      }
    }

    for (let i = 0; i < chain.length; i++) {
      if (!chain.present[i]) continue;
      const letter = chain.letters[i];
      const info = BASES[letter];
      if (!info) continue;

      const o = i * 3;
      const sc = chain.scale[i] * baseScale;
      if (sc < 0.02) continue;

      const hi = chain.highlight[i];

      if (this.showBases) {
        const target = info.ring === 'fused' ? this.purines : this.pyrimidines;
        const idx = info.ring === 'fused' ? c.pur++ : c.pyr++;
        if (idx < this.maxResidues) {
          frameMatrix(_m,
            chain.pos[o], chain.pos[o + 1], chain.pos[o + 2],
            chain.out[o], chain.out[o + 1], chain.out[o + 2],
            chain.up[o], chain.up[o + 1], chain.up[o + 2], sc);
          target.setMatrixAt(idx, _m);
          _c.setHex(info.colour);
          if (hi > 0) _c.lerp(_c.clone().offsetHSL(0, 0.1, 0.35), hi).multiplyScalar(1 + hi * 1.4);
          target.instanceColor.setXYZ(idx, _c.r, _c.g, _c.b);
        }
      }

      if (this.showSugars && c.sug < this.maxResidues) {
        frameMatrix(_m,
          chain.sugar[o], chain.sugar[o + 1], chain.sugar[o + 2],
          chain.out[o], chain.out[o + 1], chain.out[o + 2],
          chain.up[o], chain.up[o + 1], chain.up[o + 2], sc * 0.95);
        this.sugars.setMatrixAt(c.sug++, _m);
      }

      if (c.pho < this.maxResidues) {
        _m.makeScale(sc, sc, sc);
        _m.setPosition(chain.backbone[o], chain.backbone[o + 1], chain.backbone[o + 2]);
        this.phosphates.setMatrixAt(c.pho++, _m);
      }
    }
  }

  /**
   * A bond between two points. Used for hydrogen bonds within a pair and for
   * the covalent link a polymerase forms when it joins a nucleotide to the
   * chain.
   */
  addBond(from, to, { radius = 0.22, colour = 0xdfe9ff } = {}) {
    const c = this._counts;
    if (c.bond >= this.bonds.count) return;
    _a.set(from[0], from[1], from[2]);
    _b.set(to[0], to[1], to[2]);
    const len = _a.distanceTo(_b);
    if (len < 1e-4) return;
    _p.copy(_b).sub(_a).normalize();
    _q.setFromUnitVectors(_y.set(0, 1, 0), _p);
    _s.set(radius, len, radius);
    _m.compose(_a, _q, _s);
    const idx = c.bond++;
    this.bonds.setMatrixAt(idx, _m);
    _c.setHex(colour);
    this.bonds.instanceColor.setXYZ(idx, _c.r, _c.g, _c.b);
  }

  end() {
    const c = this._counts;
    this.purines.count = c.pur;
    this.pyrimidines.count = c.pyr;
    this.sugars.count = this.showSugars ? c.sug : 0;
    this.phosphates.count = c.pho;
    this.bonds.count = this.showHBonds ? c.bond : 0;

    for (const m of [this.purines, this.pyrimidines, this.sugars, this.phosphates, this.bonds]) {
      m.instanceMatrix.needsUpdate = true;
      if (m.instanceColor) m.instanceColor.needsUpdate = true;
    }
    for (let i = c.rib; i < this.ribbons.length; i++) this.ribbons[i].mesh.visible = false;
  }

  setEnvironment(env) {
    for (const m of [this.purines, this.pyrimidines, this.sugars, this.phosphates, this.bonds]) {
      m.material.envMap = env;
      m.material.needsUpdate = true;
    }
    for (const r of this.ribbons) {
      r.mesh.material.envMap = env;
      r.mesh.material.needsUpdate = true;
    }
  }
}
