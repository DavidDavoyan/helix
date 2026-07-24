/**
 * pool.js — the nucleotide pool, and how a polymerase picks from it.
 *
 * A cell does not fetch the nucleotide it needs. Free triphosphates diffuse in
 * and out of the active site constantly, and almost all of them are wrong and
 * fall straight back out; the right one is held because it can make its
 * hydrogen bonds and the site closes around it. Selection is by rejection.
 *
 * So that is what this does. Every incorporation tries one or two wrong bases
 * first, each of which is drawn into the site, fails to pair, and is flung
 * back into the pool. Only then does a matching base arrive and lock. The
 * fidelity you see on screen is the fidelity of the mechanism.
 */

import * as THREE from './lib/three.module.js';
import { Chain } from './chain.js';
import { DNA_LETTERS, BASES } from './bio.js';

const FREE = 0, TESTING = 1, REJECTED = 2, DOCKING = 3, SPENT = 4;

const _v = new THREE.Vector3();
const _w = new THREE.Vector3();

function rand(a, b) { return a + Math.random() * (b - a); }
function smootherstep(t) { t = Math.min(1, Math.max(0, t)); return t * t * t * (t * (t * 6 - 15) + 10); }

export class NucleotidePool {
  /**
   * @param {object} opts
   * @param {number} opts.count    how many free nucleotides drift about
   * @param {number[]} opts.bounds half-extents of the box they diffuse in
   */
  constructor({ count = 150, bounds = [70, 70, 150], rna = false } = {}) {
    this.count = count;
    this.bounds = bounds;
    this.rna = rna;
    this.letters = rna ? ['A', 'U', 'G', 'C'] : DNA_LETTERS;

    this.chain = new Chain(count, { kind: rna ? 'rna' : 'dna', label: 'free nucleotides' });
    this.chain.covalent = false;
    this.chain.length = count;

    this.state = new Uint8Array(count);
    this.vel = new Float32Array(count * 3);
    this.spin = new Float32Array(count * 3);   // tumbling axis * rate
    this.angle = new Float32Array(count);
    this.timer = new Float32Array(count);

    for (let i = 0; i < count; i++) {
      this.chain.letters[i] = this.letters[(Math.random() * 4) | 0];
      this.chain.present[i] = 1;
      this._respawn(i, true);
    }

    /** Active incorporation attempts. */
    this.jobs = [];
    this.sparks = new Sparks(220);
  }

  _respawn(i, anywhere = false) {
    const [bx, by, bz] = this.bounds;
    const o = i * 3;
    const p = [
      rand(-bx, bx),
      rand(-by, by),
      anywhere ? rand(-bz, bz) : rand(-bz, bz),
    ];
    // Keep the pool out of the middle, where the machinery lives, so the view
    // down the helix stays clear.
    if (Math.hypot(p[0], p[1]) < 26) {
      const s = 26 / (Math.hypot(p[0], p[1]) || 1);
      p[0] *= s * rand(1, 2.2); p[1] *= s * rand(1, 2.2);
    }
    const out = new THREE.Vector3(rand(-1, 1), rand(-1, 1), rand(-1, 1)).normalize();
    const up = new THREE.Vector3(rand(-1, 1), rand(-1, 1), rand(-1, 1));
    up.sub(out.clone().multiplyScalar(up.dot(out))).normalize();
    this.chain.placeFree(i, p, [out.x, out.y, out.z], [up.x, up.y, up.z]);
    this.state[i] = FREE;
    this.chain.scale[i] = 1;
    this.chain.highlight[i] = 0;
    this.vel[o] = rand(-6, 6); this.vel[o + 1] = rand(-6, 6); this.vel[o + 2] = rand(-6, 6);
    this.spin[o] = rand(-1, 1); this.spin[o + 1] = rand(-1, 1); this.spin[o + 2] = rand(-1, 1);
  }

  /** Nearest free nucleotide carrying `letter`, or -1. */
  _findFree(letter, near) {
    let best = -1, bestD = Infinity;
    for (let i = 0; i < this.count; i++) {
      if (this.state[i] !== FREE) continue;
      if (this.chain.letters[i] !== letter) continue;
      const o = i * 3;
      const d = (this.chain.pos[o] - near[0]) ** 2
        + (this.chain.pos[o + 1] - near[1]) ** 2
        + (this.chain.pos[o + 2] - near[2]) ** 2;
      if (d < bestD) { bestD = d; best = i; }
    }
    if (best < 0) {
      // Nothing of that letter is loose — recycle a spent one rather than
      // stalling the polymerase.
      for (let i = 0; i < this.count; i++) {
        if (this.state[i] === FREE) {
          this.chain.letters[i] = letter;
          return i;
        }
      }
    }
    return best;
  }

  /**
   * Start an incorporation.
   *
   * @param {object} req
   * @param {string} req.letter   the base that will actually be accepted
   * @param {Function} req.site   () => ({pos:[x,y,z], out:[..], up:[..]}) — read
   *                              every frame, because the polymerase is moving
   * @param {number} req.duration seconds for the whole attempt
   * @param {number} req.wrongTries how many mismatches to test first
   * @param {Function} req.onLock called once, when the right base seats
   * @returns {object} a job handle with a `done` flag
   */
  incorporate({ letter, site, duration = 0.7, wrongTries = 1, onLock = null }) {
    const s = site();
    const wrong = [];
    const others = this.letters.filter((l) => l !== letter);
    for (let k = 0; k < wrongTries; k++) wrong.push(others[(Math.random() * others.length) | 0]);

    const job = {
      letter, site, duration, onLock,
      queue: [...wrong.map((l) => ({ l, ok: false })), { l: letter, ok: true }],
      idx: -1, t: 0, current: -1, done: false, locked: false,
    };
    this.jobs.push(job);
    this._nextCandidate(job, s);
    return job;
  }

  _nextCandidate(job, s) {
    job.idx++;
    if (job.idx >= job.queue.length) { job.done = true; job.current = -1; return; }
    const spec = job.queue[job.idx];
    const i = this._findFree(spec.l, s.pos);
    if (i < 0) { job.done = true; return; }
    this.state[i] = TESTING;
    this.chain.letters[i] = spec.l;
    job.current = i;
    job.t = 0;
    // A rejected test is quick; the accepted one takes the rest of the budget.
    job.phaseTime = spec.ok ? job.duration * 0.55 : job.duration * 0.22;
    const o = i * 3;
    job.from = [this.chain.pos[o], this.chain.pos[o + 1], this.chain.pos[o + 2]];
    job.fromOut = [this.chain.out[o], this.chain.out[o + 1], this.chain.out[o + 2]];
  }

  update(dt, speed = 1) {
    const step = dt * speed;
    const [bx, by, bz] = this.bounds;

    // --- diffusion -------------------------------------------------------
    for (let i = 0; i < this.count; i++) {
      const st = this.state[i];
      if (st === TESTING || st === DOCKING || st === SPENT) continue;
      const o = i * 3;

      // Brownian kick, viscous drag. Nucleotides in water are overdamped, so
      // they jitter rather than coast.
      this.vel[o] += rand(-1, 1) * 60 * step;
      this.vel[o + 1] += rand(-1, 1) * 60 * step;
      this.vel[o + 2] += rand(-1, 1) * 60 * step;
      const drag = Math.exp(-3.2 * step);
      this.vel[o] *= drag; this.vel[o + 1] *= drag; this.vel[o + 2] *= drag;

      let px = this.chain.pos[o] + this.vel[o] * step;
      let py = this.chain.pos[o + 1] + this.vel[o + 1] * step;
      let pz = this.chain.pos[o + 2] + this.vel[o + 2] * step;

      if (Math.abs(px) > bx) { px = Math.sign(px) * bx; this.vel[o] *= -0.7; }
      if (Math.abs(py) > by) { py = Math.sign(py) * by; this.vel[o + 1] *= -0.7; }
      if (Math.abs(pz) > bz) { pz = Math.sign(pz) * bz; this.vel[o + 2] *= -0.7; }

      // Tumble: rotate the base's own frame slowly about a fixed axis.
      this.angle[i] += step * 1.4;
      const a = this.angle[i];
      _v.set(this.spin[o], this.spin[o + 1], this.spin[o + 2]).normalize();
      _w.set(Math.cos(a), Math.sin(a * 0.83), Math.sin(a * 0.61)).normalize();
      _w.sub(_v.clone().multiplyScalar(_w.dot(_v))).normalize();
      this.chain.placeFree(i, [px, py, pz], [_w.x, _w.y, _w.z], [_v.x, _v.y, _v.z]);

      if (st === REJECTED) {
        this.timer[i] -= step;
        this.chain.highlight[i] = Math.max(0, this.timer[i] / 0.4);
        if (this.timer[i] <= 0) this.state[i] = FREE;
      }
    }

    // --- active incorporations ------------------------------------------
    for (let j = this.jobs.length - 1; j >= 0; j--) {
      const job = this.jobs[j];
      if (job.done) { this.jobs.splice(j, 1); continue; }

      const i = job.current;
      if (i < 0) { this.jobs.splice(j, 1); continue; }

      job.t += step;
      const k = Math.min(1, job.t / job.phaseTime);
      const e = smootherstep(k);
      const s = job.site();
      const spec = job.queue[job.idx];

      // A candidate that will be rejected only ever gets close; the accepted
      // one goes all the way in.
      const approach = spec.ok ? 1 : 0.72;
      const off = spec.ok ? 0 : 6.5 * (1 - e * 0.35);

      const tx = s.pos[0] + s.out[0] * -off;
      const ty = s.pos[1] + s.out[1] * -off;
      const tz = s.pos[2] + s.out[2] * -off;

      const px = job.from[0] + (tx - job.from[0]) * e * approach + (spec.ok ? 0 : 0);
      const py = job.from[1] + (ty - job.from[1]) * e * approach;
      const pz = job.from[2] + (tz - job.from[2]) * e * approach;

      // Orientation swings to the site's frame as it comes in.
      const ox = job.fromOut[0] + (s.out[0] - job.fromOut[0]) * e;
      const oy = job.fromOut[1] + (s.out[1] - job.fromOut[1]) * e;
      const oz = job.fromOut[2] + (s.out[2] - job.fromOut[2]) * e;
      const l = Math.hypot(ox, oy, oz) || 1;

      this.chain.placeFree(i,
        spec.ok ? [px, py, pz] : [job.from[0] + (tx - job.from[0]) * e, job.from[1] + (ty - job.from[1]) * e, job.from[2] + (tz - job.from[2]) * e],
        [ox / l, oy / l, oz / l], s.up);
      this.chain.highlight[i] = spec.ok ? e : e * 0.5;

      if (k >= 1) {
        if (spec.ok) {
          // Seated. Hand the residue to the growing chain and release the
          // pyrophosphate — that cleavage is what pays for the bond.
          if (job.onLock) job.onLock(i);
          job.locked = true;
          this.state[i] = SPENT;
          this.chain.scale[i] = 0;
          this.sparks.burst(s.pos, 6);
          setTimeout(() => { if (this.state[i] === SPENT) this._respawn(i, true); }, 0);
          job.done = true;
        } else {
          // Wrong shape, wrong donors and acceptors. Out it goes.
          this.state[i] = REJECTED;
          this.timer[i] = 0.4;
          const o = i * 3;
          const dir = new THREE.Vector3(rand(-1, 1), rand(-1, 1), rand(-1, 1)).normalize();
          this.vel[o] = dir.x * 90; this.vel[o + 1] = dir.y * 90; this.vel[o + 2] = dir.z * 90;
          this._nextCandidate(job, s);
        }
      }
    }

    this.sparks.update(dt);
  }

  /** Abandon everything in flight — used when the mode or sequence changes. */
  reset() {
    for (const job of this.jobs) if (job.current >= 0) this._respawn(job.current, true);
    this.jobs.length = 0;
    for (let i = 0; i < this.count; i++) if (this.state[i] !== FREE) this._respawn(i, true);
    this.sparks.clear();
  }

  setVisible(v) {
    for (let i = 0; i < this.count; i++) {
      if (this.state[i] === SPENT) continue;
      this.chain.present[i] = v ? 1 : 0;
    }
    this.sparks.points.visible = v;
  }
}

/* -------------------------------------------------------------------- sparks */

/**
 * Pyrophosphate. Every nucleotide added to a chain arrives as a triphosphate
 * and leaves two of those phosphates behind; hydrolysing them is what makes
 * the reaction go. It is the cell's receipt for the bond.
 */
export class Sparks {
  constructor(max = 200) {
    this.max = max;
    this.pos = new Float32Array(max * 3);
    this.vel = new Float32Array(max * 3);
    this.life = new Float32Array(max);
    this.head = 0;

    const geo = new THREE.BufferGeometry();
    geo.setAttribute('position', new THREE.BufferAttribute(this.pos, 3));
    geo.setAttribute('alpha', new THREE.BufferAttribute(this.life, 1));
    geo.boundingSphere = new THREE.Sphere(new THREE.Vector3(), 1e4);

    const mat = new THREE.ShaderMaterial({
      uniforms: { uSize: { value: 34 }, uColour: { value: new THREE.Color(0xffc978) } },
      vertexShader: `
        attribute float alpha;
        varying float vA;
        uniform float uSize;
        void main() {
          vA = alpha;
          vec4 mv = modelViewMatrix * vec4(position, 1.0);
          gl_PointSize = uSize * alpha / max(-mv.z, 1.0) * 60.0;
          gl_Position = projectionMatrix * mv;
        }`,
      fragmentShader: `
        varying float vA;
        uniform vec3 uColour;
        void main() {
          vec2 d = gl_PointCoord - 0.5;
          float r = dot(d, d);
          if (r > 0.25) discard;
          float g = smoothstep(0.25, 0.0, r);
          gl_FragColor = vec4(uColour * (1.0 + g * 2.0), g * vA);
        }`,
      transparent: true,
      depthWrite: false,
      blending: THREE.AdditiveBlending,
    });

    this.points = new THREE.Points(geo, mat);
    this.points.frustumCulled = false;
  }

  burst(at, n = 6) {
    for (let k = 0; k < n; k++) {
      const i = this.head++ % this.max;
      const o = i * 3;
      this.pos[o] = at[0]; this.pos[o + 1] = at[1]; this.pos[o + 2] = at[2];
      const d = new THREE.Vector3(rand(-1, 1), rand(-1, 1), rand(-1, 1)).normalize().multiplyScalar(rand(14, 46));
      this.vel[o] = d.x; this.vel[o + 1] = d.y; this.vel[o + 2] = d.z;
      this.life[i] = 1;
    }
  }

  update(dt) {
    for (let i = 0; i < this.max; i++) {
      if (this.life[i] <= 0) continue;
      const o = i * 3;
      this.pos[o] += this.vel[o] * dt;
      this.pos[o + 1] += this.vel[o + 1] * dt;
      this.pos[o + 2] += this.vel[o + 2] * dt;
      const drag = Math.exp(-2.4 * dt);
      this.vel[o] *= drag; this.vel[o + 1] *= drag; this.vel[o + 2] *= drag;
      this.life[i] = Math.max(0, this.life[i] - dt * 1.5);
    }
    this.points.geometry.attributes.position.needsUpdate = true;
    this.points.geometry.attributes.alpha.needsUpdate = true;
  }

  clear() {
    this.life.fill(0);
    this.points.geometry.attributes.alpha.needsUpdate = true;
  }
}
