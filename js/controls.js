/**
 * controls.js — orbit, pan and dolly, with damping.
 *
 * three's OrbitControls lives in the addons bundle; this project vendors only
 * the core build, so here is the eighty lines of it that this actually needs.
 */

import * as THREE from './lib/three.module.js';

const _v = new THREE.Vector3();
const _off = new THREE.Vector3();

export class Orbit {
  constructor(camera, dom, { target = new THREE.Vector3(), minDistance = 20, maxDistance = 900 } = {}) {
    this.camera = camera;
    this.dom = dom;
    this.target = target.clone();
    this.minDistance = minDistance;
    this.maxDistance = maxDistance;

    _off.copy(camera.position).sub(this.target);
    this.distance = _off.length();
    this.theta = Math.atan2(_off.x, _off.z);
    this.phi = Math.acos(Math.min(1, Math.max(-1, _off.y / this.distance)));

    this.goalTheta = this.theta;
    this.goalPhi = this.phi;
    this.goalDistance = this.distance;
    this.goalTarget = this.target.clone();

    this.damping = 0.12;
    this.autoRotate = 0;      // radians per second
    this.enabled = true;

    this._pointers = new Map();
    this._lastPinch = 0;
    this._mode = null;

    const el = dom;
    el.style.touchAction = 'none';

    el.addEventListener('pointerdown', (e) => {
      if (!this.enabled) return;
      el.setPointerCapture(e.pointerId);
      this._pointers.set(e.pointerId, { x: e.clientX, y: e.clientY });
      this._mode = (e.button === 2 || e.shiftKey) ? 'pan' : 'orbit';
    });

    el.addEventListener('pointermove', (e) => {
      const p = this._pointers.get(e.pointerId);
      if (!p) return;
      const dx = e.clientX - p.x;
      const dy = e.clientY - p.y;
      p.x = e.clientX; p.y = e.clientY;

      if (this._pointers.size === 2) {
        const [a, b] = [...this._pointers.values()];
        const d = Math.hypot(a.x - b.x, a.y - b.y);
        if (this._lastPinch) this.goalDistance *= this._lastPinch / d;
        this._lastPinch = d;
        this.goalDistance = Math.min(this.maxDistance, Math.max(this.minDistance, this.goalDistance));
        return;
      }

      if (this._mode === 'pan') this._pan(dx, dy);
      else {
        this.goalTheta -= dx * 0.0055;
        this.goalPhi = Math.min(Math.PI - 0.08, Math.max(0.08, this.goalPhi - dy * 0.0055));
      }
    });

    const end = (e) => {
      this._pointers.delete(e.pointerId);
      if (this._pointers.size < 2) this._lastPinch = 0;
      if (!this._pointers.size) this._mode = null;
    };
    el.addEventListener('pointerup', end);
    el.addEventListener('pointercancel', end);
    el.addEventListener('contextmenu', (e) => e.preventDefault());

    el.addEventListener('wheel', (e) => {
      if (!this.enabled) return;
      e.preventDefault();
      this.goalDistance *= Math.exp(e.deltaY * 0.0011);
      this.goalDistance = Math.min(this.maxDistance, Math.max(this.minDistance, this.goalDistance));
    }, { passive: false });
  }

  _pan(dx, dy) {
    const scale = this.distance * Math.tan((this.camera.fov / 2) * Math.PI / 180) * 2 / this.dom.clientHeight;
    const right = _v.setFromMatrixColumn(this.camera.matrix, 0);
    this.goalTarget.addScaledVector(right, -dx * scale);
    const up = _v.setFromMatrixColumn(this.camera.matrix, 1);
    this.goalTarget.addScaledVector(up, dy * scale);
  }

  /** Ease towards a framing. Used when the scene changes. */
  frame({ target, distance, theta, phi }, snap = false) {
    if (target) this.goalTarget.set(target[0], target[1], target[2]);
    if (distance !== undefined) this.goalDistance = distance;
    if (theta !== undefined) this.goalTheta = theta;
    if (phi !== undefined) this.goalPhi = phi;
    if (snap) {
      this.target.copy(this.goalTarget);
      this.distance = this.goalDistance;
      this.theta = this.goalTheta;
      this.phi = this.goalPhi;
    }
  }

  update(dt) {
    if (this.autoRotate) this.goalTheta += this.autoRotate * dt;

    const k = 1 - Math.pow(1 - this.damping, dt * 60);
    this.theta += (this.goalTheta - this.theta) * k;
    this.phi += (this.goalPhi - this.phi) * k;
    this.distance += (this.goalDistance - this.distance) * k;
    this.target.lerp(this.goalTarget, k);

    const sp = Math.sin(this.phi);
    this.camera.position.set(
      this.target.x + this.distance * sp * Math.sin(this.theta),
      this.target.y + this.distance * Math.cos(this.phi),
      this.target.z + this.distance * sp * Math.cos(this.theta),
    );
    this.camera.lookAt(this.target);
  }
}
