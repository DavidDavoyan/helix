/**
 * geometry.js — molecular building blocks.
 *
 * The shapes here are schematic but not arbitrary: rings are built at the
 * 1.39 A aromatic bond length, purines are a genuinely fused bicycle sharing
 * an edge, and the deoxyribose is a five-membered ring. At the scale we draw
 * them this reads as "two rings versus one", which is the distinction that
 * actually matters for how a base pair works.
 */

import * as THREE from './lib/three.module.js';

const BOND = 1.39; // A, aromatic C-C

/* ------------------------------------------------------------ ring outlines */

/**
 * Vertices of a regular n-gon of the given side length, rotated so that one
 * edge is vertical and sits on the +x side of the centre. That edge is the one
 * a second ring fuses onto.
 */
function polygon(n, side, cx = 0, cy = 0) {
  const R = side / (2 * Math.sin(Math.PI / n));
  const step = (2 * Math.PI) / n;
  const start = -step / 2; // puts two vertices symmetric about +x
  const pts = [];
  for (let i = 0; i < n; i++) {
    const a = start + i * step;
    pts.push([cx + R * Math.cos(a), cy + R * Math.sin(a)]);
  }
  return { pts, R, apothem: side / (2 * Math.tan(Math.PI / n)) };
}

/**
 * Outline of a pyrimidine: one six-membered ring. The glycosidic nitrogen (N1)
 * sits at the origin and the ring extends along +x, so a base plate can simply
 * be placed at C1' and pointed at its partner.
 */
function pyrimidineOutline() {
  const hex = polygon(6, BOND);
  // Shift so the leftmost vertex — N1 — lands on the origin.
  const minX = Math.min(...hex.pts.map((p) => p[0]));
  return hex.pts.map(([x, y]) => [x - minX, y]);
}

/**
 * Outline of a purine: a five-membered ring fused to a six-membered ring along
 * a shared edge, which is the imidazole + pyrimidine bicycle of adenine and
 * guanine. N9 — where the sugar attaches — is a vertex of the five-ring on the
 * far side from the fusion, so the glycosidic bond leaves at the
 * characteristic angle rather than straight down the long axis.
 */
function purineOutline() {
  const fuseX = 0; // the shared edge sits on x = 0
  const half = BOND / 2;

  const pent = polygon(5, BOND, fuseX - polygon(5, BOND).apothem, 0);
  const hex = polygon(6, BOND, fuseX + polygon(6, BOND).apothem, 0);

  // The two shared vertices, common to both rings.
  const shared = [
    [fuseX, +half],
    [fuseX, -half],
  ];
  const isShared = (p) => shared.some((s) => Math.hypot(s[0] - p[0], s[1] - p[1]) < 1e-6);

  // Walk the outline: pentagon vertices (minus the shared edge), then the
  // hexagon's, so the result is a single simple polygon.
  const pentOuter = pent.pts.filter((p) => !isShared(p));
  const hexOuter = hex.pts.filter((p) => !isShared(p));

  const ring = [
    shared[1],
    ...hexOuter.sort((a, b) => Math.atan2(a[1] - hex.pts[0][1], a[0]) - Math.atan2(b[1], b[0])),
    shared[0],
    ...pentOuter.reverse(),
  ];

  // Order robustly by angle about the bicycle's centroid rather than trusting
  // the hand-assembled walk above.
  const cx = ring.reduce((s, p) => s + p[0], 0) / ring.length;
  const cy = ring.reduce((s, p) => s + p[1], 0) / ring.length;
  ring.sort((a, b) => Math.atan2(a[1] - cy, a[0] - cx) - Math.atan2(b[1] - cy, b[0] - cx));

  // N9: the pentagon vertex furthest from the hexagon, below the long axis.
  const n9 = pentOuter.reduce((best, p) => (p[0] < best[0] ? p : best), pentOuter[0]);
  return ring.map(([x, y]) => [x - n9[0], y - n9[1]]);
}

/* -------------------------------------------------------------- base plates */

function plateFromOutline(outline, thickness) {
  const shape = new THREE.Shape();
  shape.moveTo(outline[0][0], outline[0][1]);
  for (let i = 1; i < outline.length; i++) shape.lineTo(outline[i][0], outline[i][1]);
  shape.closePath();

  const geo = new THREE.ExtrudeGeometry(shape, {
    depth: thickness,
    bevelEnabled: true,
    bevelThickness: thickness * 0.3,
    bevelSize: 0.07,
    bevelSegments: 1,
    curveSegments: 1,
  });
  // Extrusion runs along +z; centre it so the plate straddles its base-pair
  // plane, then leave +x pointing from the glycosidic nitrogen at the partner.
  geo.translate(0, 0, -thickness / 2);
  geo.computeVertexNormals();
  return geo;
}

/**
 * The two base plates, in a local frame where the origin is the glycosidic
 * nitrogen, +x points at the partner base and +z is the helix axis.
 */
export function makeBaseGeometries(thickness = 1.25) {
  return {
    fused: plateFromOutline(purineOutline(), thickness),
    single: plateFromOutline(pyrimidineOutline(), thickness),
  };
}

/** Deoxyribose: a puckered five-ring, drawn flat. */
export function makeSugarGeometry(thickness = 1.0) {
  const pent = polygon(5, 1.45);
  return plateFromOutline(pent.pts, thickness);
}

/** Phosphate: a tetrahedral PO4, near enough a small tetrahedron. */
export function makePhosphateGeometry(r = 1.1) {
  const geo = new THREE.TetrahedronGeometry(r, 0);
  geo.computeVertexNormals();
  return geo;
}

/* -------------------------------------------------------------- ribbon tube */

/**
 * A tube of fixed topology whose vertices are rewritten every frame.
 *
 * Rebuilding a TubeGeometry per frame allocates; this does not. Frames are
 * carried along the curve by parallel transport (rotation-minimising), because
 * Frenet frames flip at inflection points and a DNA backbone is nothing but
 * inflection points.
 */
export class Ribbon {
  /**
   * @param {number} maxPoints  most spine points that will ever be supplied
   * @param {number} radial     cross-section resolution
   */
  constructor(maxPoints, radius, radial = 8, { vertexColours = false } = {}) {
    this.maxPoints = maxPoints;
    this.radius = radius;
    this.radial = radial;
    this.count = 0;

    const verts = maxPoints * radial;
    this.position = new Float32Array(verts * 3);
    this.normal = new Float32Array(verts * 3);
    this.colour = vertexColours ? new Float32Array(verts * 3) : null;

    const index = new Uint32Array((maxPoints - 1) * radial * 6);
    for (let s = 0; s < maxPoints - 1; s++) {
      for (let r = 0; r < radial; r++) {
        const a = s * radial + r;
        const b = s * radial + ((r + 1) % radial);
        const c = (s + 1) * radial + r;
        const d = (s + 1) * radial + ((r + 1) % radial);
        const o = (s * radial + r) * 6;
        index[o + 0] = a; index[o + 1] = c; index[o + 2] = b;
        index[o + 3] = b; index[o + 4] = c; index[o + 5] = d;
      }
    }

    this.geometry = new THREE.BufferGeometry();
    this.geometry.setAttribute('position', new THREE.BufferAttribute(this.position, 3));
    this.geometry.setAttribute('normal', new THREE.BufferAttribute(this.normal, 3));
    if (this.colour) this.geometry.setAttribute('color', new THREE.BufferAttribute(this.colour, 3));
    this.geometry.setIndex(new THREE.BufferAttribute(index, 1));
    this.geometry.setDrawRange(0, 0);
    this.geometry.boundingSphere = new THREE.Sphere(new THREE.Vector3(), 1e4);

    this._tangent = new THREE.Vector3();
    this._up = new THREE.Vector3();
    this._side = new THREE.Vector3();
    this._prevUp = new THREE.Vector3();
  }

  /**
   * @param {THREE.Vector3[]} pts   spine, at least two points
   * @param {number}          count how many of them to use
   * @param {number[]}        [taper] optional per-point radius multiplier
   * @param {THREE.Color[]}   [colours] optional per-point colour
   */
  update(pts, count = pts.length, taper = null, colours = null) {
    count = Math.min(count, this.maxPoints);
    this.count = count;
    if (count < 2) {
      this.geometry.setDrawRange(0, 0);
      return;
    }

    const { radial, position, normal } = this;
    const t = this._tangent;
    const up = this._up;
    const side = this._side;
    const prevUp = this._prevUp;

    // Seed an up-vector that is not parallel to the first tangent.
    t.copy(pts[1]).sub(pts[0]).normalize();
    prevUp.set(0, 0, 1);
    if (Math.abs(prevUp.dot(t)) > 0.9) prevUp.set(1, 0, 0);
    prevUp.projectOnPlane(t).normalize();

    for (let i = 0; i < count; i++) {
      // Central-difference tangent, forward/backward at the ends.
      if (i === 0) t.copy(pts[1]).sub(pts[0]);
      else if (i === count - 1) t.copy(pts[i]).sub(pts[i - 1]);
      else t.copy(pts[i + 1]).sub(pts[i - 1]);
      if (t.lengthSq() < 1e-12) t.set(0, 0, 1);
      t.normalize();

      // Parallel transport: keep the previous up, remove any component that has
      // fallen along the new tangent. No flips, no accumulated twist.
      up.copy(prevUp).projectOnPlane(t);
      if (up.lengthSq() < 1e-8) {
        up.set(0, 0, 1).projectOnPlane(t);
        if (up.lengthSq() < 1e-8) up.set(1, 0, 0).projectOnPlane(t);
      }
      up.normalize();
      prevUp.copy(up);
      side.crossVectors(t, up).normalize();

      const r = this.radius * (taper ? taper[i] : 1);
      const p = pts[i];
      const col = this.colour && colours ? colours[i] : null;
      for (let k = 0; k < radial; k++) {
        const a = (k / radial) * Math.PI * 2;
        const nx = up.x * Math.cos(a) + side.x * Math.sin(a);
        const ny = up.y * Math.cos(a) + side.y * Math.sin(a);
        const nz = up.z * Math.cos(a) + side.z * Math.sin(a);
        const o = (i * radial + k) * 3;
        position[o + 0] = p.x + nx * r;
        position[o + 1] = p.y + ny * r;
        position[o + 2] = p.z + nz * r;
        normal[o + 0] = nx;
        normal[o + 1] = ny;
        normal[o + 2] = nz;
        if (col) {
          this.colour[o + 0] = col.r;
          this.colour[o + 1] = col.g;
          this.colour[o + 2] = col.b;
        }
      }
    }

    this.geometry.attributes.position.needsUpdate = true;
    this.geometry.attributes.normal.needsUpdate = true;
    if (this.colour && colours) this.geometry.attributes.color.needsUpdate = true;
    this.geometry.setDrawRange(0, (count - 1) * radial * 6);
  }

  dispose() {
    this.geometry.dispose();
  }
}

/* ------------------------------------------------------------------- blobs */

/** Cheap 3D value noise — enough to make a protein surface look lumpy. */
function hash3(x, y, z) {
  const s = Math.sin(x * 127.1 + y * 311.7 + z * 74.7) * 43758.5453;
  return s - Math.floor(s);
}

function vnoise(x, y, z) {
  const xi = Math.floor(x), yi = Math.floor(y), zi = Math.floor(z);
  const xf = x - xi, yf = y - yi, zf = z - zi;
  const u = xf * xf * (3 - 2 * xf), v = yf * yf * (3 - 2 * yf), w = zf * zf * (3 - 2 * zf);
  const lerp = (a, b, t) => a + (b - a) * t;
  const c = (i, j, k) => hash3(xi + i, yi + j, zi + k);
  return lerp(
    lerp(lerp(c(0, 0, 0), c(1, 0, 0), u), lerp(c(0, 1, 0), c(1, 1, 0), u), v),
    lerp(lerp(c(0, 0, 1), c(1, 0, 1), u), lerp(c(0, 1, 1), c(1, 1, 1), u), v),
    w,
  );
}

function fbm(x, y, z, octaves = 3) {
  let a = 0.5, f = 1, sum = 0, norm = 0;
  for (let i = 0; i < octaves; i++) {
    sum += a * vnoise(x * f, y * f, z * f);
    norm += a;
    a *= 0.5; f *= 2.05;
  }
  return sum / norm;
}

/**
 * A globular protein: a noise-displaced icosphere, optionally squashed and
 * optionally bored through by a channel so DNA or mRNA can be threaded.
 *
 * Real enzymes are not spheres, but they are compact, lumpy and roughly
 * globular at this scale, and a channel is the one feature that has to read
 * clearly — the viewer has to see the nucleic acid go *through* the machine.
 */
export function makeProteinBlob({
  radius = 10,
  detail = 3,
  lumpiness = 0.28,
  scale = [1, 1, 1],
  seed = 0,
  channel = null, // { axis: 'x'|'y'|'z', radius: n } — carve a groove
} = {}) {
  const geo = new THREE.IcosahedronGeometry(radius, detail);
  const pos = geo.attributes.position;
  const v = new THREE.Vector3();

  for (let i = 0; i < pos.count; i++) {
    v.fromBufferAttribute(pos, i);
    const n = v.clone().normalize();
    const d = fbm(n.x * 2.1 + seed * 13.7, n.y * 2.1 + seed * 7.3, n.z * 2.1 + seed * 3.1, 3);
    let r = radius * (1 + (d - 0.5) * 2 * lumpiness);

    if (channel) {
      // Distance from the channel axis, in units of the blob radius.
      const ax = channel.axis;
      const off = ax === 'x' ? Math.hypot(n.y, n.z) : ax === 'y' ? Math.hypot(n.x, n.z) : Math.hypot(n.x, n.y);
      const bite = 1 - Math.min(1, off / (channel.radius / radius));
      r -= bite * radius * 0.55;
    }

    v.copy(n).multiplyScalar(r);
    v.x *= scale[0]; v.y *= scale[1]; v.z *= scale[2];
    pos.setXYZ(i, v.x, v.y, v.z);
  }

  geo.computeVertexNormals();
  return geo;
}

/**
 * Push a geometry's vertices along their normals by noise. Turns the smooth
 * primitives that helicases and ribosomes start life as into something with
 * the lumpy, folded look of a real protein surface.
 */
export function roughen(geo, amount = 0.12, freq = 2.4, seed = 0) {
  const pos = geo.attributes.position;
  const nor = geo.attributes.normal;
  const v = new THREE.Vector3();
  const n = new THREE.Vector3();
  for (let i = 0; i < pos.count; i++) {
    v.fromBufferAttribute(pos, i);
    n.fromBufferAttribute(nor, i);
    const d = fbm(v.x * freq * 0.06 + seed * 11.3, v.y * freq * 0.06 + seed * 5.9, v.z * freq * 0.06 + seed * 2.7, 3);
    const s = (d - 0.5) * 2 * amount;
    pos.setXYZ(i, v.x + n.x * s * 10, v.y + n.y * s * 10, v.z + n.z * s * 10);
  }
  geo.computeVertexNormals();
  return geo;
}

/**
 * A unit cylinder along +y, radius 1, height 1, origin at the base. Instanced
 * and stretched into hydrogen bonds and backbone links.
 */
export function makeUnitBond(radial = 6) {
  const geo = new THREE.CylinderGeometry(1, 1, 1, radial, 1, true);
  geo.translate(0, 0.5, 0);
  return geo;
}
