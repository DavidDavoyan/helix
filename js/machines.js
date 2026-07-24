/**
 * machines.js — the proteins.
 *
 * These are shapes, not structures: a helicase here is a lumpy ring, not the
 * DnaB hexamer. But each keeps the one feature that explains what the thing
 * does — the helicase is a ring because it encircles a single strand, the
 * polymerase has a cleft because the duplex has to sit in it, the ribosome is
 * two subunits because mRNA threads between them, and tRNA is an L because the
 * anticodon must be 7 nm from the amino acid it carries.
 *
 * Everything is a THREE.Group with a `place(pos, dir)` method. What moves them
 * is in processes.js.
 */

import * as THREE from './lib/three.module.js';
import { makeProteinBlob, roughen } from './geometry.js';
import { AMINO } from './bio.js';

const _v = new THREE.Vector3();
const _q = new THREE.Quaternion();
const _up = new THREE.Vector3(0, 0, 1);

/* ------------------------------------------------------------- appearance */

/**
 * Proteins read as soft, wet and slightly translucent — dense in the middle,
 * bright at the silhouette. That is a rim term, so the surface gets a normal
 * shell plus an additive fresnel envelope over it.
 */
function proteinMaterial(colour, opacity = 0.62) {
  return new THREE.MeshStandardMaterial({
    color: colour,
    roughness: 0.55,
    metalness: 0.05,
    transparent: true,
    opacity,
    depthWrite: false,
    envMapIntensity: 0.8,
  });
}

function fresnelShell(colour, power = 2.4, strength = 0.55) {
  return new THREE.ShaderMaterial({
    uniforms: {
      uColour: { value: new THREE.Color(colour) },
      uPower: { value: power },
      uStrength: { value: strength },
    },
    vertexShader: `
      varying vec3 vN;
      varying vec3 vV;
      void main() {
        vec4 mv = modelViewMatrix * vec4(position, 1.0);
        vN = normalize(normalMatrix * normal);
        vV = normalize(-mv.xyz);
        gl_Position = projectionMatrix * mv;
      }`,
    fragmentShader: `
      uniform vec3 uColour;
      uniform float uPower;
      uniform float uStrength;
      varying vec3 vN;
      varying vec3 vV;
      void main() {
        float f = pow(1.0 - clamp(dot(normalize(vN), normalize(vV)), 0.0, 1.0), uPower);
        gl_FragColor = vec4(uColour * f * uStrength, f * uStrength);
      }`,
    transparent: true,
    depthWrite: false,
    blending: THREE.AdditiveBlending,
    side: THREE.FrontSide,
  });
}

/** A body plus its glow envelope, as one object. */
function shelled(geo, colour, opacity, rim = 1.03) {
  const g = new THREE.Group();
  const body = new THREE.Mesh(geo, proteinMaterial(colour, opacity));
  const shell = new THREE.Mesh(geo, fresnelShell(colour));
  shell.scale.setScalar(rim);
  g.add(body, shell);
  g.userData.body = body;
  g.userData.shell = shell;
  return g;
}

/* ------------------------------------------------------------ base machine */

class Machine extends THREE.Group {
  constructor(name, caption) {
    super();
    this.name = name;
    this.caption = caption;
    this.visible = false;
  }

  /** Put the machine at `pos`, facing along `dir` (its local +z). */
  place(pos, dir) {
    this.position.set(pos[0], pos[1], pos[2]);
    if (dir) {
      _v.set(dir[0], dir[1], dir[2]).normalize();
      _q.setFromUnitVectors(_up, _v);
      this.quaternion.copy(_q);
    }
  }

  setTint(hex) {
    this.traverse((o) => {
      if (o.material && o.material.color) o.material.color.setHex(hex);
      if (o.material && o.material.uniforms && o.material.uniforms.uColour) {
        o.material.uniforms.uColour.value.setHex(hex);
      }
    });
  }
}

/* ---------------------------------------------------------------- helicase */

/**
 * Replicative helicase. A ring, because it works by threading one strand
 * through its middle and shouldering the other aside — it does not so much cut
 * the base pairs as walk between them.
 */
export function makeHelicase() {
  const m = new Machine('helicase', 'Helicase — breaking the hydrogen bonds, two at a time on A·T, three on G·C');

  const ring = roughen(new THREE.TorusGeometry(14, 6.2, 14, 26), 0.16, 2.2, 3);
  ring.scale(1, 1, 0.72);
  const body = shelled(ring, 0x6ad2ff, 0.5);
  m.add(body);

  // Two lobes on the leading face: the wedge that splits the strands.
  for (const s of [-1, 1]) {
    const lobe = shelled(makeProteinBlob({ radius: 7.4, detail: 2, lumpiness: 0.3, seed: 5 + s }), 0x53b6e8, 0.55);
    lobe.position.set(0, s * 11, 9.5);
    m.add(lobe);
  }
  return m;
}

/* -------------------------------------------------------------- polymerase */

/**
 * DNA polymerase. The classic description is a right hand: a palm holding the
 * catalytic site, fingers that close over an incoming nucleotide, and a thumb
 * gripping the duplex behind. The cleft between them is where the template
 * runs, so the cleft is the part that has to be legible.
 */
export function makePolymerase(colour = 0xffb454, scale = 1) {
  const m = new Machine('polymerase', 'DNA polymerase — adds nucleotides only to a 3′ end, so one strand runs backwards');

  const palm = shelled(makeProteinBlob({
    radius: 13 * scale, detail: 3, lumpiness: 0.26, seed: 1,
    scale: [1.15, 0.92, 1.0], channel: { axis: 'z', radius: 6.5 * scale },
  }), colour, 0.55);
  m.add(palm);

  const fingers = shelled(makeProteinBlob({ radius: 8.2 * scale, detail: 2, lumpiness: 0.3, seed: 2 }), colour, 0.5);
  fingers.position.set(0, 12 * scale, 3 * scale);
  m.add(fingers);
  m.userData.fingers = fingers;

  const thumb = shelled(makeProteinBlob({ radius: 7 * scale, detail: 2, lumpiness: 0.32, seed: 4 }), colour, 0.5);
  thumb.position.set(0, -11 * scale, -4 * scale);
  m.add(thumb);

  return m;
}

/**
 * RNA polymerase. Bigger, and it does its own unwinding — it opens a bubble of
 * about 14 base pairs, reads one strand, and pushes the transcript out through
 * a separate channel so the DNA can zip up again behind it.
 */
export function makeRnaPolymerase() {
  const m = new Machine('rnap', 'RNA polymerase — melts a 14 bp bubble, reads the template, and rewinds behind itself');

  // Kept notably see-through: the whole point of the scene is the melted
  // bubble inside it, and an opaque enzyme hides exactly what matters.
  const core = shelled(makeProteinBlob({
    radius: 17, detail: 3, lumpiness: 0.24, seed: 7,
    scale: [1.0, 1.05, 0.95], channel: { axis: 'z', radius: 8 },
  }), 0xc39bff, 0.3);
  m.add(core);

  const clamp = shelled(makeProteinBlob({ radius: 9.5, detail: 2, lumpiness: 0.28, seed: 9 }), 0xa87fe8, 0.5);
  clamp.position.set(0, 14, 4);
  m.add(clamp);

  // The RNA exit channel, pointing away from the DNA.
  const exit = shelled(roughen(new THREE.CylinderGeometry(4.2, 6.2, 14, 12), 0.1, 3, 2), 0xb98cf5, 0.45);
  exit.rotation.z = Math.PI / 2;
  exit.position.set(-16, -6, 0);
  m.add(exit);
  m.userData.exit = exit;

  return m;
}

/** Single-strand binding protein: small, and it stops the melted strands re-pairing. */
export function makeSSB() {
  const m = new Machine('ssb', 'Single-strand binding protein — keeps the open strands from snapping shut');
  const b = shelled(makeProteinBlob({ radius: 5.4, detail: 2, lumpiness: 0.34, seed: 12 }), 0x7fe3c4, 0.5);
  m.add(b);
  return m;
}

/** Primase lays the short RNA primer every polymerase needs to start from. */
export function makePrimase() {
  const m = new Machine('primase', 'Primase — polymerase cannot start a chain, so this lays a short RNA primer first');
  const b = shelled(makeProteinBlob({ radius: 8.2, detail: 2, lumpiness: 0.3, seed: 15, scale: [1, 0.85, 1] }), 0xff87c2, 0.52);
  m.add(b);
  return m;
}

/** Ligase seals the nick left between two Okazaki fragments. */
export function makeLigase() {
  const m = new Machine('ligase', 'DNA ligase — seals the nick between fragments into one continuous strand');
  const b = shelled(makeProteinBlob({ radius: 7.6, detail: 2, lumpiness: 0.28, seed: 21, scale: [1.1, 1, 0.9] }), 0x9de84f, 0.52);
  m.add(b);
  return m;
}

/* ---------------------------------------------------------------- ribosome */

/**
 * The ribosome, as two subunits with a slot between them.
 *
 * mRNA runs through that slot. The large subunit carries three tRNA sites —
 * A where a charged tRNA arrives, P where the growing chain sits, E where a
 * spent tRNA leaves — and the peptide comes out of a tunnel in its back.
 */
export function makeRibosome() {
  const m = new Machine('ribosome', 'Ribosome — three tRNA sites: arrive at A, hold the chain at P, leave from E');

  const large = shelled(makeProteinBlob({
    radius: 26, detail: 3, lumpiness: 0.2, seed: 31, scale: [1.0, 0.9, 1.05],
  }), 0x8fa4c9, 0.46);
  large.position.set(0, 17, 0);
  m.add(large);
  m.userData.large = large;

  const small = shelled(makeProteinBlob({
    radius: 19, detail: 3, lumpiness: 0.22, seed: 37, scale: [1.05, 0.75, 1.1],
  }), 0x6f86ad, 0.46);
  small.position.set(0, -14, 0);
  m.add(small);
  m.userData.small = small;

  // Local anchors, in the ribosome's own frame. +z is the direction of travel
  // along the mRNA, so E is behind, A is ahead.
  m.userData.sites = {
    E: new THREE.Vector3(0, 6, -11.5),
    P: new THREE.Vector3(0, 6, 0),
    A: new THREE.Vector3(0, 6, 11.5),
    exit: new THREE.Vector3(6, 34, -6),
  };
  return m;
}

/* -------------------------------------------------------------------- tRNA */

/**
 * Transfer RNA, in its real L shape.
 *
 * The cloverleaf everyone draws is the secondary structure; folded up it is an
 * L about 7 nm on a side, with the anticodon at one tip and the amino acid at
 * the other. That distance is the whole point — the codon is read in one place
 * and the peptide bond is made in another.
 */
export function makeTRNA(aa = 'M', anticodon = 'UAC') {
  const m = new Machine('trna', 'tRNA — anticodon at one end, amino acid at the other, 7 nm apart');
  m.userData.aa = aa;
  m.userData.anticodon = anticodon;

  // Two arms meeting at an elbow. The anticodon arm drops to the message; the
  // acceptor arm rises away from it, so the amino acid is carried up into the
  // large subunit where the peptide bond is made.
  const path = new THREE.CatmullRomCurve3([
    new THREE.Vector3(0, -24, 0),      // anticodon tip, down at the codon
    new THREE.Vector3(1.5, -15, 1),
    new THREE.Vector3(0, -6, 0),
    new THREE.Vector3(1, 1, -1),       // the elbow
    new THREE.Vector3(7, 6, 0),
    new THREE.Vector3(14, 12, 1),
    new THREE.Vector3(21, 17, 0),      // acceptor stem
  ]);
  const tube = new THREE.TubeGeometry(path, 44, 2.9, 8, false);
  const body = shelled(tube, 0xffd88a, 0.68, 1.06);
  m.add(body);

  // Anticodon: three bases at the tip, colour-coded like every other base.
  const tip = new THREE.Group();
  tip.position.set(0, -24, 0);
  m.add(tip);
  m.userData.tip = tip;

  // The amino acid it is carrying.
  const info = AMINO[aa] || AMINO.G;
  const res = new THREE.Mesh(
    new THREE.IcosahedronGeometry(5.2, 2),
    new THREE.MeshStandardMaterial({
      color: info.colour, roughness: 0.35, metalness: 0.1,
      emissive: new THREE.Color(info.colour).multiplyScalar(0.18),
    }),
  );
  res.position.set(24.5, 20, 0);
  m.add(res);
  m.userData.residue = res;

  return m;
}

/* ------------------------------------------------------------------ labels */

/**
 * A flat text sprite drawn to a canvas. Used for the tags that hang off each
 * machine so the viewer knows what they are watching without reading the panel.
 */
export function makeLabel(text, { colour = '#dbe6f5', size = 44, pad = 16 } = {}) {
  const canvas = document.createElement('canvas');
  const ctx = canvas.getContext('2d');
  const font = `600 ${size}px "Inter", "Segoe UI", system-ui, sans-serif`;
  ctx.font = font;
  const w = Math.ceil(ctx.measureText(text).width) + pad * 2;
  const h = size + pad * 2;
  canvas.width = w; canvas.height = h;

  ctx.font = font;
  ctx.fillStyle = 'rgba(8, 13, 22, 0.72)';
  ctx.beginPath();
  const r = 12;
  ctx.moveTo(r, 0); ctx.arcTo(w, 0, w, h, r); ctx.arcTo(w, h, 0, h, r);
  ctx.arcTo(0, h, 0, 0, r); ctx.arcTo(0, 0, w, 0, r); ctx.closePath(); ctx.fill();
  ctx.strokeStyle = 'rgba(140, 170, 210, 0.28)'; ctx.lineWidth = 2; ctx.stroke();
  ctx.fillStyle = colour;
  ctx.textBaseline = 'middle';
  ctx.fillText(text, pad, h / 2 + 2);

  const tex = new THREE.CanvasTexture(canvas);
  tex.colorSpace = THREE.SRGBColorSpace;
  tex.anisotropy = 4;
  const sprite = new THREE.Sprite(new THREE.SpriteMaterial({
    map: tex, transparent: true, depthTest: false, depthWrite: false,
  }));
  sprite.scale.set(w * 0.095, h * 0.095, 1);
  sprite.renderOrder = 999;
  return sprite;
}
