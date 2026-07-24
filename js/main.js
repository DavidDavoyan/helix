/**
 * main.js — boot, lighting, the frame loop, and the wiring between the UI and
 * whichever scene is running.
 */

import * as THREE from './lib/three.module.js';
import { DEFAULT_SEQUENCE, PRESETS, cleanSequence, BDNA } from './bio.js';
import { Stage, zOf } from './stage.js';
import { StructureScene, ReplicationScene } from './scene-dna.js';
import { TranscriptionScene, TranslationScene } from './scene-expression.js';
import { ProteinScene } from './scene-protein.js';
import { Orbit } from './controls.js';
import { Composer } from './post.js';
import { UI } from './ui.js';

/* ---------------------------------------------------------------- quality */

const coarse = matchMedia('(pointer: coarse)').matches;
const small = innerWidth < 900;
const LOW = coarse || small;
const CAPTURE = new URLSearchParams(location.search).has('capture');

const CFG = {
  pixelRatio: Math.min(devicePixelRatio || 1, LOW ? 1.5 : 1.9),
  quality: LOW ? 'low' : 'high',
  poolCount: LOW ? 80 : 140,
};

/* ------------------------------------------------------------------- boot */

const canvas = document.getElementById('stage');
const renderer = new THREE.WebGLRenderer({
  canvas,
  antialias: false,
  stencil: false,
  powerPreference: 'high-performance',
  preserveDrawingBuffer: CAPTURE,
});
renderer.setPixelRatio(CFG.pixelRatio);
renderer.setSize(Math.max(1, innerWidth), Math.max(1, innerHeight), false);
renderer.outputColorSpace = THREE.LinearSRGBColorSpace;  // the composer encodes

const scene = new THREE.Scene();
scene.fog = new THREE.FogExp2(0x050810, 0.0019);

const camera = new THREE.PerspectiveCamera(42, Math.max(1, innerWidth) / Math.max(1, innerHeight), 1, 4000);
camera.position.set(120, 70, 190);

const composer = new Composer(renderer, { quality: CFG.quality });
const orbit = new Orbit(camera, canvas, { minDistance: 30, maxDistance: 1100 });

/* ----------------------------------------------------------------- lights */

/**
 * A three-point rig plus an image-based environment. The environment is what
 * makes the base plates and the phosphates read as solid objects rather than
 * flat colour, and it is generated here rather than loaded so the folder stays
 * dependency-free.
 */
function buildEnvironment() {
  const env = new THREE.Scene();

  const panel = (colour, intensity, pos, scale) => {
    const m = new THREE.Mesh(
      new THREE.PlaneGeometry(1, 1),
      new THREE.MeshBasicMaterial({ color: new THREE.Color(colour).multiplyScalar(intensity), side: THREE.DoubleSide }),
    );
    m.position.set(pos[0], pos[1], pos[2]);
    m.scale.set(scale[0], scale[1], 1);
    m.lookAt(0, 0, 0);
    env.add(m);
    return m;
  };

  env.add(new THREE.Mesh(
    new THREE.SphereGeometry(60, 16, 12),
    new THREE.MeshBasicMaterial({ color: 0x0a1220, side: THREE.BackSide }),
  ));
  panel(0xbcd8ff, 2.6, [0, 42, 8], [60, 40]);      // key, cool, from above
  panel(0x3d6fd0, 1.5, [-40, -6, -22], [46, 44]);  // deep blue fill
  panel(0xff9a52, 1.0, [34, -12, 20], [30, 30]);   // warm bounce
  panel(0x63d3ff, 0.9, [10, 4, -46], [40, 26]);    // rim

  const pmrem = new THREE.PMREMGenerator(renderer);
  pmrem.compileEquirectangularShader();
  const rt = pmrem.fromScene(env, 0.03);
  pmrem.dispose();
  env.traverse((o) => { if (o.geometry) o.geometry.dispose(); });
  return rt.texture;
}

const envMap = buildEnvironment();
scene.environment = envMap;

/**
 * A backdrop, so the molecule sits in something rather than in nothing. Cool
 * at the top, a hint of warmth low down, and a soft pool of light behind the
 * subject to separate the silhouette from the background.
 */
scene.add((() => {
  const dome = new THREE.Mesh(
    new THREE.SphereGeometry(2200, 32, 24),
    new THREE.ShaderMaterial({
      uniforms: {
        uTop: { value: new THREE.Color(0x070c16) },
        uBottom: { value: new THREE.Color(0x0d1018) },
        uGlow: { value: new THREE.Color(0x14304e) },
      },
      vertexShader: `
        varying vec3 vP;
        void main() {
          vP = normalize(position);
          gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
        }`,
      fragmentShader: `
        uniform vec3 uTop, uBottom, uGlow;
        varying vec3 vP;
        void main() {
          float h = vP.y * 0.5 + 0.5;
          vec3 c = mix(uBottom, uTop, smoothstep(0.0, 1.0, h));
          // A broad pool of light behind and above the subject.
          c += uGlow * pow(max(0.0, dot(vP, normalize(vec3(-0.25, 0.55, -0.8)))), 3.0);
          gl_FragColor = vec4(c, 1.0);
        }`,
      side: THREE.BackSide,
      depthWrite: false,
      fog: false,
    }),
  );
  dome.renderOrder = -1;
  dome.frustumCulled = false;
  return dome;
})());

const key = new THREE.DirectionalLight(0xdcecff, 1.5);
key.position.set(80, 120, 90);
scene.add(key);

const fill = new THREE.DirectionalLight(0x5f8fd6, 0.7);
fill.position.set(-110, -30, -40);
scene.add(fill);

const rim = new THREE.DirectionalLight(0xffb072, 0.55);
rim.position.set(-30, 40, -140);
scene.add(rim);

scene.add(new THREE.AmbientLight(0x243349, 0.55));

/* ------------------------------------------------------------------ world */

const stage = new Stage(scene);
stage.renderer.setEnvironment(envMap);
stage.peptide.setEnvironment(envMap);

const scenes = {
  structure: new StructureScene(stage),
  replication: new ReplicationScene(stage),
  transcription: new TranscriptionScene(stage),
  translation: new TranslationScene(stage),
  protein: new ProteinScene(stage),
};

let current = null;
let currentId = null;
let playing = true;
let speed = 1;
let sequence = DEFAULT_SEQUENCE;

/* ----------------------------------------------------------------- camera */

/**
 * Frame the molecule for whichever scene has just started.
 *
 * Everything is built along +z, so theta near zero puts the camera on the axis
 * and the whole molecule collapses to a disc. Side-on is theta near pi/2.
 */
function frameScene(id, snap = false) {
  const n = stage.n;
  const span = n * BDNA.RISE;
  switch (id) {
    case 'structure':
      orbit.frame({ target: [0, 0, 0], distance: Math.max(120, span * 0.78), theta: 1.35, phi: 1.30 }, snap);
      break;
    case 'replication':
      orbit.frame({ target: [0, 0, 0], distance: Math.max(170, span * 0.95), theta: 1.48, phi: 1.24 }, snap);
      break;
    case 'transcription':
      orbit.frame({ target: [0, -4, 0], distance: Math.max(150, span * 0.85), theta: 1.72, phi: 1.26 }, snap);
      break;
    case 'translation':
      orbit.frame({ target: [0, 16, 0], distance: 165, theta: 1.52, phi: 1.30 }, snap);
      break;
    case 'protein':
      // Starts wide on the extended chain; the follow below pulls in as it
      // collapses, because a folded 36-mer is a twentieth the length.
      orbit.frame({ target: [0, 0, 0], distance: Math.max(120, span * 0.55), theta: 1.1, phi: 1.24 }, snap);
      break;
    default: break;
  }
}

function setScene(id) {
  if (current) current.exit();
  current = scenes[id];
  currentId = id;
  current.enter();
  ui.setMode(id);
  frameScene(id);
  orbit.autoRotate = ui._toggleEls.spin.input.checked ? 0.16 : 0;
}

function restart() {
  if (!current) return;
  current.enter();
  frameScene(currentId);
}

/* --------------------------------------------------------------------- UI */

const ui = new UI({
  onMode: (id) => setScene(id),
  onSequence: (seq, isNew) => {
    const clean = cleanSequence(seq);
    if (!clean.length) return;
    sequence = clean;
    stage.setSequence(clean);
    ui.setSequence(clean, { isNew });
    restart();
  },
  onPlayToggle: () => { playing = !playing; ui.setPlaying(playing); },
  onRestart: () => restart(),
  onSpeed: (v) => { speed = v; },
  // One slider, two meanings: it melts base pairs in the DNA view and
  // denatures the fold in the protein view. Both are the same idea — heat
  // against many weak bonds — so they share a control.
  onTemperature: (v) => {
    scenes.structure.temperature = v;
    scenes.protein.temperature = v;
  },
  onToggle: (k, v) => {
    const r = stage.renderer;
    if (k === 'bases') r.showBases = v;
    if (k === 'backbone') r.showBackbone = v;
    if (k === 'sugars') r.showSugars = v;
    if (k === 'hbonds') r.showHBonds = v;
    if (k === 'grooves') scenes.structure.showGrooves = v;
    if (k === 'pool') {
      stage.pool.setVisible(v && !['translation', 'structure', 'protein'].includes(currentId));
    }
    if (k === 'labels') stage.showLabels = v;
    if (k === 'spin') orbit.autoRotate = v ? 0.16 : 0;
  },
});

/* ------------------------------------------------------------------- loop */

const clock = new THREE.Clock();
const _follow = new THREE.Vector3();
let elapsed = 0;

function step(dt) {
  if (!current) return;

  if (playing) {
    current.update(dt, speed);
    stage.pool.update(dt, speed);
  } else {
    // Paused still needs geometry written, or the first frame after a mode
    // change draws nothing.
    current.update(0, 0);
  }

  // The interesting part of every process travels along the molecule, so the
  // camera tracks it rather than sitting on the midpoint watching it leave.
  const follow = (x, y, z, rate = 1.6) => {
    _follow.set(x, y, z);
    orbit.goalTarget.lerp(_follow, Math.min(1, dt * rate));
  };
  if (currentId === 'translation' && stage.machines.ribosome.visible) {
    const p = stage.machines.ribosome.position;
    follow(p.x, p.y, p.z);
  } else if (currentId === 'replication') {
    follow(0, 0, zOf(Math.max(0, current.fork), stage.n), 1.1);
  } else if (currentId === 'transcription') {
    follow(0, -4, zOf(Math.max(0, Math.min(stage.n - 1, current.pos)), stage.n), 1.1);
  } else if (currentId === 'protein' && stage.peptide.count) {
    // Track the protein and frame to its actual size, so the shot tightens as
    // the chain collapses instead of leaving it a speck in the middle.
    const p = stage.peptide;
    _follow.set(0, 0, 0);
    for (let i = 0; i < p.count; i++) _follow.add(p.pos[i]);
    _follow.multiplyScalar(1 / p.count);
    orbit.goalTarget.lerp(_follow, Math.min(1, dt * 1.8));
    if (!orbit._pointers.size) {
      const want = Math.max(78, p.radiusOfGyration() * 7.5);
      orbit.goalDistance += (want - orbit.goalDistance) * Math.min(1, dt * 0.6);
    }
  }

  stage.peptide.render();

  const r = stage.renderer;
  r.begin();
  current.draw(r);
  if (stage.pool.chain.length && currentId !== 'translation') r.addChain(stage.pool.chain);
  r.end();

  const state = current.state();
  ui.update(state, currentId, sequence);
}

function frame() {
  requestAnimationFrame(frame);
  const dt = Math.min(clock.getDelta(), 0.05);
  elapsed += dt;
  step(dt);
  orbit.update(dt);
  composer.render(scene, camera, elapsed);
}

/**
 * Never hand a zero to the renderer. A background or not-yet-composited tab
 * reports innerWidth 0, and a zero-sized render target is a WebGL error rather
 * than a no-op.
 */
function resize() {
  const w = Math.max(1, innerWidth);
  const h = Math.max(1, innerHeight);
  renderer.setSize(w, h, false);
  camera.aspect = w / h;
  camera.updateProjectionMatrix();
  composer.setSize(w, h);
}
addEventListener('resize', resize);
resize();

/* ------------------------------------------------------------------ start */

stage.setSequence(sequence);
ui.setSequence(sequence, { isNew: true });
document.getElementById('preset').value = PRESETS[0].id;
document.getElementById('seq-note').textContent = PRESETS[0].note;
setScene('structure');
frameScene('structure', true);
ui.setPlaying(playing);
ui.ready();

// Boot the loop off a timeout rather than rAF: a hidden tab never fires rAF,
// and stills still need to be pullable from it.
setTimeout(() => { frame(); }, 0);

/** Manual driving, for tuning a shot or grabbing a frame without a loop. */
window.__helix = {
  stage, scenes, orbit, camera, renderer, composer,
  setScene,
  advance(seconds, fps = 60) {
    const dt = 1 / fps;
    for (let t = 0; t < seconds; t += dt) { step(dt); orbit.update(dt); }
    elapsed += seconds;
    composer.render(scene, camera, elapsed);
  },
  render() { composer.render(scene, camera, elapsed); },
};
