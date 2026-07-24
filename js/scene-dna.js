/**
 * scene-dna.js — the molecule itself, and copying it.
 */

import * as THREE from './lib/three.module.js';
import { BDNA, BASES, complement, gcContent, meltingPoint } from './bio.js';
import { hydrogenBonds, facePair } from './chain.js';
import { Ribbon } from './geometry.js';
import { zOf, smoothstep, swingOut } from './stage.js';

const _bonds = [];

/* ========================================================== STRUCTURE ===== */

/**
 * The double helix, at rest, with a temperature control.
 *
 * The melting is the point of this scene. Turn the heat up and the molecule
 * does not come apart evenly — it opens at the ends first, and then wherever
 * the sequence is richest in A and T, because an A·T pair is held by two
 * hydrogen bonds and a G·C pair by three. Run it on the TATA box preset and
 * you can watch a promoter melt while the GC-rich flanks stay shut, which is
 * exactly why promoters are AT-rich.
 */
export class StructureScene {
  constructor(stage) {
    this.stage = stage;
    this.id = 'structure';
    this.title = 'The molecule';
    this.temperature = 0;
    this.showGrooves = false;
    this.spin = true;
    this._t = 0;

    // Two thin helical guides tracing the floor of each groove.
    this.grooves = [];
    for (const g of [
      { phase: BDNA.MINOR_GROOVE_SPAN / 2, colour: 0x3fd0c8, name: 'minor' },
      { phase: BDNA.MINOR_GROOVE_SPAN + (360 - BDNA.MINOR_GROOVE_SPAN) / 2, colour: 0xff9a3c, name: 'major' },
    ]) {
      const r = new Ribbon(280, 0.85, 6);
      const m = new THREE.Mesh(r.geometry, new THREE.MeshStandardMaterial({
        color: g.colour, emissive: new THREE.Color(g.colour).multiplyScalar(0.35),
        roughness: 0.5, metalness: 0, transparent: true, opacity: 0.85,
      }));
      m.frustumCulled = false;
      m.visible = false;
      stage.root.add(m);
      this.grooves.push({ ...g, ribbon: r, mesh: m, pts: [] });
    }
  }

  enter() {
    const s = this.stage;
    s.hideAllMachines();
    s.pool.reset();
    s.pool.setVisible(false);
    s.peptide.clear();
    for (let i = 0; i < s.n; i++) {
      s.coding.present[i] = 1;
      s.template.present[i] = 1;
    }
    s.newTop.present.fill(0);
    s.newBottom.present.fill(0);
    s.rna.present.fill(0);
    this._t = 0;
  }

  exit() {
    for (const g of this.grooves) g.mesh.visible = false;
  }

  /**
   * How hard this base pair is to open, 0..1.
   *
   * Not just this pair — DNA melts in cooperative domains, so a lone A·T
   * inside a GC block stays shut. The window average is a cheap stand-in for
   * that cooperativity. The ends get a bonus because a terminal pair is only
   * stacked on one side, and fraying from the ends is what actually happens.
   */
  _stability(i) {
    const s = this.stage;
    const n = s.n;
    const W = 4;
    let sum = 0, count = 0;
    for (let k = i - W; k <= i + W; k++) {
      if (k < 0 || k >= n) continue;
      const b = BASES[s.coding.letters[k]];
      sum += b ? b.hbonds : 2.5;
      count++;
    }
    const mean = sum / (count || 1);                     // 2 (all A·T) .. 3 (all G·C)
    const fromEnd = Math.min(i, n - 1 - i);
    const endPenalty = Math.exp(-fromEnd / 2.5) * 0.42;  // terminal pairs fray first
    return Math.max(0, (mean - 2) - endPenalty);         // 0 = melts easily
  }

  update(dt) {
    const s = this.stage;
    this._t += dt;

    const openFn = this.temperature <= 0.001 ? null : (i) => {
      const stab = this._stability(i);
      return smoothstep(stab - 0.14, stab + 0.30, this.temperature);
    };
    s.layIntact(openFn);

    // Breathing: even at rest, base pairs open and shut on their own. It is a
    // few percent of the time and it is how anything ever gets in.
    if (openFn === null) {
      for (let i = 0; i < s.n; i++) {
        const b = 0.5 + 0.5 * Math.sin(this._t * 1.7 + i * 0.9);
        s.coding.highlight[i] = s.template.highlight[i] = b * b * b * 0.10;
      }
    }

    // A groove is a feature of the intact duplex — once the strands come apart
    // there is nothing left to have a groove — so the guides fade out with the
    // first of the melting rather than tracing a helix that is no longer there.
    const grooveFade = 1 - smoothstep(0.02, 0.16, this.temperature);
    if (this.showGrooves && grooveFade > 0.01) {
      const z0 = zOf(0, s.n);
      for (const g of this.grooves) {
        g.mesh.material.opacity = 0.85 * grooveFade;
        let k = 0;
        for (let i = 0; i < s.n; i++) {
          const phase = (g.phase + i * BDNA.TWIST) * Math.PI / 180;
          if (!g.pts[k]) g.pts[k] = new THREE.Vector3();
          g.pts[k].set(
            Math.cos(phase) * (BDNA.HELIX_RADIUS + 0.9),
            Math.sin(phase) * (BDNA.HELIX_RADIUS + 0.9),
            z0 + i * BDNA.RISE,
          );
          k++;
        }
        if (k >= 2) g.ribbon.update(g.pts, k);
        g.mesh.visible = k >= 2;
      }
    } else {
      for (const g of this.grooves) g.mesh.visible = false;
    }
  }

  draw(r) {
    const s = this.stage;
    r.addChain(s.coding, { backboneColour: 0x93a9c6 });
    r.addChain(s.template, { backboneColour: 0x6f8299 });

    const openFn = this.temperature <= 0.001 ? null : (i) => {
      const stab = this._stability(i);
      return smoothstep(stab - 0.14, stab + 0.30, this.temperature);
    };

    for (let i = 0; i < s.n; i++) {
      const open = openFn ? openFn(i) : 0;
      if (open > 0.42) continue;                       // bonds are broken
      const info = BASES[s.coding.letters[i]];
      if (!info) continue;
      hydrogenBonds(s.coding, i, s.template, i, info.hbonds, _bonds);
      const fade = 1 - open / 0.42;
      for (const b of _bonds) {
        r.addBond(b.from, b.to, {
          radius: 0.22 * fade,
          colour: info.hbonds === 3 ? 0xbfe0ff : 0xffd8c0,
        });
      }
    }
  }

  state() {
    const s = this.stage;
    const gc = gcContent(s.seq);
    const melted = this.temperature > 0.001
      ? Array.from({ length: s.n }, (_, i) => (smoothstep(this._stability(i) - 0.14, this._stability(i) + 0.30, this.temperature) > 0.5 ? 1 : 0)).reduce((a, b) => a + b, 0)
      : 0;
    return {
      caption: this.temperature < 0.02
        ? 'B-DNA at rest. Two chains running in opposite directions, bases stacked in the core, backbone out in the water.'
        : melted >= s.n
          ? 'Fully denatured. Both strands are free — this is the melting step of every PCR cycle.'
          : `Melting. ${melted} of ${s.n} pairs open — the ends first, then wherever A·T runs.`,
      progress: this.temperature,
      stats: [
        ['length', `${s.n} bp`],
        ['GC content', `${(gc * 100).toFixed(0)}%`],
        ['melting point', `${meltingPoint(s.seq).toFixed(0)} °C`],
        ['turns', `${(s.n / (360 / BDNA.TWIST)).toFixed(1)}`],
        ['open pairs', `${melted}`],
      ],
      track: {
        mark: this.temperature > 0.001
          ? Array.from({ length: s.n }, (_, i) => i).filter(
            (i) => smoothstep(this._stability(i) - 0.14, this._stability(i) + 0.30, this.temperature) > 0.5)
          : [],
      },
    };
  }
}

/* ======================================================== REPLICATION ===== */

const FRAGMENT = 9;      // bp per Okazaki fragment — real ones are far longer
const DAUGHTER_OFFSET = 27;
const FORK_WIDTH = 5;

/**
 * The replication fork.
 *
 * Everything awkward about this picture follows from one fact: polymerase can
 * only add a nucleotide to a 3' end. The two template strands run in opposite
 * directions, so at a fork that travels one way, only one of them can be
 * copied continuously. The other has to be copied backwards, in pieces, each
 * piece started fresh behind the fork and run away from it until it meets the
 * last one. Those pieces are Okazaki fragments, and the nicks between them are
 * what ligase seals.
 *
 * The scene lays the two daughter duplexes out as real B-DNA around axes that
 * separate as the fork passes, so behind the fork you are looking at two
 * finished double helices, not a diagram.
 */
export class ReplicationScene {
  constructor(stage) {
    this.stage = stage;
    this.id = 'replication';
    this.title = 'Replication';
    this.reset();
  }

  reset() {
    const s = this.stage;
    this.fork = s.n - 1;
    this.leadTip = s.n - 1;
    this.leadJob = null;
    this.lagJob = null;
    this.fragments = [];
    this.pending = [];
    this.lastFragStart = s.n;
    this.current = null;      // fragment being synthesised
    this.lagFill = -1;
    this.sealing = null;
    this.done = false;
    this.events = [];
    this._t = 0;
  }

  enter() {
    const s = this.stage;
    s.hideAllMachines();
    s.pool.reset();
    s.pool.setVisible(true);
    s.peptide.clear();
    s.rna.present.fill(0);

    // New strands: the leading strand copies the coding strand, the lagging
    // strand copies the template. Both start empty.
    for (let i = 0; i < s.n; i++) {
      s.coding.present[i] = 1;
      s.template.present[i] = 1;
      s.newTop.letters[i] = complement(s.coding.letters[i]);
      s.newBottom.letters[i] = complement(s.template.letters[i]);
      s.newTop.present[i] = 0;
      s.newBottom.present[i] = 0;
    }
    s.newTop.length = s.n;
    s.newBottom.length = s.n;
    this.reset();
  }

  exit() {}

  /** 0 ahead of the fork (still duplex), 1 behind it (two daughters). */
  _open(i) {
    return smoothstep(this.fork, this.fork + FORK_WIDTH, i);
  }

  _log(text) {
    this.events.push(text);
    if (this.events.length > 5) this.events.shift();
  }

  update(dt, speed) {
    const s = this.stage;
    this._t += dt;
    const n = s.n;

    /* --- the fork advances, helicase first ------------------------------ */
    if (!this.done) {
      // G·C costs more to open than A·T. The fork actually slows in GC-rich
      // stretches, so let it.
      const at = Math.max(0, Math.min(n - 1, Math.round(this.fork)));
      const b = BASES[s.coding.letters[at]];
      const resistance = b && b.hbonds === 3 ? 0.72 : 1.0;
      this.fork -= dt * speed * 3.2 * resistance;
      if (this.fork < -FORK_WIDTH) { this.fork = -FORK_WIDTH; }
    }

    /* --- leading strand: one polymerase, never stops --------------------- */
    if (!this.leadJob && this.leadTip >= 0 && this.leadTip > this.fork) {
      const i = this.leadTip;
      this.leadJob = s.pool.incorporate({
        letter: s.newTop.letters[i],
        site: () => ({
          pos: [s.newTop.pos[i * 3], s.newTop.pos[i * 3 + 1], s.newTop.pos[i * 3 + 2]],
          out: [s.newTop.out[i * 3], s.newTop.out[i * 3 + 1], s.newTop.out[i * 3 + 2]],
          up: [s.newTop.up[i * 3], s.newTop.up[i * 3 + 1], s.newTop.up[i * 3 + 2]],
        }),
        duration: 0.30 / Math.max(0.25, speed),
        wrongTries: Math.random() < 0.55 ? 1 : 0,
        onLock: () => {
          s.newTop.present[i] = 1;
          s.newTop.highlight[i] = 1;
          this.leadTip--;
        },
      });
    }
    if (this.leadJob && this.leadJob.done) this.leadJob = null;

    /* --- lagging strand: start a fragment whenever enough template shows -- */
    if (this.lastFragStart - this.fork >= FRAGMENT && this.fork > -FORK_WIDTH) {
      const start = Math.max(0, Math.round(this.fork));
      const end = this.lastFragStart - 1;
      if (end >= start) {
        this.pending.push({ start, end, fill: start });
        this.lastFragStart = start;
        this._log(`Primase lays a primer at ${start + 1}; fragment ${this.fragments.length + this.pending.length} begins`);
      }
    }
    if (!this.current && this.pending.length) {
      this.current = this.pending.shift();
      this.lagFill = this.current.start;
    }

    if (this.current && !this.lagJob) {
      const i = this.lagFill;
      if (i > this.current.end) {
        // Fragment complete. The nick behind it is now ligase's problem.
        this.fragments.push(this.current);
        this.sealing = { at: this.current.end, t: 0 };
        this._log(`Fragment sealed — ${this.current.end - this.current.start + 1} nt joined to the strand`);
        this.current = null;
      } else {
        this.lagJob = s.pool.incorporate({
          letter: s.newBottom.letters[i],
          site: () => ({
            pos: [s.newBottom.pos[i * 3], s.newBottom.pos[i * 3 + 1], s.newBottom.pos[i * 3 + 2]],
            out: [s.newBottom.out[i * 3], s.newBottom.out[i * 3 + 1], s.newBottom.out[i * 3 + 2]],
            up: [s.newBottom.up[i * 3], s.newBottom.up[i * 3 + 1], s.newBottom.up[i * 3 + 2]],
          }),
          duration: 0.26 / Math.max(0.25, speed),
          wrongTries: Math.random() < 0.4 ? 1 : 0,
          onLock: () => {
            s.newBottom.present[i] = 1;
            s.newBottom.highlight[i] = 1;
            this.lagFill++;
          },
        });
      }
    }
    if (this.lagJob && this.lagJob.done) this.lagJob = null;

    if (this.sealing) {
      this.sealing.t += dt * speed;
      if (this.sealing.t > 1.1) this.sealing = null;
    }

    if (this.fork <= -FORK_WIDTH && this.leadTip < 0 && !this.current && !this.pending.length) {
      this.done = true;
    }

    /* --- geometry -------------------------------------------------------- */
    const z0 = zOf(0, n);
    for (let i = 0; i < n; i++) {
      const o = this._open(i);
      const dy = DAUGHTER_OFFSET * o;

      if (o < 0.02) {
        // Ahead of the fork: untouched parental duplex.
        s.coding.placeHelical(i, { z0, phaseOffset: 0 });
        s.template.placeHelical(i, { z0, phaseOffset: BDNA.MINOR_GROOVE_SPAN });
      } else {
        // Behind it: two daughter helices, each keeping its parental strand's
        // phase and taking a new partner in the place the old one left.
        s.coding.placeHelical(i, { z0, ay: dy, phaseOffset: 0, radiusScale: 1 + (1 - o) * 0.1 });
        s.template.placeHelical(i, { z0, ay: -dy, phaseOffset: BDNA.MINOR_GROOVE_SPAN });
      }
      s.newTop.placeHelical(i, { z0, ay: dy, phaseOffset: BDNA.MINOR_GROOVE_SPAN });
      s.newBottom.placeHelical(i, { z0, ay: -dy, phaseOffset: 0 });

      // Aim each base at whatever it is paired with. Ahead of the fork that is
      // the other parental strand; behind it, each template faces its new
      // partner instead. Where nothing is paired yet, the base swings out.
      if (o < 0.5) {
        facePair(s.coding, i, s.template, i);
      } else {
        if (s.newTop.present[i]) facePair(s.coding, i, s.newTop, i);
        else swingOut(s.coding, i, (o - 0.5) * 2, { ay: dy });
        if (s.newBottom.present[i]) facePair(s.template, i, s.newBottom, i);
        else swingOut(s.template, i, (o - 0.5) * 2, { ay: -dy });
      }

      const decay = Math.exp(-dt * 2.4);
      s.newTop.highlight[i] *= decay;
      s.newBottom.highlight[i] *= decay;
    }

    /* --- machines -------------------------------------------------------- */
    const m = s.machines;
    const fz = z0 + this.fork * BDNA.RISE;

    m.helicase.visible = !this.done;
    m.helicase.place([0, 0, fz], [0, 0, -1]);
    m.helicase.rotation.z = this._t * 2.2;
    s.label('helicase', 'helicase', m.helicase, [0, 2, -26]);

    const leadZ = z0 + Math.max(0, this.leadTip) * BDNA.RISE;
    m.polyLeading.visible = this.leadTip >= 0;
    m.polyLeading.place([0, DAUGHTER_OFFSET * this._open(this.leadTip) + 3, leadZ], [0, 0, -1]);
    s.label('lead', 'polymerase — leading strand', m.polyLeading, [0, 30, 0]);

    const lagI = this.current ? this.lagFill : this.fork;
    const lagZ = z0 + Math.max(0, Math.min(n - 1, lagI)) * BDNA.RISE;
    m.polyLagging.visible = !!this.current;
    m.polyLagging.place([0, -DAUGHTER_OFFSET * this._open(lagI) - 3, lagZ], [0, 0, 1]);
    s.label('lag', 'polymerase — lagging strand', m.polyLagging, [0, -32, 0]);

    m.primase.visible = !!this.current && this.lagFill - this.current.start < 3;
    if (m.primase.visible) {
      m.primase.place([12, -DAUGHTER_OFFSET * this._open(this.current.start) - 8,
        z0 + this.current.start * BDNA.RISE], [0, 0, 1]);
    }
    s.label('primase', 'primase', m.primase, [0, -24, 0]);

    m.ligase.visible = !!this.sealing;
    if (this.sealing) {
      const zz = z0 + this.sealing.at * BDNA.RISE;
      m.ligase.place([-14, -DAUGHTER_OFFSET - 6, zz], [0, 0, 1]);
      m.ligase.scale.setScalar(1 + Math.sin(this.sealing.t * 9) * 0.12);
    }
    s.label('ligase', 'ligase — sealing the nick', m.ligase, [0, -22, 0]);

    // Single-strand binding proteins sit on whatever template is exposed and
    // not yet copied.
    const exposedTop = [];
    for (let i = Math.max(0, Math.floor(this.fork)); i < n; i++) {
      if (!s.newTop.present[i]) exposedTop.push(i);
      if (exposedTop.length >= 2) break;
    }
    const exposedBottom = [];
    for (let i = Math.max(0, Math.floor(this.fork)); i < n; i++) {
      if (!s.newBottom.present[i]) exposedBottom.push(i);
      if (exposedBottom.length >= 2) break;
    }
    const spots = [
      ...exposedTop.map((i) => [i, +1]),
      ...exposedBottom.map((i) => [i, -1]),
    ];
    m.ssb.forEach((p, k) => {
      const spot = spots[k];
      p.visible = !!spot && !this.done;
      if (spot) {
        const [i, side] = spot;
        p.place([9 * side, side * (DAUGHTER_OFFSET * this._open(i) + 9), z0 + i * BDNA.RISE], [0, 0, 1]);
      }
    });
    s.label('ssb', 'single-strand binding protein', m.ssb[0], [0, 26, 34]);
  }

  draw(r) {
    const s = this.stage;
    r.addChain(s.coding, { backboneColour: 0x93a9c6 });
    r.addChain(s.template, { backboneColour: 0x93a9c6 });
    r.addChain(s.newTop, { backboneColour: 0xffb454 });
    r.addChain(s.newBottom, { backboneColour: 0xff8a5c });

    for (let i = 0; i < s.n; i++) {
      const o = this._open(i);
      const info = BASES[s.coding.letters[i]];
      if (!info) continue;

      // Parental pairs, still closed ahead of the fork.
      if (o < 0.45) {
        hydrogenBonds(s.coding, i, s.template, i, info.hbonds, _bonds);
        const fade = 1 - o / 0.45;
        for (const b of _bonds) r.addBond(b.from, b.to, { radius: 0.22 * fade, colour: 0xdfe9ff });
      }
      // New pairs, formed behind it.
      if (s.newTop.present[i]) {
        hydrogenBonds(s.coding, i, s.newTop, i, info.hbonds, _bonds);
        for (const b of _bonds) r.addBond(b.from, b.to, { radius: 0.22, colour: 0xffe3b0 });
      }
      if (s.newBottom.present[i]) {
        hydrogenBonds(s.template, i, s.newBottom, i, info.hbonds, _bonds);
        for (const b of _bonds) r.addBond(b.from, b.to, { radius: 0.22, colour: 0xffd0b0 });
      }
    }
  }

  state() {
    const s = this.stage;
    const made = s.newTop.countPresent() + s.newBottom.countPresent();
    const total = s.n * 2;
    const frag = this.fragments.length + (this.current ? 1 : 0);
    return {
      caption: this.done
        ? 'Two molecules where there was one, each keeping one old strand. That is what semi-conservative means.'
        : this.current
          ? `Leading strand runs continuously into the fork. Lagging strand is being built backwards, away from it — fragment ${frag}.`
          : 'Helicase splits the pair. Both new strands are laid 5′→3′, which is why one of them has to work in reverse.',
      progress: made / total,
      stats: [
        ['fork at', `${Math.max(0, Math.round(this.fork) + 1)} / ${s.n}`],
        ['nucleotides added', `${made} / ${total}`],
        ['Okazaki fragments', `${frag}`],
        ['leading strand', `${s.newTop.countPresent()} nt, unbroken`],
        ['lagging strand', `${s.newBottom.countPresent()} nt, in pieces`],
      ],
      events: this.events,
      track: {
        mark: [Math.max(0, Math.round(this.fork))],
        doneRange: [Math.max(0, Math.round(this.fork)) + 1, s.n - 1],
      },
    };
  }
}
