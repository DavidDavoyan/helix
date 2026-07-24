/**
 * ui.js — the panel, the readout and the sequence strip.
 *
 * The strip along the bottom is the part that matters. It is the same molecule
 * as the one in the viewport, written out as letters, with the machinery's
 * position marked on it — so when the ribosome steps three nucleotides you can
 * see which three. Clicking a base mutates it, and the readout says what that
 * did to the protein.
 */

import {
  PRESETS, BASES, DNA_LETTERS, AMINO, CODON_TABLE, cleanSequence,
  gcContent, meltingPoint, classifyMutation, transcribe,
} from './bio.js';

const $ = (id) => document.getElementById(id);

export class UI {
  constructor(handlers) {
    this.h = handlers;
    this.mutations = new Set();
    this.baseline = '';
    this._trackLetters = '';
    this._cells = [];
    this._frameCells = [];

    this.el = {
      ui: $('ui'),
      loading: $('loading'),
      modes: $('modes'),
      preset: $('preset'),
      seq: $('seq-input'),
      seqNote: $('seq-note'),
      seqLen: $('seq-len'),
      seqGc: $('seq-gc'),
      seqTm: $('seq-tm'),
      play: $('play'),
      restart: $('restart'),
      speed: $('speed'),
      speedVal: $('speed-val'),
      temp: $('temp'),
      tempVal: $('temp-val'),
      tempRow: $('temp-row'),
      toggles: $('toggles'),
      caption: $('caption'),
      progress: $('progress-bar'),
      stats: $('stats'),
      events: $('events'),
      peptideStrip: $('peptide-strip'),
      peptideList: $('peptide-list'),
      trackLabel: $('track-label'),
      trackStrip: $('track-strip'),
      trackFrame: $('track-frame'),
      mutationNote: $('mutation-note'),
      about: $('about'),
    };

    this._buildModes();
    this._buildPresets();
    this._buildToggles();
    this._wire();
  }

  /* ---------------------------------------------------------- construction */

  _buildModes() {
    const modes = [
      ['structure', 'The molecule'],
      ['replication', 'Replication'],
      ['transcription', 'Transcription'],
      ['translation', 'Translation'],
    ];
    this._modeButtons = {};
    for (const [id, label] of modes) {
      const b = document.createElement('button');
      b.type = 'button';
      b.textContent = label;
      b.setAttribute('role', 'tab');
      b.setAttribute('aria-selected', 'false');
      b.addEventListener('click', () => this.h.onMode(id));
      this.el.modes.appendChild(b);
      this._modeButtons[id] = b;
    }
  }

  _buildPresets() {
    for (const p of PRESETS) {
      const o = document.createElement('option');
      o.value = p.id;
      o.textContent = p.name;
      this.el.preset.appendChild(o);
    }
    const custom = document.createElement('option');
    custom.value = 'custom';
    custom.textContent = 'Custom';
    custom.hidden = true;
    this.el.preset.appendChild(custom);
  }

  _buildToggles() {
    const items = [
      ['bases', 'bases', true],
      ['backbone', 'backbone', true],
      ['sugars', 'sugars', true],
      ['hbonds', 'H-bonds', true],
      ['grooves', 'grooves', false],
      ['pool', 'free nucleotides', true],
      ['labels', 'labels', true],
      ['spin', 'auto-rotate', false],
    ];
    this._toggleEls = {};
    for (const [key, label, on] of items) {
      const l = document.createElement('label');
      l.className = on ? 'on' : '';
      const input = document.createElement('input');
      input.type = 'checkbox';
      input.checked = on;
      const span = document.createElement('span');
      span.textContent = label;
      l.append(input, span);
      input.addEventListener('change', () => {
        l.classList.toggle('on', input.checked);
        this.h.onToggle(key, input.checked);
      });
      this.el.toggles.appendChild(l);
      this._toggleEls[key] = { label: l, input };
    }
  }

  _wire() {
    const e = this.el;

    e.preset.addEventListener('change', () => {
      const p = PRESETS.find((x) => x.id === e.preset.value);
      if (p) {
        e.seq.value = p.seq;
        e.seqNote.textContent = p.note;
        this.h.onSequence(p.seq, true);
      }
    });

    let debounce;
    e.seq.addEventListener('input', () => {
      clearTimeout(debounce);
      debounce = setTimeout(() => {
        const clean = cleanSequence(e.seq.value);
        if (clean !== e.seq.value) {
          const at = e.seq.selectionStart;
          e.seq.value = clean;
          e.seq.setSelectionRange(Math.min(at, clean.length), Math.min(at, clean.length));
        }
        e.preset.value = 'custom';
        e.seqNote.textContent = '';
        this.h.onSequence(clean, true);
      }, 260);
    });

    e.play.addEventListener('click', () => this.h.onPlayToggle());
    e.restart.addEventListener('click', () => this.h.onRestart());

    e.speed.addEventListener('input', () => {
      const v = parseFloat(e.speed.value);
      e.speedVal.textContent = `${v.toFixed(2).replace(/0$/, '')}×`;
      this.h.onSpeed(v);
    });

    e.temp.addEventListener('input', () => {
      const v = parseFloat(e.temp.value);
      // 25 °C at rest up to 100 °C, where everything is denatured.
      e.tempVal.textContent = `${Math.round(25 + v * 75)} °C`;
      this.h.onTemperature(v);
    });

    $('about-open').addEventListener('click', () => { e.about.hidden = false; });
    $('about-close').addEventListener('click', () => { e.about.hidden = true; });
    e.about.addEventListener('click', (ev) => { if (ev.target === e.about) e.about.hidden = true; });

    window.addEventListener('keydown', (ev) => {
      if (ev.target.tagName === 'TEXTAREA' || ev.target.tagName === 'INPUT') return;
      if (ev.key === ' ') { ev.preventDefault(); this.h.onPlayToggle(); }
      if (ev.key === 'Escape') e.about.hidden = true;
      if (ev.key === 'r' || ev.key === 'R') this.h.onRestart();
    });
  }

  /* --------------------------------------------------------------- state */

  ready() {
    this.el.loading.classList.add('gone');
    this.el.ui.hidden = false;
    setTimeout(() => { this.el.loading.style.display = 'none'; }, 700);
  }

  setMode(id) {
    for (const [key, b] of Object.entries(this._modeButtons)) {
      b.setAttribute('aria-selected', String(key === id));
    }
    this.el.tempRow.style.display = id === 'structure' ? '' : 'none';
    this.el.peptideStrip.hidden = id !== 'translation';
    this._toggleEls.grooves.label.style.display = id === 'structure' ? '' : 'none';
    this.el.trackLabel.textContent = id === 'translation'
      ? 'mRNA, 5′ → 3′'
      : 'Coding strand, 5′ → 3′';
    this._trackLetters = '';   // force a rebuild
  }

  setPlaying(on) {
    this.el.play.textContent = on ? 'Pause' : 'Play';
  }

  setSequence(seq, { isNew = false } = {}) {
    if (this.el.seq.value !== seq) this.el.seq.value = seq;
    if (isNew) {
      this.baseline = seq;
      this.mutations.clear();
      this.el.mutationNote.hidden = true;
    }
    this.el.seqLen.textContent = `${seq.length} bp`;
    this.el.seqGc.textContent = `GC ${(gcContent(seq) * 100).toFixed(0)}%`;
    this.el.seqTm.textContent = `Tm ${meltingPoint(seq).toFixed(0)}°C`;
    this._trackLetters = '';
  }

  /** Repaint the readout. Called every frame; keeps DOM writes minimal. */
  update(state, sceneId, seq) {
    const e = this.el;

    if (this._lastCaption !== state.caption) {
      e.caption.textContent = state.caption;
      this._lastCaption = state.caption;
    }
    e.progress.style.width = `${Math.max(0, Math.min(1, state.progress || 0)) * 100}%`;

    const statsKey = JSON.stringify(state.stats);
    if (this._lastStats !== statsKey) {
      e.stats.replaceChildren();
      for (const [k, v] of state.stats || []) {
        const dt = document.createElement('dt');
        dt.textContent = k;
        const dd = document.createElement('dd');
        dd.textContent = v;
        e.stats.append(dt, dd);
      }
      this._lastStats = statsKey;
    }

    const evKey = (state.events || []).join('|');
    if (this._lastEvents !== evKey) {
      e.events.replaceChildren();
      for (const t of state.events || []) {
        const li = document.createElement('li');
        li.textContent = t;
        e.events.appendChild(li);
      }
      this._lastEvents = evKey;
    }

    if (sceneId === 'translation') {
      const key = (state.peptide || []).join('');
      if (this._lastPeptide !== key) {
        e.peptideList.replaceChildren();
        for (const aa of state.peptide || []) {
          const s = document.createElement('span');
          const info = AMINO[aa] || AMINO.G;
          s.textContent = info.tla;
          s.style.color = `#${info.colour.toString(16).padStart(6, '0')}`;
          e.peptideList.appendChild(s);
        }
        this._lastPeptide = key;
      }
    }

    this._track(state, sceneId, seq);
  }

  /* --------------------------------------------------------------- track */

  _track(state, sceneId, seq) {
    const letters = sceneId === 'translation' ? transcribe(seq) : seq;

    if (letters !== this._trackLetters) {
      this._trackLetters = letters;
      this.el.trackStrip.replaceChildren();
      this._cells = [];
      for (let i = 0; i < letters.length; i++) {
        const b = document.createElement('b');
        b.textContent = letters[i];
        b.className = letters[i];
        b.title = `${i + 1} · ${(BASES[letters[i]] || {}).name || ''}`;
        b.addEventListener('click', () => this._mutate(i));
        this.el.trackStrip.appendChild(b);
        this._cells.push(b);
      }
    }

    const t = state.track || {};
    const marks = new Set(t.mark || []);
    const [d0, d1] = t.doneRange || [-1, -2];

    for (let i = 0; i < this._cells.length; i++) {
      const c = this._cells[i];
      const active = marks.has(i);
      const done = i >= d0 && i <= d1;
      if (c._a !== active) { c.classList.toggle('active', active); c._a = active; }
      if (c._d !== done) { c.classList.toggle('done', done); c._d = done; }
      const mutated = this.mutations.has(i);
      if (c._m !== mutated) { c.classList.toggle('mutated', mutated); c._m = mutated; }
    }

    // Reading frame, only where there is one.
    const frame = t.frame;
    const frameKey = frame ? `${frame.start}:${frame.codon}:${letters.length}` : '';
    if (this._lastFrame !== frameKey) {
      this.el.trackFrame.replaceChildren();
      if (frame) {
        // Pad so codon boxes line up under their three bases: 3 * (17 + 1) - 1.
        const pad = document.createElement('i');
        pad.style.width = `${frame.start * 18}px`;
        pad.style.background = 'transparent';
        this.el.trackFrame.appendChild(pad);
        for (let c = 0; ; c++) {
          const at = frame.start + c * 3;
          if (at + 3 > letters.length) break;
          const codon = letters.slice(at, at + 3);
          const aa = CODON_TABLE[codon];
          const i = document.createElement('i');
          i.textContent = aa === '*' ? 'stop' : (AMINO[aa] ? AMINO[aa].tla : '?');
          if (aa === '*') { i.className = 'stop'; this.el.trackFrame.appendChild(i); break; }
          this.el.trackFrame.appendChild(i);
        }
      }
      this._lastFrame = frameKey;
      this._frameCells = [...this.el.trackFrame.querySelectorAll('i')];
    }
    if (frame) {
      this._frameCells.forEach((el, k) => {
        el.classList.toggle('current', k - 1 === frame.codon);
      });
    }
  }

  _mutate(i) {
    const before = this.el.seq.value;
    const isRna = this.el.trackLabel.textContent.startsWith('mRNA');
    const dnaIndex = i;
    if (dnaIndex >= before.length) return;

    const cur = before[dnaIndex];
    const order = DNA_LETTERS;
    const next = order[(order.indexOf(cur) + 1) % order.length];
    const after = before.slice(0, dnaIndex) + next + before.slice(dnaIndex + 1);

    if (this.baseline[dnaIndex] === next) this.mutations.delete(dnaIndex);
    else this.mutations.add(dnaIndex);

    this.el.seq.value = after;
    this.el.preset.value = 'custom';
    this.setSequence(after);
    this.h.onSequence(after, false);
    this._showMutation(after);
    void isRna;
  }

  _showMutation(after) {
    const e = this.el;
    if (!this.mutations.size || !this.baseline) { e.mutationNote.hidden = true; return; }
    const verdict = classifyMutation(this.baseline, after);
    const n = this.mutations.size;
    e.mutationNote.hidden = false;
    e.mutationNote.textContent = `${n} substitution${n > 1 ? 's' : ''} — ${verdict.label}`;
    e.mutationNote.style.color = {
      silent: '#7fe3a0', missense: '#ffcf8a', nonsense: '#ff8a94',
    }[verdict.type] || '#ffcf8a';
  }
}
