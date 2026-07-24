/**
 * post.js — a small composer: bright pass, two blur scales, then bloom,
 * tone-map, grade, vignette and grain in one composite.
 *
 * Hand-rolled rather than pulled from three's addons, so the whole thing stays
 * a static folder with one dependency.
 */

import * as THREE from './lib/three.module.js';

const VERT = `
  varying vec2 vUv;
  void main() {
    vUv = uv;
    gl_Position = vec4(position.xy, 0.0, 1.0);
  }`;

function fsQuad(material) {
  const geo = new THREE.BufferGeometry();
  geo.setAttribute('position', new THREE.BufferAttribute(new Float32Array([
    -1, -1, 0, 3, -1, 0, -1, 3, 0,
  ]), 3));
  geo.setAttribute('uv', new THREE.BufferAttribute(new Float32Array([0, 0, 2, 0, 0, 2]), 2));
  return new THREE.Mesh(geo, material);
}

export class Composer {
  constructor(renderer, { quality = 'high' } = {}) {
    this.renderer = renderer;
    this.quality = quality;
    this.scale = quality === 'low' ? 3 : 2;

    const size = renderer.getSize(new THREE.Vector2());
    const type = THREE.HalfFloatType;

    this.scene = new THREE.WebGLRenderTarget(size.x, size.y, {
      type, samples: quality === 'low' ? 0 : 4, depthBuffer: true,
    });
    this.scene.texture.colorSpace = THREE.LinearSRGBColorSpace;

    const half = (n) => Math.max(1, Math.floor(n / this.scale));
    const rt = (w, h) => {
      const t = new THREE.WebGLRenderTarget(w, h, { type, depthBuffer: false });
      t.texture.minFilter = THREE.LinearFilter;
      t.texture.magFilter = THREE.LinearFilter;
      return t;
    };
    this.brightRT = rt(half(size.x), half(size.y));
    this.blurA = rt(half(size.x), half(size.y));
    this.blurB = rt(half(size.x), half(size.y));
    this.wideA = rt(half(size.x) >> 1 || 1, half(size.y) >> 1 || 1);
    this.wideB = rt(half(size.x) >> 1 || 1, half(size.y) >> 1 || 1);

    this.bright = fsQuad(new THREE.ShaderMaterial({
      uniforms: {
        tDiffuse: { value: null },
        uThreshold: { value: 0.78 },
        uKnee: { value: 0.42 },
      },
      vertexShader: VERT,
      fragmentShader: `
        uniform sampler2D tDiffuse;
        uniform float uThreshold, uKnee;
        varying vec2 vUv;
        void main() {
          vec3 c = texture2D(tDiffuse, vUv).rgb;
          float l = dot(c, vec3(0.2126, 0.7152, 0.0722));
          // Soft knee, so highlights ease into the bloom instead of popping.
          float s = clamp((l - uThreshold + uKnee) / (2.0 * uKnee), 0.0, 1.0);
          float w = max(l - uThreshold, s * s * uKnee) / max(l, 1e-4);
          gl_FragColor = vec4(c * w, 1.0);
        }`,
      depthTest: false, depthWrite: false,
    }));

    this.blur = fsQuad(new THREE.ShaderMaterial({
      uniforms: {
        tDiffuse: { value: null },
        uDir: { value: new THREE.Vector2(1, 0) },
        uTexel: { value: new THREE.Vector2() },
      },
      vertexShader: VERT,
      fragmentShader: `
        uniform sampler2D tDiffuse;
        uniform vec2 uDir, uTexel;
        varying vec2 vUv;
        void main() {
          // Nine-tap gaussian folded into five bilinear fetches.
          vec2 o = uDir * uTexel;
          vec3 c = texture2D(tDiffuse, vUv).rgb * 0.2270270270;
          c += texture2D(tDiffuse, vUv + o * 1.3846153846).rgb * 0.3162162162;
          c += texture2D(tDiffuse, vUv - o * 1.3846153846).rgb * 0.3162162162;
          c += texture2D(tDiffuse, vUv + o * 3.2307692308).rgb * 0.0702702703;
          c += texture2D(tDiffuse, vUv - o * 3.2307692308).rgb * 0.0702702703;
          gl_FragColor = vec4(c, 1.0);
        }`,
      depthTest: false, depthWrite: false,
    }));

    this.composite = fsQuad(new THREE.ShaderMaterial({
      uniforms: {
        tScene: { value: null },
        tBloom: { value: null },
        tWide: { value: null },
        uBloom: { value: 0.62 },
        uExposure: { value: 1.16 },
        uVignette: { value: 0.7 },
        uGrain: { value: 0.012 },
        uTime: { value: 0 },
      },
      vertexShader: VERT,
      fragmentShader: `
        uniform sampler2D tScene, tBloom, tWide;
        uniform float uBloom, uExposure, uVignette, uGrain, uTime;
        varying vec2 vUv;

        // ACES filmic, Narkowicz's fit.
        vec3 aces(vec3 x) {
          const float a = 2.51, b = 0.03, c = 2.43, d = 0.59, e = 0.14;
          return clamp((x * (a * x + b)) / (x * (c * x + d) + e), 0.0, 1.0);
        }

        float hash(vec2 p) {
          return fract(sin(dot(p, vec2(127.1, 311.7))) * 43758.5453);
        }

        void main() {
          vec3 c = texture2D(tScene, vUv).rgb;
          vec3 b = texture2D(tBloom, vUv).rgb + texture2D(tWide, vUv).rgb * 0.7;
          c += b * uBloom;

          c *= uExposure;
          c = aces(c);

          // A touch of cool in the shadows, warm in the highlights.
          float l = dot(c, vec3(0.2126, 0.7152, 0.0722));
          c = mix(c * vec3(0.94, 0.98, 1.10), c * vec3(1.05, 1.01, 0.95), smoothstep(0.25, 0.9, l));

          vec2 d = vUv - 0.5;
          c *= mix(1.0, 1.0 - dot(d, d) * 1.6, uVignette);

          c += (hash(vUv * 1024.0 + fract(uTime) * 91.0) - 0.5) * uGrain;

          c = pow(max(c, 0.0), vec3(1.0 / 2.2));
          gl_FragColor = vec4(c, 1.0);
        }`,
      depthTest: false, depthWrite: false,
    }));

    this.quadScene = new THREE.Scene();
    this.quadCamera = new THREE.OrthographicCamera(-1, 1, 1, -1, 0, 1);
    this.setSize(size.x, size.y);
  }

  setSize(w, h) {
    const half = (n) => Math.max(1, Math.floor(n / this.scale));
    this.scene.setSize(w, h);
    this.brightRT.setSize(half(w), half(h));
    this.blurA.setSize(half(w), half(h));
    this.blurB.setSize(half(w), half(h));
    this.wideA.setSize(Math.max(1, half(w) >> 1), Math.max(1, half(h) >> 1));
    this.wideB.setSize(Math.max(1, half(w) >> 1), Math.max(1, half(h) >> 1));
  }

  _pass(mesh, target) {
    this.quadScene.children.length = 0;
    this.quadScene.add(mesh);
    this.renderer.setRenderTarget(target);
    this.renderer.clear();
    this.renderer.render(this.quadScene, this.quadCamera);
  }

  render(scene, camera, time) {
    const r = this.renderer;

    r.setRenderTarget(this.scene);
    r.clear();
    r.render(scene, camera);

    this.bright.material.uniforms.tDiffuse.value = this.scene.texture;
    this._pass(this.bright, this.brightRT);

    const blurInto = (src, mid, dst, w, h) => {
      const u = this.blur.material.uniforms;
      u.tDiffuse.value = src.texture;
      u.uTexel.value.set(1 / w, 1 / h);
      u.uDir.value.set(1, 0);
      this._pass(this.blur, mid);
      u.tDiffuse.value = mid.texture;
      u.uDir.value.set(0, 1);
      this._pass(this.blur, dst);
    };

    const bw = this.blurA.width, bh = this.blurA.height;
    blurInto(this.brightRT, this.blurB, this.blurA, bw, bh);
    blurInto(this.blurA, this.wideB, this.wideA, this.wideA.width, this.wideA.height);

    const u = this.composite.material.uniforms;
    u.tScene.value = this.scene.texture;
    u.tBloom.value = this.blurA.texture;
    u.tWide.value = this.wideA.texture;
    u.uTime.value = time;

    r.setRenderTarget(null);
    this._pass(this.composite, null);
  }
}
