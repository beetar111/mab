/* OTC — Optically-written Tensor Core engine.
 * One file, no dependencies. Runs in the browser (global `OTC`) and in Node (`require`).
 * Units: lengths in µm, wavelengths in µm unless named *_nm, phases in rad.
 *
 *  platform + process ─► cell design (mode solver + BPM) ─► device statistics (noisy BPM realizations)
 *  design file (blocks, weights) ─► SVD + Clements decomposition ─► realized matrix with errors ─► accuracy
 *  design file ─► layout (paths, weight patches, fiducials) ─► raster tiles ─► DMD frames / LCoS holograms ─► GDSII
 */
(function (root) {
"use strict";
const PI = Math.PI, TAU = 2 * PI;
const OTC = { version: "1.0.0" };

// ───────────────────────── numerics ─────────────────────────
function makeFFT(n) {
  const bits = Math.round(Math.log2(n)), rev = new Uint32Array(n);
  for (let i = 0; i < n; i++) { let r = 0; for (let b = 0; b < bits; b++) r |= ((i >> b) & 1) << (bits - 1 - b); rev[i] = r; }
  const cs = new Float64Array(n / 2), sn = new Float64Array(n / 2);
  for (let i = 0; i < n / 2; i++) { cs[i] = Math.cos(TAU * i / n); sn[i] = Math.sin(TAU * i / n); }
  return function (re, im, inv) {
    for (let i = 0; i < n; i++) { const j = rev[i]; if (j > i) { let t = re[i]; re[i] = re[j]; re[j] = t; t = im[i]; im[i] = im[j]; im[j] = t; } }
    for (let size = 2; size <= n; size <<= 1) {
      const half = size >> 1, step = n / size;
      for (let i = 0; i < n; i += size) for (let j = i, k = 0; j < i + half; j++, k += step) {
        const l = j + half, c = cs[k], s = inv ? sn[k] : -sn[k];
        const tr = re[l] * c - im[l] * s, ti = re[l] * s + im[l] * c;
        re[l] = re[j] - tr; im[l] = im[j] - ti; re[j] += tr; im[j] += ti;
      }
    }
    if (inv) for (let i = 0; i < n; i++) { re[i] /= n; im[i] /= n; }
  };
}
OTC.makeFFT = makeFFT;
function erf(x) { const s = x < 0 ? -1 : 1; x = Math.abs(x); const t = 1 / (1 + 0.3275911 * x); return s * (1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x)); }
function hash(i, j, seed) { let h = Math.imul(i, 374761393) ^ Math.imul(j, 668265263) ^ Math.imul(seed, 2246822519); h = Math.imul(h ^ (h >>> 13), 1274126177); h ^= h >>> 16; return (h >>> 0) / 4294967296; }
function gl(i, j, s) { const u = hash(i, j, s) + 1e-12, v = hash(i, j, s + 977); return Math.sqrt(-2 * Math.log(u)) * Math.cos(TAU * v); }
function rng(seed) { let a = seed >>> 0; return () => { a = (a + 0x6D2B79F5) | 0; let t = Math.imul(a ^ (a >>> 15), 1 | a); t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t; return ((t ^ (t >>> 14)) >>> 0) / 4294967296; }; }
function gaussR(r) { return Math.sqrt(-2 * Math.log(r() + 1e-12)) * Math.cos(TAU * r()); }
OTC.rng = rng;

// complex helpers ([re, im] pairs) and dense complex matrices (arrays of rows)
const C = {
  mul: (a, b) => [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]],
  add: (a, b) => [a[0] + b[0], a[1] + b[1]],
  sub: (a, b) => [a[0] - b[0], a[1] - b[1]],
  conj: a => [a[0], -a[1]],
  abs: a => Math.hypot(a[0], a[1]),
  abs2: a => a[0] * a[0] + a[1] * a[1],
  arg: a => Math.atan2(a[1], a[0]),
  exp: p => [Math.cos(p), Math.sin(p)],
  scale: (a, s) => [a[0] * s, a[1] * s],
  div: (a, b) => { const d = b[0] * b[0] + b[1] * b[1]; return [(a[0] * b[0] + a[1] * b[1]) / d, (a[1] * b[0] - a[0] * b[1]) / d]; },
};
OTC.C = C;
const M = {
  eye: n => Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => [i === j ? 1 : 0, 0])),
  fromReal: A => A.map(r => r.map(v => [v, 0])),
  mul: (A, B) => { const n = A.length, m = B[0].length, k = B.length, R = []; for (let i = 0; i < n; i++) { const row = []; for (let j = 0; j < m; j++) { let re = 0, im = 0; for (let t = 0; t < k; t++) { const a = A[i][t], b = B[t][j]; re += a[0] * b[0] - a[1] * b[1]; im += a[0] * b[1] + a[1] * b[0]; } row.push([re, im]); } R.push(row); } return R; },
  dagger: A => A[0].map((_, j) => A.map(r => C.conj(r[j]))),
  copy: A => A.map(r => r.map(v => [v[0], v[1]])),
};
OTC.M = M;
// apply a 2×2 block T on rows (a, a+1) from the left / on columns from the right, in place
function leftApply(U, a, T) { const n = U[0].length; for (let j = 0; j < n; j++) { const x = U[a][j], y = U[a + 1][j]; U[a][j] = C.add(C.mul(T[0][0], x), C.mul(T[0][1], y)); U[a + 1][j] = C.add(C.mul(T[1][0], x), C.mul(T[1][1], y)); } }
function rightApply(U, a, T) { const n = U.length; for (let i = 0; i < n; i++) { const x = U[i][a], y = U[i][a + 1]; U[i][a] = C.add(C.mul(x, T[0][0]), C.mul(y, T[1][0])); U[i][a + 1] = C.add(C.mul(x, T[0][1]), C.mul(y, T[1][1])); } }
const inv2 = T => [[C.conj(T[0][0]), C.conj(T[1][0])], [C.conj(T[0][1]), C.conj(T[1][1])]];   // unitary inverse

// ───────────────────────── MZI model ─────────────────────────
// Physical MZI: coupler(α1) → phase θ on upper arm → coupler(α2) → (preceded by) phase φ on upper input.
// T = g · B(α2) · diag(e^{iθ},1) · B(α1) · diag(e^{iφ},1),  B(α) = [[cos α, i sin α],[i sin α, cos α]].
// Ideal (α = π/4, g = 1): T = i e^{iθ/2} [[e^{iφ} sin(θ/2), cos(θ/2)], [e^{iφ} cos(θ/2), −sin(θ/2)]];  θ = 0 → full cross.
function mziT(theta, phi, a1 = PI / 4, a2 = PI / 4, g = 1) {
  const B = a => [[[Math.cos(a), 0], [0, Math.sin(a)]], [[0, Math.sin(a)], [Math.cos(a), 0]]];
  const b1 = B(a1), b2 = B(a2), et = C.exp(theta), ep = C.exp(phi);
  const X = [[C.mul(et, b1[0][0]), C.mul(et, b1[0][1])], [b1[1][0], b1[1][1]]];      // diag(e^{iθ},1)·B1
  const Y = [[C.add(C.mul(b2[0][0], X[0][0]), C.mul(b2[0][1], X[1][0])), C.add(C.mul(b2[0][0], X[0][1]), C.mul(b2[0][1], X[1][1]))],
             [C.add(C.mul(b2[1][0], X[0][0]), C.mul(b2[1][1], X[1][0])), C.add(C.mul(b2[1][0], X[0][1]), C.mul(b2[1][1], X[1][1]))]];
  return [[C.scale(C.mul(Y[0][0], ep), g), C.scale(Y[0][1], g)], [C.scale(C.mul(Y[1][0], ep), g), C.scale(Y[1][1], g)]];
}
OTC.mziT = mziT;
const wrap = p => ((p % TAU) + TAU) % TAU;

// ───────────────────────── Clements decomposition ─────────────────────────
// Returns {mzis:[{m, theta, phi, col}] in light-travel order, D:[complex] output phases}, U = D · Π T.
function clements(Uin) {
  const N = Uin.length, U = M.copy(Uin), right = [], left = [], eps = 1e-12;
  for (let i = 0; i < N - 1; i++) {
    if (i % 2 === 0) {
      for (let j = 0; j <= i; j++) {
        const r = N - 1 - j, c = i - j, u0 = U[r][c], u1 = U[r][c + 1];
        let th, ph;
        if (C.abs(u0) < eps) { th = PI; ph = 0; } else if (C.abs(u1) < eps) { th = 0; ph = 0; }
        else { th = 2 * Math.atan2(C.abs(u1), C.abs(u0)); ph = C.arg(u0) - C.arg(u1) - PI; }
        const T = mziT(th, ph); rightApply(U, c, inv2(T)); right.push({ m: c, theta: th, phi: ph });
      }
    } else {
      for (let j = 1; j <= i + 1; j++) {
        const x = N + j - i - 2, c = j - 1, u0 = U[x - 1][c], u1 = U[x][c];
        let th, ph;
        if (C.abs(u1) < eps) { th = PI; ph = 0; } else if (C.abs(u0) < eps) { th = 0; ph = 0; }
        else { th = 2 * Math.atan2(C.abs(u0), C.abs(u1)); ph = C.arg(u1) - C.arg(u0); }
        const T = mziT(th, ph); leftApply(U, x - 1, T); left.push({ m: x - 1, theta: th, phi: ph });
      }
    }
  }
  // move the left-hand inverses through the diagonal:  L⁻¹ · D = D' · T'
  let D = U.map((r, i) => r[i]);
  const moved = [];
  for (let k = left.length - 1; k >= 0; k--) {
    const L = left[k], a = L.m, Li = inv2(mziT(L.theta, L.phi));
    const A = [[C.mul(Li[0][0], D[a]), C.mul(Li[0][1], D[a + 1])], [C.mul(Li[1][0], D[a]), C.mul(Li[1][1], D[a + 1])]];
    const th = 2 * Math.atan2(C.abs(A[0][0]), C.abs(A[0][1])), s = Math.sin(th / 2), c = Math.cos(th / 2), K = C.mul([0, 1], C.exp(th / 2));
    let da, db, ph = 0;
    if (c < 1e-9) { da = C.div(A[0][0], K); db = C.scale(C.div(A[1][1], K), -1); }
    else if (s < 1e-9) { da = C.div(A[0][1], K); db = C.div(A[1][0], K); }
    else { da = C.scale(C.div(A[0][1], K), 1 / c); db = C.scale(C.div(A[1][1], K), -1 / s); ph = C.arg(C.div(A[0][0], C.scale(C.mul(da, K), s))); }
    D[a] = da; D[a + 1] = db; moved.push({ m: a, theta: th, phi: ph });
  }
  // light order: right list (in order), then moved list reversed... moved was built from innermost (k = last) outward:
  // U = D_fin · T'_1 · … · T'_k · R_m … R_1  → applied: R_1 … R_m, T'_k … T'_1
  const seq = right.slice().concat(moved.slice());       // moved[0] = T'_k (innermost, applied first after R)
  const next = new Array(N).fill(0);
  for (const z of seq) { const col = Math.max(next[z.m], next[z.m + 1]); z.col = col; next[z.m] = next[z.m + 1] = col + 1; z.theta = wrap(z.theta); z.phi = wrap(z.phi); }
  return { mzis: seq, D, depth: Math.max(...next) };
}
function meshMatrix(N, mzis, D, devFn) {      // forward model; devFn(z) → 2×2 realized block
  const U = M.eye(N);
  for (const z of mzis) leftApply(U, z.m, devFn ? devFn(z) : mziT(z.theta, z.phi));
  if (D) for (let i = 0; i < N; i++) for (let j = 0; j < N; j++) U[i][j] = C.mul(D[i], U[i][j]);
  return U;
}
OTC.clements = clements; OTC.meshMatrix = meshMatrix;

// ───────────────────────── real SVD (one-sided Jacobi) ─────────────────────────
function svd(W) {
  const n = W.length, A = W.map(r => r.slice()), V = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? 1 : 0)));
  for (let sweep = 0; sweep < 60; sweep++) {
    let off = 0;
    for (let p = 0; p < n - 1; p++) for (let q = p + 1; q < n; q++) {
      let a = 0, b = 0, g = 0; for (let i = 0; i < n; i++) { a += A[i][p] ** 2; b += A[i][q] ** 2; g += A[i][p] * A[i][q]; }
      if (Math.abs(g) <= 1e-15 * Math.sqrt(a * b) || g === 0) continue;
      off = Math.max(off, Math.abs(g) / Math.sqrt(a * b));
      const z = (b - a) / (2 * g), t = Math.sign(z || 1) / (Math.abs(z) + Math.sqrt(1 + z * z)), c = 1 / Math.sqrt(1 + t * t), s = c * t;
      for (let i = 0; i < n; i++) { const x = A[i][p], y = A[i][q]; A[i][p] = c * x - s * y; A[i][q] = s * x + c * y; const v = V[i][p], w = V[i][q]; V[i][p] = c * v - s * w; V[i][q] = s * v + c * w; }
    }
    if (off < 1e-13) break;
  }
  const S = [], U = Array.from({ length: n }, () => new Array(n).fill(0));
  for (let j = 0; j < n; j++) { let s = 0; for (let i = 0; i < n; i++) s += A[i][j] ** 2; s = Math.sqrt(s); S.push(s); for (let i = 0; i < n; i++) U[i][j] = s > 1e-12 ? A[i][j] / s : 0; }
  // complete U for zero singular values (Gram–Schmidt against the standard basis)
  for (let j = 0; j < n; j++) if (S[j] <= 1e-12) {
    for (let e = 0; e < n; e++) {
      const v = new Array(n).fill(0); v[e] = 1;
      for (let k = 0; k < n; k++) if (k !== j && (S[k] > 1e-12 || k < j)) { let d = 0; for (let i = 0; i < n; i++) d += v[i] * U[i][k]; for (let i = 0; i < n; i++) v[i] -= d * U[i][k]; }
      const nv = Math.hypot(...v); if (nv > 1e-6) { for (let i = 0; i < n; i++) U[i][j] = v[i] / nv; break; }
    }
  }
  return { U, S, V };   // W = U · diag(S) · Vᵀ
}
OTC.svd = svd;

// ───────────────────────── weights library ─────────────────────────
OTC.weights = function (spec, N) {
  if (Array.isArray(spec)) return spec;
  const [kind, arg] = String(spec).split(":");
  const W = Array.from({ length: N }, () => new Array(N).fill(0));
  if (kind === "dct") { for (let k = 0; k < N; k++) for (let n = 0; n < N; n++) W[k][n] = Math.sqrt((k ? 2 : 1) / N) * Math.cos(PI * (n + 0.5) * k / N); }
  else if (kind === "hadamard") { for (let i = 0; i < N; i++) for (let j = 0; j < N; j++) { let b = i & j, p = 0; while (b) { p ^= b & 1; b >>= 1; } W[i][j] = (p ? -1 : 1) / Math.sqrt(N); } }
  else if (kind === "identity") { for (let i = 0; i < N; i++) W[i][i] = 1; }
  else if (kind === "blur") { for (let i = 0; i < N; i++) for (let d = -1; d <= 1; d++) { const j = i + d; if (j >= 0 && j < N) W[i][j] = d ? 0.25 : 0.5; } }
  else { const r = rng(+(arg || 1) * 7919 + 13); for (let i = 0; i < N; i++) for (let j = 0; j < N; j++) W[i][j] = gaussR(r) / Math.sqrt(N); }
  return W;
};

// ───────────────────────── platforms & process ─────────────────────────
OTC.PLATFORMS = {
  pcm_sin: {
    name: "Sb₂Se₃ phase-change film on SiN slab", short: "PCM on SiN",
    n0: 1.80, dn: 0.03, lossDbCm: 1.0, wgLevel: 0.5, levels: 256, ld: 0.15, writeNm: 532,
    rewritable: true, volatile: false,
    stack: ["Si wafer", "4 µm thermal SiO₂", "300 nm LPCVD Si₃N₄ (unpatterned slab)", "20 nm Sb₂Se₃ (sputtered, amorphous as deposited)", "10 nm SiO₂ cap (ALD or sputter)"],
    write: "Crystallize (raise index) by heating above ≈200 °C with the projected light; graded dose gives partial crystallization = analog weight level.",
    erase: "Re-amorphize with a short (ns) melt-quench pulse pattern, or flood-erase the die.",
    notes: "Effective-index contrast of the 20 nm film on 300 nm SiN is an estimate; verify with a mode solver for your stack. Sb₂Se₃ data: Delaney et al., Adv. Funct. Mater. 2020.",
  },
  polymer: {
    name: "Holographic photopolymer film (write once)", short: "Photopolymer",
    n0: 1.50, dn: 0.02, lossDbCm: 0.5, wgLevel: 0.5, levels: 256, ld: 0.2, writeNm: 532,
    rewritable: false, volatile: false,
    stack: ["Glass or Si/SiO₂ substrate", "5 µm low-index cladding", "15 µm holographic photopolymer (laminated or spin-coated)", "Protective cover film"],
    write: "One exposure polymerizes the pattern; a UV flood then fixes it. Weights are permanent, so this makes an inference ASIC printed by light.",
    erase: "Not erasable. Print a new die for new weights.",
    notes: "Index modulation of commercial holographic photopolymers is up to ≈0.03.",
  },
  ln_fe: {
    name: "LiNbO₃:Fe photorefractive crystal", short: "LiNbO₃:Fe",
    n0: 2.14, dn: 5e-4, lossDbCm: 0.2, wgLevel: 0.5, levels: 64, ld: 0.3, writeNm: 532,
    rewritable: true, volatile: false,
    stack: ["Z-cut LiNbO₃:Fe wafer (polished)", "Optional SiO₂ buffer and electrodes for electro-optic inputs"],
    write: "Photorefractive index change; can be thermally fixed for months of retention.",
    erase: "Uniform illumination plus heating.",
    notes: "Very weak contrast, so devices are cm-scale. It is included to show the scaling limit.",
  },
  sbn: {
    name: "SBN:61 photorefractive crystal (biased)", short: "SBN:61",
    n0: 2.30, dn: 1e-3, lossDbCm: 0.5, wgLevel: 0.5, levels: 64, ld: 0.5, writeNm: 532,
    rewritable: true, volatile: true,
    stack: ["SBN:61 crystal with bias electrodes (1–3 kV/cm)"],
    write: "Photorefractive index change under bias field.",
    erase: "Uniform light; decays in the dark over hours.",
    notes: "Volatile and research-grade.",
  },
};
OTC.PROCESS_DEFAULT = {
  plate: "lcos", frames: 8, projNA: 0.4, drift: 0.01,
  lcos: { pixels: [4160, 2464], pitch_um: 3.74, levels: 256, name: "LCoS phase SLM (4K class)" },
  dmd: { pixels: [2560, 1600], pitch_um: 7.56, levels: 256, name: "DMD (binary mirrors, 8-bit gray by PWM)" },
  exposure_s: 0.5, settle_s: 0.3,
};
// speckle on the written index for each plate technology (relative rms on Δn)
OTC.speckle = proc => proc.plate === "dmd" ? 0.01 : 0.075 / Math.sqrt(Math.max(1, proc.frames));

// ───────────────────────── BPM cell engine (2-D scalar, paraxial, split-step) ─────────────────────────
const NX = 256, fftX = makeFFT(NX);
function makeProf(w, sig, dx) {
  const h = Math.min(dx / 8, 0.02), U = w / 2 + 6 * sig + 2 * dx, N = Math.ceil(2 * U / h) + 1, t = new Float64Array(N), s2 = Math.SQRT2 * Math.max(sig, 1e-3);
  for (let i = 0; i < N; i++) { const u = -U + i * h; t[i] = 0.5 * (erf((u + w / 2) / s2) - erf((u - w / 2) / s2)); }
  return u => (u <= -U || u >= U) ? 0 : t[Math.round((u + U) / h)];
}
function dnRow(g, centers, z, nz, out, aout) {
  const c = g.cfg;
  for (let j = 0; j < NX; j++) {
    let s = 0; for (const xc of centers) s += g.prof(g.xs[j] - xc);
    if (s > 1) s = 1;
    let d = c.dn * c.lvl * s;
    if (nz) d *= 1 + nz(g.xs[j], z);
    out[j] = d; if (aout) aout[j] = 0;
  }
}
function makeNoise(c, seed, dz, dx) {
  if (!seed) return null;
  const lf = Math.max(c.ellF, dx), lz = Math.max(c.ellF, dz), sf = c.sigF * Math.sqrt(c.ellF / lz), ss = c.sigS, ls = c.ellS;
  const grid = (l, lzz, sd) => {
    let rA = null, kA = null, rB = null, kB = null;
    const row = k => { if (k === kA) return rA; if (k === kB) return rB; const r = new Map(); rB = rA; kB = kA; rA = r; kA = k; return r; };
    const val = (r, k, i) => { let v = r.get(i); if (v === undefined) { v = gl(i, k, sd); r.set(i, v); } return v; };
    return (x, z) => { const fx = x / l, fz = z / lzz, ix = Math.floor(fx), iz = Math.floor(fz), tx = fx - ix, tz = fz - iz, r0 = row(iz), r1 = row(iz + 1);
      return (1 - tz) * ((1 - tx) * val(r0, iz, ix) + tx * val(r0, iz, ix + 1)) + tz * ((1 - tx) * val(r1, iz + 1, ix) + tx * val(r1, iz + 1, ix + 1)); };
  };
  const gf = grid(lf, lz, seed), gs = grid(ls, ls, seed + 7919);
  return (x, z) => 1.5 * (sf * gf(x, z) + ss * gs(x, z));
}
function makePath(secs) {
  const segs = []; let z = 0;
  for (const s of secs) { segs.push({ z0: z, z1: z + s.L, a: s.a, b: s.b }); z += s.L; }
  let idx = 0;
  return { L: z, sep(q) { if (q < segs[idx].z0) idx = 0; while (idx < segs.length - 1 && q > segs[idx].z1) idx++; const s = segs[idx]; if (s.a === s.b) return s.a; const t = Math.min(1, Math.max(0, (q - s.z0) / (s.z1 - s.z0))); return s.a + (s.b - s.a) * (1 - Math.cos(PI * t)) / 2; } };
}
function bpm(g, path, fields, z0, z1, nz, centersFn) {
  const c = g.cfg, steps = Math.max(1, Math.round((z1 - z0) / g.dz)), dz = (z1 - z0) / steps, k0 = TAU / g.lam;
  const pc = new Float64Array(NX), ps = new Float64Array(NX);
  for (let m = 0; m < NX; m++) { const mm = m < NX / 2 ? m : m - NX, kx = TAU * mm / g.X, ph = -kx * kx * dz / (2 * k0 * c.n0); pc[m] = Math.cos(ph); ps[m] = Math.sin(ph); }
  const aBg = c.lossDbCm / 4.343 / 1e4, dn = new Float64Array(NX), mc = new Float64Array(NX), ms = new Float64Array(NX), mk = new Float64Array(NX);
  for (let j = 0; j < NX; j++) mk[j] = Math.exp(-0.03 * dz * g.mask[j]);
  for (let s = 0; s < steps; s++) {
    const z = z0 + (s + 0.5) * dz, sep = path.sep(z);
    dnRow(g, centersFn ? centersFn(sep) : [sep / 2, -sep / 2], z, nz, dn);
    for (let j = 0; j < NX; j++) { const ph = k0 * dn[j] * dz / 2, a = Math.sqrt(Math.exp(-aBg * dz / 2) * mk[j]); mc[j] = a * Math.cos(ph); ms[j] = a * Math.sin(ph); }
    for (const f of fields) {
      const re = f.re, im = f.im;
      for (let j = 0; j < NX; j++) { const r = re[j], i = im[j]; re[j] = r * mc[j] - i * ms[j]; im[j] = r * ms[j] + i * mc[j]; }
      fftX(re, im, false);
      for (let m = 0; m < NX; m++) { const r = re[m], i = im[m]; re[m] = r * pc[m] - i * ps[m]; im[m] = r * ps[m] + i * pc[m]; }
      fftX(re, im, true);
      for (let j = 0; j < NX; j++) { const r = re[j], i = im[j]; re[j] = r * mc[j] - i * ms[j]; im[j] = r * ms[j] + i * mc[j]; }
    }
  }
}
const fieldFrom = f => ({ re: Float64Array.from(f), im: new Float64Array(NX) });
const overlap = (f, m, dx) => { let r = 0, i = 0; for (let j = 0; j < NX; j++) { r += f.re[j] * m[j]; i += f.im[j] * m[j]; } return [r * dx, i * dx]; };
function solveModes(dn, n0, dx, count, k0) {
  const n = dn.length, inv = 1 / (dx * dx), e2 = inv * inv, d = new Float64Array(n);
  let lo = Infinity, hi = -Infinity;
  for (let i = 0; i < n; i++) { const nn = n0 + dn[i]; d[i] = -2 * inv + k0 * k0 * nn * nn; lo = Math.min(lo, d[i]); hi = Math.max(hi, d[i]); }
  lo -= 2 * inv; hi += 2 * inv;
  const cnt = lam => { let c = 0, q = d[0] - lam; if (q < 0) c++; for (let i = 1; i < n; i++) { q = d[i] - lam - e2 / (q === 0 ? 1e-300 : q); if (q < 0) c++; } return c; };
  const out = [];
  for (let k = 0; k < count; k++) {
    let a = lo, b = hi; for (let it = 0; it < 90; it++) { const m = (a + b) / 2; if (cnt(m) > n - 1 - k) b = m; else a = m; }
    const lam = (a + b) / 2, mu = lam + 1e-9 * Math.abs(lam); let y = new Float64Array(n).fill(1);
    for (let it = 0; it < 4; it++) {
      const c = new Float64Array(n), r = new Float64Array(n); let bb = d[0] - mu; c[0] = inv / bb; r[0] = y[0] / bb;
      for (let i = 1; i < n; i++) { bb = d[i] - mu - inv * c[i - 1]; c[i] = inv / bb; r[i] = (y[i] - inv * r[i - 1]) / bb; }
      const x = new Float64Array(n); x[n - 1] = r[n - 1]; for (let i = n - 2; i >= 0; i--) x[i] = r[i] - c[i] * x[i + 1];
      let s = 0; for (let i = 0; i < n; i++) s += x[i] * x[i]; s = Math.sqrt(s * dx); for (let i = 0; i < n; i++) x[i] /= s; y = x;
    }
    let mi = 0; for (let i = 0; i < n; i++) if (Math.abs(y[i]) > Math.abs(y[mi])) mi = i; if (y[mi] < 0) for (let i = 0; i < n; i++) y[i] = -y[i];
    out.push({ beta: Math.sqrt(lam), f: y });
  }
  return out;
}
function refineMode(g, centers, f0) {
  const k0 = TAU / g.lam, n0 = g.cfg.n0, dn = new Float64Array(NX); dnRow(g, centers, 0, null, dn);
  const tau = g.dz, pk = new Float64Array(NX), pd = new Float64Array(NX); let dmax = 0; for (let j = 0; j < NX; j++) dmax = Math.max(dmax, dn[j]);
  for (let m = 0; m < NX; m++) { const mm = m < NX / 2 ? m : m - NX, kx = TAU * mm / g.X; pk[m] = Math.exp(-kx * kx * tau / (2 * k0 * n0)); }
  for (let j = 0; j < NX; j++) pd[j] = Math.exp(k0 * (dn[j] - dmax) * tau / 2);
  const re = Float64Array.from(f0), im = new Float64Array(NX);
  for (let it = 0; it < 700; it++) {
    for (let j = 0; j < NX; j++) { re[j] *= pd[j]; im[j] = 0; }
    fftX(re, im, false); for (let m = 0; m < NX; m++) { re[m] *= pk[m]; im[m] *= pk[m]; } fftX(re, im, true);
    for (let j = 0; j < NX; j++) re[j] *= pd[j];
    let s = 0; for (let j = 0; j < NX; j++) s += re[j] * re[j]; s = Math.sqrt(s * g.dx); for (let j = 0; j < NX; j++) re[j] /= s;
  }
  return re;
}
function bendLoss(g, R) {     // dB/cm, by BPM in the conformally mapped bend n(x)(1 + x/R)
  const c = g.cfg, k0 = TAU / g.lam, L = Math.min(4000, PI * R), steps = Math.max(60, Math.round(L / g.dz)), dz = L / steps;
  const pc = new Float64Array(NX), ps = new Float64Array(NX);
  for (let m = 0; m < NX; m++) { const mm = m < NX / 2 ? m : m - NX, kx = TAU * mm / g.X, ph = -kx * kx * dz / (2 * k0 * c.n0); pc[m] = Math.cos(ph); ps[m] = Math.sin(ph); }
  const dn = new Float64Array(NX); dnRow(g, [0], 0, null, dn);
  const m0 = refineMode(g, [0], solveModes(dn, c.n0, g.dx, 1, k0)[0].f), f = fieldFrom(m0);
  const mc = new Float64Array(NX), ms = new Float64Array(NX);
  for (let j = 0; j < NX; j++) { const nb = isFinite(R) ? (c.n0 + dn[j]) * (1 + g.xs[j] / R) - c.n0 : dn[j], a = Math.sqrt(Math.exp(-0.03 * dz * g.mask[j])); mc[j] = a * Math.cos(k0 * nb * dz / 2); ms[j] = a * Math.sin(k0 * nb * dz / 2); }
  const Ps = [];
  for (let s = 0; s < steps; s++) {
    for (let j = 0; j < NX; j++) { const r = f.re[j], i = f.im[j]; f.re[j] = r * mc[j] - i * ms[j]; f.im[j] = r * ms[j] + i * mc[j]; }
    fftX(f.re, f.im, false); for (let m = 0; m < NX; m++) { const r = f.re[m], i = f.im[m]; f.re[m] = r * pc[m] - i * ps[m]; f.im[m] = r * ps[m] + i * pc[m]; } fftX(f.re, f.im, true);
    for (let j = 0; j < NX; j++) { const r = f.re[j], i = f.im[j]; f.re[j] = r * mc[j] - i * ms[j]; f.im[j] = r * ms[j] + i * mc[j]; }
    let P = 0; for (let j = 0; j < NX; j++) if (Math.abs(g.xs[j]) < 2.5 * g.w) P += f.re[j] ** 2 + f.im[j] ** 2; Ps.push(P);
  }
  const sa = Math.floor(steps * 0.15), Pa = Ps[sa]; if (!(Pa > 1e-8)) return 1000;
  let sb = steps - 1; for (let s = sa + 1; s < steps; s++) if (Ps[s] < 0.01 * Pa) { sb = s; break; }
  return Math.min(1000, Math.max(0, Math.log(Pa / Ps[sb]) / ((sb - sa) * dz)) * 4.343e4);
}

// Design the MZI cell for a platform + process. Returns geometry and BPM handles.
OTC.designCell = function (platKey, proc, lam = 1.55, opt = {}) {
  const P = OTC.PLATFORMS[platKey], k0 = TAU / lam;
  const sPsf = 0.217 * (P.writeNm / 1000) / proc.projNA, blur = Math.sqrt(sPsf * sPsf + P.ld * P.ld);
  const cfg = { n0: P.n0, dn: P.dn, lvl: P.wgLevel, lossDbCm: P.lossDbCm, blur, sigF: OTC.speckle(proc) * (sPsf / blur), ellF: 2.355 * blur, sigS: proc.drift, ellS: 500 };
  const adn = Math.abs(P.dn) * P.wgLevel, NA = Math.sqrt(2 * P.n0 * adn), w = 0.8 * PI / (k0 * NA);
  const g = { platKey, P, cfg, lam, w, NA, adn, sPsf };
  g.gap = w; g.sc = 2 * w; g.D = 5 * w;                       // lane pitch = MZI arm separation
  g.X = g.D + 12 * w + 16 * blur; g.dx = g.X / NX;
  g.dz = Math.min(10, 0.03 / (k0 * Math.abs(P.dn)));
  g.xs = new Float64Array(NX); for (let j = 0; j < NX; j++) g.xs[j] = -g.X / 2 + (j + 0.5) * g.dx;
  g.mask = new Float64Array(NX); const xb = 0.8 * g.X / 2;
  for (let j = 0; j < NX; j++) { const ax = Math.abs(g.xs[j]); g.mask[j] = ax < xb ? 0 : ((ax - xb) / (g.X / 2 - xb)) ** 2; }
  g.prof = makeProf(w, blur, g.dx);
  const dn = new Float64Array(NX);
  dnRow(g, [g.D / 2], 0, null, dn);
  const mU = solveModes(dn, P.n0, g.dx, 1, k0)[0];
  g.neff = mU.beta / k0; g.modeU = refineMode(g, [g.D / 2], mU.f); g.modeL = Float64Array.from(g.modeU).reverse();
  dnRow(g, [g.sc / 2, -g.sc / 2], 0, null, dn);
  const sm = solveModes(dn, P.n0, g.dx, 2, k0); g.Lc = PI / (sm[0].beta - sm[1].beta);
  // bends: bend loss curve → minimum radius for 0.1 dB/cm; S-bends must respect both slope and radius
  const base = bendLoss(g, Infinity); g.bend = [];
  let R = 25 * g.X;
  for (let k = 0; k < 14; k++) { const t = Math.max(1e-6, bendLoss(g, R) - base); g.bend.push({ R, dBcm: t }); if (t < 0.01 && k > 2) break; R *= 1.5; }
  const use = g.bend.filter(p => p.dBcm > 0.02 && p.dBcm < 30);
  let A = 0, Bf = 0; if (use.length >= 2) { const n = use.length; let sx = 0, sy = 0, sxx = 0, sxy = 0; for (const p of use) { const y = Math.log(p.dBcm); sx += p.R; sy += y; sxx += p.R * p.R; sxy += p.R * y; } Bf = (n * sxy - sx * sy) / (n * sxx - sx * sx); A = (sy - Bf * sx) / n; }
  g.Rmin = use.length >= 2 ? (Math.log(0.1) - A) / Bf : g.bend[g.bend.length - 1].R;
  g.bendAt = Rq => use.length >= 2 ? Math.exp(A + Bf * Rq) : 0;
  const off = (g.D - g.sc) / 2, thmax = 0.15 * NA / P.n0;
  g.Ls = Math.max(PI * off / (2 * thmax), PI * Math.sqrt(g.Rmin * off / 2));
  g.Lin = Math.max(20, 10 * w);
  // calibrate 50:50 couplers (bends included, noise off)
  g.L50 = g.Lc / 2;
  const cpl = L => [{ L: g.Lin, a: g.D, b: g.D }, { L: g.Ls, a: g.D, b: g.sc }, { L, a: g.sc, b: g.sc }, { L: g.Ls, a: g.sc, b: g.D }, { L: g.Lin, a: g.D, b: g.D }];
  for (let it = 0; it < 3; it++) {
    const path = makePath(cpl(g.L50)), f = fieldFrom(g.modeU); bpm(g, path, [f], 0, path.L, null);
    const u = overlap(f, g.modeU, g.dx), l = overlap(f, g.modeL, g.dx), Pu = C.abs2(u), Pl = C.abs2(l);
    const thm = Math.asin(Math.sqrt(Pl / (Pu + Pl))), Cc = PI / (2 * g.Lc), thb = thm - Cc * g.L50;
    g.L50 = Math.max(0, (PI / 4 - thb) / Cc); g.split = Pl / (Pu + Pl);
  }
  // weight patches: 2π of phase from the headroom above the waveguide level
  g.dnW = Math.abs(P.dn) * (1 - P.wgLevel); g.Lp = lam / (g.dnW * 0.9);
  g.margin = Math.max(5, 4 * w);
  g.Lphi = g.Lp + 2 * g.margin; g.La = g.Lp + 2 * g.margin;
  g.secs = [{ L: g.Lphi, a: g.D, b: g.D }, { L: g.Ls, a: g.D, b: g.sc }, { L: g.L50, a: g.sc, b: g.sc }, { L: g.Ls, a: g.sc, b: g.D }, { L: g.La, a: g.D, b: g.D }, { L: g.Ls, a: g.D, b: g.sc }, { L: g.L50, a: g.sc, b: g.sc }, { L: g.Ls, a: g.sc, b: g.D }];
  g.Lcell = g.secs.reduce((a, s) => a + s.L, 0);
  g.zArm = g.Lphi + 2 * g.Ls + g.L50 + g.La;
  // process centering: written MZIs come out over-coupled (speckle, blur and weak arm-to-arm coupling).
  // Measure test MZIs written with the real process and correct the coupler length so the mean α is π/4.
  g.centering = [];
  const setL = L => { g.L50 = L; g.secs[2].L = L; g.secs[6].L = L; g.Lcell = g.secs.reduce((a, s) => a + s.L, 0); g.zArm = g.Lphi + 2 * g.Ls + g.L50 + g.La; };
  if (opt.center !== false) {
    let am = 0; const seeds = [11, 23, 37];
    for (const sd of seeds) { const r = OTC.realizeMZI(g, 9000 + sd * 13); am += (r.a1 + r.a2) / 2; }
    am /= seeds.length; g.centering.push({ L50: g.L50, alphaMean: am });
    setL(Math.max(0, g.L50 + (PI / 4 - am) / (PI / (2 * g.Lc))));
  }
  return g;
};
// One noisy realization of the MZI cell → extracted coupler angles, arm-phase offset, loss.
OTC.realizeMZI = function (g, seed) {
  const path = makePath(g.secs), nz = makeNoise(g.cfg, seed, g.dz, g.dx);
  const f1 = fieldFrom(g.modeU), f2 = fieldFrom(g.modeL);
  bpm(g, path, [f1, f2], 0, g.zArm, nz);
  const half = (f, up) => { const o = { re: new Float64Array(NX), im: new Float64Array(NX) }; for (let j = 0; j < NX; j++) if ((g.xs[j] > 0) === up) { o.re[j] = f.re[j]; o.im[j] = f.im[j]; } return o; };
  const F = [half(f1, true), half(f1, false), half(f2, true), half(f2, false)];
  bpm(g, path, F, g.zArm, path.L, nz);
  const ov = f => [overlap(f, g.modeU, g.dx), overlap(f, g.modeL, g.dx)];
  const a = [ov(F[0]), ov(F[2])], b = [ov(F[1]), ov(F[3])];          // a[p][o]: input p, output o, via upper arm
  const A = [[a[0][0], a[1][0]], [a[0][1], a[1][1]]], B = [[b[0][0], b[1][0]], [b[0][1], b[1][1]]];
  const o = C.abs2(A[0][0]) + C.abs2(A[0][1]) > C.abs2(A[1][0]) + C.abs2(A[1][1]) ? 0 : 1;
  const a1 = Math.acos(Math.min(1, Math.sqrt(C.abs2(A[o][0]) / (C.abs2(A[o][0]) + C.abs2(A[o][1])))));
  const p = C.abs2(A[0][0]) + C.abs2(A[1][0]) > C.abs2(A[0][1]) + C.abs2(A[1][1]) ? 0 : 1;
  const a2 = Math.acos(Math.min(1, Math.sqrt(C.abs2(A[0][p]) / (C.abs2(A[0][p]) + C.abs2(A[1][p])))));
  const dth = C.arg(B[1][0]) - C.arg(A[1][0]);
  const e = C.exp(dth), M00 = C.add(C.mul(e, A[0][0]), B[0][0]), M10 = C.add(C.mul(e, A[1][0]), B[1][0]);
  const tot = C.abs2(M00) + C.abs2(M10);
  return { a1, a2, dth: Math.atan2(Math.sin(dth), Math.cos(dth)), IL: -10 * Math.log10(tot), erCross: 10 * Math.log10(C.abs2(M10) / Math.max(1e-15, C.abs2(M00))) };
};
// Port-phase noise between devices from the written-index statistics (rad rms over one cell length).
OTC.portPhaseSigma = g => { const c = g.cfg, k0 = TAU / g.lam; return k0 * 0.9 * g.adn * Math.sqrt(c.sigF ** 2 * c.ellF * g.Lcell + c.sigS ** 2 * c.ellS * g.Lcell); };

// ───────────────────────── program + simulate a tensor core ─────────────────────────
// pipeline: array of block types from {fanin, meshV, sigma, meshU, mesh, phase, fanout}
OTC.programCore = function (core, N) {
  const W = OTC.weights(core.weights, N), types = core.pipeline;
  const hasV = types.includes("meshV"), hasU = types.includes("meshU"), hasS = types.includes("sigma");
  const out = { N, W, mode: "", blocks: {} };
  if (hasV && hasU && hasS) {
    const { U, S, V } = svd(W), smax = Math.max(...S) || 1;
    // degenerate singular values (e.g. unitary W): the SVD is not unique. Rotate each degenerate block by a random
    // orthogonal Q (U→UQ, V→VQ leaves W unchanged) so no MZI sits at an exact bar/cross state that imperfect
    // couplers cannot reach.
    const rr = rng(4242), used = new Array(N).fill(false);
    for (let i = 0; i < N; i++) {
      if (used[i]) continue; const cl = [i]; used[i] = true;
      for (let j = i + 1; j < N; j++) if (!used[j] && Math.abs(S[j] - S[i]) < 1e-9 * smax) { cl.push(j); used[j] = true; }
      if (cl.length < 2) continue;
      const m = cl.length, Q = Array.from({ length: m }, () => Array.from({ length: m }, () => gaussR(rr)));
      for (let c = 0; c < m; c++) { for (let d = 0; d < c; d++) { let dot = 0; for (let t = 0; t < m; t++) dot += Q[t][c] * Q[t][d]; for (let t = 0; t < m; t++) Q[t][c] -= dot * Q[t][d]; } let nn = 0; for (let t = 0; t < m; t++) nn += Q[t][c] ** 2; nn = Math.sqrt(nn); for (let t = 0; t < m; t++) Q[t][c] /= nn; }
      for (const X of [U, V]) for (let row = 0; row < N; row++) { const old = cl.map(c => X[row][c]); cl.forEach((c, ci) => { let v = 0; for (let t = 0; t < m; t++) v += old[t] * Q[t][ci]; X[row][c] = v; }); }
    }
    const Vt = V[0].map((_, j) => V.map(r => r[j]));
    out.scale = smax; out.S = S.map(s => s / smax);
    out.blocks.meshV = clements(M.fromReal(Vt));
    out.blocks.meshU = clements(M.fromReal(U));
    out.mode = "arbitrary (SVD: W = U·Σ·Vᵀ)";
    out.target = W.map(r => r.map(v => v / smax));
  } else if (hasU || types.includes("mesh")) {
    // closest unitary (polar decomposition) if W is not unitary
    const { U, V } = svd(W); const Q = U.map((r, i) => V.map((_, j) => r.reduce((a, _, k) => a + U[i][k] * V[j][k], 0)));
    out.blocks.meshU = clements(M.fromReal(Q)); out.scale = 1; out.S = new Array(N).fill(1);
    out.mode = "unitary only (no Σ block): the closest unitary to W is used"; out.target = Q;
  } else { out.mode = "diagonal only"; out.S = W.map((r, i) => Math.abs(r[i])); out.scale = Math.max(...out.S) || 1; out.S = out.S.map(s => s / out.scale); out.target = W.map((r, i) => r.map((v, j) => (i === j ? v / out.scale : 0))); }
  return out;
};
// Hardware-aware programming of one MZI (measured couplers a1, a2): find (θ', φ') so that the real block,
// fed with inputs that carry phase debts p, equals diag(e^{iq})·T_ideal. Returns new output debts q.
function fitMZI(Tid, p0, p1, a1, a2) {
  const tgt0 = C.abs(Tid[0][0]), tgt1 = C.abs(Tid[1][0]);
  const magErr = t => { const R = mziT(t, 0, a1, a2, 1); return Math.abs(C.abs(R[0][0]) - tgt0) + Math.abs(C.abs(R[1][0]) - tgt1); };
  const golden = (f, a, b) => { const gr = 0.6180339887; let c = b - gr * (b - a), d = a + gr * (b - a); for (let i = 0; i < 40; i++) { if (f(c) < f(d)) b = d; else a = c; c = b - gr * (b - a); d = a + gr * (b - a); } return (a + b) / 2; };
  let bt = 0, be = Infinity; for (let k = 0; k < 256; k++) { const t = TAU * k / 256, e = magErr(t); if (e < be) { be = e; bt = t; } }
  const th = golden(magErr, bt - TAU / 256, bt + TAU / 256);
  const ep = [C.exp(p0), C.exp(p1)];
  const fitCost = ph => {
    const R = mziT(th, ph, a1, a2, 1); let cost = 0; const q = [];
    for (let r = 0; r < 2; r++) {
      const Rr = [C.mul(R[r][0], ep[0]), C.mul(R[r][1], ep[1])];
      const qq = C.arg(C.add(C.mul(Rr[0], C.conj(Tid[r][0])), C.mul(Rr[1], C.conj(Tid[r][1])))); q.push(qq);
      const e = C.exp(qq); for (let c = 0; c < 2; c++) cost += C.abs2(C.sub(Rr[c], C.mul(e, Tid[r][c])));
    }
    return { cost, q };
  };
  let bp = 0, bc = Infinity; for (let k = 0; k < 128; k++) { const ph = TAU * k / 128, c = fitCost(ph).cost; if (c < bc) { bc = c; bp = ph; } }
  const ph = golden(p => fitCost(p).cost, bp - TAU / 128, bp + TAU / 128);
  const q = fitCost(ph).q;
  return { th: wrap(th), ph: wrap(ph), q0: q[0], q1: q[1] };
}
OTC.fitMZI = fitMZI;
// Build the realized N×N complex transfer of a programmed core, given device statistics.
//  uncalibrated: ideal phases written open-loop (quantized), devices carry their BPM-measured errors + port phase noise
//  calibrated:   every MZI measured; hardware-aware phases (fitMZI) written closed-loop; residual = level quantization
OTC.realizeCore = function (prog, stats, opt) {
  const N = prog.N, r = rng(opt.seed || 1), cal = opt.calibrated, levels = opt.levels || 256;
  const qStep = TAU / levels, sigRes = 0.003, sigPort = stats.portSigma;
  const q = p => Math.round(wrap(p) / qStep) * qStep;
  let k = 0; const debt = new Array(N).fill(0);
  const T = M.eye(N);
  const pick = () => stats.pool[(k++) % stats.pool.length];
  const write = (th, ph) => [q(th) + sigRes * gaussR(r), q(ph) + sigRes * gaussR(r)];
  const mzi = z => {
    const d = pick(), g = Math.pow(10, -d.IL / 20); let th, ph;
    if (cal) { const f = fitMZI(mziT(z.theta, z.phi), debt[z.m], debt[z.m + 1], d.a1, d.a2); [th, ph] = write(f.th, f.ph); debt[z.m] = f.q0; debt[z.m + 1] = f.q1; }
    else { [th, ph] = write(z.theta, z.phi); th += d.dth + sigPort * gaussR(r); ph += sigPort * gaussR(r); }
    leftApply(T, z.m, mziT(th, ph, d.a1, d.a2, g));
  };
  // pass-through lanes in each mesh column are waveguides of the same length: same material loss, own phase noise
  const gPass = Math.pow(10, -(stats.passIL || 0) / 20);
  const mesh = blk => {
    const byCol = Array.from({ length: blk.depth }, () => []); for (const z of blk.mzis) byCol[z.col].push(z);
    for (const col of byCol) {
      const used = new Set(); for (const z of col) { mzi(z); used.add(z.m); used.add(z.m + 1); }
      for (let i = 0; i < N; i++) if (!used.has(i)) { const f = C.scale(C.exp(cal ? 0 : sigPort * gaussR(r)), gPass); for (let j = 0; j < N; j++) T[i][j] = C.mul(f, T[i][j]); }
    }
  };
  if (prog.blocks.meshV) {
    mesh(prog.blocks.meshV);
    const Dv = prog.blocks.meshV.D;
    for (let i = 0; i < N; i++) {
      const d = pick(), g = Math.pow(10, -d.IL / 20), s = Math.min(1, prog.S[i]);
      let th, ph;
      if (cal) {
        let bt = 0, be = Infinity; for (let kk = 0; kk < 1024; kk++) { const t = TAU * kk / 1024, e = Math.abs(C.abs(mziT(t, 0, d.a1, d.a2, 1)[0][0]) - s); if (e < be) { be = e; bt = t; } }
        th = bt; const R00 = mziT(th, 0, d.a1, d.a2, 1)[0][0];
        ph = C.arg(Dv[i]) - C.arg(R00) - debt[i]; debt[i] = 0;          // phase of R00 ∝ e^{iφ}: cancel the debt here
        [th, ph] = write(th, ph);
      } else { th = 2 * Math.asin(s); ph = -PI / 2 - th / 2 + C.arg(Dv[i]); [th, ph] = write(th, ph); th += d.dth + sigPort * gaussR(r); ph += sigPort * gaussR(r); }
      const b = mziT(th, ph, d.a1, d.a2, g);
      for (let j = 0; j < N; j++) T[i][j] = C.mul(b[0][0], T[i][j]);
    }
  }
  if (prog.blocks.meshU) mesh(prog.blocks.meshU);
  // coherent detection: LO phase per output. Calibrated: fitted against the target (digital per-channel gain as well).
  const Du = prog.blocks.meshU ? prog.blocks.meshU.D : new Array(N).fill([1, 0]);
  const tgt = prog.target, lo = [], gain = [];
  for (let i = 0; i < N; i++) {
    let psi = -C.arg(Du[i]), gn = 1;      // output = D_U·(mesh); LO multiplies by e^{-iψ}
    if (cal) { let s = [0, 0], nn = 0; for (let j = 0; j < N; j++) { s = C.add(s, C.scale(T[i][j], tgt[i][j])); nn += tgt[i][j] ** 2; } psi = C.arg(s); const proj = C.abs(s); gn = proj > 1e-9 ? nn / proj : 1; }
    lo.push(psi); gain.push(gn);
  }
  const Teff = T.map((row, i) => row.map(v => C.scale(C.mul(v, C.exp(-lo[i])), gain[i])));
  return { T, Teff, lo, gain };
};
// Run vectors through a realized core with detector noise and converter quantization.
OTC.runCore = function (prog, real, opt) {
  const N = prog.N, r = rng((opt.seed || 1) + 99), S = opt.samples || 256, dac = opt.dacBits || 8, adc = opt.adcBits || 8, snr = opt.snrLin || 1e5;
  const qz = (v, b) => { const L = Math.pow(2, b - 1) - 1; return Math.max(-1, Math.min(1, Math.round(v * L) / L)); };
  // per-channel ADC full scale: 4σ of the expected output for uniform inputs, never more than the row's L1 norm
  const FS = prog.target.map(row => { const l2 = Math.sqrt(row.reduce((a, v) => a + v * v, 0)), l1 = row.reduce((a, v) => a + Math.abs(v), 0); return Math.max(1e-6, Math.min(l1, 4 * l2 / Math.sqrt(3))); });
  let e2 = 0, y2 = 0, ymax = 0; const errs = [], pairs = [];
  for (let s = 0; s < S; s++) {
    const x = Array.from({ length: N }, () => qz(2 * r() - 1, dac));
    for (let i = 0; i < N; i++) {
      let ideal = 0, c = [0, 0];
      for (let j = 0; j < N; j++) { ideal += prog.target[i][j] * x[j]; c = C.add(c, C.scale(real.T[i][j], x[j])); }
      let y = (c[0] * Math.cos(real.lo[i]) + c[1] * Math.sin(real.lo[i])) * real.gain[i];
      y = FS[i] * qz(y / FS[i] + gaussR(r) / Math.sqrt(snr), adc);   // receiver noise and ADC at full scale FS[i]
      const e = y - ideal; e2 += e * e; y2 += ideal * ideal; ymax = Math.max(ymax, Math.abs(ideal)); errs.push(e);
      if (pairs.length < 600) pairs.push([ideal, y]);
    }
  }
  const rmsE = Math.sqrt(e2 / (S * N)), rmsY = Math.sqrt(y2 / (S * N));
  const snrDb = 20 * Math.log10(rmsY / rmsE);
  return { rmsE, rmsY, snrDb, enob: (snrDb - 1.76) / 6.02, errs, pairs };
};
OTC.fidelity = (A, B) => {     // |Tr(A†B)|² / (Tr(A†A)·Tr(B†B))
  let t = [0, 0], na = 0, nb = 0; const n = A.length;
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) { t = C.add(t, C.mul(C.conj(A[i][j]), B[i][j])); na += C.abs2(A[i][j]); nb += C.abs2(B[i][j]); }
  return C.abs2(t) / (na * nb);
};

// ───────────────────────── system performance model ─────────────────────────
OTC.performance = function (design, g, layout, accuracy) {
  const s = design.system, N = design.N, K = design.cores.length, f = s.symbolRate_GHz * 1e9, lambdas = s.wavelengths || 1;
  const macs = K * N * N * f * lambdas, ops = 2 * macs;
  const nMziCore = design.cores.reduce((a, c) => a + (c.pipeline.includes("meshV") ? N * (N - 1) / 2 : 0) + (c.pipeline.includes("meshU") || c.pipeline.includes("mesh") ? N * (N - 1) / 2 : 0) + (c.pipeline.includes("sigma") ? N : 0), 0);
  const e = s.energy_pJ, chans = K * N * lambdas;
  const P = {
    laser: K * s.laserPerCore_mW * 1e-3 * lambdas / s.wallPlug,
    modulators: chans * e.modulator * 1e-12 * f,
    dacs: chans * e.dac * 1e-12 * f,
    receivers: chans * (e.tia + e.adc) * 1e-12 * f,
    weights: 0,
    thermal: s.tec_W,
  };
  const total = Object.values(P).reduce((a, b) => a + b, 0);
  // silicon-photonics reference: same architecture, thermo-optic phase shifters hold the weights
  const siWeights = nMziCore * 2 * s.si.Ppi_mW * 1e-3;     // mean phase π → mean power P_π per shifter, 2 shifters/MZI
  const siTotal = total - P.weights + siWeights;
  const siArea = nMziCore * s.si.mziArea_mm2;
  return {
    N, K, lambdas, macs, ops, P, total, topsW: ops / 1e12 / total, pJperMAC: total / macs * 1e12, nMzi: nMziCore,
    si: { weights: siWeights, total: siTotal, topsW: ops / 1e12 / siTotal, area_mm2: siArea },
    gpu: { name: s.gpuRef.name, topsW: s.gpuRef.tops / s.gpuRef.watts, tops: s.gpuRef.tops, watts: s.gpuRef.watts },
    die: layout ? { W: layout.die.W / 1000, H: layout.die.H / 1000, area: layout.die.W * layout.die.H / 1e6 } : null,
    latency_ns: layout ? layout.die.W * g.neff / 299.792458 / 1000 : NaN,
    enob: accuracy ? accuracy.enob : NaN,
  };
};
// SNR at the receiver for coherent (balanced homodyne, shot-noise-limited) detection
OTC.receiverSNR = function (design, insertionLossDb) {
  const s = design.system, N = design.N, B = s.symbolRate_GHz * 1e9 / 2;
  const Ps = s.laserPerCore_mW * 1e-3 * Math.pow(10, -insertionLossDb / 10) / N;   // full-scale power per output
  const R = 1.0, q = 1.602e-19, snr = 2 * R * Ps / (q * B);
  return { Ps, snr, snrDb: 10 * Math.log10(snr) };
};

// ───────────────────────── layout compiler ─────────────────────────
// Produces paths (center polylines, width, write level) and weight patches in die coordinates (µm).
// x = propagation, y = transverse (down). Cores are stacked in y.
OTC.compile = function (design, g, progs) {
  const N = design.N, P = g.D, w = g.w, wl = g.P.wgLevel, fiberPitch = design.system.fiberPitch_um || 127;
  const paths = [], patches = [], marks = [], blocks = [];
  const lv = theta => wl + (1 - wl) * (wrap(theta) / TAU);
  const rc = (x0, y0, L, dy, n = 24) => { const pts = []; for (let i = 0; i <= n; i++) { const t = i / n; pts.push([x0 + L * t, y0 + dy * (1 - Math.cos(PI * t)) / 2]); } return pts; };
  const sLen = dy => Math.max(PI * Math.sqrt(g.Rmin * Math.abs(dy) / 2), 1);   // raised-cosine S-bend whose tightest radius equals Rmin
  const coreH = Math.max(N * fiberPitch, (2 * N) * P) + 400;
  const facetConv = Math.max(100, 30 * w);
  let dieW = 0;
  design.cores.forEach((core, ci) => {
    const yc = ci * coreH + coreH / 2, lane = i => yc + (i - (N - 1) / 2) * P, fib = i => yc + (i - (N - 1) / 2) * fiberPitch;
    let x = 200;  // edge margin
    const prog = progs[ci];
    for (const b of core.pipeline) {
      const x0 = x;
      if (b === "fanin" || b === "fanout") {
        const dyMax = Math.abs(fib(0) - lane(0)), L = sLen(dyMax);
        for (let i = 0; i < N; i++) {
          const ya = b === "fanin" ? fib(i) : lane(i), yb = b === "fanin" ? lane(i) : fib(i);
          if (b === "fanin") { // facet mode converter: written level ramps up from the facet
            for (let s = 0; s < 4; s++) paths.push({ pts: [[x + s * facetConv / 4, ya], [x + (s + 1) * facetConv / 4, ya]], w, level: wl * (0.25 + 0.25 * s), kind: "converter" });
            paths.push({ pts: rc(x + facetConv, ya, L, yb - ya), w, level: wl, kind: "fan" });
          } else {
            paths.push({ pts: rc(x, ya, L, yb - ya), w, level: wl, kind: "fan" });
            for (let s = 0; s < 4; s++) paths.push({ pts: [[x + L + s * facetConv / 4, yb], [x + L + (s + 1) * facetConv / 4, yb]], w, level: wl * (1 - 0.25 * s), kind: "converter" });
          }
        }
        x += L + facetConv;
      } else if (b === "meshV" || b === "meshU" || b === "mesh") {
        const blk = prog.blocks[b === "meshV" ? "meshV" : "meshU"];
        const cols = blk.depth; const byCol = Array.from({ length: cols }, () => []);
        for (const z of blk.mzis) byCol[z.col].push(z);
        for (let c = 0; c < cols; c++) {
          const xc = x + c * g.Lcell, used = new Set();
          for (const z of byCol[c]) {
            used.add(z.m); used.add(z.m + 1);
            const yu = lane(z.m), yl = lane(z.m + 1), mid = (yu + yl) / 2, half = g.sc / 2;
            for (const [y0, sign] of [[yu, 1], [yl, -1]]) {
              const pts = [[xc, y0], [xc + g.Lphi, y0]];
              let xx = xc + g.Lphi;
              const yin = mid - sign * half;
              pts.push(...rc(xx, y0, g.Ls, yin - y0).slice(1)); xx += g.Ls;
              pts.push([xx + g.L50, yin]); xx += g.L50;
              pts.push(...rc(xx, yin, g.Ls, y0 - yin).slice(1)); xx += g.Ls;
              pts.push([xx + g.La, y0]); xx += g.La;
              pts.push(...rc(xx, y0, g.Ls, yin - y0).slice(1)); xx += g.Ls;
              pts.push([xx + g.L50, yin]); xx += g.L50;
              pts.push(...rc(xx, yin, g.Ls, y0 - yin).slice(1));
              paths.push({ pts, w, level: wl, kind: "mzi" });
            }
            patches.push({ x0: xc + g.margin, x1: xc + g.margin + g.Lp, y: yu, w, level: lv(z.phi), kind: "phi", core: ci, block: b, m: z.m, col: c, value: z.phi });
            const xa = xc + g.Lphi + 2 * g.Ls + g.L50 + g.margin;
            patches.push({ x0: xa, x1: xa + g.Lp, y: yu, w, level: lv(z.theta), kind: "theta", core: ci, block: b, m: z.m, col: c, value: z.theta });
          }
          for (let i = 0; i < N; i++) if (!used.has(i)) paths.push({ pts: [[xc, lane(i)], [xc + g.Lcell, lane(i)]], w, level: wl, kind: "pass" });
        }
        x += cols * g.Lcell;
      } else if (b === "sigma") {
        // spread lanes to 2P, MZI per lane against a dump guide, then gather back
        const sp = i => yc + (i - (N - 1) / 2) * 2 * P, L = sLen(Math.abs(sp(0) - lane(0)));
        for (let i = 0; i < N; i++) paths.push({ pts: rc(x, lane(i), L, sp(i) - lane(i)), w, level: wl, kind: "fan" });
        const xm = x + L;
        for (let i = 0; i < N; i++) {
          const yu = sp(i), yl = sp(i) + P, mid = (yu + yl) / 2, half = g.sc / 2;
          for (const [y0, sign, dump] of [[yu, 1, false], [yl, -1, true]]) {
            const yin = mid - sign * half; let xx = xm + (dump ? g.Lphi : 0);
            const pts = dump ? [[xx, y0]] : [[xm, y0], [xm + g.Lphi, y0]];
            pts.push(...rc(xx, y0, g.Ls, yin - y0).slice(1)); xx += g.Ls; pts.push([xx + g.L50, yin]); xx += g.L50;
            pts.push(...rc(xx, yin, g.Ls, y0 - yin).slice(1)); xx += g.Ls; pts.push([xx + g.La, y0]); xx += g.La;
            pts.push(...rc(xx, y0, g.Ls, yin - y0).slice(1)); xx += g.Ls; pts.push([xx + g.L50, yin]); xx += g.L50;
            pts.push(...rc(xx, yin, g.Ls, y0 - yin).slice(1)); xx += g.Ls;
            paths.push({ pts, w, level: wl, kind: dump ? "dump" : "mzi" });
            if (dump) for (let s = 0; s < 4; s++) paths.push({ pts: [[xx + s * 30, y0], [xx + (s + 1) * 30, y0]], w, level: wl * (0.75 - 0.2 * s), kind: "terminator" });
          }
          const s = Math.min(1, prog.S[i]), th = 2 * Math.asin(s), ph = -PI / 2 - th / 2 + (prog.blocks.meshV ? C.arg(prog.blocks.meshV.D[i]) : 0);
          patches.push({ x0: xm + g.margin, x1: xm + g.margin + g.Lp, y: sp(i), w, level: lv(ph), kind: "phi", core: ci, block: "sigma", m: i, col: 0, value: wrap(ph) });
          const xa = xm + g.Lphi + 2 * g.Ls + g.L50 + g.margin;
          patches.push({ x0: xa, x1: xa + g.Lp, y: sp(i), w, level: lv(th), kind: "theta", core: ci, block: "sigma", m: i, col: 0, value: wrap(th) });
        }
        const xg = xm + g.Lcell + 120;
        for (let i = 0; i < N; i++) { paths.push({ pts: [[xm + g.Lcell, sp(i)], [xg, sp(i)]], w, level: wl, kind: "pass" }); paths.push({ pts: rc(xg, sp(i), L, lane(i) - sp(i)), w, level: wl, kind: "fan" }); }
        x = xg + L;
      } else if (b === "phase") {
        for (let i = 0; i < N; i++) { paths.push({ pts: [[x, lane(i)], [x + g.Lphi, lane(i)]], w, level: wl, kind: "pass" }); patches.push({ x0: x + g.margin, x1: x + g.margin + g.Lp, y: lane(i), w, level: wl, kind: "phi", core: ci, block: "phase", m: i, col: 0, value: 0 }); }
        x += g.Lphi;
      }
      blocks.push({ core: ci, type: b, x0, x1: x, y0: yc - coreH / 2 + 100, y1: yc + coreH / 2 - 100 });
    }
    dieW = Math.max(dieW, x + 200);
  });
  const dieH = design.cores.length * coreH;
  for (const [mx, my] of [[100, 100], [dieW - 100, 100], [100, dieH - 100], [dieW - 100, dieH - 100]]) marks.push({ x: mx, y: my, size: 60, arm: 4 });
  return { paths, patches, marks, blocks, die: { W: dieW, H: dieH }, coreH, pixel: w / 4 };
};

// ───────────────────────── rasterizer, tiles, holograms ─────────────────────────
OTC.tilePlan = function (layout, proc) {
  const sl = proc.plate === "dmd" ? proc.dmd : proc.lcos, px = layout.pixel;
  const tw = sl.pixels[0] * px, th = sl.pixels[1] * px;
  const nx = Math.ceil(layout.die.W / tw), ny = Math.ceil(layout.die.H / th);
  return { px, tw, th, nx, ny, count: nx * ny, slmPix: sl.pixels, demag: sl.pitch_um / px, slm: sl };
};
OTC.raster = function (layout, x0, y0, nx, ny, px) {      // returns Uint8Array (0..255 write dose)
  const img = new Uint8Array(nx * ny), x1 = x0 + nx * px, y1 = y0 + ny * px;
  const stamp = (ax, ay, bx, by, hw, v) => {
    const minx = Math.max(x0, Math.min(ax, bx) - hw), maxx = Math.min(x1, Math.max(ax, bx) + hw), miny = Math.max(y0, Math.min(ay, by) - hw), maxy = Math.min(y1, Math.max(ay, by) + hw);
    if (minx >= maxx || miny >= maxy) return;
    const i0 = Math.floor((minx - x0) / px), i1 = Math.ceil((maxx - x0) / px), j0 = Math.floor((miny - y0) / px), j1 = Math.ceil((maxy - y0) / px);
    const dx = bx - ax, dy = by - ay, L2 = dx * dx + dy * dy || 1e-12, hw2 = hw * hw;
    for (let j = Math.max(0, j0); j < Math.min(ny, j1); j++) { const py = y0 + (j + 0.5) * px;
      for (let i = Math.max(0, i0); i < Math.min(nx, i1); i++) { const pxx = x0 + (i + 0.5) * px;
        let t = ((pxx - ax) * dx + (py - ay) * dy) / L2; t = t < 0 ? 0 : t > 1 ? 1 : t;
        const ex = ax + t * dx - pxx, ey = ay + t * dy - py;
        if (ex * ex + ey * ey <= hw2) { const k = j * nx + i; if (v > img[k]) img[k] = v; } } }
  };
  for (const p of layout.paths) {
    const v = Math.round(255 * p.level);
    let bx0 = Infinity, bx1 = -Infinity, by0 = Infinity, by1 = -Infinity; for (const q of p.pts) { bx0 = Math.min(bx0, q[0]); bx1 = Math.max(bx1, q[0]); by0 = Math.min(by0, q[1]); by1 = Math.max(by1, q[1]); }
    if (bx1 < x0 - p.w || bx0 > x1 + p.w || by1 < y0 - p.w || by0 > y1 + p.w) continue;
    for (let k = 1; k < p.pts.length; k++) stamp(p.pts[k - 1][0], p.pts[k - 1][1], p.pts[k][0], p.pts[k][1], p.w / 2, v);
  }
  for (const q of layout.patches) {
    const v = Math.round(255 * q.level), ya = q.y - q.w / 2, yb = q.y + q.w / 2;
    if (q.x1 < x0 || q.x0 > x1 || yb < y0 || ya > y1) continue;
    for (let j = Math.max(0, Math.floor((ya - y0) / px)); j < Math.min(ny, Math.ceil((yb - y0) / px)); j++) for (let i = Math.max(0, Math.floor((q.x0 - x0) / px)); i < Math.min(nx, Math.ceil((q.x1 - x0) / px)); i++) img[j * nx + i] = v;
  }
  for (const m of layout.marks) { stamp(m.x - m.size, m.y, m.x + m.size, m.y, m.arm, 255); stamp(m.x, m.y - m.size, m.x, m.y + m.size, m.arm, 255); }
  return img;
};
// Gerchberg–Saxton (noise-region variant) for a phase-only plate. target: Float64Array n×n (0..1), returns phase 0..2π and reconstruction.
OTC.hologram = function (target, n, opt = {}) {
  const fft = makeFFT(n), NN = n * n, iters = opt.iters || 25, levels = opt.levels || 256, fill = opt.fill || 1, r = rng(opt.seed || 7);
  const f2 = (re, im, inv) => { for (let k = 0; k < n; k++) fft(re.subarray(k * n, (k + 1) * n), im.subarray(k * n, (k + 1) * n), inv); const tr = new Float64Array(n), ti = new Float64Array(n); for (let c = 0; c < n; c++) { for (let k = 0; k < n; k++) { tr[k] = re[k * n + c]; ti[k] = im[k * n + c]; } fft(tr, ti, inv); for (let k = 0; k < n; k++) { re[k * n + c] = tr[k]; im[k * n + c] = ti[k]; } } };
  // target is placed in the centre "fill" fraction; the rest of the field is a free noise region
  const tS = new Float64Array(NN), win = new Uint8Array(NN), sh = i => (i + n / 2) % n, o = Math.floor((1 - fill) * n / 2);
  for (let a = 0; a < n; a++) for (let b = 0; b < n; b++) { const ra = sh(a), cb = sh(b), p = a * n + b; win[p] = ra >= o && ra < n - o && cb >= o && cb < n - o ? 1 : 0; tS[p] = Math.sqrt(target[ra * n + cb]); }
  let sT = 0; for (let p = 0; p < NN; p++) if (win[p]) sT += tS[p] * tS[p];
  const re = new Float64Array(NN), im = new Float64Array(NN), phi = new Float64Array(NN), q = TAU / levels;
  for (let p = 0; p < NN; p++) phi[p] = TAU * r();
  for (let it = 0; it <= iters; it++) {
    for (let p = 0; p < NN; p++) { re[p] = Math.cos(phi[p]); im[p] = Math.sin(phi[p]); }
    f2(re, im, false); if (it === iters) break;
    let pw = 0; for (let p = 0; p < NN; p++) if (win[p]) pw += re[p] * re[p] + im[p] * im[p];
    const s = Math.sqrt(pw / (sT || 1));
    for (let p = 0; p < NN; p++) if (win[p]) { const m = Math.hypot(re[p], im[p]) || 1e-30, a = s * tS[p]; re[p] = a * re[p] / m; im[p] = a * im[p] / m; }
    f2(re, im, true);
    for (let p = 0; p < NN; p++) { let a = Math.atan2(im[p], re[p]); if (a < 0) a += TAU; phi[p] = levels >= 256 ? a : (Math.round(a / q) % levels) * q; }
  }
  const I = new Float64Array(NN); for (let a = 0; a < n; a++) for (let b = 0; b < n; b++) { const p = a * n + b; I[sh(a) * n + sh(b)] = re[p] * re[p] + im[p] * im[p]; }
  // metrics
  let sI = 0, sTb = 0, nB = 0; for (let p = 0; p < NN; p++) if (target[p] > 0.3) { sI += I[p]; sTb += target[p]; nB++; }
  const scale = sI / (sTb || 1); let v = 0; for (let p = 0; p < NN; p++) if (target[p] > 0.3) { const e = I[p] / (scale * target[p]) - 1; v += e * e; }
  return { phase: phi, I, sigma: nB ? Math.sqrt(v / nB) : 0, scale };
};
// Fresnel (lensless) plate: add a quadratic lens phase so the slab sits at f = n·p²/λ from the plate.
OTC.fresnelPlate = function (phase, n, levels = 256) {
  const out = new Uint8Array(n * n), q = TAU / levels;
  for (let y = 0; y < n; y++) for (let x = 0; x < n; x++) { const p = y * n + x, rr = (x - n / 2) ** 2 + (y - n / 2) ** 2; let v = wrap(phase[p] - PI * rr / n); if (levels < 256) v = (Math.round(v / q) % levels) * q; out[p] = Math.round(v / TAU * 255) & 255; }
  return out;
};

// ───────────────────────── file writers ─────────────────────────
const CRC = (() => { const t = new Uint32Array(256); for (let n = 0; n < 256; n++) { let c = n; for (let k = 0; k < 8; k++) c = c & 1 ? 0xEDB88320 ^ (c >>> 1) : c >>> 1; t[n] = c >>> 0; } return t; })();
function crc32(buf) { let c = 0xFFFFFFFF; for (let i = 0; i < buf.length; i++) c = CRC[(c ^ buf[i]) & 0xFF] ^ (c >>> 8); return (c ^ 0xFFFFFFFF) >>> 0; }
// gray8 PNG; deflate(Uint8Array) must return zlib-wrapped data (Node: zlib.deflateSync; browser: CompressionStream("deflate"))
OTC.crc32 = crc32;
OTC.pngGray = async function (gray, w, h, deflate) {
  const raw = new Uint8Array((w + 1) * h); for (let y = 0; y < h; y++) { raw[y * (w + 1)] = 0; raw.set(gray.subarray(y * w, (y + 1) * w), y * (w + 1) + 1); }
  const z = await deflate(raw);
  const chunk = (type, data) => { const out = new Uint8Array(12 + data.length), dv = new DataView(out.buffer); dv.setUint32(0, data.length); for (let i = 0; i < 4; i++) out[4 + i] = type.charCodeAt(i); out.set(data, 8); dv.setUint32(8 + data.length, crc32(out.subarray(4, 8 + data.length))); return out; };
  const ihdr = new Uint8Array(13), dv = new DataView(ihdr.buffer); dv.setUint32(0, w); dv.setUint32(4, h); ihdr[8] = 8; ihdr[9] = 0;
  const parts = [new Uint8Array([137, 80, 78, 71, 13, 10, 26, 10]), chunk("IHDR", ihdr), chunk("IDAT", z), chunk("IEND", new Uint8Array(0))];
  const len = parts.reduce((a, p) => a + p.length, 0), png = new Uint8Array(len); let o = 0; for (const p of parts) { png.set(p, o); o += p.length; } return png;
};
// GDSII stream writer (1 nm database unit). Layers: 1 waveguides, 2 weight patches (datatype = dose level 0..255),
// 3 facet converters/terminators, 10 fiducials, 63 die outline.
OTC.gds = function (layout, name = "OGPU") {
  const bytes = []; const push = a => { for (const b of a) bytes.push(b); };
  const rec = (type, dt, data) => { const len = 4 + data.length; push([len >> 8, len & 255, type, dt]); push(data); };
  const i2 = (...v) => { const o = []; for (const x of v) o.push((x >> 8) & 255, x & 255); return o; };
  const i4 = (...v) => { const o = []; for (const x of v) { const u = x >>> 0; o.push((u >>> 24) & 255, (u >>> 16) & 255, (u >>> 8) & 255, u & 255); } return o; };
  const real8 = v => { const b = new Array(8).fill(0); if (v === 0) return b; let s = 0; if (v < 0) { s = 128; v = -v; } let e = 64; while (v >= 1) { v /= 16; e++; } while (v < 1 / 16) { v *= 16; e--; } b[0] = s | e; for (let i = 1; i < 8; i++) { v *= 256; const d = Math.floor(v); b[i] = d; v -= d; } return b; };
  const str = s => { const a = Array.from(s).map(c => c.charCodeAt(0) & 127); if (a.length % 2) a.push(0); return a; };
  const now = new Date(), date = i2(now.getFullYear(), now.getMonth() + 1, now.getDate(), now.getHours(), now.getMinutes(), now.getSeconds());
  const nm = v => Math.round(v * 1000);
  rec(0x00, 0x02, i2(600)); rec(0x01, 0x02, date.concat(date)); rec(0x02, 0x06, str(name.toUpperCase() + ".DB"));
  rec(0x03, 0x05, real8(1e-3).concat(real8(1e-9)));
  rec(0x05, 0x02, date.concat(date)); rec(0x06, 0x06, str(name.toUpperCase()));
  const path = (layer, dt, pts, width) => {
    for (let s = 0; s < pts.length - 1; s += 8000) {
      const seg = pts.slice(s, Math.min(pts.length, s + 8001));
      rec(0x09, 0x00, []); rec(0x0D, 0x02, i2(layer)); rec(0x0E, 0x02, i2(dt)); rec(0x21, 0x02, i2(0)); rec(0x0F, 0x03, i4(nm(width)));
      rec(0x10, 0x03, i4(...seg.flatMap(p => [nm(p[0]), nm(p[1])]))); rec(0x11, 0x00, []);
    }
  };
  const box = (layer, dt, x0, y0, x1, y1) => { rec(0x08, 0x00, []); rec(0x0D, 0x02, i2(layer)); rec(0x0E, 0x02, i2(dt)); rec(0x10, 0x03, i4(nm(x0), nm(y0), nm(x1), nm(y0), nm(x1), nm(y1), nm(x0), nm(y1), nm(x0), nm(y0))); rec(0x11, 0x00, []); };
  for (const p of layout.paths) path(p.kind === "converter" || p.kind === "terminator" ? 3 : 1, Math.round(255 * p.level), p.pts, p.w);
  for (const q of layout.patches) box(2, Math.round(255 * q.level), q.x0, q.y - q.w / 2, q.x1, q.y + q.w / 2);
  for (const m of layout.marks) { box(10, 0, m.x - m.size, m.y - m.arm, m.x + m.size, m.y + m.arm); box(10, 0, m.x - m.arm, m.y - m.size, m.x + m.arm, m.y + m.size); }
  box(63, 0, 0, 0, layout.die.W, layout.die.H);
  rec(0x07, 0x00, []); rec(0x04, 0x00, []);
  return new Uint8Array(bytes);
};
OTC.svg = function (layout, opt = {}) {
  const s = opt.scale || 0.05, W = layout.die.W * s, H = layout.die.H * s;
  const out = [`<svg xmlns="http://www.w3.org/2000/svg" width="${W.toFixed(0)}" height="${H.toFixed(0)}" viewBox="0 0 ${layout.die.W.toFixed(0)} ${layout.die.H.toFixed(0)}">`, `<rect width="100%" height="100%" fill="#0b1118"/>`];
  for (const p of layout.paths) out.push(`<polyline fill="none" stroke="rgb(${Math.round(60 + 195 * p.level)},${Math.round(120 + 80 * p.level)},255)" stroke-width="${Math.max(p.w, 1 / s)}" points="${p.pts.map(q => q[0].toFixed(1) + "," + q[1].toFixed(1)).join(" ")}"/>`);
  for (const q of layout.patches) out.push(`<rect x="${q.x0.toFixed(1)}" y="${(q.y - q.w).toFixed(1)}" width="${(q.x1 - q.x0).toFixed(1)}" height="${(2 * q.w).toFixed(1)}" fill="rgb(255,${Math.round(220 - 160 * (q.level - 0.5) * 2)},60)"/>`);
  for (const m of layout.marks) out.push(`<path d="M${m.x - m.size} ${m.y}H${m.x + m.size}M${m.x} ${m.y - m.size}V${m.y + m.size}" stroke="#fff" stroke-width="${2 * m.arm}"/>`);
  out.push("</svg>"); return out.join("\n");
};

// ───────────────────────── full flow ─────────────────────────
OTC.DEFAULT_SYSTEM = {
  symbolRate_GHz: 10, wavelengths: 1, laserPerCore_mW: 20, wallPlug: 0.2, tec_W: 2, fiberPitch_um: 127,
  facetLoss_dB: 2.0, dacBits: 8, adcBits: 8,
  energy_pJ: { modulator: 1.0, dac: 1.0, tia: 0.5, adc: 2.0 },
  si: { Ppi_mW: 20, mziArea_mm2: 0.03 },
  gpuRef: { name: "NVIDIA H100 SXM (INT8 dense)", tops: 1979, watts: 700 },
};
OTC.run = async function (design, opt = {}) {
  const log = opt.log || (() => {}), tick = opt.tick || (async () => {});
  const proc = Object.assign({}, OTC.PROCESS_DEFAULT, design.process || {});
  design.system = Object.assign({}, OTC.DEFAULT_SYSTEM, design.system || {});
  log("Designing MZI cell (mode solver, bend loss, coupler calibration)…", 0.05); await tick();
  const g = OTC.designCell(design.platform, proc, (design.wavelength_nm || 1550) / 1000);
  const R = opt.realizations || 10, pool = [];
  for (let k = 0; k < R; k++) { log(`Writing and measuring MZI realization ${k + 1}/${R} (BPM)…`, 0.1 + 0.5 * k / R); await tick(); pool.push(OTC.realizeMZI(g, 1000 + k * 37 + (design.seed || 0) * 101)); }
  const stats = { pool, portSigma: OTC.portPhaseSigma(g), passIL: g.P.lossDbCm * g.Lcell / 1e4 };
  log("Programming cores (SVD + Clements) and simulating accuracy…", 0.65); await tick();
  const progs = design.cores.map(c => OTC.programCore(c, design.N));
  const meanIL = pool.reduce((a, p) => a + p.IL, 0) / pool.length;
  const pathIL = 2 * design.system.facetLoss_dB + (2 * design.N + 1) * meanIL + 0.5;
  const rx = OTC.receiverSNR(design, pathIL);
  const acc = design.cores.map((c, i) => {
    const res = {};
    for (const cal of [false, true]) {
      const real = OTC.realizeCore(progs[i], stats, { calibrated: cal, levels: OTC.PLATFORMS[design.platform].levels, seed: 17 + i });
      const run = OTC.runCore(progs[i], real, { samples: opt.samples || 200, dacBits: design.system.dacBits, adcBits: design.system.adcBits, snrLin: rx.snr, seed: 5 + i });
      res[cal ? "cal" : "uncal"] = { real, run, fid: OTC.fidelity(M.fromReal(progs[i].target), real.Teff) };
    }
    return res;
  });
  log("Compiling layout…", 0.8); await tick();
  const layout = OTC.compile(design, g, progs);
  const tiles = OTC.tilePlan(layout, proc);
  const enob = Math.min(...acc.map(a => a.cal.run.enob));
  const perf = OTC.performance(design, g, layout, { enob });
  log("Done.", 1);
  return { design, proc, g, stats, progs, acc, layout, tiles, perf, rx, pathIL };
};

if (typeof module !== "undefined" && module.exports) module.exports = OTC; else root.OTC = OTC;
})(typeof globalThis !== "undefined" ? globalThis : this);
