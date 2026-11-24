// Core computations (no DOM). Port of misura_volumi.py.
import Delaunator from '../vendor/delaunator.js';

// ---------------------------------------------------------------- utilities
export function mulberry32(seed) {
  return function () {
    seed |= 0; seed = (seed + 0x6D2B79F5) | 0;
    let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

function median(arr) {
  const a = Float64Array.from(arr).sort();
  const n = a.length;
  if (!n) return NaN;
  return n % 2 ? a[(n - 1) / 2] : 0.5 * (a[n / 2 - 1] + a[n / 2]);
}

// Solve A x = b (A n×n as array of rows), Gaussian elimination with partial pivoting.
export function solve(A, b) {
  const n = b.length;
  const M = A.map((r, i) => [...r, b[i]]);
  for (let k = 0; k < n; k++) {
    let p = k;
    for (let i = k + 1; i < n; i++) if (Math.abs(M[i][k]) > Math.abs(M[p][k])) p = i;
    if (Math.abs(M[p][k]) < 1e-300) return null;
    [M[k], M[p]] = [M[p], M[k]];
    for (let i = k + 1; i < n; i++) {
      const f = M[i][k] / M[k][k];
      if (f === 0) continue;
      for (let j = k; j <= n; j++) M[i][j] -= f * M[k][j];
    }
  }
  const x = new Array(n).fill(0);
  for (let i = n - 1; i >= 0; i--) {
    let s = M[i][n];
    for (let j = i + 1; j < n; j++) s -= M[i][j] * x[j];
    x[i] = s / M[i][i];
  }
  return x;
}

// Least squares via normal equations: rows given by callback(i) -> array of m values, target t(i)
function lstsq(n, m, row, target) {
  const AtA = Array.from({ length: m }, () => new Array(m).fill(0));
  const Atb = new Array(m).fill(0);
  for (let i = 0; i < n; i++) {
    const r = row(i), t = target(i);
    for (let a = 0; a < m; a++) {
      Atb[a] += r[a] * t;
      for (let b = a; b < m; b++) AtA[a][b] += r[a] * r[b];
    }
  }
  for (let a = 0; a < m; a++) for (let b = 0; b < a; b++) AtA[a][b] = AtA[b][a];
  return solve(AtA, Atb);
}

// Symmetric 3×3 eigen decomposition (Jacobi). Returns {w:[ascending], V: columns as vectors V[k]}
export function eigSym3(S) {
  const a = S.map(r => r.slice());
  const v = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
  for (let sweep = 0; sweep < 50; sweep++) {
    const off = a[0][1] ** 2 + a[0][2] ** 2 + a[1][2] ** 2;
    if (off < 1e-30) break;
    for (const [p, q] of [[0, 1], [0, 2], [1, 2]]) {
      if (Math.abs(a[p][q]) < 1e-300) continue;
      const th = (a[q][q] - a[p][p]) / (2 * a[p][q]);
      const t = Math.sign(th || 1) / (Math.abs(th) + Math.sqrt(th * th + 1));
      const c = 1 / Math.sqrt(t * t + 1), s = t * c;
      for (let k = 0; k < 3; k++) {
        const akp = a[k][p], akq = a[k][q];
        a[k][p] = c * akp - s * akq; a[k][q] = s * akp + c * akq;
      }
      for (let k = 0; k < 3; k++) {
        const apk = a[p][k], aqk = a[q][k];
        a[p][k] = c * apk - s * aqk; a[q][k] = s * apk + c * aqk;
      }
      for (let k = 0; k < 3; k++) {
        const vkp = v[k][p], vkq = v[k][q];
        v[k][p] = c * vkp - s * vkq; v[k][q] = s * vkp + c * vkq;
      }
    }
  }
  const idx = [0, 1, 2].sort((i, j) => a[i][i] - a[j][j]);
  return { w: idx.map(i => a[i][i]), V: idx.map(i => [v[0][i], v[1][i], v[2][i]]) };
}

const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const norm = a => Math.hypot(a[0], a[1], a[2]);
const scale = (a, s) => [a[0] * s, a[1] * s, a[2] * s];

function det3(M) {
  return M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
    - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
    + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
}

export function rotationAtoB(a, b) {
  a = scale(a, 1 / norm(a)); b = scale(b, 1 / norm(b));
  const v = cross(a, b), s = norm(v), c = dot(a, b);
  if (s < 1e-12) {
    if (c > 0) return [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    let ax = cross(a, [1, 0, 0]);
    if (norm(ax) < 1e-6) ax = cross(a, [0, 1, 0]);
    ax = scale(ax, 1 / norm(ax));
    return [0, 1, 2].map(i => [0, 1, 2].map(j => 2 * ax[i] * ax[j] - (i === j ? 1 : 0)));
  }
  const vx = [[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]];
  const k = (1 - c) / (s * s);
  const R = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) {
    let vv = 0;
    for (let m = 0; m < 3; m++) vv += vx[i][m] * vx[m][j];
    R[i][j] += vx[i][j] + vv * k;
  }
  return R;
}

// ---------------------------------------------------------------- parsing
// chunks: async iterable of strings. Returns {pts: Float64Array(3N), nrm: Float32Array(3N)|null, n}
export async function parseAscChunks(chunks, totalBytes, onProgress) {
  let ncol = 0, cap = 1 << 20, buf = new Float64Array(cap), len = 0;
  let rest = '', read = 0, lastReport = 0;
  const handle = (line) => {
    let i = 0;
    while (i < line.length && (line.charCodeAt(i) === 32 || line.charCodeAt(i) === 9)) i++;
    const ch = line.charCodeAt(i);
    // digit, '-', '+', '.'
    if (!((ch >= 48 && ch <= 57) || ch === 45 || ch === 43 || ch === 46)) return;
    const tok = line.trim().split(/\s+/);
    if (!ncol) {
      if (tok.length < 3 || tok.some(t => isNaN(parseFloat(t)))) return;
      ncol = Math.min(tok.length, 6);
    }
    if (tok.length < ncol) return;
    if (len + ncol > cap) {
      cap *= 2;
      const nb = new Float64Array(cap); nb.set(buf); buf = nb;
    }
    for (let k = 0; k < ncol; k++) buf[len++] = +tok[k];
  };
  for await (const chunk of chunks) {
    read += chunk.length;
    const text = rest + chunk;
    let start = 0, nl;
    while ((nl = text.indexOf('\n', start)) !== -1) {
      handle(text.slice(start, nl));
      start = nl + 1;
    }
    rest = text.slice(start);
    if (onProgress && read - lastReport > 4e6) { lastReport = read; await onProgress(read / totalBytes); }
  }
  if (rest) handle(rest);
  if (!ncol) throw new Error('noNumeric');
  const n = Math.floor(len / ncol);
  const pts = new Float64Array(3 * n);
  let nrm = ncol >= 6 ? new Float32Array(3 * n) : null;
  let anyN = false;
  for (let i = 0; i < n; i++) {
    const o = i * ncol;
    pts[3 * i] = buf[o]; pts[3 * i + 1] = buf[o + 1]; pts[3 * i + 2] = buf[o + 2];
    if (nrm) {
      const a = buf[o + 3], b = buf[o + 4], c = buf[o + 5], l = Math.hypot(a, b, c);
      if (l > 1e-12) { nrm[3 * i] = a / l; nrm[3 * i + 1] = b / l; nrm[3 * i + 2] = c / l; anyN = true; }
    }
  }
  if (!anyN) nrm = null;
  return { pts, nrm, n };
}

// ---------------------------------------------------------------- 2D grid index
export class Grid2 {
  constructor(xs, ys, n, cell) {
    let minx = Infinity, miny = Infinity, maxx = -Infinity, maxy = -Infinity;
    for (let i = 0; i < n; i++) {
      const x = xs[i], y = ys[i];
      if (x < minx) minx = x; if (x > maxx) maxx = x;
      if (y < miny) miny = y; if (y > maxy) maxy = y;
    }
    const area = Math.max((maxx - minx) * (maxy - miny), 1e-12);
    cell = Math.max(cell, Math.sqrt(area / 4e6), 1e-9);
    this.cell = cell; this.minx = minx; this.miny = miny;
    this.nx = Math.max(1, Math.floor((maxx - minx) / cell) + 1);
    this.ny = Math.max(1, Math.floor((maxy - miny) / cell) + 1);
    this.xs = xs; this.ys = ys;
    const nc = this.nx * this.ny;
    const cnt = new Int32Array(nc + 1);
    const cid = new Int32Array(n);
    for (let i = 0; i < n; i++) {
      const c = Math.floor((xs[i] - minx) / cell) + this.nx * Math.floor((ys[i] - miny) / cell);
      cid[i] = c; cnt[c + 1]++;
    }
    for (let c = 0; c < nc; c++) cnt[c + 1] += cnt[c];
    const pos = cnt.slice(0, nc);
    const order = new Int32Array(n);
    for (let i = 0; i < n; i++) order[pos[cid[i]]++] = i;
    this.start = cnt; this.order = order;
  }
  // indices with (x-cx)²+(y-cy)² <= r²
  disc(cx, cy, r, filter) {
    const out = [];
    const { cell, minx, miny, nx, ny, xs, ys, start, order } = this;
    const i0 = Math.max(0, Math.floor((cx - r - minx) / cell)), i1 = Math.min(nx - 1, Math.floor((cx + r - minx) / cell));
    const j0 = Math.max(0, Math.floor((cy - r - miny) / cell)), j1 = Math.min(ny - 1, Math.floor((cy + r - miny) / cell));
    const r2 = r * r;
    for (let j = j0; j <= j1; j++) for (let i = i0; i <= i1; i++) {
      const c = i + nx * j;
      for (let k = start[c]; k < start[c + 1]; k++) {
        const p = order[k], dx = xs[p] - cx, dy = ys[p] - cy;
        if (dx * dx + dy * dy <= r2 && (!filter || filter(p))) out.push(p);
      }
    }
    return out;
  }
}

// ---------------------------------------------------------------- spacing, plane
export function estimateSpacing(pts, n, seed = 1) {
  // grid on the two axes with the largest range
  const rng = mulberry32(seed);
  const lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i < n; i++) for (let k = 0; k < 3; k++) {
    const v = pts[3 * i + k]; if (v < lo[k]) lo[k] = v; if (v > hi[k]) hi[k] = v;
  }
  const ext = [0, 1, 2].map(k => hi[k] - lo[k]);
  const ax = [0, 1, 2].sort((a, b) => ext[b] - ext[a]);
  const xs = new Float64Array(n), ys = new Float64Array(n);
  for (let i = 0; i < n; i++) { xs[i] = pts[3 * i + ax[0]]; ys[i] = pts[3 * i + ax[1]]; }
  let cell = Math.sqrt(Math.max(ext[ax[0]] * ext[ax[1]], 1e-12) / n) * 2;
  const g = new Grid2(xs, ys, n, cell);
  cell = g.cell;
  const d = [];
  const m = Math.min(3000, n);
  for (let s = 0; s < m; s++) {
    const i = Math.floor(rng() * n);
    let best = Infinity;
    for (const r of [cell, 2 * cell, 4 * cell]) {
      const ids = g.disc(xs[i], ys[i], r);
      for (const j of ids) {
        if (j === i) continue;
        const dd = Math.hypot(pts[3 * j] - pts[3 * i], pts[3 * j + 1] - pts[3 * i + 1], pts[3 * j + 2] - pts[3 * i + 2]);
        if (dd < best) best = dd;
      }
      if (best <= r) break;
    }
    if (isFinite(best) && best > 0) d.push(best);
  }
  return d.length ? median(d) : 1;
}

function planeFromCov(pts, idx) {
  let cx = 0, cy = 0, cz = 0;
  for (const i of idx) { cx += pts[3 * i]; cy += pts[3 * i + 1]; cz += pts[3 * i + 2]; }
  const m = idx.length; cx /= m; cy /= m; cz /= m;
  const S = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (const i of idx) {
    const d = [pts[3 * i] - cx, pts[3 * i + 1] - cy, pts[3 * i + 2] - cz];
    for (let a = 0; a < 3; a++) for (let b = a; b < 3; b++) S[a][b] += d[a] * d[b];
  }
  S[1][0] = S[0][1]; S[2][0] = S[0][2]; S[2][1] = S[1][2];
  const { V } = eigSym3(S);
  return { n: V[0], p0: [cx, cy, cz] };
}

export function fitPlaneRansac(pts, n, thr, iters = 500, maxSamples = 200000, seed = 0) {
  const rng = mulberry32(seed + 7);
  const ns = Math.min(n, maxSamples);
  const S = new Float64Array(3 * ns);
  for (let s = 0; s < ns; s++) {
    const i = ns === n ? s : Math.floor(rng() * n);
    S[3 * s] = pts[3 * i]; S[3 * s + 1] = pts[3 * i + 1]; S[3 * s + 2] = pts[3 * i + 2];
  }
  let bestCnt = -1, best = null;
  for (let it = 0; it < iters; it++) {
    const ia = Math.floor(rng() * ns), ib = Math.floor(rng() * ns), ic = Math.floor(rng() * ns);
    const a = [S[3 * ia], S[3 * ia + 1], S[3 * ia + 2]];
    const b = [S[3 * ib], S[3 * ib + 1], S[3 * ib + 2]];
    const c = [S[3 * ic], S[3 * ic + 1], S[3 * ic + 2]];
    let nv = cross([b[0] - a[0], b[1] - a[1], b[2] - a[2]], [c[0] - a[0], c[1] - a[1], c[2] - a[2]]);
    const nn = norm(nv);
    if (nn < 1e-12) continue;
    nv = scale(nv, 1 / nn);
    const d0 = -dot(nv, a);
    let cnt = 0;
    for (let s = 0; s < ns; s++) {
      if (Math.abs(nv[0] * S[3 * s] + nv[1] * S[3 * s + 1] + nv[2] * S[3 * s + 2] + d0) < thr) cnt++;
    }
    if (cnt > bestCnt) { bestCnt = cnt; best = { n: nv, p0: a }; }
  }
  if (!best) throw new Error('ransac');
  let { n: nv, p0 } = best;
  const inliers = () => {
    const idx = [];
    for (let i = 0; i < n; i++) {
      const r = (pts[3 * i] - p0[0]) * nv[0] + (pts[3 * i + 1] - p0[1]) * nv[1] + (pts[3 * i + 2] - p0[2]) * nv[2];
      if (Math.abs(r) < thr) idx.push(i);
    }
    return idx;
  };
  for (let k = 0; k < 4; k++) {
    const idx = inliers();
    if (idx.length < 3) break;
    ({ n: nv, p0 } = planeFromCov(pts, idx));
  }
  const idx = inliers();
  let ss = 0;
  for (const i of idx) {
    const r = (pts[3 * i] - p0[0]) * nv[0] + (pts[3 * i + 1] - p0[1]) * nv[1] + (pts[3 * i + 2] - p0[2]) * nv[2];
    ss += r * r;
  }
  return { n: nv, p0, inliers: idx, rms: idx.length ? Math.sqrt(ss / idx.length) : NaN };
}

// ---------------------------------------------------------------- fits
export function fitSphere(Q, m) {
  let mx = 0, my = 0, mz = 0;
  for (let i = 0; i < m; i++) { mx += Q[3 * i]; my += Q[3 * i + 1]; mz += Q[3 * i + 2]; }
  mx /= m; my /= m; mz /= m;
  const X = i => [Q[3 * i] - mx, Q[3 * i + 1] - my, Q[3 * i + 2] - mz];
  const sol = lstsq(m, 4, i => { const p = X(i); return [2 * p[0], 2 * p[1], 2 * p[2], 1]; },
    i => { const p = X(i); return p[0] * p[0] + p[1] * p[1] + p[2] * p[2]; });
  let c = [sol[0], sol[1], sol[2]];
  let r = Math.sqrt(Math.max(sol[3] + dot(c, c), 1e-12));
  for (let it = 0; it < 30; it++) {
    const JtJ = Array.from({ length: 4 }, () => [0, 0, 0, 0]), Jtr = [0, 0, 0, 0];
    for (let i = 0; i < m; i++) {
      const p = X(i), d = [p[0] - c[0], p[1] - c[1], p[2] - c[2]];
      const dist = Math.max(norm(d), 1e-12);
      const J = [-d[0] / dist, -d[1] / dist, -d[2] / dist, -1];
      const res = dist - r;
      for (let a = 0; a < 4; a++) { Jtr[a] -= J[a] * res; for (let b = 0; b < 4; b++) JtJ[a][b] += J[a] * J[b]; }
    }
    const st = solve(JtJ, Jtr);
    if (!st) break;
    c = [c[0] + st[0], c[1] + st[1], c[2] + st[2]]; r += st[3];
    if (Math.hypot(...st) < 1e-10) break;
  }
  let ss = 0;
  for (let i = 0; i < m; i++) {
    const p = X(i);
    ss += (Math.hypot(p[0] - c[0], p[1] - c[1], p[2] - c[2]) - Math.abs(r)) ** 2;
  }
  return { c: [c[0] + mx, c[1] + my, c[2] + mz], r: Math.abs(r), rms: Math.sqrt(ss / m) };
}

export function sphereVolumeAbovePlane(zc, r) {
  if (zc >= r) return 4 / 3 * Math.PI * r ** 3;
  if (zc <= -r) return 0;
  const hb = r - zc;
  return 4 / 3 * Math.PI * r ** 3 - Math.PI * hb * hb * (3 * r - hb) / 3;
}

export function fitEllipsoid(Q, m) {
  let mx = 0, my = 0, mz = 0;
  for (let i = 0; i < m; i++) { mx += Q[3 * i]; my += Q[3 * i + 1]; mz += Q[3 * i + 2]; }
  mx /= m; my /= m; mz /= m;
  let sc = 1e-9;
  for (let i = 0; i < m; i++) sc = Math.max(sc, Math.abs(Q[3 * i] - mx), Math.abs(Q[3 * i + 1] - my), Math.abs(Q[3 * i + 2] - mz));
  const X = i => [(Q[3 * i] - mx) / sc, (Q[3 * i + 1] - my) / sc, (Q[3 * i + 2] - mz) / sc];
  const p = lstsq(m, 9, i => { const [x, y, z] = X(i); return [x * x, y * y, z * z, 2 * x * y, 2 * x * z, 2 * y * z, 2 * x, 2 * y, 2 * z]; }, () => 1);
  if (!p || p.some(v => !isFinite(v))) return null;
  const A = [[p[0], p[3], p[4]], [p[3], p[1], p[5]], [p[4], p[5], p[2]]];
  const b = [p[6], p[7], p[8]];
  const cs = solve(A, b);
  if (!cs) return null;
  const c = cs.map(v => -v);
  const Ac = [dot(A[0], c), dot(A[1], c), dot(A[2], c)];
  const g = 1 + dot(c, Ac);
  if (g <= 0) return null;
  const { w, V } = eigSym3(A.map(r => r.map(v => v / g)));
  if (w.some(v => v <= 0)) return null;
  // M = U diag(1/sqrt(w)) * sc  (columns = semi-axis vectors)
  const M = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (let k = 0; k < 3; k++) for (let i = 0; i < 3; i++) M[i][k] = V[k][i] / Math.sqrt(w[k]) * sc;
  const cw = [c[0] * sc + mx, c[1] * sc + my, c[2] * sc + mz];
  // rms (approx geometric)
  const Minv = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (let k = 0; k < 3; k++) for (let i = 0; i < 3; i++) Minv[k][i] = V[k][i] * Math.sqrt(w[k]) / sc;
  let ss = 0;
  for (let i = 0; i < m; i++) {
    const d = [Q[3 * i] - cw[0], Q[3 * i + 1] - cw[1], Q[3 * i + 2] - cw[2]];
    const s = [dot(Minv[0], d), dot(Minv[1], d), dot(Minv[2], d)];
    const rn = Math.max(norm(s), 1e-12);
    ss += ((rn - 1) * norm(d) / rn) ** 2;
  }
  const axes = [0, 1, 2].map(k => Math.hypot(M[0][k], M[1][k], M[2][k])).sort((a, b) => b - a);
  return { c: cw, M, rms: Math.sqrt(ss / m), axes };
}

export function ellipsoidVolumeAbovePlane(c, M) {
  const mrow = M[2], nm = norm(mrow), detM = Math.abs(det3(M));
  if (nm < 1e-15) return c[2] >= 0 ? 4 / 3 * Math.PI * detM : 0;
  const k = -c[2] / nm;
  if (k <= -1) return 4 / 3 * Math.PI * detM;
  if (k >= 1) return 0;
  const h = 1 - k;
  return Math.PI * h * h * (3 - h) / 3 * detM;
}

// ---------------------------------------------------------------- surface & volume
function perpBasis(d) {
  const t = Math.abs(d[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0];
  let u = cross(d, t); u = scale(u, 1 / norm(u));
  return [u, cross(d, u)];
}

// Q: Float64Array 3n (x,y,h). Returns {V, tri: Uint32Array, foot}
export function closedVolume(Q, n, d) {
  d = scale(d, 1 / norm(d));
  if (d[2] < 0) d = scale(d, -1);
  const [u, v] = perpBasis(d);
  const coords = new Float64Array(2 * n);
  for (let i = 0; i < n; i++) {
    const p = [Q[3 * i], Q[3 * i + 1], Q[3 * i + 2]];
    coords[2 * i] = dot(p, u); coords[2 * i + 1] = dot(p, v);
  }
  const tri = Uint32Array.from(new Delaunator(coords).triangles);
  let V = 0, foot = 0;
  for (let t = 0; t < tri.length; t += 3) {
    let a = tri[t], b = tri[t + 1], c = tri[t + 2];
    const A = [Q[3 * a], Q[3 * a + 1], Q[3 * a + 2]], B = [Q[3 * b], Q[3 * b + 1], Q[3 * b + 2]], C = [Q[3 * c], Q[3 * c + 1], Q[3 * c + 2]];
    let cr = scale(cross([B[0] - A[0], B[1] - A[1], B[2] - A[2]], [C[0] - A[0], C[1] - A[1], C[2] - A[2]]), 0.5);
    if (dot(cr, d) < 0) { tri[t + 1] = c; tri[t + 2] = b; cr = scale(cr, -1); }
    V += (A[2] + B[2] + C[2]) / 3 * cr[2];
    if (A[2] > 0 || B[2] > 0 || C[2] > 0) foot += Math.abs(cr[2]);
  }
  return { V, tri, foot };
}

function boundaryVertices(tri) {
  const cnt = new Map();
  for (let t = 0; t < tri.length; t += 3) {
    for (const [a, b] of [[tri[t], tri[t + 1]], [tri[t + 1], tri[t + 2]], [tri[t + 2], tri[t]]]) {
      const k = a < b ? a * 4294967296 + b : b * 4294967296 + a;
      cnt.set(k, (cnt.get(k) || 0) + 1);
    }
  }
  const bnd = new Set();
  for (const [k, c] of cnt) if (c === 1) { bnd.add(Math.floor(k / 4294967296)); bnd.add(k % 4294967296); }
  return [...bnd];
}

// adds the base (z = 0) to obtain a closed mesh
export function closeMeshFaces(Q, tri) {
  const bnd = boundaryVertices(tri);
  if (bnd.length < 3) return tri;
  const coords = new Float64Array(2 * bnd.length);
  bnd.forEach((p, i) => { coords[2 * i] = Q[3 * p]; coords[2 * i + 1] = Q[3 * p + 1]; });
  const bt = new Delaunator(coords).triangles;
  const out = new Uint32Array(tri.length + bt.length);
  out.set(tri);
  for (let t = 0; t < bt.length; t += 3) {
    let a = bnd[bt[t]], b = bnd[bt[t + 1]], c = bnd[bt[t + 2]];
    const nz = (Q[3 * b] - Q[3 * a]) * (Q[3 * c + 1] - Q[3 * a + 1]) - (Q[3 * b + 1] - Q[3 * a + 1]) * (Q[3 * c] - Q[3 * a]);
    if (nz > 0) [b, c] = [c, b];
    out[tri.length + t] = a; out[tri.length + t + 1] = b; out[tri.length + t + 2] = c;
  }
  return out;
}

function largestCluster(Q, m, link) {
  const parent = new Int32Array(m).map((_, i) => i);
  const find = i => { while (parent[i] !== i) { parent[i] = parent[parent[i]]; i = parent[i]; } return i; };
  const xs = new Float64Array(m), ys = new Float64Array(m);
  for (let i = 0; i < m; i++) { xs[i] = Q[3 * i]; ys[i] = Q[3 * i + 1]; }
  const g = new Grid2(xs, ys, m, link);
  const l2 = link * link;
  for (let i = 0; i < m; i++) {
    for (const j of g.disc(xs[i], ys[i], link)) {
      if (j <= i) continue;
      const dz = Q[3 * i + 2] - Q[3 * j + 2], dx = xs[i] - xs[j], dy = ys[i] - ys[j];
      if (dx * dx + dy * dy + dz * dz <= l2) {
        const a = find(i), b = find(j);
        if (a !== b) parent[a] = b;
      }
    }
  }
  const cnt = new Map();
  for (let i = 0; i < m; i++) { const r = find(i); cnt.set(r, (cnt.get(r) || 0) + 1); }
  let best = -1, bc = -1;
  for (const [r, c] of cnt) if (c > bc) { bc = c; best = r; }
  const keep = new Uint8Array(m);
  for (let i = 0; i < m; i++) keep[i] = find(i) === best ? 1 : 0;
  return keep;
}

// ---------------------------------------------------------------- model
export class CloudModel {
  constructor(pts, nrm, n) {
    this.raw = pts; this.rawN = nrm; this.n = n;
    this.spacing = estimateSpacing(pts, n);
    this.flip = false;
  }

  align(thr, flip = false) {
    const { raw, rawN, n } = this;
    this.thr = thr; this.flip = flip;
    const pl = fitPlaneRansac(raw, n, thr);
    let nv = pl.n;
    if (rawN) {
      const m = [0, 0, 0];
      for (let i = 0; i < n; i++) { m[0] += rawN[3 * i]; m[1] += rawN[3 * i + 1]; m[2] += rawN[3 * i + 2]; }
      if (norm(m) > 1e-9 && dot(nv, m) < 0) nv = scale(nv, -1);
    } else {
      const isIn = new Uint8Array(n); for (const i of pl.inliers) isIn[i] = 1;
      const out = [];
      for (let i = 0; i < n; i++) if (!isIn[i]) out.push((raw[3 * i] - pl.p0[0]) * nv[0] + (raw[3 * i + 1] - pl.p0[1]) * nv[1] + (raw[3 * i + 2] - pl.p0[2]) * nv[2]);
      if (out.length && median(out) < 0) nv = scale(nv, -1);
    }
    if (flip) nv = scale(nv, -1);
    const R = rotationAtoB(nv, [0, 0, 1]);
    const o = [0, 0, 0];
    for (const i of pl.inliers) { o[0] += raw[3 * i]; o[1] += raw[3 * i + 1]; o[2] += raw[3 * i + 2]; }
    for (let k = 0; k < 3; k++) o[k] /= pl.inliers.length;
    const x = new Float64Array(n), y = new Float64Array(n), z = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      const d = [raw[3 * i] - o[0], raw[3 * i + 1] - o[1], raw[3 * i + 2] - o[2]];
      x[i] = dot(R[0], d); y[i] = dot(R[1], d); z[i] = dot(R[2], d);
    }
    let zm = 0; for (const i of pl.inliers) zm += z[i]; zm /= pl.inliers.length;
    for (let i = 0; i < n; i++) z[i] -= zm;
    let N = null;
    if (rawN) {
      N = new Float32Array(3 * n);
      for (let i = 0; i < n; i++) {
        const d = [rawN[3 * i], rawN[3 * i + 1], rawN[3 * i + 2]];
        N[3 * i] = dot(R[0], d); N[3 * i + 1] = dot(R[1], d); N[3 * i + 2] = dot(R[2], d);
      }
    }
    Object.assign(this, { x, y, z, N, R, origin: o, nInliers: pl.inliers.length, rms: pl.rms });
    this.grid = new Grid2(x, y, n, this.spacing * 4);
  }

  // aligned -> original coordinates
  toOriginal(p) {
    const R = this.R, o = this.origin;
    return [0, 1, 2].map(k => R[0][k] * p[0] + R[1][k] * p[1] + R[2][k] * p[2] + o[k]);
  }

  select(c, r, cylinder) {
    const { x, y, z } = this;
    if (cylinder) return this.grid.disc(c[0], c[1], r);
    const r2 = r * r;
    return this.grid.disc(c[0], c[1], r, p => (x[p] - c[0]) ** 2 + (y[p] - c[1]) ** 2 + (z[p] - c[2]) ** 2 <= r2);
  }

  measure(c, r, tol, { cylinder = false, keepLargest = true, ringFactor = 1.5, voxel = 0 } = {}) {
    const { x, y, z, N, spacing } = this;
    const warnings = [];
    const idx = this.select(c, r, cylinder);
    if (idx.length < 10) throw new Error('fewPoints|' + idx.length);
    let rd = 0;
    for (const p of idx) rd = Math.max(rd, Math.hypot(x[p] - c[0], y[p] - c[1]));
    rd = Math.max(rd, 3 * spacing);

    // local plane on a ring around the probe (coordinates centred on c)
    const ringR = rd * ringFactor + 3 * spacing;
    const zlim = Math.max(10 * tol, 3 * this.thr);
    const cand = this.grid.disc(c[0], c[1], ringR, p => Math.hypot(x[p] - c[0], y[p] - c[1]) > rd && Math.abs(z[p]) < zlim);
    let coef = [0, 0, 0], lrms = NaN, nring = 0;
    if (cand.length >= 20) {
      let mask = new Uint8Array(cand.length).fill(1);
      for (let it = 0; it < 6; it++) {
        const sub = cand.filter((_, k) => mask[k]);
        if (sub.length < 10) break;
        const s = lstsq(sub.length, 3, k => [x[sub[k]] - c[0], y[sub[k]] - c[1], 1], k => z[sub[k]]);
        if (!s) break;
        coef = s;
        const res = cand.map(p => z[p] - (coef[0] * (x[p] - c[0]) + coef[1] * (y[p] - c[1]) + coef[2]));
        const rm = res.filter((_, k) => mask[k]);
        const med = median(rm);
        const mad = 1.4826 * median(rm.map(v => Math.abs(v - med)));
        const lim = Math.max(3 * mad, tol);
        const nm = Uint8Array.from(res, v => (Math.abs(v) < lim ? 1 : 0));
        const same = nm.every((v, k) => v === mask[k]);
        mask = nm;
        if (same) break;
      }
      let ss = 0; nring = 0;
      cand.forEach((p, k) => {
        if (!mask[k]) return;
        const rr = z[p] - (coef[0] * (x[p] - c[0]) + coef[1] * (y[p] - c[1]) + coef[2]);
        ss += rr * rr; nring++;
      });
      lrms = nring ? Math.sqrt(ss / nring) : NaN;
    } else {
      warnings.push(['wNoRing']);
    }
    const hgt = p => z[p] - (coef[0] * (x[p] - c[0]) + coef[1] * (y[p] - c[1]) + coef[2]);

    let objIds = idx.filter(p => hgt(p) > tol);
    if (objIds.length < 10) throw new Error('noObject|' + tol);
    let nExclCl = 0;
    if (keepLargest) {
      const Qt = new Float64Array(3 * objIds.length);
      objIds.forEach((p, k) => { Qt[3 * k] = x[p]; Qt[3 * k + 1] = y[p]; Qt[3 * k + 2] = hgt(p); });
      const keep = largestCluster(Qt, objIds.length, 4 * spacing);
      const before = objIds.length;
      objIds = objIds.filter((_, k) => keep[k]);
      nExclCl = before - objIds.length;
    }
    const no = objIds.length;
    let hmax = -Infinity;
    for (const p of objIds) hmax = Math.max(hmax, hgt(p));

    // scan direction from the file "normals"
    let dvec = [0, 0, 1];
    if (N) {
      const m = [0, 0, 0];
      for (const p of objIds) { m[0] += N[3 * p]; m[1] += N[3 * p + 1]; m[2] += N[3 * p + 2]; }
      const l = norm(m);
      if (l > 1e-9 && Math.abs(m[2]) / l > 0.05) dvec = scale(m, 1 / l);
    }
    if (dvec[2] < 0) dvec = scale(dvec, -1);
    const tanT = Math.hypot(dvec[0], dvec[1]) / dvec[2];

    // plane disc around the object, wide enough to contain its shadow
    let robjMax = 0;
    for (const p of objIds) robjMax = Math.max(robjMax, Math.hypot(x[p] - c[0], y[p] - c[1]));
    const Rb = robjMax + hmax * tanT + 6 * spacing;
    const disc = this.grid.disc(c[0], c[1], Rb);
    const isObj = new Set(objIds);
    const base = [], removed = [];
    let nAbove = 0;
    for (const p of disc) {
      const h = hgt(p);
      if (h <= tol) base.push(p);
      else nAbove++;
      if (!isObj.has(p)) removed.push(x[p], y[p], h);
    }
    const nExclOut = nAbove - no;
    if (base.length < 3) warnings.push(['wFewBase']);

    // checks
    if (!cylinder) {
      let ne = 0, hmaxE = 0;
      for (const p of objIds) {
        if (Math.hypot(x[p] - c[0], y[p] - c[1], z[p] - c[2]) > r - 2 * spacing) { ne++; hmaxE = Math.max(hmaxE, hgt(p)); }
      }
      if (ne > 0.01 * no && hmaxE > 5 * tol) warnings.push(['wEdgeSphere']);
    } else {
      let hmaxE = 0;
      for (const p of objIds) if (Math.hypot(x[p] - c[0], y[p] - c[1]) > r - 2 * spacing) hmaxE = Math.max(hmaxE, hgt(p));
      if (hmaxE > 5 * tol) warnings.push(['wEdgeCyl']);
    }
    if (nExclOut > 0.05 * no) warnings.push(['wExclOut', nExclOut]);
    if (nExclCl > 0.03 * no + 20) warnings.push(['wExclCl', nExclCl]);

    // volume
    const nq = no + base.length;
    const Q = new Float64Array(3 * nq);
    objIds.forEach((p, k) => { Q[3 * k] = x[p]; Q[3 * k + 1] = y[p]; Q[3 * k + 2] = hgt(p); });
    base.forEach((p, k) => { Q[3 * (no + k)] = x[p]; Q[3 * (no + k) + 1] = y[p]; Q[3 * (no + k) + 2] = 0; });
    const cvScan = closedVolume(Q, nq, dvec);
    const cvCol = closedVolume(Q, nq, [0, 0, 1]);
    const bnd = boundaryVertices(cvScan.tri);
    if (bnd.some(p => Q[3 * p + 2] > 0)) warnings.push(['wOpen']);

    const vs = voxel > 0 ? voxel : autoVoxelSize(cvScan.V, spacing);
    const vx = voxelVolume(Q, nq, cvScan.tri, dvec, vs);
    if (vx && vx.size > vs * 1.01) warnings.push(['wVoxGrow', vx.size.toFixed(4)]);

    const sph = fitSphere(Q, no);
    const Vsph = sphereVolumeAbovePlane(sph.c[2], sph.r);
    const ell = fitEllipsoid(Q, no);
    const Vell = ell ? ellipsoidVolumeAbovePlane(ell.c, ell.M) : NaN;

    // recommended estimate: undercut?
    let fit = { name: 'sphere', V: Vsph, rms: sph.rms, zc: sph.c[2], size: sph.r };
    if (ell && ell.rms < 0.8 * sph.rms && isFinite(Vell)) fit = { name: 'ellipsoid', V: Vell, rms: ell.rms, zc: ell.c[2], size: ell.axes[2] };
    const goodFit = fit.rms < 0.05 * fit.size + 2 * tol;
    const undercut = goodFit && fit.zc > 0.15 * hmax;
    const Vbest = undercut ? fit.V : cvScan.V;

    const objOrig = new Float64Array(3 * no), objN = this.rawN ? new Float32Array(3 * no) : null;
    objIds.forEach((p, k) => {
      for (let j = 0; j < 3; j++) objOrig[3 * k + j] = this.raw[3 * p + j];
      if (objN) for (let j = 0; j < 3; j++) objN[3 * k + j] = this.rawN[3 * p + j];
    });

    return {
      center: c.slice(), r, cylinder, tol,
      Q, nq, nObj: no, nBase: base.length, tri: cvScan.tri,
      Vscan: cvScan.V, Vcol: cvCol.V, vox: vx, scanDir: dvec, Vsphere: Vsph, Vell, Vbest, undercut, fitName: fit.name,
      sphere: sph, ell, footprint: cvScan.foot, hmax,
      scanAngle: Math.atan(tanT) * 180 / Math.PI, localTilt: Math.atan(Math.hypot(coef[0], coef[1])) * 180 / Math.PI,
      localRms: lrms, nRing: nring, warnings,
      removed: Float64Array.from(removed), objOrig, objN,
    };
  }
}

// ---------------------------------------------------------------- voxel volume
// The object is the region bounded by the triangulated sheet (graph along the scan
// direction d) and the plane z = 0: a point p is inside iff z >= 0 and p·d <= sheet(p).
// Voxels whose 8 corners are all inside count fully, empty ones count zero, and
// boundary voxels are supersampled (m³ points) to get their filled fraction.
export function voxelVolume(Q, nq, tri, d, size, { sub = 4, maxCorners = 12e6 } = {}) {
  d = scale(d, 1 / norm(d));
  if (d[2] < 0) d = scale(d, -1);
  const [u, v] = perpBasis(d);
  const U = new Float64Array(nq), Vv = new Float64Array(nq), W = new Float64Array(nq);
  for (let i = 0; i < nq; i++) {
    const p = [Q[3 * i], Q[3 * i + 1], Q[3 * i + 2]];
    U[i] = dot(p, u); Vv[i] = dot(p, v); W[i] = dot(p, d);
  }
  // triangles of the raised part only (the rest lies on the plane) + bbox of the region
  const raised = [];
  let lo = [Infinity, Infinity, 0], hi = [-Infinity, -Infinity, 0];
  for (let t = 0; t < tri.length; t += 3) {
    const a = tri[t], b = tri[t + 1], c = tri[t + 2];
    if (Q[3 * a + 2] > 0 || Q[3 * b + 2] > 0 || Q[3 * c + 2] > 0) {
      raised.push(t);
      for (const p of [a, b, c]) {
        lo[0] = Math.min(lo[0], Q[3 * p]); hi[0] = Math.max(hi[0], Q[3 * p]);
        lo[1] = Math.min(lo[1], Q[3 * p + 1]); hi[1] = Math.max(hi[1], Q[3 * p + 1]);
        hi[2] = Math.max(hi[2], Q[3 * p + 2]);
      }
    }
  }
  if (!raised.length) return null;
  // grow size until the corner grid is small enough
  let s = size;
  const dims = () => [0, 1, 2].map(k => Math.max(1, Math.ceil((hi[k] - lo[k]) / s)));
  while (dims().reduce((a, b) => a * (b + 1), 1) > maxCorners) s *= 1.25;
  const [nx, ny, nz] = dims();
  // bucket grid of raised triangles in (u,v)
  let ulo = Infinity, uhi = -Infinity, vlo = Infinity, vhi = -Infinity;
  for (const t of raised) for (let k = 0; k < 3; k++) {
    const p = tri[t + k];
    ulo = Math.min(ulo, U[p]); uhi = Math.max(uhi, U[p]); vlo = Math.min(vlo, Vv[p]); vhi = Math.max(vhi, Vv[p]);
  }
  const bc = Math.max(Math.sqrt((uhi - ulo) * (vhi - vlo) / raised.length) * 1.5, 1e-9);
  const bnx = Math.max(1, Math.ceil((uhi - ulo) / bc) + 1), bny = Math.max(1, Math.ceil((vhi - vlo) / bc) + 1);
  const buckets = Array.from({ length: bnx * bny }, () => []);
  for (const t of raised) {
    const ps = [tri[t], tri[t + 1], tri[t + 2]];
    const i0 = Math.floor((Math.min(...ps.map(p => U[p])) - ulo) / bc), i1 = Math.floor((Math.max(...ps.map(p => U[p])) - ulo) / bc);
    const j0 = Math.floor((Math.min(...ps.map(p => Vv[p])) - vlo) / bc), j1 = Math.floor((Math.max(...ps.map(p => Vv[p])) - vlo) / bc);
    for (let j = j0; j <= j1; j++) for (let i = i0; i <= i1; i++) buckets[i + bnx * j].push(t);
  }
  const eps = 1e-12;
  const inside = (x, y, z) => {
    if (z < 0) return false;
    const pu = x * u[0] + y * u[1] + z * u[2], pv = x * v[0] + y * v[1] + z * v[2];
    const i = Math.floor((pu - ulo) / bc), j = Math.floor((pv - vlo) / bc);
    if (i < 0 || j < 0 || i >= bnx || j >= bny) return false;
    for (const t of buckets[i + bnx * j]) {
      const a = tri[t], b = tri[t + 1], c = tri[t + 2];
      const x1 = U[a], y1 = Vv[a], x2 = U[b], y2 = Vv[b], x3 = U[c], y3 = Vv[c];
      const den = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3);
      if (Math.abs(den) < 1e-300) continue;
      const l1 = ((y2 - y3) * (pu - x3) + (x3 - x2) * (pv - y3)) / den;
      const l2 = ((y3 - y1) * (pu - x3) + (x1 - x3) * (pv - y3)) / den;
      const l3 = 1 - l1 - l2;
      if (l1 >= -eps && l2 >= -eps && l3 >= -eps) {
        const ws = l1 * W[a] + l2 * W[b] + l3 * W[c];
        return x * d[0] + y * d[1] + z * d[2] <= ws;
      }
    }
    return false; // over the flat part of the plane or outside: not solid
  };
  // corner flags
  const cx = nx + 1, cy = ny + 1, cz = nz + 1;
  const corner = new Uint8Array(cx * cy * cz);
  for (let k = 0; k < cz; k++) for (let j = 0; j < cy; j++) for (let i = 0; i < cx; i++) {
    corner[i + cx * (j + cy * k)] = inside(lo[0] + i * s, lo[1] + j * s, k * s) ? 1 : 0;
  }
  const vox = s * s * s;
  let nFull = 0, nPart = 0, fracSum = 0;
  const occ = []; // [i, j, k, frac] of non-empty voxels, for display
  for (let k = 0; k < nz; k++) for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
    let c = 0;
    for (let dk = 0; dk < 2; dk++) for (let dj = 0; dj < 2; dj++) for (let di = 0; di < 2; di++)
      c += corner[(i + di) + cx * ((j + dj) + cy * (k + dk))];
    if (c === 0) continue;
    if (c === 8) { nFull++; occ.push(i, j, k, 1); continue; }
    let f = 0;
    for (let c3 = 0; c3 < sub; c3++) for (let b3 = 0; b3 < sub; b3++) for (let a3 = 0; a3 < sub; a3++) {
      if (inside(lo[0] + (i + (a3 + 0.5) / sub) * s, lo[1] + (j + (b3 + 0.5) / sub) * s, (k + (c3 + 0.5) / sub) * s)) f++;
    }
    f /= sub * sub * sub;
    nPart++; fracSum += f;
    if (f > 0) occ.push(i, j, k, f);
  }
  return {
    V: (nFull + fracSum) * vox, Vmin: nFull * vox, Vmax: (nFull + nPart) * vox,
    size: s, dims: [nx, ny, nz], nFull, nPart, origin: [lo[0], lo[1], 0], occ: Float32Array.from(occ),
  };
}

export function autoVoxelSize(Vscan, spacing) {
  return Math.max(Math.cbrt(Math.max(Vscan, 1e-12)) / 40, spacing / 4);
}

// ---------------------------------------------------------------- export helpers
export function ascText(pts, nrm, m) {
  const out = [];
  out.push(nrm ? 'Orig_Pos[X] Orig_Pos[Y] Orig_Pos[Z] Orig_Vec[I] Orig_Vec[J] Orig_Vec[K] ' : 'Orig_Pos[X] Orig_Pos[Y] Orig_Pos[Z] ');
  if (nrm) out.push(`SCANLINE 0  ( n = ${m} ) `);
  for (let i = 0; i < m; i++) {
    let s = `${pts[3 * i].toFixed(8)} ${pts[3 * i + 1].toFixed(8)} ${pts[3 * i + 2].toFixed(8)}`;
    if (nrm) s += ` ${nrm[3 * i].toFixed(8)} ${nrm[3 * i + 1].toFixed(8)} ${nrm[3 * i + 2].toFixed(8)}`;
    out.push(s);
  }
  return out.join('\n') + '\n';
}

export function stlBinary(Q, faces) {
  const nt = faces.length / 3;
  const buf = new ArrayBuffer(84 + nt * 50);
  const dv = new DataView(buf);
  dv.setUint32(80, nt, true);
  let o = 84;
  for (let t = 0; t < nt; t++) {
    const a = faces[3 * t], b = faces[3 * t + 1], c = faces[3 * t + 2];
    const A = [Q[3 * a], Q[3 * a + 1], Q[3 * a + 2]], B = [Q[3 * b], Q[3 * b + 1], Q[3 * b + 2]], C = [Q[3 * c], Q[3 * c + 1], Q[3 * c + 2]];
    let nv = cross([B[0] - A[0], B[1] - A[1], B[2] - A[2]], [C[0] - A[0], C[1] - A[1], C[2] - A[2]]);
    const l = norm(nv) || 1; nv = scale(nv, 1 / l);
    for (const v of [nv, A, B, C]) for (let k = 0; k < 3; k++) { dv.setFloat32(o, v[k], true); o += 4; }
    o += 2;
  }
  return buf;
}

// signed volume of a closed mesh (check)
export function meshVolume(Q, faces) {
  let V = 0;
  for (let t = 0; t < faces.length; t += 3) {
    const a = faces[t], b = faces[t + 1], c = faces[t + 2];
    V += dot([Q[3 * a], Q[3 * a + 1], Q[3 * a + 2]], cross([Q[3 * b], Q[3 * b + 1], Q[3 * b + 2]], [Q[3 * c], Q[3 * c + 1], Q[3 * c + 2]])) / 6;
  }
  return V;
}
