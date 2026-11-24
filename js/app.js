// UI for the volume measurement page
import * as THREE from '../vendor/three.module.min.js';
import { OrbitControls } from '../vendor/OrbitControls.js';
import * as core from './core.js';
import { t, setLang, getLang } from './i18n.js';

const $ = id => document.getElementById(id);
const state = { model: null, fname: '', measures: [], selRow: -1, cur: null, curIdx: -1, nSel: null };
const frame = () => new Promise(r => setTimeout(r, 16));
const f4 = v => (isFinite(v) ? v.toFixed(4) : '—');

// ---------------------------------------------------------------- colours (viridis)
const VIR = [[68, 1, 84], [72, 40, 120], [62, 73, 137], [49, 104, 142], [38, 130, 142], [31, 158, 137], [53, 183, 121], [110, 206, 88], [181, 222, 43], [253, 231, 37]];
function viridis(t01, out, o) {
  const tt = Math.min(1, Math.max(0, t01)) * (VIR.length - 1);
  const i = Math.min(VIR.length - 2, Math.floor(tt)), f = tt - i;
  for (let k = 0; k < 3; k++) out[o + k] = (VIR[i][k] + (VIR[i + 1][k] - VIR[i][k]) * f) / 255;
}

// ---------------------------------------------------------------- viewer
class Viewer {
  constructor(container) {
    this.container = container;
    this.renderer = new THREE.WebGLRenderer({ antialias: true, preserveDrawingBuffer: true });
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
    this.renderer.setClearColor(0xffffff);
    container.prepend(this.renderer.domElement);
    this.scene = new THREE.Scene();
    this.camera = new THREE.PerspectiveCamera(40, 1, 0.01, 1e5);
    this.camera.up.set(0, 0, 1);
    this.scene.add(this.camera);
    this.scene.add(new THREE.AmbientLight(0xffffff, 0.75));
    const dl = new THREE.DirectionalLight(0xffffff, 1.6);
    dl.position.set(0.3, 0.4, 1);
    this.camera.add(dl);
    this.controls = new OrbitControls(this.camera, this.renderer.domElement);
    this.controls.addEventListener('change', () => this.render());
    this.controls.zoomToCursor = true;
    new ResizeObserver(() => this.resize()).observe(container);
    this.resize();
  }
  resize() {
    const w = this.container.clientWidth, h = this.container.clientHeight;
    if (!w || !h) return;
    this.renderer.setSize(w, h, false);
    this.camera.aspect = w / h;
    this.camera.updateProjectionMatrix();
    this.render();
  }
  render() { this.renderer.render(this.scene, this.camera); }
  fit(box) {
    const c = box.getCenter(new THREE.Vector3()), s = box.getSize(new THREE.Vector3());
    const d = Math.max(s.length(), 1e-3);
    this.camera.near = d / 1000; this.camera.far = d * 100;
    this.camera.position.copy(c).add(new THREE.Vector3(0.9, -1.1, 0.95).normalize().multiplyScalar(d * 1.25));
    this.camera.updateProjectionMatrix();
    this.controls.target.copy(c);
    this.controls.update();
    this.render();
  }
}

// ---------------------------------------------------------------- i18n
function applyLang() {
  document.documentElement.lang = getLang();
  document.title = t('title');
  document.querySelectorAll('[data-i18n]').forEach(el => { el.innerHTML = t(el.dataset.i18n); });
  $('lang').value = getLang();
  updateInfo();
  updateSelCount();
  if (state.cur) renderResultText();
  if (state.cur && !$('result').hidden && rviewer) buildResultScene(false);
}

// ---------------------------------------------------------------- busy overlay
function busy(msg, frac) {
  if (msg === null) { $('busy').hidden = true; return; }
  $('busy').hidden = false;
  $('busyMsg').textContent = msg;
  $('busyBar').style.width = frac == null ? '100%' : Math.round(frac * 100) + '%';
}

// ---------------------------------------------------------------- main view
let viewer;
const scn = {};
function initMainView() {
  viewer = new Viewer($('view'));
  scn.probe = new THREE.Mesh(new THREE.SphereGeometry(1, 48, 32),
    new THREE.MeshPhongMaterial({ color: 0xe0302a, transparent: true, opacity: 0.25, depthWrite: false }));
  scn.probe.renderOrder = 2;
  const cp = [];
  for (let i = 0; i <= 128; i++) { const a = i / 128 * Math.PI * 2; cp.push(Math.cos(a), Math.sin(a), 0); }
  const cg = new THREE.BufferGeometry(); cg.setAttribute('position', new THREE.Float32BufferAttribute(cp, 3));
  scn.circle = new THREE.Line(cg, new THREE.LineBasicMaterial({ color: 0xe0302a }));
  scn.selGeom = new THREE.BufferGeometry();
  scn.sel = new THREE.Points(scn.selGeom, new THREE.PointsMaterial({ color: 0xff1a1a, size: 3.5, sizeAttenuation: false }));
  scn.sel.frustumCulled = false; scn.sel.renderOrder = 1;
  scn.axes = new THREE.AxesHelper(1);
  for (const o of [scn.probe, scn.circle, scn.sel, scn.axes]) { o.visible = false; viewer.scene.add(o); }
}

function buildCloud(resetCam) {
  const m = state.model;
  if (scn.cloud) { viewer.scene.remove(scn.cloud); scn.cloud.geometry.dispose(); }
  const stride = Math.max(1, Math.ceil(m.n / 8e6));
  const nn = Math.ceil(m.n / stride);
  const pos = new Float32Array(3 * nn), col = new Float32Array(3 * nn);
  // colour limits: percentiles of z on a sample
  const smp = [];
  for (let i = 0; i < m.n; i += Math.max(1, Math.floor(m.n / 100000))) smp.push(m.z[i]);
  smp.sort((a, b) => a - b);
  const lo = smp[Math.floor(0.01 * (smp.length - 1))], hi = Math.max(smp[Math.floor(0.998 * (smp.length - 1))], lo + 1e-6);
  let k = 0;
  for (let i = 0; i < m.n; i += stride, k++) {
    pos[3 * k] = m.x[i]; pos[3 * k + 1] = m.y[i]; pos[3 * k + 2] = m.z[i];
    viridis((m.z[i] - lo) / (hi - lo), col, 3 * k);
  }
  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.BufferAttribute(pos, 3));
  g.setAttribute('color', new THREE.BufferAttribute(col, 3));
  g.computeBoundingBox(); g.computeBoundingSphere();
  scn.cloud = new THREE.Points(g, new THREE.PointsMaterial({ size: 2, sizeAttenuation: false, vertexColors: true }));
  viewer.scene.add(scn.cloud);
  $('legend').hidden = false; $('lgMin').textContent = lo.toFixed(3); $('lgMax').textContent = hi.toFixed(3);
  const ext = g.boundingBox.getSize(new THREE.Vector3());
  scn.axes.scale.setScalar(Math.max(ext.x, ext.y) * 0.08);
  for (const o of [scn.probe, scn.circle, scn.sel, scn.axes]) o.visible = true;
  if (resetCam) viewer.fit(g.boundingBox);
}

// ---------------------------------------------------------------- probe sliders
const SL = {};
for (const k of ['x', 'y', 'z', 'r']) SL[k] = { range: $('r' + k), num: $('n' + k), min: 0, max: 1, v: 0 };
function setRange(k, min, max) {
  if (!(max > min)) max = min + 1e-6;
  Object.assign(SL[k], { min, max });
  SL[k].num.min = min; SL[k].num.max = max; SL[k].num.step = +((max - min) / 400).toPrecision(2);
  setVal(k, SL[k].v);
}
function setVal(k, v) {
  const s = SL[k];
  if (!isFinite(v)) v = s.v;
  v = Math.min(s.max, Math.max(s.min, v));
  s.v = v; s.num.value = +v.toFixed(4);
  s.range.value = Math.round((v - s.min) / (s.max - s.min) * 4000);
}
const val = k => SL[k].v;
for (const k of Object.keys(SL)) {
  SL[k].range.addEventListener('input', () => { setVal(k, SL[k].min + (SL[k].max - SL[k].min) * SL[k].range.value / 4000); probeChanged(); });
  SL[k].num.addEventListener('change', () => { setVal(k, parseFloat(SL[k].num.value)); probeChanged(); });
}
['cbRest', 'cbCyl', 'cbShow'].forEach(id => $(id).addEventListener('change', probeChanged));

function probeCenter() { return [val('x'), val('y'), $('cbRest').checked ? val('r') : val('z')]; }

let selTimer = null;
function probeChanged() {
  if (!state.model) return;
  const rest = $('cbRest').checked, cyl = $('cbCyl').checked;
  if (rest) setVal('z', val('r'));
  SL.z.range.disabled = SL.z.num.disabled = rest || cyl;
  const c = probeCenter(), r = val('r');
  scn.probe.visible = !cyl;
  scn.probe.position.set(...c); scn.probe.scale.setScalar(r);
  const rc = cyl ? r : (Math.sqrt(Math.max(r * r - c[2] * c[2], 0)) || r * 0.02);
  scn.circle.position.set(c[0], c[1], 0); scn.circle.scale.set(rc, rc, 1);
  viewer.render();
  clearTimeout(selTimer);
  selTimer = setTimeout(updateSelection, 60);
}

function updateSelection() {
  const m = state.model;
  if (!m) return;
  const idx = m.select(probeCenter(), val('r'), $('cbCyl').checked);
  state.nSel = idx.length;
  updateSelCount();
  if ($('cbShow').checked && idx.length) {
    const st = Math.max(1, Math.ceil(idx.length / 400000));
    const nn = Math.ceil(idx.length / st), pos = new Float32Array(3 * nn);
    for (let a = 0, k = 0; a < idx.length; a += st, k++) { const p = idx[a]; pos[3 * k] = m.x[p]; pos[3 * k + 1] = m.y[p]; pos[3 * k + 2] = m.z[p]; }
    scn.selGeom.setAttribute('position', new THREE.BufferAttribute(pos, 3));
    scn.sel.visible = true;
  } else scn.sel.visible = false;
  viewer.render();
}
function updateSelCount() { $('selCount').textContent = state.nSel == null ? '' : t('selCount', { n: state.nSel.toLocaleString() }); }

function nudge(dx, dy, dz, dr) {
  if (!state.model) return;
  const st = val('r') / 5;
  if (dx) setVal('x', val('x') + dx * st);
  if (dy) setVal('y', val('y') + dy * st);
  if (dz) { if ($('cbRest').checked) { $('cbRest').checked = false; setVal('z', val('r')); } setVal('z', val('z') + dz * st); }
  if (dr) setVal('r', val('r') * (dr > 0 ? 1.1 : 1 / 1.1));
  probeChanged();
}

// ---------------------------------------------------------------- picking
const mouse = { x: 0, y: 0, inside: false };
const raycaster = new THREE.Raycaster();
function pickAt(clientX, clientY) {
  if (!state.model || !scn.cloud) return;
  const rect = viewer.renderer.domElement.getBoundingClientRect();
  const ndc = new THREE.Vector2(((clientX - rect.left) / rect.width) * 2 - 1, -((clientY - rect.top) / rect.height) * 2 + 1);
  raycaster.setFromCamera(ndc, viewer.camera);
  // threshold ≈ a few pixels at the target distance
  const dist = viewer.camera.position.distanceTo(viewer.controls.target);
  const pix = 2 * dist * Math.tan(THREE.MathUtils.degToRad(viewer.camera.fov / 2)) / rect.height;
  raycaster.params.Points.threshold = Math.max(state.model.spacing * 1.5, pix * 4);
  const hit = raycaster.intersectObject(scn.cloud, false)[0];
  if (!hit) return;
  const a = scn.cloud.geometry.attributes.position, i = hit.index;
  const p = [a.getX(i), a.getY(i), a.getZ(i)];
  setVal('x', p[0]); setVal('y', p[1]);
  if (!$('cbRest').checked) setVal('z', Math.max(p[2], 0) * 0.5 + val('r') * 0.5);
  probeChanged();
}

// ---------------------------------------------------------------- loading
async function* streamChunks(stream) {
  const reader = stream.getReader();
  for (;;) { const { done, value } = await reader.read(); if (done) return; yield value; }
}

async function loadFile(file) {
  try {
    busy(t('reading', { p: 0 }), 0);
    await frame();
    const stream = file.stream().pipeThrough(new TextDecoderStream());
    const { pts, nrm, n } = await core.parseAscChunks(streamChunks(stream), file.size, async f => {
      busy(t('reading', { p: Math.min(100, Math.round(f * 100)) }), f); await frame();
    });
    busy(t('aligning'), null); await frame();
    const m = new core.CloudModel(pts, nrm, n);
    const thr = Math.max(2 * m.spacing, 1e-4);
    m.align(thr);
    state.model = m; state.fname = file.name;
    $('thr').value = +thr.toPrecision(4);
    $('tol').value = +Math.max(3 * m.rms, 0.5 * m.spacing, 1e-4).toPrecision(3);
    busy(t('building'), null); await frame();
    afterAlign(true);
    for (const id of ['grpPlane', 'grpProbe', 'grpMeasure']) $(id).disabled = false;
  } catch (e) {
    console.error(e);
    alert(t('eLoad', { 0: e.message === 'noNumeric' ? t('eNoNum') : e.message }));
  } finally { busy(null); }
}

function afterAlign(resetProbe) {
  const m = state.model;
  let lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i < m.n; i++) {
    if (m.x[i] < lo[0]) lo[0] = m.x[i]; if (m.x[i] > hi[0]) hi[0] = m.x[i];
    if (m.y[i] < lo[1]) lo[1] = m.y[i]; if (m.y[i] > hi[1]) hi[1] = m.y[i];
    if (m.z[i] < lo[2]) lo[2] = m.z[i]; if (m.z[i] > hi[2]) hi[2] = m.z[i];
  }
  state.bounds = { lo, hi };
  const ext = Math.max(hi[0] - lo[0], hi[1] - lo[1]);
  setRange('x', lo[0], hi[0]); setRange('y', lo[1], hi[1]);
  setRange('z', Math.min(lo[2], 0) - ext * 0.1, hi[2] + ext * 0.2);
  setRange('r', Math.max(m.spacing * 2, 1e-4), ext / 2);
  if (resetProbe) {
    setVal('r', Math.min(Math.max(hi[2], ext * 0.03), ext / 2));
    setVal('x', (lo[0] + hi[0]) / 2); setVal('y', (lo[1] + hi[1]) / 2); setVal('z', val('r'));
  }
  updateInfo();
  buildCloud(resetProbe);
  probeChanged();
}

function updateInfo() {
  const m = state.model;
  if (!m) { $('info').innerHTML = t('noFile'); return; }
  const b = state.bounds;
  $('info').innerHTML = t('info', {
    f: state.fname.replace(/</g, '&lt;'), n: m.n.toLocaleString(), s: m.spacing.toFixed(4),
    inl: m.nInliers.toLocaleString(), rms: m.rms.toFixed(4), nrm: m.N ? t('yes') : t('no'),
    zmin: b.lo[2].toFixed(3), zmax: b.hi[2].toFixed(3),
  });
}

async function realign(toggleFlip) {
  const m = state.model;
  if (!m) return;
  busy(t('aligning'), null); await frame();
  try {
    m.align(parseFloat($('thr').value) || m.thr, toggleFlip ? !m.flip : m.flip);
    $('tol').value = +Math.max(3 * m.rms, 0.5 * m.spacing, 1e-4).toPrecision(3);
    afterAlign(false);
  } finally { busy(null); }
}

// ---------------------------------------------------------------- measurement
async function doMeasure() {
  const m = state.model;
  if (!m) { alert(t('eFirst')); return; }
  busy(t('measuring'), null); await frame();
  let res;
  try {
    res = m.measure(probeCenter(), val('r'), parseFloat($('tol').value), { cylinder: $('cbCyl').checked, keepLargest: $('cbKeep').checked, voxel: parseFloat($('voxSize').value) || 0 });
  } catch (e) {
    busy(null);
    const [code, arg] = e.message.split('|');
    alert(code === 'fewPoints' ? t('eFew', { 0: arg }) : code === 'noObject' ? t('eNoObj', { 0: arg }) : e.message);
    return;
  }
  busy(null);
  res.file = state.fname;
  res.num = ++measureCounter;
  state.measures.push(res);
  renderTable();
  showResult(state.measures.length - 1);
}

function renderTable() {
  const tb = $('table').tBodies[0];
  tb.innerHTML = '';
  state.measures.forEach((r, i) => {
    const tr = document.createElement('tr');
    if (i === state.selRow) tr.className = 'sel';
    for (const v of [r.num ?? i + 1, r.center[0].toFixed(3), r.center[1].toFixed(3), f4(r.Vbest), f4(r.Vscan), f4(r.Vsphere), f4(r.Vell), f4(r.vox ? r.vox.V : NaN), f4(r.hmax)]) {
      const td = document.createElement('td'); td.textContent = v; tr.appendChild(td);
    }
    tr.addEventListener('click', () => { state.selRow = i; renderTable(); });
    tr.addEventListener('dblclick', () => showResult(i));
    tb.appendChild(tr);
  });
}
let measureCounter = 0;

// ---------------------------------------------------------------- result window
let rviewer;
const rs = {};
function showResult(i) {
  const res = state.measures[i];
  if (res.num == null) res.num = ++measureCounter;
  state.cur = res; state.curIdx = i;
  $('result').hidden = false;
  if (!rviewer) rviewer = new Viewer($('rview'));
  else rviewer.resize();
  renderResultText();
  buildResultScene(true);
}

function buildResultScene(resetCam) {
  const res = state.cur;
  if (!res) return;
  if (rs.group) {
    rviewer.scene.remove(rs.group);
    rs.group.traverse(o => { if (o.geometry) o.geometry.dispose(); });
  }
  const G = new THREE.Group();
  const Q = res.Q, nq = res.nq;
  const pos = Float32Array.from(Q);
  const col = new Float32Array(3 * nq);
  const hm = Math.max(res.hmax, 1e-9);
  for (let i = 0; i < nq; i++) viridis(Q[3 * i + 2] / hm, col, 3 * i);
  $('rlgMax').textContent = res.hmax.toFixed(3);
  if ($('cbMesh').checked) {
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.BufferAttribute(pos, 3));
    g.setAttribute('color', new THREE.BufferAttribute(col, 3));
    g.setIndex(new THREE.BufferAttribute(res.tri, 1));
    g.computeVertexNormals();
    G.add(new THREE.Mesh(g, new THREE.MeshPhongMaterial({ vertexColors: true, flatShading: true, side: THREE.DoubleSide, shininess: 15 })));
  }
  const og = new THREE.BufferGeometry();
  og.setAttribute('position', new THREE.BufferAttribute(pos.slice(0, 3 * res.nObj), 3));
  og.setAttribute('color', new THREE.BufferAttribute(col.slice(0, 3 * res.nObj), 3));
  G.add(new THREE.Points(og, new THREE.PointsMaterial({ size: $('cbMesh').checked ? 2 : 3, sizeAttenuation: false, vertexColors: true })));
  if ($('cbRemoved').checked && res.removed.length) {
    const rg = new THREE.BufferGeometry();
    rg.setAttribute('position', new THREE.BufferAttribute(Float32Array.from(res.removed), 3));
    G.add(new THREE.Points(rg, new THREE.PointsMaterial({ size: 2, sizeAttenuation: false, color: 0xa0a0a0 })));
  }
  if ($('cbFit').checked) {
    const sg = new THREE.SphereGeometry(res.sphere.r, 40, 20);
    sg.translate(...res.sphere.c);
    G.add(new THREE.LineSegments(new THREE.WireframeGeometry(sg), new THREE.LineBasicMaterial({ color: 0xff8800, transparent: true, opacity: 0.45 })));
    if (res.ell) {
      const eg = new THREE.SphereGeometry(1, 40, 20);
      const M = res.ell.M, c = res.ell.c;
      eg.applyMatrix4(new THREE.Matrix4().set(M[0][0], M[0][1], M[0][2], c[0], M[1][0], M[1][1], M[1][2], c[1], M[2][0], M[2][1], M[2][2], c[2], 0, 0, 0, 1));
      G.add(new THREE.LineSegments(new THREE.WireframeGeometry(eg), new THREE.LineBasicMaterial({ color: 0xcc00cc, transparent: true, opacity: 0.35 })));
    }
  }
  $('voxShown').textContent = '';
  if ($('cbVox').checked && res.vox) {
    const vx = res.vox, occ = vx.occ, nAll = occ.length / 4, LIM = 200000;
    let sel = [];
    for (let q = 0; q < nAll; q++) if (occ[4 * q + 3] >= 0.5) sel.push(q);
    let partial = false;
    if (sel.length > LIM) { sel = sel.filter(q => occ[4 * q + 3] < 1); partial = true; }
    if (sel.length > LIM) { const st = Math.ceil(sel.length / LIM); sel = sel.filter((_, k) => k % st === 0); }
    const s = vx.size;
    const im = new THREE.InstancedMesh(new THREE.BoxGeometry(s * 0.9, s * 0.9, s * 0.9),
      new THREE.MeshLambertMaterial({ color: 0xffffff }), sel.length);
    const mt = new THREE.Matrix4(), cc = new Float32Array(3), colr = new THREE.Color();
    sel.forEach((q, k) => {
      const zc = (occ[4 * q + 2] + 0.5) * s;
      mt.makeTranslation(vx.origin[0] + (occ[4 * q] + 0.5) * s, vx.origin[1] + (occ[4 * q + 1] + 0.5) * s, zc);
      im.setMatrixAt(k, mt);
      viridis(zc / hm, cc, 0); im.setColorAt(k, colr.setRGB(cc[0], cc[1], cc[2]));
    });
    G.add(im);
    $('voxShown').textContent = partial ? t('voxShownPart', { n: sel.length.toLocaleString() }) : t('voxShown', { n: sel.length.toLocaleString() });
  }
  // grid on the plane + axes
  const box = new THREE.Box3().setFromBufferAttribute(new THREE.BufferAttribute(pos, 3));
  const size = box.getSize(new THREE.Vector3()), ctr = box.getCenter(new THREE.Vector3());
  const gs = Math.max(size.x, size.y) * 1.05;
  const grid = new THREE.GridHelper(gs, 10, 0x999999, 0xcccccc);
  grid.rotation.x = Math.PI / 2; grid.position.set(ctr.x, ctr.y, -1e-4 * gs);
  G.add(grid);
  const ax = new THREE.AxesHelper(gs * 0.15); ax.position.set(box.min.x, box.min.y, 0); G.add(ax);
  rs.group = G;
  rviewer.scene.add(G);
  if (resetCam) rviewer.fit(box); else rviewer.render();
}

function renderResultText() {
  const r = state.cur;
  $('resTitle').textContent = t('resTitle', { k: r.num });
  const why = r.undercut ? t('whyUnder', { fit: t(r.fitName) }) : t('whyScan');
  const ell = r.ell ? `${f4(r.Vell)} mm³</td></tr><tr><td colspan="2" class="sub small">${t('ellDet', { a: r.ell.axes[0].toFixed(3), b: r.ell.axes[1].toFixed(3), c: r.ell.axes[2].toFixed(3), rms: r.ell.rms.toFixed(4) })}` : t('fitFail');
  const vx = r.vox;
  const vox = vx ? `<tr><td>${t('vVox', { s: vx.size.toFixed(4) })}</td><td>${f4(vx.V)} mm³</td></tr><tr><td colspan="2" class="sub small">${t('voxDet', { f: vx.nFull.toLocaleString(), p: vx.nPart.toLocaleString(), m: 4, a: vx.Vmin.toFixed(3), b: vx.Vmax.toFixed(3) })}</td></tr>` : '';
  const w = r.warnings.length ? r.warnings.map(([k, a]) => `<li class="w">${t(k, { 0: a })}</li>`).join('') : `<li>${t('none')}</li>`;
  $('resText').innerHTML = `
    <h3>${t('volTitle')}</h3>
    <div>${t('best')}:</div>
    <div class="big">${f4(r.Vbest)} mm³</div>
    <div class="small">${why}</div>
    <table>
      <tr><td>${t('vScan')}</td><td>${f4(r.Vscan)} mm³</td></tr>
      <tr><td>${t('vCol')}</td><td>${f4(r.Vcol)} mm³</td></tr>
      <tr><td>${t('vSph')}</td><td>${f4(r.Vsphere)} mm³</td></tr>
      <tr><td>${t('vEll')}</td><td>${ell}</td></tr>
      ${vox}
    </table>
    <p>${t('sphDet', { r: r.sphere.r.toFixed(4), h: r.sphere.c[2].toFixed(4), rms: r.sphere.rms.toFixed(4) })}<br>
    ${t('hmax')}: ${r.hmax.toFixed(4)} mm<br>
    ${t('foot')}: ${r.footprint.toFixed(4)} mm²<br>
    ${t('nPts', { o: r.nObj, b: r.nBase })}</p>
    <p>${t('local', { t: r.localTilt.toFixed(4), rms: f4(r.localRms), n: r.nRing })}<br>
    ${t('scanDir', { a: r.scanAngle.toFixed(2) })}<br>
    ${t('tolUsed', { t: r.tol.toFixed(4) })}</p>
    <p><b>${t('warnings')}</b></p><ul>${w}</ul>
    <p class="small">${t('note')}</p>
    <p class="small">${t('voxNote')}</p>`;
  if (vx) $('voxSize2').value = +vx.size.toPrecision(3);
}

// ---------------------------------------------------------------- exports
function download(data, name, type) {
  const blob = data instanceof Blob ? data : new Blob(Array.isArray(data) ? data : [data], { type });
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob); a.download = name;
  document.body.appendChild(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(a.href), 4000);
}
const baseName = () => (state.fname || 'cloud').replace(/\.[^.]+$/, '');

function exportCsv() {
  if (!state.measures.length) return;
  const head = ['file', 'n', 'X_aligned', 'Y_aligned', 'Z_probe', 'R_probe', 'cylinder', 'tolerance_mm',
    'V_recommended_mm3', 'criterion', 'V_surface_mm3', 'V_column_mm3', 'V_sphere_mm3', 'r_sphere_mm', 'rms_sphere_mm',
    'V_ellipsoid_mm3', 'rms_ellipsoid_mm', 'V_voxel_mm3', 'V_voxel_min_mm3', 'V_voxel_max_mm3', 'voxel_edge_mm', 'h_max_mm', 'footprint_mm2', 'object_points', 'scan_angle_deg', 'warnings'];
  const q = s => '"' + String(s).replace(/"/g, '""') + '"';
  const rows = [head.join(';')];
  state.measures.forEach((r, i) => {
    rows.push([q(r.file), r.num ?? i + 1, r.center[0].toFixed(5), r.center[1].toFixed(5), r.center[2].toFixed(5), r.r.toFixed(5),
      r.cylinder ? 1 : 0, r.tol.toFixed(5), r.Vbest.toFixed(6), r.undercut ? r.fitName + '_fit' : 'closed_surface',
      r.Vscan.toFixed(6), r.Vcol.toFixed(6), r.Vsphere.toFixed(6), r.sphere.r.toFixed(5), r.sphere.rms.toFixed(5),
      isFinite(r.Vell) ? r.Vell.toFixed(6) : '', r.ell ? r.ell.rms.toFixed(5) : '',
      r.vox ? r.vox.V.toFixed(6) : '', r.vox ? r.vox.Vmin.toFixed(6) : '', r.vox ? r.vox.Vmax.toFixed(6) : '', r.vox ? r.vox.size.toFixed(5) : '',
      r.hmax.toFixed(5), r.footprint.toFixed(5),
      r.nObj, r.scanAngle.toFixed(3), q(r.warnings.map(([k, a]) => t(k, { 0: a })).join(' | '))].join(';'));
  });
  download('﻿' + rows.join('\r\n') + '\r\n', baseName() + '_volumes.csv', 'text/csv');
}

function deleteRow() {
  if (state.selRow < 0 || state.selRow >= state.measures.length) return;
  state.measures.splice(state.selRow, 1);
  state.selRow = -1;
  renderTable();
}

async function exportAligned() {
  const m = state.model;
  if (!m) return;
  busy(t('building'), null); await frame();
  const parts = [m.N ? 'Orig_Pos[X] Orig_Pos[Y] Orig_Pos[Z] Orig_Vec[I] Orig_Vec[J] Orig_Vec[K] \n' + `SCANLINE 0  ( n = ${m.n} ) \n` : 'Orig_Pos[X] Orig_Pos[Y] Orig_Pos[Z] \n'];
  for (let s = 0; s < m.n; s += 100000) {
    const lines = [];
    for (let i = s; i < Math.min(m.n, s + 100000); i++) {
      let l = `${m.x[i].toFixed(8)} ${m.y[i].toFixed(8)} ${m.z[i].toFixed(8)}`;
      if (m.N) l += ` ${m.N[3 * i].toFixed(8)} ${m.N[3 * i + 1].toFixed(8)} ${m.N[3 * i + 2].toFixed(8)}`;
      lines.push(l);
    }
    parts.push(lines.join('\n') + '\n');
  }
  busy(null);
  download(parts, baseName() + '_aligned.asc', 'text/plain');
}

// ---------------------------------------------------------------- wiring
function wire() {
  $('lang').addEventListener('change', e => { setLang(e.target.value); applyLang(); });
  $('btnOpen').addEventListener('click', () => $('file').click());
  $('file').addEventListener('change', e => { const f = e.target.files[0]; if (f) loadFile(f); e.target.value = ''; });
  $('btnRefit').addEventListener('click', () => realign(false));
  $('btnFlip').addEventListener('click', () => realign(true));
  $('btnMeasure').addEventListener('click', doMeasure);
  $('btnCsv').addEventListener('click', exportCsv);
  $('btnDel').addEventListener('click', deleteRow);
  $('btnAligned').addEventListener('click', exportAligned);
  $('btnClose').addEventListener('click', () => { $('result').hidden = true; });
  $('btnVox').addEventListener('click', async () => {
    const r = state.cur; if (!r) return;
    const s = parseFloat($('voxSize2').value);
    if (!(s > 0)) return;
    busy(t('measuring'), null); await frame();
    try {
      r.vox = core.voxelVolume(r.Q, r.nq, r.tri, r.scanDir, s);
      r.warnings = r.warnings.filter(w => w[0] !== 'wVoxGrow');
      if (r.vox && r.vox.size > s * 1.01) r.warnings.push(['wVoxGrow', r.vox.size.toFixed(4)]);
    } finally { busy(null); }
    renderResultText(); renderTable(); buildResultScene(false);
  });
  ['cbMesh', 'cbRemoved', 'cbFit', 'cbVox'].forEach(id => $(id).addEventListener('change', () => buildResultScene(false)));
  $('btnExpPts').addEventListener('click', () => {
    const r = state.cur; if (!r) return;
    download(core.ascText(r.objOrig, r.objN, r.nObj), `${baseName()}_object${r.num}.asc`, 'text/plain');
  });
  $('btnExpStl').addEventListener('click', () => {
    const r = state.cur; if (!r) return;
    download(core.stlBinary(r.Q, core.closeMeshFaces(r.Q, r.tri)), `${baseName()}_object${r.num}.stl`, 'model/stl');
  });
  $('btnCopy').addEventListener('click', async () => {
    const txt = $('resText').innerText;
    try { await navigator.clipboard.writeText(txt); } catch (e) {
      const ta = document.createElement('textarea'); ta.value = txt; document.body.appendChild(ta); ta.select();
      try { document.execCommand('copy'); } catch (e2) { /* ignore */ }
      ta.remove();
    }
    const b = $('btnCopy'); b.textContent = t('copied'); setTimeout(() => { b.textContent = t('copy'); }, 1200);
  });

  const cv = viewer.renderer.domElement;
  cv.addEventListener('pointermove', e => { mouse.x = e.clientX; mouse.y = e.clientY; mouse.inside = true; });
  cv.addEventListener('pointerleave', () => { mouse.inside = false; });
  cv.addEventListener('dblclick', e => pickAt(e.clientX, e.clientY));
  window.addEventListener('keydown', e => {
    const tag = (e.target.tagName || '').toUpperCase();
    if (['INPUT', 'SELECT', 'TEXTAREA'].includes(tag) && e.target.type !== 'range' && e.target.type !== 'checkbox') return;
    if (!$('result').hidden) { if (e.key === 'Escape') $('result').hidden = true; return; }
    const map = {
      ArrowLeft: [-1, 0, 0, 0], ArrowRight: [1, 0, 0, 0], ArrowUp: [0, 1, 0, 0], ArrowDown: [0, -1, 0, 0],
      PageUp: [0, 0, 1, 0], PageDown: [0, 0, -1, 0], ']': [0, 0, 0, 1], '[': [0, 0, 0, -1],
    };
    if (map[e.key]) { e.preventDefault(); nudge(...map[e.key]); }
    else if ((e.key === 'p' || e.key === 'P') && mouse.inside) pickAt(mouse.x, mouse.y);
  });

  // drag & drop
  let dragDepth = 0;
  window.addEventListener('dragenter', e => { e.preventDefault(); dragDepth++; $('dropzone').hidden = false; });
  window.addEventListener('dragleave', () => { if (--dragDepth <= 0) { dragDepth = 0; $('dropzone').hidden = true; } });
  window.addEventListener('dragover', e => e.preventDefault());
  window.addEventListener('drop', e => {
    e.preventDefault(); dragDepth = 0; $('dropzone').hidden = true;
    const f = e.dataTransfer.files[0]; if (f) loadFile(f);
  });
}

try {
  initMainView();
  wire();
} catch (e) {
  console.error(e);
  alert(t('eWebgl'));
}
applyLang();
window.__mv = { state, loadFile, doMeasure, setVal, probeChanged, showResult, pickAt, viewer: () => viewer, rviewer: () => rviewer };
