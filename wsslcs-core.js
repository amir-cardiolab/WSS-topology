// wsslcs-core.js — fixed points and stable/unstable manifolds of a vector
// field on a triangulated surface.  Browser port of the Python package
// `wsslcs` (mesh.py, field.py, fixedpoints.py, manifolds.py, tracer.py).
// Pure JavaScript, no DOM: runs in a Web Worker or any JS engine.

export const TYPE_NAMES = { 0: 'unknown', 1: 'source', 2: 'sink', 3: 'saddle', 4: 'center', 6: 'attracting focus', 7: 'repelling focus' };
export const UNKNOWN = 0, SOURCE = 1, SINK = 2, SADDLE = 3, CENTER = 4, ATTRACTING_FOCUS = 6, REPELLING_FOCUS = 7;
export const UNSTABLE = 0, STABLE = 1;
const ACTIVE = 0, LEFT = 1, STUCK = 2;
const TWO_PI = 2 * Math.PI;

const hypot2 = (x, y) => Math.sqrt(x * x + y * y);
const wrapPi = a => { a = (a + Math.PI) % TWO_PI; if (a < 0) a += TWO_PI; return a - Math.PI; };

// ---------------------------------------------------------------------------
// mesh
// ---------------------------------------------------------------------------
export function buildMesh(pointsIn, trianglesIn) {
  const P = pointsIn instanceof Float64Array ? pointsIn : Float64Array.from(pointsIn);
  const N = P.length / 3;
  const triAll = trianglesIn instanceof Int32Array ? trianglesIn : Int32Array.from(trianglesIn);
  // drop degenerate triangles
  const keep = [];
  for (let t = 0; t < triAll.length / 3; t++) {
    const a = triAll[3 * t], b = triAll[3 * t + 1], c = triAll[3 * t + 2];
    if (a === b || b === c || a === c || a < 0 || b < 0 || c < 0 || a >= N || b >= N || c >= N) continue;
    const e1x = P[3 * b] - P[3 * a], e1y = P[3 * b + 1] - P[3 * a + 1], e1z = P[3 * b + 2] - P[3 * a + 2];
    const e2x = P[3 * c] - P[3 * a], e2y = P[3 * c + 1] - P[3 * a + 1], e2z = P[3 * c + 2] - P[3 * a + 2];
    const nx = e1y * e2z - e1z * e2y, ny = e1z * e2x - e1x * e2z, nz = e1x * e2y - e1y * e2x;
    if (nx * nx + ny * ny + nz * nz > 0) keep.push(t);
  }
  const M = keep.length;
  if (M === 0) throw new Error('the mesh has no valid triangles');
  const tri = new Int32Array(3 * M);
  for (let i = 0; i < M; i++) { tri[3 * i] = triAll[3 * keep[i]]; tri[3 * i + 1] = triAll[3 * keep[i] + 1]; tri[3 * i + 2] = triAll[3 * keep[i] + 2]; }
  const mesh = { points: P, triangles: tri, nPoints: N, nTris: M, warnings: [] };
  if (M !== triAll.length / 3) mesh.warnings.push(`${triAll.length / 3 - M} degenerate triangles removed`);
  buildEdges(mesh);
  buildGeometry(mesh);
  buildFrames(mesh);
  buildFans(mesh);
  return mesh;
}

function buildEdges(mesh) {
  const { triangles: tri, nPoints: N, nTris: M } = mesh;
  const map = new Map();
  const edgeV = [], edgeTris = [];
  const triEdges = new Int32Array(3 * M);
  let nonManifold = 0;
  for (let t = 0; t < M; t++) {
    for (let j = 0; j < 3; j++) {
      const a = tri[3 * t + j], b = tri[3 * t + (j + 1) % 3];
      const key = (a < b ? a : b) * N + (a < b ? b : a);
      let e = map.get(key);
      if (e === undefined) {
        e = edgeV.length / 2;
        map.set(key, e);
        edgeV.push(a < b ? a : b, a < b ? b : a);
        edgeTris.push(t, -1);
      } else if (edgeTris[2 * e + 1] < 0) {
        edgeTris[2 * e + 1] = t;
      } else {
        nonManifold++;
      }
      triEdges[3 * t + j] = e;
    }
  }
  if (nonManifold) mesh.warnings.push(`${nonManifold} non-manifold edge uses (more than two triangles)`);
  const E = edgeV.length / 2;
  const triNeighbors = new Int32Array(3 * M);
  const boundaryTriangle = new Uint8Array(M);
  const boundaryVertex = new Uint8Array(N);
  const boundaryEdge = new Uint8Array(E);
  let nBoundaryEdges = 0;
  for (let e = 0; e < E; e++) {
    if (edgeTris[2 * e + 1] < 0) { boundaryEdge[e] = 1; nBoundaryEdges++; boundaryVertex[edgeV[2 * e]] = 1; boundaryVertex[edgeV[2 * e + 1]] = 1; }
  }
  for (let t = 0; t < M; t++) {
    for (let j = 0; j < 3; j++) {
      const e = triEdges[3 * t + j];
      const t0 = edgeTris[2 * e], t1 = edgeTris[2 * e + 1];
      const nb = t0 === t ? t1 : t0;
      triNeighbors[3 * t + j] = nb;
      if (nb < 0) boundaryTriangle[t] = 1;
    }
  }
  let sumLen = 0;
  const P = mesh.points;
  for (let e = 0; e < E; e++) {
    const a = edgeV[2 * e], b = edgeV[2 * e + 1];
    sumLen += Math.sqrt((P[3 * a] - P[3 * b]) ** 2 + (P[3 * a + 1] - P[3 * b + 1]) ** 2 + (P[3 * a + 2] - P[3 * b + 2]) ** 2);
  }
  Object.assign(mesh, { edges: Int32Array.from(edgeV), edgeTris: Int32Array.from(edgeTris), nEdges: E, triEdges, triNeighbors,
    boundaryTriangle, boundaryVertex, boundaryEdge, nBoundaryEdges, meanEdgeLength: sumLen / E });
}

function buildGeometry(mesh) {
  const { points: P, triangles: tri, nPoints: N, nTris: M } = mesh;
  const normals = new Float64Array(3 * M), area = new Float64Array(M);
  const mn = [Infinity, Infinity, Infinity], mx = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i < N; i++) for (let k = 0; k < 3; k++) { const v = P[3 * i + k]; if (v < mn[k]) mn[k] = v; if (v > mx[k]) mx[k] = v; }
  const center = [(mn[0] + mx[0]) / 2, (mn[1] + mx[1]) / 2, (mn[2] + mx[2]) / 2];
  const radius = Math.sqrt((center[0] - mn[0]) ** 2 + (center[1] - mn[1]) ** 2 + (center[2] - mn[2]) ** 2);
  let sv = 0, total = 0;
  for (let t = 0; t < M; t++) {
    const a = tri[3 * t], b = tri[3 * t + 1], c = tri[3 * t + 2];
    const e1x = P[3 * b] - P[3 * a], e1y = P[3 * b + 1] - P[3 * a + 1], e1z = P[3 * b + 2] - P[3 * a + 2];
    const e2x = P[3 * c] - P[3 * a], e2y = P[3 * c + 1] - P[3 * a + 1], e2z = P[3 * c + 2] - P[3 * a + 2];
    let nx = e1y * e2z - e1z * e2y, ny = e1z * e2x - e1x * e2z, nz = e1x * e2y - e1y * e2x;
    const l = Math.sqrt(nx * nx + ny * ny + nz * nz);
    area[t] = l / 2; total += area[t];
    nx /= l; ny /= l; nz /= l;
    normals[3 * t] = nx; normals[3 * t + 1] = ny; normals[3 * t + 2] = nz;
    sv += ((center[0] - P[3 * a]) * nx + (center[1] - P[3 * a + 1]) * ny + (center[2] - P[3 * a + 2]) * nz) * area[t];
  }
  const flipped = sv > 0;
  if (flipped) for (let i = 0; i < normals.length; i++) normals[i] = -normals[i];
  const vn = new Float64Array(3 * N);
  for (let t = 0; t < M; t++) for (let j = 0; j < 3; j++) {
    const v = tri[3 * t + j];
    vn[3 * v] += normals[3 * t] * area[t]; vn[3 * v + 1] += normals[3 * t + 1] * area[t]; vn[3 * v + 2] += normals[3 * t + 2] * area[t];
  }
  for (let i = 0; i < N; i++) {
    const l = Math.sqrt(vn[3 * i] ** 2 + vn[3 * i + 1] ** 2 + vn[3 * i + 2] ** 2);
    if (l > 0) { vn[3 * i] /= l; vn[3 * i + 1] /= l; vn[3 * i + 2] /= l; }
  }
  Object.assign(mesh, { triNormals: normals, area, totalArea: total, meanArea: total / M, vertexNormals: vn,
    center, radius, bounds: [mn, mx], normalsFlipped: flipped });
}

function buildFrames(mesh) {
  const { points: P, triangles: tri, nTris: M, triNormals: nrm } = mesh;
  const LX = new Float64Array(3 * M), LY = new Float64Array(3 * M), x1 = new Float64Array(M), x2 = new Float64Array(M), y2 = new Float64Array(M);
  for (let t = 0; t < M; t++) {
    const a = tri[3 * t], b = tri[3 * t + 1], c = tri[3 * t + 2];
    const e1x = P[3 * b] - P[3 * a], e1y = P[3 * b + 1] - P[3 * a + 1], e1z = P[3 * b + 2] - P[3 * a + 2];
    const e2x = P[3 * c] - P[3 * a], e2y = P[3 * c + 1] - P[3 * a + 1], e2z = P[3 * c + 2] - P[3 * a + 2];
    const l1 = Math.sqrt(e1x * e1x + e1y * e1y + e1z * e1z);
    const lx = [e1x / l1, e1y / l1, e1z / l1];
    const nx = nrm[3 * t], ny = nrm[3 * t + 1], nz = nrm[3 * t + 2];
    let lyx = ny * lx[2] - nz * lx[1], lyy = nz * lx[0] - nx * lx[2], lyz = nx * lx[1] - ny * lx[0];
    const ll = Math.sqrt(lyx * lyx + lyy * lyy + lyz * lyz);
    lyx /= ll; lyy /= ll; lyz /= ll;
    LX[3 * t] = lx[0]; LX[3 * t + 1] = lx[1]; LX[3 * t + 2] = lx[2];
    LY[3 * t] = lyx; LY[3 * t + 1] = lyy; LY[3 * t + 2] = lyz;
    x1[t] = l1;
    x2[t] = e2x * lx[0] + e2y * lx[1] + e2z * lx[2];
    y2[t] = e2x * lyx + e2y * lyy + e2z * lyz;
  }
  Object.assign(mesh, { LX, LY, x1, x2, y2 });
}

function buildFans(mesh) {
  const { points: P, triangles: tri, nPoints: N, nTris: M, triNeighbors: nb, vertexNormals: VN } = mesh;
  // CSR of incident corners
  const counts = new Int32Array(N + 1);
  for (let i = 0; i < 3 * M; i++) counts[tri[i] + 1]++;
  for (let i = 0; i < N; i++) counts[i + 1] += counts[i];
  const vtxStart = counts;
  const fill = new Int32Array(N);
  const vtxTri = new Int32Array(3 * M), vtxCorner = new Int32Array(3 * M);
  for (let t = 0; t < M; t++) for (let j = 0; j < 3; j++) {
    const v = tri[3 * t + j];
    const k = vtxStart[v] + fill[v]++;
    vtxTri[k] = t; vtxCorner[k] = j;
  }
  // corner angles
  const ang = new Float64Array(3 * M);
  for (let t = 0; t < M; t++) for (let j = 0; j < 3; j++) {
    const v = tri[3 * t + j], a = tri[3 * t + (j + 1) % 3], b = tri[3 * t + (j + 2) % 3];
    const ax = P[3 * a] - P[3 * v], ay = P[3 * a + 1] - P[3 * v + 1], az = P[3 * a + 2] - P[3 * v + 2];
    const bx = P[3 * b] - P[3 * v], by = P[3 * b + 1] - P[3 * v + 1], bz = P[3 * b + 2] - P[3 * v + 2];
    const den = Math.max(Math.sqrt(ax * ax + ay * ay + az * az) * Math.sqrt(bx * bx + by * by + bz * bz), 1e-300);
    ang[3 * t + j] = Math.acos(Math.min(1, Math.max(-1, (ax * bx + ay * by + az * bz) / den)));
  }
  const cornerBegin = new Float64Array(3 * M);
  const cornerOther = new Int32Array(3 * M);
  for (let t = 0; t < M; t++) for (let j = 0; j < 3; j++) cornerOther[3 * t + j] = (j + 2) % 3;
  const inFan = new Uint8Array(3 * M);
  const vT = new Float64Array(3 * N), vB = new Float64Array(3 * N), vR = new Float64Array(N).fill(1), vOrient = new Float64Array(N).fill(1);
  let nonManifoldVerts = 0;
  const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const norm = a => Math.sqrt(dot(a, a));
  for (let v = 0; v < N; v++) {
    const s0 = vtxStart[v], s1 = vtxStart[v + 1], K = s1 - s0;
    const n = [VN[3 * v], VN[3 * v + 1], VN[3 * v + 2]];
    if (K === 0) {
      let T = cross(n, [1, 0, 0]); if (norm(T) < 1e-8) T = cross(n, [0, 1, 0]);
      const l = norm(T) || 1; T = [T[0] / l, T[1] / l, T[2] / l];
      const B = cross(n, T);
      vT.set(T, 3 * v); vB.set(B, 3 * v);
      continue;
    }
    let start = -1, eIn = -1;
    for (let k = 0; k < K; k++) { const t = vtxTri[s0 + k], j = vtxCorner[s0 + k]; if (nb[3 * t + (j + 2) % 3] < 0) { start = k; eIn = (j + 2) % 3; break; } }
    if (start < 0) for (let k = 0; k < K; k++) { const t = vtxTri[s0 + k], j = vtxCorner[s0 + k]; if (nb[3 * t + j] < 0) { start = k; eIn = j; break; } }
    if (start < 0) { start = 0; eIn = (vtxCorner[s0] + 2) % 3; }
    let t = vtxTri[s0 + start], j = vtxCorner[s0 + start];
    const fanT = [], fanJ = [], fanOther = [];
    const visited = new Set();
    let closed = false;
    for (;;) {
      if (visited.has(t)) break;
      visited.add(t);
      const otherIn = eIn === j ? (j + 1) % 3 : (j + 2) % 3;
      fanT.push(t); fanJ.push(j); fanOther.push(otherIn);
      const eOut = eIn !== j ? j : (j + 2) % 3;
      const tn = nb[3 * t + eOut];
      if (tn < 0) break;
      if (tn === fanT[0]) { closed = true; break; }
      const w = eOut === j ? tri[3 * t + (j + 1) % 3] : tri[3 * t + (j + 2) % 3];
      let jn = -1, wn = -1;
      for (let k = 0; k < 3; k++) { if (tri[3 * tn + k] === v) jn = k; if (tri[3 * tn + k] === w) wn = k; }
      if (jn < 0 || wn < 0) break;
      eIn = wn === (jn + 1) % 3 ? jn : (jn + 2) % 3;
      t = tn; j = jn;
    }
    if (fanT.length < K) nonManifoldVerts++;
    const Kf = fanT.length;
    const angles = fanT.map((tt, k) => ang[3 * tt + fanJ[k]]);
    const total = angles.reduce((a, b) => a + b, 0);
    const r = closed && total > 0 ? TWO_PI / total : 1.0;
    let T = null, refK = 0;
    for (let k = 0; k < Kf; k++) {
      const o = tri[3 * fanT[k] + fanOther[k]];
      const d0 = [P[3 * o] - P[3 * v], P[3 * o + 1] - P[3 * v + 1], P[3 * o + 2] - P[3 * v + 2]];
      const dn = dot(d0, n);
      const tt = [d0[0] - dn * n[0], d0[1] - dn * n[1], d0[2] - dn * n[2]];
      const nt = norm(tt);
      if (nt > 1e-12 * Math.max(norm(d0), 1e-300)) { T = [tt[0] / nt, tt[1] / nt, tt[2] / nt]; refK = k; break; }
    }
    if (T === null) {
      T = cross(n, [1, 0, 0]); if (norm(T) < 1e-8) T = cross(n, [0, 1, 0]);
      const l = norm(T) || 1; T = [T[0] / l, T[1] / l, T[2] / l];
    }
    let B = cross(n, T);
    const nB = norm(B); if (nB > 0) B = [B[0] / nB, B[1] / nB, B[2] / nB];
    const t0 = fanT[0], j0 = fanJ[0], o0 = fanOther[0];
    const third = 3 - j0 - o0;
    const tv = tri[3 * t0 + third];
    const d1 = [P[3 * tv] - P[3 * v], P[3 * tv + 1] - P[3 * v + 1], P[3 * tv + 2] - P[3 * v + 2]];
    const phi1 = Math.atan2(dot(d1, B), dot(d1, T));
    const s = phi1 >= 0 ? 1 : -1;
    const cum = new Float64Array(Kf + 1);
    for (let k = 0; k < Kf; k++) cum[k + 1] = cum[k] + angles[k] * r;
    for (let k = 0; k < Kf; k++) {
      const idx = 3 * fanT[k] + fanJ[k];
      cornerBegin[idx] = s * (cum[k] - cum[refK]);
      cornerOther[idx] = fanOther[k];
      inFan[idx] = 1;
    }
    vT.set(T, 3 * v); vB.set(B, 3 * v); vR[v] = r; vOrient[v] = s;
    if (Kf < K) {
      for (let k = 0; k < K; k++) {
        const tt = vtxTri[s0 + k], jj = vtxCorner[s0 + k];
        if (inFan[3 * tt + jj]) continue;
        const o = (jj + 2) % 3, ov = tri[3 * tt + o];
        const d = [P[3 * ov] - P[3 * v], P[3 * ov + 1] - P[3 * v + 1], P[3 * ov + 2] - P[3 * v + 2]];
        cornerBegin[3 * tt + jj] = Math.atan2(dot(d, B), dot(d, T));
        cornerOther[3 * tt + jj] = o;
      }
    }
  }
  if (nonManifoldVerts) mesh.warnings.push(`${nonManifoldVerts} non-manifold vertices`);
  // unit local direction of the reference edge of every corner
  const E0 = new Float64Array(6 * M);
  for (let t = 0; t < M; t++) for (let j = 0; j < 3; j++) {
    const v = tri[3 * t + j], o = tri[3 * t + cornerOther[3 * t + j]];
    const dx = P[3 * o] - P[3 * v], dy = P[3 * o + 1] - P[3 * v + 1], dz = P[3 * o + 2] - P[3 * v + 2];
    let ex = dx * mesh.LX[3 * t] + dy * mesh.LX[3 * t + 1] + dz * mesh.LX[3 * t + 2];
    let ey = dx * mesh.LY[3 * t] + dy * mesh.LY[3 * t + 1] + dz * mesh.LY[3 * t + 2];
    const l = hypot2(ex, ey); if (l > 0) { ex /= l; ey /= l; }
    E0[6 * t + 2 * j] = ex; E0[6 * t + 2 * j + 1] = ey;
  }
  Object.assign(mesh, { vtxStart, vtxTri, vtxCorner, cornerAngle: ang, cornerBegin, cornerOther, cornerE0: E0,
    vertexT: vT, vertexB: vB, vertexR: vR, vertexOrient: vOrient });
}

// ---- coordinate helpers -----------------------------------------------------
export function bary(mesh, t, a, b, out) {
  const a2 = b / mesh.y2[t];
  const a1 = (a - a2 * mesh.x2[t]) / mesh.x1[t];
  out = out || [0, 0, 0];
  out[0] = 1 - a1 - a2; out[1] = a1; out[2] = a2;
  return out;
}
export function localFromBary(mesh, t, alpha) {
  return [alpha[1] * mesh.x1[t] + alpha[2] * mesh.x2[t], alpha[2] * mesh.y2[t]];
}
export function localToGlobal(mesh, t, a, b) {
  const v0 = mesh.triangles[3 * t];
  const P = mesh.points, LX = mesh.LX, LY = mesh.LY;
  return [P[3 * v0] + a * LX[3 * t] + b * LY[3 * t], P[3 * v0 + 1] + a * LX[3 * t + 1] + b * LY[3 * t + 1], P[3 * v0 + 2] + a * LX[3 * t + 2] + b * LY[3 * t + 2]];
}
export function globalToLocal(mesh, t, g) {
  const v0 = mesh.triangles[3 * t];
  const P = mesh.points, LX = mesh.LX, LY = mesh.LY;
  const dx = g[0] - P[3 * v0], dy = g[1] - P[3 * v0 + 1], dz = g[2] - P[3 * v0 + 2];
  return [dx * LX[3 * t] + dy * LX[3 * t + 1] + dz * LX[3 * t + 2], dx * LY[3 * t] + dy * LY[3 * t + 1] + dz * LY[3 * t + 2]];
}
export function vectorToGlobal(mesh, t, vx, vy) {
  const LX = mesh.LX, LY = mesh.LY;
  return [vx * LX[3 * t] + vy * LY[3 * t], vx * LX[3 * t + 1] + vy * LY[3 * t + 1], vx * LX[3 * t + 2] + vy * LY[3 * t + 2]];
}

// ---------------------------------------------------------------------------
// vector field: tangential projection and transport into the triangles
// ---------------------------------------------------------------------------
export function buildField(mesh, vectors, scale = 1.0) {
  const { nPoints: N, nTris: M, triangles: tri, vertexNormals: VN, vertexT: T, vertexB: B } = mesh;
  const raw = new Float64Array(3 * N), tvec = new Float64Array(3 * N), mag = new Float64Array(N), tAngle = new Float64Array(N);
  for (let i = 0; i < N; i++) {
    const gx = vectors[3 * i] * scale, gy = vectors[3 * i + 1] * scale, gz = vectors[3 * i + 2] * scale;
    raw[3 * i] = gx; raw[3 * i + 1] = gy; raw[3 * i + 2] = gz;
    const nx = VN[3 * i], ny = VN[3 * i + 1], nz = VN[3 * i + 2];
    const d = gx * nx + gy * ny + gz * nz;
    const tx = gx - d * nx, ty = gy - d * ny, tz = gz - d * nz;
    tvec[3 * i] = tx; tvec[3 * i + 1] = ty; tvec[3 * i + 2] = tz;
    mag[i] = Math.sqrt(tx * tx + ty * ty + tz * tz);
    const a = tx * T[3 * i] + ty * T[3 * i + 1] + tz * T[3 * i + 2];
    const b = tx * B[3 * i] + ty * B[3 * i + 1] + tz * B[3 * i + 2];
    let ang = Math.atan2(b, a); if (ang < 0) ang += TWO_PI;
    tAngle[i] = ang;
  }
  const D = new Float64Array(6 * M);
  const E0 = mesh.cornerE0;
  for (let t = 0; t < M; t++) for (let j = 0; j < 3; j++) {
    const v = tri[3 * t + j];
    const phi = wrapPi(tAngle[v] - mesh.cornerBegin[3 * t + j]);
    const c = Math.cos(phi), s = Math.sin(phi);
    const ex = E0[6 * t + 2 * j], ey = E0[6 * t + 2 * j + 1];
    D[6 * t + 2 * j] = (c * ex - s * ey) * mag[v];
    D[6 * t + 2 * j + 1] = (s * ex + c * ey) * mag[v];
  }
  return { mesh, raw, tvec, mag, tAngle, D, scale };
}

export function evaluate(mesh, D, t, alpha) {
  return [alpha[0] * D[6 * t] + alpha[1] * D[6 * t + 2] + alpha[2] * D[6 * t + 4],
          alpha[0] * D[6 * t + 1] + alpha[1] * D[6 * t + 3] + alpha[2] * D[6 * t + 5]];
}

export function jacobian(mesh, D, t) {
  const dD1x = D[6 * t + 2] - D[6 * t], dD1y = D[6 * t + 3] - D[6 * t + 1];
  const dD2x = D[6 * t + 4] - D[6 * t], dD2y = D[6 * t + 5] - D[6 * t + 1];
  const x1 = mesh.x1[t], x2 = mesh.x2[t], y2 = mesh.y2[t], det = x1 * y2;
  return [[dD1x * y2 / det, (-dD1x * x2 + dD2x * x1) / det], [dD1y * y2 / det, (-dD1y * x2 + dD2y * x1) / det]];
}

// ---------------------------------------------------------------------------
// fixed points
// ---------------------------------------------------------------------------
export function poincareIndices(field, excludeBoundary = true, zeroTol = 1e-10) {
  const { mesh, D, mag } = field;
  const M = mesh.nTris, tri = mesh.triangles;
  const index = new Int8Array(M);
  for (let t = 0; t < M; t++) {
    if (excludeBoundary && mesh.boundaryTriangle[t]) continue;
    if (mag[tri[3 * t]] < zeroTol || mag[tri[3 * t + 1]] < zeroTol || mag[tri[3 * t + 2]] < zeroTol) continue;
    const a0 = Math.atan2(D[6 * t + 1], D[6 * t]), a1 = Math.atan2(D[6 * t + 3], D[6 * t + 2]), a2 = Math.atan2(D[6 * t + 5], D[6 * t + 4]);
    let total = wrapPi(a0 - a2) + wrapPi(a1 - a0) + wrapPi(a2 - a1);
    if (mesh.y2[t] < 0) total = -total;
    index[t] = Math.round(total / TWO_PI);
  }
  return index;
}

function eig2(J) {
  const a = J[0][0], b = J[0][1], c = J[1][0], d = J[1][1];
  const tr = a + d, det = a * d - b * c;
  const disc = tr * tr / 4 - det;
  if (disc >= 0) {
    const s = Math.sqrt(disc);
    const l1 = tr / 2 + s, l2 = tr / 2 - s;
    const vec = l => {
      let v;
      if (Math.abs(b) > 1e-300) v = [b, l - a];
      else if (Math.abs(c) > 1e-300) v = [l - d, c];
      else v = Math.abs(l - a) <= Math.abs(l - d) ? [1, 0] : [0, 1];
      const n = hypot2(v[0], v[1]) || 1;
      return [v[0] / n, v[1] / n];
    };
    return { real: true, values: [[l1, 0], [l2, 0]], vectors: [vec(l1), vec(l2)] };
  }
  const im = Math.sqrt(-disc);
  return { real: false, values: [[tr / 2, im], [tr / 2, -im]], vectors: null };
}

export function classify(J, relTol = 1e-8) {
  const e = eig2(J);
  const scale = Math.max(...e.values.map(v => hypot2(v[0], v[1])), 1e-300);
  let code;
  if (e.real) {
    const r1 = e.values[0][0], r2 = e.values[1][0], small = relTol * scale;
    if (Math.abs(r1) <= small || Math.abs(r2) <= small) code = UNKNOWN;
    else if (r1 > 0 && r2 > 0) code = SOURCE;
    else if (r1 < 0 && r2 < 0) code = SINK;
    else code = SADDLE;
  } else {
    const r = e.values[0][0];
    code = Math.abs(r) <= relTol * scale ? CENTER : (r < 0 ? ATTRACTING_FOCUS : REPELLING_FOCUS);
  }
  return { code, eig: e };
}

export function findFixedPoints(field, opts = {}) {
  const excludeBoundary = opts.excludeBoundary !== false;
  const zeroTol = opts.zeroTol ?? 1e-10;
  const { mesh, D } = field;
  const index = poincareIndices(field, excludeBoundary, zeroTol);
  const fps = [];
  for (let t = 0; t < mesh.nTris; t++) {
    if (index[t] === 0) continue;
    // solve [[D0x D1x D2x],[D0y D1y D2y],[1 1 1]] alpha = (0 0 1)
    const m00 = D[6 * t], m01 = D[6 * t + 2], m02 = D[6 * t + 4];
    const m10 = D[6 * t + 1], m11 = D[6 * t + 3], m12 = D[6 * t + 5];
    const det = m00 * (m11 - m12) - m01 * (m10 - m12) + m02 * (m10 - m11);
    let alpha = [1 / 3, 1 / 3, 1 / 3], located = false;
    if (Math.abs(det) > 1e-300) {
      // Cramer's rule with rhs (0,0,1)
      const a0 = (m01 * m12 - m02 * m11) / det;
      const a1 = (m02 * m10 - m00 * m12) / det;
      const a2 = (m00 * m11 - m01 * m10) / det;
      if ([a0, a1, a2].every(x => Number.isFinite(x) && x >= -1e-6 && x <= 1 + 1e-6)) {
        const c = [Math.max(a0, 0), Math.max(a1, 0), Math.max(a2, 0)];
        const s = c[0] + c[1] + c[2];
        alpha = [c[0] / s, c[1] / s, c[2] / s];
        located = true;
      }
    }
    const [la, lb] = localFromBary(mesh, t, alpha);
    const pos = localToGlobal(mesh, t, la, lb);
    const J = jacobian(mesh, D, t);
    const { code, eig } = classify(J);
    const v = evaluate(mesh, D, t, alpha);
    const fp = { id: fps.length, triangle: t, position: pos, local: [la, lb], alpha, index: index[t], type: code,
      type_name: TYPE_NAMES[code], eigenvalues: eig.values, jacobian: J, residual: hypot2(v[0], v[1]), located,
      eigvec_out_local: null, eigvec_in_local: null, eigvec_out: null, eigvec_in: null };
    if (code === SADDLE) {
      const order = eig.values[0][0] < eig.values[1][0] ? [0, 1] : [1, 0];
      fp.eigvec_in_local = eig.vectors[order[0]];
      fp.eigvec_out_local = eig.vectors[order[1]];
      fp.eigvec_in = vectorToGlobal(mesh, t, fp.eigvec_in_local[0], fp.eigvec_in_local[1]);
      fp.eigvec_out = vectorToGlobal(mesh, t, fp.eigvec_out_local[0], fp.eigvec_out_local[1]);
    }
    fps.push(fp);
  }
  // merge duplicates (a zero on a shared edge is reported by both triangles)
  const tol = 1e-6 * mesh.meanEdgeLength;
  const sorted = fps.slice().sort((p, q) => p.residual - q.residual);
  const kept = [];
  for (const fp of sorted) {
    if (kept.every(k => Math.sqrt((k.position[0] - fp.position[0]) ** 2 + (k.position[1] - fp.position[1]) ** 2 + (k.position[2] - fp.position[2]) ** 2) > tol)) kept.push(fp);
  }
  kept.sort((p, q) => p.triangle - q.triangle);
  kept.forEach((fp, i) => { fp.id = i; });
  return { fixedPoints: kept, index };
}

// ---------------------------------------------------------------------------
// tracer (single particle, arc length or time)
// ---------------------------------------------------------------------------
function vertexTransfer(mesh, D, v, dir) {
  const s0 = mesh.vtxStart[v], s1 = mesh.vtxStart[v + 1];
  let best = null, bestMargin = -Infinity;
  for (let k = s0; k < s1; k++) {
    const t = mesh.vtxTri[k], j = mesh.vtxCorner[k];
    let dx = dir * D[6 * t + 2 * j], dy = dir * D[6 * t + 2 * j + 1];
    const nd = hypot2(dx, dy);
    if (nd < 1e-300) continue;
    dx /= nd; dy /= nd;
    const xy = [[0, 0], [mesh.x1[t], 0], [mesh.x2[t], mesh.y2[t]]];
    const e1 = [xy[(j + 1) % 3][0] - xy[j][0], xy[(j + 1) % 3][1] - xy[j][1]];
    const e2 = [xy[(j + 2) % 3][0] - xy[j][0], xy[(j + 2) % 3][1] - xy[j][1]];
    const det = e1[0] * e2[1] - e1[1] * e2[0];
    if (Math.abs(det) < 1e-300) continue;
    const al = (dx * e2[1] - dy * e2[0]) / det, be = (e1[0] * dy - e1[1] * dx) / det;
    const margin = Math.min(al * hypot2(e1[0], e1[1]), be * hypot2(e2[0], e2[1]));
    if (margin > bestMargin) { bestMargin = margin; best = [t, j]; }
  }
  if (best === null) return { best: null, why: 'zero' };
  if (bestMargin < -1e-9) return { best: null, why: 'outside' };
  return { best, why: 'ok' };
}

/** Advance one particle {tri, a, b, status, prevEdge, slide} by a step h. */
export function stepParticle(mesh, D, p, h, dir, opts = {}) {
  const normalize = !!opts.normalize, rk4 = opts.integrator !== 'euler';
  const maxRounds = opts.maxRounds ?? 50, tol = 1e-9, vtol = 1e-9;
  let rem = h, used = 0, speed = 0, first = true;
  const alphaP = [0, 0, 0], alphaQ = [0, 0, 0];
  const vel = (t, a, b) => {
    bary(mesh, t, a, b, alphaQ);
    let vx = alphaQ[0] * D[6 * t] + alphaQ[1] * D[6 * t + 2] + alphaQ[2] * D[6 * t + 4];
    let vy = alphaQ[0] * D[6 * t + 1] + alphaQ[1] * D[6 * t + 3] + alphaQ[2] * D[6 * t + 5];
    if (normalize) { const n = hypot2(vx, vy); if (n > 0) { vx /= n; vy /= n; } else { vx = vy = 0; } }
    return [vx * dir, vy * dir];
  };
  const rk4Step = (t, a, b, hh, k1) => {
    k1 = k1 || vel(t, a, b);
    const k2 = vel(t, a + 0.5 * hh * k1[0], b + 0.5 * hh * k1[1]);
    const k3 = vel(t, a + 0.5 * hh * k2[0], b + 0.5 * hh * k2[1]);
    const k4 = vel(t, a + hh * k3[0], b + hh * k3[1]);
    return [a + hh / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]), b + hh / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])];
  };
  for (let round = 0; round < maxRounds && rem > 1e-14 * h; round++) {
    const t = p.tri;
    bary(mesh, t, p.a, p.b, alphaP);
    const ap = [alphaP[0], alphaP[1], alphaP[2]];
    if (first) {
      const vx = ap[0] * D[6 * t] + ap[1] * D[6 * t + 2] + ap[2] * D[6 * t + 4];
      const vy = ap[0] * D[6 * t + 1] + ap[1] * D[6 * t + 3] + ap[2] * D[6 * t + 5];
      speed = hypot2(vx, vy); first = false;
    }
    let v1 = vel(t, p.a, p.b);
    if (p.slide) {
      let l = -1;
      for (let k = 0; k < 3; k++) if (mesh.triEdges[3 * t + k] === p.prevEdge) l = k;
      if (l >= 0) {
        const xy = [[0, 0], [mesh.x1[t], 0], [mesh.x2[t], mesh.y2[t]]];
        let ex = xy[(l + 1) % 3][0] - xy[l][0], ey = xy[(l + 1) % 3][1] - xy[l][1];
        const n = hypot2(ex, ey) || 1; ex /= n; ey /= n;
        const d = v1[0] * ex + v1[1] * ey;
        v1 = [d * ex, d * ey];
      }
    }
    let q;
    if (rk4 && !p.slide) q = rk4Step(t, p.a, p.b, rem, v1);
    else q = [p.a + rem * v1[0], p.b + rem * v1[1]];
    bary(mesh, t, q[0], q[1], alphaQ);
    const aq = [alphaQ[0], alphaQ[1], alphaQ[2]];
    if (aq[0] >= -tol && aq[1] >= -tol && aq[2] >= -tol) {
      p.a = q[0]; p.b = q[1]; used += rem; rem = 0; p.slide = false;
      break;
    }
    // exit through an edge
    let s = Infinity, jstar = -1;
    for (let j = 0; j < 3; j++) {
      if (aq[j] < -tol) {
        const den = ap[j] - aq[j];
        let sj = den !== 0 ? ap[j] / den : Infinity;
        if (sj < 0) sj = 0;
        if (sj < s) { s = sj; jstar = j; }
      }
    }
    if (jstar < 0) { p.a = q[0]; p.b = q[1]; used += rem; rem = 0; break; }
    s = Math.min(Math.max(s, 0), 1);
    let ax = [ap[0] + s * (aq[0] - ap[0]), ap[1] + s * (aq[1] - ap[1]), ap[2] + s * (aq[2] - ap[2])];
    if (rk4 && !p.slide) {
      const fp = ap[jstar];
      for (let it = 0; it < 2; it++) {
        const q1 = rk4Step(t, p.a, p.b, s * rem, v1);
        bary(mesh, t, q1[0], q1[1], alphaQ);
        const den = fp - alphaQ[jstar];
        if (Math.abs(den) > 1e-300) s = Math.min(Math.max(s * fp / den, 0), 1);
      }
      const q1 = rk4Step(t, p.a, p.b, s * rem, v1);
      bary(mesh, t, q1[0], q1[1], alphaQ);
      ax = [alphaQ[0], alphaQ[1], alphaQ[2]];
    }
    ax[jstar] = 0;
    for (let k = 0; k < 3; k++) if (ax[k] < 0) ax[k] = 0;
    const sum = ax[0] + ax[1] + ax[2] || 1;
    ax = [ax[0] / sum, ax[1] / sum, ax[2] / sum];
    const dtPart = s * rem;
    used += dtPart; rem = Math.max(rem - dtPart, 0);
    const edgeLocal = (jstar + 1) % 3;
    const edgeId = mesh.triEdges[3 * t + edgeLocal];
    const nbr = mesh.triNeighbors[3 * t + edgeLocal];
    let kmax = 0; for (let k = 1; k < 3; k++) if (ax[k] > ax[kmax]) kmax = k;
    const atVertex = ax[kmax] > 1 - vtol;
    const sameEdge = edgeId === p.prevEdge && s < 1e-9 && !atVertex && !p.slide;
    const loc = localFromBary(mesh, t, ax);
    p.a = loc[0]; p.b = loc[1];
    if (sameEdge) { p.slide = true; p.prevEdge = edgeId; continue; }
    if (atVertex) {
      const v = mesh.triangles[3 * t + kmax];
      const { best, why } = vertexTransfer(mesh, D, v, dir);
      if (best === null) {
        p.status = (why === 'outside' || mesh.boundaryVertex[v]) ? LEFT : STUCK;
        rem = 0; break;
      }
      const [tn, jn] = best;
      const an = [0, 0, 0]; an[jn] = 1;
      const l2 = localFromBary(mesh, tn, an);
      p.tri = tn; p.a = l2[0]; p.b = l2[1]; p.prevEdge = -1; p.slide = false;
      continue;
    }
    if (nbr < 0) { p.status = LEFT; rem = 0; break; }
    const g = localToGlobal(mesh, t, p.a, p.b);
    const ln = globalToLocal(mesh, nbr, g);
    const an = bary(mesh, nbr, ln[0], ln[1]);
    for (let k = 0; k < 3; k++) if (an[k] < 0) an[k] = 0;
    const sn = an[0] + an[1] + an[2] || 1;
    const ln2 = localFromBary(mesh, nbr, [an[0] / sn, an[1] / sn, an[2] / sn]);
    p.tri = nbr; p.a = ln2[0]; p.b = ln2[1]; p.prevEdge = edgeId; p.slide = false;
  }
  if (rem > 1e-14 * h && p.status === ACTIVE) p.status = STUCK;
  return { used, speed };
}

// ---------------------------------------------------------------------------
// manifolds
// ---------------------------------------------------------------------------
function seedPoint(mesh, fp, dir2, eps) {
  const t = fp.triangle;
  let a = 0, b = 0;
  for (let i = 0; i < 12; i++) {
    a = fp.local[0] + eps * dir2[0]; b = fp.local[1] + eps * dir2[1];
    const al = bary(mesh, t, a, b);
    if (al[0] >= 0 && al[1] >= 0 && al[2] >= 0) return { tri: t, a, b };
    eps *= 0.5;
  }
  // fall back: the seed lies just across an edge -> use the neighbour
  a = fp.local[0] + 8 * eps * dir2[0]; b = fp.local[1] + 8 * eps * dir2[1];
  const al = bary(mesh, t, a, b);
  let jneg = 0; for (let k = 1; k < 3; k++) if (al[k] < al[jneg]) jneg = k;
  const nbr = mesh.triNeighbors[3 * t + (jneg + 1) % 3];
  if (nbr >= 0) {
    const g = localToGlobal(mesh, t, a, b);
    const ln = globalToLocal(mesh, nbr, g);
    const an = bary(mesh, nbr, ln[0], ln[1]);
    for (let k = 0; k < 3; k++) if (an[k] < 0) an[k] = 0;
    const s = an[0] + an[1] + an[2] || 1;
    const l2 = localFromBary(mesh, nbr, [an[0] / s, an[1] / s, an[2] / s]);
    return { tri: nbr, a: l2[0], b: l2[1] };
  }
  return { tri: t, a: fp.local[0], b: fp.local[1] };
}

export function computeManifolds(field, fixedPoints, opts = {}) {
  const stepFraction = opts.stepFraction ?? 0.2, maxSteps = opts.maxSteps ?? 20000;
  const perturbation = opts.perturbation ?? 0.1, captureRadius = opts.captureRadius ?? 0.5, maxCrossings = opts.maxCrossings ?? 50;
  const { mesh, D } = field;
  const maxLength = opts.maxLength || 25 * mesh.radius;
  const saddles = fixedPoints.filter(fp => fp.type === SADDLE);
  const branches = [];
  const sqArea = mesh.area.map(Math.sqrt);
  for (const fp of saddles) {
    const eps = perturbation * sqArea[fp.triangle];
    for (const [kind, evec, dir] of [[UNSTABLE, fp.eigvec_out_local, 1], [STABLE, fp.eigvec_in_local, -1]]) {
      for (const sign of [1, -1]) {
        const seed = seedPoint(mesh, fp, [sign * evec[0], sign * evec[1]], eps);
        const p = { tri: seed.tri, a: seed.a, b: seed.b, status: ACTIVE, prevEdge: -1, slide: false };
        const g = localToGlobal(mesh, p.tri, p.a, p.b);
        branches.push({ saddle: fp.id, kind, kind_name: kind === UNSTABLE ? 'unstable' : 'stable', sign, dir, p,
          pts: [g[0], g[1], g[2]], arc: [0], tris: [p.tri], spd: [], length: 0, leftStart: false, end_reason: '', end_fixed_point: -1 });
      }
    }
  }
  const fpRad = fixedPoints.map(fp => captureRadius * sqArea[fp.triangle]);
  for (let it = 0; it < maxSteps; it++) {
    let anyActive = false;
    for (const br of branches) {
      if (br.p.status !== ACTIVE) continue;
      anyActive = true;
      const h = stepFraction * sqArea[br.p.tri];
      const { used, speed } = stepParticle(mesh, D, br.p, h, br.dir, { normalize: true, integrator: 'rk4', maxRounds: maxCrossings });
      br.length += used;
      const g = localToGlobal(mesh, br.p.tri, br.p.a, br.p.b);
      br.pts.push(g[0], g[1], g[2]); br.arc.push(br.length); br.tris.push(br.p.tri); br.spd.push(speed);
      // termination
      let near = -1, nearD = Infinity, farFromOwn = true;
      for (let k = 0; k < fixedPoints.length; k++) {
        const q = fixedPoints[k].position;
        const d = Math.sqrt((g[0] - q[0]) ** 2 + (g[1] - q[1]) ** 2 + (g[2] - q[2]) ** 2);
        const own = fixedPoints[k].id === br.saddle;
        if (own && d <= 3 * fpRad[k]) farFromOwn = false;
        if (own && !br.leftStart) continue;
        if (d <= fpRad[k] && d < nearD) { nearD = d; near = k; }
      }
      if (farFromOwn) br.leftStart = true;
      if (br.p.status === LEFT) br.end_reason = 'left the domain';
      else if (br.p.status === STUCK) br.end_reason = 'zero field / vertex';
      else if (near >= 0) { br.end_reason = `reached fixed point ${fixedPoints[near].id} (${fixedPoints[near].type_name})`; br.end_fixed_point = fixedPoints[near].id; br.p.status = STUCK; }
      else if (br.length >= maxLength) { br.end_reason = 'maximum length'; br.p.status = STUCK; }
    }
    if (!anyActive) break;
  }
  for (const br of branches) {
    if (br.p.status === ACTIVE) br.end_reason = 'maximum number of steps';
    if (br.spd.length < br.arc.length) br.spd.unshift(br.spd.length ? br.spd[0] : 0);
  }
  return branches.map(br => ({ saddle: br.saddle, kind: br.kind, kind_name: br.kind_name, sign: br.sign,
    points: Float64Array.from(br.pts), arclength: Float64Array.from(br.arc), triangles: Int32Array.from(br.tris),
    speed: Float64Array.from(br.spd), length: br.length, end_reason: br.end_reason, end_fixed_point: br.end_fixed_point }));
}

// ---------------------------------------------------------------------------
// complete analysis
// ---------------------------------------------------------------------------
export function analyze(points, triangles, vectors, opts = {}) {
  const t0 = now();
  const mesh = buildMesh(points, triangles);
  const scale = opts.scale || 1.0;
  const field = buildField(mesh, vectors, scale);
  const t1 = now();
  const { fixedPoints, index } = findFixedPoints(field, opts);
  const t2 = now();
  const branches = opts.manifolds === false ? [] : computeManifolds(field, fixedPoints, opts);
  const t3 = now();
  let indexSum = 0; for (let t = 0; t < mesh.nTris; t++) indexSum += index[t];
  const magOut = new Float64Array(mesh.nPoints), vecOut = new Float64Array(3 * mesh.nPoints);
  for (let i = 0; i < mesh.nPoints; i++) { magOut[i] = field.mag[i] / scale; for (let k = 0; k < 3; k++) vecOut[3 * i + k] = field.tvec[3 * i + k] / scale; }
  return {
    n_points: mesh.nPoints, n_triangles: mesh.nTris, n_edges: mesh.nEdges, n_boundary_edges: mesh.nBoundaryEdges,
    points: mesh.points, triangles: mesh.triangles, magnitude: magOut, vectors: vecOut,
    bounds: mesh.bounds, center: mesh.center, radius: mesh.radius, mean_edge_length: mesh.meanEdgeLength,
    fixed_points: fixedPoints.map(fp => ({ id: fp.id, type: fp.type, type_name: fp.type_name, index: fp.index, triangle: fp.triangle,
      position: fp.position, eigenvalues: fp.eigenvalues.map(v => [v[0] / scale, v[1] / scale]),
      eigvec_out: fp.eigvec_out, eigvec_in: fp.eigvec_in, residual: fp.residual / scale, jacobian: fp.jacobian.map(r => r.map(x => x / scale)) })),
    branches, poincare_index_sum: indexSum, scale, warnings: mesh.warnings,
    timings: { mesh_and_field: t1 - t0, fixed_points: t2 - t1, manifolds: t3 - t2 },
  };
}

function now() { return typeof performance !== 'undefined' ? performance.now() : Date.now(); }

// ---------------------------------------------------------------------------
// exporters (legacy ASCII VTK text, CSV)
// ---------------------------------------------------------------------------
const f = x => Number(x).toPrecision(9).replace(/\.?0+$/, '').replace(/\.?0+e/, 'e');

export function fixedPointsCSV(res) {
  const rows = [['id', 'type', 'type_name', 'poincare_index', 'triangle', 'x', 'y', 'z', 'eig1_real', 'eig1_imag', 'eig2_real', 'eig2_imag', 'residual']];
  for (const fp of res.fixed_points) rows.push([fp.id, fp.type, fp.type_name, fp.index, fp.triangle, ...fp.position, fp.eigenvalues[0][0], fp.eigenvalues[0][1], fp.eigenvalues[1][0], fp.eigenvalues[1][1], fp.residual]);
  return rows.map(r => r.join(',')).join('\n') + '\n';
}

export function fixedPointsVTK(res) {
  const fps = res.fixed_points, n = fps.length;
  const L = ['# vtk DataFile Version 3.0', 'WSSLCS fixed points', 'ASCII', 'DATASET POLYDATA', `POINTS ${n} double`];
  for (const fp of fps) L.push(fp.position.map(f).join(' '));
  L.push(`VERTICES ${n} ${2 * n}`);
  for (let i = 0; i < n; i++) L.push(`1 ${i}`);
  L.push(`POINT_DATA ${n}`, 'SCALARS type int 1', 'LOOKUP_TABLE default');
  for (const fp of fps) L.push(String(fp.type));
  L.push('SCALARS poincare_index int 1', 'LOOKUP_TABLE default');
  for (const fp of fps) L.push(String(fp.index));
  L.push('SCALARS triangle int 1', 'LOOKUP_TABLE default');
  for (const fp of fps) L.push(String(fp.triangle));
  L.push('FIELD FieldData 4', `eigenvalue_real 2 ${n} double`);
  for (const fp of fps) L.push(`${f(fp.eigenvalues[0][0])} ${f(fp.eigenvalues[1][0])}`);
  L.push(`eigenvalue_imag 2 ${n} double`);
  for (const fp of fps) L.push(`${f(fp.eigenvalues[0][1])} ${f(fp.eigenvalues[1][1])}`);
  L.push(`eigvec_outgoing 3 ${n} double`);
  for (const fp of fps) L.push((fp.eigvec_out || [0, 0, 0]).map(f).join(' '));
  L.push(`eigvec_incoming 3 ${n} double`);
  for (const fp of fps) L.push((fp.eigvec_in || [0, 0, 0]).map(f).join(' '));
  return L.join('\n') + '\n';
}

export function manifoldsVTK(res) {
  const brs = res.branches;
  let nPts = 0; for (const br of brs) nPts += br.points.length / 3;
  const L = ['# vtk DataFile Version 3.0', 'WSSLCS manifolds (0 unstable, 1 stable)', 'ASCII', 'DATASET POLYDATA', `POINTS ${nPts} double`];
  for (const br of brs) for (let i = 0; i < br.points.length; i += 3) L.push(`${f(br.points[i])} ${f(br.points[i + 1])} ${f(br.points[i + 2])}`);
  L.push(`LINES ${brs.length} ${nPts + brs.length}`);
  let off = 0;
  for (const br of brs) { const m = br.points.length / 3; const ids = []; for (let i = 0; i < m; i++) ids.push(off + i); L.push(`${m} ${ids.join(' ')}`); off += m; }
  L.push(`CELL_DATA ${brs.length}`, 'SCALARS manifold int 1', 'LOOKUP_TABLE default');
  for (const br of brs) L.push(String(br.kind));
  L.push('SCALARS saddle int 1', 'LOOKUP_TABLE default');
  for (const br of brs) L.push(String(br.saddle));
  L.push('SCALARS branch int 1', 'LOOKUP_TABLE default');
  for (const br of brs) L.push(String(br.sign));
  L.push('SCALARS length double 1', 'LOOKUP_TABLE default');
  for (const br of brs) L.push(f(br.length));
  L.push(`POINT_DATA ${nPts}`, 'SCALARS arclength double 1', 'LOOKUP_TABLE default');
  for (const br of brs) for (let i = 0; i < br.arclength.length; i++) L.push(f(br.arclength[i]));
  L.push('SCALARS speed double 1', 'LOOKUP_TABLE default');
  for (const br of brs) for (let i = 0; i < br.arclength.length; i++) L.push(f((br.speed[i] ?? 0) / res.scale));
  return L.join('\n') + '\n';
}
