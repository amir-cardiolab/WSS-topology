// worker.js — runs the fixed point / manifold analysis off the main thread
import { analyze } from './wsslcs-core.js';

self.onmessage = e => {
  const { id, points, triangles, vectors, options } = e.data;
  try {
    const res = analyze(points, triangles, vectors, options);
    const transfer = [res.points.buffer, res.triangles.buffer, res.magnitude.buffer, res.vectors.buffer, res.tri_normals.buffer, res.vertex_normals.buffer];
    for (const br of res.branches) transfer.push(br.points.buffer, br.arclength.buffer, br.triangles.buffer, br.speed.buffer);
    for (const sl of res.streamlines) transfer.push(sl.points.buffer, sl.triangles.buffer, sl.speed.buffer);
    self.postMessage({ id, result: res }, transfer);
  } catch (err) {
    self.postMessage({ id, error: String(err && err.message ? err.message : err) });
  }
};
