// vtk-reader.js — reads triangulated surfaces with point vector arrays from
// legacy VTK files (ASCII or binary, PolyData / UnstructuredGrid, old and
// VTK 5.1 cell layouts) and XML files (.vtp / .vtu: ascii, inline base64 or
// appended raw/base64 data, optionally zlib-compressed).  No DOM required.
//
// readVTK(bufferOrBytes, {inflate}) -> Promise<{points: Float64Array(N*3),
//     triangles: Int32Array(M*3), pointArrays: [{name, ncomp, data}],
//     cellArrays: [...], activeVectors, warnings: []}>

const TYPE_SIZE = { Int8: 1, UInt8: 1, Int16: 2, UInt16: 2, Int32: 4, UInt32: 4, Int64: 8, UInt64: 8, Float32: 4, Float64: 8 };
const LEGACY_TYPE = { bit: 'UInt8', unsigned_char: 'UInt8', char: 'Int8', unsigned_short: 'UInt16', short: 'Int16', unsigned_int: 'UInt32', int: 'Int32',
  unsigned_long: 'UInt64', long: 'Int64', float: 'Float32', double: 'Float64', vtktypeint64: 'Int64', vtktypeuint64: 'UInt64', vtktypeint32: 'Int32', vtkidtype: 'Int32' };
const WS = new Uint8Array(256); WS[32] = WS[9] = WS[10] = WS[13] = WS[12] = 1;

export function bytesToString(u8, start = 0, end = u8.length) {
  let s = '';
  for (let i = start; i < end; i += 8192) s += String.fromCharCode.apply(null, u8.subarray(i, Math.min(i + 8192, end)));
  return s;
}

function readTyped(dv, offset, count, type, little) {
  const out = new Float64Array(count), size = TYPE_SIZE[type];
  switch (type) {
    case 'Float32': for (let i = 0; i < count; i++) out[i] = dv.getFloat32(offset + 4 * i, little); break;
    case 'Float64': for (let i = 0; i < count; i++) out[i] = dv.getFloat64(offset + 8 * i, little); break;
    case 'Int32': for (let i = 0; i < count; i++) out[i] = dv.getInt32(offset + 4 * i, little); break;
    case 'UInt32': for (let i = 0; i < count; i++) out[i] = dv.getUint32(offset + 4 * i, little); break;
    case 'Int16': for (let i = 0; i < count; i++) out[i] = dv.getInt16(offset + 2 * i, little); break;
    case 'UInt16': for (let i = 0; i < count; i++) out[i] = dv.getUint16(offset + 2 * i, little); break;
    case 'Int8': for (let i = 0; i < count; i++) out[i] = dv.getInt8(offset + i); break;
    case 'UInt8': for (let i = 0; i < count; i++) out[i] = dv.getUint8(offset + i); break;
    case 'Int64': for (let i = 0; i < count; i++) out[i] = Number(dv.getBigInt64(offset + 8 * i, little)); break;
    case 'UInt64': for (let i = 0; i < count; i++) out[i] = Number(dv.getBigUint64(offset + 8 * i, little)); break;
    default: throw new Error('unsupported data type ' + type);
  }
  return [out, offset + count * size];
}

function parseNumbersFromText(text, count) {
  const out = count != null ? new Float64Array(count) : [];
  let k = 0;
  const re = /[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?|nan|inf/gi;
  let m;
  while ((m = re.exec(text)) !== null) {
    out[k++] = parseFloat(m[0]);
    if (count != null && k >= count) break;
  }
  if (count != null && k < count) throw new Error(`expected ${count} numbers, found ${k}`);
  return count != null ? out : Float64Array.from(out);
}

// ---------------------------------------------------------------------------
// legacy format
// ---------------------------------------------------------------------------
class Cursor {
  constructor(u8) { this.u8 = u8; this.pos = 0; this.dv = new DataView(u8.buffer, u8.byteOffset, u8.byteLength); }
  eof() { return this.pos >= this.u8.length; }
  readLine() {
    const u8 = this.u8, n = u8.length;
    let i = this.pos;
    while (i < n && u8[i] !== 10) i++;
    let s = bytesToString(u8, this.pos, i);
    this.pos = Math.min(i + 1, n);
    if (s.endsWith('\r')) s = s.slice(0, -1);
    return s;
  }
  peekLine() { const p = this.pos; const s = this.readLine(); this.pos = p; return s; }
  readAscii(count) {
    const u8 = this.u8, n = u8.length, out = new Float64Array(count);
    let i = this.pos, k = 0;
    while (k < count) {
      while (i < n && WS[u8[i]]) i++;
      let j = i;
      while (j < n && !WS[u8[j]]) j++;
      if (j === i) break;
      out[k++] = parseFloat(bytesToString(u8, i, j));
      i = j;
    }
    if (k < count) throw new Error(`unexpected end of file (read ${k} of ${count} values)`);
    this.pos = i;
    return out;
  }
  readBinary(count, type) {
    const [out, next] = readTyped(this.dv, this.pos, count, type, false);   // legacy binary is big-endian
    this.pos = next;
    return out;
  }
}

function legacyCellsToTriangles(conn, offsets, types, warnings) {
  // conn/offsets in VTK 5.1 layout: offsets has ncells+1 entries
  const tris = [];
  let skipped = 0;
  const nconn = conn.length;
  for (let c = 0; c < offsets.length - 1; c++) {
    const s = offsets[c], e = offsets[c + 1], m = e - s;
    if (!(m >= 0 && m <= 1024) || s < 0 || e > nconn) throw new Error('corrupt cell connectivity (unsupported encoding?)');
    const ty = types ? types[c] : (m === 3 ? 5 : m === 4 ? 9 : 7);
    if (ty === 5 && m === 3) tris.push(conn[s], conn[s + 1], conn[s + 2]);
    else if ((ty === 9 || ty === 7) && m >= 3) { for (let k = 1; k < m - 1; k++) tris.push(conn[s], conn[s + k], conn[s + k + 1]); }
    else if (ty === 6 && m >= 3) { for (let k = 0; k < m - 2; k++) tris.push(conn[s + k], conn[s + k + 1 + (k % 2)], conn[s + k + 2 - (k % 2)]); }
    else skipped++;
  }
  if (skipped) warnings.push(`${skipped} non-triangle cells skipped`);
  return Int32Array.from(tris);
}

function oldStyleCells(flat) {
  // "n v0 v1 ... " repeated -> conn, offsets
  const conn = [], offsets = [0];
  let i = 0;
  while (i < flat.length) {
    const m = flat[i++];
    if (!(m >= 0 && m <= 1024) || i + m > flat.length) throw new Error('corrupt cell list (unsupported encoding?)');
    for (let k = 0; k < m; k++) conn.push(flat[i++]);
    offsets.push(conn.length);
  }
  return [conn, offsets];
}

function readLegacy(u8, opts = {}) {
  const log = opts.log || (() => {});
  const cur = new Cursor(u8);
  const warnings = [];
  cur.readLine();                         // # vtk DataFile Version x.y
  cur.readLine();                         // title
  const fmt = cur.readLine().trim().toUpperCase();
  const binary = fmt === 'BINARY';
  if (!binary && fmt !== 'ASCII') throw new Error(`unknown legacy file format '${fmt}'`);
  let dataset = '';
  let points = null, triangles = null;
  const pointArrays = [], cellArrays = [];
  let mode = null, modeCount = 0, activeVectors = null;
  const readValues = (count, ltype) => binary ? cur.readBinary(count, LEGACY_TYPE[ltype.toLowerCase()] || 'Float32') : cur.readAscii(count);
  const readCellArray = (nOffsetsOrCells, size) => {
    // returns [conn, offsets] (VTK 5.1 layout)
    const peek = cur.peekLine().trim().toUpperCase();
    if (peek.startsWith('OFFSETS')) {
      const ot = cur.readLine().trim().split(/\s+/)[1];
      const offsets = readValues(nOffsetsOrCells, ot);
      while (cur.peekLine().trim() === '') cur.readLine();
      const ct = cur.readLine().trim().split(/\s+/)[1];
      const conn = readValues(size, ct);
      return [Array.from(conn), Array.from(offsets)];
    }
    const flat = readValues(size, 'int');
    return oldStyleCells(flat);
  };
  let pendingCells = null, pendingTypes = null;
  while (!cur.eof()) {
    const line = cur.readLine().trim();
    if (!line) continue;
    const tok = line.split(/\s+/);
    const key = tok[0].toUpperCase();
    log('legacy keyword', line.slice(0, 60), 'at byte', cur.pos);
    if (key === 'DATASET') { dataset = tok[1].toUpperCase(); continue; }
    if (key === 'METADATA' || key === 'INFORMATION' || key === 'NAME') { continue; }
    if (key === 'DATA' && tok[1] && tok[1].includes(':')) continue;
    if (key === 'POINTS') {
      const n = parseInt(tok[1], 10);
      points = readValues(3 * n, tok[2] || 'float');
      continue;
    }
    if (key === 'POLYGONS' || key === 'TRIANGLE_STRIPS') {
      const n = parseInt(tok[1], 10), size = parseInt(tok[2], 10);
      const [conn, offsets] = readCellArray(n, size);
      const types = key === 'TRIANGLE_STRIPS' ? new Array(offsets.length - 1).fill(6) : null;
      const t = legacyCellsToTriangles(conn, offsets, types, warnings);
      triangles = triangles ? Int32Array.from([...triangles, ...t]) : t;
      continue;
    }
    if (key === 'LINES' || key === 'VERTICES') {
      const n = parseInt(tok[1], 10), size = parseInt(tok[2], 10);
      readCellArray(n, size);
      continue;
    }
    if (key === 'CELLS') {
      const n = parseInt(tok[1], 10), size = parseInt(tok[2], 10);
      pendingCells = readCellArray(n, size);
      continue;
    }
    if (key === 'CELL_TYPES') {
      const n = parseInt(tok[1], 10);
      pendingTypes = Array.from(readValues(n, 'int'));
      if (pendingCells) { triangles = legacyCellsToTriangles(pendingCells[0], pendingCells[1], pendingTypes, warnings); pendingCells = null; }
      continue;
    }
    if (key === 'POINT_DATA') { mode = 'point'; modeCount = parseInt(tok[1], 10); continue; }
    if (key === 'CELL_DATA') { mode = 'cell'; modeCount = parseInt(tok[1], 10); continue; }
    const target = mode === 'cell' ? cellArrays : pointArrays;
    if (key === 'SCALARS') {
      const name = tok[1], type = tok[2] || 'float', ncomp = tok[3] ? parseInt(tok[3], 10) : 1;
      const lt = cur.peekLine().trim().toUpperCase();
      if (lt.startsWith('LOOKUP_TABLE')) cur.readLine();
      const data = readValues(modeCount * ncomp, type);
      target.push({ name, ncomp, data });
      continue;
    }
    if (key === 'VECTORS' || key === 'NORMALS') {
      const name = tok[1], type = tok[2] || 'float';
      const data = readValues(modeCount * 3, type);
      target.push({ name, ncomp: 3, data });
      if (key === 'VECTORS' && mode === 'point' && !activeVectors) activeVectors = name;
      continue;
    }
    if (key === 'TENSORS') { target.push({ name: tok[1], ncomp: 9, data: readValues(modeCount * 9, tok[2] || 'float') }); continue; }
    if (key === 'TEXTURE_COORDINATES') { const dim = parseInt(tok[2], 10); target.push({ name: tok[1], ncomp: dim, data: readValues(modeCount * dim, tok[3] || 'float') }); continue; }
    if (key === 'COLOR_SCALARS') { const nc = parseInt(tok[2], 10); target.push({ name: tok[1], ncomp: nc, data: binary ? cur.readBinary(modeCount * nc, 'UInt8') : cur.readAscii(modeCount * nc) }); continue; }
    if (key === 'LOOKUP_TABLE') { if (tok.length > 2) { const size = parseInt(tok[2], 10); if (binary) cur.readBinary(size * 4, 'UInt8'); else cur.readAscii(size * 4); } continue; }
    if (key === 'FIELD') {
      const k = parseInt(tok[2], 10);
      for (let i = 0; i < k; i++) {
        let l = cur.readLine().trim();
        while (!cur.eof() && (!l || /^(METADATA|INFORMATION|NAME |DATA )/i.test(l))) {
          l = cur.readLine().trim();
        }
        if (!l) break;
        const a = l.split(/\s+/);
        const name = a[0], ncomp = parseInt(a[1], 10), ntup = parseInt(a[2], 10), type = a[3] || 'float';
        if (!(ncomp > 0) || !(ntup >= 0)) throw new Error('cannot parse FIELD array header: ' + l);
        const data = readValues(ncomp * ntup, type);
        if (mode === 'point' || mode === 'cell') target.push({ name, ncomp, data });
        else if (points && ntup === points.length / 3) pointArrays.push({ name, ncomp, data });
      }
      continue;
    }
    // unknown keyword: ignore
  }
  if (!points) throw new Error('no POINTS section found');
  if (!triangles || triangles.length === 0) throw new Error('no polygons / triangle cells found in the file');
  return { points, triangles, pointArrays, cellArrays, activeVectors, warnings, dataset, format: binary ? 'legacy binary' : 'legacy ascii' };
}

// ---------------------------------------------------------------------------
// XML format (.vtp / .vtu)
// ---------------------------------------------------------------------------
function attrs(tag) {
  const out = {};
  const re = /([A-Za-z_][\w.:-]*)\s*=\s*"([^"]*)"/g;
  let m;
  while ((m = re.exec(tag)) !== null) out[m[1]] = m[2];
  return out;
}

const B64 = new Int16Array(256).fill(-1);
'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/'.split('').forEach((c, i) => { B64[c.charCodeAt(0)] = i; });
B64['-'.charCodeAt(0)] = 62; B64['_'.charCodeAt(0)] = 63;

function base64Decode(str, start = 0, maxChars = Infinity) {
  // decodes from str[start] until '=' padding, end of string, or a non-base64 character
  const bytes = [];
  let buf = 0, bits = 0, i = start, n = Math.min(str.length, start + maxChars);
  for (; i < n; i++) {
    const c = str.charCodeAt(i);
    if (c === 61) break;                     // '='
    const v = c < 256 ? B64[c] : -1;
    if (v < 0) { if (c === 10 || c === 13 || c === 32 || c === 9) continue; break; }
    buf = (buf << 6) | v; bits += 6;
    if (bits >= 8) { bits -= 8; bytes.push((buf >> bits) & 255); }
  }
  return new Uint8Array(bytes);
}

async function defaultInflate(bytes) {
  if (typeof DecompressionStream === 'undefined') throw new Error('this file uses zlib compression, which this browser cannot decompress');
  const ds = new DecompressionStream('deflate');
  const writer = ds.writable.getWriter();
  writer.write(bytes); writer.close();
  const reader = ds.readable.getReader();
  const chunks = [];
  let total = 0;
  for (;;) { const { value, done } = await reader.read(); if (done) break; chunks.push(value); total += value.length; }
  const out = new Uint8Array(total);
  let o = 0; for (const c of chunks) { out.set(c, o); o += c.length; }
  return out;
}

async function decodeBlocks(getHeaderBytes, getDataBytes, headerType, compressed, little, inflate, getWhole) {
  // getHeaderBytes(nInts) -> Uint8Array with at least nInts header integers
  // getDataBytes(headerByteLength, dataByteLength) -> Uint8Array of the data
  // getWhole(nBytes) -> (base64 only) the first nBytes of the contiguous header+data stream
  const hs = TYPE_SIZE[headerType];
  const hdr = (bytes, n) => readTyped(new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength), 0, n, headerType, little)[0];
  if (!compressed) {
    const h = hdr(getHeaderBytes(1), 1);
    const nbytes = h[0];
    if (getWhole) return getWhole(hs + nbytes).subarray(hs);
    return getDataBytes(hs, nbytes);
  }
  const h0 = hdr(getHeaderBytes(3), 3);
  const nblocks = h0[0], blockSize = h0[1], lastSize = h0[2];
  if (!(nblocks > 0)) return new Uint8Array(0);
  const h = hdr(getHeaderBytes(3 + nblocks), 3 + nblocks);
  const sizes = Array.from(h.subarray(3, 3 + nblocks));
  const totalComp = sizes.reduce((a, b) => a + b, 0);
  const comp = getDataBytes(hs * (3 + nblocks), totalComp);
  const out = new Uint8Array(blockSize * (nblocks - 1) + (nblocks ? (lastSize || blockSize) : 0));
  let off = 0, oo = 0;
  for (let b = 0; b < nblocks; b++) {
    const chunk = await inflate(comp.subarray(off, off + sizes[b]));
    out.set(chunk, oo); oo += chunk.length; off += sizes[b];
  }
  return out.subarray(0, oo);
}

async function readXML(u8, opts) {
  const inflate = opts.inflate || defaultInflate;
  const log = opts.log || (() => {});
  const warnings = [];
  // locate the appended data section without decoding binary content as text
  const marker = '<AppendedData';
  let appIdx = -1;
  {
    const m0 = marker.charCodeAt(0);
    for (let i = 0; i < u8.length - marker.length; i++) {
      if (u8[i] !== m0) continue;
      let ok = true;
      for (let k = 1; k < marker.length; k++) if (u8[i + k] !== marker.charCodeAt(k)) { ok = false; break; }
      if (ok) { appIdx = i; break; }
    }
  }
  const xmlText = bytesToString(u8, 0, appIdx >= 0 ? appIdx : u8.length);
  let appEncoding = null, appStart = -1;
  if (appIdx >= 0) {
    let i = appIdx; while (i < u8.length && u8[i] !== 62) i++;           // '>'
    appEncoding = attrs(bytesToString(u8, appIdx, i + 1)).encoding || 'base64';
    let j = i + 1; while (j < u8.length && u8[j] !== 95) j++;             // '_'
    appStart = j + 1;
  }
  const root = attrs((xmlText.match(/<VTKFile\b[^>]*>/) || [''])[0]);
  const type = root.type || 'PolyData';
  const little = (root.byte_order || 'LittleEndian') !== 'BigEndian';
  const headerType = root.header_type || 'UInt32';
  const compressor = root.compressor || '';
  if (compressor && !/ZLib/i.test(compressor)) throw new Error(`unsupported compressor ${compressor} (only zlib is supported)`);
  const piece = attrs((xmlText.match(/<Piece\b[^>]*>/) || [''])[0]);
  let appText = null;
  const getAppText = () => { if (appText === null) appText = bytesToString(u8, appStart, u8.length); return appText; };

  async function decodeArray(tag, inner) {
    const a = attrs(tag);
    log('xml array', a.Name, a.type, a.format, a.offset, 'inner', inner.length);
    const dtype = a.type || 'Float32', ncomp = parseInt(a.NumberOfComponents || '1', 10), fmt = (a.format || 'ascii').toLowerCase();
    let values;
    if (fmt === 'ascii') {
      values = parseNumbersFromText(inner, null);
    } else {
      let bytes;
      const compressed = !!compressor;
      if (fmt === 'binary') {
        const text = inner.trim();
        bytes = await decodeBlocks(n => base64Decode(text, 0, Math.ceil(n * TYPE_SIZE[headerType] / 3) * 4),
          (hlen, dlen) => { const hchars = Math.ceil(hlen / 3) * 4; return base64Decode(text, hchars, Math.ceil(dlen / 3) * 4 + 4); },
          headerType, compressed, little, inflate, total => base64Decode(text, 0, Math.ceil(total / 3) * 4 + 4));
      } else if (fmt === 'appended') {
        const off = parseInt(a.offset || '0', 10);
        if (appEncoding === 'raw') {
          const base = appStart + off;
          bytes = await decodeBlocks(n => u8.subarray(base, base + n * TYPE_SIZE[headerType]),
            (hlen, dlen) => u8.subarray(base + hlen, base + hlen + dlen), headerType, compressed, little, inflate);
        } else {
          const text = getAppText();
          bytes = await decodeBlocks(n => base64Decode(text, off, Math.ceil(n * TYPE_SIZE[headerType] / 3) * 4),
            (hlen, dlen) => { const hchars = Math.ceil(hlen / 3) * 4; return base64Decode(text, off + hchars, Math.ceil(dlen / 3) * 4 + 4); },
            headerType, compressed, little, inflate, total => base64Decode(text, off, Math.ceil(total / 3) * 4 + 4));
        }
      } else throw new Error('unknown DataArray format ' + fmt);
      const count = Math.floor(bytes.byteLength / TYPE_SIZE[dtype]);
      values = readTyped(new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength), 0, count, dtype, little)[0];
    }
    return { name: a.Name || '', ncomp, data: values };
  }

  async function arraysIn(sectionRe) {
    log('section', String(sectionRe).slice(0, 30));
    const m = xmlText.match(sectionRe);
    log('section matched', !!m, m ? m[0].length : 0);
    if (!m) return { arrays: [], attrs: {} };
    const secAttrs = attrs(m[0].slice(0, m[0].indexOf('>') + 1));
    const body = m[1];
    const arrays = [];
    const re = /<DataArray\b([^>]*?)(\/>|>([\s\S]*?)<\/DataArray>)/g;
    let mm;
    while ((mm = re.exec(body)) !== null) arrays.push(await decodeArray(mm[1], mm[3] || ''));
    return { arrays, attrs: secAttrs };
  }

  const pts = await arraysIn(/<Points>([\s\S]*?)<\/Points>/);
  if (!pts.arrays.length) throw new Error('no <Points> found in the XML file');
  const points = pts.arrays[0].data;
  let triangles;
  if (/UnstructuredGrid/i.test(type)) {
    const cells = await arraysIn(/<Cells>([\s\S]*?)<\/Cells>/);
    const get = n => (cells.arrays.find(x => x.name === n) || {}).data;
    const conn = get('connectivity'), offs = get('offsets'), types = get('types');
    if (!conn || !offs) throw new Error('no cells in the unstructured grid');
    triangles = legacyCellsToTriangles(conn, [0, ...Array.from(offs)], types ? Array.from(types) : null, warnings);
  } else {
    const polys = await arraysIn(/<Polys>([\s\S]*?)<\/Polys>/);
    const get = n => (polys.arrays.find(x => x.name === n) || {}).data;
    let conn = get('connectivity'), offs = get('offsets');
    const strips = await arraysIn(/<Strips>([\s\S]*?)<\/Strips>/);
    let tri = conn && offs ? legacyCellsToTriangles(conn, [0, ...Array.from(offs)], null, warnings) : new Int32Array(0);
    const sget = n => (strips.arrays.find(x => x.name === n) || {}).data;
    if (sget('connectivity') && sget('offsets')) {
      const so = [0, ...Array.from(sget('offsets'))];
      const st = legacyCellsToTriangles(sget('connectivity'), so, new Array(so.length - 1).fill(6), warnings);
      tri = Int32Array.from([...tri, ...st]);
    }
    triangles = tri;
  }
  if (!triangles.length) throw new Error('no polygons found in the file');
  const pd = await arraysIn(/<PointData\b([^]*?)<\/PointData>/);
  // the regex above captured from the tag attributes; re-run properly to get the body
  const pdm = xmlText.match(/<PointData\b[^>]*>([\s\S]*?)<\/PointData>/);
  const pointArrays = [];
  let activeVectors = null;
  if (pdm) {
    const secAttrs = attrs((xmlText.match(/<PointData\b[^>]*>/) || [''])[0]);
    activeVectors = secAttrs.Vectors || null;
    const re = /<DataArray\b([^>]*?)(\/>|>([\s\S]*?)<\/DataArray>)/g;
    let mm;
    while ((mm = re.exec(pdm[1])) !== null) pointArrays.push(await decodeArray(mm[1], mm[3] || ''));
  }
  const cdm = xmlText.match(/<CellData\b[^>]*>([\s\S]*?)<\/CellData>/);
  const cellArrays = [];
  if (cdm) {
    const re = /<DataArray\b([^>]*?)(\/>|>([\s\S]*?)<\/DataArray>)/g;
    let mm;
    while ((mm = re.exec(cdm[1])) !== null) cellArrays.push(await decodeArray(mm[1], mm[3] || ''));
  }
  void pd;
  return { points, triangles, pointArrays, cellArrays, activeVectors, warnings, dataset: type, format: 'xml' };
}

// ---------------------------------------------------------------------------
export async function readVTK(input, opts = {}) {
  const u8 = input instanceof Uint8Array ? input : new Uint8Array(input);
  const head = bytesToString(u8, 0, Math.min(u8.length, 512));
  let result;
  if (/^\s*(<\?xml|<VTKFile)/i.test(head)) result = await readXML(u8, opts);
  else if (/^\s*#\s*vtk\s+DataFile/i.test(head)) result = readLegacy(u8, opts);
  else throw new Error('unrecognised file: expected a legacy .vtk file ("# vtk DataFile") or an XML .vtp/.vtu file');
  if (!(result.points instanceof Float64Array)) result.points = Float64Array.from(result.points);
  if (!(result.triangles instanceof Int32Array)) result.triangles = Int32Array.from(result.triangles);
  const N = result.points.length / 3;
  // cell vectors -> point vectors (area-unweighted average) when no point vectors exist
  if (!result.pointArrays.some(a => a.ncomp === 3) && result.cellArrays.some(a => a.ncomp === 3)) {
    const ca = result.cellArrays.find(a => a.ncomp === 3);
    const M = result.triangles.length / 3;
    if (ca.data.length === 3 * M) {
      const acc = new Float64Array(3 * N), cnt = new Float64Array(N);
      for (let t = 0; t < M; t++) for (let j = 0; j < 3; j++) { const v = result.triangles[3 * t + j]; cnt[v]++; for (let k = 0; k < 3; k++) acc[3 * v + k] += ca.data[3 * t + k]; }
      for (let v = 0; v < N; v++) if (cnt[v]) for (let k = 0; k < 3; k++) acc[3 * v + k] /= cnt[v];
      result.pointArrays.push({ name: ca.name, ncomp: 3, data: acc });
      result.warnings.push(`cell array '${ca.name}' averaged to the points`);
    }
  }
  return result;
}

/** Pick the vector array: by name, then the active vectors, then the first 3-component array. */
export function selectVectors(result, name) {
  const cands = result.pointArrays.filter(a => a.ncomp === 3 && a.data.length === result.points.length);
  if (!cands.length) throw new Error('the file has no 3-component point array (WSS vectors)');
  let a = name ? cands.find(x => x.name === name) : null;
  if (!a && result.activeVectors) a = cands.find(x => x.name === result.activeVectors);
  if (!a) a = cands[0];
  return { name: a.name, data: a.data, available: cands.map(x => x.name), fallback: name && a.name !== name };
}
