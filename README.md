# WSSLCS Explorer

Fixed points and stable/unstable manifolds (WSS Lagrangian coherent
structures) of a wall shear stress (WSS) vector field on a triangulated
surface, computed entirely in the browser.

**Live app: https://amir-cardiolab.github.io/WSS-topology/**

Open the page, load a surface file with a point vector array (or use the
demo data set that loads automatically), and explore the WSS topology:

* **fixed points** of the WSS field: sources, sinks, saddles, foci and
  centers, located from the Poincaré index of every triangle and the zero of
  the linear field inside it, with the eigenvalues/eigenvectors of the local
  Jacobian;
* **stable and unstable manifolds** of the saddles: the WSS trajectories
  leaving each saddle along its eigenvectors, integrated forward (unstable
  manifold, attracting WSS LCS, blue) and backward (stable manifold, repelling
  WSS LCS, red) until they reach another fixed point or leave the surface.

Results can be downloaded as legacy ASCII VTK files (readable by ParaView),
CSV and JSON.  Nothing is uploaded: the file is parsed and processed by
JavaScript in your browser.

## Input files

* legacy `.vtk` (ASCII or binary, PolyData or UnstructuredGrid, old and
  VTK 5.1 cell layouts);
* XML `.vtp` / `.vtu` (ascii, inline base64, appended raw or base64 data,
  optionally zlib-compressed).

The surface must be triangulated (quads/polygons are fan-triangulated) and
carry a 3-component point array with the WSS vectors (cell vectors are
averaged to the points).  Any array can be selected from the drop-down.

## Method

The method is the one used in

> A. Arzani, A. M. Gambaruto, G. Chen, S. C. Shadden, *Wall shear stress
> exposure time: a Lagrangian measure of near-wall stagnation and
> concentration in cardiovascular flows*, Biomechanics and Modeling in
> Mechanobiology (2017).

Vertex vectors are projected onto the tangent planes, transported into the
triangles with the discrete polar map of the vertex fans, and interpolated
linearly inside each triangle.  Manifolds are traced with RK4 in arc length
across the triangles (exact edge crossings, vertex handling).

This page is the browser version of the Python/VTK package WSSLCS, which
also computes WSS trajectories, residence time and WSS exposure time.

## Files

| file | content |
|---|---|
| `index.html` | the page (three.js rendering, controls, tables, exports) |
| `wsslcs-core.js` | mesh, field transport, fixed points, manifolds, tracer |
| `vtk-reader.js` | VTK legacy / XML reader |
| `worker.js` | runs the analysis in a Web Worker |
| `demo/P3_TAWSS-1.vtk` | demo data: time-averaged WSS on an abdominal aortic aneurysm surface |

The site is static: any web server (or GitHub Pages) can host these files.
three.js is loaded from the jsDelivr CDN and the fonts from Google Fonts.
