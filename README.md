# WSSLCS Explorer

Fixed points, stable/unstable manifolds (WSS Lagrangian coherent structures)
and the direction field of a wall shear stress (WSS) vector field on a
triangulated surface, computed and rendered entirely in the browser.

**Live app: https://amir-cardiolab.github.io/WSS-topology/**

![WSS LCS Explorer](docs/WSSLCS_explorer2.png)

*The WSS LCS Explorer on the carotid artery demo: the surface coloured by
|WSS|, evenly spaced streamlines of the WSS direction field (white), the
fixed points (saddles as yellow octahedra, a sink as a sphere, foci with a
ring) and the unstable (blue, attracting WSS LCS) and stable (red, repelling
WSS LCS) manifolds, all computed in the browser.*

Open the page, load a surface file with a point vector array (or use the
demo data set that loads automatically), and explore the WSS topology:

* **fixed points** of the WSS field: sources, sinks, saddles, foci and
  centers, located from the Poincaré index of every triangle and the zero of
  the linear field inside it, with the eigenvalues/eigenvectors of the local
  Jacobian;
* **stable and unstable manifolds** of the saddles: the WSS trajectories
  leaving each saddle along its eigenvectors, integrated forward (unstable
  manifold, attracting WSS LCS, blue) and backward (stable manifold, repelling
  WSS LCS, red) until they reach another fixed point or leave the surface;
* **streamlines** of the WSS direction field, evenly spaced over the surface,
  which show how the manifolds organise the near-wall flow.

Results can be downloaded as legacy ASCII VTK files (readable by ParaView),
CSV and JSON, and the view as PNG (screen resolution, 3× with colour bar and
legend, or with a transparent background).  Nothing is uploaded: the file is
parsed and processed by JavaScript in your browser.

## Input files

* legacy `.vtk` (ASCII or binary, PolyData or UnstructuredGrid, old and
  VTK 5.1 cell layouts);
* XML `.vtp` / `.vtu` (ascii, inline base64, appended raw or base64 data,
  optionally zlib-compressed).

The surface must be triangulated (quads/polygons are fan-triangulated) and
carry a 3-component point array with the WSS vectors (cell vectors are
averaged to the points).  Any array can be selected from the drop-down.
Triangles with inconsistent winding are re-oriented automatically.

## Visualization

* **Curves that never fight the wall.** Manifolds and streamlines are lifted
  slightly off the surface along the local normal, the wall and the curves
  carry complementary depth offsets, and every manifold is drawn over a
  depth-dependent halo in the background tone, so curves stay crisp at any
  zoom and crossings read correctly.  Optionally, the parts of a curve hidden
  behind the wall are drawn as translucent dashes (*hidden parts as dashes*).
* **Faithful colours.** The |WSS| colour map (viridis, inferno, magma,
  cividis, turbo, greys) is sampled per pixel from a lookup table, with
  linear or logarithmic scale, automatic 2–98 % or manual range, units label
  and a colour bar with round tick values; no tone mapping; diffuse shading
  with calibrated lights so that the brightest colour equals the colour-map
  colour (options: soft sheen, unlit).
* **Readable symbols.** Fixed-point type is encoded by shape (saddle
  octahedron, source-like cube, sink-like sphere, ring for foci) and colour
  (classic per-type palette, or by stability in the manifold hues), with an
  outline and a minimum size on screen; saddle eigenvectors are drawn as
  short segments; a colour-blind-safe palette (Okabe–Ito) is available.
* **Options.** Screen-space strokes with pseudo-tube shading or shaded tubes;
  arrow glyphs; ambient occlusion and depth cueing (both off by default);
  light and dark themes; adjustable line width, halo, lift and marker size.
* **Interaction.** Orbit, pan and zoom to the cursor; view presets, an
  orientation gizmo and an orthographic projection; hover tooltips and
  selection linked with the result tables; double-click a fixed point (or a
  table row) to fly to it; *Copy settings link* produces a URL that
  reproduces every setting and the camera, and every control can be preset
  from URL parameters (e.g. `?theme=light&cmap=inferno&focus=3&dist=60`).

## Method and references

The tool implements the WSS topology analysis developed in the following
papers.  Please cite them if you use it:

1. Arzani, A., Gambaruto, A. M., Chen, G. and Shadden, S. C., *Wall shear
   stress exposure time: A Lagrangian measure of near-wall stagnation and
   concentration in cardiovascular flows*, Biomechanics and Modeling in
   Mechanobiology, 16(3), 787–803, 2017.
2. Arzani, A., Shadden, S. C., *Wall shear stress fixed points in
   cardiovascular fluid mechanics*, Journal of Biomechanics, 73, 145–152,
   2018.
3. Arzani, A., Gambaruto, A. M., Chen, G. and Shadden, S. C., *Lagrangian
   wall shear stress structures and near wall transport in high Schmidt
   aneurysmal flows*, Journal of Fluid Mechanics, 790, 158–172, 2016.

The stable and unstable manifolds of the WSS fixed points are the WSS
Lagrangian coherent structures (WSS LCS) introduced in [3], which organize
near-wall transport: near-wall trajectories accumulate along the unstable
manifolds (attracting WSS LCS) and separate along the stable manifolds
(repelling WSS LCS).  Reference [2] discusses the WSS fixed points themselves
(sources, sinks, saddles, foci) and their significance in cardiovascular
flows.  Reference [1] describes the computation used here (Section 2.4):
the fixed points of the (time-averaged) WSS vector field are located in the
triangles whose Poincaré index is non-trivial, the field is linearised
around them to obtain the Jacobian, its eigenvalues and eigenvectors, and
the saddle-type fixed points are perturbed along the eigenvector of the
positive (negative) eigenvalue and integrated forward (backward) in time to
trace the unstable (stable) manifold, until the trajectory reaches another
fixed point or leaves the domain.  Reference [1] also defines the WSS
exposure time, which is computed by the desktop Python version of this
software.

Implementation notes: vertex vectors are projected onto the tangent planes,
transported into the triangles with the discrete polar map of the vertex
fans, and interpolated linearly inside each triangle; manifolds are traced
with RK4 in arc length across the triangles (exact edge crossings, vertex
handling); streamlines are seeded with a minimum separation distance
(Jobard–Lefer style) and integrated in both directions.  This page is the
browser version of the Python/VTK package WSSLCS, which also computes WSS
trajectories, residence time and WSS exposure time.

## Files

| file | content |
|---|---|
| `index.html` | the page (three.js rendering, controls, tables, exports) |
| `wsslcs-core.js` | mesh (consistent winding, normals), field transport, fixed points, manifolds, evenly spaced streamlines, VTK/CSV writers |
| `vtk-reader.js` | VTK legacy / XML reader |
| `worker.js` | runs the analysis in a Web Worker |
| `colormaps.js` | colour map tables (generated with matplotlib) |
| `docs/WSSLCS_explorer2.png` | screenshot used in this README |
| `demo/Carotid_artery_TAWSS.vtk` | demo data: time-averaged WSS on a carotid artery surface (22 899 points, 45 794 triangles, 7 fixed points) |

The site is static: any web server (or GitHub Pages) can host these files.
three.js is loaded from the jsDelivr CDN and the fonts from Google Fonts.

## Unsteady WSS fields and surface tracer transport

This page analyses one (time-averaged) WSS snapshot.  For time-resolved WSS
fields and a more comprehensive analysis of near-wall transport, use our
ParaView plugin **[WSS-topology-unsteady-Paraview](https://github.com/amir-cardiolab/WSS-topology-unsteady-Paraview)**:
surface tracers advected by the unsteady WSS field (staggered release,
residence time and WSS exposure time), fixed points and manifolds of each
time step or of the time-averaged field, fixed points tracked in time, and
the time-series metrics TAWSS, OSI, RRT, time-averaged WSS divergence and
TSVI, all directly on the surfaces loaded in ParaView.
