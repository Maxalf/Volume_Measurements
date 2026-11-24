# Volume Measurement on Plane / Misura Volumi su piano

Web app to measure the volume of sphere-like objects lying on a plane, from `.asc` point clouds
(`X Y Z I J K` per line, `SCANLINE` lines ignored). Runs entirely in the browser: nothing to install,
the file is never uploaded anywhere. Interface in English and Italian.

## How it works
1. **Open** an `.asc` file (or drag it onto the page). The dominant plane is found with RANSAC and the cloud
   is rotated/translated so the plane is z = 0 with +Z pointing towards the scanner (sign taken from the file normals,
   i.e. the scan direction). If objects end up below the plane, use **Flip side**.
2. **Probe**: move/resize the red sphere (sliders, or double-click / `P` on a point, arrows, PgUp/PgDn, `[` `]`).
3. **Isolate and measure**: the local plane is refitted on a ring around the probe, plane points are removed,
   and the object is shown in a separate window with:
   - **closed surface**: visible surface triangulated as the scanner sees it, closed on the plane, V = Σ z̄·n_z·A
     (exact without undercuts, otherwise an upper bound);
   - **voxels**: the same region split into cubes (edge auto or user-set); inner cubes count fully, boundary cubes
     by their supersampled filled fraction, with a min–max range (only inner / inner + boundary cubes);
   - **sphere / ellipsoid fits** (volume above the plane, analytic) to reconstruct parts hidden from the scanner;
   - a **recommended estimate** choosing between them.
4. Export CSV of all measurements, object points (.asc, original coordinates), closed mesh (.stl), aligned cloud.

Units are those of the file (mm → mm³).

## Publish on GitHub Pages
Upload these files to a public repository, then *Settings → Pages → Deploy from a branch → main / (root)*.

## Third-party
[three.js](https://threejs.org) (MIT) and [delaunator](https://github.com/mapbox/delaunator) (ISC), bundled in `vendor/`.
