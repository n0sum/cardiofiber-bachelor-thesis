# `roney.py`: reimplementation of method 2 of Roney et al. 2019

Roney C.H. et al., *A technique for measuring anisotropy in atrial conduction
to estimate conduction velocity and atrial fibre direction*, Computers in
Biology and Medicine 104 (2019) 278–290, doi 10.1016/j.compbiomed.2018.10.019
(open access). Section numbers below refer to that paper.

This file lists where `roney.py` follows the paper and where it departs from
it, why, and what each departure does to the comparison with the per-face
estimator. Short version:

- the **errors** (fiber angle, mean speed e^μ, anisotropy ratio) are compared only on
  faces estimated by both methods, and on those faces `roney.py` does what the
  paper describes;
- the **non-estimable area** of Roney's method is partly caused by our choices
  (no interpolation of excluded regions, our residual threshold). It is
  reported split by cause, and only the part due to the method itself is
  attributed to it.

## Input and output

`roney.py` reads only `data/reference/SCENARIO/lat.ply` (the LAT maps), never
the ground truth, and for each triplet writes
`output/SCENARIO/roney/sub_i-j-k/estimate.ply` (the estimate in format D:
fiber in 3D projected on the plane of the face and normalised;
μ = ¼ ln(λ1·λ2), α = ¼ ln(λ1/λ2), with λ1 ≥ λ2 the eigenvalues of M = P⁻¹;
NaN on the non-estimable faces) and `estimate.csv` (`face_id,k_used`). The
errors against the ground truth are computed by `load()` of
`analysis/analysis.py`, with the same formulas as for the per-face estimator.

## The method in the paper (sections 2.3, 2.4, 2.6, 3.3)

For each element:

1. a 0.5 × 0.5 cm region centred around the measurement location (sizes chosen
   in initial parameter sensitivity testing), flattened with Dijkstra geodesic
   distances along the mesh edges and the multidimensional scaling of
   Zigelman et al.;
2. twenty recording locations chosen at random (all of them if there are
   fewer; mean 19.97, range 12–20);
3. for each map, a linear least-squares planar fit t ≈ α0 + α1·X + α2·Y of the
   LAT on those points, giving one planar velocity vector per map, assigned to
   the element closest to the centre of the region;
4. assuming that the velocity vectors of the different maps approximately form
   an ellipse, a centred ellipse (orientation θ, semi-axes CV_L and CV_T) is
   fitted to their tips "following the method of Ray and Srivastava",
   minimising the residuals measured along radii of the ellipse; with 3 maps
   all three parameters are estimated, with 2 maps either θ (from a fibre
   atlas) or the anisotropy ratio is fixed;
5. regions where a wavefront collision occurs in one of the maps are excluded
   with a threshold on the planar fit residual, and fibre directions and CV
   magnitudes there are interpolated from neighbouring values with Shepard
   interpolation (sections 2.6 and 4.7; section 3.3 notes that this
   "introduces error").

The paper does not report: the value of the residual threshold or how the
residual is normalised; whether the square region is taken in 3D or in the
flattened plane, and its orientation; the parameters of the Shepard
interpolation; the optimiser of the ellipse fit and its settings (the cited
reference, Ray and Srivastava, J. Struct. Geol. 30 (2008) 1593–1602, uses a
genetic algorithm).

## Differences and their effect

| step | paper | `roney.py` | effect on the comparison |
|---|---|---|---|
| region | 0.5 × 0.5 cm square around the element | geodesic disc of radius 2.5 mm around the face centroid | same diameter, 22% smaller area, no orientation to choose; minor |
| flattening, points, planar fit | Dijkstra + MDS, 20 random points, least squares | same; one random generator with fixed seed | none |
| collision threshold | value not reported | RMS residual / (\|g\| · radius) ≤ 0.14 | decides which faces are paired and part of Roney's lost area; see below |
| excluded regions | fibre and CV interpolated (Shepard) | face not estimable for that triplet | adds to Roney's lost area; does not touch the paired errors |
| ellipse fit | minimum radial residuals | exact centred conic through the 3 points | identical where an ellipse exists |
| no ellipse through the 3 points | not reported | face not estimable | the only non-arbitrary outcome; see below |
| number of maps | 1, 2 and 3 | 3 (all 20 triplets of the 6 sites) | comparison in the paper's best case |

### Residual threshold

The residual is normalised by |g| · radius, the time the wavefront takes to
cross the disc radius, so it is dimensionless and the same threshold holds for
slow and fast tissue. 0.14 is the median value of this residual for a
wavefront that turns by 60° along a line through the centre of the disc,
sampled with 20 random points: a map is rejected if its LATs in the disc are at
least as far from a plane as such a fold. As in the paper, one rejected map is
enough to exclude the face. The value 60° is our choice; the paper gives no
number to copy. A sensitivity study of this threshold has not been done yet.

### No Shepard interpolation

An interpolated fibre direction or CV is not a measurement at that face, and
the paper does not give the interpolation parameters. Faces where at least one
map of the triplet fails the planar fit are therefore counted as
non-estimable (`k_used < 3` in the output, see `experiments/README.md`). In
the method as published those faces would have an interpolated estimate.

### Exact conic instead of a genetic algorithm

x = g/|g|²: velocity vector of the wavefront (mm/ms), the quantity Roney fits.

A centred ellipse is xᵀ P x = 1, linear in the three entries of P, so three
velocity vectors give a 3 × 3 linear system. If its solution P is positive
definite, the ellipse passes exactly through the three tips: all radial
residuals are zero, which is the minimum of the paper's objective, and the
centred ellipse through three points is unique. Any optimiser of that
objective, genetic or not, converges to it; solving the system gives it
exactly and without randomness. Then M = P⁻¹.

If P is not positive definite the conic through the three tips is a
hyperbola. In that case the radial-residual objective has **no minimum among
ellipses**: the residual keeps decreasing as the ellipse is stretched along
the fiber towards a pair of parallel lines (CV_L → ∞). In a synthetic test
(exact velocity vectors, random wavefront directions, r = 1.6–2.4) the
optimiser always stopped at the upper bound of CV_L. The estimate would be
set by the bounds of the optimiser, which the paper does not report, so these
faces are counted as non-estimable (`k_used = 3`, NaN in `estimate.ply`).

This loss is a property of the method, not of the reimplementation. The
velocity vectors of plane wavefronts lie on the pedal curve of the wavefront
ellipse, not on the ellipse itself, and for r > √2 the pedal curve is not
convex. Three tips around the concave part admit no centred ellipse even with
exact data.

### Singular-value threshold

Some 3 × 3 systems are exactly singular (parallel wavefronts); without the
threshold `np.linalg.solve` stops with "Singular matrix".
The threshold has the same value as in the per-face estimator, but the rows
here are built from velocity vectors, (mm/ms)², instead of gradients,
(ms/mm)², so the same number does not mean the same thing in the two methods.
The paper itself notes that method 2 loses accuracy where the velocity vectors
of the maps are close to collinear (section 4.2).

## What the differences do to the results

The non-estimable area of Roney's method splits in two parts, by
`k_used` in `estimate.csv`: faces where a map of the triplet fails the planar
fit (`k_used < 3`), and faces where all three maps pass but the conic is not
an ellipse or the system is singular (`k_used = 3`, NaN in the estimate).

- **Errors.** Errors are compared only on faces estimated by both methods
  (`comparison.csv`). On those faces all three maps pass the planar fit, so
  interpolation plays no role, and the exact conic is the paper's
  radial-residual minimum. The only remaining differences there are the shape
  of the region (disc instead of square, same diameter) and the residual
  threshold, which decides which faces are paired. Paired faces are those
  where Roney has three clean maps and an ellipse.
- **Non-estimable area.** The first part is due to our choices (threshold
  and no interpolation) and is not a shortcoming of the method. The second
  part is due to the method. The total should not be presented as a
  comparison of coverage between the two methods.
- **Fidelity.** On our data the median fiber angle error of `roney.py` is
  in `comparison.csv` (`output/figures/`); the paper reports 10.7–12.1° with 3 maps (table 3). The data
  differ (monodomain simulations with the Courtemanche model there, eikonal
  here), and the published errors include the interpolated regions, so this
  is only an order-of-magnitude check.

## Open points

- Sensitivity of the comparison to `fit.max_residual` (for example 0.07, 0.14,
  0.28 and no threshold).
