# Lymphatics Phase Difference Study — Session Notes

## Project Overview
CFD simulation of lymphatic vessels using the Lattice Boltzmann Method (LBM) with an Immersed Boundary/Level-Set Method (IBM/LS) for fluid-structure interaction. Simulates peristaltic pumping through a vessel with valve leaflets (lymphangion segments between valves). Solver is `Template/` (C++ / Visual Studio project `LBLS.vcxproj`, MATLAB scripts drive geometry generation and pre/post-processing). `Simulation/<name>/` folders are disposable run copies created by `run_sweep.m` / `setupSubfolders.m` — always fix bugs in `Template/`, never in a `Simulation/` copy.

**Note:** the "Current Status" content below (dated valve-extent hardcoding, speed-based motion, etc.) fully documents the *previous* session's state and has been superseded by this session's fixes — kept for history, see "This Session" below for what's now true.

---

## This Session (2026-07-22) — Fixes Applied, in Order

Each item was found by actually building and running the solver (`MSBuild` → `Template/x64/Release OpenMP/LBLS.exe`) and inspecting real output, not by code review alone.

### 1. Build didn't compile — `VESSEL_DIAMETER` undefined
- **Symptom**: `clMethod.h` referenced `c.ls.VESSEL_DIAMETER`, never defined anywhere.
- **Root cause**: leftover from an unfinished "position-based motion" feature; no `#define`, no struct member, no MATLAB writer.
- **Fix**: computed once at init in `clC::ini()` (`clVariables.h`) from the loaded geometry: `(max-min)/2` of vessel-wall (`lso==0`) node Y-coordinates. *(Later removed — see #2, its only use case was replaced by a per-node baseline.)*

### 2. Valve leaflet anchor nodes snapped into a circle on the first step
- **Symptom**: MATLAB plots of `lsc.txt` showed the leaflet root/anchor points (flat/planar at t=0) collapsing into a perfect circle within the first couple of simulation steps.
- **Root cause**: `clMethod.h`'s `solveVerlet1()`/`solveVerlet2()` (LEAFLETS branch) forced **every** `LS_CUSTOM_NODE_0/1`-flagged node onto the same absolute radius (`baseRadius = VESSEL_DIAMETER/2`) every step, ignoring that a crescent-shaped leaflet root isn't circular.
- **Fix**: target radius is now each node's own **initial** radius, read from `vls.C.c0(i)` (the untouched reference configuration), with the oscillation modulation added on top — preserves the original shape instead of erasing it. `VESSEL_DIAMETER` became unused after this and was removed (member, definition, computation site).

### 3. Valve extent (`VALVE_MIN_X`/`VALVE_MAX_X`) was wrong and the MATLAB pipeline that should've computed it was dead code
- **Symptom**: solver printed `Valve extent: X in [75.00, 300.00]` — stale hardcoded values from `ConstDef.h`.
- **Root cause**: `generateGeometry.m` *did* correctly compute the real extent, but never returned `cfg` to its caller (`function generateGeometry(cfg)`, no output arg), so `editConstDef.m`'s `LineModifier` calls for these two constants never actually fired — the whole propagation path was silently broken.
- **Fix**: compute `valveMinX`/`valveMaxX` once at init in C++ (`clC::ini()`, `clVariables.h`) directly from the loaded geometry: lowest X among `lso==1` nodes, lowest X among `lso==max(lso)` nodes. Removed the dead MATLAB computation (`generateGeometry.m`), the `LineModifier` calls (`editConstDef.m`), and the `#define`s (`ConstDef.h`).

### 4. Vessel was built at half its intended diameter
- **Symptom**: valve (leaflet) geometry visibly wider in Y than the vessel wall at t=0.
- **Root cause**: `GeometryCreatorAmir.m` passed `Vessel_diam/2` into `pipeFromGeometry.m`, which forwards that value unchanged into `generateStraightPipeGeometry.m` as `Diameter` — which *itself* does `radius = Diameter/2`. Net effect: actual vessel radius = `Vessel_diam/4`, i.e. the vessel's real diameter was half of `Vessel_diam`, while the leaflet geometry (built independently) correctly used the full `Vessel_diam`.
- **Fix**: `GeometryCreatorAmir.m:66` now passes `Vessel_diam` (not `Vessel_diam/2`). Verified numerically: with `Vessel_diam=20, Domain_offset=10`, the bug predicted `nY=nZ=30` (matched the old broken output exactly); the fix predicts and produces `nY=nZ=40`.

### 5. Peristaltic wavelength was a hardcoded constant, not tied to geometry
- **Ask**: wavelength should be `2 × intervalve (lymphangion) distance`, computed from the actual valve spacing.
- **Fix**: added `c.ls.waveLength`, computed once at init (`clVariables.h`) as `2 * (valveMaxX - valveMinX)`. Replaced all uses of the old `#define INTERVAL (150)` in the wave-motion formulas (`clMethod.h`) with `c.ls.waveLength`. Left `INTERVAL` itself defined in `ConstDef.h` since `RE_NUMBER`/`GAMMA_NUMBER` (an unrelated flow-forcing parameter) still needs it as a compile-time constant.
- **Bonus find**: `solveVerlet2()` had its own separate `OSCILLATING_VESSEL` block still using **hardcoded X-thresholds** (75/150/225/300/375 — an old fixed 3-valve layout) completely disconnected from the actual geometry. Since both `solveVerlet1()` and `solveVerlet2()` run every step and both write `vls.V(i)`, the stale block was silently overwriting the correct one. Rewrote it to match `solveVerlet1()`'s dynamic scheme.

### 6. Vessel wall radius grew without bound (quadratic in time)
- **Symptom**: tracked vessel Y-extent across saved snapshots in a live run — grew 20 → 380 over 23 dumps, with a **constant second difference (~1.48)**, the signature of quadratic growth.
- **Root cause**: plain vessel-wall nodes (`O(i)==0`) never carry the `LS_CUSTOM_NODE` flag (only leaflet nodes do), so they never went through the position-based fix from #2 — they were driven by the legacy `OSCILLATING_VESSEL` block, which prescribed a **radial velocity** (`speed = OSC_AMPLITUDE·sin(...)·sin(...)`) each step and let it integrate into position. Since `RAMP_PERIOD` is large relative to elapsed sim time, `sin(2π·Time/RAMP_PERIOD) ≈ 2π·Time/RAMP_PERIOD` early on, so the velocity itself grows ~linearly, and its integral (position) grows ~quadratically — unbounded drift, not a bounded oscillation.
- **Fix**: converted both copies of this block (`solveVerlet1()` and `solveVerlet2()`, `clMethod.h`) to the same position-based scheme as #2 — target radius computed from each node's own `c0` baseline plus bounded modulation, velocity driven as `(target-current)/dt`. Verified after the fix: radius stays in a physically reasonable range (min 4.0 / max 16.0 across a 577-frame run) instead of diverging.

### 7. Fluid was never actually being solved — root cause found and fixed
- **Symptom**: `lbro.txt`/`lbj.txt` (density/flux dumps) were **uniformly `-1`/`0` at every single frame**, across the entire run (checked frames 0, 16, 64, 200, 400, 576) — not near-zero, exactly sentinel values, meaning the LBM fluid solver's output was never populated despite `LB_SOLVER` being defined and the binary not being stale.
- **Root cause chain** (fully confirmed, traced end to end):
  1. `generateStraightPipeGeometry.m` builds the vessel as an *open* shell (no end caps) — not watertight.
  2. `pipeFromGeometry.m`'s ray-casting solid/fluid mask (`LBMl`, Möller–Trumbore test) against that open mesh comes back **uniformly 0** for the whole domain — a bug in that MATLAB ray-casting logic.
  3. `GeometryCreatorAmir.m` (line ~295) wrote this broken all-zero mask straight to `load/lbm.txt`.
  4. C++'s `clM::ini()` (`clVariables.h:3753`) does `if(load("load/lbm")==FALSE) inifill_map();` — since the file *existed* with the right size, `load()` succeeded, so the solver's own robust fallback (`inifill_map()`, a flood-fill that derives the mask correctly from the already-validated structural boundary mesh `lsbl.txt`) was **skipped**.
  5. `#define LB_SOLID 0` — so loading an all-zero mask marked **every node in the domain as solid**, leaving zero fluid nodes for the LBM solver to ever update. Density/flux stayed frozen at their pre-init sentinel values forever.
- **Fix** (option B — chosen over rewriting the fragile ray-casting): `GeometryCreatorAmir.m` no longer writes `load/lbm.txt` when `Vessel_flag==1 && No_solid_flag~=1`. Absent the file, `clM::ini()` naturally falls back to its own correct `inifill_map()`. Also deleted the two stale buggy copies already sitting in `Template/load/` and `Template/geometry creator/load/` (otherwise `generateGeometry.m`'s `copyfile` sync would've resurrected them).
- **Verified working**: fresh run's `lbro.txt`/`lbj.txt` now show 231,925 / 670,132 distinct values per 240,000/720,000-value snapshot by dump #10 (up from 1 distinct value at dump #1, i.e. flow visibly developing) — finite, no NaN/Inf, physically small density (±0.004) and flux (±0.0047) consistent with the zero-pressure-forcing (`DELTA_RO=0`) case. The dumped mask is now uniformly `1` (not a new bug — this solver couples fluid to structure via immersed-boundary forcing, not bounce-back walls, so the whole domain legitimately stays "fluid").

### 8. Second valve's leaflets were structurally undriven (not an FSI/fluid bug — a geometry-vs-region mismatch)
- **Symptom**: using the now-fixed `LivePostProcessor.m`, valve 2's leaflet (`lso==3,4`) showed uniformly near-zero velocity (dark, flat) while valve 1's (`lso==1,2`) showed strong, varying velocity — confirmed via direct `lsv.txt`/`lbj.txt` inspection (~20-40x smaller magnitude, steady across frames 20-177, not a transient still catching up).
- **Root cause, confirmed against actual geometry data**: `valveMinX`/`valveMaxX` mark where each valve's leaflet *starts* (lowest X of `lso==1` and of `lso==max(lso)` respectively). Valve 1's leaflet (X ∈ [25, 48.9]) sits entirely *inside* `[valveMinX, valveMaxX) = [25,125)`, so it received real wave modulation. Valve 2's leaflet (X ∈ [125, 148.9]) starts exactly *at* `valveMaxX` and extends *beyond* it — every one of its nodes fell into the formula's final `else` branch, forcing `radiusModulation = 0` for the entire leaflet, for all time. Ruled out boundary conditions (symmetric, confirmed via `bcPressLBx0`/`x1` and `RoX0`/`RoX1`) and the wave formula's spatial symmetry before landing on this.
- **Fix** (per user direction — remove the X-range gating entirely, let the geometry's own static flag do that job): in all four copies of the modulation logic (`solveVerlet1`/`solveVerlet2` × LEAFLETS branch / `OSCILLATING_VESSEL` branch, `clMethod.h`), removed `isOutsideValveRegion`/`x<valveMinX`/`x>=valveMaxX` clamping entirely. `radiusModulation` is now `OSC_AMPLITUDE·sin(2π·localX/waveLength)·sin(2π·Time/RAMP_PERIOD)` unconditionally, zeroed only when `T(i) & LS_STATIC_NODE` — the flag `GeometryCreatorAmir.m`'s `entrance_static_flag` block already sets on vessel-wall (`lso==0`) nodes outside the two valves' footprint (confirmed enabled: `buildConfig.m`'s `entranceStatic=true`). This flag was defined at geometry time but never actually read anywhere in `clMethod.h` before this fix. Leaflet nodes never get this flag (the static-marking loop only touches `O==1`/vessel wall), so both valves' leaflets now get real, symmetric modulation based on their own position — valve 2's now ranges over roughly the second half of one wavelength (opposite phase from valve 1), which is physically the expected traveling-wave behavior between two valves.
- Also fixed in passing: `ConstDef.h`'s `LB_P_X0` (inlet pressure) was hardcoded to `1./3.`, ignoring `DELTA_RO_WALL` (the line that included it was commented out), while `LB_P_X1` correctly used it. Didn't affect this run (`DELTA_RO=0` makes both forms equal) but would break any nonzero-pressure-difference sweep run. Restored to `(1./3. + DELTA_RO_WALL)`.
- Also fixed: `LivePostProcessor.m`'s `leaflet_mask = lso >= 2` silently excluded `lso==1` (valve 1's first leaflet) — changed to `lso ~= 0` so both valves show all their leaflets.

### 9. `LivePostProcessor.m` showed nothing
- **First bug**: `lbro.txt`/`lbj.txt` aren't flat `nX*nY*nZ` vectors — `clArray3D::save2file()` (C++) writes `nX` blocks (one per X) of `nY` rows × `nZ` (or `nZ*3` for vectors) columns, blank-line-separated between blocks. A plain `load()`+`reshape()` either errors (dimension mismatch) or silently scrambles the data (MATLAB `reshape` is column-major, the save format is row-major). **Fix**: added `Template/loadMultiBlock3D.m` and `Template/loadMultiBlock3Dv3.m` — parse the actual block format via `textscan` (which skips blank lines) + `reshape`/`permute` into the correct axis order. Same parsing convention already existed independently in `loadfile3d.m`/`loadfile3d_xyz.m` (manual triple-loop reindexing) — cross-checked, matches.
- **Second bug (this round)**: even after the loader fix, the script could still silently do nothing for several reasons that weren't diagnosable from the UI alone — wrong `pwd` when launched (silently polls an empty/nonexistent `results/` forever, zero output), a live race between a `results/N/` folder being created and all its files finishing their atomic rename into it (frame gets permanently skipped, not retried), and `streamline` being a fragile choice for a mostly-tiny velocity field. **Fix**: rewrote `Template/LivePostProcessor.m` — validates `load/`/`results/` exist at startup with a clear error, prints periodic heartbeat status while idle (so "nothing new yet" is visibly different from "the script died"), only advances past a frame *after* it renders successfully (retries transient failures instead of skipping), replaced `streamline` with `quiver` for the midplane velocity, added explicit/adaptive color limits so a small-magnitude density field doesn't wash out to a flat color, and downsampled the velocity quiver grid.

### 10. `LivePostProcessor.m` right-hand subplot converted to 3D
- **Ask**: view the valve/leaflet nodes in 3D (X/Y/Z) instead of the flattened X/Y-only 2D scatter.
- **Change**: `Template/LivePostProcessor.m`'s second subplot (`ax2`) now uses `scatter3`/`quiver3` (X, Y, Z + 3D velocity vectors) instead of `scatter`/`quiver`. Still colored/scaled by velocity magnitude, same `lso~=0` leaflet mask from fix #8.
- **View angle persistence**: `view(ax2, 3)`, `DataAspectRatio`, and axis labels are now set **once at startup**, not inside the per-frame loop. `cla(ax2)` (used every update) only clears plotted objects, not axis/camera properties, so if you manually rotate the 3D view with the mouse, that angle now survives across live updates instead of snapping back to the default every 0.5s.

### 11. `LivePostProcessor.m` merged into a single combined 3D view
- **Ask**: combine the two side-by-side subplots into one figure.
- **Change**: replaced the two `subplot(1,2,·)` axes with **one** 3D axes showing everything together:
  - Y-midplane density, drawn as a semi-transparent colored `surf` slice positioned at the actual midplane Y-coordinate (`midY_coord = midY-1`, converting the 1-indexed MATLAB slice index to the same 0-indexed physical coordinate `lsc` uses) — this is what drives the axes colormap (`jet`) and colorbar.
  - Fluid velocity on that slice: black `quiver3` arrows (now full 3D — added the Y-component, which the old 2D version silently dropped).
  - Leaflet/valve nodes: magenta `scatter3`, with **marker size** (not color) encoding velocity magnitude — color was kept reserved for density so structure markers can never visually blend into the jet-colormapped fluid slice (two independent color scales in one axes isn't natively supported without a fragile dual-colorbar hack, so size was used as the second channel instead) — plus matching magenta `quiver3` velocity arrows (arrow length still separately scaled to stay visible regardless of true magnitude, as before).
  - View angle, `DataAspectRatio`, axis labels, and the `jet` colormap are set once at startup (same reasoning as #10 — `cla` doesn't touch those, so manual rotation persists across live updates).

### 12. `LivePostProcessor.m` now shows the vessel wall itself
- **Ask**: see where the vessel wall (`lso==0`) is and which way it's moving, at what velocity.
- **Change**: added a third layer to the combined 3D view, styled distinctly from both the fluid slice (`jet` colormap) and the leaflets (magenta):
  - **Shadow**: `scatter3` of wall nodes (downsampled to ~1200 points for a clean tube outline), small gray dots (`[0.55 0.55 0.6]`), `MarkerFaceAlpha=0.25` — dense enough to read as a translucent tube silhouette, transparent enough to not obscure the fluid slice or leaflets behind it.
  - **Velocity**: `quiver3` arrows in steel blue (`[0.2 0.45 0.75]`), downsampled much more sparsely (~150 arrows) than the shadow — enough to see the peristaltic wave's spatial pattern (direction reverses between the compressing/expanding halves of the cycle) without turning the tube into a solid mass of arrows. Arrow length is independently scaled to the wall's own max velocity in view (same scheme as the leaflet arrows, kept separate since wall and leaflet speeds are generally very different magnitudes).
  - Added a legend (`findobj(ax,'-regexp','DisplayName','.+')`, rebuilt every frame since `cla()` deletes the underlying handles each update) so all four now-present categories — fluid velocity, vessel wall shadow, vessel wall velocity, leaflet/valve — are distinguishable at a glance instead of relying on remembering the color code.

### 13. Fix #8's X-range removal got partially reverted — valve 2 was silently undriven again
- **Investigation**: asked to check whether valve 2's leaflet is "even flexible" using a live post-fix run (using the #12 visualization). Isolated `LS_CUSTOM_NODE_0`-flagged anchor nodes specifically (via `lst.txt`, checking bit 256) and tracked one anchor node's *actual radius* (not velocity — see why below) over 44 frames at a matched, non-degenerate X on each valve (X=41 for valve 1, X=141 for valve 2, both away from the sine's zero-crossings). Valve 1's anchor smoothly displaced (`dr`: 0 → 1.49 over 44 frames, a slow ramp consistent with `RAMP_PERIOD`). Valve 2's anchor stayed at **exactly** `dr=0.0000` every single frame — not smaller, literally undriven.
- **Side-finding, not a bug**: anchor-node *velocity* (`lsv.txt`) reads near-zero for both valves regardless — this is an artifact of the dual-half-step scheme: `solveVerlet1` moves the anchor's position exactly to that step's target, then `solveVerlet2` (same, not-yet-advanced `c.t.Time`) recomputes the same target against a position that's already there, so its saved velocity is ≈0 by construction. Position/displacement is the correct signal to check for anchor nodes, not velocity.
- **Root cause**: `Template/source/clMethod.h`'s two LEAFLETS branches (`solveVerlet1`/`solveVerlet2`, the ones that actually drive valve anchors) had reverted to `if(x >= valveMinX && x < valveMaxX) { modulation = ... }` — i.e. fix #8's X-range removal had been undone in these two spots specifically (the original explanatory comment from fix #8 was still sitting directly above the reverted code, now contradicting it). The two `OSCILLATING_VESSEL` branches (vessel wall motion) were untouched and still correct. Cause of the revert unconfirmed — possibly a manual edit made directly against a `Simulation/` copy's source and not reconciled back.
- **Fix**: removed the reintroduced X-range check from both LEAFLETS branches again, back to the unconditional formula gated only by `LS_STATIC_NODE`, matching the two `OSCILLATING_VESSEL` branches. Rebuilt clean.
- **Action needed**: the currently-running `Simulation/pressureDiff_0` process was compiled from the reverted code (confirmed: its `source/clMethod.h` has the same regression) — its results still show valve 2 undriven. It needs a fresh run from the now-corrected `Template/` to actually test whether valve 2 is flexible.

### 14. New batch flow-field animation exporter (`Template/FlowAnimationExporter.m`)
- **Ask**: "I want an animation of the flow field, I had something for this — what was it?"
- **What they had**: `AnimePlotterMainFlow_FIXED.m` + `computeFrameDataFlow.m` — a per-frame PNG-then-GIF exporter with selectable modes (`Velocity_contour`, `Pressure_contour`, `Vorticity_contour`, `Stress_contour`, etc.). Found it in a state that wouldn't produce a useful flow animation as-is: `Stress_contour=1` active (not a flow mode), the frame loop was `for subfolderIdx = 15%200:nFrames` (the `%` makes everything after it a comment, so it only rendered frame 15), hardcoded `xlim`/`ylim`/`zlim` and normalization constants (e.g. `A = (2*pi*10*3*5*25)^(1/3)/100000`) from an older/different domain size, and — the real blocker — `computeFrameDataFlow.m` has its own hardcoded peristaltic wave formula (`r = 10 - 0.5*3*sin(...)`, valve positions assumed at X≈50/175) for masking the vessel interior, completely disconnected from the actual dynamic `valveMinX`/`valveMaxX`/`waveLength` this session wired up in `clMethod.h`. Patching it would mean fixing all of the above plus re-deriving that mask formula.
- **Decision** (user's choice, offered as an alternative to patching the old script): build a new exporter reusing `LivePostProcessor.m`'s already-correct, already-tested combined 3D rendering (density slice + fluid velocity + vessel wall shadow/velocity + leaflet/valve nodes/velocity) instead of the old per-mode contour system.
- **New file**: `Template/FlowAnimationExporter.m`. Same rendering block as `LivePostProcessor.m`, but instead of polling for the latest frame in a `while true` loop, it enumerates every saved `results/N` folder up front (configurable `frameStride` to skip frames for speed), renders each with an *offscreen* figure (`Visible','off'`, faster for batch work), saves a PNG per frame via `exportgraphics`, then stitches all PNGs into `FlowAnimation.gif` (same `rgb2ind`+`imwrite` approach as the old script's GIF step, which was fine and got reused as-is). Run from the simulation folder, same as `LivePostProcessor.m`.
- **Parallelized** (follow-up ask: "as fast as physically able"): rendering is embarrassingly parallel (each frame only reads its own `results/N` and writes its own PNG), so the per-frame loop is now `parfor` across a pool sized to `feature('numcores')` (physical cores — avoids hyperthreading oversubscription for this CPU-bound rendering work). Note the *old* `AnimePlotterMainFlow_FIXED.m` called `parpool()` at the top too but only ever used a plain `for` after that, so it was never actually parallel despite starting a pool.
  - Since each `parfor` iteration runs on a separate MATLAB worker *process*, figures/axes can't be shared or reused across iterations the way `LivePostProcessor.m` reuses one via `cla()` — every iteration now builds its own offscreen figure from scratch and closes it when done.
  - `savedFiles`/`renderOK` are pre-sized (`cell(nFrames,1)` / `false(nFrames,1)`) and written by fixed index (`savedFiles{fi}=...`) rather than grown with `end+1`, since `parfor` requires sliced, fixed-size output assignment — a failed frame just leaves its slot empty/false instead of breaking indexing for every frame after it.
  - GIF stitching stays a plain sequential loop after the `parfor` completes — one output file, must be written in frame order, and it's I/O-bound (fast) next to rendering anyway, so there's nothing to gain from parallelizing it.

### 15. Connected this codebase to the actual paper it implements
- **Context**: user shared the full LaTeX source of "Pumping dynamics of a single lymphangion" (Poorghani, Karimi, Dixon, Alexeev) and asked how this solver relates to it.
- **Finding**: this *is* the paper's computational framework, not just something similar — confirmed via multiple concrete, verifiable matches:
  - `LBLS.cpp`'s header (`// (c) Alexander Alexeev, 2006`) — Alexeev is the paper's senior author, and the paper's own Methodology section says the model "builds upon a prior computational framework developed by our group."
  - `clVariablesLS` (which this session's notes had been calling "level-set") is actually **Lattice Spring Method (LSM)** per the paper — a correction to earlier assumptions in these notes. `clVariablesLB` = the paper's LBM (D3Q19, single relaxation time, Zou-He pressure BCs = `bcPressLBx0()`/`bcPressLBx1()`).
  - Exact matching constants: paper states $R=10$, $\rho=1$, $\mu=1/6$, $40\times40$ cross-section mesh (their "M3", validated and used for all runs), valve length $\ell=1.25$ — all match `Vessel_diam=20` ($\to R=10$), `MU_NUMBER=1/6`, `nY=nZ=40`, `Leaflet_length=1.25` exactly.
  - **Smoking gun**: `AnimePlotterMainFlow_FIXED.m`'s `Velocity_contour` block has the annotation string *"Supplementary Animation 1... $\mathcal{L}=5$, $\Gamma=0.14$, $K=0.065$, $\Delta P=0$"* — a verbatim match to the paper's Figure 3 caption. That script is what generated Figure 3 / Supp. Animations 1-2.
  - `DELTA_RO` ↔ paper's $\Delta P$; `Simulation/pressureDiff_0` (`DELTA_RO=0`) is the paper's $\Delta P=0$ baseline used in most figures.
- **Caveat noted**: the `nX=150` test domain used throughout this session's debugging is smaller than the paper's $\mathcal{L}=5$ baseline (`nX=225`, mesh M3) — a scaled-down config for fast iteration, not literally reproducing a paper figure until this session's fix (#16 below).

### 16. `buildConfig.m` rebuilt around the paper's non-dimensional parameters
- **Ask**: enter the paper's non-dimensional numbers directly ($\mathcal{L}$, $\Gamma$, $K$, $\Delta P$) instead of raw LBM values, and have the config generate the matching raw simulation automatically. Specifically requested for the Figure 3 case.
- **Real bug found and fixed along the way**: `GeometryCreatorAmir.m` had `EIv=1.0; EEv=1.0;` hardcoded immediately after the function signature, silently discarding whatever valve stiffness was passed in from the caller. This meant **valve bending stiffness could never actually be controlled** through `buildConfig.m`/`cfg`, regardless of what anyone set there — matching the paper's $K$ was impossible until this was removed. `generateGeometry.m` updated to actually pass `cfg.geometry.valveBendStiff`/`valveExtStiff` through (wall stiffness `EEW` intentionally left fixed — wall motion is fully prescribed/position-based, not stiffness-driven, so it isn't tied to any of the paper's non-dimensional groups).
- **Derivations implemented** (all in `buildConfig.m`, from five paper inputs `L_nd`($\mathcal{L}$), `Gamma`($\Gamma$), `K`, `dP`($\Delta P$), plus the paper's fixed `A_nd`=0.15 and `ell_nd`=1.25):
  - Period: Womersley number $\Gamma=R\sqrt{2\pi/(\nu T)} \Rightarrow T=2\pi R^2/(\nu\Gamma^2)$, then `cfg.sim.period = 2T` (accounts for `RAMP_PERIOD = CYCLE_PERIOD/RAMP_PER_CYCLE(=2)` in `ConstDef.h`, so `RAMP_PERIOD` ends up exactly `T`).
  - Amplitude: $A=2R\mathcal{A}$ (raw), `cfg.sim.oscAmplitude = A/2` (accounts for `amplitudeScaleFactor=2` in `editConstDef.m`, which doubles it back to `A` when written as `OSC_AMPLITUDE`).
  - Valve spacing: `valveMaxX-valveMinX` (the geometry's actual valve-to-valve span) = `spacingA * vesselDiameter * leafletLength`, and we want that to equal $L=2R\mathcal{L}$, so `spacingA = L_nd/ell_nd`.
  - Valve bending stiffness: $K=D_b/(2\pi R^3 P_0) \Rightarrow D_b$, via $\Delta V=2\pi RAL$, $Q_0=2\Delta V/T$, $P_0=8\mu LQ_0/(\pi R^4)$. Then $D_b \to$ raw `EIv`: worked backward through `leafletCrescent.m`'s `DEw = EIw/((1-\nu_p^2)(NZr-1)L_{eq}\sin(\pi/3))` and confirmed algebraically that `C0w` (the actual bending spring constant used) equals the paper's $k_b=4D_b/(3\sqrt3)$ exactly when `DEw=D_b` — so `EIw` is just `D_b` scaled by the mesh's own geometric factors ($NZr$, $L_{eq}$, computed the same way `GeometryCreatorAmir.m` does, from `Vessel_diam` and the fixed `Leq_min=1.5`, $\nu_p=1/3$).
  - Stretching stiffness (`EEv`): the paper states this has a minor effect on valve deformation and gives no matching non-dimensional group, so it's kept equal to the derived `EIv` (same convention the deprecated `SimulationCreator.m` used).
- **Verified independently**: re-implemented the exact same formulas outside MATLAB (bash/awk) and got identical numbers to the `.m` file's logic — for the Figure 3 defaults ($\mathcal{L}=5$, $\Gamma=0.14$, $K=0.065$): $T\approx192{,}342$, `spacingA=4`, `EIv\approx6.04`, `cfg.sim.period\approx384{,}685`, `cfg.sim.oscAmplitude=1.5`.
- **Left as free/non-paper-derived parameters** (documented in comments): `spacingB` (immobilized end-extension length — paper doesn't specify one) and vessel wall stiffness (not stiffness-driven).
- **Not run yet** — no MATLAB available in this session to execute it; the derivation is verified arithmetically but not against an actual simulation result.

### 17. PACE (Georgia Tech HPC) deployment
- **Ask**: run the paper-matched config on PACE.
- Confirmed the codebase already has an established PACE workflow: `Template/Makefile` (Intel `icpc`, `-qopenmp`) and `Template/run.sbatch` (SLURM, account `gts-aalexeev3`, email `apoorghani3@gatech.edu` already filled in, `module load intel/20.0.4`).
- Checked the C++ source for Windows-only calls that would break a Linux/icpc build (`_mkdir`, path separators) — everything is already properly guarded with `#if defined(_WIN32)`/`#else` (e.g. `clOutput.h`/`clOutput.cpp`'s `MakeDir`/`ini()` correctly fall back to POSIX `mkdir()`), so no portability fixes were needed.
- **No direct PACE access from this session** (no SSH/credentials) — provided the transfer/build/submit command sequence (`rsync` → `module load` + `make` → `sbatch run.sbatch`) instead of running it. Flagged three things worth checking before submitting: the leftover job name (`-JOnlyContrTest`), whether the walltime (`-t1-6`) is enough for ~385,000 timesteps at PACE's actual steps/second (unknown without a first run), and whether `OMP_NUM_THREADS=5` is using the full node or should be increased.

---

## This Session (2026-07-27 – 2026-07-30) — PACE Deployment, Visualization Overhaul, Reviewer-Response Figure

### 18. Actually got the paper-matched config running on PACE
- **Account fix**: `run.sbatch`'s `-Agts-aalexeev3` was wrong for this allocation. First tried `-Aaalexeev3-afrench` (rejected — "Invalid account"), then confirmed via `sacctmgr -p show assoc user=apoorghani3` that the real account is `gts-aalexeev3-afrench` (needs the `gts-` prefix after all).
- **Compiler module fix**: `module load intel/20.0.4` no longer exists on this PACE node (RHEL9, Lmod). `module avail intel` showed `intel/2021.9.0` is the version that still provides the classic `icpc` (needed since `Makefile` targets `icpc`/`-qopenmp`, not the newer `icpx` in `intel-oneapi-compilers`).
- **Link error fix**: `-fast` (in `Makefile`'s `CCFLAGS`/`LINKFLAGS`) implies `-static` on Linux, and this RHEL9 node has no static `libstdc++.a` installed → `cannot find -lstdc++` at link time. Fixed by replacing `-fast` with its non-static equivalent: `-ipo -O3 -no-prec-div -fp-model fast=2` (same optimization level, dynamically linked).
- **Job name / walltime / threads**: renamed the stale `-JOnlyContrTest` to something descriptive per run; walltime bumped `1-6` → `2-0` (still turned out insufficient — see #21); replaced the hardcoded `OMP_NUM_THREADS=5` with `--cpus-per-task=24` + `export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK` so thread count always matches whatever's actually requested.
- **Result**: `lbls_run` (job 11518716, the `Simulation/pressureDiff_0` / L_nd=5 / paper Fig. 3 config) submitted and confirmed genuinely computing (CPU time climbing in lockstep with wall-clock across repeated `ps` checks, ~860-960% CPU / ~9 of the 24 requested cores actually utilized — worth revisiting if throughput matters, OpenMP isn't scaling to the full core count).

### 19. Set up direct SSH access (key-based, no more copy-pasting terminal output)
- Generated a dedicated passphrase-less ed25519 keypair (`~/.ssh/pace_ed25519`) specifically for this automation — user appended the public key to PACE's `~/.ssh/authorized_keys` via their own already-Duo-authenticated session (one manual step, unavoidable — Duo/passwords can't be scripted).
- Added an `~/.ssh/config` `Host pace` alias. From here on, `ssh pace '<command>'` runs directly without relaying through the user's terminal.
- **Caveat documented for the user**: this key has full, unrestricted shell access to their PACE account — revocable anytime by deleting the line from `authorized_keys`.

### 20. Generated a second run: shorter vessel (L_nd=4)
- **Ask**: sweep a shorter lymphangion length than the paper's tested minimum (L_nd=5).
- `Template/buildConfig.m`: `paper.L_nd` 5→4 (documented as an extrapolation below the paper's validated 5-9 range), `cfg.simulationFolderName` changed to `Simulation_shortVessel_L4` (a distinct top-level folder, specifically to avoid `run_sweep.m` silently overwriting the existing `Simulation/pressureDiff_0` that mirrors the live PACE job).
- Ran `run_sweep.m` → generated `Simulation_shortVessel_L4/pressureDiff_0` locally (`nX=155`, vs. the wider L_nd=5 config — confirms the geometry actually got shorter). Shipped to PACE via `tar | ssh | tar` (rsync isn't installed in this shell environment — used as a one-shot substitute), built cleanly with the now-fixed `Makefile`/module, submitted as job `lbls_L4` (job 11518778).

### 21. Both PACE jobs hit their walltime limit before finishing
- `sacct` showed both `lbls_run` and `lbls_L4` ended in `State: TIMEOUT` at almost exactly 2 days elapsed — killed by the `-t2-0` walltime cap set in #18, not because they completed the ~385,000-step run. **Not yet resubmitted with a longer walltime or a checkpoint/restart** — see Known Issues.

### 22. Fluid visualization: mask out the domain padding around the vessel
- **Ask**: zero out density/velocity outside the vessel's actual wall bounds, so only true lumen fluid is visualized (matching the paper's own masked contour style).
- New `Template/computeVesselRadiusProfile.m`: derives per-X wall center/radius **directly from that frame's actual wall-node positions** (`lso==0`), not by re-deriving the C++ peristaltic wave formula in MATLAB — reimplementing that formula independently is exactly what caused the old `computeFrameDataFlow.m`'s hardcoded, now-stale `r = 10-0.5*3*sin(...)` mask (tied to a since-changed geometry). `interp1(...,'extrap')` fills every lattice-X column, since the wall mesh's own ring spacing (every ~2 units) is coarser than the fluid lattice's 1-unit spacing — an earlier version that required an *exact* `round(X)` match left gaps, which rendered as leftover un-masked "always inside" stripes.
- New `Template/computeVesselRadiusMask.m`: builds the `nX x nZ` inside/outside mask from that profile for the Y-midplane slice.
- Wired into `LivePostProcessor.m`/`FlowAnimationExporter.m`: masked-out fluid points set to **NaN** (not 0) — `surf`/`quiver3` skip NaN-valued faces/arrows entirely, so padding renders as nothing instead of a flat zero-colored region.

### 23. Rendering artifact: striping turned out to be two separate bugs
- First render attempt still showed cyan stripes outside the vessel despite the mask being numerically correct (verified: mask row-sums matched the expected radius-based band width, e.g. 20/40 at the wide ends, 14/40 at the narrow valve constriction — no bug there).
- **Real cause**: `FaceAlpha < 1` (had been 0.85) on a `surf` at a steep oblique 3D angle is a known MATLAB alpha-blending depth-sort artifact, made worse once the translucent wall surface (#24) was added on top. Fixed with two changes: density `surf` set to fully opaque (`FaceAlpha=1.0`), and `ax.SortMethod = 'childorder'` (composites strictly in plotting order instead of per-pixel depth-sorting transparent objects) set once at axes creation in both scripts.

### 24. Point-cloud → real triangulated surfaces for wall/valve, with valve trimming
- **Ask**: stop rendering the vessel wall and leaflets as `scatter3` point clouds; use actual surfaces, and trim the portion of the valve that pokes outside the vessel's current bounds. Referenced the old (now largely superseded) `AnimePlotterMainFlow.m`/`computeFrameDataFlow.m` for the general approach (triangle-list → patch, per-node radial trim).
- New `Template/buildTrimmedStructureSurfaces.m`: builds `patch`-ready `Vertices`/`Faces` from `load/lsbl.txt` (static Nx3, 0-indexed triangle connectivity — confirmed 7684 triangles, matches wall+valve triangle counts summing correctly with no loss) split into wall (`lso==0`) vs. valve triangles via an explicit boolean mask (not an assumed contiguous sort order, unlike the old script). Leaflet node positions used for the surface are radially clamped to `computeVesselRadiusProfile`'s current wall radius *before* triangulating (a rendering-only adjustment — real node velocities/positions used elsewhere are untouched), so valve tips that have flexed past the wall don't visually poke through it.
- Wired into both visualization scripts, replacing the old `scatter3` wall-shadow + leaflet-dot blocks.

### 25. Wall/valve display toggles + vectors restricted to the midplane
- **Ask**: disable the wall surface, disable the valve surface, only show velocity vectors on the Y-midplane (not the separate wall-velocity/leaflet-velocity arrow blocks that existed before).
- Added `showWall`/`showValve` config flags (both `false` by default) gating the patch calls in both scripts; entirely removed the wall-velocity and leaflet-velocity `quiver3` blocks (not just gated — an absolute simplification per the ask, independent of the wall/valve toggles) so only the fluid velocity quiver on the midplane remains.

### 26. Frame-save format: `.fig` then back to `.png`
- First changed `FlowAnimationExporter.m` to `savefig(...)` (MATLAB `.fig`, one per frame) per an explicit ask — necessarily dropped the old GIF-stitching step in the process, since `imread`/`rgb2ind` can't operate on `.fig` files.
- Then reverted to `exportgraphics(..., 'Resolution', 150)` → `.png` per a follow-up ask ("now that i have the files, edit the code so it saves them as regular image files") — GIF stitching was **not** restored (wasn't asked for); output is per-frame PNGs only.

### 27. Reviewer-response figure: inlet valve delayed closure / transient backflow
- **Context**: Physics of Fluids reviewer (Referee #3, MS POF26-AR-02859R) asked for a dedicated visualization of the delayed-valve-closure/transient-backflow phase (~τ=0.8-0.9 per the paper's own Fig. 4a for the inlet valve) — not covered by the existing Figure 3's four snapshots (τ=0, 0.25, 0.5, 0.75).
- **Data-driven confirmation (not assumed from the paper)**: tracked, across the `Simulation/pressureDiff_0` run's first physical period (100 dumps), (a) the inlet valve's leaflet-tip minimum distance from the vessel axis (a direct "how open is the gap" measure, using the *actual* leaflet geometry — not the wall radius, which was tried first and is nearly the wrong signal since it only reflects the wall's own mild ±20% breathing, not the leaflet's much larger swing) and (b) fluid velocity at the valve gap's centerline. Found: the gap collapses to near-zero (closed) by τ≈0.53 with a sharp transient backflow spike at that instant, then — while the leaflet gap stays visually flat/shut through τ≈0.96 — the centerline velocity doesn't stay at zero: it drifts negative (backflow) through τ≈0.55-0.76, crosses zero, then **sustains and escalates a small positive (forward) leakage flow through τ≈0.77-0.96**, immediately before the valve mechanically reopens at τ≈0.97. That escalating-leakage window lines up with the reviewer's cited τ≈0.8-0.9.
- New `Template/DelayedClosureFigure.m`: renders 3 frames (dumps 82/86/90, τ≈0.81/0.85/0.89) zoomed to the inlet valve region (X∈[25,75]), velocity-**magnitude** contour (not density — matches what the reviewer specifically asked to see) + midplane velocity vectors + wall/valve surfaces shown for context (re-enabled just for this script's own local `showWall`/`showValve`-equivalent hardcoded `true`, without changing the general scripts' now-`false` defaults from #25). Uses **one shared color scale across all three frames** (computed as a `globalMax` pre-pass before rendering) — an initial version let each frame auto-scale its own colorbar, which would have made the escalating-jet story partly a color-scaling artifact rather than a genuine comparison.
- Saves both `.fig` (editable) and `.png` (300 DPI, manuscript-ready) per frame to `Visualization/DelayedClosure_frame0{82,86,90}.{fig,png}`.
- **Flagged to user, not yet acted on**: the reviewer's 28-Jun-2026 revision deadline has already passed (today is 2026-07-30) — worth confirming status with Dr. Alexeev before finalizing the revision. Also flagged as open questions: whether this 3D-oblique rendering style should be made to match the manuscript's existing Figure 3 style (a different, "lofted tube" rendering from the old `AnimePlotterMainFlow_FIXED.m`) more closely, and that the reviewer's alternative asks (vortex structures, leaflet WSS) haven't been addressed — only the velocity-vector visualization has.

---

## Current Solver Behavior (verified)
- Vessel diameter, valve extent, and peristaltic wavelength are all computed **once at init from the actual loaded geometry** — no more hardcoded/stale constants for these.
- Vessel wall and leaflet motion are both **position-based** (target radius computed directly each step from each node's own reference configuration), not velocity-integration — bounded by construction, can't drift.
- The LBM fluid solver is confirmed running and coupled to the moving structure (real, finite, developing density/flux fields).
- `LivePostProcessor.m` is self-diagnosing; if it still shows nothing, the Command Window will say why.

## Known Issues / Next Steps
1. Watch the fluid field over a longer run — max wall radius was still slowly trending upward (10 → ~16 by dump 576 of the pre-fluid-fix run) rather than settling into a clean bounded oscillation; worth confirming this is genuine physical response vs. a residual amplitude-calibration issue once a full post-fix run completes.
2. `pipeFromGeometry.m`'s ray-casting mask (`LBMl`) is still broken (option A from #7, not pursued since B unblocked the solver) — fine to leave unused/unwritten indefinitely, but if anything ever needs that MATLAB-side mask independently of the C++ solver, it'll need a real fix (likely needs proper end-cap triangles on the tube mesh, or a different inside/outside test entirely).
3. Current default run config: `Vessel_diam=20`, `Domain_offset=10` → domain `nY=nZ=40`, `nX=150`; 2 valves at X≈25/125; `DELTA_RO=0` (zero pressure-driven forcing — flow comes only from wall motion). Note this describes `Template/`'s baseline — the two active PACE runs (`Simulation/pressureDiff_0` at L_nd=5, `Simulation_shortVessel_L4/pressureDiff_0` at L_nd=4) use the paper-matched derivation from `buildConfig.m` instead (see session #16/#20), which computes different actual dimensions.
4. **Both PACE jobs need resubmitting** — `lbls_run` and `lbls_L4` both hit `TIMEOUT` at the `-t2-0` walltime cap (session #21) before finishing their ~385,000-step run. Either raise the walltime further and resubmit from scratch, or check whether the solver supports restarting from its last `results/N` checkpoint (`c_t_flContinue` in `par.m`/`clConstants.h` suggests a continuation mechanism exists but wasn't investigated this session).
5. OpenMP isn't scaling to the full 24 requested cores on PACE (~9 cores' worth of CPU time actually used per session #18) — worth profiling if faster turnaround matters, since PACE allocation hours are being spent on mostly-idle requested cores.
6. Reviewer-response figure (`DelayedClosureFigure.m`, session #27) only covers the velocity-vector ask — Referee #3's alternative suggestions (vortex structures in the wake, leaflet-surface wall shear stress) haven't been visualized, and the 28-Jun-2026 revision deadline has already passed as of this session.

## Key Files (cumulative across sessions)
```
Template/
├── source/
│   ├── clMethod.h                     — position-based wall/leaflet motion (solveVerlet1/2), waveLength usage
│   ├── clConstants.h / .cpp           — valveMinX/valveMaxX/waveLength members
│   ├── clVariables.h                  — clC::ini(): valve extent + waveLength computed from geometry
│   └── ConstDef.h                     — removed VALVE_MIN_X/MAX_X defines (now computed, not hardcoded)
├── geometry creator/
│   └── GeometryCreatorAmir.m          — Vessel_diam fix (was /2'd twice), stopped writing broken lbm.txt
├── generateGeometry.m                 — removed dead valve-extent computation
├── editConstDef.m                     — removed dead LineModifier calls for valve extent
├── buildConfig.m                      — paper non-dimensional parameters → derived raw LBM config (L_nd currently 4)
├── Makefile                           — PACE/RHEL9 fixes: -fast → -ipo -O3 -no-prec-div -fp-model fast=2 (drops implicit -static)
├── run.sbatch                         — PACE fixes: account gts-aalexeev3-afrench, intel/2021.9.0, $SLURM_CPUS_PER_TASK threads
├── computeVesselRadiusProfile.m       — NEW: per-X wall center/radius from actual wall-node positions (this frame)
├── computeVesselRadiusMask.m          — NEW: Y-midplane inside/outside fluid mask, built on the profile above
├── buildTrimmedStructureSurfaces.m    — NEW: wall/valve patch geometry from lsbl connectivity, valve radially trimmed
├── LivePostProcessor.m                — masked fluid, real surfaces (not scatter3), showWall/showValve toggles, midplane-only vectors
├── FlowAnimationExporter.m            — same visualization fixes as above, batch/parallel, saves PNG per frame
├── DelayedClosureFigure.m             — NEW: reviewer-response figure, inlet valve delayed-closure window, shared color scale
├── loadMultiBlock3D.m                 — correct scalar-field block-format loader
└── loadMultiBlock3Dv3.m               — correct vector-field block-format loader
```

---

## This Session (2026-08-12) — Flow Visualization Pipeline & Paper Figure

Goal: produce a publication figure of inlet-valve flow for *Pumping dynamics of single lymphangion*
(`OneDrive/Attachments/Pumping_dynamics_of_single_lymphangion_BD.docx`), in response to
reviewer/advisor requests. Along the way the post-processing was rebuilt and **two real bugs
were found and fixed** — both had already produced figures that were wrong.

### 1. Post-processors reduced to flow-only, then extended to all time steps
- `LivePostProcessor.m` and `FlowAnimationExporter.m` previously drew a density slice + wall +
  leaflet surfaces. Stripped to velocity vectors only, then streamlines.
- `LivePostProcessor.m` used to render **only the latest dump**; it now plays through every
  dump in ascending order and then keeps polling for new ones.
- `FlowAnimationExporter.m` gained `plotMode = 'vectors' | 'streamlines'`, with separate output
  folders per mode (`Visualization/` vs `VisualizationStreamlines/`) so one never overwrites the
  other — `makeAnimation.m` globs `*.png`, so mixing them in one folder would interleave frames.

### 2. `exportgraphics` produced variable-size frames → broke video encoding
- **Symptom**: all 182 PNGs written, then `VideoWriter` failed with "properties of the video is invalid".
- **Root cause**: `exportgraphics` crops to drawn content. With the axes hidden every frame came out
  a different size, and frame 1 (t=0, fluid genuinely at rest, no arrows) cropped to **475×26** —
  just the title. `makeAnimation` pinned all frames to frame 1's size.
- **Fix**: `print -dpng` instead, which honours the figure size regardless of content. `makeAnimation.m`
  additionally picks the **modal** frame size from an `imfinfo` header scan and resizes outliers,
  rather than trusting frame 1.

### 3. Arrow/streamline legibility — dynamic range is ~3 decades
- Measured: within one frame the median speed is ~5–15% of that frame's peak, and the peak itself
  swings ~5× over a cycle. A single linear length scale leaves typical arrows sub-pixel.
- **Vectors**: `drawSpeedQuiver.m` — gamma-compressed length (`arrowGamma = 0.35`) renormalised
  per frame (`arrowNorm = 'frame'`, 95th percentile), with the frame's reference speed captioned
  since length is then comparable *within* a frame but not between frames.
- **Streamlines**: `drawStreamlines.m` — evenly-spaced placement (Jobard & Lefer occupancy raster);
  raw seed grids gave 827 overlapping lines and a solid black mass. Density is set by `dSep` alone.
  - Integrated **forward and backward** from each seed, or a recirculation loop never closes.
  - `maxTurn = 270°` cap: a vortex core here is ~3 cells across, so an uncapped line laps it ~15
    times and paints the core solid. Accumulates **signed** turning (meanders cancel; only sustained
    rotation trips it).
  - Arrowheads are explicit filled triangles, not `quiver3` heads — `quiver3` sizes its head as a
    fraction of the shaft, which on a 2-cell direction marker is a couple of pixels.

### 4. BUG (fixed): `lambda_ci` was computing the *inverse* of a vortex criterion
- **Symptom**: none visible — the vortex panels looked plausible.
- **Root cause**: `gradient()` returns the **column-direction** derivative first. Arrays here are
  laid out (rows = X, cols = Z), so `[duxdx, duxdz] = gradient(ux)` assigned ∂ux/∂Z to `duxdx`.
  This transposes the velocity gradient tensor; the discriminant becomes `ω² − 4a²` instead of
  `4a² + s² − ω²`, i.e. it fires where **strain** beats rotation and stays silent on real vortices.
- **Verified** against analytic fields — solid-body rotation returned 0.000 (missed), pure strain
  returned 1.000 (false positive).
- **Fix**: `swirlingStrength.m`, one shared implementation with a built-in self-test:
  `matlab -batch "swirlingStrength('selftest')"` → rotation 1, shear 0, strain 0.
- **Impact**: the buggy/corrected per-frame maxima correlate at 0.97 (strain and rotation are
  co-located at the valve-gap shear layer), so frame *rankings* survived — but any earlier claim
  resting on λ_ci is void. The recirculation conclusion (below) does **not** depend on it.

### 5. BUG (fixed): τ → dump mapping was wrong; figures were generated at the wrong frames
- **Root cause**: `DelayedClosureFigure.m` carries the comment
  `targetDumps = [82, 86, 90]; % tau ~= 0.81, 0.85, 0.89`, implying `dump = 100τ + 1`.
  **That comment is wrong.** Amir supplied the correct anchors: τ = 0.8 → dump 131, τ = 0.9 → 141.
- **Correct mapping**: `dump = round(100·τ) + 51`, i.e. `τ = (dump − 51)/100`. One cycle = dumps 51–151.
- **`DelayedClosureFigure.m` is still targeting the wrong frames** — its phase labels are off by 50
  dumps. Not fixed this session; fix before reusing it.

### 6. Recirculation analysis — two persistent zones, one per valve
Method deliberately avoided relying on any single metric:
- **Counterflow test** (gradient-free, unaffected by the λ_ci bug): do forward and reverse axial
  flow coexist at the *same* axial station? A global oscillation reverses the whole cross-section
  at once; only a recirculation cell has both signs side by side. → **99% of frames**, sharply
  localised at X 40–69 (inlet valve) and X 140–169 (outlet valve), essentially nothing in the open
  tube between. Peak: 65% of frames at X = 44.
- **Streamlines** independently show closed spiral cores at the same two locations.
- Reverse-flow speed reaches 1.24× the mean forward speed on average, peaking at 3.8×.
- Naïve "27% of vessel area is reversed" is misleading — the net flux flips sign through the run,
  so much of that is global oscillation, not recirculation.

### 7. Symmetry — mirrored, not just symmetrically seeded
- The vessel is axisymmetric; asymmetric-looking streamlines read as a solver artifact.
- Symmetric **seeding** is not sufficient: the even-spacing rule lets whichever line is processed
  first claim the shared corridor, so mirror-image seeds still truncate differently.
- **Fix**: `drawStreamlines` `mirrorZ` — integrate the upper half only, reflect it. Exact by construction.
- Axis taken from the wall nodes, not the domain centre: **Z = 19.500, spread across phases 4e-05**,
  which is the check that makes mirroring legitimate rather than cosmetic.
- `mirrorField.m` handles scalar fields, with an explicit **parity** argument. Vorticity is **odd**
  under the reflection (ux even, uz odd, ∂/∂z̃ = −∂/∂z) — copying it across unchanged would draw a
  co-rotating pair, which a symmetric flow cannot produce.

### 8. Vorticity normalised by the contraction period (Alexeev's request)
- `ω* = ω·T`, dimensionless — reads as radians of rotation per contraction cycle.
- **T = 192,300 lattice steps** = 100 dumps × 1923 steps/dump.
- **Trap**: `par.m` has `c_t_oscPeriod1 = 384684` steps = **200.04 dumps, twice the cycle**.
  Normalising by that would put every value out by a factor of 2. Confirmed 100 dumps three ways:
  Amir's τ anchors, `oscPeriod1/2 / Tdumpstep` = 100.02, and autocorrelation of the mean wall
  radius (peak at lag 97; FFT ~91, coarse because the run covers only ~1.8 cycles).
- `cyclePeriodSteps.m` derives T from the run's own `par.m` and warns on a >5% mismatch.

### 9. Paper figure
`PaperFigure_ValveFlow.m` → `Snapshots/Fig_ValveFlow.png` (7.0 × 1.0 in, 600 dpi) and `.pdf`.
Style taken from the manuscript, not from MATLAB defaults:
- Serif (Times New Roman) — manuscript Figures 3–6 are serif.
- `(a)`–`(e)` panel labels in parentheses; thin panel frames as in Figure 2's snapshot grid;
  colorbar label above the bar as in Figure 2.
- Sized in **inches** to a 7.0 in text width so it is placed at 100% and type lands at its stated
  point size. Rescaling in Word is what makes fonts inconsistent between figures.
- Content: 1×5 row, τ = 0.75…0.95, streamlines above the vessel axis / normalised vorticity below.
  Legitimate only because of the symmetry check in item 7.
- **The manuscript states ideal valves close at τ = 0.75**, so this sequence starts exactly at ideal
  closure and documents the delayed-closure window — ties directly to the Figure 3 discussion.
- τ = 0.75 is genuinely quiescent (max speed 4.6e-07 vs 1.2e-04 at τ = 0.85, **~265× slower**), which
  is why its vorticity half is blank. Note: streamlines normalise direction and discard magnitude,
  so panel (a) looks as substantial as the others — worth a caption note.

#### 9a. Text rendering — LaTeX, not MATLAB's `tex` (revised 2026-08-12, later)
Two label changes requested after first review of the PDF. **Plot content was not touched** —
only text objects.
- **Panel labels moved to top CENTRE** of each panel (was top-left): `text(ax, 0.5, 1.0, ...)` in
  normalised axis units with `'HorizontalAlignment','center'`, so placement is identical on every
  panel regardless of data ranges.
- **τ now rendered with `'Interpreter','latex'`**, string `$\mathrm{(a)}\;\tau = 0.75$`.
  MATLAB's default `'tex'` interpreter draws Greek from its own symbol font, which is a visibly
  different glyph from the Computer Modern τ that `Main_Rev7.tex` produces — the mismatch is what
  makes a figure look imported. `\mathrm` keeps the panel letter upright while τ stays italic,
  matching the manuscript's own Figure 3a.
- The `ωT` colorbar label was switched to LaTeX **as well**, for consistency: leaving it on `tex`
  would have put two different Greek fonts in one figure (a MATLAB ω beside a CM τ), which is more
  conspicuous than either alone.
- **Gotcha**: under `'latex'` MATLAB **ignores `FontName`**. The `fontName = 'Times New Roman'`
  config now only affects the colorbar's numeric tick labels (5.5 / 0 / −5.5), which were left as
  Times deliberately — changing them would alter the figure beyond text. If full consistency is
  wanted later, it is one line: `cb.TickLabelInterpreter = 'latex'`.

#### 9b. Where the figure lives / how to regenerate
```
Simulation/pressureDiff_0/Snapshots/Fig_ValveFlow.pdf   <- the one for Overleaf (vector text)
Simulation/pressureDiff_0/Snapshots/Fig_ValveFlow.png   <- 600 dpi raster, slides/preview
```
Regenerate with, from `Simulation/pressureDiff_0`:
```
matlab -batch "PaperFigure_ValveFlow"
```
The Overleaf `\includegraphics` currently points at `Final Physics of Fluids/Fig_ValveFlow.pdf`,
a placeholder path never checked against the project — see manuscript open item 1.

### Files added this session (all in `Template/`, mirrored into `Simulation/pressureDiff_0/`)
| File | Purpose |
|---|---|
| `drawStreamlines.m` | evenly-spaced streamlines, bidirectional, turn-capped, triangle arrowheads, mirroring |
| `drawSpeedQuiver.m` | velocity arrows with gamma-compressed / per-frame-normalised length |
| `overlayStructure.m` | wall + leaflet midplane cross-section (wall = envelope per X; leaflets ordered along their own principal axis, one `lso` group at a time) |
| `swirlingStrength.m` | signed λ_ci, **with self-test** |
| `inPlaneVorticity.m` | ω = ∂uz/∂x − ∂ux/∂z, **with self-test** |
| `mirrorField.m` | reflect a scalar field about the axis, explicit even/odd parity |
| `cyclePeriodSteps.m` | contraction period in lattice steps, from `par.m` |
| `estimateSpeedLimits.m` | global speed range for colour scales |
| `makeAnimation.m` | PNG folder → MP4 + GIF, headless (`print`-based, no `getframe`) |
| `InletValvePanels.m` | per-timestep inlet-valve panels (streamlines + λ_ci), all dumps |
| `StreamlineSnapshots.m` | curated τ snapshots + vertical strips |
| `PanelFigures.m` | 1×5 row figures: streamlines / vorticity / split |
| `PaperFigure_ValveFlow.m` | **the manuscript figure** |

### Outputs in `Simulation/pressureDiff_0/`
```
Snapshots/Fig_ValveFlow.png|.pdf      <- the paper figure
Snapshots/Row_{Streamlines,Vorticity,Split}.png
Snapshots/Streamlines_tau*.png, Vorticity_tau*.png, Vortex_tau*.png
Snapshots/{Streamlines,Vorticity}_AllPhases.png
VisualizationInletStream/ , VisualizationInletVortex/   (182 panels each)
VisualizationStreamlines/ , Visualization/              (182 frames each, full vessel)
{Inlet,}{Stream,Vortex,Flow}Animation.mp4|.gif
RecirculationMap.png    (NOTE: its swirl panel used the pre-fix λ_ci — regenerate if used)
```

### Open items / handoff
1. ~~**Move `Fig_ValveFlow` into the manuscript**~~ — **DONE** in the Overleaf half of this session
   (see the next section). Placed after Figure 4 as `\label{fig:DelayedClosureFlow}`, not as
   "Figure 7" as originally guessed here. Draft caption written at the time:
   > **Figure 7:** Flow through the inlet valve during the delayed-closure window for Γ = 0.14 and
   > ∆P = 0. Panels show the mid-plane at (a) τ = 0.75, (b) τ = 0.80, (c) τ = 0.85, (d) τ = 0.90,
   > and (e) τ = 0.95, where ideal valves close instantaneously at τ = 0.75. In each panel the upper
   > half shows streamlines and the lower half shows the vorticity ω normalised by the contraction
   > period T; the flow is symmetric about the vessel axis (dotted line). Grey lines mark the vessel
   > wall and black lines the valve leaflets.

   ~~Verify Γ and ∆P~~ — **DONE**, see "Parameter verification" below: Γ = 0.14, ∆P = 0 and L̄ = 5
   all confirmed from the run's own data. **K = 0.065 is still unverified.**
2. ~~**Rebuttal letter** — needs the reviewer comments~~ — the comments **do** exist:
   `SinglePaperFigures/Reviewer Comment.pdf` (Referee 3, MS #POF26-AR-02859R). The letter itself is
   still **not drafted**. This note originally said the comments "were never in this repo" — that
   was wrong; they were saved during the Overleaf half of the session.
3. `DelayedClosureFigure.m` still uses the wrong τ mapping (item 5). **Its phase labels are off by
   50 dumps.** If any existing manuscript figure came from it, that is a substantive error, not a
   cosmetic one — check before resubmission.
4. `RecirculationMap.png` middle panel predates the λ_ci fix.
5. Vertical `Vorticity_*.png` strips from `StreamlineSnapshots.m` are **not** period-normalised;
   only `PanelFigures.m` and `PaperFigure_ValveFlow.m` are.
6. Viscosity / Γ convention — see the ⚠ OPEN QUESTION block below. Affects every Γ in the paper.

---

## This Session (2026-08-12), continued — Manuscript Revision in Overleaf

Separate from the MATLAB/PACE work above: this half of the session worked directly in the
paper's Overleaf project (`Main_Rev7.tex`, "Pumping dynamics of a single lymphangion", MS
#POF26-AR-02859R) to close out the two outstanding Referee 3 comments from the round-2 decision
letter (`Reviewer Comment.pdf`, saved to `SinglePaperFigures/` the same day). All other comments
(1, 2, 4, 6) were already addressed in a prior round.

### Outstanding comments identified from the round-2 letter
- **Comment 3** (Referee 3): the four-snapshot Figure 3 (τ=0, 0.25, 0.5, 0.75) does not cover the
  delayed-closure/backflow window the paper itself identifies at τ≈0.8–0.9 for the inlet valve.
  Reviewer asked for a dedicated figure or panel showing streamlines, wake vortex structures, or
  leaflet WSS specifically during that phase.
- **Comment 5**: the Discussion's single-lymphangion-limitations paragraph acknowledged that
  downstream compliance and upstream pressure pulses "may influence" delayed closure without
  stating a direction — reviewer asked for at least a qualitative amplify/attenuate call.

### False start: six-panel Figure 3, then reverted
- First attempt extended Figure 3 itself to six panels (τ=0, 0.25, 0.5, 0.75, 0.8, 0.9) using a
  PowerPoint-authored `SnapshotFigureAfterReview.png`/`.pdf`, with a matching sentence added to the
  in-text walkthrough.
- **Two PowerPoint export bugs hit along the way** (both real, both fixed by the user in
  PowerPoint, not by me): (1) "Save as Picture" flattens the whole slide to a raster PNG, losing
  the vector text/labels that every other figure in the paper has — fixed by exporting via
  `File → Export → Create PDF/XPS` instead. (2) The exported PDF was clipped to the slide's
  canvas bounds, silently dropping the bottom row (τ=0.9) even though it was visible while
  editing — fixed by resizing the slide/group so the full 6-row content fit inside the slide
  frame before re-exporting.
- **Reverted per the user's supervisor's direction**: rather than growing Figure 3, the delayed-
  closure visualization was moved to its own dedicated figure (matching the reviewer's explicit
  "add a dedicated figure, **or** an additional panel" wording — the supervisor chose the former).
  Figure 3 was restored to its original four-panel scope (`SnapshotFigure.pdf`, τ=0/0.25/0.5/0.75,
  streamlines kept as a general improvement), and the six-panel-specific sentence was removed from
  the in-text walkthrough.

### Comment 5 — resolved with a qualitative direction
Added to the Discussion, replacing the direction-less "may influence" sentence:
> upstream pressure pulses arriving during valve closure reinforce the transient driving pressure
> and amplify the delayed closure/backflow; downstream compliance damps outlet-side pressure
> transients and likely attenuates it.

### Comment 3 — resolved using the MATLAB-side `Fig_ValveFlow` figure (see item 27/29 above)
The supervisor approved using the actual data-driven figure produced earlier this session by the
MATLAB/PACE pipeline (`PaperFigure_ValveFlow.m` → `Fig_ValveFlow.png`/`.pdf`, five panels at
τ=0.75/0.80/0.85/0.90/0.95, streamlines above the axis / period-normalized vorticity below,
mirrored for symmetry) rather than a new hand-built PowerPoint figure. Placement decided with the
user: **after Figure 4** (`fig:elastic-delayed`, which establishes *when* the delay occurs) and
**before** the flow-rate figure (`fig:ideal_elastic_flowrate`, which quantifies the backflow it
causes) — so the new figure sits in between as the visual explanation of the mechanism.
- New `\begin{figure*}` block added, `\label{fig:DelayedClosureFlow}`, referencing
  `Final Physics of Fluids/Fig_ValveFlow.pdf` — **filename not yet confirmed against what's
  actually in the Overleaf project**, needs checking.
- Caption and a new in-text paragraph written to match the manuscript's existing prose style/depth
  (no other section's writing was altered) — describes quiescent flow at τ=0.75, counter-rotating
  recirculation cells forming as the gap narrows through τ=0.80–0.90 (the mechanism behind the
  Figure 5a backflow), and weakening recirculation by τ=0.95.
- **Γ=0.14, K=0.065, ΔP=0 used in the caption are the paper's usual baseline values but were not
  independently re-verified against this specific run's actual config** — the MATLAB-side notes
  (item 29 above) flag that these were inferred from the `pressureDiff_0` folder name, not read
  from a parameter file. Confirm before submitting.
- All new manuscript text wrapped in `\hl{}` (green highlight via `soul`), following the paper's
  existing convention for marking round-2 changes in the "Marked-Up Manuscript" upload, with
  `% Reviewer 3 Comment N` / `%Start` / `%End` comment markers matching the pattern already used
  for the Comment 2/5 edit.

### Parameter verification for `Fig_ValveFlow` (done 2026-08-12, closes MATLAB-side item 1)
Measured from this run's own data, not inferred from the folder name:
- **ΔP = 0** — confirmed: the sweep is `linspace(dP,dP,1)`, a single run, folder `pressureDiff_0`.
- **L̄ = 5** — confirmed: valve roots at X = 38 and X = 138 → L = 100, 2R = 20.
  (Note `buildConfig.m` currently has `L_nd = 4`; it was edited **after** this run was generated, so
  it does not describe this run. Cite the measurement, not `buildConfig`.)
- **Γ = 0.14** — confirmed *as configured*: run period T = 192,300 steps vs the 192,342 that
  `buildConfig.m` derives for Γ = 0.14 assuming ν = 1/6. **0.02% agreement** — not a coincidence.
- **K = 0.065** — NOT verified. Would need the leaflet bending stiffness compared against P0.

### ⚠ OPEN QUESTION — viscosity convention, affects every Γ in the paper
Found while verifying Γ above. Two mutually inconsistent viscosities in the same run:
- `ConstDef.h`: `MU_NUMBER = 1/6`, and `buildConfig.m` derives the period from **ν = 1/6** → Γ = 0.1400
- The solver actually ran with `lblamda = -1.4724` → τ = 0.679164 → **ν = (τ−½)/3 = 0.0597** → Γ = 0.2339

That is a factor of **2.79 in ν, 1.67 in Γ**. The period being tuned so precisely to
Γ = 0.14-with-ν = 1/6 says ν = 1/6 is the design intent, and `MU_NUMBER` only feeds unrelated
legacy formulas (`MAGNETIC_FORCE`, `KSI_NUMBER`), which suggests this solver's τ→ν relation is
simply not the textbook BGK one. **But that has not been established.**

This is **not specific to the new figure** — every Γ in the manuscript comes from the same
derivation. If the standard relation does hold here, all Γ labels are off by 1.67×.
Settle it with a Poiseuille calibration run, or by asking whoever knows the solver's convention,
**before resubmission**. Keeping Γ = 0.14 in the new caption is right either way — it is
consistent with the rest of the paper.

### Still open (manuscript side)
1. Confirm the `Fig_ValveFlow.pdf` filename/path matches what's actually uploaded to the Overleaf
   project (a placeholder path was used based on the MATLAB output convention).
   **Source of truth:** `Simulation/pressureDiff_0/Snapshots/Fig_ValveFlow.pdf` — re-upload it,
   the file was regenerated after the label changes in item 9a above.
2. ~~Verify Γ and ΔP~~ — **DONE**, see "Parameter verification" above. Γ = 0.14, ΔP = 0, L̄ = 5 all
   confirmed from the run's own data. **K = 0.065 remains unverified** — either verify it or drop
   it from the caption.
3. Point-by-point response letter to Referee 3 — **not yet drafted**, and unblocked: the reviewer
   text is in `SinglePaperFigures/Reviewer Comment.pdf` / `Referee 3 Review Attachment 1`.
4. Alt-text Word doc (`Pumping dynamics of a single lymphangion alt text PoF.txt`) needs an entry
   added for the new figure — AIP Publishing requires alt text for all figures on resubmission.
   Draft, from the figure's actual content: *"Five panels showing flow through the lymphatic inlet
   valve at five successive phases of the contraction cycle. In each panel streamlines are drawn in
   the upper half of the vessel and colour-mapped vorticity in the lower half; the valve leaflets
   converge toward the vessel axis, and a counter-rotating vortex pair develops at the leaflet tips
   as the gap narrows."*
5. Clean "Article File" (accept all `\hl` highlights) and "Marked-Up Manuscript" (keep them) both
   need to be exported from Overleaf once the figure filename/params above are confirmed.
6. **Sequencing note**: items 1 and 5 depend on the figure being final. It now is (item 9a), unless
   the colorbar tick numbers are also switched to LaTeX — the one remaining cosmetic choice.

---

---

## This Session (2026-08-13) — Valve Under-Opening Investigation & Faster Run

### Why this is happening
The supervisor concluded that **the valve opening in the current simulations does not match
the case presented in the paper**. The original paper's simulation code and results are
**LOST** — they cannot be inspected or re-run, so the only way to reconcile is to reproduce
the paper's behaviour from the current code.

Amir's hypothesis: the reduced opening is caused by **slower flow not pushing the leaflets
open hard enough**, so a faster run was launched to test it.

A **new Claude Code session will do the two-run comparison.** Everything it needs is below.

### The two runs

| | OLD (`SimulationPaper/pressureDiff_0`) | NEW (`Simulation/pressureDiff_0`) |
|---|---|---|
| Gamma requested | 0.14 | **0.28** |
| Gamma actual | 0.234 | 0.468 |
| cycle T | 192,342 steps | **48,085 steps** (4x faster) |
| `oscPeriod1` | 384,684 | 96,170 |
| `Tdumpstep` | 1923 | 480 |
| dumps per cycle | 100 | 100 (same) |
| total dumps | 182 (of 200) | 200 expected |
| Tstop | 384,685 (2 cycles) | 96,171 (2 cycles) |
| runtime | — | ~15 h |
| geometry | L=5, R=10, 175x40x40 | identical |
| LB_TAU_NUMBER | 0.6791635268 | unchanged |

**SUPERSEDED for the NEW run** - the table above describes the Gamma=0.28 attempt that was
cancelled and deleted. The run actually going now is the corrected one described under
"RESOLVED" below (Gamma = 0.2 true, K effective 0.0325, ~35 h).

**Amir has explicitly deprioritised the absolute Gamma value** — the goal of the new run is
simply "faster than the previous one". Do not spend the new session on Gamma labelling.
(The 1.67x Gamma error is still documented above under the viscosity OPEN QUESTION.)

### Sanity check already done — dumps 1 and 2, at MATCHED cycle phase
Both runs put dump N at tau ~ N/100, so dump 2 is the same phase in each.

| | NEW | OLD | ratio |
|---|---|---|---|
| wall radius (dump 1 / 2) | 10.000 / 10.067 | 10.000 / 10.067 | **identical** |
| valve gap (dump 1 / 2) | 1.219 / 1.194 | 1.219 / 1.194 | **identical** |
| mean speed (dump 2) | 5.109e-05 | 1.464e-05 | **3.5x** |
| max abs(rho-1) (dump 2) | 4.750e-04 | 3.014e-05 | **15.8x** |

- Identical wall radius and valve gap at matched dump confirms the geometry rebuilt correctly
  and the two runs are phase-aligned the same way. It also implies the **tau = (dump-51)/100
  mapping carries over to the new run** (worth re-confirming once ~100 dumps exist).
- Mean speed 3.5x against a predicted 4x (U ~ Gamma_req^2). As expected.
- Pressure 15.8x ~ 16 = (0.28/0.14)^4, the signature of an **inertia-dominated** pressure
  field (unsteady term ~ rho*(U/T)*L ~ Gamma^4, vs viscous ~ Gamma^2). Correct behaviour for
  a higher-Womersley case.

### MEASURED: the valve really is under-opening (supervisor confirmed, quantitatively)
Minimum 3D leaflet-to-leaflet separation over all 182 dumps of the OLD run, with
DeltaV = 2*pi*R*A*L = 18,850 so DeltaV^(1/3) = 26.61:

```
  MAX opening  = 2.194 cells  at dump 129 (tau = 0.78)
  MIN opening  = 0.040 cells  at dump 158 (tau = 1.07)   (i.e. fully closed)
  delta_max / DeltaV^(1/3) = 0.0825
  delta_max / 2R           = 0.110
```

**The paper's Figure 3a peaks at delta/DeltaV^(1/3) ~ 0.29 (Gamma=0.1) to ~0.41 (Gamma=0.2).**
So this run opens **3.5x to 5x less than the paper**. In physical terms the paper's valve
opens to ~55% of the vessel diameter; this one reaches 11%. The supervisor's read is correct
and now has a number attached.

### STRONGEST CANDIDATE CAUSE — the viscosity bug makes the leaflets 2.79x too STIFF
This is probably more important than the flow-speed hypothesis, and it is a direct
consequence of the same `mu = 1/6` error documented earlier:

- `buildConfig.m` sets the leaflet bending rigidity from `Db = K * 2*pi*R^3 * P0`,
  with `P0 = 8*mu*L*Q0/(pi*R^4)`.
- `Q0 = 2*DeltaV/T` depends only on geometry and the period, both of which are written to the
  config — so **P0 is proportional to mu alone**.
- With mu assumed 1/6 but actually 0.0597, the real P0 is 2.79x SMALLER than assumed, so the
  achieved dimensionless stiffness is **K_actual = K_requested x 2.79**.

| K requested | K actually simulated |
|---|---|
| 0.065 (paper's flexible case) | **0.181** |
| 0.325 (paper's stiff case) | 0.907 |

So every run intended as the flexible K=0.065 case has been running at K~0.18 — **2.8x
stiffer than intended**, and a third of the way to the paper's stiff case. Stiffer leaflets
open less. This alone is a strong candidate for the under-opening.

### CAUTION on the faster-flow hypothesis
By the paper's own trend, doubling Gamma (0.1 -> 0.2) only raises delta/DeltaV^(1/3) from
0.29 to 0.41 — a factor of **1.4**. Closing a **3.5x** gap by flow speed alone would need far
more than the 2x Gamma increase in the new run. The faster run is still worth having as a
data point, but on the paper's own scaling it should NOT be expected to reproduce the paper's
opening. **Fixing the stiffness error is the more likely route**, and it costs nothing to
test: it is the same one-line `buildConfig` fix.

### FLAG — valve timing does not match the paper's phase convention
Using tau = (dump-51)/100, the measured OLD-run valve is:
- **fully closed** (gap 0.040) for tau ~ 0.06 to 0.46
- **open**, peaking at 2.194, for tau ~ 0.54 to 1.05, with the **maximum at tau = 0.78**

But the paper states ideal valves open at tau = 0.25 and close at **tau = 0.75**, and the
manuscript's delayed-closure window (tau ~0.77-0.96) is supposed to be *after* closure.
**At tau = 0.80 this run's valve is at its widest opening, not closed.**

The measured open window [0.54, 1.05] maps onto the ideal [0.25, 0.75] under a shift of
about **+0.29**, i.e. `dump = 100*tau + 80` rather than `+51`. Two possibilities, unresolved:
1. the tau offset is wrong by ~29 dumps, or
2. the valve timing genuinely differs from the paper.

**This matters for `Fig_ValveFlow`** (the new manuscript figure): it is captioned as the
post-closure delayed-closure window, but under the measured behaviour those five phases
straddle peak opening. Resolve before the figure is submitted.

### Tools available for the comparison session
All in `Template/` and mirrored into each run folder:
- `InletSplitPanels.m` — per-dump streamline/vorticity panel, all dumps + animation
- `PaperFigure_ValveFlow.m` — the 5-phase paper figure
- `ValveSurface3D.m` — valve as a 3D surface, cropped to the valve, prints min leaflet
  separation per phase (**this is the valve-opening measurement**)
- `PanelFigures.m`, `StreamlineSnapshots.m`, `InletValvePanels.m`
- helpers: `drawStreamlines`, `drawSpeedQuiver`, `overlayStructure`, `swirlingStrength`,
  `inPlaneVorticity`, `mirrorField`, `cyclePeriodSteps`, `makeAnimation`, `midplaneSlice`

Scratch analysis scripts used to produce the numbers above are in this session's scratchpad
(`valveOpening.m`, `compareRuns.m`, `verifyParams.m`, `vortexAudit.m`) — re-create if needed,
they are short.

### RESOLVED — the viscosity bug is fixed, and a corrected run is now going

The Gamma=0.28 run described above was **cancelled and its folder deleted** before producing
results. Amir's position: he does not care about the absolute Gamma value, he wants a run that
reproduces the paper's valve opening. The following was then applied.

#### 1. `source/ConstDef.h` — `LB_TAU_NUMBER` 0.6791635268 -> **1.0**
tau = 1.0 gives nu = cs^2 (tau - 1/2) = **1/6**, which is what `buildConfig.m` had always
assumed. This single change removes BOTH errors at once:
- `Gamma_actual = 1.67 * Gamma_requested`  -> now 1:1
- `K_actual = 2.79 * K_requested`          -> now 1:1 (the leaflets stop being 2.8x too stiff)

tau = 1.0 is also the numerically best-behaved BGK relaxation time. **Requires a rebuild.**

#### 2. `buildConfig.m` — viscosity is READ, never assumed
New local function `readLbTau()` parses `LB_TAU_NUMBER` out of `source/ConstDef.h`
(comment-aware, so a historical value in a `//` line cannot win) and derives
nu = (tau - 1/2)/3. The hardcoded `mu = 1/6` is gone.

It now also **asserts** what the configuration will actually produce, and errors at build time
rather than after a multi-hour run:
```
Gamma requested 0.2000 -> T = 94248 steps -> Gamma achieved 0.2000
K     requested 0.0650 -> D_b = 0.69333   -> K     achieved 0.0650
Expected peak valve opening delta/DeltaV^(1/3) ~ 0.41 (paper Fig 3a trend)
```

#### 3. `buildConfig.m` — new `paper.valveSoften` knob (ABSOLUTE leaflet softening)
Divides the raw bending rigidity EIv directly, leaving period / Reynolds / geometry untouched.
There is no "absolute softness" that is not also a change in K — deflection is load/stiffness,
so K IS the governing ratio — but this knob changes only the leaflet, which is what was wanted.
The effective K is printed so a softened run is never mistaken for a paper K value.

| valveSoften | EIv | effective K | predicted delta |
|---|---|---|---|
| 1 | 12.33 | 0.0650 | 0.21 |
| **2 (set)** | **6.16** | **0.0325** | **0.42** |
| 3 | 4.11 | 0.0217 | 0.63 |

**Why 2:** two independent framings agree. (a) EIv 6.16 is essentially the OLD under-opening
run's absolute stiffness of 6.04 — same leaflets, but the corrected viscosity now loads them
2.79x harder. (b) delta ~ 1/K anchored on the measured 0.0825 at K=0.181 predicts 0.42,
against the paper's 0.41.

**Caveat:** delta ~ 1/K is an unvalidated linear model, and effective K = 0.0325 is HALF the
paper's stated 0.065. So this run targets the paper's *opening*, not its *stated stiffness*.
Record that distinction if it feeds the manuscript.

### The corrected run — `Simulation/pressureDiff_0`, launched 2026-08-13 ~18:50

Verified after `run_sweep.m` + rebuild:

| | value | how verified |
|---|---|---|
| **nu** | **0.166667** | `lbtau.txt` = 1, `lblamda.txt` = -1, **written by the solver at startup** |
| **Gamma** | **0.2000** | recomputed from runtime `par.m` + solver's own tau |
| K nominal / effective | 0.065 / **0.0325** | `valveSoften = 2` |
| EIv | 6.163 | `lsak` min = 0.066716, exactly half the unsoftened 0.133432 |
| geometry | 175x40x40, R=10.000, valve roots X=38/138, L/2R=**5.000** | identical to the paper run |
| CYCLE_PERIOD | 188,496 -> cycle T = 94,248 | |
| dumps | 200 total, **100.1 per cycle** | `Tdumpstep` = 942 |
| cycles | 2.00 | `Tstop` = 188,496 |
| exe | Release OpenMP 18:50, newer than ConstDef.h 18:49 | |

**The definitive bug check is `lbtau.txt`**, because the solver writes it at startup from the
value actually compiled in — not from what the header says:
```
cat Simulation/pressureDiff_0/lbtau.txt
   1         -> fixed
   0.679164  -> the old bug
```

**Runtime: ~35 h** (`R34:52` steady), not the ~29 h originally estimated. The run is ~5,200
steps/h vs the cancelled run's ~6,400. Most likely the softened leaflets costing more work per
step in the FSI coupling. Expect completion early on 2026-08-15.

**Health at 0.11%:** M 280000.0 -> 279999.9 (negligible drift), R ~ 0.000 (stable density),
N3372/**0** detached nodes, S (max leaflet strain) 0.006 -> 0.007 climbing gently through
ramp-up. `avgl` (inlet flux) rising as expected — ~0.6x the cancelled Gamma=0.468 run at
matched cycle phase, against ~0.51 predicted from Q ~ 1/T. Nothing anomalous.

**LEFTOVER RISK:** `x64/Debug/LBLS.exe` (dated 07-13) is still present and was compiled with
the OLD `LB_TAU_NUMBER`. Launching it — stale shortcut, VS defaulting to Debug — silently
reinstates the bug, and nothing in the console output would reveal it (the M/R/S status line
does not print nu). Delete or rebuild it.

### Plan for the comparison session
1. Wait for the run (~35 h, 200 dumps). A spot check at ~25% (~9 h in, past the first cycle)
   is worthwhile — if delta is already far below expectation, that saves a day.
2. Measure peak valve opening with `ValveSurface3D.m` (prints minimum leaflet separation per
   phase). Divide by DeltaV^(1/3) = 26.61 and compare against the baseline **0.0825**:
   - **~0.42** -> stiffness diagnosis confirmed, paper opening reproduced
   - **~0.21** -> delta ~ 1/K roughly holds but something else also suppresses the opening
   - **~0.08** -> softening is not the mechanism; look elsewhere
3. Re-derive the tau offset for this run from the wall-radius trace. Dumps per cycle are 100
   in both runs, so only the offset is in question — but see the timing FLAG above, which is
   still unresolved and affects `Fig_ValveFlow`'s caption.
4. Compare like-for-like against `SimulationPaper/pressureDiff_0` at matched cycle phase, the
   way `compareRuns.m` did for dumps 1-2.
## This Session (2026-09-23) — Root cause of the valve under-opening: the fluid never saw the wall move

Context: the French Collaboration paper (ex vivo / in vitro / in silico) needs its simulation curves
regenerated, because its current ones (`French Collaboration/LymphaticsData*.xlsx`, Fig1-3) came from
the lost code. The plan: run one sample case, compare with the lost-code data, and if it matches, run long sims on PACE.

### Also found after the 08-13 notes (undocumented work from 08-13/14)
- `Template_v2/`: stiffness magnitudes now live in ConstDef.h (`LS_*STIFFNESS_1/2`), the mesh stores
  topology weights, and `VerifyV2Equivalence.m` checks that v1 and v2 give the same physics. It has a direct knob, `paper.leafletBendStiff`.
- `Simulation_v2/pressureDiff_0`: `results1/` holds 56 dumps at k_b = 0.8 and `results/` holds 8 dumps at k_b = **0.08**.
  Its compiled ConstDef.h has 0.08, but its buildConfig.m says 0.8, so the two are out of sync. It stopped on 08-14.

### Finding 1 — the opening metric did not match the paper's definition
`ValveSurface3D.m`/`CompareValveOpening.m` measure the MINIMUM leaflet separation. The paper and the old
data use the MAXIMUM tip-to-tip distance (the widest point of the orifice). Leaflet free edge = boundary nodes of the
leaflet's `lsbl` triangles (every triangle is listed twice, so a boundary edge appears 2x, not 1x) without the 256/512 flag,
which gives 13 nodes per leaflet. With the paper metric, SimulationPaper peaks at 4.28 cells = **0.16** (not 0.0825),
against the old data's 11.5 cells = 0.40.

### Finding 2 — opening does not depend on stiffness
`Simulation_v2/results1` (k_b 0.8, 5x softer) reaches a peak inlet tip-to-tip of 4.32, the same as baseline 4.28.
Fluid force on leaflet nodes (`lsfb`) is ~1e-6 per node, about 300x below a load of a P0 pressure drop. The leaflets were
being moved almost entirely by their anchors, not by the flow.

### Finding 3 (ROOT CAUSE) — prescribed-node velocity was zero where the fluid reads it
- `solve()` accumulates `vls.V.avr` (the wall velocity used by the LB bounce-back) **after `solveVerlet2`**.
- Position-based motion (fix #6/#2) set `V = (target-current)/dt`. `solveVerlet1` puts the node on its
  target, so in `solveVerlet2` that expression is ~0 (the saved wall `lsv` is 1e-16).
- Result: the fluid saw a wall that jumped position each step but never moved. The fluid next to the wall moved at ~2e-6,
  against a wall speed of ~1e-4. Q/Q0 peaked at ~0.05, where the old data has ~1.5 (30x). There was no pressure behind the valves, and
  none of the valve-opening tuning (tau, K, valveSoften, v2 stiffness) could have worked.
- **Fix (Template_v2/source/clMethod.h, all 4 blocks; backup `clMethod.h.bak_before_wallvel`)**: prescribed nodes now
  carry the analytic rate `dr/dt = OSC_AMPLITUDE*sin(2*pi*x/lambda)*(2*pi/RAMP_PERIOD)*cos(2*pi*t/RAMP_PERIOD)` along
  e_r. In Verlet1 it is added to the (target-current)/dt tracking term, and in Verlet2 it is used alone.

### Finding 4 — geometry convention of the old data
The old sheets normalise by DeltaV^(1/3) = 28.67, so DeltaV = 23562 and **L = 125** (their "L = 5" means L/l = 5, l = 25).
Revived runs used L = 100. The sample case therefore uses `paper.L_nd = 6.25`.

### Still open, not yet fixed
- Wall/anchor shape `sin(2*pi*x/(2L))` turns NEGATIVE past the outlet valve root, so the outlet valve's wall and anchors
  move in ANTI-phase (up to ~2 cells). The paper clamps the outlet valve to a stationary extension and uses
  `r = R - A/2 (1-cos(2*pi*x/L)) sin(wt)`. If the sample run still disagrees, fix this next.
- The phase convention is unchanged: the solver starts by expanding, and tau_paper = (dump-51)/100.

### Sample run — `Simulation_wallvel/pressureDiff_0`, launched 2026-09-23 16:05 (local)
Gamma 0.2 (T = 94,248), K 0.065 derived (leaflet LS_ANGLE_STIFFNESS 12.51), L = 125 (nX = 200), tau = 1, 2 cycles.
Compare with `python Template_v2/CompareOldData.py <run> [stride]` (it writes compare_old.csv/png in the run folder). Baseline
SimulationPaper gives Q/Q0 <= 0.05 and delta/DeltaV^(1/3) <= 0.16. The old targets are ~1.5 and ~0.40.

### Multithreading (2026-09-23) — 4.5x faster, results identical to serial
The solver used only ~3 of 20 logical cores. Profiling with timers around each phase gave **0.73 of 0.78 s per step in
`bcFindBoundaryNodes`**, a serial loop over all 7,684 surface triangles that searches for the lattice links each one cuts. Actually applying the
bounce-back took 0.015 s, and collision/propagation took 0.03 s.
Changes (Template_v2/source; backups `*.bak_before_lsomp`, `clMethod.h.bak_before_bcomp`):
1. `BC_USE_OPENMP` (new, in ConstDef.h): the boundary search is now two-phase. The read-only geometric search runs in parallel and
   records each triangle's cut links. Then `bcSetBoundaryNode` is applied serially in the original triangle order, since the first
   triangle to reach a link claims it (Bflag) and the apply step accumulates into shared arrays. Serial order means the result is exact.
2. `LS_USE_OPENMP` enabled. Its old branch was broken in three ways, all fixed:
   - `solveLSspring`'s OpenMP branch had its force accumulation commented out, so it would have dropped ALL spring forces.
     It is rewritten in gather form (each node sums its own springs, using the lower-index row for K/L0/Ds as the serial loop does).
   - `solveLSangle`'s OpenMP branch added into `vls.F` (which is overwritten later) instead of `vls.F.s`, so it would have dropped all bending.
   - `solveVerlet2` used an invalid `reduction(> : maxV)`. maxV is now taken in a serial pass after the loop.
   - `fill_s`/`fill_sf` stay serial (they accumulate flux corrections into shared LB cells) under their own `LS_FILL_USE_OPENMP`.
Verified by running serial and parallel builds on the same case for 300 steps: lsc and lsfb are bit-identical, and lsv/lbj/lbro differ by <=1e-10.
Speed: 70 s -> 15.5 s per 100 steps (~5,100 -> ~23,000 steps/h). One Gamma=0.2 cycle now takes ~4 h locally.

### Finding 5 (ROOT CAUSE #2) — the vessel wall was porous: scrambled and one-sided surface mesh
Even with the wall velocity fixed, the fluid next to the wall moved at ~2% of wall speed and inflow stayed ~50x too small.
Instrumenting the bounce-back showed the correct wall velocity reaching `bcSetNodeBC`, but only ~65k cut links per step
instead of the ~130k a two-sided shell should give. The cause was in the MATLAB wall mesh (`lsbl`):
1. `generateStraightPipeGeometry.m` built the nodes with `TH(:)` (column-major, so AXIS-fastest) while the connectivity assumes
   node `(iz-1)*ntheta+it` (ring-fastest). Every wall triangle joined the wrong nodes: 4,202/6,000 had zero area, some were huge,
   and the triangulated wall covered only 58% of the cylinder. The wall spring/bending network (built from the same faces in
   `pipeFromGeometry.m`) is scrambled the same way. That doesn't matter while the wall motion is prescribed, but it would if the wall were ever made passive.
2. The wall triangles were listed in ONE winding only (mostly outward normals). The solver bounces back only links that arrive on
   the side the normal faces, so fluid INSIDE the vessel passed through the wall. The leaflets (leafletCrescent) and the solver's own
   `clLS::testflatBL` both list thin sheets in both windings.
Fix (`Template_v2/geometry creator`, backups `*.bak_before_meshfix`): nodes are generated ring-major (`THr = TH.'`), and
`pipeFromGeometry.m` appends the reversed-winding copy of every wall triangle (and duplicates INOUT to match).
Check: 12,000 wall triangles, every one of area 2.111, total 25,338 = 2 x cylinder area, 6,000 out + 6,000 in.
Result (first dumps): fluid u_r at the wall is 0.95-1.03 x wall speed, net inflow is ~0.8 Q0 (was ~0.03), and the outlet valve closes/inlet
valve opens during expansion. The old leaky run is kept as `Simulation_wallvel_leakywall` for reference.
Cost: twice the triangles gives ~2 min/dump, ~28k steps/h, ~3.4 h per Gamma=0.2 cycle.
**Implication: every run made with the revived MATLAB generator (SimulationPaper, Simulation_v2, the L4 PACE run, and
Fig_ValveFlow's data) had a porous wall. Their flow fields are not physical.**

### First full cycle with both fixes (legacy sin wall law) — pumping works, but the valves are kinematic
`Simulation_wallvel` cycle 1: Q peaks at about +-1 Q0 (right sign and phase), but there is large backflow (inlet -0.93 Q0 during contraction,
outlet -0.70 Q0 during filling), so net Q is only ~0.14 Q0 against the old ~0.93. delta stays in 0.04-0.145 and follows the WALL POSITION,
not the flow: the anchors carry the leaflets.

### Finding 6 — wall law differed from the paper (`WALL_LAW_PAPER`, clMethod.h `wallMod`/`wallRate`)
The legacy law `A sin(2 pi x/2L) sin(wt)` moved the outlet-valve region in ANTI-phase and moved the inlet anchors at nearly full amplitude.
The paper uses `r = R - A/2 (1-cos(2 pi x/L)) sin(wt)` on [0,L] only: the outlet valve sits on a stationary extension, and tau = t/T has contraction first.
Implemented behind `#define WALL_LAW_PAPER` (on in Template_v2/ConstDef.h; backup `clMethod.h.bak_before_paperwall`).
`CompareOldData.py <run> <stride> paper` uses tau = t/T. Result (`Simulation_paperwall`, first 0.27 cycle): the inlet seals within 0.1 cycle
(backflow ~0) and the outlet now passes the ejection flow (up to 0.6 Q0), BUT delta_out only goes 0.120 -> 0.138 (old data ~0.40), because the leaflets barely bend.

### Finding 7 (ROOT CAUSE #3) — leaflets were 15x too stiff for the stated K
`buildConfig.m` computed the paper's k_b = 4 D_b/(3 sqrt 3) correctly, then multiplied it by v1's uniform gains (`gainAK = 15`,
`gainK = 7.5`, kept on purpose so v2 would reproduce v1 exactly). The solver applies weight*constant with no further factor, so the leaflet
bending constant applied was 12.5 against the paper's 0.834 for K = 0.065. **Every "K = 0.065" run so far was really K ~ 0.98.**
Fix: leaflet constants no longer carry the gains (`paper.legacyLeafletGain = false`, backup `buildConfig.m.bak_before_gainfix`).
For Gamma 0.2 and K 0.065 it now gives LS_ANGLE_STIFFNESS_1/2 = 0.83395 and LS_STIFFNESS_1/2 = 2.50185.
Runs going (10 threads each, `NUM_THREADS 10`): `Simulation_paperK` (paper-correct), `Simulation_paperwall_soft0.1` (k_b 1.25, 1.5x paper).
Reference: `Simulation_paperwall` (k_b 12.5, 28 dumps).

### Results at the paper-correct stiffness (`Simulation_paperK`, k_b 0.834) — better, but a mass leak remains
- Outlet opening delta_out peaks at 0.214 (0.138 at the 15x stiffness). The old data give ~0.40. The inlet seals (backflow ~0).
- 1.5x stiffer (`soft0.1`) gives 0.209, so the opening is NOT stiffness-limited any more. The leaflets are not stretch-limited either
  (strain 0.2%) and are nowhere near the wall. They simply see a small load, because the FLOW is ~half of what the kinematics require:
  Q_out ~0.73 Q0 at tau 0.06, against dV/dt ~1.5 Q0 (and the old data's 1.45).
- Control-volume budget over x in [30,193], dumps 5-9: the vessel volume falls at 0.68/step, outflow is 0.37, compression is 0.06,
  and **~0.25/step (37%) is unaccounted for, consistently** -> mass leaves through the wall region.

### Finding 8 — moving thin membrane with fluid on both sides: swept nodes carry mass across the wall
- The mask (`lbm.txt`) is uniformly 1: inside and outside the vessel are the same fluid. `inifill_map` seeds fluid on the side of each
  triangle's normal. With the two-sided wall (Finding 5), and with the leaflet triangles that extend outside the wall (anchors at r up to 14.4),
  the whole domain is seeded as fluid.
- The mask is STATIC. Only volumetric solids (`fill_s`, Sf tetrahedra) update it each step. As the wall sweeps across lattice nodes, those
  nodes switch sides WITH their fluid, so interior mass crosses the membrane instead of being pushed out axially.
- The solver HAS machinery for this, but all of it is disabled in this lineage:
  - the "fix for capsules" (inside/outside map ids BF vs BF+1 by normal side, in bcSetBoundaryNode, commented out);
  - `bcTestCrossSolid` (detects a node crossing the membrane, commented out);
  - `bcSetNewF` (refills a crossed node), which starts with `return;` and the comment "Wenbin said it is inaccurate, and it also seems to have made code not run on cluster".
- An exterior that is simply SOLID won't work either: the lumen expands by up to A = 3 cells into the static solid.
- **Open decision (asked Amir 2026-09-24):**
  (a) implement explicit mass redistribution for swept nodes, verified with a valve-less contracting-tube mass-balance test;
  (b) re-enable and repair the original capsule/crossing/refill path;
  (c) ask Alexeev/Karimi how the lost paper code handled the moving membrane.

### Mass-conservation study (2026-09-24) — `Template_v2/MassBudget.py` + `PlotMassBudget.py` -> `<run>/MassBudget.png`
Run: `Simulation_paperK` (2 cycles). The per-slice lumen budget is leak = -(d/dt[pi R^2 (1+rho')] + dQ/dx).
- Cycle 2: ~46% of the volume the wall squeezes out passes THROUGH THE WALL instead of out of the ends.
  It sloshes (out during contraction, back in during filling), so over the cycle the net is ~0 and the cycle-mean Q(x) is flat
  at ~0.26 Q0 (old data ~0.93). Inlet and outlet planes balance per cycle (0.52 vs 0.49 dV).
- Where, as |leak| per cycle in units of dV: inlet-valve footprint 0.40, outlet-valve footprint 0.32, lymphangion middle 0.16,
  static tube ends 0.05. Per unit length the valve footprints leak ~10x more than the middle. **The outlet footprint's wall is
  static under WALL_LAW_PAPER, so the leak there comes from the leaflet-wall junction, not wall motion.** (Leaflet anchor rows
  sit at r 9.5-14.4, so the leaflet sheet pierces the wall.)
- Static tube sections: a pressure-driven permeability, outward flux = 4.7e-4 x (rho_in - rho_out).
- Global: the domain mass swings by about +-4000 per cycle. The internal source (whole cross-section) is small except at the valve footprints and
  at the Zou-He planes (x = 0-2, 197-199), so the mismatch with the boundary flux in panel (f) is mostly pressure-BC plane accounting.
- NEXT: inspect the leaflet-wall junction (links that cross both a leaflet and a wall triangle, Bflag first-claim, anchor
  triangles outside the wall), and the half-way bounce-back test (`Simulation_test_hwbb`, set up but not run).

### Finding 9 (ROOT CAUSE #4, FIXED 2026-09-24) — the interpolated bounce-back leaked through the wall
- Link audit (`Template_v2/LinkLeakAudit.py`, debug dump of Bflag at step 300): all 161,112 lattice links that cross the wall
  or a leaflet ARE flagged. There are no holes in the surfaces.
- So the leak was in how the flagged links were treated. Bouzidi interpolation (for q < 0.5 it reads the population at the node
  BEHIND ii0) takes no account of a second surface between those nodes, which happens at the leaflet-wall junctions. It also
  runs uncorrected: `LB_NO_MASS_CORRECTION` is defined.
- Tested over the same window (t/T 0-0.12, `Template_v2/LeakByRegion.py`), wall leak as a share of the squeezed volume:
  interpolated BB 35% | mass correction re-enabled 2.3% | **half-way BB (`LB_HW_BB`) 1.6%**. With half-way BB the valve footprints leak 0.002 dV (was 0.05).
- **Fix: `#define LB_HW_BB` in Template_v2/source/ConstDef.h** (backup `ConstDef.h.bak_before_hwbb`). The outflow rises ~45% at the same phase
  (Q_out 1.05 Q0 at tau 0.09, was 0.73). The residual gap to the kinematic value is compression, which Amir says to ignore.
- The outlet opening is unchanged (delta_out 0.219 vs 0.214), so under-opening of the valves is a SEPARATE issue from the leak.
- Figure: `LeakFix_before_after.png`. The run continuing to 2 cycles is `Simulation_test_hwbb` (10 threads, started 09:3x).
- Stopped: `Simulation_test_masscor` (14 dumps kept).

### Valve under-opening after the leak fix (2026-09-24) — what it is NOT
All runs use the leak-fixed config (LB_HW_BB, WALL_LAW_PAPER, K = 0.065 -> k_b 0.83395, k_s 2.50185). Metric: delta_out at tau 0.09-0.12.
Baseline `Simulation_test_hwbb`: 0.219 (old data ~0.39-0.40). Flow now matches the old data (Q_out peak 1.45 vs 1.47).
| test | delta_out | effect |
|---|---|---|
| bending x1/3 (`test_bend0.33`) | 0.226 | +3% |
| stretching x0.1 (`test_stretch0.1`) | 0.254 | +16% |
| clamp band 76 -> 40 anchors (`test_band2rows`, TrimLeafletAnchors.py band 1.0) | 0.222 | none |
| clamp band 76 -> 26 anchors (`test_band1row`, band 0.5) | 0.222 | none |
| leaflet density 1000 -> 100 / 30 (`LS_MASS_1/2` 0.1 / 0.03) | 0.219 / 0.219 | none |
| leaflet density 10 / 3 / 1 (`LS_MASS_1/2` 0.01 / 0.003 / 0.001) | - | UNSTABLE (explicit FSI added-mass), crash in minutes |
- lsk/lsak audit: load/ weights x ConstDef constants = exactly the intended k_b/k_s. Written once at startup, never modified at runtime
  (clK::ini / clAK::ini only), no slaved pairs, no rigid bodies, no static/frozen leaflet nodes (bit definitions match MATLAB).
- Leaflet mass: GeometryCreatorAmir.m has Ro0 = 1000 against fluid 1 (lattice units), so nodes weigh 240-490. This doesn't matter (quasi-static).
- Free edge: 23.53 long over a 20.00 commissure chord (3.53 excess). An inextensible edge alone would allow delta ~0.36,
  but at peak the edge strains only 1.2% (2.2% with stretching x0.1). The opening is limited by the in-plane (membrane) constraint
  of a sheet clamped along a curved line: it cannot bow without stretching. Stretching is the only lever that responded.
- Metric caveat: after trimming anchors, CompareOldData's free-edge detection picks up the new cut boundary -> use the original crescent
  nodes (matched by position) for trimmed meshes.

## Session 2026-09-24 .. 26 — commissure fix, continuity, hinge bending  (READ THIS FIRST in a new session)

### State at hand-off (2026-09-25 ~23:55)
**UPDATE 2026-09-26 00:36: the running case is now `Simulation_hinge_wallsoft_edgestretch0.1` (see Change 15); `Simulation_hinge_wallsoft` and `..._edgestretch` (0.2) were stopped.**
- RUNNING (detached from Claude, survives session end): `Simulation_hinge_wallsoft/pressureDiff_0` (LBLS.exe, started 23:50, 2 cycles,
  ~2.5 min/dump on 20 threads, ~8 h total). Launched via its `run_detached.cmd` with Win32_Process.Create, so it is NOT a child of the
  Claude session. Stop it with: `taskkill /IM LBLS.exe /F`.
- RUNNING (detached): the live valve dashboard updater `Template_v2/LiveValveDashboard.py` -> `LiveDashboard/dashboard.png|html`
  (log `LiveDashboard.log`); it exits by itself when the solver stops. To view it, serve the folder:
  `python -m http.server 8765 --directory LiveDashboard`, then open http://localhost:8765/dashboard.html
  (`French Collaboration/.claude/launch.json` has this as configuration "live-dashboard").
- 2026-09-26 00:03: updater restarted (now PID 34640, log `LiveDashboard/updater.log|.err`) with a tabbed large-image view:
  it also writes `inlet.png`, `outlet.png`, `history.png` + `tabs.html` -> open http://localhost:8765/tabs.html
  (buttons: Inlet valve / Outlet valve / Opening & flow; images refresh every 20 s). The old dashboard.html still works.
  Pre-change copy of the script: `Template_v2/LiveValveDashboard.py.bak`.
- Closing the laptop lid sleeps the machine (the run pauses and continues on wake) unless the Windows lid action on AC is "Do nothing"
  (`powercfg /setacvalueindex SCHEME_CURRENT SUB_BUTTONS LIDACTION 0` + `powercfg /setactive SCHEME_CURRENT`).
- There are NO checkpoints: dumps hold rho, J, lsc, lsv only (no distribution functions), so a stopped run can only restart from t = 0.

### Current solver = Template_v2 (switches in source/ConstDef.h)
`LB_HW_BB` (half-way bounce-back; fixes the wall leak), `WALL_LAW_PAPER`, `LS_USE_OPENMP` + `BC_USE_OPENMP` (4.5x faster),
`LS_HINGE_BENDING` (new, below), `LS_HINGE_WALL_FACTOR 1.0` / `LS_HINGE_WALL_RAMP 3.0` (wall softening, OFF by default),
`LB_NO_MASS_CORRECTION` (still defined). Geometry: the commissure fix in `geometry creator` (below).
buildConfig: L_nd 6.25, Gamma 0.2, K 0.065 (k_b 0.834, k_s 2.50); `cfg.simulationFolderName` is currently 'Simulation_commissure_v2'
-> change it before run_sweep, or that reference folder gets overwritten.

### Finding 10 — the commissure gap (FIXED; Amir confirmed it was the leak cause)
In `leafletCrescent.m` the crescent cut DELETED the two outermost grid rows at the free end (and shifted rows 2 / Nrows-1), so the free edge
stopped at |y| = 9.1 < R. `GeometryCreatorAmir.m` also placed each leaflet 0.61 off the mid-plane. Result: an open wedge at each commissure;
~90% of the closed-valve backflow went through it (`Template_v2/CommissureGapFlow.py`).
Fix (backups `*.bak_before_commissure`): outer rows are trimmed to the crescent like every other row, and `snapToMidplane()` moves each leaflet's
point nearest the mid-plane onto it. The upper and lower corner nodes now meet on the wall at (|y| = R, z = 0).
Flow through the commissure regions: -0.020 -> 0.000 (outlet), -0.004 -> ~0 (inlet). Closed inlet delta 0.008 (was 0.040).
The closed outlet keeps a thin central slit (delta 0.036, ~ -0.035 Q0): the outlet leaflets do not fully meet (the under-deflection issue).
Run `Simulation_commissure_v2` (147 dumps, 1.46 cycles) = the TRIPLET-bending reference used below.
Figures: `GapFix_compare.png`, `GapFix_zoom.png`, `CommissureFix_rest.png`, `SideView_widest_open.png`.

### Finding 11 — continuity: the moving wall creates mass (NOT fixed; Amir: do not worry about compression)
`Template_v2/ContinuityCheck.py` over the last full period of commissure_v2: IN 0.673 dV, OUT 0.742 dV, storage -0.002, so +0.067 dV per cycle
appears through the wall, almost all along the MOVING lymphangion wall (x 64-162) and only during expansion (tau 0.35-0.65). The flow through the whole
cross-section also grows along x, so this is mass CREATED, not an exchange with the outside. Figure `Continuity_commissure.png`.
Test `Simulation_test_masscorr_hwbb` (mass correction re-enabled, 47 dumps): the expansion-phase source is IDENTICAL and flow drops ~20%,
so global mass correction is NOT the fix. A targeted correction for nodes crossing the moving wall would be needed (not done).

### Finding 12 — the valve under-opening is NOT from bad lsk/lsak, the clamp band or mass
- lsk/lsak audit: load weights x ConstDef constants = exactly the intended k_b/k_s; set once in clK::ini/clAK::ini, never modified at runtime.
- Bending x1/3: +3%. Stretching x0.1: +16%. Clamp band 76 -> 40/26 anchors (`TrimLeafletAnchors.py`): 0%. Density 1000 -> 100/30: 0%.
  Density 10/3/1 crashes: an explicit-FSI added-mass instability that then runs out of memory (0xc0000409, "paging file too small").
- Opening at K 0.065: inlet 0.30, outlet 0.23 (old data ~0.40). Timing and flows match the old data well (Q_out peak 1.45 vs 1.47).
- Amir's position (2026-09-25): geometry is correct; the remaining difference is HOW the lost code modelled the valve stiffness.
- How stiffness is modelled now (explained in chat): in-plane = Hooke springs k_s on the 6-neighbour triangular lattice (E2D = 2 k_s / sqrt(3),
  nu = 1/3, damping Ds = 0.3); bending = k_b (1 + cos theta) on straight node triplets (D_b = 3 sqrt(3)/4 k_b); k_s is tied to k_b by the old
  "EEv = EIv" rule, giving E2D = 2.89, i.e. ~180x softer in-plane than a real 0.15-thick plate.

### Change 13 — leaflet bending by dihedral HINGES (Amir's request, 2026-09-25)
`LS_HINGE_BENDING` in `clMethod.h` (buildHinges / solveLShinges; backup `clMethod.h.bak_before_hinge`). Every interior edge shared by two leaflet
triangles is a hinge (tips i, l; edge j-k). E = 0.5 kh (theta - theta0)^2, theta = the signed dihedral angle between the triangle normals, theta0 from
the reference mesh. kh = 1.5 * LS_ANGLE_STIFFNESS_1/2, which gives the same D_b (sheet: D_b = sqrt(3)/2 kh; triplets: D_b = 3 sqrt(3)/4 k_b), so
K keeps its meaning. Leaflet straight triplets are skipped; the wall keeps triplets. 1768 hinges.
Verified in `Template_v2/HingeBendingCheck.py`: gradient vs finite differences 1e-9 (the first version had a wrong sign on the edge-node terms,
now fixed); cylinder-bending energy on the real leaflet mesh: hinge/plate 0.94-0.99 (triplet 0.87-1.00).
Run `Simulation_hinge` (11 dumps, stopped): early delta/Q within a few % of the triplet reference.

### Change 14 — wall-proximity hinge softening (Amir's request, 2026-09-25; RUNNING)
kh *= f, f = LS_HINGE_WALL_FACTOR + (1 - FACTOR) * min(dw / LS_HINGE_WALL_RAMP, 1), where dw = distance of the hinge-edge midpoint from the
wall in the reference mesh (backup `clMethod.h.bak_before_wallsoft`). Running case `Simulation_hinge_wallsoft`: FACTOR 0.2, RAMP 3.0
(the solver reports factor min 0.200, mean 0.597). These values were chosen by Claude (Amir gave none) - revisit.
NEXT: compare delta_in/out (peak inlet ~tau 0.39, outlet ~tau 1.0) with commissure_v2 (0.30 / 0.23) and the old data (~0.40):
`python Template_v2/CompareOldData.py Simulation_hinge_wallsoft/pressureDiff_0 1 paper`.

### Change 15 — free-edge softening of the leaflet in-plane springs (Amir's request, 2026-09-26; RUNNING)
- 2026-09-26 00:10: `Simulation_hinge_wallsoft` STOPPED by Amir's request at dump ~10 (t/T ~0.07); superseded by the run below.
- `LS_STRETCH_EDGE_FACTOR` / `LS_STRETCH_EDGE_RAMP` in ConstDef.h (template default 1.0 = off). `clMethod.h` buildEdgeSoftening(), called once
  at the first solveLSspring (backups `clMethod.h|ConstDef.h.bak_before_edgestretch`): every leaflet spring K *= f,
  f = FACTOR + (1 - FACTOR) * min(de / RAMP, 1), de = distance of the spring midpoint from the leaflet's free edge (reference mesh).
  Free edge = leaflet boundary edges with at least one end outside the clamp band (lst 256/512) = 12 segments / 13 nodes per leaflet.
  Damping Ds unchanged. Note `lsk.txt` is written BEFORE the scaling (still shows unscaled weights); the scaled factors are in `lsedgesoft.txt`
  (i, j, de, f) in the run folder.
- Amir's values: FACTOR 0.2 (5x softer at the free edge), RAMP 24 lu (= whole leaflet: 24.1 lu is the farthest node from the free edge),
  hinge wall softening kept (0.2 over 3 lu). Solver: 48 free-edge segments, 1976 leaflet springs (all 494 per leaflet), factor min 0.200, mean 0.584.
  Checked against an independent Python calculation from the load files: identical to print precision (5e-5).
- Run `Simulation_hinge_wallsoft_edgestretch/pressureDiff_0`, started 2026-09-26 00:14 (detached, Win32_Process.Create of run_detached.cmd,
  PID 6896, 20 threads, ~8 h for 2 cycles). Live dashboard updater restarted on it (PID 41972): http://localhost:8765/tabs.html
- 2026-09-26 ~00:33: that run (FACTOR 0.2) STOPPED by Amir at 11 dumps; he wants the edge 10x softer.
- Run `Simulation_hinge_wallsoft_edgestretch0.1/pressureDiff_0`: identical except LS_STRETCH_EDGE_FACTOR 0.1 (ramp 24 lu, hinge wall
  softening 0.2/3 kept). Started 2026-09-26 00:36 (detached; the first Win32_Process.Create right after MSBuild silently did nothing -
  relaunch worked). Solver: factor min 0.100, mean 0.532 over 1976 springs. Dashboard updater restarted on it.
- RESULT (finished 07:06, 200 dumps, 2 cycles, err.log empty; `compare_old.csv|png` in the run folder):
  inlet peak 0.256 @ tau 0.47 (both cycles) vs ref 0.302 @ 0.39 -> 15% LOWER and later; outlet peak 0.268 (cycle 2) vs ref 0.237 -> +13%;
  Q_out max 1.52 (old 1.47), Q_in max 1.31 (old 1.49). Old data: in 0.42, out 0.40. Still well short; the edge softening helps the outlet
  about as much as uniform stretch x0.1 did (+16%) but reduces the inlet peak.
- NEXT: decide with Amir (e.g. ramp shorter than 24 lu, uniform stretch softening, or other stiffness conventions).

### Change 16 — half period (Amir's request, 2026-09-26; RUNNING)
- NOTE on T: `CYCLE_PERIOD` is 2T (run length); the wall period is RAMP_PERIOD = CYCLE_PERIOD/RAMP_PER_CYCLE(2). Full-T runs: T = 94248
  (= 2 pi R^2/(nu Gamma^2), nu = 1/6, Gamma 0.2). Wall: R 10, A 3, r = 7..13 at mid-lymphangion, L 125 (roots x 38/163), 40x40 section.
- `Simulation_edgestretch0.1_halfT`: copy of `Simulation_hinge_wallsoft_edgestretch0.1` with ONLY `CYCLE_PERIOD 94248` -> T = 47124,
  dumps every 471 steps (T/100), 200 dumps = 2 cycles (~3.3 h). Raw stiffness constants unchanged, so effectively Gamma -> 0.2*sqrt(2) = 0.283
  and the leaflets are relatively softer against the (2x larger) flow forces. Started 09:23, PID 35808. Live dashboard: half T vs full T
  (updater PID 33848) at http://localhost:8765/tabs.html.

### Finding 17 — the whole fluid breathes: the open ends do not hold the pressure (2026-09-26; NOT fixed yet)
`Template_v2/CompressibilityCheck.py <out.png> <run> <label> [...]` -> `Compressibility_fullT_vs_halfT.png`.
- Mach is small (max 0.018), yet rho' swings about -0.065..+0.082 over each cycle, UNIFORMLY over the whole domain, including x = 0 and x = 199
  (dump 100: rho' 0.062 at both ends and 0.077 mid). Local gradients are small next to this global swing.
- Mass budget (full T): compression storage d/dt Int A rho' reaches 1.2 Q0 against |dV/dt| max 1.6 Q0; the mass stored by compression swings
  0.31 stroke volumes over a cycle -> about a third of the wall's stroke goes into squeezing the fluid instead of through the valves.
  Likely the same thing as Finding 11 ("moving wall creates mass").
- Why the ends float: `lsinout.txt` is all zeros, so INLET_P_BC/OUTLET_P_BC (bcSetNodeBC) never activate. The x-ends use bcPressLBx0/x1
  (Zou-He, rho target 1 + RoX0*invel = 1). In bcPressLBx0 (clMethod.h ~2416) `if(Bflag & LF[k]) return;` exits the WHOLE function at the first
  node whose link was already set by bounce-back (the tube-edge nodes come first in the j/z loop), so almost no end nodes get the pressure BC.
  Should probably be `continue`. Also `refsum = ...` should be `refsum += ...`. Same pattern to check in bcPressLBx1. Not yet confirmed by a test run.
- 2026-09-26 09:43: tried that fix (return->continue, refsum += with C[m], both ends) in `Simulation_edgestretch0.1_pressfix`: the run wrote
  nothing in 6 min (stuck or blew up at start; stdout was buffered so no clue where). Stopped, and on Amir's request the fix was REVERTED:
  Template_v2 and the pressfix case's clMethod.h are again identical to the pre-fix solver. Still in place: unbuffered stdout
  (setvbuf in LBLS.cpp), CompressibilityCheck.py live mode, Compressibility tab in LiveValveDashboard.py. Half-T run stopped at 15 dumps.

### Change 18 — leaflet base stiffness x0.1 (Amir: "these are the correct stiffness values"), 2026-09-26
- `Simulation_edgestretch0.1_base0.1`: copy of the edge-stretch x0.1 full-T run with ONLY LS_STIFFNESS_1/2 2.501851166 -> 0.2501851166 and
  LS_ANGLE_STIFFNESS_1/2 0.8339503888 -> 0.08339503888 (wall, softenings, T unchanged; pressure-BC fix NOT included). Started 09:56,
  stopped by Amir at 13 dumps (tau 0..0.12). compare_old.csv in the run folder.
- Non-dimensional numbers (paper definitions, Main_Rev7 / Downloads PDF): Gamma 0.200, Re 0.0365, A/2R 0.15, L/2R 6.25, ell 1.25 all in the
  paper's range; K = D_b/(2 pi R^3 P0) = 0.0065 (full), 0.0013 at the wall -> 10x BELOW the paper's K = 0.065. P0 0.002653, Q0 0.5, D_b 0.1083.
  Stretching (no paper group): E2D 0.289 (full), 0.029 at the free edge; E/(P0 R) = 10.9 / 1.1.
- vs paper data (LymphaticsDataNew.xlsx, K 0.065, Gamma 0.2) at the same tau: outlet opening 0.376 at tau 0.12 vs paper 0.384 (full-stiffness
  run 0.249 at the same tau). Caveat: the sim starts from rest (cycle 1) while the paper's curve is the periodic state, so tau < ~0.1 mixes in
  start-up (paper outlet already 0.40 at tau 0; sim Q_out rises from 0). Inlet: sim keeps a ~0.025 slit where the paper is fully shut (0.000),
  and sim Q_in is -0.4..-0.05 Q0 (backflow) where the paper's is ~0.

### Change 19 — closed-valve gaps: view fix + commissure ties (Amir chose A+B), 2026-09-26
- Full-stiffness run `Simulation_hinge_wallsoft_edgestretch0.1` moved to the Recycle Bin on Amir's request (Amir: the x0.1 base stiffness is correct).
- Diagnosis (base0.1 run, dump 13, inlet closed): centre seals (free edges at z = +-0.02); slits from |y| ~6 outward, 0.32 lu at |y| 6.6 and
  0.68 at |y| 8, closing at the corner node. The white slivers next to the wall in the end-on view were an artefact: the circle was drawn at one x,
  but the leaflets are clamped up to x ~64 where the wall has contracted to r 9.25. Real open area 3.3 lu^2 (1.2% of lumen); flow through the
  |y|>6 band ~0.2% Q0. Figure `Simulation_edgestretch0.1_base0.1/pressureDiff_0/closed_gap_check.png`.
- A: `ValveCrossSection.py` endon() draws the wall at the NARROWEST radius along the valve's x-span (wallRmin); used by all end-on panels.
- B: `clMethod.h` solveLSties()/buildCommissureTies() (backups *.bak_before_committie): zero-length spring K (LS_COMMISSURE_TIE_K, 0 = off)
  + damping D on relative velocity between the outermost unclamped free-edge nodes of opposite leaflets (1-2, 3-4), both sides -> 4 ties at
  |y-yc| 8.33, rest separation 2.09.
- Run `Simulation_base0.1_commtie` = base0.1 run + K 5, D 5. Started 10:40. Ties hold (separation 0.026 at dump 3). Peak strain spiked to 1.13
  during the initial pull-in and sits at ~0.38, but ONLY in the 16 springs touching tie nodes (median strain ~0); a local distortion.
  Live: tabs.html (tie vs no tie) + compressibility tab.
- Amir stopped it at 7 dumps: the permanent tie pinches the OPEN orifice (outlet open area -30..37%, delta_out 0.261 vs 0.298 at tau 0.06).
  The tie node is already the last free-edge node (1.7 lu from the wall), so a permanent tie cannot be made smaller on this mesh.
  Untied pair separation: 0.67 lu closed vs 2.7-3.9 lu open -> CLOSING-ONLY ties: LS_COMMISSURE_TIE_OPEN (lu); tie scaled by
  w = clamp((OPEN - dc)/(OPEN/2), 0, 1), dc = distance between the valve's two central free-edge nodes (full once the centre has shut, off
  once it opens past OPEN). Starts inactive (rest centre gap 3.8), so no pull-in snap.
- Run `Simulation_base0.1_commtie_closing` = commtie + OPEN 0.5. Started 11:02; early strain 0.027 (no snap). Live vs the untied run.
  Result to dump 10: ties engage at tau 0.04 (pair gap 0.039); closed-inlet open area 1.5 vs 4.4 lu^2 untied (-65%); outlet delta and open
  area IDENTICAL to the untied run (no blocking). Stopped by Amir at 16 dumps. (29-min stall between dumps 10-11 = machine asleep;
  OneDrive syncing the dumps also eats CPU.)

### Other things done this session
- Deleted the one byte-identical duplicate `French Collaboration/Fig1 (1).png` (other look-alike files differ; kept).
- Audit of `French Collaboration/main3.tex`: several captions/claims inconsistent with the figures (see chat); nothing edited yet.

### Tools added (all in Template_v2/; usage in each docstring)
CompareOldData.py (delta/Q vs the lost-code data; 'paper' phase arg) | AnimateRun.py + makeMp4.m (animations) | MassBudget.py + PlotMassBudget.py +
LeakByRegion.py + PlotLeakCompare.py (mass budget / wall leak) | ContinuityCheck.py (IN vs OUT over the last period) | BudgetMonitor.py (live residual) |
LinkLeakAudit.py (bounce-back link audit; needs a debug dump) | TrimLeafletAnchors.py | ValveCrossSection.py (end-on views; anim/grid) |
SideViewOpen.py | CommissureGapFlow.py | HingeBendingCheck.py | LiveValveMonitor.py | LiveDashboard.py / LiveValveDashboard.py |
old_data_by_dump.csv (old data interpolated per dump).

### Open items / next steps
1. Evaluate `Simulation_hinge_wallsoft` against the triplet reference and the old data; tune FACTOR/RAMP with Amir.
2. The lost code's stiffness convention is unknown; the LymphaticsData.xlsx 'Kb' column (e.g. 106.1 at L/l=2) may encode its raw values - decode it.
3. Moving-wall mass source (Finding 11): a targeted swept-node correction, if it matters.
4. Checkpoint/restart (save F + LS state every N dumps) before running on PACE (both earlier PACE jobs hit the 48 h walltime).
5. The PoF reviewer figure `Fig_ValveFlow` was made from a run with the old porous wall - regenerate it before resubmission.

## Contact Info
User: Amir (amirp11.ap@gmail.com)
Folder: `C:\Users\Amir\OneDrive - Georgia Institute of Technology\Lymphatics\PhaseDifferenceStudy\`
