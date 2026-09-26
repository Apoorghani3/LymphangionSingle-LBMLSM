# GeometryCreatorAmir.m — Complete Documentation

## Overview
`GeometryCreatorAmir()` generates 3D LBM geometry for lymphatic vessel simulations with leaflets. It creates:
- **Vessel geometry**: cylindrical pipe with specified diameter and length
- **Leaflet geometry**: crescent-shaped valve leaflets at specified positions
- **Static regions**: pinned vessel wall nodes at entrance/exit
- **Output files**: saved to `load/` in ASCII format for C++ solver

---

## Function Signature

```matlab
GeometryCreatorAmir(flShow, leaflets_flag, Vessel_flag, Leaflet_length, 
                    entrance_static_flag, num_valve, Domain_offset, 
                    Vessel_diam, Vessel_Length, a, b, Crescent_depth)
```

### Input Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `flShow` | bool | 1 → display 3D plot after generation |
| `leaflets_flag` | bool | 1 → create leaflet geometry |
| `Vessel_flag` | bool | 1 → create vessel geometry |
| `Leaflet_length` | double | leaflet length / vessel diameter ratio |
| `entrance_static_flag` | bool | 1 → pin vessel nodes at entrance/exit |
| `num_valve` | int | number of valve pairs |
| `Domain_offset` | int | LBM domain padding (lattice nodes) |
| `Vessel_diam` | double | vessel diameter (lattice units) |
| `Vessel_Length` | double | vessel length (lattice units) |
| **`a` (spacingA)** | double | **inter-valve gap multiplier** |
| **`b` (spacingB)** | double | **entrance/exit gap multiplier** |
| `Crescent_depth` | double | crescent cut depth for leaflets |

---

## Section 1: Vessel Geometry

**Lines 58-84** — Creates cylindrical pipe via `pipeFromGeometry()`:

- Domain dimensions: `nX × nY × nZ` (lattice nodes)
- All vessel nodes initialized with `O = 1` (object ID for vessel)
- Centerline computed: `Cent_cap = [nX/2, nY/2, nZ/2]`
- If vessel disabled: default `nX=140`, centered vessel position

**Output arrays after this section:**
- `Cx, Cy, Cz` — node coordinates
- `N, M, K, L0` — connectivity and spring data
- `T` — node type bitmask
- `O` — object ID (initially all 1 for vessel)

---

## Section 2: Leaflet Geometry

**Lines 89-222** — Creates and positions valve leaflets.

### Grid Sizing (Lines 91-112)

```matlab
Leq_min = 1.5;           % minimum element size
nZr = floor(Vessel_diam / Leq_min);
Leq = (2/√3) * Vessel_diam / (nZr - 1);
nXr = round((nZr - 1) * Leaflet_length);  % leaflet X extent
```

- `nZr`: circumferential resolution (leaflet height in Z)
- `nXr`: axial extent of leaflet in X direction
- `Angle3D`: lookup table maps `Leaflet_length` → tilt angles (pre-computed for 1–5)

### **Valve Positioning — Spacing-Based (Lines 129-149)**

#### **NEW FIX** ✓

Valves now positioned using `spacingA` and `spacingB` parameters:

```matlab
valve_width_x = nXr;                          % leaflet width in X
space_scale = Vessel_diam * Leaflet_length;   % scaling factor

first_valve_min_x = b/2 * space_scale;        % entrance gap
inter_valve_spacing = a * space_scale;        % gap between valves

for leafletc = 1:num_valve
    valve_min_x = first_valve_min_x + (leafletc - 1) * inter_valve_spacing;
    XcenterR = round(valve_min_x);            % center of this valve
```

#### Layout Interpretation:
```
[  spacingB/2    ] [ Valve 1 ] [  spacingA  ] [ Valve 2 ] [  spacingB/2  ]
|<-- entrance -->|           |<-- inter -->|           |<-- exit -->|
```

**Scaling rationalization:**
- `buildConfig` computes `Vessel_Length ∝ (spacingA × (nValves-1) + spacingB) × Vessel_diam × Leaflet_length`
- Spacing parameters must scale with `Vessel_diam × Leaflet_length` for consistency

### Leaflet Shape

Each valve pair consists of **two crescent leaflets** (top and bottom):

**Leaflet 1 (top, lines 150-164):**
- Generated via `leafletCrescent()` with angle `Angle3D`
- Positioned at X = `XcenterR`, Y/Z centered in vessel
- Placed in upper half of vessel (offset by `+Vessel_diam/2` in Z)

**Leaflet 2 (bottom, lines 166-181):**
- Same geometry but mirrored (angle reversed: `Angle3D(2) → -Angle3D(2)`)
- Placed in lower half of vessel (offset by `-Vessel_diam/2` in Z)

**Boundary tagging:**
- Nodes outside vessel radius → marked with `Tcustom0` bitmask (repulsive forces)

### Leaflet Appending (Lines 183-193)

Each leaflet appended to global arrays via `appendLeaflet()`:
- Object IDs increment: first valve leaflets get `O = 2, 3`; second valve `O = 4, 5`; etc.
- All coordinates, connectivity, and flags aggregated

### Optional Cylindrical Wrap (Lines 203-219)
- Disabled by default (`Leaflet_alter_flag = 0`)
- If enabled: wraps leaflet nodes to cylindrical surface for smooth curvature

---

## Section 3: Static Entrance/Exit Nodes

**Lines 224-253** — Pin vessel wall nodes outside valve region.

### **NEW FIX** ✓ — Linked to Valve Extents

```matlab
if entrance_static_flag == 1
    maxO = max(O);
    
    % First valve: O=2 and O=3
    first_valve_indices = find(O == 2 | O == 3);
    if ~isempty(first_valve_indices)
        firstValveMinX = min(Cx(first_valve_indices));
    else
        firstValveMinX = 0;
    end
    
    % Last valve: O=maxO and O=maxO-1
    last_valve_indices = find(O == maxO | O == maxO - 1);
    if ~isempty(last_valve_indices)
        lastValveMaxX = max(Cx(last_valve_indices));
    else
        lastValveMaxX = max(Cx);
    end
    
    % Mark static: O=1 vessel nodes outside [firstValveMinX, lastValveMaxX]
    for i = 1:length(O)
        if O(i) == 1 && (Cx(i) < firstValveMinX || Cx(i) > lastValveMaxX)
            T(i) = T(i) + Tstatic;
        end
    end
end
```

#### Behavior:
- **First valve bounds** determined by `min(X)` of `O=2` and `O=3` leaflets
- **Last valve bounds** determined by `max(X)` of last two leaflets
- **Static flag** (`Tstatic = 2^4 = 16`) added to type `T` for entrance/exit vessel nodes
- All valve leaflets remain **fluid-driven** (not pinned)

#### Advantage of dynamic linking:
- Static regions automatically adjust if `spacingA`, `spacingB`, or `Vessel_Length` change
- No hardcoded X values needed
- Always matches actual valve positions

---

## Section 4: Data Output

**Lines 263-301** — Saves geometry arrays to ASCII files in `load/` directory.

### Generated Files:

| File | Content | Description |
|------|---------|-------------|
| `lsn.txt` | `N - 1` | Node neighbor connectivity (1-indexed) |
| `lsm.txt` | `M` | Node mass |
| `lsk.txt` | `K` | Spring constants |
| `lsc.txt` | `[Cx, Cy, Cz]` | Node coordinates (3 columns) |
| `lsbn.txt` | `BN - 1` | Boundary node indices |
| `lsbs.txt` | `BS` | Boundary spring stiffness |
| `lsbf.txt` | `BF = BS + 1` | Boundary force flag |
| `lbsize.txt` | `[nX; nY; nZ]` | Domain dimensions |
| `lst.txt` | `T` | Node type bitmasks |
| `lso.txt` | `O - 1` | Object IDs (0-indexed: 0=vessel, 1,2=1st valve leaflets, etc.) |
| `lspo.txt` | `LSPO` | Pressure offset points |
| `lspn.txt` | `SN - 1` | Spring node pairs |
| `lspnc.txt` | `LSN` | Spring node coefficients |
| `lsbl.txt` | `lsbl - 1` | Surface triangle indices |
| `lsinout.txt` | `INOUT` | In/out flags for triangles |
| `lsan.txt` | `AN - 1` | Bending triplet indices |
| `lsak.txt` | `AK` | Bending spring constants |
| `lbm.txt` | `lbm (3D)` | Lattice Boltzmann map (1=fluid, 0=solid) |

---

## Type Bitmasks

**Lines 32-38** — Node type flags (bitwise-OR combination):

```matlab
Tflat          = 2^1  = 2    % flat element
Tstatic        = 2^4  = 16   % pinned/static
Twallrepulsive = 2^6  = 64   % wall repulsion
Touter         = 2^7  = 128  % outer boundary
Tcustom0       = 2^8  = 256  % custom marker 0 (leaflet outside radius)
Tcustom1       = 2^9  = 512  % custom marker 1 (custom marker 1)
```

A node's type `T(i)` can combine multiple flags. Example:
- Vessel wall node: `T(i) = 0` (not static, not outer, etc.)
- Static entrance vessel: `T(i) = 16` (Tstatic)
- Leaflet outside vessel: `T(i) = 256` (Tcustom0)

---

## Object IDs (O array)

| Value | Meaning |
|-------|---------|
| 0 | (reserved/unused) |
| 1 | Vessel wall |
| 2 | First valve, leaflet 1 (top) |
| 3 | First valve, leaflet 2 (bottom) |
| 4 | Second valve, leaflet 1 (top) |
| 5 | Second valve, leaflet 2 (bottom) |
| ... | (pattern continues for additional valves) |

**Note**: `appendLeaflet()` increments object ID after each valve pair.

---

## Visualization (Lines 305-324)

If `flShow = 1`, displays 3D scatter plot color-coded by node type/object:

- **Black filled circles** — static nodes (Tstatic)
- **Black X** — outer boundary nodes
- **Red filled circles** — leaflet nodes outside vessel radius (Tcustom0)
- **Blue circles** — vessel and certain leaflet nodes
- **Cyan circles** — other leaflet nodes
- **Green circles** — specific boundary marker

---

## Key Parameters in buildConfig.m

These settings drive GeometryCreatorAmir:

```matlab
cfg.geometry.spacingA = 10;       % inter-valve gap (multiplier)
cfg.geometry.spacingB = 10;       % entrance/exit gap (multiplier)
cfg.geometry.vesselDiameter = 20; % [lattice units]
cfg.geometry.leafletLength = 1.25; % ratio to vessel diameter

% Derived:
cfg.geometry.vesselLength = Vessel_diam * Leaflet_length 
                          * ((num_valves - 1) * spacingA + spacingB)
```

---

## Workflow: From Config to Geometry

1. **buildConfig.m** defines spacing parameters and computes vessel length
2. **generateGeometry.m** calls `GeometryCreatorAmir()` with these values
3. **GeometryCreatorAmir** generates geometry using:
   - `spacingA` and `spacingB` for valve positioning
   - Dynamic valve extent detection for static region linking
4. **Geometry files saved** to `load/` (synced to all simulation subfolders)
5. **C++ solver** reads from `load/` files for simulation

---

## Recent Changes Summary

### ✓ Valve Positioning Fix
- **Before**: Valves evenly distributed via `XcenterR = round(nX * leafletc / (num_valve + 1))` — ignored `a` and `b`
- **After**: Valves positioned using spacingA and spacingB:
  - Entrance gap: `b/2 × Vessel_diam × Leaflet_length`
  - Inter-valve gap: `a × Vessel_diam × Leaflet_length`
  - Scaling factor: `space_scale = Vessel_diam × Leaflet_length`

### ✓ Static Entrance/Exit Linking
- **Before**: Hardcoded static regions at X < 25 and X > (max(X) - 25)
- **After**: Dynamic linking to actual valve extents:
  - First valve min X from `O=2,3` nodes
  - Last valve max X from `O=maxO, maxO-1` nodes
  - Static vessel nodes automatically positioned outside this range

---

## Usage Notes

- Always run `generateGeometry(cfg)` after modifying spacing or vessel parameters
- The `load/` folder is synchronized to simulation subfolders by `setupSubfolders.m`
- To visualize geometry: set `cfg.geometry.showGeometry = true`
- For debugging: check `lso.txt` (object IDs) and `lsc.txt` (coordinates) to verify valve positions
- Valve widths depend on `Leaflet_length` and `Vessel_diam` — nXr ≈ (Vessel_diam / 1.5) × Leaflet_length

