# Vessel Wall Motion: Speed-Based → Position-Based

## Summary of Changes

The vessel wall peristaltic motion in `clMethod.h` has been refactored from **speed-based** (velocity-driven) to **position-based** (location-driven). This decouples the motion from time-stepping and provides more stable, direct control over node positions.

---

## Old Approach (Speed-Based)

### Logic
```cpp
double speed = OSC_AMPLITUDE * sin(2π × localX/INTERVAL) * sin(2π × Time/RAMP_PERIOD);
vls.V(i).y = speed * (y - midy) / radius;
vls.V(i).z = speed * (z - midz) / radius;
```

### Issues
- Velocity is calculated and applied each timestep
- Motion amplitude depends on the speed formula **and** how the solver integrates velocity
- Decouples the intended motion from actual position reached
- Less predictable in terms of actual radial displacement

---

## New Approach (Position-Based)

### Logic
```cpp
// 1. Calculate what the radius SHOULD be at this location and time
double radiusModulation = OSC_AMPLITUDE * sin(2π × localX/INTERVAL) * sin(2π × Time/RAMP_PERIOD);
double targetRadius = baseRadius + radiusModulation;

// 2. Scale current position to match target radius
double currentRadius = sqrt((y - midy)² + (z - midz)²);
double scaleFactor = targetRadius / currentRadius;

double targetY = midy + (y - midy) × scaleFactor;
double targetZ = midz + (z - midz) × scaleFactor;

// 3. Calculate velocity needed to reach target in one timestep
vls.V(i).y = (targetY - C(i).y) / dTls;
vls.V(i).z = (targetZ - C(i).z) / dTls;
```

### Benefits
- **Position is explicit**: We directly specify where the node should be
- **Decoupled from speed**: The formula controls position, not velocity
- **Deterministic**: Same input (location, time) → same target position
- **Stable**: No accumulation of velocity errors
- **Adaptive**: Velocity automatically adjusts to reach target in current timestep

---

## Key Variables

| Variable | Meaning |
|----------|---------|
| `radiusModulation` | How much to add/subtract from base radius (replaces `speed`) |
| `baseRadius` | Vessel diameter / 2 (equilibrium radius) |
| `targetRadius` | `baseRadius + radiusModulation` (desired radius at this instant) |
| `currentRadius` | Actual radial distance of node from vessel centerline |
| `scaleFactor` | `targetRadius / currentRadius` — how to scale node's current position |
| `targetY, targetZ` | Where the node **should be** at this instant |

---

## Motion Formula (Unchanged)

The **wave formula itself is identical**:

```
radiusModulation = OSC_AMPLITUDE × sin(2π × localX / INTERVAL) × sin(2π × Time / RAMP_PERIOD)
```

- `OSC_AMPLITUDE` — peak modulation (e.g., ±3 units)
- `localX` — distance from valve region start
- `INTERVAL` — spatial wavelength
- `Time` — current simulation time
- `RAMP_PERIOD` — temporal period (one full cycle)

The only difference: this formula now directly controls **radius**, not velocity.

---

## Static Regions (Unchanged)

Entrance and exit regions remain stationary:
- Outside valve region: `radiusModulation = 0` → `targetRadius = baseRadius` → no motion
- This is enforced through the location-based logic (same as before)

---

## Implementation Details

### Functions Modified
1. `solveVerlet1()` — lines 395–443 (approx.)
2. `solveVerlet2()` — lines 639–694 (approx.)

### Edge Cases
- **Near-zero current radius** (node at vessel center): Skip motion calculation to avoid division by zero
- **Leaflet nodes** (O > 1): Same logic applies — they experience the same radial modulation if they fall within valve region
- **Outside valve buffer** (10-unit margins): Nodes treated as stationary

---

## Simulation Impact

### What Changes
- Vessel wall contracts/expands **to a specific radius** (not just at a certain speed)
- Motion is **reproducible** — same location + time = same position
- Energy dissipation is absorbed by the FSI coupling (fluid pressure), not hidden in velocity integration

### What Stays the Same
- The peristaltic wave **frequency** and **amplitude** remain as configured
- **Entrance/exit regions** still stationary
- **Valve region** still contains the active wave
- **Overall pumping effect** — fluid moves forward due to the wave

---

## Testing / Validation

To verify the change works correctly:

1. Run a short simulation (e.g., 1–2 cardiac cycles)
2. Monitor vessel wall node positions and velocities:
   - In valve region: should see smooth radial oscillation
   - Outside valve region: should remain at constant radius
   - Velocity should be smooth (no spikes at region boundaries)
3. Check fluid flow: pump should still move fluid forward through valves
4. Compare to earlier speed-based results: final flow rate should be similar or improved

---

## Future Refinements

Potential enhancements (not in this update):
- Smooth transition at valve region boundaries (add buffer zone with ramping)
- Separate amplitude control for different valve pairs
- Non-uniform wavelength (e.g., shorter wavelength near valves)
- Damping factor to control how fast nodes move to target position

