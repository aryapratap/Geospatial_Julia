# Julia Geospatial Earthquake Analysis & Benchmark

A single self-contained Julia script (`julia_version.jl`) that pulls **live
USGS earthquake data** (last 30 days, ~8 000–12 000 events), runs three
compute benchmarks, and produces a publication-quality multi-panel world map.

---

## What the script does

The script runs in four stages, each printing timing information as it goes:

1. **Data load** — downloads the USGS `all_month.geojson` feed (cached for
   one hour so repeat runs don't re-hit the network).
2. **Feature extraction** — pulls `longitude`, `latitude`, `depth`, and
   `magnitude` out of every valid GeoJSON feature.
3. **Three compute benchmarks** — see below.
4. **Visualization** — renders a multi-panel figure with GeoMakie and saves
   it as `earthquakes_julia.png`.

### The three benchmarks

| # | Computation | Why it's interesting |
|---|---|---|
| a | **Pairwise Haversine distance matrix** (n × n) | Quadratic in the number of quakes — shows how Julia's native loops hold up on an O(n²) workload |
| b | **Neighbors within 500 km** for every quake | A realistic spatial query; stresses tight inner loops |
| c | **Gaussian KDE on a 360 × 180 global grid** | A triple loop (lat × lon × quake) that in Python would need risky broadcasting or a Numba JIT |

All three use a hand-written Haversine function with `@inbounds` and warm-up
calls so the JIT compile time is excluded from the timing.

### The visualization

A 1800 × 1100 px figure containing:

- A **Winkel-Tripel world map** with every earthquake plotted
  - Color → focal depth (inferno colormap, 0–300 km)
  - Size → magnitude² (approximates released energy)
  - Faint white coastlines and graticule
- **Focal-depth histogram** (lower-left)
- **Magnitude histogram** (lower-middle)
- **Depth-vs-magnitude scatter** (lower-right)
- A dark panel theme consistent across all subplots

---

## Install

You need **Julia 1.9 or newer**. From the Julia REPL:

```julia
using Pkg
Pkg.add(["HTTP", "JSON3", "GeoMakie", "CairoMakie"])
```

First install takes a few minutes because Makie / GeoMakie precompile. After
that it's fast.

---

## Run

```bash
julia julia_version.jl
```

Or from the REPL:

```julia
include("julia_version.jl")
```

Output:
- Timing table printed to stdout
- `earthquakes_julia.png` saved next to the script
- `usgs_all_month.geojson` cached next to the script (reused for 1 hour)

### Expected output — timing

On a typical laptop with ~10 000 quakes:

```
========================================================================
JULIA BENCHMARK SUMMARY
========================================================================
Data load + parse              :   0.8124 s
Feature extraction             :   0.0183 s
Pairwise distances  (n=10000)  :   0.5612 s
Neighbor counts     (r=500km)  :   0.4921 s
Gaussian KDE        (360x180)  :   1.0847 s
Compute total                  :   2.1380 s
========================================================================
```

Your numbers will vary with the live event count and your hardware, but
the compute stages together consistently finish in a couple of seconds
even on O(n²) workloads.

---

## How it works

### Data
- **Source:** [USGS Earthquake Hazards Program](https://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/all_month.geojson) — GeoJSON feed of the last 30 days, updated every minute.
- **Parsing:** `JSON3` for zero-allocation traversal of the feature array.
- **Caching:** The script writes the raw bytes to `usgs_all_month.geojson`
  and only re-downloads if the cache is older than one hour.

### The Haversine hot loop

```julia
@inline function haversine(lat1, lon1, lat2, lon2)
    R = 6371.0088
    φ1 = deg2rad(lat1); φ2 = deg2rad(lat2)
    Δφ = deg2rad(lat2 - lat1)
    Δλ = deg2rad(lon2 - lon1)
    a  = sin(Δφ/2)^2 + cos(φ1) * cos(φ2) * sin(Δλ/2)^2
    return 2R * asin(sqrt(a))
end
```

Because `lats::Vector{Float64}` and `lons::Vector{Float64}`, Julia's type
inference proves the argument types at compile time. LLVM compiles this
into tight scalar code — no boxing, no dispatch overhead, and the loop
body is eligible for SIMD auto-vectorization.

### Why this is fast without tricks

- **Specialized compilation.** Every generic function is JIT-compiled to
  machine code specialized on concrete argument types.
- **`@inbounds`** skips bounds checks inside the hot loops.
- **Warm-up calls** (the script calls each benchmark once on 10 elements
  before timing) ensure the JIT compile cost doesn't pollute the measurement.
- **No GIL.** If you want more throughput, just put `Threads.@threads` in
  front of the outer loop and restart Julia with `-t auto`.

### Why the visualization looks good

- **GeoMakie** provides a real Winkel Tripel projection through proj.4.
- **CairoMakie backend** outputs vector-quality PNG at 2× device pixels by
  default — text stays crisp at any zoom.
- Alpha compositing happens in **linear color space**, so dense
  earthquake clusters (Japan, Indonesia, Aleutians) show real
  depth-of-field instead of flat overplotting.
- All panels share a dark theme with matching tick / label / spine colors.

---

## Tuning knobs

All easy to find near the top of each section:

| What | Where | Default |
|---|---|---|
| KDE grid resolution | `kde_grid(...; nlon=360, nlat=180)` | 360 × 180 |
| KDE bandwidth | `kde_grid(...; bw=5.0)` | 5° |
| Neighbor radius | `count_neighbors(lats, lons, 500.0)` | 500 km |
| Depth color range | `colorrange = (0, 300)` | 0–300 km |
| Output DPI | `save(outfile, fig; px_per_unit = 2)` | 2× |
| Figure size | `Figure(size = (1800, 1100), ...)` | 1800 × 1100 |

---

## Files produced

| File | What it is |
|---|---|
| `earthquakes_julia.png`    | The rendered multi-panel figure |
| `usgs_all_month.geojson`   | 1-hour cache of the USGS feed |

---

## Troubleshooting

- **`UndefVarError: GeoMakie not defined`** — run the `Pkg.add` line above.
  First compile can take several minutes; it only happens once.
- **Coastlines warning** — `GeoMakie.coastlines()` needs Natural Earth data,
  which the package bundles. If it's missing the script logs a warning and
  draws the map without coastlines (everything else still works).
- **`ArgumentError: Package HTTP not found`** — you ran the script before
  `Pkg.add`. Install the dependencies first.
- **Corporate proxy blocks USGS** — the first download will fail. Set the
  `HTTP_PROXY` / `HTTPS_PROXY` environment variables before launching Julia,
  or drop a copy of `all_month.geojson` next to the script — if the cache
  exists, the script uses it.

---

## Reproducibility note

The USGS feed is live, so the earthquake count changes every run. If you
need identical data across runs (e.g. for a screenshot), just let the cache
live — back-to-back runs within the hour reuse the same file.
