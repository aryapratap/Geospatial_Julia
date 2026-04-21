#=
    julia_version.jl
    ----------------
    Geospatial analysis + visualization benchmark in Julia.

    What it does:
      1. Downloads the last 30 days of earthquakes from USGS (GeoJSON feed)
      2. Extracts lat / lon / depth / magnitude
      3. Runs three compute benchmarks:
           (a) full pairwise Haversine distance matrix
           (b) neighbor-count within 500 km for every quake
           (c) 2-D Gaussian kernel density estimate on a 360x180 grid
      4. Renders a high-quality multi-panel figure with GeoMakie
         (world map color-coded by focal depth, size by magnitude,
         plus depth / magnitude / scatter sub-panels).

    First run:   using Pkg; Pkg.add(["HTTP","JSON3","GeoMakie","CairoMakie"])
    Then:        julia --project julia_version.jl
=#

using HTTP, JSON3
using GeoMakie, CairoMakie
using Printf, Statistics

const BAR = "=" ^ 72
println(BAR)
println("JULIA: Geospatial Earthquake Analysis & Benchmark")
println(BAR)

# ---------------------------------------------------------------- data load
println("\n[1/4] Loading USGS earthquake data (last 30 days)...")
const URL = "https://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/all_month.geojson"

# tiny cache to avoid repeated network hits during iteration
const CACHE = joinpath(@__DIR__, "usgs_all_month.geojson")

t0 = time()
raw = if isfile(CACHE) && (time() - mtime(CACHE) < 3600)
    read(CACHE)
else
    body = HTTP.get(URL).body
    open(CACHE, "w") do io; write(io, body); end
    body
end
data = JSON3.read(raw)
t_load = time() - t0
@printf("   Download + parse              : %8.4f s\n", t_load)

# ---------------------------------------------------------------- extraction
println("\n[2/4] Extracting earthquake features...")
t0 = time()
features = data.features

lons   = Float64[]
lats   = Float64[]
depths = Float64[]
mags   = Float64[]

for f in features
    c = f.geometry.coordinates
    length(c) >= 3 || continue
    mag = f.properties.mag
    mag === nothing && continue
    push!(lons,   Float64(c[1]))
    push!(lats,   Float64(c[2]))
    push!(depths, Float64(c[3]))
    push!(mags,   Float64(mag))
end
n = length(lats)
t_extract = time() - t0
@printf("   Extracted %-5d quakes         : %8.4f s\n", n, t_extract)

# ---------------------------------------------------------------- benchmarks
println("\n[3/4] Running compute benchmarks...")

@inline function haversine(lat1, lon1, lat2, lon2)
    R = 6371.0088
    φ1 = deg2rad(lat1); φ2 = deg2rad(lat2)
    Δφ = deg2rad(lat2 - lat1)
    Δλ = deg2rad(lon2 - lon1)
    a  = sin(Δφ/2)^2 + cos(φ1) * cos(φ2) * sin(Δλ/2)^2
    return 2R * asin(sqrt(a))
end

# --- (a) pairwise distance matrix --------------------------------
function pairwise_distances(lats, lons)
    n = length(lats)
    D = zeros(Float64, n, n)
    @inbounds for i in 1:n-1
        lati = lats[i]; loni = lons[i]
        for j in i+1:n
            d = haversine(lati, loni, lats[j], lons[j])
            D[i, j] = d
            D[j, i] = d
        end
    end
    return D
end

pairwise_distances(lats[1:10], lons[1:10])             # JIT warmup
t0 = time()
D = pairwise_distances(lats, lons)
t_dist = time() - t0
@printf("   Pairwise distances (%d x %d)  : %8.4f s\n", n, n, t_dist)

# --- (b) neighbor counting --------------------------------------
function count_neighbors(lats, lons, radius_km)
    n = length(lats)
    counts = zeros(Int, n)
    @inbounds for i in 1:n
        c = 0
        lati = lats[i]; loni = lons[i]
        for j in 1:n
            i == j && continue
            if haversine(lati, loni, lats[j], lons[j]) <= radius_km
                c += 1
            end
        end
        counts[i] = c
    end
    return counts
end

count_neighbors(lats[1:10], lons[1:10], 500.0)          # warmup
t0 = time()
neighbor_counts = count_neighbors(lats, lons, 500.0)
t_neigh = time() - t0
@printf("   Neighbors within 500 km        : %8.4f s\n", t_neigh)

# --- (c) Gaussian KDE on a grid ---------------------------------
function kde_grid(lats, lons, weights; nlon=360, nlat=180, bw=5.0)
    lon_grid = range(-180, 180, length=nlon)
    lat_grid = range(-90,  90,  length=nlat)
    Z = zeros(Float64, nlat, nlon)
    inv2b2 = 1.0 / (2 * bw^2)
    @inbounds for j in 1:nlon
        glon = lon_grid[j]
        for i in 1:nlat
            glat = lat_grid[i]
            s = 0.0
            for k in eachindex(lats)
                dlat = glat - lats[k]
                dlon = glon - lons[k]
                s += weights[k] * exp(-(dlat*dlat + dlon*dlon) * inv2b2)
            end
            Z[i, j] = s
        end
    end
    return lon_grid, lat_grid, Z
end

kde_grid(lats[1:10], lons[1:10], mags[1:10]; nlon=36, nlat=18)   # warmup
t0 = time()
lon_grid, lat_grid, density = kde_grid(lats, lons, mags)
t_kde = time() - t0
@printf("   Gaussian KDE (360 x 180 grid)  : %8.4f s\n", t_kde)

# ---------------------------------------------------------------- visualize
println("\n[4/4] Rendering visualization with GeoMakie...")

BG     = RGBf(0.045, 0.05, 0.09)
PANEL  = RGBf(0.07,  0.085, 0.14)
AXCOL  = RGBf(0.75,  0.78,  0.85)
TEXT   = :white

fig = Figure(
    size = (1800, 1100),
    backgroundcolor = BG,
    fontsize = 14,
)

# ---- title band --------------------------------------------------
Label(fig[0, 1:4],
      "Global Seismicity — Last 30 Days";
      fontsize = 30, color = TEXT, font = :bold,
      halign = :left, padding = (20, 0, 0, 0))
Label(fig[0, 1:4, Bottom()],
      @sprintf("USGS feed · %d earthquakes · rendered in Julia / GeoMakie", n);
      fontsize = 14, color = AXCOL,
      halign = :left, padding = (20, 0, 2, 0))

# ---- main world map ----------------------------------------------
ga = GeoAxis(fig[1, 1:3];
    dest = "+proj=wintri",
    xgridcolor = RGBAf(1, 1, 1, 0.07),
    ygridcolor = RGBAf(1, 1, 1, 0.07),
    xticklabelcolor = AXCOL,
    yticklabelcolor = AXCOL,
)

# coastlines as a thin pale stroke
try
    lines!(ga, GeoMakie.coastlines();
           color = RGBAf(1, 1, 1, 0.35), linewidth = 0.5)
catch err
    @warn "Coastlines not available offline, skipping." exception=err
end

# marker size scales with magnitude squared (energy ~ 10^(1.5 M))
msize = @. 3.5 + (mags ^ 2) * 0.8

sc = scatter!(ga, lons, lats;
    color       = depths,
    colormap    = :inferno,
    colorrange  = (0, 300),
    markersize  = msize,
    strokewidth = 0.15,
    strokecolor = RGBAf(1, 1, 1, 0.55),
    alpha       = 0.9,
)

Colorbar(fig[1, 4], sc;
    label = "Focal depth (km)",
    labelcolor = TEXT, labelsize = 14,
    ticklabelcolor = AXCOL,
    height = Relative(0.72),
    width  = 18,
)

# ---- lower panels ------------------------------------------------
function styled_axis(pos; title, xlabel, ylabel="")
    ax = Axis(fig[pos...];
        title = title, titlecolor = TEXT, titlesize = 16,
        xlabel = xlabel, ylabel = ylabel,
        xlabelcolor = AXCOL, ylabelcolor = AXCOL,
        xticklabelcolor = AXCOL, yticklabelcolor = AXCOL,
        xtickcolor = AXCOL, ytickcolor = AXCOL,
        backgroundcolor = PANEL,
        xgridcolor = RGBAf(1, 1, 1, 0.06),
        ygridcolor = RGBAf(1, 1, 1, 0.06),
        topspinecolor = RGBAf(1,1,1,0.15),
        bottomspinecolor = RGBAf(1,1,1,0.15),
        leftspinecolor = RGBAf(1,1,1,0.15),
        rightspinecolor = RGBAf(1,1,1,0.15),
    )
    return ax
end

ax_d = styled_axis((2, 1); title="Focal depth distribution",
                   xlabel="Depth (km)", ylabel="Count")
hist!(ax_d, depths; bins = 60,
      color = :orangered, strokecolor = RGBAf(1,1,1,0.5), strokewidth = 0.3)

ax_m = styled_axis((2, 2); title="Magnitude distribution",
                   xlabel="Magnitude", ylabel="Count")
hist!(ax_m, mags; bins = 40, color = :gold,
      strokecolor = RGBAf(1,1,1,0.5), strokewidth = 0.3)

ax_s = styled_axis((2, 3); title="Depth vs Magnitude",
                   xlabel="Magnitude", ylabel="Depth (km)")
scatter!(ax_s, mags, depths; color = depths, colormap = :inferno,
         markersize = 5, alpha = 0.75)

# column widths: give map the bulk of the space
colsize!(fig.layout, 4, Fixed(60))
rowsize!(fig.layout, 1, Relative(0.70))
rowsize!(fig.layout, 2, Relative(0.25))

outfile = joinpath(@__DIR__, "earthquakes_julia.png")
save(outfile, fig; px_per_unit = 2)
println("   Saved: ", outfile)

# ---------------------------------------------------------------- summary
println("\n" * BAR)
println("JULIA BENCHMARK SUMMARY")
println(BAR)
@printf("Data load + parse              : %8.4f s\n", t_load)
@printf("Feature extraction             : %8.4f s\n", t_extract)
@printf("Pairwise distances  (n=%d)  : %8.4f s\n", n, t_dist)
@printf("Neighbor counts     (r=500km)  : %8.4f s\n", t_neigh)
@printf("Gaussian KDE        (360x180)  : %8.4f s\n", t_kde)
@printf("Compute total                  : %8.4f s\n",
        t_dist + t_neigh + t_kde)
println(BAR)
