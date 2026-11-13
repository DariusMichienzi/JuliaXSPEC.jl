# How accurately the kerrz line-profile grids interpolate.
#
#   julia --project=builds/kerrz scripts/kerrz_interp_test.jl [geometries...] [options]
#
#   geometries: any of lamp, ring, disc, umbrella (default: all four)
#
# For every pair of neighbouring grid points along an axis, with all other
# parameters at their defaults, this compares the profile the grid gives
# halfway between them with the exact kerrz profile there. Halfway is in the
# axis's interpolation coordinate (the geometric mean on a log axis), which is
# where linear interpolation is usually worst.
#
# The error is max |interpolated - exact| as a percentage of the exact
# profile's peak, after smoothing both with a Gaussian of width sigma in g.
# sigma = 0 is kerrz's own resolution (g spacing 0.002). As a guide,
# sigma_g = sigma_E / E: a CCD (about 60 eV sigma at 6.4 keV) is about 0.01.
# The emissivity profiles carry about 1 % Monte Carlo noise, so errors of a
# percent or two are at the noise floor.
#
# Options:
#   --axes=inc,a            only these axes (default: all; names as in XSPEC)
#   --sigma=0,0.005,0.02    smoothing widths in g, one column each
#   --stride=n              test every n-th interval (default 1: all)
#   --threshold=3           flag intervals worse than this many percent,
#                           judged at the first non-zero sigma
#   --out=results.csv       also write every result to a CSV file
#
# Each interval needs one exact profile at a value off the grid, i.e. new kerrz
# runs: about a second for a lamppost or a new inclination, a few seconds for
# the disc-like and umbrella emissivities. A fine inclination axis therefore
# takes a few minutes per geometry; --stride or --axes shortens a run. Results
# are cached as usual, so a repeat run is quick.

using JuliaXSPEC, KerrzXSPEC
using Printf

const K = KerrzXSPEC
const PARAMS = (lamp = K.LAMP_PARAMETERS, ring = K.RING_PARAMETERS,
                disc = K.DISC_PARAMETERS, umbrella = K.UMB_PARAMETERS)

function option(name::AbstractString, default::AbstractString)
    for arg in ARGS
        startswith(arg, "--$(name)=") && return String(split(arg, '='; limit = 2)[2])
    end
    return default
end

const SELECTED = let chosen = [a for a in ARGS if !startswith(a, "--")]
    isempty(chosen) ? [string(k) for k in keys(K.GEOMETRIES)] : chosen
end
const AXIS_FILTER = let s = option("axes", "")
    isempty(s) ? nothing : split(s, ',')
end
const SIGMAS = parse.(Float64, split(option("sigma", "0,0.005,0.02"), ','))
const INTERVAL_STRIDE = parse(Int, option("stride", "1"))
const THRESHOLD = parse(Float64, option("threshold", "3"))
const OUTFILE = option("out", "")
const JUDGE = something(findfirst(>(0), SIGMAS), 1)   # which column the threshold applies to

"Trapezoidal area on G_GRID."
area(L) = sum(diff(K.G_GRID) .* (L[1:end-1] .+ L[2:end]) ./ 2)
normalise(L) = (L = max.(L, 0.0); L ./ area(L))

"Gaussian smoothing in g, of width sigma; sigma = 0 returns L unchanged."
function smooth(L::Vector{Float64}, sigma::Float64)
    sigma <= 0 && return L
    dg = K.G_GRID[2] - K.G_GRID[1]
    half = ceil(Int, 5 * sigma / dg)
    offsets = -half:half
    w = [exp(-0.5 * (k * dg / sigma)^2) for k in offsets]
    w ./= sum(w)
    n = length(L)
    out = zeros(n)
    for i in 1:n, (j, k) in enumerate(offsets)
        m = i + k
        1 <= m <= n && (out[i] += w[j] * L[m])
    end
    return out
end

"Max |a - b| as a percentage of the peak of b, after smoothing both."
function percent_error(a, b, sigma)
    sa = smooth(a, sigma)
    sb = smooth(b, sigma)
    return 100 * maximum(abs.(sa .- sb)) / maximum(sb)
end

midpoint(axis, lo, hi) = axis.scale == :log ? sqrt(lo * hi) : (lo + hi) / 2

sigma_header() = join([@sprintf("σ=%-6g", s) for s in SIGMAS], " ")

results = NamedTuple[]
worst_per_axis = NamedTuple[]

for gname in SELECTED
    sym = Symbol(gname)
    haskey(K.GEOMETRIES, sym) || error("unknown geometry $(gname): use lamp, ring, disc or umbrella")
    profiles = K.GEOMETRIES[sym]
    defaults = [p.default for p in PARAMS[sym]]
    axis_names = [string(n) for n in profiles.names]
    println("\n== $(gname)  (defaults: ", join(["$(n)=$(d)" for (n, d) in zip(axis_names, defaults)], ", "), ")")
    @printf("  %-6s %-19s  %s   time\n", "axis", "interval", sigma_header())

    for (k, axis) in enumerate(profiles.axes)
        axname = axis_names[k]
        AXIS_FILTER === nothing || axname in AXIS_FILTER || continue
        v = axis.values
        length(v) >= 2 || continue
        worst = (err = -Inf, lo = NaN, hi = NaN)
        for j in 1:INTERVAL_STRIDE:length(v)-1
            point = copy(defaults)
            point[k] = midpoint(axis, v[j], v[j+1])
            local interp, exact
            t = @elapsed begin
                interp = normalise(profiles(point...))
                exact = normalise(evaluate_exact(profiles, point...))
            end
            errs = [percent_error(interp, exact, s) for s in SIGMAS]
            flag = errs[JUDGE] > THRESHOLD ? "  <-- above $(THRESHOLD) %" : ""
            @printf("  %-6s %8.4g – %-8.4g  %s  %5.1f s%s\n", axname, v[j], v[j+1],
                    join([@sprintf("%7.2f ", e) for e in errs], " "), t, flag)
            push!(results, (geometry = gname, axis = axname, lo = v[j], hi = v[j+1],
                            mid = point[k], errors = errs))
            errs[JUDGE] > worst.err && (worst = (err = errs[JUDGE], lo = v[j], hi = v[j+1]))
        end
        push!(worst_per_axis, merge((geometry = gname, axis = axname), worst))
    end
end

println("\n== Summary: worst interval per axis at σ = $(SIGMAS[JUDGE]) (threshold $(THRESHOLD) %)")
for s in worst_per_axis
    verdict = s.err > THRESHOLD ? "add grid points" : "ok"
    @printf("  %-9s %-6s %6.2f %%  between %g and %g   %s\n", s.geometry, s.axis, s.err, s.lo, s.hi, verdict)
end

if !isempty(OUTFILE)
    open(OUTFILE, "w") do io
        println(io, "geometry,axis,lo,hi,mid,", join(["err_sigma_$(s)" for s in SIGMAS], ","))
        for r in results
            println(io, join([r.geometry, r.axis, r.lo, r.hi, r.mid, r.errors...], ","))
        end
    end
    println("\nWrote every interval to $(OUTFILE)")
end
