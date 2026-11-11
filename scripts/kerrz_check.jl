# Sanity checks for the kerrz corona models, outside XSPEC.
#
#   julia --project=builds/kerrz scripts/kerrz_check.jl [lamp ring disc umbrella]
#
# For each geometry: the kernel at the default parameters (area, mean g and
# where it is non-zero), then the grid-interpolated profile at a point between
# grid nodes against the exact kerrz profile there. That difference is the
# interpolation error of the grid; if it is too large for your data, add grid
# points along the axis that matters (see KerrzXSPEC/src/geometries.jl).

using JuliaXSPEC, KerrzXSPEC
using Printf

const K = KerrzXSPEC
const PARAMS = (lamp = K.LAMP_PARAMETERS, ring = K.RING_PARAMETERS,
                disc = K.DISC_PARAMETERS, umbrella = K.UMB_PARAMETERS)

"A point between grid nodes: halfway (in the axis coordinate) above each default."
function off_grid_point(profiles, defaults)
    map(enumerate(profiles.axes)) do (i, axis)
        v = axis.values
        j = clamp(searchsortedfirst(v, defaults[i]), 1, length(v) - 1)
        axis.scale == :log ? sqrt(v[j] * v[j+1]) : (v[j] + v[j+1]) / 2
    end |> Tuple
end

describe(L) = (g = K.G_GRID; nz = findall(>(0), L);
               @sprintf("mean g %.4f, non-zero for %.3f ≤ g ≤ %.3f", sum(g .* L) / sum(L), g[first(nz)], g[last(nz)]))

selected = isempty(ARGS) ? collect(string.(keys(K.GEOMETRIES))) : ARGS
for name in selected
    sym = Symbol(name)
    profiles = K.GEOMETRIES[sym]
    defaults = Tuple(p.default for p in PARAMS[sym])
    println("== $name")
    t = @elapsed kernel = K.profile_kernel(profiles, defaults)
    area = sum(diff(kernel.g) .* (kernel.L[1:end-1] .+ kernel.L[2:end]) ./ 2)
    @printf("  defaults %s: area %.6f, %s (%.1f s)\n", defaults, area, describe(kernel.L), t)

    point = off_grid_point(profiles, defaults)
    t = @elapsed begin
        interpolated = profiles(point...)
        exact = evaluate_exact(profiles, point...)
    end
    a = interpolated ./ sum(interpolated)
    b = exact ./ sum(exact)
    @printf("  between nodes %s: max |interpolated - exact| = %.1f %% of the peak (%.1f s)\n",
            map(x -> round(x; sigdigits = 4), point), 100 * maximum(abs.(a .- b)) / maximum(b), t)
end
