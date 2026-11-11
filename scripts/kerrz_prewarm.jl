# Calculate every grid corner of one or more kerrz geometries ahead of time, so
# that fits never wait for kerrz. Safe to stop and restart: finished corners
# are on disk and are skipped.
#
#   julia --project=builds/kerrz scripts/kerrz_prewarm.jl lamp ring disc umbrella
#
# Set KERRZ_NTHREADS to control kerrz's threads. Rough single-thread totals
# with the default grids: lamp a few minutes, ring about half an hour, disc and
# umbrella a couple of hours each; several threads cut the emissivity part.

using JuliaXSPEC, KerrzXSPEC
using Printf

isempty(ARGS) && error("name the geometries to fill: lamp, ring, disc and/or umbrella")
for name in ARGS
    profiles = KerrzXSPEC.GEOMETRIES[Symbol(name)]
    nodes = [axis.values for axis in profiles.axes]
    total = prod(length, nodes)
    println("== $name: $total corners")
    start = time()
    for (k, point) in enumerate(Iterators.product(nodes...))
        profiles(point...)   # at a grid node only that corner is used
        if k % 100 == 0 || k == total
            elapsed = time() - start
            @printf("  %d / %d corners, %.0f s elapsed, about %.0f s left\n",
                    k, total, elapsed, elapsed / k * (total - k))
        end
        # Keep memory flat on long runs; the corners stay on disk.
        k % 2000 == 0 && empty_cache!(profiles)
    end
end
