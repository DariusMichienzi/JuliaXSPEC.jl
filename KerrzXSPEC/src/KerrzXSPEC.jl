"""
    KerrzXSPEC

Relativistic line-profile kernels from the kerrz command-line tool for four
corona geometries (lamppost, ring, disc-like, umbrella), registered as
JuliaXSPEC models.

`kzlamp`, `kzring`, `kzdisc` and `kzumb` are convolutions: in XSPEC each blurs
whatever rest-frame spectrum it multiplies. `kzlampx`, `kzringx`, `kzdiscx` and
`kzumbx` blur xillverD-5 on the table's own energy grid.

Line profiles are calculated by kerrz on a grid of parameters when a fit first
needs them and remembered on disk; see `geometries.jl` for the grids and
`cli.jl` for how kerrz is run. kerrz must be on the PATH (or `KERRZ_BIN` set)
wherever XSPEC runs.
"""
module KerrzXSPEC

using Printf
using JuliaXSPEC
using JuliaXSPECModels

# ---------------------------------------------------------------- settings
# Changing any of these changes the cache names, so old profiles are not mixed
# with new ones.

"Outer edge of the disc in GM/c² (inner edge: the ISCO)."
const R_OUT = 1000.0

"Photon index of the coronal spectrum, used for the corona-to-disc energy shift in the emissivity."
const EMISSIVITY_PHOTON_INDEX = 2.0

"Extra arguments for every `kerrz emissivity` run, e.g. `[\"--velocity\", \"corotate\"]`."
const EMISSIVITY_EXTRA_ARGS = String[]

# Photons traced per emissivity profile; for the extended coronae, per ring of
# the corona. These are a balance of noise against time: 20 000 per ring for
# the disc-like corona agrees with 200 000 to about 1 % in the line profile.
const LAMP_NPHOTONS = 20_000
const RING_NPHOTONS = 100_000
const DISC_NPHOTONS = 20_000
const UMB_NPHOTONS = 20_000

"Rings composing the disc-like and umbrella coronae."
const DISC_NRINGS = 10
const UMB_NRINGS = 10

"""
The g grid of every profile: kerrz's own lineprof grid (1000 points over
0 ≤ g ≤ 2) without g = 0, which a RedshiftKernel cannot hold.
"""
const G_GRID = collect(range(0.0, 2.0; length = 1000))[2:end]

include("cli.jl")
include("geometries.jl")

function __init__()
    for model in (kzlamp, kzring, kzdisc, kzumb, kzlampx, kzringx, kzdiscx, kzumbx)
        register!(model)
    end
end

end # module
