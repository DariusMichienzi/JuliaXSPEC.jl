# The four kerrz corona geometries. For each one:
#
#   *_PARAMETERS  the XSPEC parameters: spin, inclination, then the corona's own
#   *_PROFILES    a GridInterpolator of line profiles over the same parameters,
#                 filled in on demand and remembered on disk
#
# and two XSPEC models built from them by `kerrz_models`:
#
#   kz<name>      a convolution model: blurs whatever it multiplies in XSPEC,
#                 e.g. `kzring * atable{my_restframe.fits}`
#   kz<name>x     xillverD-5 blurred on the table's own energy grid
#
# Grid points are deliberately explicit vectors so they are easy to edit. Each
# parameter's limits are the ends of its grid, and every default is a grid
# point, so the first evaluation at the defaults needs a single corner.
# Changing a grid invalidates only that geometry's cached profiles; the kerrz
# table and emissivity files keyed by exact values are reused where they match.

# ---------------------------------------------------------------- shared axes

# Spin exactly 0 makes kerrz 0.4.1 crash (table) or hang (emissivity), so the
# grid starts at 0.1, the same lower limit as the Isambard table model.
const SPIN_NODES = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.98, 0.998]
const INC_NODES = vcat(collect(5.0:1.0:30.0), collect(30.5:0.5:65.0), collect(65.2:0.2:85.0))

spin_parameter() = Parameter("a", 0.9; min = first(SPIN_NODES), max = last(SPIN_NODES), delta = 0.01)
inc_parameter() = Parameter("inc", 30.0; unit = "degrees", min = first(INC_NODES), max = last(INC_NODES), delta = 1.0)

"Cache name for one geometry: changes whenever a setting that alters its profiles changes."
function cache_name(geometry::AbstractString, nphotons::Integer)
    extra = isempty(EMISSIVITY_EXTRA_ARGS) ? "" : "-x" * string(hash(EMISSIVITY_EXTRA_ARGS); base = 16)
    return "kerrz-$(geometry)-v1-n$(nphotons)-G$(fmt(EMISSIVITY_PHOTON_INDEX))-rout$(fmt(R_OUT))$(extra)"
end

"The interpolated line profile at `values` as a unit-area kernel in g."
profile_kernel(profiles, values::Tuple) = RedshiftKernel(G_GRID, max.(profiles(values...), 0.0))

"""
    kerrz_models(name, parameters, profiles, what) -> (convolution, xillver_blurred)

The two XSPEC models for one geometry: `name` (a convolution) and `name * "x"`
(xillverD-5 blurred on its own grid).
"""
function kerrz_models(name::AbstractString, parameters::Vector{Parameter}, profiles, what::AbstractString)
    conv = ConvolutionModel(
        name, parameters;
        description = "kerrz $(what) blur of the components it multiplies",
    ) do spectrum, p
        convolve(spectrum, profile_kernel(profiles, values(p)); method = :fft)
    end
    blurred = Blurred(
        name * "x", JuliaXSPECModels.xillver_table, parameters, JuliaXSPECModels.XILLVER_PARAMETERS;
        description = "xillverD-5 blurred by the kerrz $(what), on the table energy grid",
        method = :fft,
    ) do kernel
        profile_kernel(profiles, values(kernel))
    end
    return conv, blurred
end

# ---------------------------------------------------------------- lamppost

const LAMP_H = [2.5, 3.0, 4.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0, 50.0, 75.0, 100.0]

const LAMP_PARAMETERS = [
    spin_parameter(),
    inc_parameter(),
    Parameter("h", 10.0; unit = "GM/c^2", min = first(LAMP_H), max = last(LAMP_H), delta = 0.1),
]

const LAMP_PROFILES = GridInterpolator(
    (a = SPIN_NODES, inc = INC_NODES, h = GridAxis(LAMP_H; scale = :log));
    cache = cache_name("lamp", LAMP_NPHOTONS),
) do a, inc, h
    em = emissivity_file(a, "--lamppost", "h:$(fmt(h)),vr:0", LAMP_NPHOTONS)
    line_profile(a, inc, em)
end

# ---------------------------------------------------------------- ring

const RING_H = [2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0, 50.0]
const RING_X = [1.5, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0]

const RING_PARAMETERS = [
    spin_parameter(),
    inc_parameter(),
    Parameter("h", 5.0; unit = "GM/c^2", min = first(RING_H), max = last(RING_H), delta = 0.1),
    Parameter("x", 4.0; unit = "GM/c^2", min = first(RING_X), max = last(RING_X), delta = 0.1),
]

const RING_PROFILES = GridInterpolator(
    (a = SPIN_NODES, inc = INC_NODES, h = GridAxis(RING_H; scale = :log), x = GridAxis(RING_X; scale = :log));
    cache = cache_name("ring", RING_NPHOTONS),
) do a, inc, h, x
    em = emissivity_file(a, "--ring-like", "h:$(fmt(h)),x:$(fmt(x))", RING_NPHOTONS)
    line_profile(a, inc, em)
end

# ---------------------------------------------------------------- disc-like
# Parameterised by inner radius and radial extent rather than inner and outer
# radius, so that every grid corner is a valid corona (rout > rin).

const DISC_H = [2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0]
const DISC_RIN = [1.5, 2.0, 3.0, 5.0, 8.0]
const DISC_DR = [1.0, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 14.0, 20.0]

const DISC_PARAMETERS = [
    spin_parameter(),
    inc_parameter(),
    Parameter("h", 5.0; unit = "GM/c^2", min = first(DISC_H), max = last(DISC_H), delta = 0.1),
    Parameter("rin", 2.0; unit = "GM/c^2", min = first(DISC_RIN), max = last(DISC_RIN), delta = 0.1),
    Parameter("dr", 5.0; unit = "GM/c^2", min = first(DISC_DR), max = last(DISC_DR), delta = 0.1),
]

const DISC_PROFILES = GridInterpolator(
    (a = SPIN_NODES, inc = INC_NODES, h = GridAxis(DISC_H; scale = :log),
     rin = GridAxis(DISC_RIN; scale = :log), dr = GridAxis(DISC_DR; scale = :log));
    cache = cache_name("disc", DISC_NPHOTONS),
) do a, inc, h, rin, dr
    props = "h:$(fmt(h)),rin:$(fmt(rin)),rout:$(fmt(rin + dr)),nr:$(DISC_NRINGS)"
    em = emissivity_file(a, "--disc-like", props, DISC_NPHOTONS)
    line_profile(a, inc, em)
end

# ---------------------------------------------------------------- umbrella
# Parameterised by the inner opening angle and the angular width, so that the
# outer angle (inner + width) is always the larger.

const UMB_R = [2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 15.0, 20.0]
const UMB_INNER = [0.0, 10.0, 20.0, 30.0, 45.0]
const UMB_WIDTH = [5.0, 10.0, 20.0, 30.0, 45.0]

const UMB_PARAMETERS = [
    spin_parameter(),
    inc_parameter(),
    Parameter("r", 5.0; unit = "GM/c^2", min = first(UMB_R), max = last(UMB_R), delta = 0.1),
    Parameter("inner", 10.0; unit = "degrees", min = first(UMB_INNER), max = last(UMB_INNER), delta = 1.0),
    Parameter("width", 30.0; unit = "degrees", min = first(UMB_WIDTH), max = last(UMB_WIDTH), delta = 1.0),
]

const UMB_PROFILES = GridInterpolator(
    (a = SPIN_NODES, inc = INC_NODES, r = GridAxis(UMB_R; scale = :log),
     inner = UMB_INNER, width = UMB_WIDTH);
    cache = cache_name("umbrella", UMB_NPHOTONS),
) do a, inc, r, inner, width
    props = "r:$(fmt(r)),inner:$(fmt(inner)),angle:$(fmt(inner + width)),nr:$(UMB_NRINGS)"
    em = emissivity_file(a, "--umbrella", props, UMB_NPHOTONS)
    line_profile(a, inc, em)
end

# ---------------------------------------------------------------- models

const kzlamp, kzlampx = kerrz_models("kzlamp", LAMP_PARAMETERS, LAMP_PROFILES, "lamppost")
const kzring, kzringx = kerrz_models("kzring", RING_PARAMETERS, RING_PROFILES, "ring corona")
const kzdisc, kzdiscx = kerrz_models("kzdisc", DISC_PARAMETERS, DISC_PROFILES, "disc-like corona")
const kzumb, kzumbx = kerrz_models("kzumb", UMB_PARAMETERS, UMB_PROFILES, "umbrella corona")

"Every geometry's profile cache, by name; used by the check and prewarm scripts."
const GEOMETRIES = (
    lamp = LAMP_PROFILES,
    ring = RING_PROFILES,
    disc = DISC_PROFILES,
    umbrella = UMB_PROFILES,
)
