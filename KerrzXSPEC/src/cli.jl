# Running the kerrz command-line tool and keeping what it writes.
#
# A line profile is built in three kerrz steps, and the two expensive ones are
# remembered as files under `kerrz_files_dir()`:
#
#   kerrz table       Cunningham transfer functions, one file per (spin, inc)
#   kerrz emissivity  the corona's emissivity profile, one file per (spin, corona)
#   kerrz lineprof    integrates the two (about 0.05 s); not stored here, because
#                     each geometry's GridInterpolator keeps the profiles.
#
# So a new grid corner usually costs one emissivity run at most: the tables are
# shared by every corona, and an emissivity file by every inclination.

"Numbers as passed to kerrz and used in file names: short, and the same every time."
fmt(x::Real) = @sprintf("%.6g", x)

"""
    kerrz_binary() -> String

The kerrz executable: `KERRZ_BIN` if set, otherwise `kerrz` or `kerrzcli` (the
name `pip install kerrz-cli` uses) from the PATH.
"""
function kerrz_binary()
    explicit = get(ENV, "KERRZ_BIN", "")
    if !isempty(explicit)
        isfile(explicit) || error("KERRZ_BIN is set to $(explicit), which does not exist")
        return explicit
    end
    for name in ("kerrz", "kerrzcli")
        path = Sys.which(name)
        path === nothing || return path
    end
    error("kerrz not found: put `kerrz` on the PATH or set KERRZ_BIN to its full path")
end

"`--nthreads n` if `KERRZ_NTHREADS` is set; otherwise kerrz picks the thread count."
function nthreads_args()
    n = get(ENV, "KERRZ_NTHREADS", "")
    return isempty(n) ? String[] : ["--nthreads", n]
end

"Seconds a single kerrz run may take before it is stopped (`KERRZ_TIMEOUT`, default 1800)."
kerrz_timeout() = parse(Float64, get(ENV, "KERRZ_TIMEOUT", "1800"))

"Where kerrz tables and emissivity files are kept, beside the grid caches."
function kerrz_files_dir()
    dir = joinpath(cache_directory(), "kerrz-files")
    mkpath(dir)
    return dir
end

"The end of a kerrz log, without its progress-bar escape codes."
function log_tail(path::AbstractString; nchars::Integer = 1500)
    isfile(path) || return ""
    text = read(path, String)
    text = replace(text, r"\e\[[0-9;]*[A-Za-z]" => "", '\r' => '\n')
    return last(strip(text), nchars)
end

function describe_exit(proc)
    proc.termsignal != 0 && return "killed by signal $(proc.termsignal)"
    return "exit code $(proc.exitcode)"
end

"""
    run_kerrz(args; dir)

Run `kerrz args...` in `dir`, with its output going to a log file. Throws with
the end of that log if kerrz fails, and stops it after `kerrz_timeout()` s
(kerrz 0.4.1 hangs rather than failing for some inputs, e.g. spin exactly 0).
"""
function run_kerrz(args::Vector{String}; dir::AbstractString)
    exe = kerrz_binary()
    logfile = joinpath(dir, "kerrz-$(getpid())-$(time_ns()).log")
    cmd = Cmd(`$exe $args`; dir = String(dir))
    try
        open(logfile, "w") do io
            proc = run(pipeline(cmd; stdout = io, stderr = io); wait = false)
            if timedwait(() -> process_exited(proc), kerrz_timeout(); pollint = 0.05) === :timed_out
                kill(proc)
                wait(proc)
                error("kerrz $(first(args)) did not finish within $(kerrz_timeout()) s and was stopped " *
                      "(KERRZ_TIMEOUT sets this limit). Command: $(cmd)")
            end
            success(proc) ||
                error("kerrz $(first(args)) failed ($(describe_exit(proc))). Command: $(cmd)\n" * log_tail(logfile))
        end
    finally
        rm(logfile; force = true)
    end
    return nothing
end

"""
    cached_kerrz_file(make, name) -> path

The file `name` in `kerrz_files_dir()`, creating it first with
`make(dir, tmpname)` if needed. `make` writes to a temporary name that is
renamed into place, so an interrupted run never leaves a partial file behind.
"""
function cached_kerrz_file(make, name::AbstractString)
    dir = kerrz_files_dir()
    path = joinpath(dir, name)
    isfile(path) && return path
    tmp = joinpath(dir, "partial-$(getpid())-$(time_ns())-$(name)")
    try
        make(dir, tmp)
        isfile(tmp) || error("kerrz did not write $(tmp)")
        mv(tmp, path; force = true)
    finally
        rm(tmp; force = true)
    end
    return path
end

"""
    transfer_table(spin, inc) -> path

The kerrz transfer-function table (`kerrz table`) for this spin and observer
inclination in degrees, out to `R_OUT`. About a second to compute.
"""
function transfer_table(spin::Real, inc::Real)
    name = "tf_a$(fmt(spin))_i$(fmt(inc))_rout$(fmt(R_OUT)).fits"
    return cached_kerrz_file(name) do dir, out
        run_kerrz(["table", "--spin", fmt(spin), "--incl", fmt(inc), "--rout", fmt(R_OUT),
                   "-o", out, nthreads_args()...]; dir)
    end
end

"""
    emissivity_file(spin, corona_flag, properties, nphotons) -> path

The emissivity profile (`kerrz emissivity`) of one corona, e.g.
`emissivity_file(0.9, "--ring-like", "h:5,x:4", 100_000)`. For the extended
coronae `nphotons` is per ring of the corona.
"""
function emissivity_file(spin::Real, corona_flag::AbstractString, properties::AbstractString, nphotons::Integer)
    tag = replace(properties, ':' => "", ',' => "_")
    extra = isempty(EMISSIVITY_EXTRA_ARGS) ? "" : "_x" * string(hash(EMISSIVITY_EXTRA_ARGS); base = 16)
    name = "em_a$(fmt(spin))_$(lstrip(corona_flag, '-'))_$(tag)_n$(nphotons)" *
           "_G$(fmt(EMISSIVITY_PHOTON_INDEX))$(extra).fits"
    return cached_kerrz_file(name) do dir, out
        run_kerrz(["emissivity", "--spin", fmt(spin), String(corona_flag), String(properties),
                   "--nphotons", string(nphotons), "--photon-index", fmt(EMISSIVITY_PHOTON_INDEX),
                   EMISSIVITY_EXTRA_ARGS..., "-o", out, nthreads_args()...]; dir)
    end
end

"Read the two comma-separated columns (g, flux) that `kerrz lineprof` writes."
function read_lineprofile(path::AbstractString)
    g = Float64[]
    L = Float64[]
    for line in eachline(path)
        parts = split(strip(line), ',')
        length(parts) == 2 || continue
        x = tryparse(Float64, strip(parts[1]))
        y = tryparse(Float64, strip(parts[2]))
        (x === nothing || y === nothing) && continue
        push!(g, x)
        push!(L, y)
    end
    length(g) >= 2 || error("no line profile could be read from $(path)")
    return g, L
end

"Linear interpolation of (g, L) onto `grid`, zero outside the tabulated range."
function resample(g::Vector{Float64}, L::Vector{Float64}, grid::AbstractVector{Float64})
    out = zeros(length(grid))
    for (k, x) in pairs(grid)
        (x < g[1] || x > g[end]) && continue
        j = clamp(searchsortedlast(g, x), 1, length(g) - 1)
        t = (x - g[j]) / (g[j+1] - g[j])
        out[k] = max((1 - t) * L[j] + t * L[j+1], 0.0)
    end
    return out
end

"""
    line_profile(spin, inc, emissivity_path) -> Vector{Float64}

The line profile `L(g)` on `G_GRID` for a disc from the ISCO to `R_OUT`,
illuminated as described by the emissivity file. Not normalised here;
`RedshiftKernel` scales it to unit area.
"""
function line_profile(spin::Real, inc::Real, emissivity_path::AbstractString)
    table = transfer_table(spin, inc)
    dir = kerrz_files_dir()
    out = joinpath(dir, "lineprofile-$(getpid())-$(time_ns()).csv")
    try
        run_kerrz(["lineprof", "--using-table", table, "--emissivity", String(emissivity_path),
                   "-o", out, nthreads_args()...]; dir)
        g, L = read_lineprofile(out)
        profile = resample(g, L, G_GRID)
        sum(profile) > 0 || error("kerrz returned an empty line profile for spin=$(spin), inc=$(inc), " *
                                  "emissivity $(basename(emissivity_path))")
        return profile
    finally
        rm(out; force = true)
    end
end
