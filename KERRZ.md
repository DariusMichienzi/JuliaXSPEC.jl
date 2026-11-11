# kerrz corona models for XSPEC

`KerrzXSPEC` turns the kerrz command-line tool into XSPEC models for four
corona geometries. Each geometry gives two models:

| Corona     | Convolution (blurs any spectrum) | xillverD-5 blurred | Corona parameters (after `a`, `inc`)              |
|------------|----------------------------------|--------------------|---------------------------------------------------|
| Lamppost   | `kzlamp`                         | `kzlampx`          | `h` height                                        |
| Ring       | `kzring`                         | `kzringx`          | `h` height, `x` ring radius                       |
| Disc-like  | `kzdisc`                         | `kzdiscx`          | `h` height, `rin` inner radius, `dr` radial extent |
| Umbrella   | `kzumb`                          | `kzumbx`           | `r` offset radius, `inner` angle, `width` (outer = inner + width) |

Every model starts with spin `a` (0.1–0.998) and observer inclination `inc`
(5–85°). Lengths are in GM/c², angles in degrees. The disc runs from the ISCO
to 1000 GM/c².

## How it works

A line profile is calculated by kerrz in three steps: a transfer-function table
for (spin, inc) with `kerrz table`, the corona's emissivity profile with
`kerrz emissivity`, and the integration of the two with
`kerrz lineprof --using-table … --emissivity-profile …`. The first two are kept
as files and shared between grid points; the third takes about 0.05 s.

Profiles are cached on a grid (`KerrzXSPEC/src/geometries.jl`) and
interpolated between grid points, exactly like `jllamp`/`jlring`. A fit only
calculates the grid corners it visits, and every corner is remembered on disk
for later sessions. `scripts/kerrz_prewarm.jl` fills a whole grid in advance.

## Requirements

kerrz must be runnable wherever XSPEC runs: on the PATH as `kerrz` (or
`kerrzcli`, the name `pip install kerrz-cli` uses), or named by `KERRZ_BIN`.
Built and checked against kerrz 0.4.1. Older versions may lack the extended
coronae or `lineprof --using-table`; check `kerrz --help emissivity` and
`kerrz --help lineprof`.

Environment variables, set before starting XSPEC:

| Variable                | Meaning                                                    |
|-------------------------|------------------------------------------------------------|
| `KERRZ_BIN`             | Full path to the kerrz executable                          |
| `KERRZ_NTHREADS`        | Threads per kerrz run (default: kerrz's own choice)        |
| `KERRZ_TIMEOUT`         | Seconds before a kerrz run is stopped (default 1800)       |
| `JULIAXSPEC_CACHE_DIR`  | Where profiles and kerrz files are kept                    |
| `JULIAXSPEC_TABLE_DIR`  | Folder holding `xillverD-5.fits`, for the `…x` models      |

## Build

```bash
./build-julia.sh kerrz
./build-xspec.sh kerrz
```

then in XSPEC `lmod juliaxspec <repo>/build/kerrz/xspec`.

## Checking the profiles (outside XSPEC)

```bash
julia --project=builds/kerrz scripts/kerrz_check.jl            # all four
julia --project=builds/kerrz scripts/kerrz_prewarm.jl lamp ring # fill grids
```

`kerrz_check.jl` reports each kernel at the defaults and the interpolation
error halfway between grid nodes. If that error matters for your data, add
nodes to the relevant axis in `geometries.jl`; only that geometry's profiles
are recalculated.

## Using the convolution models

A convolution model only sees the energies on XSPEC's grid, so widen the grid
first or photons that should be shifted into the band are missing:

```
XSPEC12> energies 0.05 1000 5000 log
XSPEC12> model kzring * atable{/path/to/my_restframe_table.fits}
```

In the `…x` models the kernel parameters come first, then xillver's Gamma,
A_Fe, logXi, Dens and Incl. xillver's own inclination is therefore parameter 8
for `kzlampx`, 9 for `kzringx`, and 10 for `kzdiscx` and `kzumbx`. To make it
follow the disc inclination, tie it to parameter 2, e.g. `newpar 8 = 2` for
`kzlampx`. xillverD-5 only covers Incl ≥ 18.2°.

## Settings that change the physics

At the top of `KerrzXSPEC/src/KerrzXSPEC.jl`: the disc's outer radius, the
photon index used for the emissivity, extra `kerrz emissivity` arguments (for
example a velocity profile for the corona), and photons per emissivity run.
Changing any of them changes the cache names, so old and new profiles are never
mixed.

## Limitations

- One rest-frame spectrum is blurred by one kernel, so ionisation, density and
  emission angle are the same over the whole disc.
- Interpolating between profiles at neighbouring grid points can smear sharp
  features (the blue edge, the ISCO cut-off) when the grid is coarse. Check
  with `kerrz_check.jl` and refine the grid where needed.
- Line profiles are normalised to unit area, so the reflection strength of a
  geometry is not carried; the model norm absorbs it.
