# Compile the reference models and the kerrz corona models into one library.
# Run via ./build-julia.sh kerrz.
#
# Output, kept apart from the core and Gradus builds:
#   build/kerrz/lib/libjuliaxspec_models.{dylib,so}
#   build/kerrz/include/julia_init.h
#   build/kerrz/xspec/model.dat
#   build/kerrz/xspec/juliaxspec_wrappers.c

using PackageCompiler
using JuliaXSPEC
using JuliaXSPECModels
using KerrzXSPEC

const ROOT = dirname(@__DIR__)
include(joinpath(ROOT, "scripts", "stub_tbbmalloc_proxy.jl"))

const DEST = joinpath(ROOT, "build", "kerrz")

println("Compiling the reference models and the kerrz corona models into build/kerrz/")
println("This takes several minutes. If kerrz is on the PATH, one line profile per")
println("geometry is calculated at the default parameters and remembered.")
create_library(
    joinpath(ROOT, "KerrzXSPEC"),
    DEST;
    lib_name = "juliaxspec_models",
    force = true,
    include_transitive_dependencies = true,
    include_lazy_artifacts = true,
    precompile_execution_file = joinpath(ROOT, "scripts", "precompile_kerrz.jl"),
)

stub_tbbmalloc_proxy!(DEST)

xspec_dir = joinpath(DEST, "xspec")
mkpath(xspec_dir)
println("Writing $(xspec_dir)/model.dat for: ",
        join([m.name for m in registered_models()], ", "))
write_xspec_package(xspec_dir; julia_init_header = "../include/julia_init.h")
println("Done. Next: ./build-xspec.sh kerrz")
