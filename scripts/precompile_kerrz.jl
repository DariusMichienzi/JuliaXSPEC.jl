# Exercised while building the kerrz library, so the code paths XSPEC uses are
# compiled ahead of time. Like precompile_models.jl, but every model is tried
# separately: a missing xillverD-5.fits or kerrz executable skips that model
# instead of failing the build.

using JuliaXSPEC, JuliaXSPECModels, KerrzXSPEC

const NEEDS_TABLE = ("jltable", "jltableblur", "kzlampx", "kzringx", "kzdiscx", "kzumbx")

edges = collect(logrange(0.1, 50.0, 1001))
for model in registered_models()
    try
        model.name in NEEDS_TABLE && resolve_table_path("xillverD-5.fits")
        values = parameter_defaults(model)
        if model isa ConvolutionModel
            evaluate(model, BinnedSpectrum(edges, gaussian_per_bin(edges, 6.4, 0.1)), values)
        else
            evaluate(model, edges, values)
        end
        flux = model isa ConvolutionModel ? gaussian_per_bin(edges, 6.4, 0.1) : zeros(length(edges) - 1)
        flux_error = zeros(length(edges) - 1)
        GC.@preserve edges values flux flux_error begin
            JuliaXSPEC.juliaxspec_evaluate(
                Base.unsafe_convert(Ptr{Cchar}, model.name), pointer(edges), Cint(length(edges) - 1),
                pointer(values), Cint(1), pointer(flux), pointer(flux_error), Base.unsafe_convert(Ptr{Cchar}, ""),
            )
        end
    catch err
        println("precompile: skipping $(model.name): ", first(sprint(showerror, err), 300))
    end
end
