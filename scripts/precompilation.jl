import Pkg; Pkg.activate("..")

using PackageCompiler
PackageCompiler.create_sysimage(:TESParetoHourly; sysimage_path="../TESSysImage.so",
                               precompile_execution_file="precompile_example.jl")
