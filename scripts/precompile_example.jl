using DrWatson
@quickactivate "TESParetoHourly"

using TESParetoHourly
using JuMP
using ArgParse
using Printf

function parse_commandline(ARGS)
    s = ArgParseSettings(description = "This script finds pareto optimal intermediate points")

    @add_arg_table! s begin
        "--threads", "-t"
            help = "Number of threads"
            arg_type = Int
            default = 1
        "metYear"
            help = "Year for which meteorological and concentration data is used"
            arg_type = Int
            required = true
        "treeAge"
            help = "Age of trees that are planted"
            arg_type = Int
            required = true
        "sampleNum"
            help = "Index of sample being collected"
            arg_type = Int
            required = true
        "totSample"
            help = "Total number of ssamples taken for pareto points"
            arg_type = Int
            required = true
        "airShedArea"
            help = "Area of air shed in sq. km"
            arg_type = Float64
            default = 15.0
        "--datadir", "-d"
            help = "Path where data will be read from and results will be written"
            arg_type = String
            default = datadir()
        "--timelim","-l"
            help = "Time limit for each optimization run in hrs"
            arg_type = Float64
            default = 4.0
        "--tol", "-g"
            help = "MIP Tolerance Gap"
            arg_type = Float64
            default = 0.0
        "--private", "-p"
            help = "Flag to calculate lowest private cost"
            action = :store_true
        "--mode", "-m"
            help = "Solution mode. 1 for TES, 2 for TechOnly, 3 for EcoOnly"
            arg_type = Int
            default = 1
    end
    return parse_args(ARGS, s)
end


ARGS = [2009, 20, 1, 20] .|> string
parsed_args = parse_commandline(ARGS)
metYear = parsed_args["metYear"]
treeAge = parsed_args["treeAge"]
modeType = parsed_args["mode"]
tol = parsed_args["tol"]
timelim = parsed_args["timelim"] * 3600
Aₜ = parsed_args["airShedArea"]

@assert (2006 ≤ metYear ≤ 2015) "Meteorological data is available only for year 2006 to 2015"
@assert (1 ≤ treeAge ≤ 40) "Age of reforested trees should be between 1 to 40 years."
@assert (modeType in 1:3) "Only three modes are possible. Please select any integer from 1 to 3"
@assert (0.0 ≤ tol ≤ 1.0) "Tolerance can only be a fraction between 0 and 1"
@assert timelim > 0 "Time limit needs to be positive"

# Normal Tradeoff set up part
m = create_base_model(metYear, treeAge, parsed_args["datadir"], Aₜ;
                      threads = parsed_args["threads"], MIPTol = tol, TimeLim = timelim)
# Name the mode and set Aᵣ and Sₜ to appropiate values
mode_txt = if modeType==1
    "TES"
elseif modeType==2
    JuMP.fix.(m[:Aᵣ], 0.0; force = true)
    "Tech"
else modeType==3
    JuMP.fix(m[:Sₜ], 0.0; force = true)
    "Eco"
end

foldername = join(["pareto", metYear, treeAge, Aₜ, mode_txt], "_")
folderpath = joinpath(parsed_args["datadir"], foldername)
mkpath(folderpath)
m[:simType] = joinpath(parsed_args["datadir"], foldername, "extreme")
objValue = base_tradeoff(m, parsed_args["private"])

open(joinpath(folderpath, "tradeoff.csv"), "a") do io
    join(io, objValue, ",")
    print(io, "\n")
end

# # ϵ-constraint part
# m = create_base_model(metYear, treeAge, parsed_args["datadir"], Aₜ;
#                       threads = parsed_args["threads"], MIPTol = tol, TimeLim = timelim)
# # Name the mode and set Aᵣ and Sₜ to appropiate values
# mode_txt = if modeType==1
#     "TES"
# elseif modeType==2
#     JuMP.fix.(m[:Aᵣ], 0.0; force = true)
#     "Tech"
# else modeType==3
#     JuMP.fix(m[:Sₜ], 0.0; force = true)
#     "Eco"
# end

# foldername = join(["pareto", metYear, treeAge, Aₜ, mode_txt], "_")
# folderpath = joinpath(parsed_args["datadir"], foldername)
# mkpath(folderpath)
# m[:simType] = joinpath(parsed_args["datadir"], foldername)
# objValue = eps_con_opt(m, parsed_args["sampleNum"], parsed_args["totSample"])

# open(joinpath(folderpath, "pareto.csv"), "a") do io
#     print(io, "$(parsed_args["sampleNum"]),")
#     join(io, objValue, ",")
#     print(io, "\n")
# end
