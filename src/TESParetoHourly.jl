module TESParetoHourly

export create_base_model, base_tradeoff, eps_con_opt
export AQIOzone1hr, AQIOzone8hr, AQINO2
export ConcOzone1hr, ConcOzone8hr, ConcNO2


import Gurobi
using CSV, Printf
using BSON
using JuMP
using Sobol
using Lazy: @>
using RollingFunctions: runmean
using Statistics: mean

include("model_builder.jl")
include("aqiConverter.jl")

"""

    base_tradeoff(m::JuMP.Model, pvtCost::Bool)
Calculates the extreme pareto point in one direction. If `pvtCost` is `true`, find minimal private cost, otherwise will find minimal public cost
"""
function base_tradeoff(m::JuMP.Model, pvtCost::Bool; saveat::Real=1200)
    #Deciding order of Objectives
    objType, firstObj, secondObj = begin
        if pvtCost
            "pvt.bson", m[:Z_pvt], m[:Z_pub]
        else
            "pub.bson", m[:Z_pub], m[:Z_pvt]
        end
    end
    m[:simType] = join([m[:simType], objType], "_")
    # Set up objectives and optimize
    @objective(m, Min, firstObj)
    optimize_with_save!(m; saveat=saveat)
    set_upper_bound(firstObj, firstObj |> JuMP.value)
    @objective(m, Min, secondObj)
    optimize_with_save!(m; saveat=saveat)
    set_upper_bound(secondObj, secondObj |> JuMP.value)
    # Fixing done to aid variation minimization problem
    fix(m[:Aᵣ][end], m[:Aᵣ][end] |> value; force = true)
    fix(m[:Sₜ], m[:Sₜ] |> value; force = true)
    set_optimizer_attribute(m, "MIPFocus", 1) #Force it to find improving feasible solution # Dont care about optimality here
    @objective(m, Min, m[:Z_flow])
    optimize_with_save!(m; interrupt = false, saveat = saveat) # This optimization problem stays way too long in presolve stage
    return (getindex.(m, [:Z_pvt, :Z_pub]) .|> value)
end #base_tradeoff

"""

    eps_con_opt(m::JuMP.Model, sampleNum::Integer, totSample::Integer; saveat::Real=1200.0)
Given a constraint on public cost, the method tries to find optimal in private cost
"""
# function eps_con_opt(m::JuMP.Model, pubConstraint::Real; saveat::Real=1200.0)
function eps_con_opt(m::JuMP.Model, sampleNum::Integer, totSample::Integer; saveat::Real=1200.0)
    @assert sampleNum ≤ totSample "The sample number should be smaller than total sample size."
    tradeoffTable = load_start_point(m) #warmstart
    # Create sequence and find random sample
    s = SobolSeq(tradeoffTable[2]...)
    println("Running sim for sample $(sampleNum) of $(totSample)")
    toSkipN = 1 << (totSample |> log2 |> floor |> Int)
    for i = 1:(toSkipN + sampleNum - 1)
        next!(s)
    end
    # skip(s, toSkipN + sampleNum - 1) # Skip initial samples as prescribed in Joe and Kuo 2003
    pubConstraint = next!(s)[1]
    println("Public Constraint is $(pubConstraint)")

    fname = join(["sol", @sprintf("%02d", sampleNum), "of", @sprintf("%02d", totSample)], "_") * ".bson"
    m[:simType] = joinpath(m[:simType], fname)

    @objective(m, Min, m[:Z_pvt])
    set_upper_bound(m[:Z_pub], pubConstraint)
    optimize_with_save!(m; saveat=saveat)
    set_upper_bound(m[:Z_pvt], value(m[:Z_pvt]))
    # Fixing done to aid variation minimization problem
    fix(m[:Aᵣ][end], m[:Aᵣ][end] |> value; force = true)
    fix(m[:Sₜ], m[:Sₜ] |> value; force = true)
    set_optimizer_attribute(m, "MIPFocus", 1) #Force it to find improving feasible solution # Dont care about optimality here
    @objective(m, Min, m[:Z_flow])
    optimize_with_save!(m; interrupt = false, saveat = saveat) # This optimization problem stays way too long in presolve stage
    return (getindex.(m, [:Z_pvt, :Z_pub]) .|> value)
end #eps_con_opt


"""

    optimize_with_save!(m::JuMP.Model; interrupt::Bool = true, saveat::Real=1200.0)
    Optimizes a model while saving all the variables at `saveat` seconds. Save path is assumed to be present at `m[:simType]`. If `interrupt` is false, doesn't interupt the optimization until original timelimit runs out (falls back to behavior of `JuMP.optimize! |> save_vars`)
"""
function optimize_with_save!(m::JuMP.Model; interrupt::Bool = true, saveat::Real = 1200.0)
    # Setting up optimizer with saving in between
    timeLim = get_optimizer_attribute(m, "TimeLimit")
    timeLim_list = fill(saveat, ceil(Int, timeLim/saveat))
    timeLim_list[end] -= sum(timeLim_list) - timeLim
    timeLim_list = interrupt ? timeLim_list : [timeLim]
    #The next part ensures that atleast one feasible solution to the problem is stored
    if !interrupt
        MIPTol = get_optimizer_attribute(m, "MIPGap")
        set_optimizer_attribute(m, "MIPGap", 1.0) # Just want a feasible solution to start with
        optimize!(m)
        save_vars(m)
        set_optimizer_attribute(m, "MIPGap", MIPTol) # Reset it to original tolerance
    end
    for tl in timeLim_list
        set_optimizer_attribute(m, "TimeLimit", tl)
        optimize!(m)
        save_vars(m)
        if has_values(m) && (m |> relative_gap |> iszero)
            break
        end
    end
    has_values(m) || error("Model was not solved correctly")
    set_optimizer_attribute(m, "TimeLimit", timeLim) # Reset time value to original value
end #optimize_with_save!

"""

    save_vars(m::JuMP.Model) -> Nothing
Save variables at path `m[:simType]`
"""
function save_vars(m::JuMP.Model)::Nothing
    objs = getindex.(m, [:Z_pvt, :Z_pub])
    if termination_status(m) == MOI.OPTIMAL
        varsVal = m |> all_variables .|> JuMP.value
        objVal = objs .|> JuMP.value
        BSON.@save m[:simType] varsVal objVal
    elseif termination_status(m) == MOI.TIME_LIMIT && has_values(m)
        varsVal = m |> all_variables .|> JuMP.value
        objVal = objs .|> JuMP.value
        BSON.@save m[:simType] varsVal objVal
    end
    return nothing
end

"""

    load_start_point(m::Model) -> tradeoffTable
This function helps warm start the solver. Use this function only when tradeoff table is available i.e. for ϵ-constraint optimization only. It also assigns bounds to objectives using tradeoff_table.
"""
function load_start_point(m::Model)
    objs = getindex.(m, [:Z_pvt, :Z_pub])
    BSON.@load joinpath(m[:simType], "warmstart.bson") varsVal varsName
    tradeoffTable = CSV.read(joinpath(m[:simType], "tradeoff.csv"); header = false) |> Matrix
    @assert all((m |> all_variables .|> name) .== varsName) "Stored variable names don't match"
    set_start_value.(m |> all_variables, varsVal)
    set_lower_bound.(objs, (@> tradeoffTable minimum(dims = 1) dropdims(dims = 1)))
    set_upper_bound.(objs, (@> tradeoffTable maximum(dims = 1) dropdims(dims = 1)))
    return extrema(tradeoffTable; dims = 1)
end #load_start_point

end # Module
