const Fₗ = 0.2*2.7e4*71*3.6*0.16/(2*96485) #ton/hr of Chlorine Production Rate
const Fc = 2.57657084 #Changeover point

"""

`create_base_model(y::Int64, tree_age::Int64, fpath::String, Aₜ::Real=15.0)`

Creates an yearly model at hourly scale. `y` represent the weather data year and `tree_age` is the size of tree planted. The model does not have a predefined objective. Although, one can use variables `[Z_pvt, Z_pub, Z_flow]` for optimization.

# Arguments
- `threads::Int64`=4 Threads provided to gurobi optimizer
- `MIPTol::Float64`=0.0 Relative Optimality gap
- `TimeLim::Real`=14400 Timelimit for optimization run in secs
"""
function create_base_model(y::Int64, tree_age::Int64, fpath::String, Aₜ::Real=15.0; threads::Int64=4,
                        MIPTol::Float64=0.0, TimeLim::Real=14400)
    GRB_ENV = Gurobi.Env();
    @assert 2006 ≤ y ≤ 2015 "The system is only configured for years 2006 to 2015"
    @assert 1 ≤ tree_age ≤ 20 "Tree age can only vary between 1 to 20 years old"
    tr = ["ae","lo","rm","sm","wa"]
    BSON.@load joinpath(fpath, "met_files", "NOx$(y).bson") C₀ Cᵁ
    BSON.@load joinpath(fpath, "met_files", "h$(y).bson") H
    BSON.@load joinpath(fpath, "met_files", "VdNOx$(y).bson") VNO2 VO3
    gf = BSON.load(joinpath(fpath, "met_files", "grow_norm.bson"))[:g][tree_age, :]
    #= gf[Not(4)] .= 0.0 #Setting SM as the only species option =#
    C₀ = C₀ ./ 1e3 #Converting in kg/km³
    Cᵁ = Cᵁ ./ 1e3 #Converting in kg/km³
    H = H ./ 1e3 #Converting in km
    VNO2 .= VNO2 .* gf'
    VO3 .= VO3 .* gf'
    VNO2 .= VNO2 .* (VNO2 .≥ 1e-4)
    VO3 .= VO3 .* (VO3 .≥ 1e-4)
    Tₙ = 365*24
    Tₛ = 90*24 #Summer Start
    Tₑ = 273*24 #Summer End
    pols = ["NO2", "O3"]
    # pols = ["NO2", "SO2", "PM10", "PM2.5"]
    # pols = ["NO2", "SO2", "PM10", "PM2.5", "O3","SO224hr", "O38hr"]
    # E2m = [1.88, 2.22, 0.63, 0.59] #kg/MWhr
    # E2m = [1.88, 2.22, 0.63, 0.59]/393.11 #kg/kgCoal
    E2m = 1.88/393.11 #kg/kg
    #E2m_summ = [0.1, 0, 19.62, 2.22, 2.22, 0, 0.63, 0.59] #kg/MWhr
    Cal_Coal = 30274e-3 #MJ/kg #13000 BTU/lb
    ηₒ = 0.1*ones(Tₙ) # NO₂ to O₃ production efficiency
    # ηₒ[Tₛ:Tₑ] .= 10
    ηₒ[filter(t ->  8 ≤ rem(t, 24) ≤ 19, Tₛ:Tₑ)] .= 10
    ηₙ = ones(Tₙ)
    # ηₙ[Tₛ:Tₑ] .= 0.1
    ηₙ[filter(t ->  8 ≤ rem(t, 24) ≤ 19, Tₛ:Tₑ)] .= 0.1
    Sₜᵁ = 20 #MW
    Sₜᴸ = 6 #MW
    E_α1 = 1.018095
    E_β1 = 3.961352
    E_α2 = -10.90876
    E_β2 = 8.5903199
    C_α = 58.37186
    C_β = 16.060
    O_β = 218.7694

    # Creating different metrics of initial concentration
    C₀D1m = @> C₀ begin
        reshape(24, :, 2)
        maximum(dims=1)
        dropdims(dims=1)
    end
    C₀R = (C₀ |> eachcol .|> (x->runmean(x, 8)) |> (x->reduce(hcat, x)))
    C₀D24a = @> C₀ begin
        reshape(24, :, 2)
        mean(dims=1)
        dropdims(dims=1)
    end

    # Creating lower bound for different metrics of final concentration
    Cᶠᴸ = similar(C₀)
    Cᶠᴸ[:,1] .= C₀[:,1] .* H ./ (H + VNO2[:,5])
    Cᶠᴸ[:,2] .= C₀[:,2] .* H ./ (H + VO3[:,5])
    CᶠᴸD1m = @> Cᶠᴸ begin
        reshape(24, :, 2)
        maximum(dims=1)
        dropdims(dims=1)
    end
    CᶠᴸR = (Cᶠᴸ |> eachcol .|> (x->runmean(x, 8)) |> (x->reduce(hcat, x)))
    CᶠᴸD24a = @> Cᶠᴸ begin
        reshape(24, :, 2)
        mean(dims=1)
        dropdims(dims=1)
    end

    # Creating upper bound for different metrics of final concentration
    CᵁD1m = @> Cᵁ begin
        reshape(24, :, 3)
        maximum(dims=1)
        dropdims(dims=1)
    end
    CᵁR = Cᵁ |> eachcol .|> (x->runmean(x, 8)) |> (x->reduce(hcat, x))
    CᵁD24a = @> Cᵁ begin
        reshape(24, :, 3)
        mean(dims=1)
        dropdims(dims=1)
    end

    # E_α = 0.0640375
    # E_β = 4.4334503
    # C_α = [58.37186, 1161.2, 0]
    # C_β = [16.060, 81.82195, 1.1151]
    # O_β = [218.7694, 158.506, 353]

    m = JuMP.Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GRB_ENV),
        "NonConvex" => 2, "Presolve" => -1, "MIPGap" => MIPTol, "TimeLimit" => TimeLim,
        "Threads" => threads, "ImproveStartGap" => 0.1)) #"OptimalityTol" => 1e-2, 
    # Find better feasible solutions at 10% gap
    # NO2 is pollutant 1
    # O3 is pollutant 2
    # O3 8hr is pollutant 3

    @variables(m, begin
            0 ≤ F[1:Tₙ] ≤ 3*Fₗ #Chlorine Flow Rate
            yF[1:2, 1:Tₙ], Bin #Energy segment MILP
            yPwr[1:Tₙ÷24], Bin #Batch Variable to turn on power or off
            0 ≤ ΔF[1:Tₙ] ≤ 0.25*2Fₗ #Change in flowrate
            F1[1:Tₙ] ≥ 0 #Flowrate in lower energy curve
            F2[1:Tₙ] ≥ 0 #Flowrate in higher energy curve
            Fweek[1:20] ≥ 0 #Cumulative flow over 20 weeks of 18 days
            Eₜ[1:Tₙ] ≥ 0 #Energy split for each pollutant towards tech
            Eₑ[1:Tₙ] ≥ 0 #Energy split for each pollutant towards Ecology
            Dₑ[1:2,1:Tₙ] ≥ 0 #Demand of ecosystem service for each pollutant
            Sₑ[1:2,1:Tₙ] ≥ 0 #Supply of ecosystem service for each pollutant
            Cᶠᴸ[t,i] ≤ Cᶠ[i=1:2,t=1:Tₙ] ≤ Cᵁ[t,i] #Final Concentration in area
            CᶠᴸR[t,i] ≤ CᶠR[i=1:2,t=1:Tₙ] ≤ CᵁR[t,i] #Final Concentration in area
            CᶠᴸD1m[t,i] ≤ CᶠD1m[i=1:2,t=1:(Tₙ÷24)] ≤ CᵁD1m[t,i]
            CᶠᴸD24a[t,i] ≤ CᶠD24a[i=1:2,t=1:(Tₙ÷24)] ≤ CᵁD24a[t,i]
            (CᶠᴸD1m[t,i] - C₀D1m[t,i]) ≤ ΔCD1m[i=1:2,t=1:(Tₙ÷24)] ≤ (CᵁD1m[t,i] - C₀D1m[t,i]) #Final - Intial to be minimized
            (CᶠᴸD24a[t,i] - C₀D24a[t,i]) ≤ ΔCD24a[i=1:2,t=1:(Tₙ÷24)] ≤ (CᵁD24a[t,i] - C₀D24a[t,i]) #Final - Intial to be minimized
            0 ≤ Aᵣ[1:5] ≤ Aₜ #Total Area under restoration
            Sₜ ≥ 0 #Size of technological unit
            Yₜ, Bin #Purchasing of technological unit
            yPCU[1:Tₙ÷24], Bin #Turning on or off the Tech PCU for a batch
            En[1:Tₙ] ≥ 0 #Energy requirement for chlorine production in MWhr
            S_Pwr ≥ 0 #Size of power generating boiler unit in MW
            C_Coal[1:Tₙ] ≥ 0 #Cost of Coal every time step
            C_Tech ≥ 0 #Capital cost of Tech
            O_Tech ≥ 0 #Operational cost of Tech
            Z_pvt ≥ 0 # Total Cost of Production
            Z_pub # Public Health Cost
            Z_flow ≥ 0 # Variance in flow
        end)

    #Chlorine Production Model and technological PCU
    @constraints(m, begin
            fLwrCurve1, yF[1, :] .* Fₗ .≤ F1 #Lower Part of Curve activating
            fLwrCurve2, F1 .≤ yF[1, :] .* Fc
            fhighCurve1, Fc .* yF[2, :] .≤ F2 #Upper part of curve activating
            fhighCurve2, F2 .≤ 3*Fₗ .* yF[2, :]
            fOnOff1[t=1:Tₙ], F[t] .≥ Fₗ .* yPwr[(t-1)÷24 + 1] #Is the generator on or off
            fOnOff2[t=1:Tₙ], F[t] .≤ 3*Fₗ .* yPwr[(t-1)÷24 + 1]
            PwrOnOff[t=1:Tₙ], sum(yF[:, t]) .== yPwr[(t-1)÷24 + 1] #Select power curve only if generator is on
            FPeriodTot, F .== F1 + F2 #Flow is sum of both the curves
            FweeklySum[t=1:20], Fweek[t]*100 .== sum(F[((t-1)*438+1):(t*438)]) #Accumalate flow over weeks and divide by 100
            Energy_Req, En .== E_α1.*yF[1, :] + F1.*E_β1 + E_α2.*yF[2, :] + F2.*E_β2 #Energy Req as linear function of chlorine production rate
            Tot_Prod, sum(Fweek) == Fc*Tₙ/100 #Average Chlorine production target needs to be met
            Utility_Size, S_Pwr .≥ En #Utility generator should be able to meet demand
            Tech_Eco_Split, En .== Eₜ .+ Eₑ #Split b\w Tech and Eco
            Tech_PCU_Size, Eₜ .≤ Sₜ # Max processing capacity is less than install cap
            Tech_Exist_LB, Sₜ .≥ Sₜᴸ .* Yₜ #Lower bound of semi-cont S
            Tech_Exist_UB, Sₜ .≤ Sₜᵁ .* Yₜ #Upper bound of semi-cont S
            Tech_PCU_Switch1[t=1:Tₙ], Eₜ[t] .≥ Sₜᴸ .* yPCU[(t-1)÷24 + 1] #Turn on or off SCR for each batch period
            Tech_PCU_Switch2[t=1:Tₙ], Eₜ[t] .≤ Sₜᵁ .* yPCU[(t-1)÷24 + 1]
            Ramp_Up_Limit[t=2:Tₙ], F[t] - F[t-1] ≤ ΔF[t] + (2 - yPwr[(t-1)÷24 + 1] - yPwr[(t-2)÷24 + 1])*Fₗ*0.5 #This equations hold only when the power equipment state is on for both the time steps considered. A recently turned on equipment will only operate at lowest flow rate
            Ramp_Dwn_Limit[t=2:Tₙ], F[t] - F[t-1] ≥ -ΔF[t] - (2 - yPwr[(t-1)÷24 + 1] - yPwr[(t-2)÷24 + 1])*Fₗ*0.5
            CoalCost, C_Coal .== En .* 3.6e3/Cal_Coal*0.06916
    #         Ramp_Up_Limit[t=2:Tₙ], F[t] - F[t-1] ≤ Fₗ/10.0
    #         Ramp_Dwn_Limit[t=2:Tₙ], F[t] - F[t-1] ≥ -Fₗ/10.0
            end)

    #AQI Constraints
    @constraints(m, begin
            Conc_Final[i=1:2, t=1:Tₙ], C₀[t,i]*Aₜ*H[t] + Dₑ[i,t] - Sₑ[i,t] == Cᶠ[i,t]*Aₜ*H[t]
            Conc_Roll_8hr[t=8:Tₙ], CᶠR[:,t] .== (@> Cᶠ[:,(t-7):t] mean(dims=2) dropdims(dims=2))
            Conc_Run_8hr[t=1:7], CᶠR[:,t] .== (@> Cᶠ[:,1:t] mean(dims=2) dropdims(dims=2))
            Conc_UB, Cᶠ .≤ Cᵁ[:,1:2]' #Maintain AQI
            Conc_UB_O38hr, CᶠR[2,:] .≤ Cᵁ[:,3] #Maintain AQI
            Conc_D1_max[d=1:365], CᶠD1m[:,d] .≥ Cᶠ[:, (24d-23):(24d)]
            Conc_D24_mean, CᶠD24a .== (@> Cᶠ reshape(2, 24, :) mean(dims=2) dropdims(dims=2))
            Conc_Delta_D1max, ΔCD1m .== CᶠD1m - C₀D1m'
            Conc_Delta_D24mean, ΔCD24a .==  CᶠD24a - C₀D24a'
            end)

    #Ecological Supply and Demand
    #Ecological Constraints
    @constraints(m, begin
            Demand_Eco_NO2, Dₑ[1,:] .== ηₙ .* Eₑ .* 3.6e3/Cal_Coal*E2m # Calculating Emmisions per hour NOx
            Demand_Eco_O3,  Dₑ[2,:] .== ηₒ .* Eₑ .* 3.6e3/Cal_Coal*E2m # Calculating Emmisions per hour Ozone Summer
            Supply_EcoNO2, Sₑ[1,:] .== Cᶠ[1,:] .* (VNO2 * Aᵣ) #Supply based on area of reforestation # Sₑ should be in gms
            Supply_EcoO3,  Sₑ[2,:] .== Cᶠ[2,:] .* (VO3 * Aᵣ) #Supply based on area of reforestation # Sₑ should be in gms
            Tot_Area, Aₜ .≥ sum(Aᵣ) #Total Area of Restoration
            end)

    # Economics of the system
    @constraints(m, begin
            Cost_Tech, C_Tech == C_α*Yₜ + C_β*Sₜ
            Oper_Tech, O_Tech == O_β*sum(Eₜ)/1000 #Converting the cost into 1k USD
            Tot_Cost, Z_pvt == C_Tech + O_Tech/Tₙ + 75*sum(Aᵣ)*0.095 + sum(C_Coal)/1e3 #+ S_Pwr*1200*0.085
            Health_impact, 2*Z_pub == mean(ΔCD1m[2,:])*1120 + mean(ΔCD24a[2,:])*1480
            Flow_Variance, Z_flow == sum(ΔF)
            end)

    #Creating starting points and a dummy soluton to warm start optimization
    JuMP.fix.(Aᵣ[1:4], 0.0; force = true)
    # @objective(m, Min, Z_pub)
    # JuMP.optimize!(m)
    return m
end
