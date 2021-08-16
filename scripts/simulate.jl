using DrWatson
@quickactivate "TESParetoHourly"
using TESParetoHourly
using JLD2
using RollingFunctions: runmean, rollmean
using DataFrames

function calc_rawdemand(F)
  Fc = 2.57657084 #Changeover point
  Fₗ = 0.2*2.7e4*71*3.6*0.16/(2*96485) #ton/hr of Chlorine Production Rate
  E_α1 = 1.018095
  E_β1 = 3.961352
  E_α2 = -10.90876
  E_β2 = 8.5903199
  E2m = 1.88/393.11 #kg/kg
  calCoal = 30274e-3 #MJ/kg #13000 BTU/lb
  energy = if F < Fₗ
    0.0
  elseif F < Fc
    E_α1 + E_β1*F
  else
    E_α2 + E_β2*F
  end
  demand = energy * 3.6e3/calCoal*E2m
  return demand
end

function calc_finalconc(F, w, Aᵣ)
  Aₜ = 15.0 #km2
  C₀ = w[1:2]
  Vd = w[3:4]
  H = w[5]
  η = w[6:7]
  Dₑ = η .* calc_rawdemand(F)
  Cᶠ = (C₀.*Aₜ.*H .+ Dₑ) ./ (Aᵣ.*Vd .+ Aₜ*H)
  return Cᶠ
end

function calc_violations(F, W, Aᵣ)
  C₀ = W[1:2, :]
  C₀R = runmean(C₀[2, :], 8)
  C₀All = [C₀; C₀R']
  CᵁAll = ceil.(Int, mapreduce(ci2cu, hcat, eachcol(C₀All)))
  Cᶠ = mapreduce(hcat, F, eachcol(W)) do f, w
    calc_finalconc(f, w, Aᵣ)
  end
  CᶠAll = ceil.(Int, [Cᶠ; runmean(Cᶠ[2, :], 8)'])
  violations = vec(reduce(|, CᶠAll .> CᵁAll; dims = 1))
  df = DataFrame(
                F = F, bad = violations,
                C_NO2_0 = C₀[1, :], C_O3_0 = C₀[2, :],
                C_NO2_f = Cᶠ[1, :], C_O3_f = Cᶠ[2, :],
               )
  return df
end

function ci2cu(Conc::AbstractVector; AQI_diff = 10)
  mapreduce(ismissing, |, Conc) && return similar(Conc, Missing)
  pollutants = reshape(["NO2", "O3", "O38hr"], size(Conc))
  AQI_In = TESParetoHourly.c2a.(Conc, pollutants) |> maximum
  AQI_Out = AQI_In ≤ 50 ? min(50, AQI_In + AQI_diff) : AQI_In
  Conc_out = TESParetoHourly.a2c.(AQI_Out, pollutants)
  Conc_out = max.(Conc, Conc_out)
end
