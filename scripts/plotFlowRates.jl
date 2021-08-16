using DrWatson
@quickactivate "TESParetoHourly"
using TESParetoHourly
using JuMP: all_variables, start_value, set_start_value
using RollingFunctions: runmean, rollmean
using JLD2
using Printf, Statistics
using StatsPlots; pyplot()
using Dates: monthabbr
using Lazy: @>, @>>, @as

pal = palette(:default)

function generate_monthwise_plot(year)
  m = create_base_model(2009, 20, datadir(), 15)
  JLD2.@load datadir("sol_years", "sol2_$(year)_20_15.0_TES.jld2") varsVal objVal
  set_start_value.(m |> all_variables, varsVal)
  F = round.(start_value.(m[:F]); digits = 3)
  pvtCost, socCost = objVal
  area = sum(start_value.(m[:Aᵣ]))
  println("Objectives for ", year, " ", objVal)
  println("For year ", year, " ", "area is ", area)
  F = F + 0.5 .* iszero.(F)
  Fmonth = @> F begin
    reshape(:, 12)
    mean(dims = 1)
    dropdims(dims = 1)
    repeat(inner = 730)
  end
  xticksLab = (rollmean(range(1, 8760; length = 13), 2), monthabbr.(1:12))
  p1 = plot(
            F; α = 0.4, label = "Hourly Production Rate",
            ylabel = "F Cl2 (ton/hr)", guidefontsize = 12, size = (625, 180),
            legend = true, color = pal[1],
           )
  plot!(
        p1, Fmonth; linewidth = 2, linestyle = :dash, color = pal[1],
        label = "Mean Monthly Production Rate", xticks = xticksLab,
        yticks = ([0.5:1:3.5;], string.([0.0; 1.5:1:3.5])), grid = :y,
        tickfontsize = 12, legendfontsize = 12,
       )
  hline!(
         p1, [2.5766]; linewidth = 1.5, color = :black, legend = true,
         linestyle = :dot, label = "Mean Annual Production Rate",
        )
  p2 = twinx(p1)
  plot!(
        p2, 1; legend = false, grid = :x, showaxis = :x,
        xticks = (range(1, 8760; length = 13), fill("", 13)),
        gridstyle = :dash,
       )
  savefig(plotsdir("compPlots", "p$(year).svg"))
  return
end

# generate_monthwise_plot.(2008:2010)
