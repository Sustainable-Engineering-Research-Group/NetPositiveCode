"""
This script is used for plotting annual data on monthwise scale
"""
using Plots
using BSON: @load
using Dates: monthabbr
pyplot()

fname = ARGS[1]

@load fname monthrate badprod badAQIdays

p1 = bar(monthrate;
    fillalpha = 0.3,
    linecolor = :steelblue,
    bar_width = 1,
    legend = false,
    ylims = [1.1, 3.5],
    ylabel = "Mean Monthly Production Rate",
    grid = false,
    tickfontsize = 12,
#     legendfontsize = 12,
    xticks = (1:12, "\$" .* monthabbr.(1:12) .* "\$"))
hline!(p1, [2.5766];
    linewidth = 1.5,
    linestyle = :dashdot)
#     label = "Annual Average Target production rate")
p2 = twinx(p1)
plot!(p2, [], label = nothing)
plot!(p2, [];
    linestyle = :dashdot,
    label = "Target production rate")
plot!(p2, [badprod badAQIdays];
    label = ["Production Shutdown" "Bad AQI Days"],
    ylabel = "Number of Days",
    tickfontsize = 12,
    legendfontsize = 12,
    linewidth = 3,
    linealpha = 0.7,
    grid = false,
    framestyle = :box,
    xticks = nothing,
    ylims = [-0.25, 12.25])
replace!(fname, "data" => "plots")
savefig(p1, fname * ".svg")
