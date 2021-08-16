using Plots
using Dates: monthabbr
using RollingFunctions
using BSON

monthTicks = (rollmean(range(1,365; length = 13), 2),
              "\$" .* monthabbr.(1:12) .* "\$")
monthLines = (range(1, 365; length = 13), fill("", 13))
begin
  p1 = plot(S;
    label = "Daily Health Impact Cost",
    α = 0.8,
    xlabel = "Month", ylabel = "Cost (in kUSD)",
    legendfontsize = 11, tickfontsize = 11,
    xticks = monthTicks, grid = :y)
  plot!(p1, cumsum(S) ./ collect(1:365);
    label = "Cummulative Health Impact Cost",
    lw = 5, α = 0.6)
    ylims!(p1, -8e3, 4.5e3)
  p2 = twinx(p1)
  plot!(p2, 1; xticks = monthLines, grid = :x,
    showaxis = :x, legend = false)
end
