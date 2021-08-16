using BSON, Plots; pyplot()
using Dates: monthabbr
using RollingFunctions

monthTicks = (rollmean(range(1, 8760; length = 13), 2),
              "\$" .* monthabbr.(1:12) .* "\$")
monthLines = (range(1, 8760; length = 13), fill("", 13))

BSON.@load "data/met_files/NOx2009.bson" C₀
BSON.@load "data/SIinfo/zero_pub.bson" Cf

C0 = C₀[:, 2] ./ 1e3
Cf = Cf[:, 2]

p1 = plot(Cf - C0;
          size = (560, 220), framestyle = :zerolines,
          label = "Cf- C0",
          α = 0.8, xlabel = "Month",
          ylabel = "3Conc. (in ug/m3)",
          legendfontsize = 11, tickfontsize = 11,
          xticks = monthTicks, grid = :y
          )
p2 = twinx(p1)
plot!(p2, 1; xticks = monthLines, grid = :x,
      showaxis = :x, legend = false)
savefig("plots/zeroConc.svg")
