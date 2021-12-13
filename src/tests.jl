using Plots

n3bia(x) = 2/3 + x
n3bib(x) = 8/3 - 4/x
n3bic(x) = 2/3 + 4/x

x = 1:1:10

scatter(n3bia.(x))
scatter!(n3bib.(x))
scatter!(n3bic.(x))
xlims!((-5,5))
ylims!((-50,50))