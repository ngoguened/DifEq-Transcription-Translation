using DifferentialEquations
using Plots

function rates(constants)
    gene, transcript, protein = constants
    gene_rate = -d*protein
    transcript_rate = c*gene - v*transcript
    protein_rate = k*transcript - u*protein
    return [gene_rate, transcript_rate, protein_rate]
end

c = k = 0.1
u = v = .05
d = 0

time_span = (0,100)
intial_values = [2,0,0]

problem = ODEProblem(rates, intial_values, time_span)
solution = solve(problem)
plot(solution)
savefig("transcriptionTranslationDE.png") 