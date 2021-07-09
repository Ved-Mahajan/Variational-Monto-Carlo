using CSV
using DataFrames
using LinearAlgebra
using Plots


function ansatz(r1,r2,α)
    norm_r1 = norm(r1)
    norm_r2 = norm(r2)
    r12 = norm(r1-r2)
    wf = exp(-2*norm_r1)*exp(-2*norm_r2)*exp(r12/(2*(1+α*r12)))
    return wf
end

function prob_density(r1,r2,α)
    return ansatz(r1,r2,α)^2
end

function E_local(r1,r2,α)
    norm_r1 = norm(r1)
    norm_r2 = norm(r2)
    r12 = norm(r1-r2)
    dot_product = dot(r1/norm_r1-r2/norm_r2,r1-r2)
    energy = -4+dot_product/(r12*(1+α*r12)^2)-1/(r12*(1+α*r12)^3)-1/(4*(1+α*r12)^4)+1/r12
    return energy
end

function metropolis(N, α)

    L = 1
    r1 = rand(3)*2*L .- L
    r2 = rand(3)*2*L .- L #random number from -L to L
    E = 0
    E2 = 0
    Eln_average = 0
    ln_average = 0
    rejection_ratio = 0
    step = 0
    max_steps = -1


    for i in 1:N
        chose = rand()
        step = step + 1
        if chose < 0.5
            r1_trial = r1 + 0.5*(rand(3)*2*L .- L)
            r2_trial = r2
        else
            r2_trial = r2 + 0.5*(rand(3)*2*L .- L)
            r1_trial = r1
        end
        if prob_density(r1_trial,r2_trial,α) >= prob_density(r1,r2,α)
            r1 = r1_trial
            r2 = r2_trial
        else
            dummy = rand()
            if dummy < prob_density(r1_trial,r2_trial,α)/prob_density(r1,r2,α)
                r1 = r1_trial
                r2 = r2_trial
            else
                rejection_ratio += 1/N
            end
        end

        if step > max_steps
            E += E_local(r1,r2,α)/(N-max_steps)
            E2 += E_local(r1,r2,α)^2/(N-max_steps)
            r12 = norm(r1-r2)
            Eln_average += (E_local(r1,r2,α)*-r12^2/(2*(1+α*r12)^2))/(N-max_steps)
            ln_average += -r12^2/(2*(1+α*r12)^2)/(N-max_steps)
        end
    end
    return E, E2, Eln_average, ln_average, rejection_ratio
end


function iterate_alpha()
	α = 0

	alpha_iterations = 15
	N_metropolis = 5000
	random_walkers = 200
	γ = 0.5

	energy_plot = []
	alpha_plot = []
	variance_plot = []

	for i in 1:alpha_iterations
	    E = 0
	    E2 = 0
	    dE_dalpha = 0
	    Eln = 0
	    ln = 0
	    rejection_ratio = 0

	    for j in 1:random_walkers
	        E_met, E2_met, Eln_met, ln_met, rejections_met = metropolis(N_metropolis,α)
	        E += E_met/random_walkers
	        E2 += E2_met/random_walkers
	        Eln += Eln_met/random_walkers
	        ln += ln_met/random_walkers
	        rejection_ratio += rejections_met/random_walkers
	    end

		α = α + 0.05


	    energy_plot = append!(energy_plot, E)
	    alpha_plot = append!(alpha_plot, α)
	    variance_plot = append!(variance_plot, E2-E^2)
	end
	return energy_plot,alpha_plot,variance_plot
end

energy,alpha,variance = iterate_alpha();
CSV.write(datadir("sims/helium_atom.csv"),(Energy = energy, Variance = variance, Alpha = alpha))


energy_plot = plot()
plot!(energy_plot,energy,
	label = false,
	title = "Energy expenctation Plot",
	linecolor = :blue,
	ylabel = "<E>",
	xlabel = "timestep"
)
safesave(plotsdir("energy_helium_atom.png"),energy_plot)

variance_plot = plot()
plot!(variance_plot,variance,
	title = "Variance in energy",
	label = false,
	linecolor = :blue,
	ylabel = "Var(E)",
	xlabel = "timestep"
)
safesave(plotsdir("variance_helium_atom.png"),variance_plot)

alpha_plot = plot()
plot!(alpha_plot,alpha,
	title = "Variation Parameter",
	label = false,
	linecolor = :blue,
	ylabel = "α",
	xlabel = "timestep"
)
safesave(plotsdir("alpha_helium_atom.png"),alpha_plot)


alpha_vs_energy_plot = plot()
scatter!(alpha_vs_energy_plot,alpha,
	energy,
	title = "Energy vs α",
	xlabel = "α",
	ylabel = "<E>",
	label = false,
	markercolor = :blue
 )
safesave(plotsdir("alpha_vs_energy_helium_atom.png"),alpha_vs_energy_plot)
