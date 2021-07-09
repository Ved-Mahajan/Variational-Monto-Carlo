using CSV
using DataFrames
using Plots

function ansatz(x,α)
	return exp(-α*x^2)
end

function probability_density(x,α)
	return (ansatz(x,α)^2)/(sqrt(pi/α))
end

function E_local(x,α)
	return α + (x^2)*(1/2 - 2*α^2)
end

function metropolis(N,α)
	L = 3/sqrt(2*α)
	x = rand(-L:0.000001:L)
	E = 0
	E2 = 0
	Eln_avg = 0
	ln_avg = 0
	rejection_ratio = 0
	#start the loop
	for i in 1:N
		x_t = x + 0.4*(2*rand()*L - L)
		if probability_density(x_t,α) >= probability_density(x,α)
			x = x_t
		else
			dummy = rand()
			if dummy < probability_density(x_t,α)/probability_density(x,α)
				x = x_t
			else
				rejection_ratio += 1/N
			end
		end

		E += E_local(x,α)/N
		E2 += E_local(x,α)^2/N
		Eln_avg += (-1*E_local(x,α)*(x^2))/N
		ln_avg += (-x^2)/N
	end
	return E,E2,Eln_avg,ln_avg,rejection_ratio
end

function iterate_alpha()
	α = 0.1
	alpha_iter = 20
	N = 500
	random_walkers = 200
	γ = 0.9
	energy = Vector{Float64}(undef,alpha_iter)
	alpha = Vector{Float64}(undef,alpha_iter)
	variance = Vector{Float64}(undef,alpha_iter)
	En_analytical = Vector{Float64}(undef,alpha_iter)
	Var_analytical = Vector{Float64}(undef,alpha_iter)
	for i in 1:alpha_iter
		E = 0
		E2 = 0
		dE_dalpha = 0
		Eln = 0
		ln = 0
		rejection_ratio = 0
		for j in 1:random_walkers
			E_met, E2_met, Eln_met, ln_met, rejections_met = metropolis(N,α)
			E += E_met/random_walkers
			E2 += E2_met/random_walkers
			Eln += Eln_met/random_walkers
			ln += ln_met/random_walkers
			rejection_ratio += rejections_met/random_walkers
		end
		E_analytical = α/2 + 1/8/α
		var_analytical = (1-4*α^2)^2/32/α^2

		# Define the next α
		# dE_dalpha = 2*(Eln-E*ln)
		α = α + 0.05

		# Update the arrays for plotting
		energy[i] = E
		alpha[i] = α
		variance[i] = E2 - E^2
		En_analytical[i] = E_analytical
		Var_analytical[i] = var_analytical
	end
	return energy,alpha,variance,En_analytical,Var_analytical
end

energy,alpha,variance,E_analytical,var_analytical = iterate_alpha()
CSV.write(datadir("sims/harmonic_oscillator.csv"),(Energy = energy, Variance = variance, Alpha = alpha))

energy_plot = plot()
plot!(energy_plot,energy,
	label = "calculated expectation of energy",
	title = "Energy Plots",
	xlabel = "timestep",
	ylabel = "<E>",
	linecolor = :blue,
)
plot!(energy_plot,E_analytical,label = "Analytical Energy")

safesave(plotsdir()*"/energy_ho.png",energy_plot)

variance_plot = plot()
plot(variance_plot,variance,
	title = "Variance in energy",
	label = false,
	ylabel = "Var(E)",
	xlabel = "timestep",
	linecolor = :blue
     )

safesave(plotsdir()*"/variance_ho.png",variance_plot)

energy_vs_alpha = plot()
scatter(energy_vs_alpha,alpha,energy,
	label = "Calculated energy",
	title = "Energy vs α",
	markercolor = :blue,
    xlabel = "α",
	ylabel = "<E>"
)

scatter!(energy_vs_alpha,alpha,
	E_analytical,
	label = "Calculated energy",
	title = "Energy vs α",
	markercolor = :red
)
safesave(plotsdir()*"/energy_vs_alpha_ho.png",energy_vs_alpha)
