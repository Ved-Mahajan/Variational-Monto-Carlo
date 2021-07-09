using CSV
using DataFrames
using LinearAlgebra
using Plots

function ansatz(x,α)
	x = norm(x)
	return exp(-α*x)
end

function probability_density(x,α)
	# note that it is not normalised becuase it does not matter!
	return ansatz(x,α)^2
end

function E_local(x,α)
	x = norm(x)
	return -1/x - (α/2)*(α - 2/x)
end

function metropolis(N,α)
	L = 3/sqrt(2*α)
	x = rand(3)*2*L .- L
	E = 0
	E2 = 0
	Eln_avg = 0
	ln_avg = 0
	rejection_ratio = 0
	#start the loop
	for i in 1:N
		x_t = x + 0.2*(rand(3)*2*L .- L)
		if probability_density(x_t,α) >= probability_density(x,α)
			x = x_t
		else
			rand_number = rand()
			if rand_number < probability_density(x_t,α)/probability_density(x,α)
				x = x_t
			else
				rejection_ratio += 1/N
			end
		end

		E += E_local(x,α)/N
		E2 += E_local(x,α)^2/N
		Eln_avg += -1 * E_local(x,α) * (norm(x)/N)
        ln_avg += -1 * norm(x)/N

	end
	return E,E2,Eln_avg,ln_avg,rejection_ratio
end

function iterate_alpha()
	α = 0.2
	alpha_iter = 30
	N = 500
	random_walkers = 200
	γ = 0.5
	energy = Vector{Float64}(undef,alpha_iter)
	alpha = Vector{Float64}(undef,alpha_iter)
	variance = Vector{Float64}(undef,alpha_iter)
	for i in 1:alpha_iter
		E = 0
		E2 = 0
		dE_dalpha = 0
		Eln = 0
		ln = 0
		rejection_ratio = 0
		for j in 1:random_walkers
			E, E2, Eln, ln, rejections = metropolis(N,α)
			E += E/random_walkers
			E2 += E2/random_walkers
			Eln += Eln/random_walkers
			ln += ln/random_walkers
			# rejection_ratio += rejections_met/random_walkers
		end

		# Define the next α
		dE_dalpha = 2*(Eln - E*ln)
		α = α - γ*dE_dalpha
		# α = α + 0.05

		# Update the arrays for plotting
		energy[i] = E
		alpha[i] = α
		variance[i] = E2 - E^2

	end
	return energy,alpha,variance
end

energy,alpha,variance = iterate_alpha()
CSV.write(datadir("sims/hydrogen_atom.csv"),(Energy = energy, Variance = variance, Alpha = alpha))


energy_plot = plot()
plot!(energy_plot,energy,
	label = false,
	title = "Energy expenctation Plot",
	linecolor = :blue,
	ylabel = "<E>",
	xlabel = "timestep"
 )

safesave(plotsdir("energy_hydrogen_atom.png"),energy_plot)

variance_plot = plot()
plot!(variance_plot,variance,
	title = "Variance in energy",
	label = false,
	linecolor = :blue,
	ylabel = "Var(E)",
	xlabel = "timestep"
     )
safesave(plotsdir("variance_hydrogen_atom.png"),variance_plot)


alpha_plot = plot()
plot!(alpha_plot,alpha,
	title = "Variation Parameter",
	label = false,
	linecolor = :blue,
	ylabel = "α",
	xlabel = "timestep"
)

safesave(plotsdir("alpha_hydrogen_atom.png"),alpha_plot)

alpha_vs_energy_plot = plot()
scatter!(alpha_vs_energy_plot,alpha,
	energy,
	title = "Energy vs α",
	xlabel = "α",
	ylabel = "<E>",
	label = false,
	markercolor = :blue
)
safesave(plotsdir("alpha_vs_energy_hydrogen_atom.png"),alpha_vs_energy_plot)
