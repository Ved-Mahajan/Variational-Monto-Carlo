### A Pluto.jl notebook ###
# v0.14.3

using Markdown
using InteractiveUtils

# ╔═╡ 028bce7d-e731-47e9-9468-84185721f0ab
using Pkg; Pkg.add(["Plots","LinearAlgebra"])

# ╔═╡ c912668e-a13f-11eb-0230-8be9b4bd6d15
using Plots,LinearAlgebra

# ╔═╡ f04b3d61-d68c-4349-ab98-7136acc75447
using CSV,DataFrames

# ╔═╡ aa7fc637-2e9d-4c14-a284-50294ba7cbe5
md"""
# Define all the functions
"""

# ╔═╡ 89e3a6ae-4f5c-4882-81fd-215379d609dd
function ansatz(x,α)
	x = norm(x)
	return exp(-α*x)
end;

# ╔═╡ 1197b58f-7614-4e7c-adeb-96af732b7916
function probability_density(x,α)
	# note that it is not normalised becuase it does not matter!
	return ansatz(x,α)^2
end;

# ╔═╡ db71e022-2374-45bc-98af-3fd5e7256784
function E_local(x,α)
	x = norm(x)
	return -1/x - (α/2)*(α - 2/x)
end;

# ╔═╡ c59f8de5-73db-495f-9a7a-a4b9bf0b58bc
md" # Function for running the metropolis algorithm."

# ╔═╡ 7e5edaff-8597-4f8e-9104-047fff790e99
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
end;

# ╔═╡ 3d1f4ba0-6eec-4d3b-941f-9ad981b1773e
md"## Iterations for finding the α which minimizes the energy expectation value."

# ╔═╡ bf659597-2c99-4af6-8bfd-ba9adbba3a3c
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
end;

# ╔═╡ d5f7f53c-294a-4c46-a3ed-12a8f3a26724
md"# Running the simulation."

# ╔═╡ cee4be07-8a7c-4ebb-a3a3-16ba865bce64
energy,alpha,variance = iterate_alpha();

# ╔═╡ 575bb7ee-49bc-4d38-b473-086e32684261
md" # Plots"

# ╔═╡ 1411e1d3-86e9-4d7b-af64-133aee70636c
plot(energy,
	label = false,
	title = "Energy expenctation Plot",
	linecolor = :blue,
	ylabel = "<E>",
	xlabel = "timestep"
)

# ╔═╡ 1c094b3f-b76f-4b8c-8aa9-efc8967352fc
plot(variance,
	title = "Variance in energy",
	label = false,
	linecolor = :blue,
	ylabel = "Var(E)",
	xlabel = "timestep"
)

# ╔═╡ f0ad854a-5d32-42f5-81e7-9e850081a64e
plot(alpha,
	title = "Variation Parameter",
	label = false,
	linecolor = :blue,
	ylabel = "α",
	xlabel = "timestep"
)

# ╔═╡ a1e1fc43-87d3-4f1e-9d3a-1549aa115aef
scatter(alpha,
	energy,
	title = "Energy vs α",
	xlabel = "α",
	ylabel = "<E>",
	label = false,
	markercolor = :blue
)

# ╔═╡ Cell order:
# ╠═028bce7d-e731-47e9-9468-84185721f0ab
# ╠═c912668e-a13f-11eb-0230-8be9b4bd6d15
# ╠═f04b3d61-d68c-4349-ab98-7136acc75447
# ╟─aa7fc637-2e9d-4c14-a284-50294ba7cbe5
# ╠═89e3a6ae-4f5c-4882-81fd-215379d609dd
# ╠═1197b58f-7614-4e7c-adeb-96af732b7916
# ╠═db71e022-2374-45bc-98af-3fd5e7256784
# ╟─c59f8de5-73db-495f-9a7a-a4b9bf0b58bc
# ╠═7e5edaff-8597-4f8e-9104-047fff790e99
# ╟─3d1f4ba0-6eec-4d3b-941f-9ad981b1773e
# ╠═bf659597-2c99-4af6-8bfd-ba9adbba3a3c
# ╟─d5f7f53c-294a-4c46-a3ed-12a8f3a26724
# ╠═cee4be07-8a7c-4ebb-a3a3-16ba865bce64
# ╟─575bb7ee-49bc-4d38-b473-086e32684261
# ╠═1411e1d3-86e9-4d7b-af64-133aee70636c
# ╠═1c094b3f-b76f-4b8c-8aa9-efc8967352fc
# ╠═f0ad854a-5d32-42f5-81e7-9e850081a64e
# ╠═a1e1fc43-87d3-4f1e-9d3a-1549aa115aef
