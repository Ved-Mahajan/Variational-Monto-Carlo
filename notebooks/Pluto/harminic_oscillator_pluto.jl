### A Pluto.jl notebook ###
# v0.14.3

using Markdown
using InteractiveUtils

# ╔═╡ c4f5040e-7a5c-4635-9402-9008b0f41a71
using Pkg; Pkg.add("Plots")

# ╔═╡ 96b47390-a131-11eb-0f4c-7782f9731bd1
using Plots

# ╔═╡ 763f61ce-189f-4e6b-b422-0a79defd2ca4
using CSV,DataFrames

# ╔═╡ 02fcc95b-2349-4867-a0f3-716a4904b922
md"""
# Define all the functions
"""

# ╔═╡ bb81b5b9-5c35-41e2-b23d-f76e0a1486a9
function ansatz(x,α)
	return exp(-α*x^2)
end;

# ╔═╡ 9f41b093-123d-44f3-84db-ec987cf62628
function probability_density(x,α)
	return (ansatz(x,α)^2)/(sqrt(pi/α))
end;

# ╔═╡ 4b4252e5-4694-4bdd-ae9a-166bd66f41a1
function E_local(x,α)
	return α + (x^2)*(1/2 - 2*α^2)
end;

# ╔═╡ d8ef5c0c-736b-4746-84f2-f5dde310fec8
md" # Function for running the metropolis algorithm."

# ╔═╡ 281d664f-2d2b-4b98-960d-bbc1af5cf1a9
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
end;

# ╔═╡ f74a56fe-595e-4eed-9025-d9ecbc760f45
md"## Iterations for finding the α which minimizes the energy expectation value."

# ╔═╡ bc50f526-55fe-4059-9b12-8c568aea3009
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
end;

# ╔═╡ 9ec4cea0-bcc3-45f2-97e5-33d116da3fb1
md"# Running the simulation."

# ╔═╡ a4b16f08-0cc9-47fd-8836-188efe891b84
energy,alpha,variance,E_analytical,var_analytical = iterate_alpha();

# ╔═╡ 9b375035-df0d-415b-8975-dbb4fcec0381
md" # Plots"

# ╔═╡ 2276a325-9649-4eee-98c9-138dcb457b8f
begin
	plot(energy,
		label = "calculated expectation of energy",
		title = "Energy Plots",
		xlabel = "timestep",
		ylabel = "<E>",
		linecolor = :blue,
	)
	plot!(E_analytical,label = "Analytical Energy")
end

# ╔═╡ 285e1d2b-3e11-4c45-8646-884ce3d06c77
plot(variance,
	title = "Variance in energy",
	label = false,
	ylabel = "Var(E)",
	xlabel = "timestep",
	linecolor = :blue
)

# ╔═╡ 6ecec340-0769-4881-bee2-eedd9fdc4eb4
begin
	scatter(alpha,energy,
		label = "Calculated energy",
		title = "Energy vs α",
		markercolor = :blue,
		xlabel = "α",
		ylabel = "<E>"
	)

	scatter!(alpha,
		E_analytical,
		label = "Calculated energy",
		title = "Energy vs α",
		markercolor = :red
	)
end

# ╔═╡ Cell order:
# ╠═c4f5040e-7a5c-4635-9402-9008b0f41a71
# ╠═96b47390-a131-11eb-0f4c-7782f9731bd1
# ╠═763f61ce-189f-4e6b-b422-0a79defd2ca4
# ╟─02fcc95b-2349-4867-a0f3-716a4904b922
# ╠═bb81b5b9-5c35-41e2-b23d-f76e0a1486a9
# ╠═9f41b093-123d-44f3-84db-ec987cf62628
# ╠═4b4252e5-4694-4bdd-ae9a-166bd66f41a1
# ╟─d8ef5c0c-736b-4746-84f2-f5dde310fec8
# ╠═281d664f-2d2b-4b98-960d-bbc1af5cf1a9
# ╟─f74a56fe-595e-4eed-9025-d9ecbc760f45
# ╠═bc50f526-55fe-4059-9b12-8c568aea3009
# ╟─9ec4cea0-bcc3-45f2-97e5-33d116da3fb1
# ╠═a4b16f08-0cc9-47fd-8836-188efe891b84
# ╟─9b375035-df0d-415b-8975-dbb4fcec0381
# ╠═2276a325-9649-4eee-98c9-138dcb457b8f
# ╠═285e1d2b-3e11-4c45-8646-884ce3d06c77
# ╠═6ecec340-0769-4881-bee2-eedd9fdc4eb4
