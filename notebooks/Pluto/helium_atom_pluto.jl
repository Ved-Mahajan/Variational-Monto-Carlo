### A Pluto.jl notebook ###
# v0.14.3

using Markdown
using InteractiveUtils

# ╔═╡ dba34531-0a3d-4598-a4ba-52cf81203caf
using Pkg; Pkg.add(["Plots","LinearAlgebra"])

# ╔═╡ ed657f60-a5f3-11eb-2974-8714af0f16dc
using LinearAlgebra,Plots

# ╔═╡ 97916012-f5e5-4957-99c3-4f91d7b86b9d
function ansatz(r1,r2,α)
    norm_r1 = norm(r1)
    norm_r2 = norm(r2)
    r12 = norm(r1-r2)
    wf = exp(-2*norm_r1)*exp(-2*norm_r2)*exp(r12/(2*(1+α*r12)))
    return wf
end;

# ╔═╡ 9a7c2034-7b47-4c70-b469-60d0e1bcf159
function prob_density(r1,r2,α)
    return ansatz(r1,r2,α)^2
end;

# ╔═╡ 9f80362c-5e59-4e73-b878-ddba6429974a
function E_local(r1,r2,α)
    norm_r1 = norm(r1)
    norm_r2 = norm(r2)
    r12 = norm(r1-r2)
    dot_product = dot(r1/norm_r1-r2/norm_r2,r1-r2)
    energy = -4+dot_product/(r12*(1+α*r12)^2)-1/(r12*(1+α*r12)^3)-1/(4*(1+α*r12)^4)+1/r12
    return energy
end;

# ╔═╡ 3c677d47-d02c-4328-9e12-7ad2ce540ade
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
end;

# ╔═╡ 4ad41744-94c3-4624-b5b3-f5bf679a7957

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
end;

# ╔═╡ 82e68fc3-f949-43f8-8704-f05e19dec3aa
energy,alpha,variance = iterate_alpha();

# ╔═╡ 6f50f51d-1ec8-4031-af1a-0d6cd69171a0
plot(energy,
	label = false,
	title = "Energy expenctation Plot",
	linecolor = :blue,
	ylabel = "<E>",
	xlabel = "timestep"
)

# ╔═╡ d797e924-293c-496b-8ac8-f06b59b9ad89
plot(variance,
	title = "Variance in energy",
	label = false,
	linecolor = :blue,
	ylabel = "Var(E)",
	xlabel = "timestep"
)

# ╔═╡ 03294208-6397-4602-9157-58ca24154003
plot(alpha,
	title = "Variation Parameter",
	label = false,
	linecolor = :blue,
	ylabel = "α",
	xlabel = "timestep"
)

# ╔═╡ 53b76a30-942b-4191-982a-f06fb38bd52f
scatter(alpha,
	energy,
	title = "Energy vs α",
	xlabel = "α",
	ylabel = "<E>",
	label = false,
	markercolor = :blue
)

# ╔═╡ Cell order:
# ╠═dba34531-0a3d-4598-a4ba-52cf81203caf
# ╠═ed657f60-a5f3-11eb-2974-8714af0f16dc
# ╠═97916012-f5e5-4957-99c3-4f91d7b86b9d
# ╠═9a7c2034-7b47-4c70-b469-60d0e1bcf159
# ╠═9f80362c-5e59-4e73-b878-ddba6429974a
# ╠═3c677d47-d02c-4328-9e12-7ad2ce540ade
# ╠═4ad41744-94c3-4624-b5b3-f5bf679a7957
# ╠═82e68fc3-f949-43f8-8704-f05e19dec3aa
# ╠═6f50f51d-1ec8-4031-af1a-0d6cd69171a0
# ╠═d797e924-293c-496b-8ac8-f06b59b9ad89
# ╠═03294208-6397-4602-9157-58ca24154003
# ╠═53b76a30-942b-4191-982a-f06fb38bd52f
