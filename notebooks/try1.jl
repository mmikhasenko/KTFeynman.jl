### A Pluto.jl notebook ###
# v0.17.2

using Markdown
using InteractiveUtils

# ╔═╡ ac94273e-6f9a-11ec-07c9-ed5690c52aa8
begin
	cd(joinpath(@__DIR__, ".."))
	# 
	using Pkg
	Pkg.activate(".")
	Pkg.instantiate()
end

# ╔═╡ 78c3ad9f-58a2-4c4f-b0bd-2735376fe07c
begin
	using SymbolicTensors
	using SymPy
end

# ╔═╡ cc98cf46-721d-409e-9341-792f0849dc36
begin
	lorenz = TensorIndexType("lorenz","L",4)
	# lorenz.set_metric([1,-1,-1,-1])
end

# ╔═╡ aba8a13b-b096-40b0-be9a-66c190e66947
lorenz

# ╔═╡ 1caab6ba-87db-4f4c-a959-d10f4d3f232f
typeof([SymPy.Sym("lorenz")])

# ╔═╡ dc8d9780-f16e-4dff-8f03-f562efabb796
begin
	@indices lorenz μ ν σ ρ η i
	μ,ν,σ,ρ,η
end

# ╔═╡ 668f202f-154c-4fe8-bf8e-3c46cd86bbe7
begin
	δ = lorenz.delta
	g = lorenz.metric
	ϵ = lorenz.epsilon
end

# ╔═╡ 9fd00605-f6c8-4986-89bc-1e2ed8209917
lorenz

# ╔═╡ 79dd5294-016f-4b64-b57b-2ab731844f07
lorenz.epsilon

# ╔═╡ 7f08dd0f-9099-4e9f-a0b8-710eda39625f


# ╔═╡ d2d1f03d-562a-49dd-9e69-f660cc1cf199
@heads lorenz  pX p3

# ╔═╡ 8b1f9660-8df1-4814-b3b8-866b41436b10
begin
	# @macroexpand(@heads [lorenz] pX p3 p1 p2)
	pX = TensorHead("pX", [lorenz])
	p3 = TensorHead("p3", [lorenz])
	p1 = TensorHead("p1", [lorenz])
	p2 = TensorHead("p2", [lorenz])
end

# ╔═╡ 35b4c06e-4e81-410b-8f12-08b5d216693e
let
	lorenz = TensorIndexType("lorenz","L",4)
	@heads [lorenz] pX p3 p1 p2
end

# ╔═╡ 691504b0-e169-4a95-9358-26a03a6e6d13
let
	lorenz = TensorIndexType("lorenz","L",4)
	@indices lorenz i j
	p = TensorHead("p", [lorenz])
	p(i)*p(-i)
end

# ╔═╡ 70ccf0f2-1f1f-416c-a61e-152ea85d6258
pX(i)*pX(-i)

# ╔═╡ 67daf4a4-b47f-4a82-96ef-3a47d0b4cee9
lorenz.epsilon(μ,ν,σ,ρ)

# ╔═╡ 564f76f0-8381-45df-99f2-c198a4652270
let
	lorenz = TensorIndexType("lorenz","L",4)
	@indices lorenz μ ν σ ρ
	lorenz.epsilon(μ,ν,σ,ρ)
end

# ╔═╡ Cell order:
# ╠═ac94273e-6f9a-11ec-07c9-ed5690c52aa8
# ╠═78c3ad9f-58a2-4c4f-b0bd-2735376fe07c
# ╠═cc98cf46-721d-409e-9341-792f0849dc36
# ╠═aba8a13b-b096-40b0-be9a-66c190e66947
# ╠═1caab6ba-87db-4f4c-a959-d10f4d3f232f
# ╠═dc8d9780-f16e-4dff-8f03-f562efabb796
# ╠═668f202f-154c-4fe8-bf8e-3c46cd86bbe7
# ╠═9fd00605-f6c8-4986-89bc-1e2ed8209917
# ╠═79dd5294-016f-4b64-b57b-2ab731844f07
# ╠═7f08dd0f-9099-4e9f-a0b8-710eda39625f
# ╠═d2d1f03d-562a-49dd-9e69-f660cc1cf199
# ╠═8b1f9660-8df1-4814-b3b8-866b41436b10
# ╠═35b4c06e-4e81-410b-8f12-08b5d216693e
# ╠═691504b0-e169-4a95-9358-26a03a6e6d13
# ╠═70ccf0f2-1f1f-416c-a61e-152ea85d6258
# ╠═67daf4a4-b47f-4a82-96ef-3a47d0b4cee9
# ╠═564f76f0-8381-45df-99f2-c198a4652270
