### A Pluto.jl notebook ###
# v0.17.2

using Markdown
using InteractiveUtils

# ╔═╡ 4e7dd940-751b-11ec-1132-377cc9370086
begin
	using SymPy
	import PyCall
	#
	PyCall.pyimport_conda("sympy.physics.quantum.spin", "sympy")
	PyCall.pyimport_conda("sympy.physics.wigner",       "sympy")
	
	import_from(sympy.physics.quantum.spin, (:WignerD,), typ=:Any)
	import_from(sympy.physics.wigner)
	import_from(sympy.physics.quantum.spin)
end

# ╔═╡ 17947eb8-1889-4a2f-9f4f-159b98ef1e35
md"""
# Crossing matrix of $X \to 3\pi$
"""

# ╔═╡ 7fe8a151-5407-465c-9f02-fd58e825b8d7
wignerd(j,m1,m2,θ) = WignerD(j,m1,m2, 0, θ, 0).doit()

# ╔═╡ 2164c234-eaa7-4413-a0c4-ffe11555b521
begin
	mX, m1, m2, m3 = @vars m_X m_1 m_2 m_3 positive=true
	s, t = @vars s t positive=true
	mX, m1, m2, m3, s, t
end

# ╔═╡ fbaf7a48-9f0d-4922-8e9a-2d675f06eb33
λXs, λs, λXt, λt = @vars lambda_Xs lambda_s lambda_Xt lambda_t positive=true

# ╔═╡ 7a82b4fd-5872-453d-a8bd-36f850633887
begin
	# crossing angle
	ω, = @vars omega
	n, ϕ = @vars n phi
	# 
	cosω = n/sqrt(λXt*λXs)
	sinω = 2*sqrt(phi)*mX/sqrt(λXt*λXs)
	# 
	cos(ω)=>cosω, sin(ω)=>sinω
end

# ╔═╡ 215261eb-412a-45ba-aaf1-9d8f53c9373e
begin
	K(λ, Yj, Kallen) = (2sqrt(ϕ))^λ * sqrt(Kallen)^Yj
	Ks(λ, Yj) = K(λ, Yj, λXs)
	Kt(λ, Yj) = K(λ, Yj, λXt)
end

# ╔═╡ 8f8f7af3-0fa8-443b-a1df-0e44a4451264
function paritywignerd(J,λ,λ′,fac)
    wd1 = wignerd(J, λ,λ′, ω)
    wd2 = wignerd(J, λ, -λ′, ω)
    return wd1 + fac * (-1)^λ′ * wd2
end

# ╔═╡ 3643c393-adfb-4bfc-81a1-3bc55f9eb18c
ℂij(J,λ,λ′,Yⱼs, Yⱼt,fac) = 
	paritywignerd(J,λ,λ′,fac) / Ks(λ, Yⱼs) * Kt(λ′, Yⱼt)

# ╔═╡ 0a1a9546-0b65-41f5-b45c-fe3deef0bfa9
Yj(j,YX) = abs(j - YX) - j

# ╔═╡ 30541e5c-4480-4831-a115-ed920dc5d820
const ηπ = -1

# ╔═╡ 4194d764-d62b-4d5e-8775-5d993056745d
YX(J,ηX) = J - (1 + ηX) / Sym(2)

# ╔═╡ c2d0691d-3d8d-4cb8-977f-600af3725ed6
function ℂ(J;ηX,js,jt)
    Yₓ = YX(J,ηX)
    Yⱼs, Yⱼt = Yj(js,Yₓ), Yj(jt,Yₓ)
    ηs = ηπ^3 * ηX
	λmin = ηX==1 ? 1 : 0
	# 
    m = [ℂij(J, λ, λ′, Yⱼs, Yⱼt, ηs) for λ in λmin:js, λ′ in λmin:jt]
    return m
end

# ╔═╡ e551a603-dce1-403f-8ca7-678ff54fa187
replacesincos(e) = simplify(
	subs(
		sympy.expand_trig(e),
		cos(ω) => cosω, sin(ω) => sinω))

# ╔═╡ 65e5acfd-c3f4-455f-abe5-0adea9186d07
cosmetics(e) = simplify(subs(e, λXs*λXt=>n^2 + 4mX^2*ϕ))

# ╔═╡ a509c206-8a86-4c68-9c13-e9ed1259e51d
md"""
## Examples
"""

# ╔═╡ 3046f5d2-64fe-4c1a-9de4-6d87535e5bf0
md"""
### $J^{PC}=1^{--}\,\rho \pi$ 
"""

# ╔═╡ bf42f57e-d81e-4d78-992e-09f6a6eb3b6b
ℂϕ = replacesincos.(ℂ(1; ηX=1, js=1, jt=1))

# ╔═╡ 03096361-5378-4618-a1c4-2630092de4c9
md"""
### $J^{PC}=1^{++}\,f_2 \pi$
"""

# ╔═╡ ceaef4b3-d410-4054-8bfc-c73208073c0a
ℂa1 = replacesincos.(ℂ(1; ηX=-1, js=1, jt=1))

# ╔═╡ dcfe1ec0-7d82-4e4e-828d-7d30e25aadba
md"""
### $J^{PC}=2^{-+}\,f_2 \pi$
"""

# ╔═╡ 7139fd4b-3e6c-4c2b-92ac-89188f6e6214
ℂπ2 = cosmetics.(replacesincos.(ℂ(2; ηX=-1, js=2, jt=2) .* λXt^2))

# ╔═╡ d79f9672-cb49-4021-820b-f767d62263c6
md"""
### $J^{PC}=3^{++} \rho_3 \pi$
"""

# ╔═╡ eb90abab-603e-4a59-9af4-40f8db08f249
ℂa3 = cosmetics.(replacesincos.(ℂ(3; ηX=-1, js=3, jt=3) .* λXt^3))

# ╔═╡ 92ee9179-a2c0-4ba9-9555-2e9583d72156
md"""
## Cases to be checked: j<J
"""

# ╔═╡ eb2bacaa-5f88-41cc-a5b3-2667beda18b2
md"""
### $J^{PC}=1^{++}\,\sigma \pi \leftrightarrow \rho \pi$
"""

# ╔═╡ 87febb4c-7553-4cf7-80ee-d822e8938aa9
replacesincos.(ℂ(1; ηX=-1, js=0, jt=0))

# ╔═╡ a03c5169-f81f-4fca-b7ed-f415094693d9
replacesincos.(ℂ(1; ηX=-1, js=1, jt=0))

# ╔═╡ 6e0a0576-54fa-494d-8d95-2714466be9f1
replacesincos.(ℂ(1; ηX=-1, js=0, jt=1))

# ╔═╡ 69bd1da4-0479-4034-8fb4-8ac8048d0f64
md"""
### $J^{PC}=2^{-+}\,\rho \pi$
"""

# ╔═╡ a52f082e-7a33-46e7-ba97-36d09785c4e1
replacesincos.(ℂ(2; ηX=-1, js=1, jt=1))

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"

[compat]
PyCall = "~1.93.0"
SymPy = "~1.1.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.0"
manifest_format = "2.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "926870acb6cbcf029396f2f2de030282b6bc1941"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.11.4"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[deps.CommonEq]]
git-tree-sha1 = "d1beba82ceee6dc0fce8cb6b80bf600bbde66381"
uuid = "3709ef60-1bee-4518-9f2f-acd86f176c50"
version = "0.2.0"

[[deps.CommonSolve]]
git-tree-sha1 = "68a0743f578349ada8bc911a5cbd5a2ef6ed6d1f"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.0"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "44c37b4636bc54afac5c574d2d02b625349d6582"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.41.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Conda]]
deps = ["Downloads", "JSON", "VersionParsing"]
git-tree-sha1 = "6cdc8832ba11c7695f494c9d9a1c31e90959ce0f"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.6.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "a7254c0acd8e62f1ac75ad24d5db43f5f19f3c65"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.2"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "22df5b96feef82434b07327e2d3c770a9b21e023"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a8f4f279b6fa3c3c4f1adadd78a621b13a506bce"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.9"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "e5718a00af0ab9756305a0392832c8952c7426c1"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.6"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "92f91ba9e5941fc781fecf5494ac1da87bdac775"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "2cf929d64681236a2e074ffafb8d568733d2e6af"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.PyCall]]
deps = ["Conda", "Dates", "Libdl", "LinearAlgebra", "MacroTools", "Serialization", "VersionParsing"]
git-tree-sha1 = "71fd4022ecd0c6d20180e23ff1b3e05a143959c2"
uuid = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
version = "1.93.0"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "e08890d19787ec25029113e88c34ec20cac1c91e"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.0.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.SymPy]]
deps = ["CommonEq", "CommonSolve", "Latexify", "LinearAlgebra", "Markdown", "PyCall", "RecipesBase", "SpecialFunctions"]
git-tree-sha1 = "8f8d948ed59ae681551d184b93a256d0d5dd4eae"
uuid = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
version = "1.1.2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.VersionParsing]]
git-tree-sha1 = "e575cf85535c7c3292b4d89d89cc29e8c3098e47"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.2.1"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╟─17947eb8-1889-4a2f-9f4f-159b98ef1e35
# ╠═4e7dd940-751b-11ec-1132-377cc9370086
# ╠═7fe8a151-5407-465c-9f02-fd58e825b8d7
# ╠═2164c234-eaa7-4413-a0c4-ffe11555b521
# ╠═fbaf7a48-9f0d-4922-8e9a-2d675f06eb33
# ╠═7a82b4fd-5872-453d-a8bd-36f850633887
# ╠═215261eb-412a-45ba-aaf1-9d8f53c9373e
# ╠═8f8f7af3-0fa8-443b-a1df-0e44a4451264
# ╠═3643c393-adfb-4bfc-81a1-3bc55f9eb18c
# ╠═0a1a9546-0b65-41f5-b45c-fe3deef0bfa9
# ╠═30541e5c-4480-4831-a115-ed920dc5d820
# ╠═4194d764-d62b-4d5e-8775-5d993056745d
# ╠═c2d0691d-3d8d-4cb8-977f-600af3725ed6
# ╠═e551a603-dce1-403f-8ca7-678ff54fa187
# ╠═65e5acfd-c3f4-455f-abe5-0adea9186d07
# ╠═a509c206-8a86-4c68-9c13-e9ed1259e51d
# ╟─3046f5d2-64fe-4c1a-9de4-6d87535e5bf0
# ╠═bf42f57e-d81e-4d78-992e-09f6a6eb3b6b
# ╟─03096361-5378-4618-a1c4-2630092de4c9
# ╠═ceaef4b3-d410-4054-8bfc-c73208073c0a
# ╟─dcfe1ec0-7d82-4e4e-828d-7d30e25aadba
# ╠═7139fd4b-3e6c-4c2b-92ac-89188f6e6214
# ╟─d79f9672-cb49-4021-820b-f767d62263c6
# ╠═eb90abab-603e-4a59-9af4-40f8db08f249
# ╟─92ee9179-a2c0-4ba9-9555-2e9583d72156
# ╟─eb2bacaa-5f88-41cc-a5b3-2667beda18b2
# ╠═87febb4c-7553-4cf7-80ee-d822e8938aa9
# ╠═a03c5169-f81f-4fca-b7ed-f415094693d9
# ╠═6e0a0576-54fa-494d-8d95-2714466be9f1
# ╟─69bd1da4-0479-4034-8fb4-8ac8048d0f64
# ╠═a52f082e-7a33-46e7-ba97-36d09785c4e1
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
