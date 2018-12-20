using CSV, DataFrames, StatsBase, Statistics, Plots, LinearAlgebra, Grep

df_corr_file = "/Users/matt/OneDrive/UCSF/biostats/df_corr.csv"
corr = CSV.File(df_corr_file, delim=',', header=1, allowmissing=:auto) |> DataFrame
corr_array = convert(Array, corr[2:size(corr, 2)])
idx = corr[1]

N = 113725 #how many samples
S = 1000 #how many permutations
perm = Array{Float64}(undef, S, N)
for p in 1:S
    perm[p,:] = sample(corr_f2, N)
end

mean!(Array{Float64}(undef, S, 1), perm)

expression_file = "/Users/matt/OneDrive/UCSF/biostats/average_expression.csv"
exp_df = CSV.File(expression_file, delim=',', header=1, allowmissing=:auto) |> DataFrame
exp_array = convert(Array, exp_df)

n_row = size(exp_array)[1]
n_col = size(exp_array)[2]

r = UpperTriangular(corspearman(transpose(convert(Array{Float64,2}, exp_array[:,2:n_col]))))

pr_array_flat = filter(x -> (x != 0.0 && x != 1.0), pr_array)
h = fit(Histogram, convert(Array{Float32}, pr_array_flat))
plot(h)
iteractome_file = "/Users/matt/OneDrive/UCSF/biostats/DataS1_interactome.tsv"
it_df = CSV.File(iteractome_file, delim='\t', header=0, allowmissing=:auto) |> DataFrame
it_array = convert(Array, it_df)
idx = exp_array[:,1]
#create new r_it matrix
nodes = convert(Array{Int32,2}, it_array[:,1:2])
grep("100847079|7318", nodes)
using Match

r_it_idx = []
for potential_nodes in idx[1:460]
    multival_nodes = split(potential_nodes, "|")
    for n in multival_nodes
        if !Match.ismatch(n, unique(nodes))
            append!(r_it_idx, multival_nodes)
        end
    end
    #println(split(potential_nodes, "|"))
end
idx_fixed = convert(Array{String}, collect(skipmissing(idx)))
pr_array = []
for i in 1:size(nodes)[1]
    aa , bb = nodes[i,:]
    aa_loc = findall(x->x==string(aa), collect(idx_fixed))
    bb_loc = findall(x->x==string(bb), collect(idx_fixed))
    if (!isempty(aa_loc) && !isempty(bb_loc))
        if !(aa_loc == bb_loc)
            if isequal(vec(r[aa_loc, bb_loc])[1], 0)
                pr = r[bb_loc, aa_loc]
                append!(pr_array, pr)
            end
            pr = r[aa_loc, bb_loc]
            append!(pr_array, pr)
        end
    end
end

corspearman([24,1420,342],[134,1353,99425])
plot()
pr_array_flat = filter(x -> (x != 0.0), pr_array)
h = fit(Histogram, convert(Array{Float32}, pr_array_flat))
plot(h)
mean(pr_array_flat)

pr_abs_array_flat = []
for pr in convert(Array{Float32}, pr_array_flat)
    pr = abs(pr)
    append!(pr_abs_array_flat, pr)
end

mean(pr_abs_array_flat) #0.3685038f0

#simulation
sim_mean = Array{Float64}(undef, 1000)
for ii in 1:1000
    sim = filter(x -> (x != 0.0 && x != 200.0), sample(r_copy,(98915*2)))
    sim_abs = Array{Float64}(undef, size(sim)[1])
    for i in 1:size(sim)[1]
        npr = abs(sim[i])
        sim_abs[i] = npr
    end
    sim_mean[ii] = mean(sim_abs)
end
has = fit(Histogram, convert(Array{Float32}, sim_mean))
theme(:default)
plot(h, background_color = "black", dpi = 200)
savefig( "/Users/matt/OneDrive/UCSF/biostats/mean_null.png")
