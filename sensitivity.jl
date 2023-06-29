# load the include -
include("Include.jl")

# build performance function -
function performance(κ, model::BSTModel, visit_df::DataFrame, i::Int64)

    # main simulation -
    SF = 1e9

    # setup static -
    sfa = model.static_factors_array
    sfa[1] = visit_df[i,:TFPI]       # 1 TFPI
    sfa[2] = visit_df[i,:AT]         # 2 AT
    sfa[3] = (5e-12) * SF            # 3 TF
    sfa[6] = 0.005                   # 6 TRAUMA
    sfa[8] = 0.00075                 # 8 HEM
     
    # setup dynamic -
    # grab the multiplier from the data -
    ℳ = model.number_of_dynamic_states
    xₒ = zeros(ℳ)
    xₒ[1] = visit_df[i, :II]         # 1 FII 
    xₒ[2] = visit_df[i, :VII]        # 2 FVII 
    xₒ[3] = visit_df[i, :V]          # 3 FV
    xₒ[4] = visit_df[i, :X]          # 4 FX
    # xₒ[5] = visit_df[i, :VIII]       # 5 FVIII
    # xₒ[6] = visit_df[i, :IX]         # 6 FIX
    xₒ[7] = visit_df[i, :XI]           # 7 FXI
    xₒ[8] = visit_df[i, :XII]          # 8 FXII 
    xₒ[9] = (1e-14)*SF                 # 9 FIIa
    # xₒ[10] = 50.0                      # 10 FVIIa
    xₒ[19] = visit_df[i, :PLT]         # 19 PL
    model.initial_condition_array = xₒ
    
    #get the parameters -
    model.α = κ[1:10]

    # set G values -
    G = model.G;

    AT_idx = findfirst(x->x=="AT",model.total_species_list)
    TFPI_idx = findfirst(x->x=="TFPI",model.total_species_list)
    FIIa_idx = findfirst(x->x=="FIIa",model.total_species_list)
    AP_idx = findfirst(x->x=="AP",model.total_species_list)
    PL_idx = findfirst(x->x=="PL",model.total_species_list)
    FVIIa_idx = findfirst(x->x=="FVIIa",model.total_species_list)
    FXa_idx = findfirst(x->x=="FXa",model.total_species_list)
    FVa_idx = findfirst(x->x=="FVa",model.total_species_list)
    HEM_idx = findfirst(x->x=="HEM",model.total_species_list)
 
    # r1 -
   G[TFPI_idx,1] = -1*κ[11];

   # r2 -
   G[AP_idx,2] = κ[12];
   G[HEM_idx,2] = κ[13];
   
   # r4 -
   G[AP_idx,4] = κ[14];
  
   # r5 -
   G[AP_idx,5] = κ[15];
   G[FVIIa_idx,5] = κ[16];
   G[HEM_idx,5] = κ[17];

   # r6 -
   G[FXa_idx,6] = κ[18];
   G[FVa_idx,6] = κ[19];
   G[HEM_idx,6] = κ[20];
  
   # r9 -
   G[AT_idx,9] = κ[21];

   # r10 -
   G[FIIa_idx,10] = κ[22];

    # run the model -
    global (T,U) = evaluate(model,tspan=(0.0,120.0))

    # test -
    return integrate(T,U[:,9])    # AUC
end

# build the model structure -
path_to_model_file = joinpath(pwd(), "model", "Feedback.bst")

# build the default model structure -
model = build(path_to_model_file)

# load the training data -
_PATH_TO_DATA = joinpath(pwd(),"data")
path_to_training_data = joinpath(_PATH_TO_DATA, "Training-Composition-Transformed-w-Labels.csv")
training_df = CSV.read(path_to_training_data, DataFrame)

# which visit?
visit = 4;

# let's filter visit 4s since we look to train using that visit
visit_df = filter(:Visit => x->(x==visit), training_df) 

# size of training set -
(R,C) = size(visit_df)

# a = [0.8, 1.0, 1.0, 1.0, 0.95, 1.0, 1.0, 1.0, 0.2, 1.0]

#update G -
# G = model.G
# g = [0.65, 0.01, 0.9, 0.25, 0.2, 0.15, 0.5, 0.95, 0.9, 0.3, 0.045, 0.01] # look at sample_ensemble.jl for specific G values

# setup sensitivity analysis -
# load a pset -
pset_filename = "PSET-Actual-P6.csv"
pset_df = CSV.read(joinpath(_PATH_TO_ACTUAL_ENSEMBLE, pset_filename), DataFrame)

κ = pset_df[2:end,:parameters]

# create lower and upper bound array -
NP = length(pset_df[!,:parameters])
L = zeros(NP-1)
U = zeros(NP-1)
for pᵢ ∈ 1:(NP - 1)
    L[pᵢ] = 0.1*pset_df[pᵢ + 1,:parameters]
    U[pᵢ] = 10.0*pset_df[pᵢ + 1,:parameters]
end

# samples = 1000;
# bootreps = 100;

# initialize -
# sampler = SobolSample()
    
# generate a sampler -
# (A,B) = QuasiMonteCarlo.generate_design_matrices(samples,L,U,sampler)

# setup call to Sobol method -
F(κ) =  performance(κ, model, visit_df, 6)
m = gsa(F,Morris(num_trajectory=10000),[[L[i],U[i]] for i in 1:(NP-1)]);

# dump -
# results_array = hcat(m.ST, m.S1, m.ST_Conf_Int, m.S1_Conf_Int)
means = m.means;
means_star = m.means_star;
variances = m.variances;
results_array = vcat(means,means_star,variances)

# write -
# CSV.write(joinpath(pwd(),"sobol","Sensitivity-Sobol-$(samples)-second-order-boot-$(bootreps).csv"), Tables.table(m.S2))
# CSV.write(joinpath(pwd(),"sobol","Sensitivity-Sobol-$(samples)-second-order-CI-boot-$(bootreps).csv"), Tables.table(m.S2_Conf_Int))
# CSV.write(joinpath(pwd(),"sobol","Sensitivity-Sobol-$(samples)-boot-$(bootreps).csv"), Tables.table(results_array), header = vcat("Total_order", "First_order", "Total_order_CI", "First_order_CI"))
# dump sensitivity data to disk -
CSV.write(joinpath(pwd(),"data","Sensitivity-Morris-test-10000-H.csv"), Tables.table(transpose(results_array)), header = vcat("mean", "mean_star","variance"))