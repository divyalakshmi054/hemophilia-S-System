# load the includes -
include("Include.jl")

function objective_function(parameter_guess_array,model,exp_df)
    
    # parameter update -
    model.α = parameter_guess_array[1:10]

    # update G -
    G = model.G
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
   G[TFPI_idx,1] = -1*parameter_guess_array[11];

   # r2 -
   G[AP_idx,2] = parameter_guess_array[12];
   G[HEM_idx,2] = parameter_guess_array[13];

   # r3 -
   # G[HEM_idx,3] = 0.2;

   # r4 -
   G[AP_idx,4] = parameter_guess_array[14];
  
   # r5 -
   G[AP_idx,5] = parameter_guess_array[15];
   G[FVIIa_idx,5] = parameter_guess_array[16];
   G[HEM_idx,5] = parameter_guess_array[17];

   # r6 -
   G[FXa_idx,6] = parameter_guess_array[18];
   G[FVa_idx,6] = parameter_guess_array[19];
   G[HEM_idx,6] = parameter_guess_array[20];
  
   # r9 -
   G[AT_idx,9] = parameter_guess_array[21];

   # r10 -
   G[FIIa_idx,10] = parameter_guess_array[22];

   # put it back -
   model.G = G;

   # solve -
   (T,U) = evaluate(model,tspan=(0.0,120.0))

   # compute simulation error -
   t = exp_df[:,:Time]

   # simulated values -
   FIIa_itp = LinearInterpolation(T,U[:,9])
   FIIa_sim = FIIa_itp[t]

   # get experimental data -
   FIIa_exp = exp_df[:,:FIIa]

   # compute errors -
   error_vector = (FIIa_exp .- FIIa_sim)
   error_term_array = transpose(error_vector)*error_vector

   # total error -
   error_total = error_term_array

   # return -
   return error_total
end

function learn_routine(i,model,visit_df,exp_df; pₒ::Union{Nothing,Array{Float64,1}} = nothing)

    SF = 1e9

    # setup static -
    sfa = model.static_factors_array
    sfa[1] = visit_df[i,:TFPI]       # 1 TFPI
    sfa[2] = visit_df[i,:AT]         # 2 AT
    sfa[3] = (5e-12) * SF            # 3 TF
    sfa[6] = 0.005                   # 6 TRAUMA
    sfa[8] = 0.00075                 # 8 HEM

    # grab the multiplier from the data -
    ℳ = model.number_of_dynamic_states
    xₒ = zeros(ℳ)
    xₒ[1] = visit_df[i, :II]         # 1 FII 
    xₒ[2] = visit_df[i, :VII]        # 2 FVII 
    xₒ[3] = visit_df[i, :V]          # 3 FV
    xₒ[4] = visit_df[i, :X]          # 4 FX
    # xₒ[5] = visit_df[i, :VIII]     # 5 FVIII
    # xₒ[6] = visit_df[i, :IX]       # 6 FIX
    xₒ[7] = visit_df[i, :XI]         # 7 FXI
    xₒ[8] = visit_df[i, :XII]        # 8 FXII 
    xₒ[9] = (1e-14)*SF               # 9 FIIa
    xₒ[10] = 25.0                    # 10 FVIIa --- CHANGE THIS FOR DIFFERENT CONDITIONS
    xₒ[19] = visit_df[i, :PLT]       # 19 PL
    model.initial_condition_array = xₒ

    # set up obj fn -
    OF(p) = objective_function(p,model,exp_df)

    parameter_guess_array = 
    [0.8    0.01    10.0; # 1
    1.0     0.01    10.0; # 2
    1.0     0.01    10.0; # 3
    1.0     0.01    10.0; # 4
    0.95    0.01    10.0; # 5
    1.0     0.01    10.0; # 6
    1.0     0.01    10.0; # 7
    1.0     0.01    10.0; # 8
    0.2     0.01    10.0; # 9
    1.0     0.01    10.0; # 10
    0.65    0.01    10.0; # 11
    0.01    0.001   10.0; # 12
    0.9     0.01    10.0; # 13
    0.25    0.01    10.0; # 14
    0.2     0.01    10.0; # 15
    0.15    0.01    10.0; # 16
    0.5     0.01    10.0; # 17
    0.95    0.01    10.0; # 18
    0.9     0.01    10.0; # 19
    0.3     0.01    10.0; # 20
    0.045   0.01    10.0; # 21
    0.01    0.005   10.0  # 22
    ];

   # if (isnothing(pₒ) == true)
   #     P = length(parameter_guess_array[:,1])
   #     σ = 0.1 # move up to 10%
   #     pₒ = parameter_guess_array[:,1].*(1 .+ σ*rand(-1:1,P))
   # end

    inner_optimizer = NelderMead() 
    options = Optim.Options(time_limit = 600, show_trace = true, show_every = 10, iterations = 100)
    result = optimize(OF, parameter_guess_array[:,2],parameter_guess_array[:,3],parameter_guess_array[:,1],Fminbox(inner_optimizer),options)

    return result;
end

# load exp data & model -
exp_df = CSV.read(joinpath(_PATH_TO_DATA,"thrombin-hemophilia-fVIIa-25-nM-100-fII.csv"),DataFrame)
model_file = joinpath(_PATH_TO_MODEL,"Feedback.bst")

# build the model -
model = build(model_file)

# load the training data -
data_file = joinpath(_PATH_TO_DATA,"Training-Composition-Transformed-w-Labels.csv")
training_df = CSV.read(data_file,DataFrame)

# which visit?
visit = 4;

visit_df = filter(:Visit => x->(x==visit), training_df) 

# size of training set -
(R,C) = size(visit_df)

p_previous = nothing;
number_of_parameters = 22;
ensemble_archive = zeros(number_of_parameters+1)

for i ∈ 3:R
    result = learn_routine(i,model,visit_df,exp_df; pₒ=nothing)
    p_best = Optim.minimizer(result)

    # error -
    fitness = Optim.minimum(result)
    global ensemble_archive[1] = fitness
    
    # cache the parameters -
    for k ∈ 1:number_of_parameters
        global ensemble_archive[k+1] = p_best[k]
    end

    # run with best parameters -
    model.α = p_best[1:10]

    # update G -
    G = model.G
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
   G[TFPI_idx,1] = -1*p_best[11];

   # r2 -
   G[AP_idx,2] = p_best[12];
   G[HEM_idx,2] = p_best[13];

   # r3 -
   # G[HEM_idx,3] = 0.2;

   # r4 -
   G[AP_idx,4] = p_best[14];
  
   # r5 -
   G[AP_idx,5] = p_best[15];
   G[FVIIa_idx,5] = p_best[16];
   G[HEM_idx,5] = p_best[17];

   # r6 -
   G[FXa_idx,6] = p_best[18];
   G[FVa_idx,6] = p_best[19];
   G[HEM_idx,6] = p_best[20];
  
   # r9 -
   G[AT_idx,9] = p_best[21];

   # r10 -
   G[FIIa_idx,10] = p_best[22];

   # put it back -
   model.G = G;

   # solve -
   (T,U) = evaluate(model,tspan=(0.0,120.0))

   # dump p_best to disk -
   CSV.write(joinpath(_PATH_TO_ACTUAL_ENSEMBLE,"PSET-Actual-25-100-P$(i).csv"),Tables.table(ensemble_archive),header = ["parameters"]) # 50 rFVIIa 100% fII

   # dump SIM to disk -
   data = [T U]

    path_to_sim_data = joinpath(_PATH_TO_ACTUAL_ENSEMBLE, "SIM-Hemophilia-25-100-$(i).csv")
    CSV.write(path_to_sim_data, Tables.table(hcat(data),header=vcat("Time",model.list_of_dynamic_species)))

    global p_previous = nothing;
end
