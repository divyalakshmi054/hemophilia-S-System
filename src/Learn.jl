# include -
include("Simulation.jl")

function learn_optim(index::Int, model::BSTModel, training_df::DataFrame, Y; 
    pₒ::Union{Nothing,Array{Float64,1}} = nothing)

    # main simulation loop -
    SF = 1e9

    # setup static -
    sfa = model.static_factors_array
    sfa[1] = training_df[index,:TFPI]       # 1 TFPI
    sfa[2] = training_df[index,:AT]         # 2 AT
    sfa[3] = (5e-12) * SF               # 3 TF
    sfa[6] = 0.005                      # 6 TRAUMA

    # grab the multiplier from the data -
    ℳ = model.number_of_dynamic_states
    xₒ = zeros(ℳ)
    xₒ[1] = training_df[index, :II]         # 1 FII 
    xₒ[2] = training_df[index, :VII]        # 2 FVII 
    xₒ[3] = training_df[index, :V]          # 3 FV
    xₒ[4] = training_df[index, :X]          # 4 FX
    xₒ[5] = training_df[index, :VIII]       # 5 FVIII
    xₒ[6] = training_df[index, :IX]         # 6 FIX
    xₒ[7] = training_df[index, :XI]         # 7 FXI
    xₒ[8] = training_df[index, :XII]        # 8 FXII 
    xₒ[9] = (1e-14)*SF                      # 9 FIIa
    model.initial_condition_array = xₒ

    # setup initial parameter values and bounds array -
    κ = [
            
            # default: hand fit set -
            0.70    0.001 10.0   ; # 1
            1.0     0.001 10.0   ; # 2
            1.0     0.001 10.0   ; # 3
            1.0     0.001 10.0   ; # 4
            1.0     0.001 10.0   ; # 5
            1.0     0.001 10.0   ; # 6
            1.0     0.001 10.0   ; # 7
            1.0     0.001 10.0   ; # 8
            0.20    0.001 10.0   ; # 9
            1.0     0.001 10.0   ; # 10
            0.65    0.001 10.0   ; # 11
            0.01    0.001 10.0   ; # 12
            0.25    0.001 10.0   ; # 13
            0.20    0.001 10.0   ; # 14
            0.90    0.001 10.0   ; # 15
            0.95    0.001 10.0   ; # 16
            0.90    0.001 10.0   ; # 17
            0.045   0.001 10.0   ; # 18
            0.01    0.001 10.0   ; # 19
        ];

    # set default set as the start -
    if (isnothing(pₒ) == true)
        P = length(κ[:,1])
        σ = 0.1 # move up to 10%
        pₒ = κ[:,1].*(1 .+ σ*rand(-1:1,P))
    end

    # setup the objective function -
    inner_optimizer = NelderMead()
    obj_function(p) =  loss_scalar(p, Y, model)
    results = optimize(obj_function, κ[:,2], κ[:,3], pₒ, Fminbox(inner_optimizer),
        Optim.Options(time_limit = 600, show_trace = true, show_every = 10, iterations=100))
    
    # grab the best parameters -
    p_best = Optim.minimizer(results)
    
    # run the sim w/the best parameters -
    # 1 - 9 : α vector
    model.α = p_best[1:10]
    G = model.G

    AT_idx = findfirst(x->x=="AT",model.total_species_list)
    TFPI_idx = findfirst(x->x=="TFPI",model.total_species_list)
    FIIa_idx = findfirst(x->x=="FIIa",model.total_species_list)
    AP_idx = findfirst(x->x=="AP",model.total_species_list)
    PL_idx = findfirst(x->x=="PL",model.total_species_list)
    FVIIa_idx = findfirst(x->x=="FVIIa",model.total_species_list)
    FXa_idx = findfirst(x->x=="FXa",model.total_species_list)
    FVa_idx = findfirst(x->x=="FVa",model.total_species_list)

    
   # r1 -
   G[TFPI_idx,1] = -1*p_best[11];

   # r2 -
   G[AP_idx,2] = p_best[12];

   # r4 -
   G[AP_idx,4] = p_best[13];

   # r5 -
   G[AP_idx,5] = p_best[14];
   G[FVIIa_idx,5] = p_best[15];

   # r6 -
   G[FXa_idx,6] = p_best[16];
   G[FVa_idx,6] = p_best[17];

   # r9 -
   G[AT_idx,9] = p_best[18];

   # r10 -
   G[FIIa_idx,10] = p_best[19];

    # run the model -
    (T,U) = evaluate(model,tspan=(0.0,120.0))
    data = [T U]
    Op = hcat(data)
    Yₘ = model_output_vector(T, Op[:,10]) # properties of the Thrombin curve 
    
    return (p_best, T, U, Yₘ)
end