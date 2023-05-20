function loss_scalar(κ::Array{Float64,1}, Y::Array{Float64,1},  model::BSTModel)

    # get the G matrix -
    G = model.G

    # κ map
    # 1 - 10 : α vector
    model.α = κ[1:10]

    AT_idx = findfirst(x->x=="AT",model.total_species_list)
    TFPI_idx = findfirst(x->x=="TFPI",model.total_species_list)
    FIIa_idx = findfirst(x->x=="FIIa",model.total_species_list)
    AP_idx = findfirst(x->x=="AP",model.total_species_list)
    PL_idx = findfirst(x->x=="PL",model.total_species_list)
    FVIIa_idx = findfirst(x->x=="FVIIa",model.total_species_list)
    FXa_idx = findfirst(x->x=="FXa",model.total_species_list)
    FVa_idx = findfirst(x->x=="FVa",model.total_species_list)

    
   # r1 -
   G[TFPI_idx,1] = -1*κ[11];

   # r2 -
   G[AP_idx,2] = κ[12];

   # r4 -
   G[AP_idx,4] = κ[13];

   # r5 -
   G[AP_idx,5] = κ[14];
   G[FVIIa_idx,5] = κ[15];

   # r6 -
   G[FXa_idx,6] = κ[16];
   G[FVa_idx,6] = κ[17];

   # r9 -
   G[AT_idx,9] = κ[18];

   # r10 -
   G[FIIa_idx,10] = κ[19];


    # run the model -
    (T,U) = evaluate(model,tspan=(0.0,120.0))
    data = [T U]
    Op = hcat(data)
    Yₘ = model_output_vector(T, Op[:,10]) # properties of the Thrombin curve 
    ϵ = norm((Y .- Yₘ).^2)

    # @info ϵ

    # return -
    return ϵ
end