# include -
include("Include.jl")
size_t = 12001;

# visit_1_FIIa = Array{Float64}(undef,(size_t,1))
# sd_visit_1_FIIa = Array{Float64}(undef,(size_t,1))
visit_4_FIIa = Array{Float64}(undef,(size_t,1))
sd_visit_4_FIIa = Array{Float64}(undef,(size_t,1))

SIM_FIIa = Array{Float64}(undef,(size_t,12))
Z = zeros(size_t)
fVIIa = 25;
prothrombin = 150;

T_df = CSV.read(joinpath(_PATH_TO_TMP_CONSTRUCTION,"SIM-Hemophilia-1.csv"),DataFrame)
Time = T_df[1:size_t,:Time]
for i ∈ 1:12
    FIIa_df = CSV.read(joinpath(_PATH_TO_TMP_CONSTRUCTION,"SIM-Hemophilia-$(i).csv"),DataFrame)
    Temp_Thrombin = FIIa_df[1:size_t,:FIIa]
    SIM_FIIa[:,i] = Temp_Thrombin
end

visit_4_FIIa = mean(SIM_FIIa,dims = 2)
sd_visit_4_FIIa = std(SIM_FIIa,dims = 2)
U = sd_visit_4_FIIa
plot(Time,visit_4_FIIa,color=colorant"#0077BB",label="",lw=2,fillrange=(visit_4_FIIa-U,U+visit_4_FIIa),fillalpha=0.4)
xlabel!("Time (min)")
ylabel!("Activated Thrombin FIIa (nM)")
plot!([-5],[-5],xlim=(0.0,120.0),line=:scatter,markerstrokecolor=colorant"#0077BB",color=colorant"#0077BB",label = "",foreground_color_legend = nothing)

exp_df = CSV.read(joinpath(_PATH_TO_DATA,"thrombin-hemophilia-exp.csv"),DataFrame)
scatter!(exp_df[:,:Time],exp_df[:,:FIIa],color=colorant"#0077BB",markerstrokecolor=colorant"#0077BB",label="",yerr = 0.1.*exp_df[:,:FIIa])

#for i ∈ 1:10
#    FIIa_df = CSV.read(joinpath(_PATH_TO_TMP_CONSTRUCTION,"SIM-visit-1-TF-$(i).csv"),DataFrame)
#    Temp_Thrombin = FIIa_df[1:size_t,:FIIa]
#    SIM_FIIa[:,i] = Temp_Thrombin
#end

#visit_1_FIIa = mean(SIM_FIIa,dims = 2)
#sd_visit_1_FIIa = std(SIM_FIIa,dims = 2)
#U = sd_visit_1_FIIa
#plot!(Time,visit_1_FIIa,color=colorant"#EE7733",label="",lw=2,fillrange=(visit_1_FIIa-U,U+visit_1_FIIa),fillalpha=0.4)
#y = plot!([-5],[-5],xlim=(0.0,120.0),line=:scatter,color=colorant"#FFFFFF",markerstrokecolor = colorant"#FFFFFF",label = "rfVIIa: $(fVIIa) nM",foreground_color_legend = nothing)
y = plot!([-5],[-5],xlim=(0.0,120.0),line=:scatter,color=colorant"#FFFFFF",markerstrokecolor = colorant"#FFFFFF",label = "Hemophilia",foreground_color_legend = nothing)
savefig(y,joinpath(pwd(),"figs","Actual-Hemophilia.pdf"))