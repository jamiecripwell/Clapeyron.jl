using Clapeyron, Metaheuristics, Plots

function logger(st)
    if st.iteration % 1 == 0
        sol = st.best_sol
        it = st.iteration
        println("Iteration: ",st.iteration)
        println("Best objective: ",st.best_sol.f)
        println("Best solution: ",st.best_sol.x)
    end
end

# Define the model
model = SAFTVRMie(["methanol","hexane"])

# Define parameters to be fit
toestimate = [
    Dict(
        :param => :epsilon,
        :indices => (2,2),
        :lower => 250.,
        :upper => 450.,
        :guess => 300.
    ),
    Dict(
        :param => :sigma,
        :indices => (2,2),
        :factor => 1e-10,
        :lower => 3.4,
        :upper => 4.2,
        :guess => 3.7
    )
    ,
    Dict(
        :param => :segment,
        :indices => 2,
        :lower => 1.5,
        :upper => 3.0,
        :guess => 1.
    ),
    Dict(
        :param => :lambda_r,
        :indices => (2,2),
        :lower => 12.,
        :upper => 20.,
        :guess => 16.
    ),
    Dict(
        :param => :epsilon,
        :indices => (1,2),
        :lower => 250.,
        :upper => 400.,
        :guess => 350.
    )
]


# Define property estimation functions
function saturation_P_and_rho(model::EoSModel,T)
    sat = saturation_pressure(model,T)
    return sat[1], 1/sat[2]
end
     
function bubble_point(model::EoSModel,T,x)
    bub = bubble_pressure(model,T,[x,1-x])
    return bub[1], bub[4][1]
end

method = ECA(;options=Options(f_tol_rel=1e-2));

# Construct estimator
estimator,objective,initial,upper,lower = Estimation(model,toestimate,["C:/Users/cripwell/OneDrive - Stellenbosch University/Documents/Research/Clapeyron/Clapeyron.jl/Projects/data/hex_saturation_pressure_liquid_density.csv",
                                                                       "C:/Users/cripwell/OneDrive - Stellenbosch University/Documents/Research/Clapeyron/Clapeyron.jl/Projects/data/methanol_hex_343K.csv"]);

# Perform optimization
params, model = optimize(objective, estimator, method; verbose = true, logger = logger)

export_model(model);

## Plot results
Nexp = length(estimator.data)
    # plt = plot(grid=:off,framestyle=:box,foreground_color_legend = nothing,legend_font=font(12))

for i in 1:Nexp
    x = estimator.data[i].inputs
    y_exp = estimator.data[i].outputs[1]
    prop = estimator.data[i].method
    species = estimator.data[i].species

    idx_r = zeros(length(model))
    for i in species
        idx_r += model.components .== i
    end
    model_r = index_reduction(model,idx_r)[1]
    if length(x) == 1
        x = [x[1][i] for i in 1:length(y_exp)]
        y = prop.(model_r,x)
    elseif length(x) == 2
        x1 = [x[1][i] for i in 1:length(y_exp)]
        x = [x[2][i] for i in 1:length(y_exp)]
        y = prop.(model_r,x1,x)
    end

    if length(y[1]) == 1
        plt = plot(grid=:off,framestyle=:box,foreground_color_legend = nothing,legend_font=font(12))
        plot!(plt, x, y_exp, label="Experimental", seriestype=:scatter, color=:black, markersize=3)
        plot!(plt, x, y, label="Model", color=:red)
        savefig(plt,"fit_$(species[1])_$(prop).png")
    elseif length(y[1]) >= 2
        for j in 1:length(y[1])
            _y = [y[k][j] for k in 1:length(y)]
            plt = plot(grid=:off,framestyle=:box,foreground_color_legend = nothing,legend_font=font(12))
            plot!(plt, x, estimator.data[i].outputs[j], label="Experimental", seriestype=:scatter, color=:black, markersize=3)
            plot!(plt, x, _y, label="Model", color=:red)
            savefig(plt,"fit_$(species[1])_$(prop)_$(j).png")
        end
    end
end