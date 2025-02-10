using Clapeyron, Metaheuristics, PyCall, LaTeXStrings, Plots, CSV, DataFrames
import PyPlot as plt

pyplot()
# Enable LaTeX globally
PyPlot.matplotlib[:rc]("text", usetex=true);

methanol_VRMie_reg = SAFTVRMie(["methanol", "hexane"]);

toestimate = [
    Dict(
        :param => :epsilon,
        :lower => 50.,
        :upper => 500.,
        :guess => 167.7
    ),
    Dict(
        :param => :sigma,
        :factor => 1e-10,
        :lower => 3.0,
        :upper => 4.5,
        :guess => 3.31
    ),
    Dict(
        :param => :segment,
        :lower => 1.,
        :upper => 5.,
        :guess => 1.53
    ),
    Dict(
        :param => :lambda_r,
        :lower => 6.0,
        :upper => 30.0,
        :guess => 8.64
    ),
    Dict(
        :param => :epsilon_assoc,
        :lower => 1000.,
        :upper => 3500.,
        :guess => 2852.1
    ),
    Dict(
        :param => :bondvol,
        :lower => 1.0e-29,
        :upper => 2.0e-28,
        :guess => 1.0657e-28
    )
]

function saturation_p_rhol(model_multi::EoSModel,T)
    params_vector = fieldnames(typeof(model_multi.params))

    user_locations_string = ";"
    for i ∈ 1:length(params_vector)
            param_name = params_vector[i]
            if getfield(model_multi.params, params_vector[i]).values isa AbstractArray{Float64}
                if param_name == :sigma
                    param_value = getfield(model_multi.params, param_name).values[1].*1e10
                    user_locations_string *= "$param_name=$param_value,"
                else
                    param_value = getfield(model_multi.params, param_name).values[1]
                    user_locations_string *= "$param_name=$param_value,"
                end
            else
                param_value = getfield(model_multi.params, param_name).values.values[1]
                comp_name = model_multi.components[1]
                user_locations_string *= """$param_name=Dict((("$comp_name","e"),("$comp_name","H")) => $param_value),"""
            end
    end
    user_locations_string = user_locations_string[1:end-1]

     
    model_type_name = nameof(typeof(model_multi));
    # model_type = eval(Symbol(model_type_name))
    comps = [model_multi.components[1]]
    # model = model_type(comps)

     # Construct the full expression to be evaluated
     expression = "$model_type_name($comps; userlocations=($user_locations_string))"

     # Evaluate the expression
     model = eval(Meta.parse(expression))

    # println(model_multi.params.sigma.values)
    # println(model.params.sigma.values)

    sat = saturation_pressure(model,T)
    return sat[1], 1/sat[2]
end

function saturation_p_rho(model::EoSModel,T)
    sat = saturation_pressure(model,T)
    return sat[1], 1/sat[2]
end

function bubble_point(model::EoSModel,T,x)
    bub = bubble_pressure(model,T,[x,1-x])
    return bub[1], bub[4][1]
end

method = ECA();

# estimator,objective,initial,upper,lower = Estimation(ethanol_VRMie_reg,toestimate,["data/ethanol_sat.csv","data/ethanol_hept_333K.csv"]);
# estimator,objective,initial,upper,lower = Estimation(ethanol_VRMie_reg,toestimate,["data/ethanol_hept_333K.csv"]);
estimator,objective,initial,upper,lower = Estimation(methanol_VRMie_reg,toestimate,["data/methanol_sat.csv","data/methanol_hex_343K.csv"]);
# estimator,objective,initial,upper,lower = Estimation(methanol_VRMie_reg,toestimate,["data/methanol_sat.csv"]);

# params, ethanol_VRMie_reg = optimize(objective, estimator, method)
params, methanol_VRMie_reg = optimize(objective, estimator, method)

export_model(methanol_VRMie_reg);
