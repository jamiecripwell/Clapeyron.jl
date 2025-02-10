using Clapeyron, Metaheuristics, PyCall
import PyPlot as plt

model = SAFTVRMie(["ethanol"])

toestimate = [
    Dict(
        :param => :epsilon,
        :lower => 50.,
        :upper => 500.,
        :guess => 238.5
    ),
    Dict(
        :param => :sigma,
        :factor => 1e-10,
        :lower => 3.0,
        :upper => 4.5,
        :guess => 3.80
    ),
    Dict(
        :param => :segment,
        :lower => 1.,
        :upper => 5.,
        :guess => 3.48
    )
];

function saturation_p_rhol(model::EoSModel,T)
    sat = saturation_pressure(model,T)
    return sat[1], 1/sat[2]
end

method = ECA();

estimator,objective,initial,upper,lower = Estimation(model,toestimate,["Projects\\data\\ethanol_sat.csv"])

params, model = optimize(objective, estimator, method)