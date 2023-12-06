function mean_ionic_activity_coefficient(model::ESElectrolyteModel,salts,T,m,zsolvent=[1.])
    isolvent = model.icomponents[model.charge.==0]
    iions = model.icomponents[model.charge.!=0]
    ions = model.components[model.charge.!=0]

    method = FugBubblePressure(nonvolatiles=ions)

    ν = salt_stoichiometry(model,salts)

    z0 = molality_to_composition(model,salts,ones(length(m)).*1e-10,zsolvent)
    z = molality_to_composition(model,salts,m,zsolvent)

    (p0,vl0,vv0,y0) = bubble_pressure(model,T,z0,method)
    φ0 = VT_fugacity_coefficient(model,vl0,T,z0)[iions]

    (p,vl,vv,y) = bubble_pressure(model,T,z,method)
    φ = VT_fugacity_coefficient(model,vl,T,z)[iions]

    γim = φ./φ0.*sum(z[isolvent])/sum(z)
    γsm = (prod(γim'.^ν,dims=2)).^(1 ./sum(ν,dims=2))
    return γsm
end

function mean_ionic_activity_coefficient(model::ESElectrolyteModel,salts,p,T,m,zsolvent=[1.])
    isolvent = model.icomponents[model.charge.==0]
    iions = model.icomponents[model.charge.!=0]
    ions = model.components[model.charge.!=0]

    ν = salt_stoichiometry(model,salts)

    z0 = molality_to_composition(model,salts,ones(length(m)).*1e-10,zsolvent)
    z = molality_to_composition(model,salts,m,zsolvent)

    φ0 = fugacity_coefficient(model,p,T,z0)[iions]

    φ = fugacity_coefficient(model,p,T,z)[iions]

    γim = φ./φ0.*sum(z[isolvent])/sum(z)
    γsm = (prod(γim'.^ν,dims=2)).^(1 ./sum(ν,dims=2))
    return γsm
end

function osmotic_coefficient(model::ESElectrolyteModel,salts,T,m,zsolvent=[1.])
    isolvent = model.icomponents[model.charge.==0]
    iions = model.icomponents[model.charge.!=0]
    ions = model.components[model.charge.!=0]

    method = FugBubblePressure(nonvolatiles=ions)

    ν = salt_stoichiometry(model,salts)

    z0 = molality_to_composition(model,salts,ones(length(m)).*1e-10,zsolvent)
    z = molality_to_composition(model,salts,m,zsolvent)

    (p0,vl0,vv0,y0) = bubble_pressure(model,T,z0,method)
    φ0 = VT_fugacity_coefficient(model,vl0,T,z0)[isolvent]

    (p,vl,vv,y) = bubble_pressure(model,T,z,method)
    φ = VT_fugacity_coefficient(model,vl,T,z)[isolvent]

    asolv = φ./φ0.*z[isolvent]/sum(z)
    Mw = model.neutralmodel.params.Mw.values[isolvent].*1e-3
    # println(sum(ν.*m).*Mw)
    println(φ./φ0)
    return -1 ./(sum(ν.*m).*Mw).*log.(asolv)
end

export mean_ionic_activity_coefficient, osmotic_coefficient