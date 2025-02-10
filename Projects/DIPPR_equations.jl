using CSV
using DataFrames

function DIPPR_calc(comp::string, T::Float64, property::string, Tc::Float64 = 273.15)
    # Calculate the property of interest for a given component and temperature
    # comp: component name
    # T: temperature in K
    # property: property type, one of: 
    #   "P_sat", "rho_sat", "H_vap", 
    #   "ig_cp", "l_cp", "2vc"
    # Tc: critical temperature in K, required for some properties
    # Returns the property value


    # Load the CSV file
    file_path = "/c:/Users/cripwell/OneDrive - Stellenbosch University/Documents/Research/Clapeyron/Clapeyron.jl/Projects/DIPPR_corrs.csv"
    data = CSV.read(file_path, DataFrame)

    # Determine the property type
    prop_type = ""
    if property == "P_sat"
        prop_type = "vapour pressure"
    elseif property == "rho_sat"
        prop_type = "liquid density"
    elseif property == "H_vap"
        prop_type = "heat of vaporization"
    elseif property == "ig_cp"
        prop_type = "ideal gas isobaric heat capacity"
    elseif property == "l_cp"
        prop_type = "liquid isobaric heat capacity"
    elseif property == "2vc"
        prop_type = "second virial coefficients"
    else
        error("Invalid property type")
    end

    # Find the row corresponding to the component and property type
    row = data[(data[:, 1] .== comp) .& (data[:, 2] .== prop_type), :]

    if nrow(row) == 0
        error("Component or property type not found in the CSV file")
    end

    # Store the equation number and coefficients A to E
    equation_num = row[1, 3]
    A, B, C, D, E = row[1, 6:10]

    if equation_num == 101
        # vapour pressure - Pa
        return exp(A + B./T + C.*log(T) + D.*T.^E);
    elseif equation_num == 100
        # liquid isobaric heat capacity - J/K.mol
        # TODO: Check the units of the output
        return A + B.*T + C.*√T + D.*T.^3 + E.*T.^4;
    elseif equation_num == 106
        # heat of vaporization - J/mol
        Tr = T/Tc;
        return A .* (1 - Tr)^(B + C.*Tr + D.*Tr.^2 + E.*Tr.^3) ./ 1000;
    elseif equation_num == 116
        # molar liquid density - mol/m3
        Tr = T/Tc;
        return (A + (1 - Tr) + B.*(1 - Tr).^0.35 + C.*(1 - Tr).^(2/3) + D.*(1 - Tr) + E.*(1 - Tr).^(4/3)) .* 1000;
    elseif equation_num == 107
        # ideal gas isobaric heat capacity - J/K.mol
        return (A + B.*(C./T)/(sinh(C./T)) + D.*(E./T)/(cosh(E./T))) ./ 1000;
    elseif equation_num == 105
        # saturated molar liquid density - mol/m3
        return (A ./ (B.^(1 + (1 - T./C).^D))) .* 1000;
end