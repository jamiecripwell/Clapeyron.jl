using XLSX
using CSV
using DataFrames
using Statistics

function DIPPR_const_props(comp::String)
    # Load the constant properties of a given component
    # comp: component name
    # Returns a dictionary of constant properties

    # Load the Excel file
    xlsx_path = "C:/Users/cripwell/OneDrive - Stellenbosch University/Documents/Research/Clapeyron/Clapeyron.jl/Projects/Chem_Info.xlsx"
    xlsx = XLSX.readxlsx(xlsx_path)

    # Access the "Chem_Info" tab
    chem_info_sheet = DataFrame(XLSX.readtable(xlsx_path, "Chem_Info"))
    # Check if the last character of the component name is a number
    if isnumeric(comp[end])
        # Find the ChemID corresponding to the component CAS number
        chem_id_row = chem_info_sheet[chem_info_sheet[:, "CASN"] .== comp, :]
    else
        # Find the ChemID corresponding to the component name (case insensitive)
        chem_id_row = chem_info_sheet[uppercase.(chem_info_sheet[:, "Name"]) .== uppercase(comp), :]
    end
    if nrow(chem_id_row) == 0
        error("Component not found in the Chem_Info tab")
    end
    chem_id = chem_id_row[1, 1]

    # Access the "Const_Values" tab
    const_values_sheet = DataFrame(XLSX.readtable(xlsx_path, "Const_Values"))

    properties_dict = Dict{String, Float64}()
    molar_props = ["VC","HFUS"]

    for const_prop ∈ ["MW","TC","PC","VC","MP","HFUS","ACEN","DM"]
        # Find the row corresponding to the ChemID and property
        const_row = const_values_sheet[(const_values_sheet[:, "ChemID"] .== chem_id) .& (const_values_sheet[:, "PropertyID"] .== const_prop), :]
        if nrow(const_row) == 0
            # If no matching rows are found, set the property value to NaN
            properties_dict[const_prop] = NaN
        elseif nrow(const_row) > 1            
            if nrow(const_row[.!ismissing.(const_row[:, "NoteID"]), :]) == 0
                # If none of the rows have an entry in the "NoteID" column, use the average of the "Const_Value" column
                properties_dict[const_prop] = Statistics.mean(const_row[:, "Const_Value"])
            elseif nrow(const_row[.!ismissing.(const_row[:, "NoteID"]), :]) > 1
                # If multiple rows have an entry in the "NoteID" column, use the average of the "Const_Value" column
                properties_dict[const_prop] = Statistics.mean(const_row[:, "Const_Value"])
            else
                # If only a single row is contains an entry in the "NoteID" column, use that row
                const_row = const_row[.!ismissing.(const_row[:, "NoteID"]), :]
                properties_dict[const_prop] = const_row[1, "Const_Value"]
            end
        else
            # If a single row is found, store the property value
            properties_dict[const_prop] = const_row[1, "Const_Value"]
        end
        if const_prop ∈ molar_props
            # Convert molar properties from kmol to mol basis
            properties_dict[const_prop] /= 1000
        elseif const_prop == "DM"
            # Convert dipole moment from Coulomb-meters to Debye
            properties_dict[const_prop] /= 3.33564e-30
        end
    end
    # Return the dictionary of constant properties
    return properties_dict
end


function DIPPR_calc(comp::String, T::Vector{Float64}, property::String, properties_dict::Dict{String, Float64} = Dict{String, Float64}())
    # Calculate the property of interest for a given component and temperature
    # comp: component name
    # T: vector of temperatures in K
    # property: property type, one of: 
    #   "VP" - vapour pressure
    #   "SVP" - solid vapour pressure
    #   "LDN" - saturated molar liquid density
    #   "SDN" - solid molar density
    #   "HVP" - heat of vaporization
    #   "ICP" - ideal gas isobaric heat capacity
    #   "LCP" - liquid isobaric heat capacity
    #   "SCP" - solid heat capacity
    #   "LVS" - liquid viscosity
    #   "VVS" - vapour viscosity
    #   "LTC" - liquid thermal conductivity
    #   "VTC" - vapour thermal conductivity
    #   "ST" - surface tension
    #   "SVR" - second virial coefficients
    # Tc: critical temperature in K, required for some properties
    # properties_dict: optional dictionary of properties
    # Returns the property value and the dictionary of properties

    # Check if the property is valid
    VALID_PROPERTIES = ["VP", "SVP", "LDN", "SDN", "HVP", "ICP", "LCP", "SCP", "LVS", "VVS", "LTC", "VTC", "ST", "SVR"]

    if !(property in VALID_PROPERTIES)
        error("Invalid property. Valid properties are: $(join(VALID_PROPERTIES, ", "))")
    end
    
    # Load the constant properties if not provided
    if isempty(properties_dict)
        properties_dict = DIPPR_const_props(comp)
    end

    # Load the Excel file
    xlsx_path = "C:/Users/cripwell/OneDrive - Stellenbosch University/Documents/Research/Clapeyron/Clapeyron.jl/Projects/Chem_Info.xlsx"
    xlsx = XLSX.readxlsx(xlsx_path)

    # Access the "Chem_Info" tab
    chem_info_sheet = DataFrame(XLSX.readtable(xlsx_path, "Chem_Info"))
    # Check if the last character of the component name is a number
    if isnumeric(comp[end])
        # Find the ChemID corresponding to the component CAS number
        chem_id_row = chem_info_sheet[chem_info_sheet[:, "CASN"] .== comp, :]
    else
        # Find the ChemID corresponding to the component name (case insensitive)
        chem_id_row = chem_info_sheet[uppercase.(chem_info_sheet[:, "Name"]) .== uppercase(comp), :]
    end
    if nrow(chem_id_row) == 0
        error("Component not found in the Chem_Info tab")
    end
    chem_id = chem_id_row[1, 1]

    # Access the "Tdep_Set_Info" tab
    tdep_set_info_sheet = DataFrame(XLSX.readtable(xlsx_path, "Tdep_Set_Info"))

    # Find the rows corresponding to the ChemID and property
    subset_rows = tdep_set_info_sheet[(tdep_set_info_sheet[:, "ChemID"] .== chem_id) .& (tdep_set_info_sheet[:, "PropertyID"] .== property), :]
    if nrow(subset_rows) == 0
        error("No matching rows found in the Tdep_Set_Info tab for the given ChemID and property")
    end

    # Find the row with a non-empty "Coeff_SetID"
    coeff_set_id_row = subset_rows[.!ismissing.(subset_rows[:, "Coeff_SetID"]), :]
    if nrow(coeff_set_id_row) == 0
        error("No Coeff_Set_ID found for the given ChemID and property")
    end

    coeff_set_id = coeff_set_id_row[1, "Coeff_SetID"]

    # Access the "Coefficients" and "Tdep_Coefficients" tabs
    coefficients_sheet      = DataFrame(XLSX.readtable(xlsx_path, "Coefficients"))
    tdep_coefficients_sheet = DataFrame(XLSX.readtable(xlsx_path, "Tdep_Coefficients"))

    # Find the EqnID and coefficients corresponding to the Coeff_SetID
    equation_id_row = tdep_coefficients_sheet[tdep_coefficients_sheet[:, "Coeff_SetID"] .== coeff_set_id, :]
    equation_id = equation_id_row[1, "EqnID"]

    coefficients_rows = coefficients_sheet[coefficients_sheet[:, "Coeff_SetID"] .== coeff_set_id, :]
    coeffs = zeros(5)
    for i ∈ 1:nrow(coefficients_rows)
        coeffs[i] = coefficients_rows[i, "Coeff_Value"]
    end

    if equation_id == 101
        # vapour pressure - Pa
        return exp.(coeffs[1] .+ coeffs[2] ./ T .+ coeffs[3] .* log.(T) .+ coeffs[4] .* T .^ coeffs[5]), properties_dict
    elseif equation_id == 100
        # liquid isobaric heat capacity - J/K.mol
        return (coeffs[1] .+ coeffs[2] .* T .+ coeffs[3] .* T .^ 2 .+ coeffs[4] .* T .^ 3 .+ coeffs[5] .* T .^ 4) ./ 1000, properties_dict
    elseif equation_id == 106
        # heat of vaporization - J/mol
        Tr = T ./ properties_dict["TC"]
        return coeffs[1] .* (1 .- Tr) .^ (coeffs[2] .+ coeffs[3] .* Tr .+ coeffs[4] .* Tr .^ 2 .+ coeffs[5] .* Tr .^ 3) ./ 1000, properties_dict
    elseif equation_id == 116
        # molar liquid density - mol/m3
        Tr = T ./ properties_dict["TC"]
        return (coeffs[1] .+ (1 .- Tr) .+ coeffs[2] .* (1 .- Tr) .^ 0.35 .+ coeffs[3] .* (1 .- Tr) .^ (2/3) .+ coeffs[4] .* (1 .- Tr) .+ coeffs[5] .* (1 .- Tr) .^ (4/3)) .* 1000, properties_dict
    elseif equation_id == 107
        # ideal gas isobaric heat capacity - J/K.mol
        return (coeffs[1] .+ coeffs[2] .* (coeffs[3] ./ T) ./ sinh.(coeffs[3] ./ T) .+ coeffs[4] .* (coeffs[5] ./ T) ./ cosh.(coeffs[5] ./ T)) ./ 1000, properties_dict
    elseif equation_id == 105
        # saturated molar liquid density - mol/m3
        return (coeffs[1] ./ (coeffs[2] .^ (1 .+ (1 .- T ./ coeffs[3]) .^ coeffs[4]))) .* 1000, properties_dict
    end
end

function DIPPR_calc(component::String, T::Float64, property::String)
    DIPPR_calc(component, [T], property)  # Wrap T in a vector and call the original function
end


# Test the function
# DIPPR_calc("74-82-8",298.,"VP")
# LDN, DIPPR_props = DIPPR_calc("methanol",[298.,303.],"LDN")
# println(LDN)
# println(DIPPR_props)
# print([Psat, rholiq, Hvap, ig_cp, l_cp])