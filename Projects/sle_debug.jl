using Clapeyron, GCIdentifier, ChemicalIdentifiers

oxalic_acid_smiles = search_chemical("oxalic acid").smiles
lactic_acid_smiles = search_chemical("lactic acid").smiles
glycine_smiles = search_chemical("glycine").smiles

(~, oxalic_acid_groups,    oxalic_acid_connectivity)    = get_groups_from_smiles(oxalic_acid_smiles,    SAFTgammaMieGroups, connectivity=true)
(~, lactic_acid_groups,    lactic_acid_connectivity)    = get_groups_from_smiles(lactic_acid_smiles,    SAFTgammaMieGroups, connectivity=true)
(~, glycine_groups,        glycine_connectivity)        = get_groups_from_smiles(glycine_smiles,       SAFTgammaMieGroups, connectivity=true)

ox_gly = CompositeModel([("oxalic acid"=>oxalic_acid_groups),("glycine"=>glycine_groups)];fluid=SAFTgammaMie,solid=SolidHfus)
la_gly = CompositeModel([("lactic acid"=>lactic_acid_groups),("glycine"=>glycine_groups)];fluid=SAFTgammaMie,solid=SolidHfus)
models = [ox_gly, la_gly]

M = length(models)
TE = zeros(M,1)
sol = zeros(10,2,M)
T = zeros(10,2,M)

for k in 1:1
    
    model = models[k]

    # Use the melting temperatures as the starting points
    Tm = model.solid.params.Tm.values

    # Obtain the eutetic temperature as the end point
    TE[k] = eutectic_point(model)[1]

    # Consider each side of the SLE curve separately
    for i in 1:2
        T[:,i,k] = LinRange(TE[k],Tm[i]*0.9999,10)
        for j in 1:length(T[:,i,k])
            sol[j,i,k] = sle_solubility(model,1e5,T[j,i,k],ones(2);solute=[model.components[i]])[1]
        end
        
    end
    
end