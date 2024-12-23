function qt_f0_p!(K,z,p,ps,β0)
    K .= ps ./ p
    return rachfordrice(K,z) - β0
end


function qt_flash_x0(model,β,T,z,method::FlashMethod)
    if method.p0 == nothing
        if 0 <= β <= 0.01
            x = z ./ sum(z)
            p,vl,vv,y = __x0_bubble_pressure(model,T,x)
            y ./= sum(y)
            
            return FlashResult(p,T,SA[x,y],SA[1.0-β,1.0*β],SA[vl,vv],sort = false)
        elseif 0.99 <= β <= 1.0
            y = z ./ sum(z)
            p,vl,vv,x = __x0_dew_pressure(model,T,y)
            x ./= sum(x)
            return FlashResult(p,T,SA[x,y],SA[1.0-β,1.0*β],SA[vl,vv],sort = false)
        else
            pures = split_model(model)
            sat = extended_saturation_pressure.(pures,T)
            ps = first.(sat)    
            K = similar(ps)
            pmin,pmax = extrema(ps)
            p0 = β*pmin + (1-β)*pmax
            fp(p) = qt_f0_p!(K,z,p,ps,β)
            prob = Roots.ZeroProblem(fp,p0)
            p = Roots.solve(prob)
        end
    else
        p = method.p0
    end
    res =  pt_flash_x0(model,p,T,z,method;k0 = :suggest)
    return res
end

function qt_flash(model::EoSModel,β,T,z;kwargs...)
    method = init_preferred_method(qt_flash,model,kwargs)
    return qt_flash(model,β,T,z,method)
end

function init_preferred_method(method::typeof(qt_flash),model::EoSModel,kwargs) 
    GeneralizedXYFlash(;kwargs...)
end

function qt_flash(model,β,T,z,method::FlashMethod)
    check_arraysize(model,z)
    if supports_reduction(method)
        model_r,idx_r = index_reduction(model,z)
        z_r = z[idx_r]
        method_r = index_reduction(method,idx_r)
    else
        model_r,idx_r = model,trues(length(model))
        method_r,z_r = method,z
    end
    if length(model_r) == 1
        result1 = βflash_pure(model,temperature,T,βv,z)
        return index_expansion(result1,idx_r)
    end
    
    result = qt_flash_impl(model_r,β,T,z_r,method_r)
    if !issorted(result.volumes)
        #this is in case we catch a bad result.
        result = FlashResult(result)
    end
    ∑β = sum(result.fractions)
    result.fractions ./= ∑β
    result.fractions .*= sum(z)
    return index_expansion(result,idx_r)
end

function qt_flash_impl(model,β,T,z,method::GeneralizedXYFlash)
    flash0 = qt_flash_x0(model,β,T,z,method)
    isone(numphases(flash0)) && return flash0
    spec = FlashSpecifications(Vfrac(2),β,temperature,T)
    return xy_flash(model,spec,z,flash0,method)
end

function bubble_pressure_impl(model::EoSModel,T,z,method::GeneralizedXYFlash)
    result = Clapeyron.qt_flash(model,0,T,z,method)
    x1,x2 = result.compositions
    v1,v2 = result.volumes
    if x1 ≈ z
        y = x2
        vl,vv = v1,v2
    else
        y = x1
        vl,vv = v2,v1
    end
    return pressure(result),vl,vv,y
end

function dew_pressure_impl(model::EoSModel,T,z,method::GeneralizedXYFlash)
    result = Clapeyron.qt_flash(model,1,T,z,method)
    x1,x2 = result.compositions
    v1,v2 = result.volumes
    if x1 ≈ z
        x = x2
        vl,vv = v2,v1
    else
        x = x1
        vl,vv = v1,v2
    end
    return pressure(result),vl,vv,x
end

export qt_flash
