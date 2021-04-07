function rh = comp_rh(type, srfc, par)

    if any(strcmp(type, {'era5', 'erai', 'era5c'}))
        e = calc_esat(srfc.d2m, par.frz); % calculate actual vapor pressure
        esat = calc_esat(srfc.t2m, par.frz); % saturation vapor pressure
        rh = e./esat;
    end

end