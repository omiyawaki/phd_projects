function ma_out = calc_ma_hurs(ma_in, plev, par)
    if isnan(ma_in.ps) | isnan(ma_in.tas) | isnan(ma_in.hurs)
        ma_out = nan(size(plev)); % if any of the initial values are nan, ouput nan profile
    else
        if isfield(ma_in, 'pinit')
            pa_cd(1) = ma_in.pinit;
        else
            pa_cd(1) = ma_in.ps;
        end
        ta_cd(1) = ma_in.tas; % set initial temperature
        esat_cd(1) = calc_esat(ta_cd(1), par.frz); % calculate initial saturation vapor pressure
        qsat_cd(1) = calc_q(pa_cd(1), esat_cd(1), par); % calculate initial specific humidity
        rh(1) = ma_in.hurs/100; % set initial relative humidity

        r(1) = par.eps * rh(1) * esat_cd(1) / (pa_cd(1) - rh(1) * esat_cd(1));
        deriv_dry(1) = dry(pa_cd(1), ta_cd(1), par);
        [pa_d, ta_pa_d]=ode45(@(pa, T) dry(pa, T, par), par.pa_span, ta_cd(1));
        i = 2;
        while rh(i-1) < 1
            idx_lcl = i;
            pa_cd(i) = pa_cd(i-1) + par.dpa;
            ta_cd(i, 1) = ta_cd(i-1) + par.dpa * deriv_dry(i-1);
            esat_cd(i, 1) = calc_esat(ta_cd(i), par.frz);
            qsat_cd(i, 1) = calc_q(pa_cd(i), esat_cd(i), par);
            r(i, 1) = r(i-1);
            rh(i, 1) = r(i) * pa_cd(i) / (esat_cd(i) * (r(i) + par.eps));
            deriv_dry(i, 1) = dry(pa_cd(i), ta_cd(i), par);
            i = i + 1;
        end
        pa_lcl = pa_cd(end);
        if strcmp(par.ma_type, 'reversible')
            par.r_t = r(end); % moisture at LCL is the total moisture of rising parcel
            [pa_cs, ta_cs] = ode45(@(pa, T) sat_rev(pa, T, par.frz, par), [pa_lcl, par.pa_span(end)], ta_cd(end));
        elseif strcmp(par.ma_type, 'pseudo')
            [pa_cs, ta_cs] = ode45(@(pa, T) sat_pseudo(pa, T, par.frz, par), [pa_lcl, par.pa_span(end)], ta_cd(end));
        elseif strcmp(par.ma_type, 'std')
            [pa_cs, ta_cs] = ode45(@(pa, T) sat_std(pa, T, par.frz, par), [pa_lcl, par.pa_span(end)], ta_cd(end));
        else
            error(sprintf('The chosen moist adiabat type, %s, is not available. Choose between reversible, pseudo, or std.', par.ma_type));
        end
        for i = 1:length(pa_cs)
            if strcmp(par.ma_type, 'reversible')
                [deriv_cs(i, 1), qsat_cs(i, 1)] = sat_rev(pa_cs(i), ta_cs(i), par.frz, par);
            elseif strcmp(par.ma_type, 'pseudo')
                [deriv_cs(i, 1), qsat_cs(i, 1)] = sat_pseudo(pa_cs(i), ta_cs(i), par.frz, par);
            elseif strcmp(par.ma_type, 'std')
                [deriv_cs(i, 1), qsat_cs(i, 1)] = sat_std(pa_cs(i), ta_cs(i), par.frz, par);
            else
                error(sprintf('The chosen moist adiabat type, %s, is not available. Choose between reversible, pseudo, or std.', par.ma_type));
            end
            deriv_csd(i, 1) = dry(pa_cs(i), ta_cs(i), par);
        end
        ta_s(:, 1) = [ta_cd; ta_cs(2:end)];
        pa_s(:, 1) = [pa_cd(:); pa_cs(2:end)];
        qsat_s(:, 1) = [qsat_cd(:, 1); qsat_cs(2:end)];
        dtadpa_s(:, 1) = [deriv_dry; deriv_cs(2:end)];

        ta_s = interp1(pa_s, ta_s, plev, 'spline', nan);
        qsat_s = interp1(pa_s, qsat_s, plev, 'spline', nan);
        dtadpa_s = interp1(pa_s, dtadpa_s, plev, 'spline', nan);

        % OUTPUT
        ma_out = ta_s;
    end

end % compute moist adiabat based on RH
