function ma_out = calc_maz_dew_pa(ma_in, plev, par, type)
    if any(strcmp(type, {'erai', 'era5', 'era5c'}))
        ps = ma_in.sp; tas = ma_in.t2m; tdas = ma_in.d2m;
    elseif strcmp(type, 'echam')
        ps = ma_in.aps; tas = ma_in.temp2; tdas = ma_in.dew2;
    end

    if isnan(ps) | isnan(tas) | isnan(tdas)
        ma_out = nan(size(z_int)); % if any of the initial values are nan, ouput nan profile
    else
        % integrate temperature profiles over height
        z_cd(1) = ma_in.zs;

    if isfield(ma_in, 'pinit')
        pa_cd(1) = ma_in.pinit;
    else
        pa_cd(1) = ps;
    end
        ta_cd(1) = tas;
        esat_cd(1) = calc_esat(ta_cd(1), par.frz);
        qsat_cd(1) = calc_q(pa_cd(1), esat_cd(1), par);

        % convert dew point temperature to relative humidity
        e_cd(1) = calc_esat(tdas, par.frz); % calculate initial vapor pressure
        rh(1) = e_cd(1)/esat_cd(1); % set initial relative humidity

        r(1) = par.eps * rh(1) * esat_cd(1) / (pa_cd(1) - rh(1) * esat_cd(1));
        dtdz_dpdz(1,:) = dry_z([ta_cd(1), pa_cd(1)], par);
        T_p0 = [ta_cd(1), pa_cd(1)];
        [z_d, ta_pa_d]=ode45(@(z, T_p) dry_z(T_p, par), par.z_span, T_p0);
        i = 2;
        while rh(i-1) < 1
            idx_lcl = i;
            z_cd(i, 1) = z_cd(i-1) + par.dz;
            pa_cd(i, 1) = pa_cd(i-1) + par.dz * dtdz_dpdz(i-1, 2);
            ta_cd(i, 1) = ta_cd(i-1) + par.dz * dtdz_dpdz(i-1, 1);
            esat_cd(i, 1) = calc_esat(ta_cd(i), par.frz);
            qsat_cd(i, 1) = calc_q(pa_cd(i), esat_cd(i), par);
            r(i, 1) = r(i-1);
            rh(i, 1) = r(i) * pa_cd(i) / (esat_cd(i) * (r(i) + par.eps));
            dtdz_dpdz(i,:) = dry_z([ta_cd(i), pa_cd(i)], par);
            i = i + 1;
        end
        z_lcl = z_cd(end);
        pa_lcl = pa_cd(end);
        ta_pa_cd = [ta_cd, pa_cd];
        deriv_dry(:,1) = dtdz_dpdz(:,1);
        if strcmp(par.ma_type, 'reversible')
            par.r_t = r(end); % moisture at LCL is the total moisture of rising parcel
            [z_cs, ta_pa_cs] = ode45(@(z, T_p) sat_rev_z(T_p, par.frz, par), [z_lcl, par.z_span(end)], ta_pa_cd(end,:));
            % h1 = interp1(ta_pa_cs(:,1), z_cs, 240); % middle of troposphere
        elseif strcmp(par.ma_type, 'pseudo')
            [z_cs, ta_pa_cs] = ode45(@(z, T_p) sat_pseudo_z(T_p, par.frz, par), [z_lcl, par.z_span(end)], ta_pa_cd(end,:));
            % h1 = interp1(ta_pa_cs(:,1), z_cs, 240); % middle of troposphere
        elseif strcmp(par.ma_type, 'std')
            [z_cs, ta_pa_cs] = ode45(@(z, T_p) sat_z(T_p, par.frz, par), [z_lcl, par.z_span(end)], ta_pa_cd(end,:));
            % h1 = interp1(ta_pa_cs(:,1), z_cs, 240); % middle of troposphere
        else
            error(sprintf('The chosen moist adiabat type, %s, is not available. Choose between reversible, pseudo, or std.', zr.ma_type));
        end
        for i = 1:length(z_cs)
            if strcmp(par.ma_type, 'reversible')
                [dt_dp_cs(i, :), qsat_cs(i, 1)] = sat_rev_z(ta_pa_cs(i,:), par.frz, par);
            elseif strcmp(par.ma_type, 'pseudo')
                [dt_dp_cs(i, :), qsat_cs(i, 1)] = sat_pseudo_z(ta_pa_cs(i,:), par.frz, par);
            elseif strcmp(par.ma_type, 'std')
                [dt_dp_cs(i, :), qsat_cs(i, 1)] = sat_z(ta_pa_cs(i,:), par.frz, par);
            else
                error(sprintf('The chosen moist adiabat type, %s, is not available. Choose between reversible, pseudo, or std.', par.ma_type));
            end
            dt_dp_csd(i, :) = dry_z([ta_pa_cs(i,:)], par);
        end
        deriv_cs(:,1) = dt_dp_cs(:,1);
        deriv_csd(:,1) = dt_dp_csd(:,1);
        z_s(:, 1) = [z_cd; z_cs(2:end)];
        ta_s(:, 1) = real([ta_cd; ta_pa_cs(2:end,1)]);
        pa_s(:, 1) = [pa_cd(:); ta_pa_cs(2:end,2)];
        qsat_s(:, 1) = [qsat_cd(:, 1); qsat_cs(2:end)];
        dtadz_s(:, 1) = [deriv_dry; deriv_cs(2:end)];

        nanfilter = isnan(real(pa_s));
        pa_s(nanfilter) = [];
        dtadz_s(nanfilter) = [];

        dtadz_s = -1e3*interp1(real(pa_s), real(dtadz_s), plev, 'spline', nan); % output in K/km at standard pressure grid

        % output temperature
        ma_out = dtadz_s;
    end

end
