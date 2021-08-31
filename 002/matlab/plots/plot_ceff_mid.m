function plot_ceff_mid(type, par)
    
    make_dirs(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/flux_z.mat', prefix_proc)); % load lat x mon RCAE_ALT data
    load(sprintf('%s/ta_mon_lat.mat', prefix_proc)); % load lat x mon RCAE_ALT data

    tas = tasi.lo(:,:,1); % extract just surface temp
    clear tasi

    % tmon = ([1:12] - 0.5) /12 * 365 * 86400; % express monthly time in terms of seconds
    tmon = linspace(0,365,12) * 86400; % express monthly time in terms of seconds
    teval = linspace(0,365,100) * 86400; 
    omega = 2*pi/(365*86400);
    par.rho = 1000; % density of water kg m^-3
    % par.rho = 997; % density of water kg m^-3
    par.cw = 4000; % specific heat capacity of water J kg^-1 K^-1

    lat_center = [50 -50];
    par.lat_bound = 10;
    
    for f = {'mse_old'}; fw = f{1};
        for l = par.land_list; land = l{1};
            if strcmp(land, 'lo'); land_text = 'Land + Ocean';
            elseif strcmp(land, 'l'); land_text = 'Land';
            elseif strcmp(land, 'o'); land_text = 'Ocean';
            end

            for lb = 1:length(lat_center); lat_c = lat_center(lb);

                [lat, clat, clat_mon, ~] = make_midlatitude_lat(lat_c, par);

                folder = sprintf('%s/ceff/%s/%s/0_lat_%g_%g', plotdir, fw, land, lat(1), lat(end));
                if ~exist(folder, 'dir'); mkdir(folder); end;
                
                % interpolate at select latitude
                q = interp1(grid.dim2.lat, -flux_z.lo.swsfc, lat);
                q = nansum(clat_mon.*q,1) / nansum(clat);
                q_mon = squeeze(q);

                t = interp1(grid.dim2.lat, tas, lat);
                t = nansum(clat_mon.*t,1) / nansum(clat);
                t_mon = squeeze(t);

                % fit to sinusoid with a specified period
                ft = fittype(@(a0, A, phi, t) a0 + A*cos(omega*t+phi), 'independent', 't');

                % fit
                qfit = fit(tmon', q', ft, 'startpoint', [1 1 1], 'lower', [0 0 0]);
                tfit = fit(tmon', t', ft, 'startpoint', [1 1 1], 'lower', [0 0 0]);

                % compute effective heat capacity
                tphi = wrapTo2Pi(tfit.phi);
                qphi = wrapTo2Pi(qfit.phi);
                tA = abs(tfit.A);
                qA = abs(qfit.A);
                ceff = abs(sin(tphi - qphi))*(qA/tA)/omega;
                deff = ceff/(par.rho*par.cw);
                disp(sprintf('Ceff = %g (equiv to %g m of water) averaged for lat = %g to %g', ceff, deff, lat(1), lat(end)))

                figure(); clf; hold all; box on;
                plot(tmon, q, '-k')
                plot(teval, feval(qfit, teval), '-r')
                xlabel('time (s)'); ylabel('$Q$ (W m$^{-2}$)');
                axis('tight');
                % make_title_type_lat(type, lat(1), lat(end), titlepar);
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                set(gca, 'fontsize', par.fs, 'xminortick', 'on', 'xlim', [-inf inf])
                print(sprintf('%s/q_fitcheck.png', folder), '-dpng', '-r300');
                close;

                figure(); clf; hold all; box on;
                plot(tmon, t, '-k')
                plot(teval, feval(tfit, teval), '-r')
                xlabel('time (s)'); ylabel('$T_s$ (K)');
                axis('tight');
                % make_title_type_lat(type, lat(1), lat(end), titlepar);
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                set(gca, 'fontsize', par.fs, 'xminortick', 'on', 'xlim', [-inf inf])
                print(sprintf('%s/t_fitcheck.png', folder), '-dpng', '-r300');
                close;

            end % lat bounds

        end % land/ocean
    end % framework

end
