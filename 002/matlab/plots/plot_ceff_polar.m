function plot_ceff_polar(type, par)
    
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

    lat_list = [80 -80];
    
    for f = {'mse_old'}; fw = f{1};
        for l = par.land_list; land = l{1};
            if strcmp(land, 'lo'); land_text = 'Land + Ocean';
            elseif strcmp(land, 'l'); land_text = 'Land';
            elseif strcmp(land, 'o'); land_text = 'Ocean';
            end

            for lb = 1:length(lat_list); par.lat_bound = lat_list(lb);

                [lat, clat, clat_mon, ~] = make_polar_lat(par);

                folder = sprintf('%s/ceff/%s/%s/0_lat_%g_%g', plotdir, fw, land, lat(1), lat(end));
                if ~exist(folder, 'dir'); mkdir(folder); end;
                
                % interpolate at select latitude
                q = interp1(grid.dim2.lat, -flux_z.lo.swsfc, lat);
                q = nansum(clat_mon.*q,1) / nansum(clat);
                q_mon = squeeze(q);

                f = interp1(grid.dim2.lat, -flux_z.lo.fsfc.(fw), lat);
                f = nansum(clat_mon.*f,1) / nansum(clat);
                f_mon = squeeze(f);

                t = interp1(grid.dim2.lat, tas, lat);
                t = nansum(clat_mon.*t,1) / nansum(clat);
                t_mon = squeeze(t);

                % fit to sinusoid with a specified period
                ft = fittype(@(a0, A, phi, t) a0 + A*cos(omega*t+phi), 'independent', 't');

                % fit
                qfit = fit(tmon', q', ft, 'startpoint', [1 1 1], 'lower', [0 0 0]);
                ffit = fit(tmon', f', ft, 'startpoint', [1 1 1], 'lower', [0 0 0]);
                tfit = fit(tmon', t', ft, 'startpoint', [1 1 1], 'lower', [0 0 0]);

                % compute effective heat capacity
                tphi = wrapTo2Pi(tfit.phi);
                fphi = wrapTo2Pi(ffit.phi);
                qphi = wrapTo2Pi(qfit.phi);
                tA = abs(tfit.A);
                fA = abs(ffit.A);
                qA = abs(qfit.A);

                % ceff = abs(sin(tphi - fphi))*(fA/tA)/omega;
                ceff = abs(sin(tphi - qphi))*(qA/tA)/omega;

                deff = ceff/(par.rho*par.cw);
                disp(sprintf('Ceff = %g (equiv to %g m of water) averaged for lat = %g to %g', ceff, deff, lat(1), lat(end)))

                % save data to file
                save(sprintf('%s/ceff.mat', folder), 'ceff', 'deff');

                figure(); clf; hold all; box on;
                plot(tmon, f, '-k')
                plot(teval, feval(ffit, teval), '-r')
                xlabel('time (s)'); ylabel('$F_{\mathrm{SFC}}$ (W m$^{-2}$)');
                axis('tight');
                % make_title_type_lat(type, lat(1), lat(end), titlepar);
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                set(gca, 'fontsize', par.fs, 'xminortick', 'on', 'xlim', [-inf inf])
                print(sprintf('%s/f_fitcheck.png', folder), '-dpng', '-r300');
                close;

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

                % compare with ECHAM simulations
                figure(); clf; hold all; box on;
                line([10, 50], ceff*[1 1], 'color', 'k', 'linestyle', '--')
                d.rp000126 = 50;
                d.rp000148 = 45;
                d.rp000134 = 40;
                d.rp000146 = 35;
                d.rp000130 = 30;
                d.rp000144 = 25;
                d.rp000132 = 20;
                d.rp000140 = 15;
                d.rp000124 = 10;
                for i = 1:length(par.echam_clims); echam_clim = par.echam_clims{i};
                    % load data
                    parecham = par;
                    parecham.lat_interp = 'native';
                    parecham.echam.clim = echam_clim;
                    plotdirecham = make_plotdir('echam', parecham);
                    folderecham = sprintf('%s/ceff/%s/%s/0_lat_%g_%g', plotdirecham, fw, land, lat(1), lat(end));
                    load(sprintf('%s/ceff.mat', folderecham));

                    plot(d.(echam_clim), ceff, '*k');
                end
                xlabel('$d$ (m)'); ylabel('$C$ (J m$^{-2}$ K$^{-1}$)');
                axis('tight');
                % make_title_type_lat(type, lat(1), lat(end), titlepar);
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                set(gca, 'fontsize', par.fs, 'xminortick', 'on', 'xlim', [-inf inf])
                print(sprintf('%s/ceff_echamcomp.png', folder), '-dpng', '-r300');
                close;

            end % lat bounds

        end % land/ocean
    end % framework

end
