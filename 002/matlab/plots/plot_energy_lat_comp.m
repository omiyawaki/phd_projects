function plot_energy_lat_comp(type, par)
% latitude vs energy flux line plots, comparable to Hartmann (2016)
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);

    tmp = load(sprintf('%s/flux_zt.mat', prefix_proc)); flux = tmp.flux_zt; flux_std = tmp.flux_zt_std; lat = tmp.lat;
    load(sprintf('%s/vh.mat', prefix_proc));
    load(sprintf('%s/vh_mon.mat', prefix_proc));
    
    make_dirs(type, par)
    par.lat_interp = 'native';
    
    [flux_era5c, vh_era5c, vh_mon_era5c, lat_era5c, plotdir_era5c] = load_flux('era5c', par); % load data
    [flux_merra, vh_merra, vh_mon_merra, lat_merra, plotdir_merra] = load_flux('merra2', par); % load data
    [flux_jra55, vh_jra55, vh_mon_jra55, lat_jra55, plotdir_jra55] = load_flux('jra55', par); % load data
    
    for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
        fw = 'mse_old';

        % for l = {'lo', 'l', 'o'}; land = l{1};
        for l = {'lo'}; land = l{1};
            if strcmp(land, 'lo'); land_text = '';
            elseif strcmp(land, 'l'); land_text = 'Land';
            elseif strcmp(land, 'o'); land_text = 'Ocean';
            end
            
            [lh, sh] = rename_stf(type, flux, land, time);
            [lh_std, sh_std] = rename_stf(type, flux_std, land, time);
            [lh_era5c, sh_era5c] = rename_stf('era5c', flux_era5c, land, time);
            [lh_merra, sh_merra] = rename_stf('merra2', flux_merra, land, time);
            [lh_jra55, sh_jra55] = rename_stf('jra55', flux_jra55, land, time);
            
            ra = flux.(land).(time).ra.(fw);
            ra_std = flux_std.(land).(time).ra.(fw);
            ra_era5c = flux_era5c.(land).(time).ra.(fw);
            ra_merra = flux_merra.(land).(time).ra.(fw);
            ra_jra55 = flux_jra55.(land).(time).ra.(fw);

            r1z = flux.(land).(time).r1z.(fw);
            r1z_std = flux_std.(land).(time).r1z.(fw);
            r1z_era5c = flux_era5c.(land).(time).res.(fw)./flux_era5c.(land).(time).ra.(fw);
            r1z_merra = flux_merra.(land).(time).res.(fw)./flux_merra.(land).(time).ra.(fw);
            r1z_jra55 = flux_jra55.(land).(time).res.(fw)./flux_jra55.(land).(time).ra.(fw);

            figure(); clf; hold all; box on;
            ylim_lo = -0.6;
            ylim_up = 1.5;
            rcemax = par.ep;
            vertices = [-90 ylim_lo; 90 ylim_lo; 90 rcemax; -90 rcemax];
            patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
            raemin = par.ga;
            vertices = [-90 raemin; 90 raemin; 90 ylim_up; -90 ylim_up];
            patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
            cmip = plot(lat, r1z, 'k');
            lat2 = [lat, fliplr(lat)];
            inbtw = [r1z+r1z_std, fliplr(r1z-r1z_std)];
            fill(lat2, inbtw, 'k', 'facealpha', 0.2, 'edgealpha', 0.2);
            era5 = plot(lat_era5c, r1z_era5c, 'color', par.blue);
            merra = plot(lat_merra, r1z_merra, 'color', par.orange);
            jra55 = plot(lat_jra55, r1z_jra55, 'color', par.green);
            title(sprintf('%s', upper(time)));
            xlabel('latitude (deg)');
            ylabel('$R_1$ (unitless)');
            legend([cmip, era5, merra, jra55], 'CMIP5', 'ERA5', 'MERRA2', 'JRA55', 'location', 'eastoutside');
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide)
            set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on', 'ylim', [ylim_lo ylim_up])
            print(sprintf('%s/energy-flux-comp/%s/%s/r1z', plotdir, land, time), '-dpng', '-r300');
            close;

            figure(); clf; hold all; box on;
            line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
            cmip=plot(lat, sh, 'color', 'k');
            lat_err = [lat, fliplr(lat)];
            sh_err = [sh+sh_std, fliplr(sh-sh_std)];
            fill(lat_err, sh_err, 'k', 'facealpha', 0.2, 'edgealpha', 0.2);
            era5=plot(lat_era5c, sh_era5c, 'color', par.blue);
            merra=plot(lat_merra, sh_merra, 'color', par.orange);
            jra55=plot(lat_jra55, sh_jra55, 'color', par.green);
            title(sprintf('%s', upper(time)));
            xlabel('latitude (deg)'); ylabel('SH (Wm$^{-2}$)');
            legend([cmip era5 merra jra55], 'CMIP5', 'ERA5', 'MERRA2', 'JRA55', 'location', 'eastoutside');
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide)
            set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
            if strcmp(fw, 'mse') & strcmp(time, 'ann'); set(gca, 'ylim', [-inf inf]); end
            print(sprintf('%s/energy-flux-comp/%s/%s/sensible', plotdir, land, time), '-dpng', '-r300');
            close;

            figure(); clf; hold all; box on;
            line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
            cmip=plot(lat, lh, 'color', 'k');
            lat_err = [lat, fliplr(lat)];
            lh_err = [lh+lh_std, fliplr(lh-lh_std)];
            fill(lat_err, lh_err, 'k', 'facealpha', 0.2, 'edgealpha', 0.2);
            era5=plot(lat_era5c, lh_era5c, 'color', par.blue);
            merra=plot(lat_merra, lh_merra, 'color', par.orange);
            jra55=plot(lat_jra55, lh_jra55, 'color', par.green);
            title(sprintf('%s', upper(time)));
            xlabel('latitude (deg)'); ylabel('LH (Wm$^{-2}$)');
            legend([cmip era5 merra jra55], 'CMIP5', 'ERA5', 'MERRA2', 'JRA55', 'location', 'eastoutside');
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide)
            set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
            if strcmp(fw, 'mse') & strcmp(time, 'ann'); set(gca, 'ylim', [-inf inf]); end
            print(sprintf('%s/energy-flux-comp/%s/%s/latent', plotdir, land, time), '-dpng', '-r300');
            close;

            figure(); clf; hold all; box on;
            line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
            cmip=plot(lat, ra, 'color', 'k');
            lat_err = [lat, fliplr(lat)];
            ra_err = [ra+ra_std, fliplr(ra-ra_std)];
            fill(lat_err, ra_err, 'k', 'facealpha', 0.2, 'edgealpha', 0.2);
            era5=plot(lat_era5c, ra_era5c, 'color', par.blue);
            merra=plot(lat_merra, ra_merra, 'color', par.orange);
            jra55=plot(lat_jra55, ra_jra55, 'color', par.green);
            title(sprintf('%s', upper(time)));
            xlabel('latitude (deg)'); ylabel('$R_a$ (Wm$^{-2}$)');
            legend([cmip era5 merra jra55], 'CMIP5', 'ERA5', 'MERRA2', 'JRA55', 'location', 'eastoutside');
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_verywide)
            set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
            if strcmp(fw, 'mse') & strcmp(time, 'ann'); set(gca, 'ylim', [-inf inf]); end
            print(sprintf('%s/energy-flux-comp/%s/%s/ra', plotdir, land, time), '-dpng', '-r300');
            close;
        
        end % land

    end % time

end
