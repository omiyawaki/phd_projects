function plot_energy_lat_comp(type, par)
% latitude vs energy flux line plots, comparable to Hartmann (2016)
    [flux, vh, vh_mon, lat, plotdir] = load_flux(type, par); % load data
    make_dirs(type, par)
    par.lat_interp = 'native';
    [flux_era5, vh_era5, vh_mon_era5, lat_era5, plotdir_era5] = load_flux('era5c', par); % load data
    [flux_jra55, vh_jra55, vh_mon_jra55, lat_jra55, plotdir_jra55] = load_flux('jra55', par); % load data

    for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
        fw = 'mse';

        % for l = {'lo', 'l', 'o'}; land = l{1};
        for l = {'lo'}; land = l{1};
        if strcmp(land, 'lo'); land_text = '';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end

        figure(); clf; hold all; box on;
        line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
        sh=plot(lat,interp1(lat_era5,-flux_era5.(land).(time).sshf,lat)-flux.(land).(time).hfss, '-', 'color', par.orange);
        lh=plot(lat,interp1(lat_era5,-flux_era5.(land).(time).slhf,lat)-flux.(land).(time).hfls, '-', 'color', par.blue);
        ra=plot(lat,interp1(lat_era5,flux_era5.(land).(time).ra.mse,lat)-flux.(land).(time).ra.mse, '-', 'color', par.gray);
        divfm=plot(lat,interp1(lat_era5,flux_era5.(land).(time).res.mse,lat)-flux.(land).(time).res.mse, '-', 'color', par.maroon);
        title(sprintf('$\\mathrm{ERA5 - CMIP5\\,mmm}$, %s', upper(time)));
        xlabel('latitude (deg)'); ylabel('Energy flux (Wm$^{-2}$)');
        legend([divfm, ra, lh, sh], '$\Delta(\nabla\cdot F_m)$', '$\Delta R_a$', '$\Delta$LH', '$\Delta$SH', 'location', 'eastoutside');
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
        set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
        if strcmp(fw, 'mse') & strcmp(time, 'ann'); set(gca, 'ylim', [-inf inf]); end
        print(sprintf('%s/energy-flux-comp/%s/%s/era5-diff-all', plotdir, land, time), '-dpng', '-r300');
        close;

        figure(); clf; hold all; box on;
        line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
        sh=plot(lat,interp1(lat_era5,-flux_era5.(land).(time).sshf,lat)-flux.(land).(time).hfss, '-', 'color', par.orange);
        lh=plot(lat,interp1(lat_era5,-flux_era5.(land).(time).slhf,lat)-flux.(land).(time).hfls, '-', 'color', par.blue);
        ra=plot(lat,interp1(lat_era5,flux_era5.(land).(time).ra.mse,lat)-flux.(land).(time).ra.mse, '-', 'color', par.gray);
        divfm=plot(lat,interp1(lat_era5,flux_era5.(land).(time).res.mse_ac,lat)-flux.(land).(time).res.mse, '-', 'color', par.maroon);
        title(sprintf('$\\mathrm{ERA5AC - CMIP5\\,mmm}$, %s', upper(time)));
        xlabel('latitude (deg)'); ylabel('Energy flux (Wm$^{-2}$)');
        legend([divfm, ra, lh, sh], '$\Delta(\nabla\cdot F_m)$', '$\Delta R_a$', '$\Delta$LH', '$\Delta$SH', 'location', 'eastoutside');
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
        set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
        if strcmp(fw, 'mse') & strcmp(time, 'ann'); set(gca, 'ylim', [-inf inf]); end
        print(sprintf('%s/energy-flux-comp/%s/%s/era5ac-diff-all', plotdir, land, time), '-dpng', '-r300');
        close;

        figure(); clf; hold all; box on;
        line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
        sh=plot(lat,interp1(lat_jra55,flux_jra55.(land).(time).hfss,lat)-flux.(land).(time).hfss, '-', 'color', par.orange);
        lh=plot(lat,interp1(lat_jra55,flux_jra55.(land).(time).hfls,lat)-flux.(land).(time).hfls, '-', 'color', par.blue);
        ra=plot(lat,interp1(lat_jra55,flux_jra55.(land).(time).ra.mse,lat)-flux.(land).(time).ra.mse, '-', 'color', par.gray);
        divfm=plot(lat,interp1(lat_jra55,flux_jra55.(land).(time).res.mse,lat)-flux.(land).(time).res.mse, '-', 'color', par.maroon);
        title(sprintf('$\\mathrm{JRA55 - CMIP5\\,mmm}$, %s', upper(time)));
        xlabel('latitude (deg)'); ylabel('Energy flux (Wm$^{-2}$)');
        legend([divfm, ra, lh, sh], '$\Delta(\nabla\cdot F_m)$', '$\Delta R_a$', '$\Delta$LH', '$\Delta$SH', 'location', 'eastoutside');
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
        set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
        if strcmp(fw, 'mse') & strcmp(time, 'ann'); set(gca, 'ylim', [-inf inf]); end
        print(sprintf('%s/energy-flux-comp/%s/%s/jra55-diff-all', plotdir, land, time), '-dpng', '-r300');
        close;

        figure(); clf; hold all; box on;
        % ylim_lo = min(r1z); if isnan(ylim_lo)|ylim_lo==0; ylim_lo = -1; end;
        % ylim_up = max(r1z); if isnan(ylim_up)|ylim_up==0; ylim_up = 1; end
        ylim_lo = -0.6;
        ylim_up = 1.5;
        rcemax = par.ep;
        vertices = [-90 ylim_lo; 90 ylim_lo; 90 rcemax; -90 rcemax];
        patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
        raemin = par.ga;
        vertices = [-90 raemin; 90 raemin; 90 ylim_up; -90 ylim_up];
        patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
        line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
        cmip=plot(lat,flux.(land).(time).res.mse./flux.(land).(time).ra.mse, 'k');
        era=plot(lat_era5,flux_era5.(land).(time).res.mse./flux_era5.(land).(time).ra.mse, '--k');
        xlabel('latitude (deg)'); ylabel('$R_1$ (unitless)');
        legend([cmip,era],'CMIP5 mmm','ERA5', 'location', 'southoutside');
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
        set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
        if strcmp(fw, 'mse') & strcmp(time, 'ann'); set(gca, 'ylim', [-inf inf]); end
        print(sprintf('%s/energy-flux-comp/%s/%s/era5_r1z', plotdir, land, time), '-dpng', '-r300');
        close;

        figure(); clf; hold all; box on;
        line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
        cmip=plot(lat,flux.(land).(time).res.mse, 'color', par.maroon);
        era=plot(lat_era5,flux_era5.(land).(time).res.mse, '--', 'color', par.maroon);
        xlabel('latitude (deg)'); ylabel('$\nabla\cdot F_m$ (Wm$^{-2}$)');
        legend([cmip,era],'CMIP5 mmm','ERA5', 'location', 'southoutside');
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
        set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
        if strcmp(fw, 'mse') & strcmp(time, 'ann'); set(gca, 'ylim', [-inf inf]); end
        print(sprintf('%s/energy-flux-comp/%s/%s/era5_divfm', plotdir, land, time), '-dpng', '-r300');
        close;

        figure(); clf; hold all; box on;
        line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
        cmip=plot(lat,flux.(land).(time).res.mse, 'color', par.maroon);
        era=plot(lat_era5,flux_era5.(land).(time).res.mse_ac, '--', 'color', par.maroon);
        xlabel('latitude (deg)'); ylabel('$\nabla\cdot F_m$ (Wm$^{-2}$)');
        legend([cmip,era],'CMIP5 mmm','ERA5AC', 'location', 'southoutside');
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
        set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
        if strcmp(fw, 'mse') & strcmp(time, 'ann'); set(gca, 'ylim', [-inf inf]); end
        print(sprintf('%s/energy-flux-comp/%s/%s/era5_divfm-ac', plotdir, land, time), '-dpng', '-r300');
        close;

        figure(); clf; hold all; box on;
        line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
        cmip=plot(lat,flux.(land).(time).res.mse, 'color', par.maroon);
        era=plot(lat_era5,flux_era5.(land).(time).res.mse_sc, '--', 'color', par.maroon);
        xlabel('latitude (deg)'); ylabel('$\nabla\cdot F_m$ (Wm$^{-2}$)');
        legend([cmip,era],'CMIP5 mmm','ERA5SC', 'location', 'southoutside');
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
        set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
        if strcmp(fw, 'mse') & strcmp(time, 'ann'); set(gca, 'ylim', [-inf inf]); end
        print(sprintf('%s/energy-flux-comp/%s/%s/era5_divfm-sc', plotdir, land, time), '-dpng', '-r300');
        close;

        figure(); clf; hold all; box on;
        line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
        cmip=plot(lat,flux.(land).(time).ra.mse, 'color', par.gray);
        era=plot(lat_era5,flux_era5.(land).(time).ra.mse, '--', 'color', par.gray);
        xlabel('latitude (deg)'); ylabel('$R_a$ (Wm$^{-2}$)');
        legend([cmip,era],'CMIP5 mmm','ERA5', 'location', 'southoutside');
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
        set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
        if strcmp(fw, 'mse') & strcmp(time, 'ann'); set(gca, 'ylim', [-inf inf]); end
        print(sprintf('%s/energy-flux-comp/%s/%s/era5_ra', plotdir, land, time), '-dpng', '-r300');
        close;

        figure(); clf; hold all; box on;
        line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
        cmip=plot(lat,flux.(land).(time).hfls, 'color', par.blue);
        era=plot(lat_era5,-flux_era5.(land).(time).slhf, '--', 'color', par.blue);
        xlabel('latitude (deg)'); ylabel('LH (Wm$^{-2}$)');
        legend([cmip,era],'CMIP5 mmm','ERA5', 'location', 'southoutside');
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
        set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
        if strcmp(fw, 'mse') & strcmp(time, 'ann'); set(gca, 'ylim', [-inf inf]); end
        print(sprintf('%s/energy-flux-comp/%s/%s/era5_lh', plotdir, land, time), '-dpng', '-r300');
        close;

        figure(); clf; hold all; box on;
        line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
        cmip=plot(lat,flux.(land).(time).hfss, 'color', par.orange);
        era=plot(lat_era5,-flux_era5.(land).(time).sshf, '--', 'color', par.orange);
        xlabel('latitude (deg)'); ylabel('SH (Wm$^{-2}$)');
        legend([cmip,era],'CMIP5 mmm','ERA5', 'location', 'southoutside');
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
        set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
        if strcmp(fw, 'mse') & strcmp(time, 'ann'); set(gca, 'ylim', [-inf inf]); end
        print(sprintf('%s/energy-flux-comp/%s/%s/era5_sh', plotdir, land, time), '-dpng', '-r300');
        close;

        % northward M/DSE transport
        figure(); clf; hold all; box on;
        line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
        line([0 0], [min(vh.(land).(time).mse) max(vh.(land).(time).mse)]*10^-15, 'linewidth', 0.5, 'color', 'k');
        if strcmp(fw, 'dse'); plot(lat, vh.(land).(time).mse*10^-15, '--', 'color', par.maroon);
        else; plot(lat, vh.(land).(time).mse*10^-15, 'color', par.maroon);
        end
        xlabel('latitude (deg)'); ylabel('PW')
        if strcmp(fw, 'db13s'); title(sprintf('Northward %s Transport, %s', 'DB13*', upper(time)));
        else
                if any(strcmp(type, {'erai', 'era5', 'era5c', 'merra2', 'jra55'}))
                title(sprintf('%s, Northward %s Transport, %s', upper(type), upper(mse), upper(time)));
                elseif strcmp(type, 'gcm')
                if contains(par.model, 'mmm')
                        title(sprintf('CMIP5 %s, Northward %s Transport, %s', par.gcm.clim, upper('mse'), upper(time)));
                else
                        title(sprintf('%s, Northward %s Transport, %s', par.model, upper('mse'), upper(time)));
                end
                end
        end;
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
        set(gca, 'fontsize', par.fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
        if strcmp(fw, 'dse'); set(gca, 'ytick', [-5:5]); end;
        print(sprintf('%s/transport/%s/%s/%s', plotdir, land, time, fw), '-dpng', '-r300');
        close;
        % if ~contains(fw, {'db'})
        % % MSE/DSE transport plotted together
        %     figure(); clf; hold all;
        %     line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
        %     line([0 0], [min([vh.(land).(time).mse; vh.(land).(time).dse; vh.(land).(time).mse-vh.(land).(time).dse]) max([vh.(land).(time).mse; vh.(land).(time).dse; vh.(land).(time).mse-vh.(land).(time).dse])]*10^-15, 'linewidth', 0.5, 'color', 'k');
        %     h_mse = plot(lat, vh.(land).(time).mse*10^-15, 'color', par.maroon);
        %     h_dse = plot(lat, vh.(land).(time).dse*10^-15, '--', 'color', par.maroon);
        %     h_lh = plot(lat, (vh.(land).(time).mse-vh.(land).(time).dse)*10^-15, ':', 'color', par.maroon);
        %     legend([h_mse h_dse h_lh], '$F_m$', '$F_s$', '$F_{m}-F_{s}$', 'location', 'eastoutside');
        %     xlabel('latitude (deg)'); ylabel('PW')
        %     title(sprintf('Northward Energy Transport, %s', upper(time)));
        %     axis('tight');
        %     set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
        %     set(gca, 'fontsize', par.fs, 'xlim', [-90 90], 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
        %     print(sprintf('%s/transport/%s/%s/all', plotdir, land, time), '-dpng', '-r300');
        %     close;
        % end
        end % land

    end % time

end