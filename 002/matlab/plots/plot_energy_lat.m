function plot_energy_lat(type, par)
% latitude vs energy flux_zt line plots, comparable to Hartmann (2016)
    
    make_dirs(type, par)
    
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);

    load(sprintf('%s/flux_zt.mat', prefix_proc));
    load(sprintf('%s/vh.mat', prefix_proc));
    load(sprintf('%s/vh_mon.mat', prefix_proc));

    for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
        f_vec = assign_fw(type, par);
        for f = f_vec; fw = f{1};
            % for l = {'lo', 'l', 'o'}; land = l{1};
            for l = {'lo'}; land = l{1};
                if strcmp(land, 'lo'); land_text = '';
                elseif strcmp(land, 'l'); land_text = 'Land';
                elseif strcmp(land, 'o'); land_text = 'Ocean';
                end

                if strcmp(fw, 'mse2')
                    figure(); clf; hold all; box on;
                    line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
                    plot(lat, flux_zt.(land).(time).lw, 'color', par.gray); text(0, 0.85*interp1(lat,flux_zt.(land).(time).lw,0), '\boldmath{$\mathrm{LW}$}', 'color', par.gray);
                    if contains(fw, 'db')
                        plot(lat,flux_zt.(land).(time).res.(fw), 'color', par.maroon); text(-42, 1.2*interp1(lat,flux_zt.(land).(time).res.(fw),-42), '\boldmath{$\nabla\cdot F_m - \mathrm{SW}$}', 'color', par.maroon);
                    else
                        plot(lat,flux_zt.(land).(time).res.(fw), 'color', par.maroon); text(-42, 1.2*interp1(lat,flux_zt.(land).(time).res.(fw),-42), '\boldmath{$\nabla\cdot F_m - \mathrm{SW}$}', 'color', par.maroon);
                    end
                    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
                        if contains(fw, 'db')
                            plot(lat, flux_zt.(land).(time).stf.(fw), 'color', par.orange); text(20, 1.3*interp1(lat, flux_zt.(land).(time).stf.(fw), 20)-25, sprintf('\\boldmath{$\\mathrm{LH+SH}$}'));
                            plot(lat, flux_zt.(land).(time).stf.(fw), '--', 'color', par.blue);
                        elseif strcmp(fw, 'div')
                            plot(lat, flux_zt.(land).(time).stf.(fw), 'color', par.orange); text(20, 1.3*interp1(lat, flux_zt.(land).(time).stf.(fw), 20)-25, sprintf('\\boldmath{$\\mathrm{LH+SH}$}'));
                            plot(lat, flux_zt.(land).(time).stf.(fw), '--', 'color', par.blue);
                        else
                            if strcmp(fw, 'mse'); plot(lat, -flux_zt.(land).(time).slhf, 'color', par.blue); text(20, interp1(lat, -1.2*flux_zt.(land).(time).slhf,20), '\boldmath{$\mathrm{LH}$}', 'color', par.blue);
                            elseif strcmp(fw, 'dse'); plot(lat, par.L*(flux_zt.(land).(time).cp+flux_zt.(land).(time).lsp), '--', 'color', par.blue); text(15, 1.5*interp1(lat, par.L*(flux_zt.(land).(time).cp+flux_zt.(land).(time).lsp),15), '\boldmath{$LP$}', 'color', par.blue); end
                            plot(lat, -flux_zt.(land).(time).sshf, 'color', par.orange); text(80, interp1(lat, flux_zt.(land).(time).sshf, 80)-25, '\boldmath{$\mathrm{SH}$}', 'color', par.orange);
                        end
                    elseif strcmp(type, 'hahn')
                        if contains(fw, 'mse'); plot(lat, flux_zt.(land).(time).LHFLX, 'color', par.blue); text(20, 1.2*interp1(lat, flux_zt.(land).(time).LHFLX,20), '\boldmath{$\mathrm{LH}$}', 'color', par.blue);
                        elseif strcmp(fw, 'dse'); plot(lat, par.L*(flux_zt.(land).(time).PRECC+flux_zt.(land).(time).PRECL+flux_zt.(land).(time).PRECSC+flux_zt.(land).(time).PRECSL), '--', 'color', par.blue); text(15, 1.5*interp1(lat, par.L*(flux_zt.(land).(time).PRECC+flux_zt.(land).(time).PRECL+flux_zt.(land).(time).PRECSC+flux_zt.(land).(time).PRECSL),15), '\boldmath{$LP$}', 'color', par.blue); end
                        plot(lat, flux_zt.(land).(time).SHFLX, 'color', par.orange); text(80, interp1(lat, flux_zt.(land).(time).SHFLX, 80)-25, '\boldmath{$\mathrm{SH}$}', 'color', par.orange);
                    elseif strcmp(type, 'merra2')
                        if contains(fw, 'mse'); plot(lat, flux_zt.(land).(time).EFLUX, 'color', par.blue); text(20, 1.2*interp1(lat, flux_zt.(land).(time).EFLUX,20), '\boldmath{$\mathrm{LH}$}', 'color', par.blue);
                        elseif strcmp(fw, 'dse'); plot(lat, par.L*flux_zt.(land).(time).PRECTOT, '--', 'color', par.blue); text(15, 1.5*interp1(lat, par.L*flux_zt.(land).(time).PRECTOT,15), '\boldmath{$LP$}', 'color', par.blue); end
                        plot(lat, flux_zt.(land).(time).HFLUX, 'color', par.orange); text(80, interp1(lat, flux_zt.(land).(time).HFLUX, 80)-25, '\boldmath{$\mathrm{SH}$}', 'color', par.orange);
                    elseif any(strcmp(type, {'gcm', 'jra55', 'rea'}))
                        if contains(fw, 'mse'); plot(lat, flux_zt.(land).(time).hfls, 'color', par.blue); text(20, 1.2*interp1(lat, flux_zt.(land).(time).hfls,20), '\boldmath{$\mathrm{LH}$}', 'color', par.blue);
                        elseif strcmp(fw, 'dse'); plot(lat, par.L*flux_zt.(land).(time).pr, '--', 'color', par.blue); text(15, 1.5*interp1(lat, par.L*flux_zt.(land).(time).pr,15), '\boldmath{$LP$}', 'color', par.blue); end
                        plot(lat, flux_zt.(land).(time).hfss, 'color', par.orange); text(80, interp1(lat, flux_zt.(land).(time).hfss, 80)-25, '\boldmath{$\mathrm{SH}$}', 'color', par.orange);
                    elseif strcmp(type, 'echam')
                        if strcmp(fw, 'mse'); plot(lat, -flux_zt.(land).(time).ahfl, 'color', par.blue); text(20, 1.2*interp1(lat, -flux_zt.(land).(time).ahfl,20), '\boldmath{$\mathrm{LH}$}', 'color', par.blue);
                        elseif strcmp(fw, 'dse'); plot(lat, par.L*(flux_zt.(land).(time).aprc+flux_zt.(land).(time).aprl), '--', 'color', par.blue); text(15, 1.5*interp1(lat, par.L*(flux_zt.(land).(time).aprc+flux_zt.(land).(time).aprl),15), '\boldmath{$LP$}', 'color', par.blue);
                        end
                        plot(lat, -flux_zt.(land).(time).ahfs, 'color', par.orange); text(80, interp1(lat, -flux_zt.(land).(time).ahfs, 80)-25, '\boldmath{$\mathrm{SH}$}', 'color', par.orange);
                    end
                    make_title_type_time(type, time, par);
                    xlabel('latitude (deg)'); ylabel('energy flux_zt (Wm$^{-2}$)');
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
                    if strcmp(fw, 'mse') & strcmp(time, 'ann'); set(gca, 'ylim', [-inf inf]); end
                    print(sprintf('%s/energy-flux_zt/%s/%s/%s-all', plotdir, land, time, fw), '-dpng', '-r300');
                    close;
                else
                    figure(); clf; hold all; box on;
                    line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
                    plot(lat, flux_zt.(land).(time).ra.(fw), 'color', par.gray); text(0, 0.75*interp1(lat,flux_zt.(land).(time).ra.(fw),0), '\boldmath{$R_a$}', 'color', par.gray);
                    % if contains(type, 'era')
                    %     plot(lat, flux_zt.(land).(time).tend, 'color', 'k'); text(0, 0.75*interp1(lat,flux_zt.(land).(time).tend,0), '\boldmath{$\partial_t h$}', 'color', 'k');
                    % end
                    if strcmp(fw, 'dse'); plot(lat,flux_zt.(land).(time).res.dse, '--', 'color', par.maroon); text(-30, 2*interp1(lat,flux_zt.(land).(time).res.dse,-30), '\boldmath{$\nabla\cdot F_s$}', 'color', par.maroon);
                    else; plot(lat,flux_zt.(land).(time).res.(fw), 'color', par.maroon); text(-42, 2*interp1(lat,flux_zt.(land).(time).res.(fw),-42), '\boldmath{$\nabla\cdot F_m$}', 'color', par.maroon); end;
                    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
                        if contains(fw, 'db')
                            plot(lat, flux_zt.(land).(time).stf.(fw), 'color', par.orange); text(20, 1.3*interp1(lat, flux_zt.(land).(time).stf.(fw), 20)-25, sprintf('\\boldmath{$\\mathrm{LH+SH}$}'));
                            plot(lat, flux_zt.(land).(time).stf.(fw), '--', 'color', par.blue);
                        elseif contains(fw, 'div')
                            plot(lat, flux_zt.(land).(time).stf.(fw), 'color', par.orange); text(20, 1.3*interp1(lat, flux_zt.(land).(time).stf.(fw), 20)-25, sprintf('\\boldmath{$\\mathrm{LH+SH}$}'));
                            plot(lat, flux_zt.(land).(time).stf.(fw), '--', 'color', par.blue);
                        else
                            if strcmp(fw, 'mse'); plot(lat, -flux_zt.(land).(time).slhf, 'color', par.blue); text(20, interp1(lat, -1.2*flux_zt.(land).(time).slhf,20), '\boldmath{$\mathrm{LH}$}', 'color', par.blue);
                            elseif strcmp(fw, 'dse'); plot(lat, par.L*(flux_zt.(land).(time).cp+flux_zt.(land).(time).lsp), '--', 'color', par.blue); text(15, 1.5*interp1(lat, par.L*(flux_zt.(land).(time).cp+flux_zt.(land).(time).lsp),15), '\boldmath{$LP$}', 'color', par.blue);
                            end
                            plot(lat, -flux_zt.(land).(time).sshf, 'color', par.orange); text(80, interp1(lat, flux_zt.(land).(time).sshf, 80)-25, '\boldmath{$\mathrm{SH}$}', 'color', par.orange);
                        end
                    elseif strcmp(type, 'hahn')
                        if strcmp(fw, 'mse'); plot(lat, flux_zt.(land).(time).LHFLX, 'color', par.blue); text(20, 1.2*interp1(lat, flux_zt.(land).(time).LHFLX,20), '\boldmath{$\mathrm{LH}$}', 'color', par.blue);
                        elseif strcmp(fw, 'dse'); plot(lat, par.L*(flux_zt.(land).(time).PRECC+flux_zt.(land).(time).PRECL+flux_zt.(land).(time).PRECSC+flux_zt.(land).(time).PRECSL), '--', 'color', par.blue); text(15, 1.5*interp1(lat, par.L*(flux_zt.(land).(time).PRECC+flux_zt.(land).(time).PRECL+flux_zt.(land).(time).PRECSC+flux_zt.(land).(time).PRECSL),15), '\boldmath{$LP$}', 'color', par.blue);
                        end
                        plot(lat, flux_zt.(land).(time).SHFLX, 'color', par.orange); text(80, interp1(lat, flux_zt.(land).(time).SHFLX, 80)-25, '\boldmath{$\mathrm{SH}$}', 'color', par.orange);
                    elseif strcmp(type, 'merra2')
                        if strcmp(fw, 'mse'); plot(lat, flux_zt.(land).(time).EFLUX, 'color', par.blue); text(20, 1.2*interp1(lat, flux_zt.(land).(time).EFLUX,20), '\boldmath{$\mathrm{LH}$}', 'color', par.blue);
                        elseif strcmp(fw, 'dse'); plot(lat, par.L*flux_zt.(land).(time).PRECTOT, '--', 'color', par.blue); text(15, 1.5*interp1(lat, par.L*flux_zt.(land).(time).PRECTOT,15), '\boldmath{$LP$}', 'color', par.blue);
                        end
                        plot(lat, flux_zt.(land).(time).HFLUX, 'color', par.orange); text(80, interp1(lat, flux_zt.(land).(time).HFLUX, 80)-25, '\boldmath{$\mathrm{SH}$}', 'color', par.orange);
                    elseif any(strcmp(type, {'gcm', 'jra55', 'rea'}))
                        if strcmp(fw, 'mse'); plot(lat, flux_zt.(land).(time).hfls, 'color', par.blue); text(20, 1.2*interp1(lat, flux_zt.(land).(time).hfls,20), '\boldmath{$\mathrm{LH}$}', 'color', par.blue);
                        elseif strcmp(fw, 'dse'); plot(lat, par.L*flux_zt.(land).(time).pr, '--', 'color', par.blue); text(15, 1.5*interp1(lat, par.L*flux_zt.(land).(time).pr,15), '\boldmath{$LP$}', 'color', par.blue);
                        end
                        plot(lat, flux_zt.(land).(time).hfss, 'color', par.orange); text(80, interp1(lat, flux_zt.(land).(time).hfss, 80)-25, '\boldmath{$\mathrm{SH}$}', 'color', par.orange);
                    elseif strcmp(type, 'echam')
                        if strcmp(fw, 'mse'); plot(lat, -flux_zt.(land).(time).ahfl, 'color', par.blue); text(20, 1.2*interp1(lat, -flux_zt.(land).(time).ahfl,20), '\boldmath{$\mathrm{LH}$}', 'color', par.blue);
                        elseif strcmp(fw, 'dse'); plot(lat, par.L*(flux_zt.(land).(time).aprc+flux_zt.(land).(time).aprl), '--', 'color', par.blue); text(15, 1.5*interp1(lat, par.L*(flux_zt.(land).(time).aprc+flux_zt.(land).(time).aprl),15), '\boldmath{$LP$}', 'color', par.blue);
                        end
                        plot(lat, -flux_zt.(land).(time).ahfs, 'color', par.orange); text(80, interp1(lat, -flux_zt.(land).(time).ahfs, 80)-25, '\boldmath{$\mathrm{SH}$}', 'color', par.orange);
                    end
                    make_title_type_time(type, time, par);
                    xlabel('latitude (deg)'); ylabel('energy flux (Wm$^{-2}$)');
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
                    if strcmp(fw, 'mse') & strcmp(time, 'ann'); set(gca, 'ylim', [-inf inf]); end
                    print(sprintf('%s/energy-flux/%s/%s/%s-all', plotdir, land, time, fw), '-dpng', '-r300');
                    close;
                end

                if strcmp(fw, 'mse')
                    figure(); clf; hold all; box on;
                    line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
                    plot(lat, flux_zt.(land).(time).ra.(fw), 'color', par.gray); text(0, 0.75*interp1(lat,flux_zt.(land).(time).ra.(fw),0), '\boldmath{$R_a$}', 'color', par.gray);
                    plot(lat, flux_zt.(land).(time).sw, 'color', par.yellow); text(0, 0.75*interp1(lat,flux_zt.(land).(time).sw,0), '\boldmath{$\mathrm{SWABS}$}', 'color', par.yellow);
                    plot(lat, flux_zt.(land).(time).lw, 'color', par.green); text(0, 0.75*interp1(lat,flux_zt.(land).(time).lw,0), '\boldmath{$\mathrm{LWCOOL}$}', 'color', par.green);
                    make_title_type_time(type, time, par);
                    xlabel('latitude (deg)'); ylabel('energy flux (Wm$^{-2}$)');
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
                    if strcmp(fw, 'mse') & strcmp(time, 'ann'); set(gca, 'ylim', [-inf inf]); end
                    print(sprintf('%s/energy-flux/%s/%s/ra', plotdir, land, time), '-dpng', '-r300');
                    close;

                    figure(); clf; hold all; box on;
                    line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
                    plot(lat, flux_zt.(land).(time).lw, 'color', par.green); text(0, 0.75*interp1(lat,flux_zt.(land).(time).lw,0), '\boldmath{$\mathrm{LWCOOL}$}', 'color', par.green);
                    plot(lat, flux_zt.(land).(time).olr, '--', 'color', par.green); text(0, 0.95*interp1(lat,flux_zt.(land).(time).olr,0), '\boldmath{$\mathrm{OLR}$}', 'color', par.green);
                    plot(lat, flux_zt.(land).(time).lwsfc, ':', 'color', par.green); text(0, 0.75*interp1(lat,flux_zt.(land).(time).lwsfc,0), '\boldmath{$\mathrm{LW_{SFC}}$}', 'color', par.green);
                    make_title_type_time(type, time, par);
                    xlabel('latitude (deg)'); ylabel('energy flux (Wm$^{-2}$)');
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
                    if strcmp(fw, 'mse') & strcmp(time, 'ann'); set(gca, 'ylim', [-inf inf]); end
                    print(sprintf('%s/energy-flux/%s/%s/lwcool', plotdir, land, time), '-dpng', '-r300');
                    close;

                end

                if strcmp(type, 'gcm') & strcmp(par.gcm.clim, 'hist-pi')
                    figure(); clf; hold all; box on;
                    line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
                    if strcmp(fw, 'dse'); plot(lat,flux_zt.(land).(time).r1.dse, '--k');
                    else; plot(lat,flux_zt.(land).(time).r1.(fw), '-k'); end
                    make_title_type_time(type, time, par);
                    xlabel('latitude (deg)');
                    ylabel('$R_1^{\mathrm{historical}}-R_1^{\mathrm{piControl}}$ (unitless)')
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
                    set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
                    print(sprintf('%s/energy-flux/%s/%s/%s-r1', plotdir, land, time, fw), '-dpng', '-r300');
                    close;

                else
                    figure(); clf; hold all; box on;
                    ylim_lo = min(flux_zt.(land).(time).r1.(fw)); if isnan(ylim_lo)|ylim_lo==0; ylim_lo = -1; end;
                    ylim_up = max(flux_zt.(land).(time).r1.(fw)); if isnan(ylim_up)|ylim_up==0; ylim_up = 1; end
                    rcemax = par.ep;
                    vertices = [-90 ylim_lo; 90 ylim_lo; 90 rcemax; -90 rcemax];
                    patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
                    raemin = par.ga;
                    vertices = [-90 raemin; 90 raemin; 90 ylim_up; -90 ylim_up];
                    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                    line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
                    % line([-90 90], [1 1]*par.ep, 'linewidth', 0.5, 'color', par.orange);
                    % if ~strcmp(fw, 'mse2'); line([-90 90], -[1 1]*par.ep, 'linewidth', 0.5, 'color', par.orange); end
                    % line([-90 90], [1 1]*par.ga, 'linewidth', 0.5, 'color', par.blue);
                    if strcmp(fw, 'dse'); plot(lat,flux_zt.(land).(time).r1.dse, '--k');
                    else; plot(lat,flux_zt.(land).(time).r1.(fw), '-k'); end
                    make_title_type_time(type, time, par);
                    xlabel('latitude (deg)');
                    if strcmp(fw, 'mse2'); ylabel('$R_1^*$ (unitless)');
                    else ylabel('$R_1$ (unitless)'); end
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
                    set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
                    print(sprintf('%s/energy-flux/%s/%s/%s-r1', plotdir, land, time, fw), '-dpng', '-r300');
                    close;

                    % R1Z
                    figure(); clf; hold all; box on;
                    [r1z r1zu r1zl] = deal(nan(size(lat)));
                    if strcmp(type, 'rea')
                        r1z=flux_zt.(land).(time).r1z.(fw);
                        r1z_min = flux_zt_min.(land).(time).r1z.(fw);
                        r1z_max = flux_zt_max.(land).(time).r1z.(fw);
                        r1zu = r1z_min;
                        r1zl = r1z_max;
                    elseif (any(strcmp(type, {'gcm'})) & strcmp(par.model, 'mmm'))
                        r1z=flux_zt.(land).(time).r1z.(fw);
                        r1z_25 = flux_zt_25.(land).(time).r1z.(fw);
                        r1z_75 = flux_zt_75.(land).(time).r1z.(fw);
                        r1zu = r1z_25;
                        r1zl = r1z_75;
                    else
                        r1z=flux_zt.(land).(time).res.(fw)./flux_zt.(land).(time).ra.(fw);
                    end
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
                    % line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
                    % line([-90 90], [1 1]*par.ep, 'linewidth', 0.5, 'color', par.orange);
                    % if ~strcmp(fw, 'mse2'); line([-90 90], -[1 1]*par.ep, 'linewidth', 0.5, 'color', par.orange); end
                    % line([-90 90], [1 1]*par.ga, 'linewidth', 0.5, 'color', par.blue);
                    if strcmp(type, 'rea') | (any(strcmp(type, {'gcm'})) & strcmp(par.model, 'mmm'))
                        lat2 = [lat, fliplr(lat)];
                        inbtw = [r1zu, fliplr(r1zl)];
                        fill(lat2, inbtw, 'k', 'facealpha', 0.2, 'edgealpha', 0.2);
                    end
                    if strcmp(fw, 'dse'); plot(lat,r1z, '--k');
                    else; plot(lat,r1z, '-k'); end
                    make_title_type_time(type, time, par);
                    xlabel('latitude (deg)');
                    if strcmp(fw, 'mse2'); ylabel('$R_1^*$ (unitless)');
                    else ylabel('$R_1$ (unitless)'); end
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
                    set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on', 'ylim', [-0.6 1.5])
                    print(sprintf('%s/energy-flux/%s/%s/%s-r1z', plotdir, land, time, fw), '-dpng', '-r300');
                    close;

                    % % R2Z
                    % figure(); clf; hold all; box on;
                    % [r2z r2zu r2zl] = deal(nan(size(lat)));
                    % if strcmp(type, 'rea') | (any(strcmp(type, {'gcm'})) & strcmp(par.model, 'mmm'))
                    %     r2z=flux_zt.(land).(time).r2z.(fw);
                    %     r2z_std = flux_zt_std.(land).(time).r2z.(fw);
                    %     r2zu=r2z + r2z_std;
                    %     r2zl=r2z - r2z_std;
                    % else
                    %     r2z=flux_zt.(land).(time).stf.(fw)./flux_zt.(land).(time).ra.(fw);
                    % end
                    % % ylim_lo = min(r2z); if isnan(ylim_lo)|ylim_lo==0; ylim_lo = -1; end;
                    % % ylim_up = max(r2z); if isnan(ylim_up)|ylim_up==0; ylim_up = 1; end
                    % ylim_lo = -0.6-1;
                    % ylim_up = 1.5-1;
                    % rcemax = par.ep-1;
                    % vertices = [-90 ylim_lo; 90 ylim_lo; 90 rcemax; -90 rcemax];
                    % patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
                    % raemin = par.ga-1;
                    % vertices = [-90 raemin; 90 raemin; 90 ylim_up; -90 ylim_up];
                    % patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
                    % % line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
                    % % line([-90 90], [1 1]*par.ep, 'linewidth', 0.5, 'color', par.orange);
                    % % if ~strcmp(fw, 'mse2'); line([-90 90], -[1 1]*par.ep, 'linewidth', 0.5, 'color', par.orange); end
                    % % line([-90 90], [1 1]*par.ga, 'linewidth', 0.5, 'color', par.blue);
                    % if strcmp(type, 'rea') | (any(strcmp(type, {'gcm'})) & strcmp(par.model, 'mmm'))
                    %     lat2 = [lat, fliplr(lat)];
                    %     inbtw = [r2zu, fliplr(r2zl)];
                    %     fill(lat2, inbtw, 'k', 'facealpha', 0.2, 'edgealpha', 0.2);
                    % end
                    % if strcmp(fw, 'dse'); plot(lat,r2z, '--k');
                    % else; plot(lat,r2z, '-k'); end
                    % make_title_type_time(type, time, par);
                    % xlabel('latitude (deg)');
                    % if strcmp(fw, 'mse2'); ylabel('$R_2^*$ (unitless)');
                    % else ylabel('$R_2$ (unitless)'); end
                    % axis('tight');
                    % set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
                    % set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on', 'ylim', [-0.6-1 1.5-1])
                    % print(sprintf('%s/energy-flux/%s/%s/%s-r2z', plotdir, land, time, fw), '-dpng', '-r300');
                    % close;

                end

            % northward M/DSE transport
                figure(); clf; hold all; box on;
                line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
                line([0 0], [min(vh.(land).(time).(fw)) max(vh.(land).(time).(fw))]*10^-15, 'linewidth', 0.5, 'color', 'k');
                if strcmp(fw, 'dse'); plot(lat, vh.(land).(time).(fw)*10^-15, '--', 'color', par.maroon);
                else; plot(lat, vh.(land).(time).(fw)*10^-15, 'color', par.maroon);
                end
                xlabel('latitude (deg)'); ylabel('PW')
                if strcmp(fw, 'db13s'); title(sprintf('Northward %s Transport, %s', 'DB13*', upper(time)));
                else
                    if any(strcmp(type, {'erai', 'era5', 'merra2', 'jra55'}))
                        if contains(fw, 'mse')
                            title(sprintf('%s, Northward MSE Transport, %s', upper(type), upper(time)));
                        elseif contains(fw, 'dse')
                            title(sprintf('%s, Northward DSE Transport, %s', upper(type), upper(time)));
                        end
                    elseif strcmp(type, 'era5c')
                        if contains(fw, 'mse')
                            title(sprintf('%s, Northward MSE Transport, %s', upper('era5'), upper(time)));
                        elseif contains(fw, 'dse')
                            title(sprintf('%s, Northward DSE Transport, %s', upper('era5'), upper(time)));
                        end
                    elseif strcmp(type, 'gcm')
                        if contains(par.model, 'mmm')
                            title(sprintf('CMIP5 %s, Northward %s Transport, %s', par.gcm.clim, upper(fw), upper(time)));
                        else
                            title(sprintf('%s, Northward %s Transport, %s', par.model, upper(fw), upper(time)));
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
        end % end mse/dse loop

    end % time

    for f = f_vec; fw = f{1};
        % for l = {'lo', 'l', 'o'}; land = l{1};
        for l = {'lo'}; land = l{1};
            % northward M/DSE transport, mon x lat
            [mesh_lat, mesh_mon] = meshgrid(1:12, lat);
            figure(); clf; hold all; box on;
            cmp = colCog(20);
            colormap(cmp);
            imagesc([1 12], [lat(1) lat(end)], vh_mon.(land).(fw)*1e-15);
            cb = colorbar('ticks', [-5:1:5], 'ticklabelinterpreter', 'latex');
            ylabel(cb, '$F_m$ (PW)', 'interpreter', 'latex');
            % [C, h] = contour(mesh_lat, mesh_mon, vh_mon.(land).(fw)*10^-15, -5:1:5);
            % clabel(C, h, [-4:2:4], 'fontsize', 6, 'interpreter', 'latex');
            caxis([-5 5]);
            xlabel('Month'); ylabel('latitude (deg)');
            if contains(fw, 'mse')
                title(sprintf('Northward MSE Transport (PW)'));
            elseif contains(fw, 'dse')
                title(sprintf('Northward MSE Transport (PW)'));
            end
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
            set(gca, 'fontsize', par.fs, 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on')
            print(sprintf('%s/transport/%s/all/%s', plotdir, land, fw), '-dpng', '-r300');
            close;
        end
    end

end
