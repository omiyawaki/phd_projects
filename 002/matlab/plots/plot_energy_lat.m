function plot_energy_lat(type, par)
% latitude vs energy flux line plots, comparable to Hartmann (2016)
    [flux, vh, vh_mon, lat, plotdir] = load_flux(type, par); % load data
    make_dirs(type, par)

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
                    plot(lat, flux.(land).(time).lw, 'color', par.gray); text(0, 0.85*interp1(lat,flux.(land).(time).lw,0), '\boldmath{$\mathrm{LW}$}', 'color', par.gray);
                    if contains(fw, 'db')
                        plot(lat,flux.(land).(time).res.(fw), 'color', par.maroon); text(-42, 1.2*interp1(lat,flux.(land).(time).res.(fw),-42), '\boldmath{$\nabla\cdot F_m - \mathrm{SW}$}', 'color', par.maroon);
                    else
                        plot(lat,flux.(land).(time).res.(fw), 'color', par.maroon); text(-42, 1.2*interp1(lat,flux.(land).(time).res.(fw),-42), '\boldmath{$\nabla\cdot F_m - \mathrm{SW}$}', 'color', par.maroon);
                    end
                    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
                        if contains(fw, 'db')
                            plot(lat, flux.(land).(time).stf.(fw), 'color', par.orange); text(20, 1.3*interp1(lat, flux.(land).(time).stf.(fw), 20)-25, sprintf('\\boldmath{$\\mathrm{LH+SH}$}'));
                            plot(lat, flux.(land).(time).stf.(fw), '--', 'color', par.blue);
                        elseif strcmp(fw, 'div')
                            plot(lat, flux.(land).(time).stf.(fw), 'color', par.orange); text(20, 1.3*interp1(lat, flux.(land).(time).stf.(fw), 20)-25, sprintf('\\boldmath{$\\mathrm{LH+SH}$}'));
                            plot(lat, flux.(land).(time).stf.(fw), '--', 'color', par.blue);
                        else
                            if strcmp(fw, 'mse'); plot(lat, -flux.(land).(time).slhf, 'color', par.blue); text(20, interp1(lat, -1.2*flux.(land).(time).slhf,20), '\boldmath{$\mathrm{LH}$}', 'color', par.blue);
                            elseif strcmp(fw, 'dse'); plot(lat, par.L*(flux.(land).(time).cp+flux.(land).(time).lsp), '--', 'color', par.blue); text(15, 1.5*interp1(lat, par.L*(flux.(land).(time).cp+flux.(land).(time).lsp),15), '\boldmath{$LP$}', 'color', par.blue); end
                            plot(lat, -flux.(land).(time).sshf, 'color', par.orange); text(80, interp1(lat, flux.(land).(time).sshf, 80)-25, '\boldmath{$\mathrm{SH}$}', 'color', par.orange);
                        end
                    elseif strcmp(type, 'merra2')
                        if contains(fw, 'mse'); plot(lat, flux.(land).(time).EFLUX, 'color', par.blue); text(20, 1.2*interp1(lat, flux.(land).(time).hfls,20), '\boldmath{$\mathrm{LH}$}', 'color', par.blue);
                        elseif strcmp(fw, 'dse'); plot(lat, par.L*flux.(land).(time).PRECTOT, '--', 'color', par.blue); text(15, 1.5*interp1(lat, par.L*flux.(land).(time).pr,15), '\boldmath{$LP$}', 'color', par.blue); end
                        plot(lat, flux.(land).(time).HFLUX, 'color', par.orange); text(80, interp1(lat, flux.(land).(time).hfss, 80)-25, '\boldmath{$\mathrm{SH}$}', 'color', par.orange);
                    elseif any(strcmp(type, {'gcm', 'jra55'}))
                        if contains(fw, 'mse'); plot(lat, flux.(land).(time).hfls, 'color', par.blue); text(20, 1.2*interp1(lat, flux.(land).(time).hfls,20), '\boldmath{$\mathrm{LH}$}', 'color', par.blue);
                        elseif strcmp(fw, 'dse'); plot(lat, par.L*flux.(land).(time).pr, '--', 'color', par.blue); text(15, 1.5*interp1(lat, par.L*flux.(land).(time).pr,15), '\boldmath{$LP$}', 'color', par.blue); end
                        plot(lat, flux.(land).(time).hfss, 'color', par.orange); text(80, interp1(lat, flux.(land).(time).hfss, 80)-25, '\boldmath{$\mathrm{SH}$}', 'color', par.orange);
                    elseif strcmp(type, 'echam')
                        if strcmp(fw, 'mse'); plot(lat, -flux.(land).(time).ahfl, 'color', par.blue); text(20, 1.2*interp1(lat, -flux.(land).(time).ahfl,20), '\boldmath{$\mathrm{LH}$}', 'color', par.blue);
                        elseif strcmp(fw, 'dse'); plot(lat, par.L*(flux.(land).(time).aprc+flux.(land).(time).aprl), '--', 'color', par.blue); text(15, 1.5*interp1(lat, par.L*(flux.(land).(time).aprc+flux.(land).(time).aprl),15), '\boldmath{$LP$}', 'color', par.blue);
                        end
                        plot(lat, -flux.(land).(time).ahfs, 'color', par.orange); text(80, interp1(lat, -flux.(land).(time).ahfs, 80)-25, '\boldmath{$\mathrm{SH}$}', 'color', par.orange);
                    end
                    make_title_type_time(type, time, par);
                    xlabel('latitude (deg)'); ylabel('energy flux (Wm$^{-2}$)');
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
                    if strcmp(fw, 'mse') & strcmp(time, 'ann'); set(gca, 'ylim', [-inf inf]); end
                    print(sprintf('%s/energy-flux/%s/%s/%s-all', plotdir, land, time, fw), '-dpng', '-r300');
                    close;
                else
                    figure(); clf; hold all; box on;
                    line([-90 90], [0 0], 'linewidth', 0.5, 'color', 'k');
                    plot(lat, flux.(land).(time).ra.(fw), 'color', par.gray); text(0, 0.75*interp1(lat,flux.(land).(time).ra.(fw),0), '\boldmath{$R_a$}', 'color', par.gray);
                    if strcmp(fw, 'dse'); plot(lat,flux.(land).(time).res.dse, '--', 'color', par.maroon); text(-30, 2*interp1(lat,flux.(land).(time).res.dse,-30), '\boldmath{$\nabla\cdot F_s$}', 'color', par.maroon);
                    else; plot(lat,flux.(land).(time).res.(fw), 'color', par.maroon); text(-42, 2*interp1(lat,flux.(land).(time).res.(fw),-42), '\boldmath{$\nabla\cdot F_m$}', 'color', par.maroon); end;
                    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
                        if contains(fw, 'db')
                            plot(lat, flux.(land).(time).stf.(fw), 'color', par.orange); text(20, 1.3*interp1(lat, flux.(land).(time).stf.(fw), 20)-25, sprintf('\\boldmath{$\\mathrm{LH+SH}$}'));
                            plot(lat, flux.(land).(time).stf.(fw), '--', 'color', par.blue);
                        elseif contains(fw, 'div')
                            plot(lat, flux.(land).(time).stf.(fw), 'color', par.orange); text(20, 1.3*interp1(lat, flux.(land).(time).stf.(fw), 20)-25, sprintf('\\boldmath{$\\mathrm{LH+SH}$}'));
                            plot(lat, flux.(land).(time).stf.(fw), '--', 'color', par.blue);
                        else
                            if strcmp(fw, 'mse'); plot(lat, -flux.(land).(time).slhf, 'color', par.blue); text(20, interp1(lat, -1.2*flux.(land).(time).slhf,20), '\boldmath{$\mathrm{LH}$}', 'color', par.blue);
                            elseif strcmp(fw, 'dse'); plot(lat, par.L*(flux.(land).(time).cp+flux.(land).(time).lsp), '--', 'color', par.blue); text(15, 1.5*interp1(lat, par.L*(flux.(land).(time).cp+flux.(land).(time).lsp),15), '\boldmath{$LP$}', 'color', par.blue);
                            end
                            plot(lat, -flux.(land).(time).sshf, 'color', par.orange); text(80, interp1(lat, flux.(land).(time).sshf, 80)-25, '\boldmath{$\mathrm{SH}$}', 'color', par.orange);
                        end
                    elseif strcmp(type, 'merra2')
                        if strcmp(fw, 'mse'); plot(lat, flux.(land).(time).EFLUX, 'color', par.blue); text(20, 1.2*interp1(lat, flux.(land).(time).EFLUX,20), '\boldmath{$\mathrm{LH}$}', 'color', par.blue);
                        elseif strcmp(fw, 'dse'); plot(lat, par.L*flux.(land).(time).PRECTOT, '--', 'color', par.blue); text(15, 1.5*interp1(lat, par.L*flux.(land).(time).PRECTOT,15), '\boldmath{$LP$}', 'color', par.blue);
                        end
                        plot(lat, flux.(land).(time).HFLUX, 'color', par.orange); text(80, interp1(lat, flux.(land).(time).HFLUX, 80)-25, '\boldmath{$\mathrm{SH}$}', 'color', par.orange);
                    elseif any(strcmp(type, {'gcm', 'jra55'}))
                        if strcmp(fw, 'mse'); plot(lat, flux.(land).(time).hfls, 'color', par.blue); text(20, 1.2*interp1(lat, flux.(land).(time).hfls,20), '\boldmath{$\mathrm{LH}$}', 'color', par.blue);
                        elseif strcmp(fw, 'dse'); plot(lat, par.L*flux.(land).(time).pr, '--', 'color', par.blue); text(15, 1.5*interp1(lat, par.L*flux.(land).(time).pr,15), '\boldmath{$LP$}', 'color', par.blue);
                        end
                        plot(lat, flux.(land).(time).hfss, 'color', par.orange); text(80, interp1(lat, flux.(land).(time).hfss, 80)-25, '\boldmath{$\mathrm{SH}$}', 'color', par.orange);
                    elseif strcmp(type, 'echam')
                        if strcmp(fw, 'mse'); plot(lat, -flux.(land).(time).ahfl, 'color', par.blue); text(20, 1.2*interp1(lat, -flux.(land).(time).ahfl,20), '\boldmath{$\mathrm{LH}$}', 'color', par.blue);
                        elseif strcmp(fw, 'dse'); plot(lat, par.L*(flux.(land).(time).aprc+flux.(land).(time).aprl), '--', 'color', par.blue); text(15, 1.5*interp1(lat, par.L*(flux.(land).(time).aprc+flux.(land).(time).aprl),15), '\boldmath{$LP$}', 'color', par.blue);
                        end
                        plot(lat, -flux.(land).(time).ahfs, 'color', par.orange); text(80, interp1(lat, -flux.(land).(time).ahfs, 80)-25, '\boldmath{$\mathrm{SH}$}', 'color', par.orange);
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
                    plot(lat, flux.(land).(time).ra.(fw), 'color', par.gray); text(0, 0.75*interp1(lat,flux.(land).(time).ra.(fw),0), '\boldmath{$R_a$}', 'color', par.gray);
                    plot(lat, flux.(land).(time).sw, 'color', par.yellow); text(0, 0.75*interp1(lat,flux.(land).(time).sw,0), '\boldmath{$\mathrm{SWABS}$}', 'color', par.yellow);
                    plot(lat, flux.(land).(time).lw, 'color', par.green); text(0, 0.75*interp1(lat,flux.(land).(time).lw,0), '\boldmath{$\mathrm{LWCOOL}$}', 'color', par.green);
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
                    plot(lat, flux.(land).(time).lw, 'color', par.green); text(0, 0.75*interp1(lat,flux.(land).(time).lw,0), '\boldmath{$\mathrm{LWCOOL}$}', 'color', par.green);
                    plot(lat, flux.(land).(time).olr, '--', 'color', par.green); text(0, 0.95*interp1(lat,flux.(land).(time).olr,0), '\boldmath{$\mathrm{OLR}$}', 'color', par.green);
                    plot(lat, flux.(land).(time).lwsfc, ':', 'color', par.green); text(0, 0.75*interp1(lat,flux.(land).(time).lwsfc,0), '\boldmath{$\mathrm{LW_{SFC}}$}', 'color', par.green);
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
                    if strcmp(fw, 'dse'); plot(lat,flux.(land).(time).r1.dse, '--k');
                    else; plot(lat,flux.(land).(time).r1.(fw), '-k'); end
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
                    ylim_lo = min(flux.(land).(time).r1.(fw)); if isnan(ylim_lo)|ylim_lo==0; ylim_lo = -1; end;
                    ylim_up = max(flux.(land).(time).r1.(fw)); if isnan(ylim_up)|ylim_up==0; ylim_up = 1; end
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
                    if strcmp(fw, 'dse'); plot(lat,flux.(land).(time).r1.dse, '--k');
                    else; plot(lat,flux.(land).(time).r1.(fw), '-k'); end
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
                    r1z=flux.(land).(time).res.(fw)./flux.(land).(time).ra.(fw);
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
