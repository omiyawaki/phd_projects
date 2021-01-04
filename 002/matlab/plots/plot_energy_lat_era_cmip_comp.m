function plot_energy_lat_era_cmip_comp(type, par)
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
                    elseif strcmp(type, 'gcm')
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
                    if any(strcmp(type, {'era5', 'era5c', 'erai'}));
                        if strcmp(fw, 'db13s'); title(sprintf('%s, %s, %s, %s', upper(type), 'DB13*', upper(time), land_text));
                        else; title(sprintf('%s, %s, %s, %s', upper(type), upper(fw), upper(time), land_text)); end;
                    elseif strcmp(type, 'merra2'); title(sprintf('%s, %s, %s, %s', upper(type), upper(fw), upper(time), land_text));
                    elseif strcmp(type, 'gcm');
                        if contains(par.model, 'mmm')
                            title(sprintf('CMIP5 %s, %s, %s', par.gcm.clim, upper(fw), upper(time)));
                        else
                            title(sprintf('%s, %s, %s, %s', par.model, upper(fw), upper(time), land_text));
                        end
                    elseif strcmp(type, 'echam'); title(sprintf('%s, %s, %s, %s', upper(type), upper(fw), upper(time), land_text)); end
                    xlabel('latitude (deg)'); ylabel('energy flux (Wm$^{-2}$)');
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on')
                    if strcmp(fw, 'mse') & strcmp(time, 'ann'); set(gca, 'ylim', [-inf inf]); end
                    print(sprintf('%s/energy-flux/%s/%s/%s-all', plotdir, land, time, fw), '-dpng', '-r300');
                    close;
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
            title(sprintf('Northward %s Transport (PW)', upper(fw)));
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
            set(gca, 'fontsize', par.fs, 'xtick', [1:12], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on')
            print(sprintf('%s/transport/%s/all/%s', plotdir, land, fw), '-dpng', '-r300');
            close;
        end
    end

end
