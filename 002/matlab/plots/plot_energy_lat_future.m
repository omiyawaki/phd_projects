function plot_energy_lat_future(type, par)
% latitude vs energy flux_zt line plots, comparable to Hartmann (2016)
    
    make_dirs(type, par)
    
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);

    load(sprintf('%s/flux_zt.mat', prefix_proc));

    par2 = par;
    par2.gcm.clim = 'rcp85';
    par2.gcm.yr_span = '207001-209912';
    prefix_proc2 = make_prefix_proc(type, par2);
    rcp=load(sprintf('%s/flux_zt.mat', prefix_proc2));

    for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
        f_vec = assign_fw(type, par);
        for f = f_vec; fw = f{1};
            % for l = {'lo', 'l', 'o'}; land = l{1};
            for l = {'lo'}; land = l{1};
                if strcmp(land, 'lo'); land_text = '';
                elseif strcmp(land, 'l'); land_text = 'Land';
                elseif strcmp(land, 'o'); land_text = 'Ocean';
                end

                % R1Z
                figure(); clf; hold all; box on;
                [r1z r1zu r1zl] = deal(nan(size(lat)));
                % if strcmp(type, 'rea')
                %     r1z=flux_zt.(land).(time).r1z.(fw);
                %     r1z_min = flux_zt_min.(land).(time).r1z.(fw);
                %     r1z_max = flux_zt_max.(land).(time).r1z.(fw);
                %     r1zu = r1z_min;
                %     r1zl = r1z_max;
                % elseif (any(strcmp(type, {'gcm'})) & strcmp(par.model, 'mmm'))
                %     r1z=flux_zt.(land).(time).r1z.(fw);
                %     r1z_25 = flux_zt_25.(land).(time).r1z.(fw);
                %     r1z_75 = flux_zt_75.(land).(time).r1z.(fw);
                %     r1zu = r1z_25;
                %     r1zl = r1z_75;
                % else
                    r1z=flux_zt.(land).(time).res.(fw)./flux_zt.(land).(time).ra.(fw);
                    r1z2=rcp.flux_zt.(land).(time).res.(fw)./rcp.flux_zt.(land).(time).ra.(fw);
                % end

                % ylim_lo = min(r1z); if isnan(ylim_lo)|ylim_lo==0; ylim_lo = -1; end;
                % ylim_up = max(r1z); if isnan(ylim_up)|ylim_up==0; ylim_up = 1; end
                ylim_lo = -0.6;
                ylim_up = 1.8;
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
                % if strcmp(type, 'rea') | (any(strcmp(type, {'gcm'})) & strcmp(par.model, 'mmm'))
                %     lat2 = [lat, fliplr(lat)];
                %     inbtw = [r1zu, fliplr(r1zl)];
                %     fill(lat2, inbtw, 'k', 'facealpha', 0.2, 'edgealpha', 0.2);
                % end
                h1=plot(lat,r1z, '-k');
                h2=plot(lat,r1z2, '--k');
                legend([h1 h2],'Historical', 'RCP8.5', 'location', 'north');
                % make_title_type_time(type, time, par);
                title(sprintf('CMIP5, %s', upper(time)));
                xlabel('Latitude (deg)');
                if strcmp(fw, 'mse2'); ylabel('$R_1^*$ (unitless)');
                else ylabel('$R_1$ (unitless)'); end
                axis('tight');
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
                set(gca, 'fontsize', par.fs, 'xtick', [-90:30:90], 'xminortick', 'on', 'yminortick', 'on', 'ylim', [ylim_lo ylim_up])
                print(sprintf('%s/energy-flux/%s/%s/%s-r1z-rcp', plotdir, land, time, fw), '-dpng', '-r300');
                if par.make_tikz
                    matlab2tikz(sprintf('%s/energy-flux/%s/%s/%s-r1z-rcp.tex', plotdir, land, time, fw));
                end
                close;

            end % land
        end % fw
    end % time

end
