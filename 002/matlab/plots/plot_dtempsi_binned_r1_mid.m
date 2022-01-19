function plot_ga_malr_binned_r1_lat(type, par)
    titlepar = par;
    titlepar.gcm.clim = 'RCP8.5$-$historical';

    
    make_dirs(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/flux_z.mat', prefix_proc)); % load lat x mon RCAE_ALT data
    load(sprintf('%s/dtempsi_mon_lat.mat', prefix_proc)); % load lat x mon RCAE_ALT data
    load(sprintf('%s/dmasi_mon_lat.mat', prefix_proc)); % load lat x mon RCAE_ALT data
    
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
                clat_mon_lev = repmat(clat_mon, [1,1,size(dtempsi.(land),3)]);

                folder = sprintf('%s/dtempsi_binned_r1/%s/%s/0_lat_%g_%g', plotdir, fw, land, lat(1), lat(end));
                if ~exist(folder, 'dir'); mkdir(folder); end;
                
                % interpolate at select latitude
                r1 = flux_z.(land).res.(fw)./flux_z.(land).ra.(fw);

                r1 = interp1(grid.dim2.lat, r1, lat);
                r1 = nansum(clat_mon.*r1,1) / nansum(clat);
                r1_mon = squeeze(r1);

                dta = interp1(grid.dim3.lat, dtempsi.(land), lat);
                dta = nansum(clat_mon_lev.*dta,1) / nansum(clat);
                dta_mon = squeeze(dta);

                dma = interp1(grid.dim3.lat, dmasi.(land), lat);
                dma = nansum(clat_mon_lev.*dma,1) / nansum(clat);
                dma_mon = squeeze(dma);

                %%%%%%%% DTEMP %%%%%%%
                %% Months 1--6
                figure(); clf; hold all; box on;
                r1_bins_fine = linspace(par.r1_bins_45(1), par.r1_bins_45(end), 101);
                cmp = flip(parula(length(r1_bins_fine)-1));
                for mon = 1:6
                    [~, idxbin] = min(abs(r1(mon) - r1_bins_fine));
                    plot(dta_mon(mon,:), grid.dim3.si, 'color', cmp(idxbin,:), 'linewidth', 0.5);
                    text(dta_mon(mon,mon*3)-0.01*abs(dta_mon(mon,mon*3)), grid.dim3.si(mon*3), num2str(mon), 'color', cmp(idxbin, :), 'fontsize', 6);

                    vavg_dev_rce = nanmean(interp1(grid.dim3.si, dta_mon(mon,:), linspace(0.9, 0.3, 101)));
                end
                xlabel('$\Delta T$ (K)'); ylabel('$\sigma$ (unitless)');
                axis('tight');
                % caxis([min(par.r1_bins_45) max(par.r1_bins_45)]);
                % c = colorbar('ticks', [min(par.r1_bins_45):0.2:max(par.r1_bins_45)], 'ticklabels', strtrim(cellstr(num2str(flip([min(par.r1_bins_45):0.2:max(par.r1_bins_45)])', '%.1f'))'), 'ticklabelinterpreter', 'latex', 'ydir', 'reverse');
                % if ~strcmp(fw, 'mse_old')
                %     ylabel(c, '$R_1$ (unitless)', 'interpreter', 'latex');
                % else
                %     ylabel(c, '$R_1^*$ (unitless)', 'interpreter', 'latex');
                % end
                make_title_type_lat(type, lat(1), lat(end), titlepar);
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.3 1], 'xminortick', 'on', 'xlim', [2 6])
                print(sprintf('%s/dtempsi_r1_mon_1_6.png', folder), '-dpng', '-r300');
                close;
                
                %% Months 7--12
                figure(); clf; hold all; box on;
                r1_bins_fine = linspace(par.r1_bins_45(1), par.r1_bins_45(end), 101);
                cmp = flip(parula(length(r1_bins_fine)-1));
                %cmp = parula(length(par.r1_bins_45)-1);
                for mon = 7:12
                    [~, idxbin] = min(abs(r1(mon) - r1_bins_fine));
                    plot(dta_mon(mon,:), grid.dim3.si, 'color', cmp(idxbin,:), 'linewidth', 0.5);
                    text(dta_mon(mon,(mon-6)*10)-0.01*abs(dta_mon(mon,(mon-6)*10)), grid.dim3.si((mon-6)*10), num2str(mon), 'color', cmp(idxbin, :), 'fontsize', 6);

                    vavg_dev_rce = nanmean(interp1(grid.dim3.si, dta_mon(mon,:), linspace(0.9, 0.3, 101)));
                end
                xlabel('$\Delta T$ (K)'); ylabel('$\sigma$ (unitless)');
                axis('tight');
                % caxis([min(par.r1_bins_45) max(par.r1_bins_45)]);
                % c = colorbar('ticks', [min(par.r1_bins_45):0.2:max(par.r1_bins_45)], 'ticklabels', strtrim(cellstr(num2str(flip([min(par.r1_bins_45):0.2:max(par.r1_bins_45)])', '%.1f'))'), 'ticklabelinterpreter', 'latex', 'ydir', 'reverse');
                % if ~strcmp(fw, 'mse_old')
                %     ylabel(c, '$R_1$ (unitless)', 'interpreter', 'latex');
                % else
                %     ylabel(c, '$R_1^*$ (unitless)', 'interpreter', 'latex');
                % end
                make_title_type_lat(type, lat(1), lat(end), titlepar);
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.3 1], 'xminortick', 'on', 'xlim', [2 6])
                print(sprintf('%s/dtempsi_r1_mon_7_12.png', folder), '-dpng', '-r300');
                close;

                %%%%%%%% DTEMP WITH MA %%%%%%%
                %% Months 1--6
                figure(); clf; hold all; box on;
                r1_bins_fine = linspace(par.r1_bins_45(1), par.r1_bins_45(end), 101);
                cmp = flip(parula(length(r1_bins_fine)-1));
                for mon = 1:6
                    [~, idxbin] = min(abs(r1(mon) - r1_bins_fine));
                    plot(dma_mon(mon,:)-dta_mon(mon,:), grid.dim3.si, 'color', cmp(idxbin,:), 'linewidth', 0.5);
                    % text(dta_mon(mon,mon*3)-0.01*abs(dta_mon(mon,mon*3)), grid.dim3.si(mon*3), num2str(mon), 'color', cmp(idxbin, :), 'fontsize', 6);

                    vavg_dev_rce = nanmean(interp1(grid.dim3.si, dma_mon(mon,:)-dta_mon(mon,:), linspace(0.9, 0.3, 101)));
                end
                xlabel('$\Delta T_m - \Delta T$ (K)'); ylabel('$\sigma$ (unitless)');
                axis('tight');
                % caxis([min(par.r1_bins_45) max(par.r1_bins_45)]);
                % c = colorbar('ticks', [min(par.r1_bins_45):0.2:max(par.r1_bins_45)], 'ticklabels', strtrim(cellstr(num2str(flip([min(par.r1_bins_45):0.2:max(par.r1_bins_45)])', '%.1f'))'), 'ticklabelinterpreter', 'latex', 'ydir', 'reverse');
                % if ~strcmp(fw, 'mse_old')
                %     ylabel(c, '$R_1$ (unitless)', 'interpreter', 'latex');
                % else
                %     ylabel(c, '$R_1^*$ (unitless)', 'interpreter', 'latex');
                % end
                make_title_type_lat(type, lat(1), lat(end), titlepar);
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.3 1], 'xminortick', 'on', 'xlim', [-3 3])
                print(sprintf('%s/ddevsi_r1_mon_1_6.png', folder), '-dpng', '-r300');
                close;
                
                %% Months 7--12
                figure(); clf; hold all; box on;
                r1_bins_fine = linspace(par.r1_bins_45(1), par.r1_bins_45(end), 101);
                cmp = flip(parula(length(r1_bins_fine)-1));
                %cmp = parula(length(par.r1_bins_45)-1);
                for mon = 7:12
                    [~, idxbin] = min(abs(r1(mon) - r1_bins_fine));
                    plot(dma_mon(mon,:)-dta_mon(mon,:), grid.dim3.si, 'color', cmp(idxbin,:), 'linewidth', 0.5);
                    % text(dta_mon(mon,(mon-6)*10)-0.01*abs(dta_mon(mon,(mon-6)*10)), grid.dim3.si((mon-6)*10), num2str(mon), 'color', cmp(idxbin, :), 'fontsize', 6);

                    vavg_dev_rce = nanmean(interp1(grid.dim3.si, dta_mon(mon,:), linspace(0.9, 0.3, 101)));
                end
                xlabel('$\Delta T_m - \Delta T$ (K)'); ylabel('$\sigma$ (unitless)');
                axis('tight');
                % caxis([min(par.r1_bins_45) max(par.r1_bins_45)]);
                % c = colorbar('ticks', [min(par.r1_bins_45):0.2:max(par.r1_bins_45)], 'ticklabels', strtrim(cellstr(num2str(flip([min(par.r1_bins_45):0.2:max(par.r1_bins_45)])', '%.1f'))'), 'ticklabelinterpreter', 'latex', 'ydir', 'reverse');
                % if ~strcmp(fw, 'mse_old')
                %     ylabel(c, '$R_1$ (unitless)', 'interpreter', 'latex');
                % else
                %     ylabel(c, '$R_1^*$ (unitless)', 'interpreter', 'latex');
                % end
                make_title_type_lat(type, lat(1), lat(end), titlepar);
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.3 1], 'xminortick', 'on', 'xlim', [-3 3])
                print(sprintf('%s/ddevsi_r1_mon_7_12.png', folder), '-dpng', '-r300');
                close;

                %%%%%%%% NORM DTEMP WITH MA %%%%%%%
                %% Months 1--6
                figure(); clf; hold all; box on;
                r1_bins_fine = linspace(par.r1_bins_45(1), par.r1_bins_45(end), 101);
                cmp = flip(parula(length(r1_bins_fine)-1));
                dt_norm = 100 * ( dma_mon-dta_mon ) ./ dma_mon;
                for mon = 1:6
                    [~, idxbin] = min(abs(r1(mon) - r1_bins_fine));
                    plot(dt_norm(mon,:), grid.dim3.si, 'color', cmp(idxbin,:), 'linewidth', 0.5);
                    % text(dta_mon(mon,mon*3)-0.01*abs(dta_mon(mon,mon*3)), grid.dim3.si(mon*3), num2str(mon), 'color', cmp(idxbin, :), 'fontsize', 6);

                    vavg_dev_rce = nanmean(interp1(grid.dim3.si, dma_mon(mon,:)-dta_mon(mon,:), linspace(0.9, 0.3, 101)));
                end
                xlabel('$(\Delta T_m - \Delta T)/\Delta T_m$ (\%)'); ylabel('$\sigma$ (unitless)');
                axis('tight');
                % caxis([min(par.r1_bins_45) max(par.r1_bins_45)]);
                % c = colorbar('ticks', [min(par.r1_bins_45):0.2:max(par.r1_bins_45)], 'ticklabels', strtrim(cellstr(num2str(flip([min(par.r1_bins_45):0.2:max(par.r1_bins_45)])', '%.1f'))'), 'ticklabelinterpreter', 'latex', 'ydir', 'reverse');
                % if ~strcmp(fw, 'mse_old')
                %     ylabel(c, '$R_1$ (unitless)', 'interpreter', 'latex');
                % else
                %     ylabel(c, '$R_1^*$ (unitless)', 'interpreter', 'latex');
                % end
                make_title_type_lat(type, lat(1), lat(end), titlepar);
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.3 1], 'xminortick', 'on', 'xlim', [-50 50])
                print(sprintf('%s/ddevnormsi_r1_mon_1_6.png', folder), '-dpng', '-r300');
                close;
                
                %% Months 7--12
                figure(); clf; hold all; box on;
                r1_bins_fine = linspace(par.r1_bins_45(1), par.r1_bins_45(end), 101);
                cmp = flip(parula(length(r1_bins_fine)-1));
                %cmp = parula(length(par.r1_bins_45)-1);
                for mon = 7:12
                    [~, idxbin] = min(abs(r1(mon) - r1_bins_fine));
                    plot(dt_norm(mon,:), grid.dim3.si, 'color', cmp(idxbin,:), 'linewidth', 0.5);
                    % text(dta_mon(mon,(mon-6)*10)-0.01*abs(dta_mon(mon,(mon-6)*10)), grid.dim3.si((mon-6)*10), num2str(mon), 'color', cmp(idxbin, :), 'fontsize', 6);

                    vavg_dev_rce = nanmean(interp1(grid.dim3.si, dta_mon(mon,:), linspace(0.9, 0.3, 101)));
                end
                xlabel('$(\Delta T_m - \Delta T)/\Delta T_m$ (\%)'); ylabel('$\sigma$ (unitless)');
                axis('tight');
                % caxis([min(par.r1_bins_45) max(par.r1_bins_45)]);
                % c = colorbar('ticks', [min(par.r1_bins_45):0.2:max(par.r1_bins_45)], 'ticklabels', strtrim(cellstr(num2str(flip([min(par.r1_bins_45):0.2:max(par.r1_bins_45)])', '%.1f'))'), 'ticklabelinterpreter', 'latex', 'ydir', 'reverse');
                % if ~strcmp(fw, 'mse_old')
                %     ylabel(c, '$R_1$ (unitless)', 'interpreter', 'latex');
                % else
                %     ylabel(c, '$R_1^*$ (unitless)', 'interpreter', 'latex');
                % end
                make_title_type_lat(type, lat(1), lat(end), titlepar);
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.3 1], 'xminortick', 'on', 'xlim', [-50 50])
                print(sprintf('%s/ddevnormsi_r1_mon_7_12.png', folder), '-dpng', '-r300');
                close;

                %%%%%%%% NORMED DTEMP %%%%%%%
                %% Months 1--6
                figure(); clf; hold all; box on;
                r1_bins_fine = linspace(par.r1_bins_45(1), par.r1_bins_45(end), 101);
                cmp = flip(parula(length(r1_bins_fine)-1));
                for mon = 1:6
                    [~, idxbin] = min(abs(r1(mon) - r1_bins_fine));
                    plot(dta_mon(mon,:)/dta_mon(mon,1), grid.dim3.si, 'color', cmp(idxbin,:), 'linewidth', 0.5);
                    text(dta_mon(mon,mon*10)-0.01*abs(dta_mon(mon,mon*10)), grid.dim3.si(mon*10), num2str(mon), 'color', cmp(idxbin, :), 'fontsize', 6);

                    vavg_dev_rce = nanmean(interp1(grid.dim3.si, dta_mon(mon,:), linspace(0.9, 0.3, 101)));
                end
                xlabel('$\Delta T$ (K)'); ylabel('$\sigma$ (unitless)');
                axis('tight');
                % caxis([min(par.r1_bins_45) max(par.r1_bins_45)]);
                % c = colorbar('ticks', [min(par.r1_bins_45):0.2:max(par.r1_bins_45)], 'ticklabels', strtrim(cellstr(num2str(flip([min(par.r1_bins_45):0.2:max(par.r1_bins_45)])', '%.1f'))'), 'ticklabelinterpreter', 'latex', 'ydir', 'reverse');
                % if ~strcmp(fw, 'mse_old')
                %     ylabel(c, '$R_1$ (unitless)', 'interpreter', 'latex');
                % else
                %     ylabel(c, '$R_1^*$ (unitless)', 'interpreter', 'latex');
                % end
                make_title_type_lat(type, lat(1), lat(end), titlepar);
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.3 1], 'xminortick', 'on', 'xlim', [3 6])
                print(sprintf('%s/dtempsi_norm_r1_mon_1_6.png', folder), '-dpng', '-r300');
                close;
                
                %% Months 7--12
                figure(); clf; hold all; box on;
                r1_bins_fine = linspace(par.r1_bins_45(1), par.r1_bins_45(end), 101);
                cmp = flip(parula(length(r1_bins_fine)-1));
                %cmp = parula(length(par.r1_bins_45)-1);
                for mon = 7:12
                    [~, idxbin] = min(abs(r1(mon) - r1_bins_fine));
                    plot(dta_mon(mon,:)/dta_mon(mon,1), grid.dim3.si, 'color', cmp(idxbin,:), 'linewidth', 0.5);
                    text(dta_mon(mon,(mon-6)*10)-0.01*abs(dta_mon(mon,(mon-6)*10)), grid.dim3.si((mon-6)*10), num2str(mon), 'color', cmp(idxbin, :), 'fontsize', 6);

                    vavg_dev_rce = nanmean(interp1(grid.dim3.si, dta_mon(mon,:), linspace(0.9, 0.3, 101)));
                end
                xlabel('$\Delta T$ (K)'); ylabel('$\sigma$ (unitless)');
                axis('tight');
                % caxis([min(par.r1_bins_45) max(par.r1_bins_45)]);
                % c = colorbar('ticks', [min(par.r1_bins_45):0.2:max(par.r1_bins_45)], 'ticklabels', strtrim(cellstr(num2str(flip([min(par.r1_bins_45):0.2:max(par.r1_bins_45)])', '%.1f'))'), 'ticklabelinterpreter', 'latex', 'ydir', 'reverse');
                % if ~strcmp(fw, 'mse_old')
                %     ylabel(c, '$R_1$ (unitless)', 'interpreter', 'latex');
                % else
                %     ylabel(c, '$R_1^*$ (unitless)', 'interpreter', 'latex');
                % end
                make_title_type_lat(type, lat(1), lat(end), titlepar);
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.3 1], 'xminortick', 'on', 'xlim', [3 6])
                print(sprintf('%s/dtempsi_norm_r1_mon_7_12.png', folder), '-dpng', '-r300');
                close;

            end % lat bounds

        end % land/ocean
    end % framework

end
