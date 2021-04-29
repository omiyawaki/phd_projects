function plot_ga_malr_binned_r1_lat(type, par)
    make_dirs(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/flux_z.mat', prefix_proc)); % load lat x mon RCAE_ALT data
    load(sprintf('%s/ga_diff_mon_lat.mat', prefix_proc)); % load lat x mon RCAE_ALT data
    load(sprintf('%s/ga_frac_mon_lat.mat', prefix_proc)); % load lat x mon RCAE_ALT data
    load(sprintf('%s/gad_frac_mon_lat.mat', prefix_proc)); % load lat x mon RCAE_ALT data

    lat_list = [-85 85 45 -45];
    
    for f = {'mse_old'}; fw = f{1};
        for l = par.land_list; land = l{1};
            if strcmp(land, 'lo'); land_text = 'Land + Ocean';
            elseif strcmp(land, 'l'); land_text = 'Land';
            elseif strcmp(land, 'o'); land_text = 'Ocean';
            end

            for lb = 1:length(lat_list); lat_eval = lat_list(lb);

                folder = sprintf('%s/ga_diff_binned_r1/%s/%s/0_lat_%g', plotdir, fw, land, lat_eval);
                if ~exist(folder, 'dir'); mkdir(folder); end;
                
                % interpolate at select latitude
                res = interp1(grid.dim2.lat, flux_z.(land).res.(fw), lat_eval);
                ra = interp1(grid.dim2.lat, flux_z.(land).ra.(fw), lat_eval);
                ga = 100 * interp1(grid.dim3.lat, ga_diff.(land), lat_eval); % multiply by 100 to express in percentage
                ga_mon = squeeze(ga);
                ga_fr = 100 * interp1(grid.dim3.lat, ga_frac.(land), lat_eval);
                ga_fr_mon = squeeze(ga_fr);
                gad_fr = 100 * interp1(grid.dim3.lat, gad_frac.(land), lat_eval);
                gad_fr_mon = squeeze(gad_fr);
            
                % reshape array to vector
                r1 = res(:)./ra(:);
                ga = reshape(ga, [], length(grid.dim3.si));
                ga_fr = reshape(ga_fr, [], length(grid.dim3.si));
                gad_fr = reshape(gad_fr, [], length(grid.dim3.si));

                id_bins = discretize(r1, par.r1_bins_45);

                [ga_area ga_fr_area gad_fr_area] = deal(nan([length(par.r1_bins_45)-1, length(grid.dim3.si)]));

                for bin = 1:length(par.r1_bins_45)-1
                    idx_bins = find(id_bins==bin);

                    ga_area(bin,:) = nanmean(ga(idx_bins,:),1);
                    ga_fr_area(bin,:) = nanmean(ga_fr(idx_bins,:),1);
                    gad_fr_area(bin,:) = nanmean(gad_fr(idx_bins,:),1);

                end

                % %%%%%%%% GA DIFF BINNED BY R1 %%%%%%%%
                % [~,idx09]=min(abs(par.r1_bins_45-0.9));
                % figure(); clf; hold all; box on;
                % line([0 0], [0 1], 'color', 'k', 'linewidth', 0.5);
                % cmp = flip(parula(length(par.r1_bins_45)-1));
                % %cmp = parula(length(par.r1_bins_45)-1);
                % for bin = 1:length(par.r1_bins_45)-1
                %     if bin==idx09
                %         plot(ga_area(bin,:), grid.dim3.si, 'color', cmp(bin,:), 'linewidth', 1.5);
                %     else
                %         plot(ga_area(bin,:), grid.dim3.si, 'color', cmp(bin,:), 'linewidth', 0.5);
                %     end
                % end
                % xlabel('$\Gamma_m - \Gamma$ (K km$^{-1}$)'); ylabel('$\sigma$ (unitless)');
                % axis('tight');
                % caxis([min(par.r1_bins_45) max(par.r1_bins_45)]);
                % c = colorbar('ticks', [min(par.r1_bins_45):0.2:max(par.r1_bins_45)], 'ticklabels', strtrim(cellstr(num2str(flip([min(par.r1_bins_45):0.2:max(par.r1_bins_45)])', '%.1f'))'), 'ticklabelinterpreter', 'latex', 'ydir', 'reverse');
                % if ~strcmp(fw, 'mse_old')
                %     ylabel(c, '$R_1$ (unitless)', 'interpreter', 'latex');
                % else
                %     ylabel(c, '$R_1^*$ (unitless)', 'interpreter', 'latex');
                % end
                % make_title_type_lat_pt(type, lat_eval, par);
                % set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                % set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.3 1], 'xminortick', 'on')
                % print(sprintf('%s/ga_diff_r1_all.png', folder), '-dpng', '-r300');
                % close;

                % %%%%%%%%%% GA DIFF %%%%%%%%%%
                % %% Months 1--6
                % figure(); clf; hold all; box on;
                % line([0 0], [0 1], 'color', 'k', 'linewidth', 0.5);
                % r1_bins_fine = linspace(par.r1_bins_45(1), par.r1_bins_45(end), 101);
                % cmp = flip(parula(length(r1_bins_fine)-1));
                % %cmp = parula(length(par.r1_bins_45)-1);
                % for mon = 1:6
                %     [~, idxbin] = min(abs(r1(mon) - r1_bins_fine));
                %     plot(ga_mon(mon,:), grid.dim3.si, 'color', cmp(idxbin,:), 'linewidth', 0.5);
                %     text(ga_mon(mon,mon*10)-0.1, grid.dim3.si(mon*10), num2str(mon), 'color', cmp(idxbin, :), 'fontsize', 6);
                % end
                % xlabel('$\Gamma_m - \Gamma$ (K km$^{-1}$)'); ylabel('$\sigma$ (unitless)');
                % axis('tight');
                % caxis([min(par.r1_bins_45) max(par.r1_bins_45)]);
                % c = colorbar('ticks', [min(par.r1_bins_45):0.2:max(par.r1_bins_45)], 'ticklabels', strtrim(cellstr(num2str(flip([min(par.r1_bins_45):0.2:max(par.r1_bins_45)])', '%.1f'))'), 'ticklabelinterpreter', 'latex', 'ydir', 'reverse');
                % if ~strcmp(fw, 'mse_old')
                %     ylabel(c, '$R_1$ (unitless)', 'interpreter', 'latex');
                % else
                %     ylabel(c, '$R_1^*$ (unitless)', 'interpreter', 'latex');
                % end
                % make_title_type_lat_pt(type, lat_eval, par);
                % set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                % set(gca, 'fontsize', par.fs, 'xlim', [-2 6], 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.3 1], 'xminortick', 'on')
                % print(sprintf('%s/ga_diff_r1_mon_1_6.png', folder), '-dpng', '-r300');
                % close;

                % %% Months 7--12
                % figure(); clf; hold all; box on;
                % line([0 0], [0 1], 'color', 'k', 'linewidth', 0.5);
                % r1_bins_fine = linspace(par.r1_bins_45(1), par.r1_bins_45(end), 101);
                % cmp = flip(parula(length(r1_bins_fine)-1));
                % %cmp = parula(length(par.r1_bins_45)-1);
                % for mon = 7:12
                %     [~, idxbin] = min(abs(r1(mon) - r1_bins_fine));
                %     plot(ga_mon(mon,:), grid.dim3.si, 'color', cmp(idxbin,:), 'linewidth', 0.5);
                %     text(ga_mon(mon,(mon-6)*10)-0.1, grid.dim3.si((mon-6)*10), num2str(mon), 'color', cmp(idxbin, :), 'fontsize', 6);
                % end
                % xlabel('$\Gamma_m - \Gamma$ (K km$^{-1}$)'); ylabel('$\sigma$ (unitless)');
                % axis('tight');
                % caxis([min(par.r1_bins_45) max(par.r1_bins_45)]);
                % c = colorbar('ticks', [min(par.r1_bins_45):0.2:max(par.r1_bins_45)], 'ticklabels', strtrim(cellstr(num2str(flip([min(par.r1_bins_45):0.2:max(par.r1_bins_45)])', '%.1f'))'), 'ticklabelinterpreter', 'latex', 'ydir', 'reverse');
                % if ~strcmp(fw, 'mse_old')
                %     ylabel(c, '$R_1$ (unitless)', 'interpreter', 'latex');
                % else
                %     ylabel(c, '$R_1^*$ (unitless)', 'interpreter', 'latex');
                % end
                % make_title_type_lat_pt(type, lat_eval, par);
                % set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                % set(gca, 'fontsize', par.fs, 'xlim', [-2 6], 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.3 1], 'xminortick', 'on')
                % print(sprintf('%s/ga_diff_r1_mon_7_12.png', folder), '-dpng', '-r300');
                % close;

                %%%%%%%% GA FRAC %%%%%%%
                %% Months 1--6
                figure(); clf; hold all; box on;
                line([0 0], [0 1], 'color', 'k', 'linewidth', 0.5);
                line(100*[1 1], [0 1], 'linestyle', '--', 'color', 'k', 'linewidth', 0.5);
                r1_bins_fine = linspace(par.r1_bins_45(1), par.r1_bins_45(end), 101);
                cmp = flip(parula(length(r1_bins_fine)-1));
                %cmp = parula(length(par.r1_bins_45)-1);
                for mon = 1:6
                    [~, idxbin] = min(abs(r1(mon) - r1_bins_fine));
                    plot(ga_fr_mon(mon,:), grid.dim3.si, 'color', cmp(idxbin,:), 'linewidth', 0.5);
                    if abs(lat_eval) < 60
                        text(ga_fr_mon(mon,mon*10)-0.1*abs(ga_fr_mon(mon,mon*10)), grid.dim3.si(mon*10), num2str(mon), 'color', cmp(idxbin, :), 'fontsize', 6);
                    else
                        text(ga_fr_mon(mon,mon*3)-0.1*abs(ga_fr_mon(mon,mon*3)), grid.dim3.si(mon*3), num2str(mon), 'color', cmp(idxbin, :), 'fontsize', 6);
                    end

                    vavg_dev_rce = nanmean(interp1(grid.dim3.si, ga_fr_mon(mon,:), linspace(0.9, 0.3, 101)));
                    disp(sprintf('The vertically integrated lapse rate deviation from MALR during Month %g at %g deg latitude is %g%%.', mon, lat_eval, vavg_dev_rce))
                end
                xlabel('$\frac{\Gamma_m - \Gamma}{\Gamma_m}$ (\%)'); ylabel('$\sigma$ (unitless)');
                axis('tight');
                % caxis([min(par.r1_bins_45) max(par.r1_bins_45)]);
                % c = colorbar('ticks', [min(par.r1_bins_45):0.2:max(par.r1_bins_45)], 'ticklabels', strtrim(cellstr(num2str(flip([min(par.r1_bins_45):0.2:max(par.r1_bins_45)])', '%.1f'))'), 'ticklabelinterpreter', 'latex', 'ydir', 'reverse');
                % if ~strcmp(fw, 'mse_old')
                %     ylabel(c, '$R_1$ (unitless)', 'interpreter', 'latex');
                % else
                %     ylabel(c, '$R_1^*$ (unitless)', 'interpreter', 'latex');
                % end
                make_title_type_lat_pt(type, lat_eval, par);
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.3 1], 'xminortick', 'on', 'xlim', [-50 150])
                print(sprintf('%s/ga_fr_r1_mon_1_6.png', folder), '-dpng', '-r300');
                close;

                %% Months 7--12
                figure(); clf; hold all; box on;
                line([0 0], [0 1], 'color', 'k', 'linewidth', 0.5);
                line(100*[1 1], [0 1], 'linestyle', '--', 'color', 'k', 'linewidth', 0.5);
                r1_bins_fine = linspace(par.r1_bins_45(1), par.r1_bins_45(end), 101);
                cmp = flip(parula(length(r1_bins_fine)-1));
                %cmp = parula(length(par.r1_bins_45)-1);
                for mon = 7:12
                    [~, idxbin] = min(abs(r1(mon) - r1_bins_fine));
                    plot(ga_fr_mon(mon,:), grid.dim3.si, 'color', cmp(idxbin,:), 'linewidth', 0.5);
                    if abs(lat_eval) < 60
                        text(ga_fr_mon(mon,(mon-6)*10)-0.2*abs(ga_fr_mon(mon,(mon-6)*10)), grid.dim3.si((mon-6)*10), num2str(mon), 'color', cmp(idxbin, :), 'fontsize', 6);
                    else
                        text(ga_fr_mon(mon,(mon-6)*3)-0.2*abs(ga_fr_mon(mon,(mon-6)*3)), grid.dim3.si((mon-6)*3), num2str(mon), 'color', cmp(idxbin, :), 'fontsize', 6);
                    end

                    vavg_dev_rce = nanmean(interp1(grid.dim3.si, ga_fr_mon(mon,:), linspace(0.9, 0.3, 101)));
                    disp(sprintf('The vertically integrated lapse rate deviation from MALR during Month %g at %g deg latitude is %g%%.', mon, lat_eval, vavg_dev_rce))
                end
                xlabel('$\frac{\Gamma_m - \Gamma}{\Gamma_m}$ (\%)'); ylabel('$\sigma$ (unitless)');
                axis('tight');
                % caxis([min(par.r1_bins_45) max(par.r1_bins_45)]);
                % c = colorbar('ticks', [min(par.r1_bins_45):0.2:max(par.r1_bins_45)], 'ticklabels', strtrim(cellstr(num2str(flip([min(par.r1_bins_45):0.2:max(par.r1_bins_45)])', '%.1f'))'), 'ticklabelinterpreter', 'latex', 'ydir', 'reverse');
                % if ~strcmp(fw, 'mse_old')
                %     ylabel(c, '$R_1$ (unitless)', 'interpreter', 'latex');
                % else
                %     ylabel(c, '$R_1^*$ (unitless)', 'interpreter', 'latex');
                % end
                make_title_type_lat_pt(type, lat_eval, par);
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.3 1], 'xminortick', 'on', 'xlim', [-50 150])
                print(sprintf('%s/ga_fr_r1_mon_7_12.png', folder), '-dpng', '-r300');
                close;

            end % lat bounds

        end % land/ocean
    end % framework

end
