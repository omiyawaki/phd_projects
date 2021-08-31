function plot_ga_malr_binned_r1_lat(type, par)
    titlepar = par;
    titlepar.gcm.clim = 'abrupt4$\times$CO2$-$historical';

    
    make_dirs(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/flux_z.mat', prefix_proc)); % load lat x mon RCAE_ALT data
    load(sprintf('%s/dtempsi_mon_lat.mat', prefix_proc)); % load lat x mon RCAE_ALT data
    
    lat_list = [-85 85 45 -45];
    
    for f = {'mse_old'}; fw = f{1};
        for l = par.land_list; land = l{1};
            if strcmp(land, 'lo'); land_text = 'Land + Ocean';
            elseif strcmp(land, 'l'); land_text = 'Land';
            elseif strcmp(land, 'o'); land_text = 'Ocean';
            end

            for lb = 1:length(lat_list); lat_eval = lat_list(lb);

                folder = sprintf('%s/dtempsi_binned_r1/%s/%s/0_lat_%g', plotdir, fw, land, lat_eval);
                if ~exist(folder, 'dir'); mkdir(folder); end;
                
                % interpolate at select latitude
                res = interp1(grid.dim2.lat, flux_z.(land).res.(fw), lat_eval);
                ra = interp1(grid.dim2.lat, flux_z.(land).ra.(fw), lat_eval);
                dta = interp1(grid.dim3.lat, dtempsi.(land), lat_eval);
                dta_mon = squeeze(dta);
            
                % reshape array to vector
                r1 = res(:)./ra(:);
                dta = reshape(dta, [], length(grid.dim3.si));

                id_bins = discretize(r1, par.r1_bins_45);

                [dta_area] = deal(nan([length(par.r1_bins_45)-1, length(grid.dim3.si)]));

                for bin = 1:length(par.r1_bins_45)-1
                    idx_bins = find(id_bins==bin);

                    dta_area(bin,:) = nanmean(dta(idx_bins,:),1);

                end

                %%%%%%%% DTEMP %%%%%%%
                %% Months 1--6
                figure(); clf; hold all; box on;
                r1_bins_fine = linspace(par.r1_bins_45(1), par.r1_bins_45(end), 101);
                cmp = flip(parula(length(r1_bins_fine)-1));
                for mon = 1:6
                    [~, idxbin] = min(abs(r1(mon) - r1_bins_fine));
                    plot(dta_mon(mon,:), grid.dim3.si, 'color', cmp(idxbin,:), 'linewidth', 0.5);
                    if abs(lat_eval) < 60
                        text(dta_mon(mon,mon*10)-0.01*abs(dta_mon(mon,mon*10)), grid.dim3.si(mon*10), num2str(mon), 'color', cmp(idxbin, :), 'fontsize', 6);
                    else
                        text(dta_mon(mon,mon*3)-0.01*abs(dta_mon(mon,mon*3)), grid.dim3.si(mon*3), num2str(mon), 'color', cmp(idxbin, :), 'fontsize', 6);
                    end

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
                % make_title_type_lat_pt(type, lat_eval, titlepar);
                title(sprintf('$\\phi = %g^\\circ$', lat_eval))
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                if abs(lat_eval) < 60
                    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.3 1], 'xminortick', 'on', 'xlim', [2 6])
                else
                    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.3 1], 'xminortick', 'on', 'xlim', [0 18])
                end
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
                    if abs(lat_eval) < 60
                        text(dta_mon(mon,(mon-6)*10)-0.01*abs(dta_mon(mon,(mon-6)*10)), grid.dim3.si((mon-6)*10), num2str(mon), 'color', cmp(idxbin, :), 'fontsize', 6);
                    else
                        text(dta_mon(mon,(mon-6)*3)-0.01*abs(dta_mon(mon,(mon-6)*3)), grid.dim3.si((mon-6)*3), num2str(mon), 'color', cmp(idxbin, :), 'fontsize', 6);
                    end

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
                % make_title_type_lat_pt(type, lat_eval, titlepar);
                title(sprintf('$\\phi = %g^\\circ$', lat_eval))
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                if abs(lat_eval) < 60
                    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.3 1], 'xminortick', 'on', 'xlim', [2 6])
                else
                    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.3 1], 'xminortick', 'on', 'xlim', [0 18])
                end
                print(sprintf('%s/dtempsi_r1_mon_7_12.png', folder), '-dpng', '-r300');
                close;

                %%%%%%%% NORMED DTEMP %%%%%%%
                %% Months 1--6
                figure(); clf; hold all; box on;
                r1_bins_fine = linspace(par.r1_bins_45(1), par.r1_bins_45(end), 101);
                cmp = flip(parula(length(r1_bins_fine)-1));
                for mon = 1:6
                    [~, idxbin] = min(abs(r1(mon) - r1_bins_fine));
                    plot(dta_mon(mon,:)/dta_mon(mon,1), grid.dim3.si, 'color', cmp(idxbin,:), 'linewidth', 0.5);
                    if abs(lat_eval) < 60
                        text(dta_mon(mon,mon*10)-0.01*abs(dta_mon(mon,mon*10)), grid.dim3.si(mon*10), num2str(mon), 'color', cmp(idxbin, :), 'fontsize', 6);
                    else
                        text(dta_mon(mon,mon*3)-0.01*abs(dta_mon(mon,mon*3)), grid.dim3.si(mon*3), num2str(mon), 'color', cmp(idxbin, :), 'fontsize', 6);
                    end

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
                make_title_type_lat_pt(type, lat_eval, titlepar);
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                if abs(lat_eval) < 60
                    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.3 1], 'xminortick', 'on', 'xlim', [3 6])
                else
                    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.3 1], 'xminortick', 'on', 'xlim', [0 18])
                end
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
                    if abs(lat_eval) < 60
                        text(dta_mon(mon,(mon-6)*10)-0.01*abs(dta_mon(mon,(mon-6)*10)), grid.dim3.si((mon-6)*10), num2str(mon), 'color', cmp(idxbin, :), 'fontsize', 6);
                    else
                        text(dta_mon(mon,(mon-6)*3)-0.01*abs(dta_mon(mon,(mon-6)*3)), grid.dim3.si((mon-6)*3), num2str(mon), 'color', cmp(idxbin, :), 'fontsize', 6);
                    end

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
                make_title_type_lat_pt(type, lat_eval, titlepar);
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                if abs(lat_eval) < 60
                    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.3 1], 'xminortick', 'on', 'xlim', [3 6])
                else
                    set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.3 1], 'xminortick', 'on', 'xlim', [0 18])
                end
                print(sprintf('%s/dtempsi_norm_r1_mon_7_12.png', folder), '-dpng', '-r300');
                close;

            end % lat bounds

        end % land/ocean
    end % framework

end
