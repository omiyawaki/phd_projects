function plot_temp_binned_r1_hl(type, par)
    make_dirs(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/flux_z.mat', prefix_proc)); % load lat x mon RCAE_ALT data
    load(sprintf('%s/ta_mon_lat.mat', prefix_proc));
    load(sprintf('%s/ma_mon_lat.mat', prefix_proc));

    lat_bound_list = [80];
    
    for f = {'mse', 'mse_old'}; fw = f{1};
        % for l = {'lo', 'l', 'o'}; land = l{1};
        for l = {'lo'}; land = l{1};
            if strcmp(land, 'lo'); land_text = 'Land + Ocean';
            elseif strcmp(land, 'l'); land_text = 'Land';
            elseif strcmp(land, 'o'); land_text = 'Ocean';
            end

            for lb = 1:length(lat_bound_list); lat_bound = lat_bound_list(lb);
                dlat = 0.25; % step size for standard lat grid
                if lat_bound>0; lat_pole = 90; lat = lat_bound:dlat:lat_pole; monlabel=par.monlabel; shiftby=0;
                else lat_pole = -90; lat = lat_bound:-dlat:lat_pole; monlabel=par.monlabelsh; shiftby=6; end;
                clat = cosd(lat); % cosine of latitude for cosine weighting

                folder = sprintf('%s/temp_binned_r1/%s/%s/0_poleward_of_lat_%g', plotdir, fw, land, lat_bound);
                if ~exist(folder, 'dir'); mkdir(folder); end;
                
                % interpolate to high latitude bounds
                res = interp1(grid.dim2.lat, flux_z.(land).res.(fw), lat);
                ra = interp1(grid.dim2.lat, flux_z.(land).ra.(fw), lat);
                ta = interp1(grid.dim3.lat, tasi.(land), lat);
                ma = interp1(grid.dim3.lat, masi.(land), lat);
            
                % reshape array to vector
                lat = repmat(lat', [1, 12]);
                lat = lat(:);
                clat = cosd(lat);
                clat_vert = repmat(clat, [1, length(grid.dim3.si)]);
                r1 = res(:)./ra(:);
                ta = reshape(ta, [], length(grid.dim3.si));
                ma = reshape(ma, [], length(grid.dim3.si));
                
                id_bins = discretize(r1, par.r1_bins_hl);

                [ta_area, ma_area] = deal(nan([length(par.r1_bins_hl)-1, length(grid.dim3.si)]));

                for bin = 1:length(par.r1_bins_hl)-1
                    idx_bins = find(id_bins==bin);

                    ta_area(bin,:) = nansum(ta(idx_bins,:).*clat_vert(idx_bins),1)/nansum(clat(idx_bins));
                    ma_area(bin,:) = nansum(ma(idx_bins,:).*clat_vert(idx_bins),1)/nansum(clat(idx_bins));

                    % moist adiabats have nans below the initialization level so the nansums are 0 there. Make these spurious zeros nans.
                    if ~strcmp(par.ma_init, 'surf')
                        ma_area(bin,grid.dim3.si>par.ma_init) = nan;
                    end
                end

                [~,idx09]=min(abs(par.r1_bins_hl-0.9));
                figure(); clf; hold all; box on;
                cmp = flip(parula(length(par.r1_bins_hl)-1));
                %cmp = parula(length(par.r1_bins_hl)-1);
                for bin = 1:length(par.r1_bins_hl)-1
                    if bin==idx09
                        plot(ta_area(bin,:), grid.dim3.si, 'color', cmp(bin,:), 'linewidth', 1.5);
                    else
                        plot(ta_area(bin,:), grid.dim3.si, 'color', cmp(bin,:), 'linewidth', 0.5);
                    end
                end
                xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
                axis('tight');
                caxis([min(par.r1_bins_hl) max(par.r1_bins_hl)]);
                c = colorbar('ticks', [min(par.r1_bins_hl):0.2:max(par.r1_bins_hl)], 'ticklabels', strtrim(cellstr(num2str(flip([min(par.r1_bins_hl):0.2:max(par.r1_bins_hl)])', '%.1f'))'), 'ticklabelinterpreter', 'latex', 'ydir', 'reverse');
                if ~strcmp(fw, 'mse_old')
                    ylabel(c, '$R_1$ (unitless)', 'interpreter', 'latex');
                else
                    ylabel(c, '$R_1^*$ (unitless)', 'interpreter', 'latex');
                end
                make_title_type(type, par);
                set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                set(gca, 'fontsize', par.fs, 'xlim', [210 275], 'xtick', [200:10:300], 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.2 1], 'xminortick', 'on')
                print(sprintf('%s/temp_r1_all.png', folder), '-dpng', '-r300');
                close;
                
                % for bin = 1:length(par.r1_bins)-1
                    % figure(); clf; hold all; box on;
                    % plot(ta_area, grid.dim3.si, 'k');
                    % if par.r1_bins(bin) < 0.7
                    %     plot(ma_area, grid.dim3.si, ':k');
                    % end
                    % xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
                    % title(sprintf('$%g \\le R_1 < %g$', round(par.r1_bins(bin), 1), round(par.r1_bins(bin+1), 1)));
                    % % legend([h_rce, h_rcae, h_rae], 'RCE', 'RCAE', 'RAE', 'location', 'southwest');
                    % axis('tight');
                    % set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
                    % set(gca, 'fontsize', par.fs, 'xlim', [210 300], 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.2 1], 'xminortick', 'on')
                    % print(sprintf('%s/temp_binned_r1/%s/%s/temp_r1_%g_to_%g.png', plotdir, fw, land, round(par.r1_bins(bin),1), round(par.r1_bins(bin+1),1)), '-dpng', '-r300');
                    % close;
                % end
                
            end % lat bounds

        end % land/ocean
    end % framework

end
