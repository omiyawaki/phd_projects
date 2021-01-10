function plot_temp_binned_r1(type, par)
    make_dirs(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/flux_z.mat', prefix_proc)); % load lat x mon RCAE_ALT data
    load(sprintf('%s/ta_mon_lat.mat', prefix_proc));
    load(sprintf('%s/ma_mon_lat.mat', prefix_proc));

    for f = {'mse'}; fw = f{1};
        % for l = {'lo', 'l', 'o'}; land = l{1};
        for l = {'lo'}; land = l{1};
            if strcmp(land, 'lo'); land_text = 'Land + Ocean';
            elseif strcmp(land, 'l'); land_text = 'Land';
            elseif strcmp(land, 'o'); land_text = 'Ocean';
            end

            r1 = flux_z.(land).res.(fw)(:)./flux_z.(land).ra.(fw)(:);
            ta = reshape(tasi.(land), [], length(grid.dim3.si));
            ma = reshape(masi.(land), [], length(grid.dim3.si));
            lat = repmat(grid.dim3.lat', [1, 12]);
            lat = lat(:);
            clat = cosd(lat);
            clat_vert = repmat(clat, [1, length(grid.dim3.si)]);

            id_bins = discretize(r1, par.r1_bins);

            [ta_area, ma_area] = deal(nan([length(par.r1_bins)-1, length(grid.dim3.si)]));

            for bin = 1:length(par.r1_bins)-1
                idx_bins = find(id_bins==bin);

                ta_area(bin,:) = nansum(ta(idx_bins,:).*clat_vert(idx_bins),1)/nansum(clat(idx_bins));
                ma_area(bin,:) = nansum(ma(idx_bins,:).*clat_vert(idx_bins),1)/nansum(clat(idx_bins));

                % moist adiabats have nans below the initialization level so the nansums are 0 there. Make these spurious zeros nans.
                if ~strcmp(par.ma_init, 'surf')
                    ma_area(bin,grid.dim3.si>par.ma_init) = nan;
                end
            end

            [~,idx09]=min(abs(par.r1_bins-0.85));
            [~,idx01]=min(abs(par.r1_bins-0.05));
            figure(); clf; hold all; box on;
            cmp = flip(parula(length(par.r1_bins)-1));
            for bin = 1:length(par.r1_bins)-1
                if bin==idx09 | bin==idx01
                    plot(ta_area(bin,:), grid.dim3.si, 'k', 'color', cmp(bin,:), 'linewidth', 1.5);
                else
                    plot(ta_area(bin,:), grid.dim3.si, 'k', 'color', cmp(bin,:), 'linewidth', 0.5);
                end
            end
            % text(5+ta_area(1,50), grid.dim3.si(50), sprintf('$%g \\le R_1 < %g$',par.r1_bins(1),par.r1_bins(2)), 'fontsize',6, 'rotation', -55);
            % text(-5+ta_area(end,50), grid.dim3.si(50), sprintf('$%g \\le R_1 < %g$',par.r1_bins(end-1),par.r1_bins(end)), 'fontsize',6, 'rotation', -55);
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            axis('tight');
            c = colorbar('ticks', linspace(0,1,ceil(length(par.r1_bins)/2)+1), 'ticklabels', strtrim(cellstr(num2str(flip([-0.6:0.2:1.4])', '%.1f'))'), 'ticklabelinterpreter', 'latex', 'ydir', 'reverse');
            ylabel(c, '$R_1$ (unitless)', 'interpreter', 'latex');
            make_title_type(type, par);
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
            set(gca, 'fontsize', par.fs, 'xlim', [200 300], 'xtick', [200:20:300], 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.2 1], 'xminortick', 'on')
            print(sprintf('%s/temp_binned_r1/%s/%s/temp_r1_all.png', plotdir, fw, land), '-dpng', '-r300');
            close;

            figure(); clf; hold all; box on;
            cmp = flip(parula(length(par.r1_bins)-1));
            for bin = idx09-1:idx09+1
                plot(ta_area(bin,:), grid.dim3.si, 'k', 'color', cmp(bin,:));
            end
            % plot(ta_area(idx09,:), grid.dim3.si, 'k', 'color', cmp(idx09,:));
            % text(5+ta_area(1,50), grid.dim3.si(50), sprintf('$%g \\le R_1 < %g$',par.r1_bins(1),par.r1_bins(2)), 'fontsize',6, 'rotation', -55);
            % text(-5+ta_area(end,50), grid.dim3.si(50), sprintf('$%g \\le R_1 < %g$',par.r1_bins(end-1),par.r1_bins(end)), 'fontsize',6, 'rotation', -55);
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            legend(sprintf('$%g \\le R_1 < %g$',par.r1_bins(idx09-1),par.r1_bins(idx09)),...
                   sprintf('$%g \\le R_1 < %g$',par.r1_bins(idx09),par.r1_bins(idx09+1)),...
                   sprintf('$%g \\le R_1 < %g$',par.r1_bins(idx09+1),par.r1_bins(idx09+2)),...
                   'location', 'southoutside');
            axis('tight');
            make_title_type(type, par);
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
            set(gca, 'fontsize', par.fs, 'xtick', [200:20:300], 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.2 1], 'xminortick', 'on')
            print(sprintf('%s/temp_binned_r1/%s/%s/temp_r1_0-9.png', plotdir, fw, land), '-dpng', '-r300');
            close;

            figure(); clf; hold all; box on;
            cmp = flip(parula(length(par.r1_bins)-1));
            for bin = idx01-1:idx01+1
                plot(ta_area(bin,:), grid.dim3.si, 'color', cmp(bin,:));
                plot(ma_area(bin,:), grid.dim3.si, ':', 'color', cmp(bin,:));
            end
            % plot(ta_area(idx09,:), grid.dim3.si, 'k', 'color', cmp(idx09,:));
            % text(5+ta_area(1,50), grid.dim3.si(50), sprintf('$%g \\le R_1 < %g$',par.r1_bins(1),par.r1_bins(2)), 'fontsize',6, 'rotation', -55);
            % text(-5+ta_area(end,50), grid.dim3.si(50), sprintf('$%g \\le R_1 < %g$',par.r1_bins(end-1),par.r1_bins(end)), 'fontsize',6, 'rotation', -55);
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            legend(sprintf('$%g \\le R_1 < %g$',par.r1_bins(idx01-1),par.r1_bins(idx01)),...
                   sprintf('$%g \\le R_1 < %g$',par.r1_bins(idx01),par.r1_bins(idx01+1)),...
                   sprintf('$%g \\le R_1 < %g$',par.r1_bins(idx01+1),par.r1_bins(idx01+2)),...
                   'location', 'southoutside');
            axis('tight');
            make_title_type(type, par);
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
            set(gca, 'fontsize', par.fs, 'xtick', [200:20:300], 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.2 1], 'xminortick', 'on')
            print(sprintf('%s/temp_binned_r1/%s/%s/temp_r1_0-1.png', plotdir, fw, land), '-dpng', '-r300');
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

        end % land/ocean
    end % framework

end
