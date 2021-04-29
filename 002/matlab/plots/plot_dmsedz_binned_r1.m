function plot_dmsedz_binned_r1(type, par)
    make_dirs(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/flux_z.mat', prefix_proc)); % load lat x mon RCAE_ALT data
    load(sprintf('%s/dmsedz_mon_lat.mat', prefix_proc));

    for f = par.(type).fw; fw = f{1};
        for l = par.land_list; land = l{1};
            if strcmp(land, 'lo'); land_text = 'Land + Ocean';
            elseif strcmp(land, 'l'); land_text = 'Land';
            elseif strcmp(land, 'o'); land_text = 'Ocean';
            end

            r1 = flux_z.(land).res.(fw)(:)./flux_z.(land).ra.(fw)(:);
            eptlr = reshape(dmsedz.(land), [], length(grid.dim3.si)); % multiply by 100 to express in percentage
            lat = repmat(grid.dim3.lat', [1, 12]);
            lat = lat(:);
            clat = cosd(lat);
            clat_vert = repmat(clat, [1, length(grid.dim3.si)]);

            id_bins = discretize(r1, par.r1_bins);

            [eptlr_area] = deal(nan([length(par.r1_bins)-1, length(grid.dim3.si)]));

            for bin = 1:length(par.r1_bins)-1
                idx_bins = find(id_bins==bin);

                eptlr_area(bin,:) = nansum(eptlr(idx_bins,:).*clat_vert(idx_bins),1)/nansum(clat(idx_bins));

            end

            [~,idx09]=min(abs(par.r1_bins-0.85));
            [~,idx03]=min(abs(par.r1_bins-0.25));
            figure(); clf; hold all; box on;
            line([0 0], [0 1], 'color', 'k', 'linewidth', 0.5);
            line(100*[1 1], [0 1], 'linestyle', '--', 'color', 'k', 'linewidth', 0.5);
            cmp = flip(parula(length(par.r1_bins)-1));
            for bin = 1:length(par.r1_bins)-1
                vavg_dev = nanmean(interp1(grid.dim3.si, eptlr_area(bin,:), linspace(0.9, 0.3, 101)));
                disp(sprintf('For R1 = %g, mse vertically integrated lapse rate deviation from MALR is %g%%.', 1/2*(par.r1_bins(bin)+par.r1_bins(bin+1)), vavg_dev))
                if bin==idx09 | bin==idx03
                    plot(eptlr_area(bin,:), grid.dim3.si, 'k', 'color', cmp(bin,:), 'linewidth', 1.5);
                    % compute vertically-averaged deviation from sigma = 0.9 to 0.3
                else
                    plot(eptlr_area(bin,:), grid.dim3.si, 'k', 'color', cmp(bin,:), 'linewidth', 0.5);
                end
            end
            % draw arrow and label as inversion 
            arrows(105, 0.85, 30, 0, 'Cartesian', [6e-4, 0.1, 0.05, 1e-5]);
            text(105, 0.8, 'Inversion', 'fontsize', 7);
            % text(5+eptlr_area(1,50), grid.dim3.si(50), sprintf('$%g \\le R_1 < %g$',par.r1_bins(1),par.r1_bins(2)), 'fontsize',6, 'roeptlrtion', -55);
            % text(-5+eptlr_area(end,50), grid.dim3.si(50), sprintf('$%g \\le R_1 < %g$',par.r1_bins(end-1),par.r1_bins(end)), 'fontsize',6, 'roeptlrtion', -55);
            xlabel('$\frac{\partial h}{\partial z}$ (kJ km$^{-1}$)'); ylabel('$\sigma$ (unitless)');
            axis('tight');
            caxis([min(par.r1_bins) max(par.r1_bins)]);
            c = colorbar('ticks', [min(par.r1_bins):0.2:max(par.r1_bins)], 'ticklabels', strtrim(cellstr(num2str(flip([-0.6:0.2:1.4])', '%.1f'))'), 'ticklabelinterpreter', 'latex', 'ydir', 'reverse');
            % c = colorbar('ticks', linspace(0,1,ceil(length(par.r1_bins)/2)+1), 'ticklabels', strtrim(cellstr(num2str(flip([-0.6:0.2:1.4])', '%.1f'))'), 'ticklabelinterpreter', 'latex', 'ydir', 'reverse');
            ylabel(c, '$R_1$ (unitless)', 'interpreter', 'latex');
            make_title_type(type, par);
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
            set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.3 1], 'xminortick', 'on', 'xlim', [-10 10], 'yminortick', 'on')
            print(sprintf('%s/mselr_binned_r1/%s/%s/r1_all.png', plotdir, fw, land), '-dpng', '-r300');
            close;

            [~,idx09]=min(abs(par.r1_bins-0.85));
            [~,idx03]=min(abs(par.r1_bins-0.25));
            figure(); clf; hold all; box on;
            line([0 0], [0 1], 'color', 'k', 'linewidth', 0.5);
            line(100*[1 1], [0 1], 'linestyle', '--', 'color', 'k', 'linewidth', 0.5);
            cmp = flip(parula(length(par.r1_bins)-1));
            for bin = 1:length(par.r1_bins)-1
                if bin==idx09 | bin==idx03
                    plot(eptlr_area(bin,:), grid.dim3.si, 'k', 'color', cmp(bin,:), 'linewidth', 1.5);
                else
                    plot(eptlr_area(bin,:), grid.dim3.si, 'k', 'color', cmp(bin,:), 'linewidth', 0.3);
                end
            end
            % text(5+eptlr_area(1,50), grid.dim3.si(50), sprintf('$%g \\le R_1 < %g$',par.r1_bins(1),par.r1_bins(2)), 'fontsize',6, 'roeptlrtion', -55);
            % text(-5+eptlr_area(end,50), grid.dim3.si(50), sprintf('$%g \\le R_1 < %g$',par.r1_bins(end-1),par.r1_bins(end)), 'fontsize',6, 'roeptlrtion', -55);
            xlabel('$\frac{\partial h}{\partial z}$ (kJ km$^{-1}$)'); ylabel('$\sigma$ (unitless)');
            axis('tight');
            make_title_type(type, par);
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_wide)
            set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.3 1], 'xminortick', 'on', 'xlim', [-10 10], 'yminortick', 'on')
            print(sprintf('%s/mselr_binned_r1/%s/%s/r1_all_nocb.png', plotdir, fw, land), '-dpng', '-r300');
            close;

            figure(); clf; hold all; box on;
            axis off
            cmp = flip(parula(length(par.r1_bins)-1));
            c = colorbar('ticks', linspace(0,1,ceil(length(par.r1_bins)/2)+1), 'ticklabels', strtrim(cellstr(num2str(flip([-0.6:0.2:1.4])', '%.1f'))'), 'ticklabelinterpreter', 'latex', 'ydir', 'reverse', 'location', 'north');
            ylabel(c, '$R_1$ (unitless)', 'interpreter', 'latex');
            set(gcf, 'paperunits', 'inches', 'paperposition', [0 0 18/3 2/3]) 
            print(sprintf('%s/mselr_binned_r1/%s/%s/r1_all_largecb.png', plotdir, fw, land), '-dpng', '-r300');
            close;

        end % land/ocean
    end % framework

end
