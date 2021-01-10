function plot_dtdz_binned_r1(type, par)
    make_dirs(type, par)
    % load data
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.(type).yr_span, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.(type).yr_span);
    elseif strcmp(type, 'gcm')
        par.plotdir = sprintf('./figures/%s/%s/%s/%s', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    elseif strcmp(type, 'echam')
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/%s/flux_z.mat', prefix_proc, par.lat_interp)); % load lat x mon RCAE_ALT data
    load(sprintf('%s/%s/dtdz_mon_lat.mat', prefix_proc, par.lat_interp));
    % load(sprintf('%s/%s/ma_mon_lat.mat', prefix_proc, par.lat_interp));
    return

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

            figure(); clf; hold all; box on;
            cmp = flip(parula(length(par.r1_bins)-1));
            for bin = 1:length(par.r1_bins)-1
                plot(ta_area(bin,:), grid.dim3.si, 'k', 'color', cmp(bin,:));
            end
            text(5+ta_area(1,50), grid.dim3.si(50), sprintf('$R_1=%g$',(par.r1_bins(1)+par.r1_bins(2))/2), 'fontsize',6, 'rotation', -55);
            text(-5+ta_area(end,50), grid.dim3.si(50), sprintf('$R_1=%g$',(par.r1_bins(end)+par.r1_bins(end-1))/2), 'fontsize',6, 'rotation', -55);
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            axis('tight');
            c = colorbar('ticks', linspace(0,1,ceil(length(par.r1_bins)/2)), 'ticklabels', strtrim(cellstr(num2str(flip([-0.6:0.2:1.4])', '%.1f'))'), 'ticklabelinterpreter', 'latex', 'ydir', 'reverse');
            ylabel(c, '$R_1$ (unitless)', 'interpreter', 'latex');
            if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
                title(sprintf('%s', upper(type)));
            elseif strcmp(type, 'gcm')
                if contains(par.model, 'mmm')
                    title(sprintf('CMIP5 %s', par.gcm.clim));
                else
                    title(sprintf('%s', par.model));
                end
            elseif strcmp(type, 'echam')
                title(sprintf('%s', upper(type)));
            end
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
            set(gca, 'fontsize', par.fs, 'xlim', [200 300], 'xtick', [200:20:300], 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.2 1], 'xminortick', 'on')
            print(sprintf('%s/temp_binned_r1/%s/%s/temp_r1_all.png', par.plotdir, fw, land), '-dpng', '-r300');
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
                % print(sprintf('%s/temp_binned_r1/%s/%s/temp_r1_%g_to_%g.png', par.plotdir, fw, land, round(par.r1_bins(bin),1), round(par.r1_bins(bin+1),1)), '-dpng', '-r300');
                % close;
            % end

        end % land/ocean
    end % framework

end
