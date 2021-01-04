function plot_tempz_ann(type, par)

    % load data
    [~, ~, ~, lat, par] = load_flux(type, par);
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.(type).yr_span);
        % load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
        % load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_mon_lat.mat', type, par.lat_interp));
        % load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ma_mon_lat.mat', type, par.lat_interp));
    elseif strcmp(type, 'gcm')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
        % load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.model, par.gcm.clim));
        % load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/ta_lon_lat.mat', type, par.model, par.gcm.clim, par.lat_interp));
        % load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/ma_lon_lat.mat', type, par.model, par.gcm.clim, par.lat_interp));
    elseif strcmp(type, 'echam')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
        % load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/grid.mat', type, par.echam.clim));
        % load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/ta_mon_lat.mat', type, par.echam.clim, par.lat_interp));
        % load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/ma_mon_lat.mat', type, par.echam.clim, par.lat_interp));
    end
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/%s/flux_z.mat', prefix_proc, par.lat_interp)); % load lat x mon RCAE_ALT data
    load(sprintf('%s/%s/ta_mon_lat.mat', prefix_proc, par.lat_interp));
    load(sprintf('%s/%s/ma_mon_lat.mat', prefix_proc, par.lat_interp));

    make_dirs_ep(type, par)

    time = 'ann';
    crit = 'def';

    for f = {'mse'}; fw = f{1};
        % for l = {'lo', 'l', 'o'}; land = l{1};
        for l = {'lo'}; land = l{1};
            if strcmp(land, 'lo'); land_text = 'Land + Ocean';
            elseif strcmp(land, 'l'); land_text = 'Land';
            elseif strcmp(land, 'o'); land_text = 'Ocean';
            end

            % take annual mean
            r1_ann = nanmean(flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), 2);
            ta_ann = squeeze(nanmean(tasi.(land), 2));
            ma_ann = squeeze(nanmean(masi.(land), 2));

            % locate rce, rcae, and rae for NH and SH
            idx_rce_nh = find(r1_ann<=par.ep & grid.dim3.lat>0);
            idx_rcae_nh = find(r1_ann>par.ep & r1_ann<par.ga & grid.dim3.lat>0);
            idx_rae_nh = find(r1_ann>=par.ga & grid.dim3.lat>0);

            idx_rce_sh = find(r1_ann<=par.ep & grid.dim3.lat<0);
            idx_rcae_sh = find(r1_ann>par.ep & r1_ann<par.ga & grid.dim3.lat<0);
            idx_rae_sh = find(r1_ann>=par.ga & grid.dim3.lat<0);

            % take area averaged temperature profile
            ta_rce_nh = nansum( ta_ann(idx_rce_nh, :) .* repmat(cosd(grid.dim3.lat(idx_rce_nh)), [1 size(ta_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rce_nh)));
            ma_rce_nh = nansum( ma_ann(idx_rce_nh, :) .* repmat(cosd(grid.dim3.lat(idx_rce_nh)), [1 size(ma_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rce_nh)));
            ta_rcae_nh = nansum( ta_ann(idx_rcae_nh, :) .* repmat(cosd(grid.dim3.lat(idx_rcae_nh)), [1 size(ta_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rcae_nh)));
            ma_rcae_nh = nansum( ma_ann(idx_rcae_nh, :) .* repmat(cosd(grid.dim3.lat(idx_rcae_nh)), [1 size(ma_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rcae_nh)));
            ta_rae_nh = nansum( ta_ann(idx_rae_nh, :) .* repmat(cosd(grid.dim3.lat(idx_rae_nh)), [1 size(ta_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rae_nh)));

            ta_rce_sh = nansum( ta_ann(idx_rce_sh, :) .* repmat(cosd(grid.dim3.lat(idx_rce_sh)), [1 size(ta_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rce_sh)));
            ma_rce_sh = nansum( ma_ann(idx_rce_sh, :) .* repmat(cosd(grid.dim3.lat(idx_rce_sh)), [1 size(ma_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rce_sh)));
            ta_rcae_sh = nansum( ta_ann(idx_rcae_sh, :) .* repmat(cosd(grid.dim3.lat(idx_rcae_sh)), [1 size(ta_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rcae_sh)));
            ma_rcae_sh = nansum( ma_ann(idx_rcae_sh, :) .* repmat(cosd(grid.dim3.lat(idx_rcae_sh)), [1 size(ma_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rcae_sh)));
            ta_rae_sh = nansum( ta_ann(idx_rae_sh, :) .* repmat(cosd(grid.dim3.lat(idx_rae_sh)), [1 size(ta_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rae_sh)));

            % moist adiabats have nans below the initialization level so the nansums are 0 there. Make these spurious zeros nans.
            if ~strcmp(par.ma_init, 'surf')
                ma_rce_nh(grid.dim3.si>par.ma_init) = nan;
                ma_rcae_nh(grid.dim3.si>par.ma_init) = nan;
                ma_rce_sh(grid.dim3.si>par.ma_init) = nan;
                ma_rcae_sh(grid.dim3.si>par.ma_init) = nan;
            end

            % ALL NH compared with moist adiabat
            figure(); clf; hold all; box on;
            h_rce = plot(ta_rce_nh, grid.dim3.si, 'color', par.orange);
            h_rce_ma_si = plot(ma_rce_nh, grid.dim3.si, ':', 'color', par.orange);
            h_rcae = plot(ta_rcae_nh, grid.dim3.si, 'color', 0.25*[1 1 1]);
            h_rcae_ma_si = plot(ma_rcae_nh, grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
            h_rae = plot(ta_rae_nh, grid.dim3.si, 'color', par.blue);
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            title(sprintf('NH, %s', upper(time)));
            % legend([h_rce, h_rcae, h_rae], 'RCE', 'RCAE', 'RAE', 'location', 'southwest');
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
            set(gca, 'fontsize', par.fs, 'xlim', [210 300], 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.2 1], 'xminortick', 'on')
            print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/all_nh_ann', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
            close;

            % ALL SH compared with moist adiabat
            figure(); clf; hold all; box on;
            h_rce = plot(ta_rce_sh, grid.dim3.si, 'color', par.orange);
            h_rce_ma_si = plot(ma_rce_sh, grid.dim3.si, ':', 'color', par.orange);
            h_rcae = plot(ta_rcae_sh, grid.dim3.si, 'color', 0.25*[1 1 1]);
            h_rcae_ma_si = plot(ma_rcae_sh, grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
            h_rae = plot(ta_rae_sh, grid.dim3.si, 'color', par.blue);
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            title(sprintf('SH, %s', upper(time)));
            % legend([h_rce, h_rcae, h_rae], 'RCE', 'RCAE', 'RAE', 'location', 'southwest');
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
            set(gca, 'fontsize', par.fs, 'xlim', [210 300], 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.2 1], 'xminortick', 'on')
            print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/all_sh_ann', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
            close;

            % NARROW ALL NH compared with moist adiabat
            figure(); clf; hold all; box on;
            h_rce = plot(ta_rce_nh, grid.dim3.si, 'color', par.orange);
            h_rce_ma_si = plot(ma_rce_nh, grid.dim3.si, ':', 'color', par.orange);
            h_rcae = plot(ta_rcae_nh, grid.dim3.si, 'color', 0.25*[1 1 1]);
            h_rcae_ma_si = plot(ma_rcae_nh, grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
            h_rae = plot(ta_rae_nh, grid.dim3.si, 'color', par.blue);
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            title(sprintf('NH, %s', upper(time)));
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_vert)
            set(gca, 'fontsize', par.fs, 'xlim', [190 300], 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.2 1], 'xminortick', 'on')
            print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/all_nh_ann_vert', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
            close;

            % NARROW ALL SH compared with moist adiabat
            figure(); clf; hold all; box on;
            h_rce = plot(ta_rce_sh, grid.dim3.si, 'color', par.orange);
            h_rce_ma_si = plot(ma_rce_sh, grid.dim3.si, ':', 'color', par.orange);
            h_rcae = plot(ta_rcae_sh, grid.dim3.si, 'color', 0.25*[1 1 1]);
            h_rcae_ma_si = plot(ma_rcae_sh, grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
            h_rae = plot(ta_rae_sh, grid.dim3.si, 'color', par.blue);
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            title(sprintf('SH, %s', upper(time)));
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_vert)
            set(gca, 'fontsize', par.fs, 'xlim', [190 300], 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.2 1], 'xminortick', 'on')
            print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp/all_sh_ann_vert', par.plotdir, par.ep, par.ga, fw, crit, land, time), '-dpng', '-r300');
            close;

        end % land
    end % MSE/DSE framework
    % Legend for RCE and RAE separated into NH and SH
end
