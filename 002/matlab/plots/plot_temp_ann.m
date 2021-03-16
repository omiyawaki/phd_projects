function plot_temp_ann(type, par)

    % load data
    % [~, ~, ~, lat, plotdir] = load_flux(type, par);

    % if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
    %     prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
    %     prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.(type).yr_span);
    %     % load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
    %     % load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_mon_lat.mat', type, par.lat_interp));
    %     % load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ma_mon_lat.mat', type, par.lat_interp));
    % elseif strcmp(type, 'gcm')
    %     prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
    %     prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s', type, par.model, par.gcm.clim);
    %     % load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.model, par.gcm.clim));
    %     % load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/ta_lon_lat.mat', type, par.model, par.gcm.clim, par.lat_interp));
    %     % load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/ma_lon_lat.mat', type, par.model, par.gcm.clim, par.lat_interp));
    % elseif strcmp(type, 'echam')
    %     prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
    %     prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.echam.clim);
    %     % load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/grid.mat', type, par.echam.clim));
    %     % load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/ta_mon_lat.mat', type, par.echam.clim, par.lat_interp));
    %     % load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/ma_mon_lat.mat', type, par.echam.clim, par.lat_interp));
    % end

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/flux_z.mat', prefix_proc)); % load lat x mon RCAE_ALT data
    load(sprintf('%s/ta_mon_lat.mat', prefix_proc));
    load(sprintf('%s/ma_mon_lat_%s.mat', prefix_proc, num2str(par.ma_init)));

    if strcmp(type, 'gcm') & contains(par.model, 'mmm')
        grid.dim3.lat = grid.dim3.lat';
    end

    make_dirs_ep(type, par)

    time = 'ann';
    crit = 'def';

    for f = {'mse_old'}; fw = f{1};
        % for l = {'lo', 'l', 'o'}; land = l{1};
        for l = {'lo'}; land = l{1};
            if strcmp(land, 'lo'); land_text = 'Land + Ocean';
            elseif strcmp(land, 'l'); land_text = 'Land';
            elseif strcmp(land, 'o'); land_text = 'Ocean';
            end

            % take annual mean
            r1_ann = nanmean(flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), 2);
            ta_ann = squeeze(nanmean(tasi.(land), 2));
            ta_std_ann = squeeze(nanmean(tasi_std.(land), 2));
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

            ta_std_rce_nh = nansum( ta_std_ann(idx_rce_nh, :) .* repmat(cosd(grid.dim3.lat(idx_rce_nh)), [1 size(ta_std_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rce_nh)));
            ta_std_rcae_nh = nansum( ta_std_ann(idx_rcae_nh, :) .* repmat(cosd(grid.dim3.lat(idx_rcae_nh)), [1 size(ta_std_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rcae_nh)));
            ta_std_rae_nh = nansum( ta_std_ann(idx_rae_nh, :) .* repmat(cosd(grid.dim3.lat(idx_rae_nh)), [1 size(ta_std_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rae_nh)));

            ta_rce_sh = nansum( ta_ann(idx_rce_sh, :) .* repmat(cosd(grid.dim3.lat(idx_rce_sh)), [1 size(ta_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rce_sh)));
            ma_rce_sh = nansum( ma_ann(idx_rce_sh, :) .* repmat(cosd(grid.dim3.lat(idx_rce_sh)), [1 size(ma_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rce_sh)));
            ta_rcae_sh = nansum( ta_ann(idx_rcae_sh, :) .* repmat(cosd(grid.dim3.lat(idx_rcae_sh)), [1 size(ta_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rcae_sh)));
            ma_rcae_sh = nansum( ma_ann(idx_rcae_sh, :) .* repmat(cosd(grid.dim3.lat(idx_rcae_sh)), [1 size(ma_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rcae_sh)));
            ta_rae_sh = nansum( ta_ann(idx_rae_sh, :) .* repmat(cosd(grid.dim3.lat(idx_rae_sh)), [1 size(ta_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rae_sh)));

            ta_std_rce_sh = nansum( ta_std_ann(idx_rce_sh, :) .* repmat(cosd(grid.dim3.lat(idx_rce_sh)), [1 size(ta_std_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rce_sh)));
            ta_std_rcae_sh = nansum( ta_std_ann(idx_rcae_sh, :) .* repmat(cosd(grid.dim3.lat(idx_rcae_sh)), [1 size(ta_std_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rcae_sh)));
            ta_std_rae_sh = nansum( ta_std_ann(idx_rae_sh, :) .* repmat(cosd(grid.dim3.lat(idx_rae_sh)), [1 size(ta_std_ann,2)]) )/nansum(cosd(grid.dim3.lat(idx_rae_sh)));

            % moist adiabats have nans below the initialization level so the nansums are 0 there. Make these spurious zeros nans.
            if ~strcmp(par.ma_init, 'surf')
                ma_rce_nh(grid.dim3.si>par.ma_init) = nan;
                ma_rcae_nh(grid.dim3.si>par.ma_init) = nan;
                ma_rce_sh(grid.dim3.si>par.ma_init) = nan;
                ma_rcae_sh(grid.dim3.si>par.ma_init) = nan;
            end

            % ALL NH compared with moist adiabat
            figure(); clf; hold all; box on;
            if strcmp(type, 'rea') | (strcmp(type, 'gcm') & strcmp(par.model, 'mmm'))
                ta_rce_nh_u = ta_rce_nh + ta_std_rce_nh;
                ta_rce_nh_l = ta_rce_nh - ta_std_rce_nh;
                ta_rcae_nh_u = ta_rcae_nh + ta_std_rcae_nh;
                ta_rcae_nh_l = ta_rcae_nh - ta_std_rcae_nh;
                ta_rae_nh_u = ta_rae_nh + ta_std_rae_nh;
                ta_rae_nh_l = ta_rae_nh - ta_std_rae_nh;

                ta_rce_nh_2 = [ta_rce_nh_l'; flipud(ta_rce_nh_u')];
                ta_rcae_nh_2 = [ta_rcae_nh_l'; flipud(ta_rcae_nh_u')];
                ta_rae_nh_2 = [ta_rae_nh_l'; flipud(ta_rae_nh_u')];
                si2 = [grid.dim3.si'; flipud(grid.dim3.si')];

                fill(ta_rae_nh_2, si2, par.blue, 'facealpha', 0.2, 'edgealpha', 0.2);
                fill(ta_rcae_nh_2, si2, 0.25*[1 1 1], 'facealpha', 0.2, 'edgealpha', 0.2);
                fill(ta_rce_nh_2, si2, par.orange, 'facealpha', 0.2, 'edgealpha', 0.2);
            end
            h_rce = plot(ta_rce_nh, grid.dim3.si, 'color', par.orange);
            h_rce_ma_si = plot(ma_rce_nh, grid.dim3.si, ':', 'color', par.orange);
            h_rcae = plot(ta_rcae_nh, grid.dim3.si, 'color', 0.25*[1 1 1]);
            h_rcae_ma_si = plot(ma_rcae_nh, grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
            h_rae = plot(ta_rae_nh, grid.dim3.si, 'color', par.blue);
            text(ta_rce_nh(60)+5, grid.dim3.si(60), '{RCE}', 'color', par.orange);
            text(ta_rcae_nh(45), grid.dim3.si(45), '{RCAE}', 'color', 0.25*[1 1 1]);
            text(ta_rae_nh(30)-15, grid.dim3.si(30), '{RAE}', 'color', par.blue);
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            title(sprintf('NH, %s', upper(time)));
            % legend([h_rce, h_rcae, h_rae], 'RCE', 'RCAE', 'RAE', 'location', 'southwest');
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
            set(gca, 'fontsize', par.fs, 'xlim', [210 300], 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.2 1], 'xminortick', 'on')
            print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp_%s/all_nh_ann', plotdir, par.ep, par.ga, fw, crit, land, time, num2str(par.ma_init)), '-dpng', '-r300');
            close;

            % ALL SH compared with moist adiabat
            figure(); clf; hold all; box on;
            if strcmp(type, 'rea') | (strcmp(type, 'gcm') & strcmp(par.model, 'mmm'))
                ta_rce_sh_u = ta_rce_sh + ta_std_rce_sh;
                ta_rce_sh_l = ta_rce_sh - ta_std_rce_sh;
                ta_rcae_sh_u = ta_rcae_sh + ta_std_rcae_sh;
                ta_rcae_sh_l = ta_rcae_sh - ta_std_rcae_sh;
                ta_rae_sh_u = ta_rae_sh + ta_std_rae_sh;
                ta_rae_sh_l = ta_rae_sh - ta_std_rae_sh;

                ta_rce_sh_2 = [ta_rce_sh_l'; flipud(ta_rce_sh_u')];
                ta_rcae_sh_2 = [ta_rcae_sh_l'; flipud(ta_rcae_sh_u')];
                ta_rae_sh_2 = [ta_rae_sh_l'; flipud(ta_rae_sh_u')];
                si2 = [grid.dim3.si'; flipud(grid.dim3.si')];

                fill(ta_rae_sh_2, si2, par.blue, 'facealpha', 0.2, 'edgealpha', 0.2);
                fill(ta_rcae_sh_2, si2, 0.25*[1 1 1], 'facealpha', 0.2, 'edgealpha', 0.2);
                fill(ta_rce_sh_2, si2, par.orange, 'facealpha', 0.2, 'edgealpha', 0.2);
            end
            h_rce = plot(ta_rce_sh, grid.dim3.si, 'color', par.orange);
            h_rce_ma_si = plot(ma_rce_sh, grid.dim3.si, ':', 'color', par.orange);
            h_rcae = plot(ta_rcae_sh, grid.dim3.si, 'color', 0.25*[1 1 1]);
            h_rcae_ma_si = plot(ma_rcae_sh, grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
            h_rae = plot(ta_rae_sh, grid.dim3.si, 'color', par.blue);
            text(ta_rce_sh(60)+5, grid.dim3.si(60), '{RCE}', 'color', par.orange);
            text(ta_rcae_sh(45), grid.dim3.si(45), '{RCAE}', 'color', 0.25*[1 1 1]);
            text(ta_rae_sh(30)-15, grid.dim3.si(30), '{RAE}', 'color', par.blue);
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            title(sprintf('SH, %s', upper(time)));
            % legend([h_rce, h_rcae, h_rae], 'RCE', 'RCAE', 'RAE', 'location', 'southwest');
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
            set(gca, 'fontsize', par.fs, 'xlim', [210 300], 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.2 1], 'xminortick', 'on')
            print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp_%s/all_sh_ann', plotdir, par.ep, par.ga, fw, crit, land, time, num2str(par.ma_init)), '-dpng', '-r300');
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
            set(gca, 'fontsize', par.fs, 'xlim', [210 300], 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.2 1], 'xminortick', 'on')
            print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp_%s/all_nh_ann_vert', plotdir, par.ep, par.ga, fw, crit, land, time, num2str(par.ma_init)), '-dpng', '-r300');
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
            set(gca, 'fontsize', par.fs, 'xlim', [210 300], 'ydir', 'reverse', 'yscale', 'linear', 'ytick', [0:0.1:1], 'ylim', [0.2 1], 'xminortick', 'on')
            print(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp_%s/all_sh_ann_vert', plotdir, par.ep, par.ga, fw, crit, land, time, num2str(par.ma_init)), '-dpng', '-r300');
            close;

        end % land
    end % MSE/DSE framework
    % Legend for RCE and RAE separated into NH and SH
end
