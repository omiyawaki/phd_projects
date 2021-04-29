function plot_sice(type, par)
% plot sicedow depth
    make_dirs(type, par)

    % load data
    % [~, ~, ~, lat, par] = load_flux(type, par);

    prefix = make_prefix(type, par);
    plotdir = make_plotdir(type, par);
    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/sice.mat', prefix)); % read sea ice
    load(sprintf('%s/albedo.mat', prefix)); % read albedo
    load(sprintf('%s/sftlf.mat', prefix)); % read land fraction to plot coastline
    load(sprintf('%s/srfc.mat', prefix)); % read tas
    tas = rename_tas(type, srfc);
    ts = rename_ts(type, srfc);
    load(sprintf('%s/stf.mat', prefix)); % read lh
    lh = rename_lh(type, stf);

    if length(size(sftlf)) > 2
        sftlf = squeeze(nanmean(sftlf,3));
    end

    folder_prefix = sprintf('%s/sice/', plotdir);
    folder_list = {'nh_hl', 'sh_hl'};
    for fs = 1:length(folder_list); folder_suffix = folder_list{fs};
        foldername = sprintf('%s%s', folder_prefix, folder_suffix);
        if ~exist(foldername, 'dir')
            mkdir(foldername);
        end
    end

    % NH
    idx_nh = grid.dim2.lat >= 45;
    idx_sh = grid.dim2.lat <= -45;

    [mesh_ll_lat, mesh_ll_lon] = meshgrid(grid.dim2.lon, grid.dim2.lat);
    [mesh_ll_lat_nh, mesh_ll_lon_nh] = meshgrid(grid.dim2.lon, grid.dim2.lat(idx_nh));
    [mesh_ll_lat_sh, mesh_ll_lon_sh] = meshgrid(grid.dim2.lon, grid.dim2.lat(idx_sh));

    for mon = 1:12
        mon_str = make_mon_str(mon);
        [xmin, xmax] = make_lon_bounds(type);

        %%%%%%%% STEREOGRAPHIC PROJECTION %%%%%%%%%%%

        %%%%%% NORTH %%%%%%%%%

        % SEA ICE
        var_text = 'Sea ice fraction (1)';
        figure(); clf; hold all; box on;
        colormap('parula')
        contourfpsn(mesh_ll_lon_nh, mesh_ll_lat_nh, squeeze(sice(:,idx_nh,mon))', 0:0.1:1, 'linecolor', 'none', 'km');
        contourpsn(mesh_ll_lon_nh, mesh_ll_lat_nh, sftlf(:,idx_nh)', 0.5*[1 1], 'k', 'km');
        contourpsn(mesh_ll_lon_nh, mesh_ll_lat_nh, squeeze(sice(:,idx_nh,mon))', 0.5*[1 1], 'r', 'km');
        contourpsn(mesh_ll_lon_nh, mesh_ll_lat_nh, mesh_ll_lon_nh, 80*[1 1], ':', 'linecolor', 'k', 'km');
        caxis([0 1]);
        cb = colorbar('limits', [0 1], 'ytick', [0:0.2:1], 'location', 'eastoutside');
        cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        ylabel(cb, var_text);
        xlabel('eastings (km)');
        ylabel('northings (km)');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_gold);
        set(gca, 'xtick', [-5000:2500:5000], 'xminortick', 'on', 'ytick', [-5000:2500:5000], 'yminortick', 'on');
        print(sprintf('%s/sice/nh_hl/stereo_sice_ll_%02d', plotdir, mon), '-dpng', '-r300');

        % ALBEDO
        var_text = 'Albedo (1)';
        figure(); clf; hold all; box on;
        colormap('parula')
        contourfpsn(mesh_ll_lon_nh, mesh_ll_lat_nh, squeeze(albedo(:,idx_nh,mon))', 0:0.1:1, 'linecolor', 'none', 'km');
        contourpsn(mesh_ll_lon_nh, mesh_ll_lat_nh, sftlf(:,idx_nh)', 0.6*[1 1], 'k', 'km');
        contourpsn(mesh_ll_lon_nh, mesh_ll_lat_nh, squeeze(albedo(:,idx_nh,mon))', 0.5*[1 1], 'r', 'km');
        contourpsn(mesh_ll_lon_nh, mesh_ll_lat_nh, mesh_ll_lon_nh, 80*[1 1], ':', 'linecolor', 'k', 'km');
        caxis([0 1]);
        cb = colorbar('limits', [0 1], 'ytick', [0:0.2:1], 'location', 'eastoutside');
        cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        ylabel(cb, var_text);
        xlabel('eastings (km)');
        ylabel('northings (km)');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_gold);
        set(gca, 'xtick', [-5000:2500:5000], 'xminortick', 'on', 'ytick', [-5000:2500:5000], 'yminortick', 'on');
        print(sprintf('%s/sice/nh_hl/stereo_albedo_ll_%02d', plotdir, mon), '-dpng', '-r300');

        % TS
        var_text = '$T_{s}$ (K)';
        figure(); clf; hold all; box on;
        colormap('parula')
        contourfpsn(mesh_ll_lon_nh, mesh_ll_lat_nh, squeeze(ts(:,idx_nh,mon))', 200:2.5:300,'linecolor', 'none', 'km');
        contourpsn(mesh_ll_lon_nh, mesh_ll_lat_nh, sftlf(:,idx_nh)', 0.5*[1 1], 'k', 'km');
        contourpsn(mesh_ll_lon_nh, mesh_ll_lat_nh, squeeze(ts(:,idx_nh,mon))', 271.15*[1 1], 'r', 'km');
        contourpsn(mesh_ll_lon_nh, mesh_ll_lat_nh, mesh_ll_lon_nh, 80*[1 1], ':', 'linecolor', 'k', 'km');
        caxis([200 270]);
        cb = colorbar('limits', [200 270], 'ytick', [200:10:270], 'location', 'eastoutside');
        cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        ylabel(cb, var_text);
        xlabel('eastings (km)');
        ylabel('northings (km)');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_gold);
        set(gca, 'xtick', [-5000:2500:5000], 'xminortick', 'on', 'ytick', [-5000:2500:5000], 'yminortick', 'on');
        print(sprintf('%s/sice/nh_hl/stereo_ts_ll_%02d', plotdir, mon), '-dpng', '-r300');

        % LH
        var_text = 'LH (W m$^{-2}$)';
        figure(); clf; hold all; box on;
        cmp = colCog(40);
        colormap(cmp)
        contourfpsn(mesh_ll_lon_nh, mesh_ll_lat_nh, squeeze(lh(:,idx_nh,mon))', 0:2.5:50, 'linecolor', 'none', 'km');
        contourpsn(mesh_ll_lon_nh, mesh_ll_lat_nh, sftlf(:,idx_nh)', 0.5*[1 1], 'k', 'km');
        contourpsn(mesh_ll_lon_nh, mesh_ll_lat_nh, mesh_ll_lon_nh, 80*[1 1], ':', 'linecolor', 'k', 'km');
        % contourpsn(mesh_ll_lon_nh, mesh_ll_lat_nh, squeeze(lh(:,idx_nh,mon))', 273.15*[1 1], 'r', 'km');
        caxis([-50 50]);
        cb = colorbar('limits', [0 50], 'ytick', [0:10:50], 'location', 'eastoutside');
        cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        ylabel(cb, var_text);
        xlabel('eastings (km)');
        ylabel('northings (km)');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_gold);
        set(gca, 'xtick', [-5000:2500:5000], 'xminortick', 'on', 'ytick', [-5000:2500:5000], 'yminortick', 'on');
        print(sprintf('%s/sice/nh_hl/stereo_lh_ll_%02d', plotdir, mon), '-dpng', '-r300');


        %%%%%% SOUTH %%%%%%%%%

        % SEA ICE
        var_text = 'Sea ice fraction (1)';
        figure(); clf; hold all; box on;
        colormap('parula')
        contourfps(mesh_ll_lon_sh, mesh_ll_lat_sh, squeeze(sice(:,idx_sh,mon))', 0:0.1:1, 'linecolor', 'none', 'km');
        contourps(mesh_ll_lon_sh, mesh_ll_lat_sh, sftlf(:,idx_sh)', 0.5*[1 1], 'k', 'km');
        contourps(mesh_ll_lon_sh, mesh_ll_lat_sh, squeeze(sice(:,idx_sh,mon))', 0.5*[1 1], 'r', 'km');
        contourps(mesh_ll_lon_sh, mesh_ll_lat_sh, mesh_ll_lon_sh, -80*[1 1], ':', 'linecolor', 'k', 'km');
        caxis([0 1]);
        cb = colorbar('limits', [0 1], 'ytick', [0:0.2:1], 'location', 'eastoutside');
        cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        ylabel(cb, var_text);
        xlabel('eastings (km)');
        ylabel('northings (km)');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_gold);
        set(gca, 'xtick', [-5000:2500:5000], 'xminortick', 'on', 'ytick', [-5000:2500:5000], 'yminortick', 'on');
        print(sprintf('%s/sice/sh_hl/stereo_sice_ll_%02d', plotdir, mon), '-dpng', '-r300');

        % ALBEDO
        var_text = 'Albedo (1)';
        figure(); clf; hold all; box on;
        colormap('parula')
        contourfps(mesh_ll_lon_sh, mesh_ll_lat_sh, squeeze(albedo(:,idx_sh,mon))', 0:0.1:1, 'linecolor', 'none', 'km');
        contourps(mesh_ll_lon_sh, mesh_ll_lat_sh, sftlf(:,idx_sh)', 0.6*[1 1], 'k', 'km');
        contourps(mesh_ll_lon_sh, mesh_ll_lat_sh, squeeze(albedo(:,idx_sh,mon))', 0.5*[1 1], 'r', 'km');
        contourps(mesh_ll_lon_sh, mesh_ll_lat_sh, mesh_ll_lon_sh, -80*[1 1], ':', 'linecolor', 'k', 'km');
        caxis([0 1]);
        cb = colorbar('limits', [0 1], 'ytick', [0:0.2:1], 'location', 'eastoutside');
        cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        ylabel(cb, var_text);
        xlabel('eastings (km)');
        ylabel('northings (km)');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_gold);
        set(gca, 'xtick', [-5000:2500:5000], 'xminortick', 'on', 'ytick', [-5000:2500:5000], 'yminortick', 'on');
        print(sprintf('%s/sice/sh_hl/stereo_albedo_ll_%02d', plotdir, mon), '-dpng', '-r300');

        % TS
        var_text = '$T_{s}$ (K)';
        figure(); clf; hold all; box on;
        colormap('parula')
        contourfps(mesh_ll_lon_sh, mesh_ll_lat_sh, squeeze(ts(:,idx_sh,mon))', 200:2.5:300,'linecolor', 'none', 'km');
        contourps(mesh_ll_lon_sh, mesh_ll_lat_sh, sftlf(:,idx_sh)', 0.5*[1 1], 'k', 'km');
        contourps(mesh_ll_lon_sh, mesh_ll_lat_sh, squeeze(ts(:,idx_sh,mon))', 271.15*[1 1], 'r', 'km');
        contourps(mesh_ll_lon_sh, mesh_ll_lat_sh, mesh_ll_lon_sh, -80*[1 1], ':', 'linecolor', 'k', 'km');
        caxis([200 270]);
        cb = colorbar('limits', [200 270], 'ytick', [200:10:270], 'location', 'eastoutside');
        cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        ylabel(cb, var_text);
        xlabel('eastings (km)');
        ylabel('northings (km)');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_gold);
        set(gca, 'xtick', [-5000:2500:5000], 'xminortick', 'on', 'ytick', [-5000:2500:5000], 'yminortick', 'on');
        print(sprintf('%s/sice/sh_hl/stereo_ts_ll_%02d', plotdir, mon), '-dpng', '-r300');

        % LH
        var_text = 'LH (W m$^{-2}$)';
        figure(); clf; hold all; box on;
        cmp = colCog(40);
        colormap(cmp)
        contourfps(mesh_ll_lon_sh, mesh_ll_lat_sh, squeeze(lh(:,idx_sh,mon))', 0:2.5:50, 'linecolor', 'none', 'km');
        contourps(mesh_ll_lon_sh, mesh_ll_lat_sh, sftlf(:,idx_sh)', 0.5*[1 1], 'k', 'km');
        contourps(mesh_ll_lon_sh, mesh_ll_lat_sh, mesh_ll_lon_sh, -80*[1 1], ':', 'linecolor', 'k', 'km');
        % contourpsn(mesh_ll_lon_sh, mesh_ll_lat_sh, squeeze(lh(:,idx_sh,mon))', 273.15*[1 1], 'r', 'km');
        caxis([-50 50]);
        cb = colorbar('limits', [0 50], 'ytick', [0:10:50], 'location', 'eastoutside');
        cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        ylabel(cb, var_text);
        xlabel('eastings (km)');
        ylabel('northings (km)');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_gold);
        set(gca, 'xtick', [-5000:2500:5000], 'xminortick', 'on', 'ytick', [-5000:2500:5000], 'yminortick', 'on');
        print(sprintf('%s/sice/sh_hl/stereo_lh_ll_%02d', plotdir, mon), '-dpng', '-r300');


        % %%%%%%%%%% GLOBAL %%%%%%%%%
        % var_text = 'Sea ice fraction (1)';
        % figure(); clf; hold all; box on;
        % % cmp = colCog(40);
        % colormap('gray');
        % contourf(mesh_ll_lat, mesh_ll_lon, squeeze(sice(:,:,mon))', [0:0.1:1], 'linecolor', 'none');
        % caxis([0 1]);
        % cb = colorbar('limits', [0 1], 'ytick', [0:0.2:1], 'location', 'eastoutside');
        % cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        % ylabel(cb, var_text);
        % contour(mesh_ll_lat, mesh_ll_lon, sftlf', 0.5*[1 1], 'w', 'linewidth', 1);
        % contour(mesh_ll_lat, mesh_ll_lon, sftlf', 0.5*[1 1], 'k');
        % make_title_type_mon(type, mon_str, par);
        % xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
        % set(gca, 'xlim', [xmin xmax], 'xtick', [0:60:360], 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
        % print(sprintf('%s/sice/sice_ll_%02d', plotdir, mon), '-dpng', '-r300');
        % close;

        % %%%%%%%%%%%%%%%%%%%%% NH HL %%%%%%%%%%%%%%%

        % % PLOT SEA ICE
        % figure(); clf; hold all; box on;
        % % cmp = colCog(40);
        % % colormap(cmp);
        % colormap('gray');
        % contourf(mesh_ll_lat, mesh_ll_lon, squeeze(sice(:,:,mon))', [0:0.1:1], 'linecolor', 'none');
        % caxis([0 1]);
        % cb = colorbar('limits', [0 1], 'ytick', [0:0.2:1], 'location', 'eastoutside');
        % cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        % ylabel(cb, var_text);
        % % contour(mesh_ll_lat, mesh_ll_lon, sftlf', 0.5*[1 1], 'w', 'linewidth', 1);
        % contour(mesh_ll_lat, mesh_ll_lon, sftlf', 0.5*[1 1], 'k');
        % make_title_type_mon(type, mon_str, par);
        % xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
        % set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_superwide)
        % set(gca, 'xlim', [xmin xmax], 'xtick', [0:60:360], 'ylim', [60 90], 'ytick', [60:10:90], 'yminortick', 'on', 'tickdir', 'out');
        % print(sprintf('%s/sice/nh_hl/sice_ll_%02d', plotdir, mon), '-dpng', '-r300');
        % close;

        % % PLOT TAS
        % var_text = '$T_{2\,\mathrm{m}}$ (K)';
        % figure(); clf; hold all; box on;
        % colormap('parula');
        % contourf(mesh_ll_lat, mesh_ll_lon, squeeze(tas(:,:,mon))', [0:5:300], 'linecolor', 'none');
        % contour(mesh_ll_lat, mesh_ll_lon, squeeze(tas(:,:,mon))', 273.15*[1 1], 'linecolor', 'r', 'linewidth', 1);
        % caxis([200 300]);
        % cb = colorbar('limits', [240 300], 'ytick', [240:10:300], 'location', 'eastoutside');
        % cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        % ylabel(cb, var_text);
        % contour(mesh_ll_lat, mesh_ll_lon, sftlf', 0.5*[1 1], 'k');
        % make_title_type_mon(type, mon_str, par);
        % xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
        % set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_superwide)
        % set(gca, 'xlim', [xmin xmax], 'xtick', [0:60:360], 'ylim', [60 90], 'ytick', [60:10:90], 'yminortick', 'on', 'tickdir', 'out');
        % print(sprintf('%s/sice/nh_hl/tas_ll_%02d', plotdir, mon), '-dpng', '-r300');
        % close;

        % % PLOT LH
        % var_text = 'LH (W m$^{-2}$)';
        % figure(); clf; hold all; box on;
        % cmp = colCog(40);
        % colormap(cmp);
        % contourf(mesh_ll_lat, mesh_ll_lon, squeeze(lh(:,:,mon))', [0:5:50], 'linecolor', 'none');
        % contour(mesh_ll_lat, mesh_ll_lon, squeeze(lh(:,:,mon))', 10*[1 1], 'linecolor', 'r', 'linewidth', 1);
        % caxis([-50 50]);
        % cb = colorbar('limits', [0 50], 'ytick', [0:10:50], 'location', 'eastoutside');
        % cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
        % ylabel(cb, var_text);
        % contour(mesh_ll_lat, mesh_ll_lon, sftlf', 0.5*[1 1], 'k');
        % make_title_type_mon(type, mon_str, par);
        % xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
        % set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_superwide)
        % set(gca, 'xlim', [xmin xmax], 'xtick', [0:60:360], 'ylim', [60 90], 'ytick', [60:10:90], 'yminortick', 'on', 'tickdir', 'out');
        % print(sprintf('%s/sice/nh_hl/lh_ll_%02d', plotdir, mon), '-dpng', '-r300');
        % close;

    end

    sice = squeeze(nanmean(sice,1)); % zonal mean
    [mesh_lat, mesh_mon] = meshgrid([1:12], grid.dim2.lat);

    % mon x lat of siceow depth
    figure(); clf; hold all; box on;
    cmp = colCog(20);
    colormap(flip(cmp));
    [C,h] = contourf(mesh_lat, mesh_mon, sice, [0:0.1:1], 'linecolor', 'none');
    % clabel(C, h, 0:0.5:10, 'fontsize', 6, 'interpreter', 'latex', 'color', 0.75*[1 1 1]);
    caxis([-1 1]);
    ylabel('Latitude (deg)');
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        title(sprintf('%s', upper(type)));
    elseif strcmp(type, 'echam')
        title(sprintf('%s, %s', upper(type), par.echam.(par.echam.clim)));
    elseif strcmp(type, 'gcm')
        title(sprintf('%s', par.model));
    end
    cb = colorbar('limits', [0 1], 'ytick', [0:0.2:1], 'location', 'eastoutside');
    cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    ylabel(cb, sprintf('Sea ice fraction (1)'));
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', par.monlabelnh, 'ylim', [-90 90], 'ytick', [-90:30:90], 'yminortick', 'on', 'tickdir', 'out');
    print(sprintf('%s/sice/sice_mon_lat', plotdir), '-dpng', '-r300');
    close;

end
