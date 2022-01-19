function plot_r1_land_ocean(type, par)
    
    sh = 0; % use SH?
    if sh; shiftby=6; else; shiftby=0; end;
    if sh; monlabel=par.monlabelsh; else; monlabel=par.monlabelnh; end

    make_dirs(type, par)
    plotdir = make_plotdir(type, par);
    
    [~, ~, ~, lat, ~] = load_flux(type, par);

    lat_mid = linspace(40,60,50);
    if sh; lat_mid = -lat_mid; end;

    clat = cosd(lat_mid);

    % load echam data
    par.echam_clims = {"rp000046",... % 50 m
                       "rp000135",... % 40 m
                        "rp000145",... % 25 m
                       "rp000141",... % 15 m
                        "rp000086",... % 5 m
                        "rp000172"}; % 3 m

    r1_mid_echam={};
    n_clims = length(par.echam_clims);
    pare = par;
    pare.lat_interp = 'native';
    for i = 1:n_clims; clim=par.echam_clims{i};
        pare.echam.clim = clim;
        % load data
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', 'echam', clim);
        prefix_proc = make_prefix_proc('echam', pare);
        tmp=load(sprintf('%s/grid.mat', prefix)); grid=tmp.grid; clear tmp; % read grid data
        tmp=load(sprintf('%s/flux_z.mat', prefix_proc)); flux=tmp.flux_z; clear tmp; % read clear sky albedoedo data

        r1 = flux.lo.res.mse_old./flux.lo.ra.mse_old;
    
        % area average over AA
        r1_i = interp1(grid.dim2.lat, r1, lat_mid);
        r1_mid_echam.(clim) = nansum(clat'.*r1_i) / nansum(clat);
        
        clear grid flux r1 r1_i

    end                                    

    % load GCM data
    prefix=make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    load(sprintf('%s/grid.mat', prefix));
    load(sprintf('%s/flux_z.mat', prefix_proc));
    % load(sprintf('%s/si_bl_0.9/ga_malr_bl_diff_poleward_of_lat_80.mat', prefix_proc));

    if strcmp(type, 'rea')
        parera = par;
        parera.lat_interp = 'native';
        prefixera=make_prefix('era5c', parera);
        temp=load(sprintf('%s/grid.mat', prefixera)); grid_era = temp.grid; clear temp;
        temp=load(sprintf('%s/sftlf.mat', prefixera)); sftlf_era = temp.sftlf; clear temp;
        % zonal mean land fraction
        sftlf = squeeze(nanmean(sftlf_era, 1));
    else
        sftlf = squeeze(nanmean(sftlf, 1));
    end

    % r1 land and ocean
    r1 = flux_z.lo.res.mse_old./flux_z.lo.ra.mse_old;
    r1_land = flux_z.l.res.mse_old./flux_z.lo.ra.mse_old;
    r1_ocean = flux_z.o.res.mse_old./flux_z.lo.ra.mse_old;

    % r1s land and ocean
    r1s = flux_z.lo.r1.mse_old;
    r1s_land = flux_z.l.r1.mse_old;
    r1s_ocean = flux_z.o.r1.mse_old;

    % r1ss land and ocean
    r1ss_land = flux_z.l.res.mse_old./flux_z.l.ra.mse_old;
    r1ss_ocean = flux_z.o.res.mse_old./flux_z.o.ra.mse_old;

    % r1ss decomp land and ocean
    ar1ss_land = nanmean(r1ss_land,2);
    dr1ss_land = r1ss_land - ar1ss_land;
    ar1ss_ocean = nanmean(r1ss_ocean,2);
    dr1ss_ocean = r1ss_ocean - ar1ss_ocean;
    ares_land = nanmean(flux_z.l.res.mse_old,2);
    dres_land = flux_z.l.res.mse_old - ares_land;
    ara_land = nanmean(flux_z.l.ra.mse_old,2);
    dra_land = flux_z.l.ra.mse_old - ara_land;
    ares_ocean = nanmean(flux_z.o.res.mse_old,2);
    dres_ocean = flux_z.o.res.mse_old - ares_ocean;
    ara_ocean = nanmean(flux_z.o.ra.mse_old,2);
    dra_ocean = flux_z.o.ra.mse_old - ara_ocean;

    dyn_land = dres_land ./ ara_land;
    rad_land = - ares_land ./ (ara_land).^2 .* dra_land;
    hot_land = dr1ss_land - (dyn_land + rad_land);
    dyn_ocean = dres_ocean ./ ara_ocean;
    rad_ocean = - ares_ocean ./ (ara_ocean).^2 .* dra_ocean;
    hot_ocean = dr1ss_ocean - (dyn_ocean + rad_ocean);

    % area average over AA
    if strcmp(type, 'rea')
        sftlf_i = interp1(grid_era.dim2.lat, sftlf, lat_mid);
        sftlf_mid = (nansum(clat'.*sftlf_i,1) / nansum(clat));
    else
        sftlf_i = interp1(grid.dim2.lat, sftlf, lat_mid);
        sftlf_mid = nansum(clat.*sftlf_i) / nansum(clat);
    end
    
    r1_i = interp1(grid.dim2.lat, r1, lat_mid);
    r1_mid = nansum(clat'.*r1_i) / nansum(clat);

    r1s_i = interp1(grid.dim2.lat, r1s, lat_mid);
    r1s_mid = nansum(clat'.*r1s_i) / nansum(clat);

    r1_land_i = interp1(grid.dim2.lat, r1_land, lat_mid);
    r1_land_mid = nansum(clat'.*r1_land_i) / nansum(clat);
    r1_ocean_i = interp1(grid.dim2.lat, r1_ocean, lat_mid);
    r1_ocean_mid = nansum(clat'.*r1_ocean_i) / nansum(clat);

    r1s_land_i = interp1(grid.dim2.lat, r1s_land, lat_mid);
    r1s_land_mid = nansum(clat'.*r1s_land_i) / nansum(clat);
    r1s_ocean_i = interp1(grid.dim2.lat, r1s_ocean, lat_mid);
    r1s_ocean_mid = nansum(clat'.*r1s_ocean_i) / nansum(clat);

    r1ss_land_i = interp1(grid.dim2.lat, r1ss_land, lat_mid);
    r1ss_land_mid = nansum(clat'.*r1ss_land_i) / nansum(clat);
    r1ss_ocean_i = interp1(grid.dim2.lat, r1ss_ocean, lat_mid);
    r1ss_ocean_mid = nansum(clat'.*r1ss_ocean_i) / nansum(clat);

    ar1ss_land_i = interp1(grid.dim2.lat, ar1ss_land, lat_mid);
    ar1ss_land_mid = nansum(clat.*ar1ss_land_i) / nansum(clat);
    ar1ss_ocean_i = interp1(grid.dim2.lat, ar1ss_ocean, lat_mid);
    ar1ss_ocean_mid = nansum(clat.*ar1ss_ocean_i) / nansum(clat);

    dr1ss_land_i = interp1(grid.dim2.lat, dr1ss_land, lat_mid);
    dr1ss_land_mid = nansum(clat'.*dr1ss_land_i) / nansum(clat);
    dr1ss_ocean_i = interp1(grid.dim2.lat, dr1ss_ocean, lat_mid);
    dr1ss_ocean_mid = nansum(clat'.*dr1ss_ocean_i) / nansum(clat);

    dyn_land_i = interp1(grid.dim2.lat, dyn_land, lat_mid);
    dyn_land_mid = nansum(clat'.*dyn_land_i) / nansum(clat);
    rad_land_i = interp1(grid.dim2.lat, rad_land, lat_mid);
    rad_land_mid = nansum(clat'.*rad_land_i) / nansum(clat);
    hot_land_i = interp1(grid.dim2.lat, hot_land, lat_mid);
    hot_land_mid = nansum(clat'.*hot_land_i) / nansum(clat);

    dyn_ocean_i = interp1(grid.dim2.lat, dyn_ocean, lat_mid);
    dyn_ocean_mid = nansum(clat'.*dyn_ocean_i) / nansum(clat);
    rad_ocean_i = interp1(grid.dim2.lat, rad_ocean, lat_mid);
    rad_ocean_mid = nansum(clat'.*rad_ocean_i) / nansum(clat);
    hot_ocean_i = interp1(grid.dim2.lat, hot_ocean, lat_mid);
    hot_ocean_mid = nansum(clat'.*hot_ocean_i) / nansum(clat);

    % compute residual
    r1_res_mid = r1_mid - r1_land_mid - r1_ocean_mid;
    r1s_res_mid = r1_mid - r1s_land_mid - r1s_ocean_mid;
    r1ss_res_mid = r1_mid - sftlf_mid.*r1ss_land_mid - (1-sftlf_mid).*r1ss_ocean_mid;

    r1s_selfres_mid = r1s_mid - r1s_land_mid - r1s_ocean_mid;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % mon x lat of r1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(); clf; hold all; box on;
    rcemax=0.1;
    raemin=0.9;
    ylim_lo =-0.8;
    ylim_up =1.1;
    vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
    vertices = [1 ylim_lo; 12 ylim_lo ;12 rcemax; 1 rcemax];
    patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
    h_res=plot(1:12, circshift(r1_res_mid,shiftby), '-.k');
    h=plot(1:12, circshift(r1_mid,shiftby), '-k');
    h_land=plot(1:12, circshift(r1_land_mid,shiftby), '-', 'color', par.green);
    h_ocean=plot(1:12, circshift(r1_ocean_mid,shiftby), '-', 'color', par.blue);
    leg = legend([h, h_res, h_land, h_ocean], '$R_1$', '$R_1-R_{1,L}-R_{1,O}$', '$R_{1,L}$', '$R_{1,O}$', 'location', 'north', 'numcolumns',2);
    leg.ItemTokenSize = [15,1];
    ylabel('$R_1$ (unitless)');
    make_title_type(type, par);
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    print(sprintf('%s/r1_lo/r1_mon_mid', plotdir), '-dpng', '-r300');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % mon x lat of r1s
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(); clf; hold all; box on;
    rcemax=0.1;
    raemin=0.9;
    ylim_lo =-0.8;
    ylim_up =1.1;
    vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
    vertices = [1 ylim_lo; 12 ylim_lo ;12 rcemax; 1 rcemax];
    patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
    h_res=plot(1:12, circshift(r1s_res_mid,shiftby), '-.k');
    h=plot(1:12, circshift(r1_mid,shiftby), '-k');
    h_land=plot(1:12, circshift(r1s_land_mid,shiftby), '-', 'color', par.green);
    h_ocean=plot(1:12, circshift(r1s_ocean_mid,shiftby), '-', 'color', par.blue);
    leg = legend([h, h_res, h_land, h_ocean], '$R_1$', '$R_1-R_{1,L}^*-R_{1,O}^*$', '$R_{1,L}^*$', '$R_{1,O}^*$', 'location', 'north', 'numcolumns',2);
    leg.ItemTokenSize = [15,1];
    ylabel('$R_1$ (unitless)');
    make_title_type(type, par);
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    print(sprintf('%s/r1_lo/r1s_mon_mid', plotdir), '-dpng', '-r300');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % mon x lat of r1s (selfres)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(); clf; hold all; box on;
    rcemax=0.1;
    raemin=0.9;
    ylim_lo =-0.8;
    ylim_up =1.1;
    vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
    vertices = [1 ylim_lo; 12 ylim_lo ;12 rcemax; 1 rcemax];
    patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
    h_res=plot(1:12, circshift(r1s_selfres_mid,shiftby), '-.k');
    h=plot(1:12, circshift(r1s_mid,shiftby), '-k');
    h_land=plot(1:12, circshift(r1s_land_mid,shiftby), '-', 'color', par.green);
    h_ocean=plot(1:12, circshift(r1s_ocean_mid,shiftby), '-', 'color', par.blue);
    leg = legend([h, h_res, h_land, h_ocean], '$R_1^*$', '$R_1^*-R_{1,L}^*-R_{1,O}^*$', '$R_{1,L}^*$', '$R_{1,O}^*$', 'location', 'north', 'numcolumns',2);
    leg.ItemTokenSize = [15,1];
    ylabel('$R_1$ (unitless)');
    make_title_type(type, par);
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    print(sprintf('%s/r1_lo/r1s_selfres_mon_mid', plotdir), '-dpng', '-r300');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % mon x lat of r1ss
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(); clf; hold all; box on;
    rcemax=0.1;
    raemin=0.9;
    ylim_lo =-0.8;
    ylim_up =1.1;
    vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
    vertices = [1 ylim_lo; 12 ylim_lo ;12 rcemax; 1 rcemax];
    patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
    h_res=plot(1:12, circshift(r1ss_res_mid,shiftby), '-.k');
    h=plot(1:12, circshift(r1_mid,shiftby), '-k');
    h_land=plot(1:12, circshift(sftlf_mid.*r1ss_land_mid,shiftby), '--', 'color', 'k');
    h_ocean=plot(1:12, circshift((1-sftlf_mid).*r1ss_ocean_mid,shiftby), ':', 'color', 'k');
    % leg = legend([h, h_res, h_land, h_ocean], '$R_1$', 'Residual', '$fR_{1,\mathrm{Land}}$', '$(1-f)R_{1,\mathrm{Ocean}}$', 'location', 'north', 'numcolumns',2);
    % leg.ItemTokenSize = [15,1];
    ylabel('$R_1$ (unitless)');
    % make_title_type(type, par);
    title('Reanalysis mean, NH Midlatitudes')
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    print(sprintf('%s/r1_lo/r1ss_mon_mid', plotdir), '-dpng', '-r300');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % mon x lat DECOMPOSITION of r1ss OVER LAND
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(); clf; hold all; box on;
    rcemax=0.1;
    raemin=0.9;
    ylim_lo =-0.8;
    ylim_up =1.1;
    colororder({'k', 'k'});
    yyaxis left
    vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
    vertices = [1 ylim_lo; 12 ylim_lo ;12 rcemax; 1 rcemax];
    patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
    line([1, 12], nanmean(sftlf_mid)*ar1ss_land_mid*[1, 1], 'color', 'k', 'linewidth', 0.5)
    h_r1=plot(1:12, circshift(sftlf_mid.*r1ss_land_mid,shiftby), '-', 'color', 'k');
    ylabel('$R_1$ (unitless)');
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
    % h_res=plot(1:12, circshift(r1ss_res_mid,shiftby), '-.k');
    yyaxis right
    ylimr_lo = -nanmean(sftlf_mid)*ar1ss_land_mid+ylim_lo;
    ylimr_up = -nanmean(sftlf_mid)*ar1ss_land_mid+ylim_up;
    h_hot=plot(1:12, circshift(sftlf_mid.*hot_land_mid,shiftby), '-.', 'color', 'k');
    h_r1=plot(1:12, circshift(sftlf_mid.*dr1ss_land_mid,shiftby), '-', 'color', 'k');
    h_dyn=plot(1:12, circshift(sftlf_mid.*dyn_land_mid,shiftby), '-', 'color', par.maroon);
    h_rad=plot(1:12, circshift(sftlf_mid.*rad_land_mid,shiftby), '-', 'color', par.gray);
    % leg = legend([h_r1, h_res, h_land, h_ocean], '$R_1$', 'Residual', '$fR_{1,\mathrm{Land}}$', '$(1-f)R_{1,\mathrm{Ocean}}$', 'location', 'north', 'numcolumns',2);
    leg.ItemTokenSize = [15,1];
    ylabel('$\Delta R_1$ (unitless)');
    % make_title_type(type, par);
    title('NH Midlatitudes (Land only)');
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylimr_lo ylimr_up], 'yminortick', 'on', 'tickdir', 'out');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    print(sprintf('%s/r1_lo/r1ss_land_decomp_mon_mid', plotdir), '-dpng', '-r300');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % mon x lat DECOMPOSITION of r1ss OVER ocean
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(); clf; hold all; box on;
    rcemax=0.1;
    raemin=0.9;
    ylim_lo =-0.8;
    ylim_up =1.1;
    colororder({'k', 'k'});
    yyaxis left
    vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
    vertices = [1 ylim_lo; 12 ylim_lo ;12 rcemax; 1 rcemax];
    patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
    % h_res=plot(1:12, circshift(r1ss_res_mid,shiftby), '-.k');
    line([1, 12], nanmean(1-sftlf_mid)*ar1ss_ocean_mid*[1, 1], 'color', 'k', 'linewidth', 0.5)
    h_r1=plot(1:12, circshift((1-sftlf_mid).*r1ss_ocean_mid,shiftby), '-', 'color', 'k');
    ylabel('$R_1$ (unitless)');
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
    yyaxis right
    ylimr_lo = -nanmean(1-sftlf_mid)*ar1ss_ocean_mid+ylim_lo;
    ylimr_up = -nanmean(1-sftlf_mid)*ar1ss_ocean_mid+ylim_up;
    h_hot=plot(1:12, circshift((1-sftlf_mid).*hot_ocean_mid,shiftby), '-.', 'color', 'k');
    h_r1=plot(1:12, circshift((1-sftlf_mid).*dr1ss_ocean_mid,shiftby), '-', 'color', 'k');
    h_dyn=plot(1:12, circshift((1-sftlf_mid).*dyn_ocean_mid,shiftby), '-', 'color', par.maroon);
    h_rad=plot(1:12, circshift((1-sftlf_mid).*rad_ocean_mid,shiftby), '-', 'color', par.gray);
    % leg = legend([h_r1, h_res, h_ocean, h_ocean], '$R_1$', 'Residual', '$fR_{1,\mathrm{ocean}}$', '$(1-f)R_{1,\mathrm{Ocean}}$', 'location', 'north', 'numcolumns',2);
    leg.ItemTokenSize = [15,1];
    ylabel('$\Delta R_1$ (unitless)');
    % make_title_type(type, par);
    title('NH Midlatitudes (Ocean only)');
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylimr_lo ylimr_up], 'yminortick', 'on', 'tickdir', 'out');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    print(sprintf('%s/r1_lo/r1ss_ocean_decomp_mon_mid', plotdir), '-dpng', '-r300');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % r1, r1s, and r1ss comparison (LAND)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(); clf; hold all; box on;
    rcemax=0.1;
    raemin=0.9;
    ylim_lo =-0.8;
    ylim_up =1.1;
    vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
    vertices = [1 ylim_lo; 12 ylim_lo ;12 rcemax; 1 rcemax];
    patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
    hr1=plot(1:12, circshift(r1_land_mid,shiftby), '-', 'color', par.green);
    hr1s=plot(1:12, circshift(r1s_land_mid,shiftby), ':', 'color', par.green);
    hr1ss=plot(1:12, circshift(sftlf_mid.*r1ss_land_mid,shiftby), '--', 'color', par.green);
    haqua=plot(1:12, circshift(sftlf_mid.*r1_mid_echam.rp000172,shiftby), '-k');
    % hres=plot(1:12, circshift(r1s_mid-r1_mid,shiftby), '-.k');
    leg=legend([hr1, hr1s, hr1ss, haqua], '$R_{1,L}$', '$R_{1,L}^*$', '$fR_{1,L}^{**}$', '$f R_{1,\mathrm{3\,m\,AQUA}}$', 'location', 'north', 'numcolumns',3);
    leg.ItemTokenSize = [25,1];
    ylabel('$R_1$ (unitless)');
    make_title_type(type, par);
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    print(sprintf('%s/r1_lo/r1_landcomp_mon_mid', plotdir), '-dpng', '-r300');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % r1, r1s, and r1ss comparison (OCEAN)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(); clf; hold all; box on;
    rcemax=0.1;
    raemin=0.9;
    ylim_lo =-0.8;
    ylim_up =1.1;
    vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
    vertices = [1 ylim_lo; 12 ylim_lo ;12 rcemax; 1 rcemax];
    patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
    hr1=plot(1:12, circshift(r1_ocean_mid,shiftby), '-', 'color', par.blue);
    hr1s=plot(1:12, circshift(r1s_ocean_mid,shiftby), ':', 'color', par.blue);
    hr1ss=plot(1:12, circshift((1-sftlf_mid).*r1ss_ocean_mid,shiftby), '--', 'color', par.blue);
    haqua=plot(1:12, circshift((1-sftlf_mid).*r1_mid_echam.rp000141,shiftby), '-k');
    % hres=plot(1:12, circshift(r1s_mid-r1_mid,shiftby), '-.k');
    leg=legend([hr1, hr1s, hr1ss, haqua], '$R_{1,O}$', '$R_{1,O}^*$', '$(1-f)R_{1,O}^{**}$', '$(1-f) R_{1,\mathrm{15\,m\,AQUA}}$', 'location', 'north', 'numcolumns',3);
    leg.ItemTokenSize = [25,1];
    ylabel('$R_1$ (unitless)');
    make_title_type(type, par);
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    print(sprintf('%s/r1_lo/r1_oceancomp_mon_mid', plotdir), '-dpng', '-r300');

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % leg only
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % axis off;
    % axis([100 101 100 101])
    % title('');
    % legend(flip(leg), '5 m', '15 m', '25 m', '40 m', '50 m', 'location', 'northwest', 'numcolumns', 5, 'orientation', 'horizontal');
    % set(gcf, 'paperunits', 'inches', 'paperposition', [0 0 6.2 0.5])
    % print(sprintf('%s/alb/mld_legend', plotdir), '-dpng', '-r300');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % best fit land r1 with AQUA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(); clf; hold all; box on;
    rcemax=0.1;
    raemin=0.9;
    ylim_lo =-1.1;
    ylim_up =1.1;
    vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
    vertices = [1 ylim_lo; 12 ylim_lo ;12 rcemax; 1 rcemax];
    patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
    h_land=plot(1:12, circshift(r1_land_mid,shiftby), '-', 'color', par.green);
    h_aqua=plot(1:12, circshift(r1_mid_echam.rp000172,shiftby), '--', 'color', par.green);
    legend([h_land h_aqua], 'Land', '3 m AQUA', 'location', 'north');
    ylabel('$R_1$ (unitless)');
    make_title_type(type, par);
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    print(sprintf('%s/r1_lo/r1_mon_mid_land_aqua', plotdir), '-dpng', '-r300');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % best fit land r1 with AQUA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(); clf; hold all; box on;
    rcemax=0.1;
    raemin=0.9;
    ylim_lo =-0.8;
    ylim_up =1.1;
    vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
    vertices = [1 ylim_lo; 12 ylim_lo ;12 rcemax; 1 rcemax];
    patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
    h_land=plot(1:12, circshift(sftlf_mid.*r1ss_land_mid,shiftby), '--', 'color', 'k');
    h_aqua=plot(1:12, circshift(r1_mid_echam.rp000086,shiftby), '-', 'color', par.green);
    legend([h_land h_aqua], '$R_{1,\mathrm{Land}}$', '5 m AQUA', 'location', 'north');
    ylabel('$R_1$ (unitless)');
    make_title_type(type, par);
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    print(sprintf('%s/r1_lo/r1_mon_mid_landcomp_aqua', plotdir), '-dpng', '-r300');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % best fit ocean r1 with AQUA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(); clf; hold all; box on;
    rcemax=0.1;
    raemin=0.9;
    ylim_lo =-0.8;
    ylim_up =1.1;
    vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
    vertices = [1 ylim_lo; 12 ylim_lo ;12 rcemax; 1 rcemax];
    patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
    h_ocean=plot(1:12, circshift((1-sftlf_mid).*r1ss_ocean_mid,shiftby), ':', 'color', 'k');
    h_aqua=plot(1:12, circshift(r1_mid_echam.rp000135,shiftby), '-', 'color', par.blue);
    legend([h_ocean h_aqua], '$R_{1,\mathrm{Ocean}}$', '40 m AQUA', 'location', 'north');
    ylabel('$R_1$ (unitless)');
    make_title_type(type, par);
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    print(sprintf('%s/r1_lo/r1_mon_mid_oceancomp_aqua', plotdir), '-dpng', '-r300');

    close;

end
