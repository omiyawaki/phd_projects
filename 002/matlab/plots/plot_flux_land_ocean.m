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
    % load(sprintf('%s/sftlf.mat', prefix));
    load(sprintf('%s/flux_z.mat', prefix_proc));
    % load(sprintf('%s/si_bl_0.9/ga_malr_bl_diff_poleward_of_lat_80.mat', prefix_proc));

    % % zonal mean land fraction
    % sftlf = squeeze(nanmean(sftlf, 1));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ra land and ocean
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ra = flux_z.lo.ra.mse_old;
    ra_land = flux_z.l.ra.mse_old;
    ra_ocean = flux_z.o.ra.mse_old;

    ra_i = interp1(grid.dim2.lat, ra, lat_mid);
    ra_mid = nansum(clat'.*ra_i) / nansum(clat);
    ra_land_i = interp1(grid.dim2.lat, ra_land, lat_mid);
    ra_land_mid = nansum(clat'.*ra_land_i) / nansum(clat);
    ra_ocean_i = interp1(grid.dim2.lat, ra_ocean, lat_mid);
    ra_ocean_mid = nansum(clat'.*ra_ocean_i) / nansum(clat);

    % compute residual
    ra_res_mid = ra_mid - ra_land_mid - ra_ocean_mid;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % mon x lat of ra
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ylim_lo = -150;
    ylim_up = 50;
    figure(); clf; hold all; box on;
    % h_res=plot(1:12, circshift(ra_res_mid,shiftby), '-.k');
    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
    h=plot(1:12, circshift(ra_mid,shiftby), '-k');
    h_land=plot(1:12, circshift(ra_land_mid,shiftby), '-', 'color', par.green);
    h_ocean=plot(1:12, circshift(ra_ocean_mid,shiftby), '-', 'color', par.blue);
    % legend([h, h_res, h_land, h_ocean], 'Total', 'Residual', 'Land', 'Ocean', 'location', 'north', 'numcolumns',2);
    leg=legend([h, h_land, h_ocean], 'Total', 'Land', 'Ocean', 'location', 'north', 'numcolumns',3);
    leg.ItemTokenSize = [10,1];
    ylabel('$[R_a]$ (W m$^{-2}$)');
    make_title_type(type, par);
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    print(sprintf('%s/r1_lo/ra_mon_mid', plotdir), '-dpng', '-r300');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % div land and ocean
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    div = flux_z.lo.res.mse_old;
    div_land = flux_z.l.res.mse_old;
    div_ocean = flux_z.o.res.mse_old;

    div_i = interp1(grid.dim2.lat, div, lat_mid);
    div_mid = nansum(clat'.*div_i) / nansum(clat);
    div_land_i = interp1(grid.dim2.lat, div_land, lat_mid);
    div_land_mid = nansum(clat'.*div_land_i) / nansum(clat);
    div_ocean_i = interp1(grid.dim2.lat, div_ocean, lat_mid);
    div_ocean_mid = nansum(clat'.*div_ocean_i) / nansum(clat);

    % compute residual
    div_res_mid = div_mid - div_land_mid - div_ocean_mid;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % mon x lat of div
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ylim_lo = -150;
    ylim_up = 50;
    figure(); clf; hold all; box on;
    % h_res=plot(1:12, circshift(div_res_mid,shiftby), '-.k');
    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
    h=plot(1:12, circshift(div_mid,shiftby), '-k');
    h_land=plot(1:12, circshift(div_land_mid,shiftby), '-', 'color', par.green);
    h_ocean=plot(1:12, circshift(div_ocean_mid,shiftby), '-', 'color', par.blue);
    % legend([h, h_res, h_land, h_ocean], 'Total', 'Residual', 'Land', 'Ocean', 'location', 'south', 'numcolumns',2);
    leg=legend([h, h_land, h_ocean], 'Total', 'Land', 'Ocean', 'location', 'south', 'numcolumns',3);
    leg.ItemTokenSize = [10,1];
    ylabel('$[\partial_t m + \nabla\cdot (vm)]$ (W m$^{-2}$)');
    make_title_type(type, par);
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    print(sprintf('%s/r1_lo/div_mon_mid', plotdir), '-dpng', '-r300');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % stf land and ocean
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    stf = flux_z.lo.stf.mse_old;
    stf_land = flux_z.l.stf.mse_old;
    stf_ocean = flux_z.o.stf.mse_old;

    stf_i = interp1(grid.dim2.lat, stf, lat_mid);
    stf_mid = nansum(clat'.*stf_i) / nansum(clat);
    stf_land_i = interp1(grid.dim2.lat, stf_land, lat_mid);
    stf_land_mid = nansum(clat'.*stf_land_i) / nansum(clat);
    stf_ocean_i = interp1(grid.dim2.lat, stf_ocean, lat_mid);
    stf_ocean_mid = nansum(clat'.*stf_ocean_i) / nansum(clat);

    % compute residual
    stf_res_mid = stf_mid - stf_land_mid - stf_ocean_mid;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % mon x lat of stf
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ylim_lo = -50;
    ylim_up = 100;
    figure(); clf; hold all; box on;
    % h_res=plot(1:12, circshift(stf_res_mid,shiftby), '-.k');
    line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
    h=plot(1:12, circshift(stf_mid,shiftby), '-k');
    h_land=plot(1:12, circshift(stf_land_mid,shiftby), '-', 'color', par.green);
    h_ocean=plot(1:12, circshift(stf_ocean_mid,shiftby), '-', 'color', par.blue);
    % legend([h, h_res, h_land, h_ocean], 'Total', 'Residual', 'Land', 'Ocean', 'location', 'south', 'numcolumns',2);
    leg=legend([h, h_land, h_ocean], 'Total', 'Land', 'Ocean', 'location', 'south', 'numcolumns',3);
    leg.ItemTokenSize = [10,1];
    ylabel('$[\mathrm{LH+SH}]$ (W m$^{-2}$)');
    make_title_type(type, par);
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    print(sprintf('%s/r1_lo/stf_mon_mid', plotdir), '-dpng', '-r300');


end
