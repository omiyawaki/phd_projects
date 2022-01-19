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
    load(sprintf('%s/si_bl_%g/ga_malr_diff_si_mon_lat_%g.mat', prefix_proc, par.si_bl, par.si_up));

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


    % area average over AA
    if strcmp(type, 'rea')
        sftlf_i = interp1(grid_era.dim2.lat, sftlf, lat_mid);
        sftlf_mid = (nansum(clat'.*sftlf_i,1) / nansum(clat));
    else
        sftlf_i = interp1(grid.dim2.lat, sftlf, lat_mid);
        sftlf_mid = nansum(clat.*sftlf_i) / nansum(clat);
    end
    
    ga_i = interp1(grid.dim2.lat, ga_malr_diff.lo, lat_mid);
    ga_mid = nansum(clat'.*ga_i) / nansum(clat);

    ga_land_i = interp1(grid.dim2.lat, ga_malr_diff.l, lat_mid);
    ga_land_mid = nansum(clat'.*ga_land_i) / nansum(clat);
    ga_ocean_i = interp1(grid.dim2.lat, ga_malr_diff.o, lat_mid);
    ga_ocean_mid = nansum(clat'.*ga_ocean_i) / nansum(clat);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % mon x lat of ga
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(); clf; hold all; box on;
    % rcemax=0.1;
    % raemin=0.9;
    % ylim_lo =-0.8;
    % ylim_up =1.1;
    % vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
    % patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
    % vertices = [1 ylim_lo; 12 ylim_lo ;12 rcemax; 1 rcemax];
    % patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
    h=plot(1:12, circshift(ga_mid,shiftby), '-k');
    h_land=plot(1:12, circshift(ga_land_mid,shiftby), '--', 'color', 'k');
    h_ocean=plot(1:12, circshift(ga_ocean_mid,shiftby), ':', 'color', 'k');
    leg = legend([h, h_land, h_ocean], 'Zonal mean', 'Land', 'Ocean', 'location', 'south', 'numcolumns',3);
    leg.ItemTokenSize = [15,1];
    ylabel('$[(\Gamma - \Gamma_m)/\Gamma_m]_{0.7}^{0.3}$ (\%)');
    make_title_type(type, par);
    % set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    print(sprintf('%s/ga_lo/ga_mon_mid', plotdir), '-dpng', '-r300');

    close;

end
