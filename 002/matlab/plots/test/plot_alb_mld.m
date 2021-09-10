function plot_alb_mld(type, par)
    
    sh = 0; % use SH?
    if sh; shiftby=6; else; shiftby=0; end;
    if sh; monlabel=par.monlabelsh; else; monlabel=par.monlabelnh; end

% plot albedoow depth
    make_dirs(type, par)
    plotdir = make_plotdir(type, par);
    
    if ~strcmp(type, 'echam')
        error('This analysis only works for ECHAM')
    end

    % configurations to loop over
    % par.echam_clims = {            "rp000126",... % 50 m
    %                                  "rp000148",... % 45 m
    %                                  "rp000134",... % 40 m
    %                                  "rp000146",... % 35 m
    %                                  "rp000130",... % 30 m
    %                                  "rp000144"}; % 25 m

    par.echam_clims = {            "rp000126",... % 50 m
                                 "rp000134",... % 40 m
                                     "rp000144"}; % 25 m

    par.echam_clims_all = {            "rp000126",... % 50 m
                                     "rp000148",... % 45 m
                                 "rp000134",... % 40 m
                                     "rp000146",... % 35 m
                                     "rp000130",... % 30 m
                                     "rp000144"}; % 25 m

    colors = {par.blue, par.maroon, par.orange, par.darkbrown, par.brown, par.yellow};

    [~, ~, ~, lat, ~] = load_flux(type, par);

    lat_aa = linspace(80,90,50);
    if sh; lat_aa = -lat_aa; end;

    clat = cosd(lat_aa);
        
    albedo_aa={}; r1_aa={};
    n_clims = length(par.echam_clims);
    n_clims_all = length(par.echam_clims_all);
    for i = 1:n_clims_all; clim=par.echam_clims_all{i};
        par.echam.clim = clim;
        % load data
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, clim);
        prefix_proc = make_prefix_proc(type, par);
        % /project2/tas1/miyawaki/projects/002/data/proc/echam/rp000144/native/si_bl_0.9/ga_malr_bl_diff_poleward_of_lat_80.mat
        tmp=load(sprintf('%s/grid.mat', prefix)); grid=tmp.grid; clear tmp; % read grid data
        tmp=load(sprintf('%s/albedo.mat', prefix)); albedo=tmp.albedo; clear tmp; % read surface albedo
        tmp=load(sprintf('%s/sice.mat', prefix)); sice=tmp.sice; clear tmp; % read sea ice fraction (by grid size)
        tmp=load(sprintf('%s/siced.mat', prefix)); siced=tmp.siced; clear tmp; % read sea iced fraction (by grid size)
        tmp=load(sprintf('%s/flux_z.mat', prefix_proc)); flux=tmp.flux_z; clear tmp; % read clear sky albedoedo data
        tmp=load(sprintf('%s/si_bl_0.9/ga_malr_bl_diff_poleward_of_lat_80.mat', prefix_proc));
        ga_frac_aa.(clim)=tmp.ga_frac_lat.lo; clear tmp; % read boundary layer lapse rate deviation 

        r1 = flux.lo.res.mse_old./flux.lo.ra.mse_old;
        albedo = squeeze(nanmean(albedo,1)); % zonal mean
        sice = squeeze(nanmean(sice,1)); % zonal mean
        siced = squeeze(nanmean(siced,1)); % zonal mean
    
        % area average over AA
        albedo_i = interp1(grid.dim2.lat, albedo, lat_aa);
        albedo_aa.(clim) = nansum(clat'.*albedo_i) / nansum(clat);
        sice_i = interp1(grid.dim2.lat, sice, lat_aa);
        sice_aa.(clim) = nansum(clat'.*sice_i) / nansum(clat);
        siced_i = interp1(grid.dim2.lat, siced, lat_aa);
        siced_aa.(clim) = nansum(clat'.*siced_i) / nansum(clat);
        r1_i = interp1(grid.dim2.lat, r1, lat_aa);
        r1_aa.(clim) = nansum(clat'.*r1_i) / nansum(clat);
        
        clear grid albedo flux r1 albedo_i r1_i

    end                                    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % albedo seasonality overlay
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(); clf; hold all; box on;
    % cmp = colCog(20);
    for i = 1:n_clims_all; clim=par.echam_clims_all{i};
        % leg(i)=plot(1:12, circshift(albedo_aa.(clim),shiftby), '-', 'color', 1/2*(1/2+i/n_clims)*[1 1 1]);
        leg(i)=plot(1:12, circshift(albedo_aa.(clim),shiftby), '-', 'color', colors{i});
    end
    ylabel('Surface albedo (unitless)');
    title('ECHAM w/ ice')
    % legend(leg, '50 m', '45 m', '40 m', '40 m', '30 m', '25 m', 'location', 'southwest')
    % legend(leg, '50 m', '45 m', '40 m', '35 m', '30 m', '25 m', 'location', 'southwest')
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [0 1], 'yminortick', 'on', 'tickdir', 'out');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    print(sprintf('%s/alb/albedo_mon_icemld', plotdir), '-dpng', '-r300');
    close;
    clear leg

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LR deviation seasonality overlay
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(); clf; hold all; box on;
    % cmp = colCog(20);
    invmin=100;
    ylim_lo =-100;
    ylim_up =500;
    vertices = [1 invmin; 12 invmin; 12 ylim_up; 1 ylim_up];
    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
    for i = 1:n_clims; clim=par.echam_clims{i};
        leg(i)=plot(1:12, circshift(ga_frac_aa.(clim),shiftby), '-');
    end
    ylabel('$\left\langle(\Gamma_m - \Gamma)/\Gamma_m\right\rangle_{1.0}^{0.9}$ (\%)');
    title('ECHAM w/ ice')
    % legend(leg, '50 m', '40 m', '25 m', 'location', 'southwest')
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    print(sprintf('%s/alb/gafrac_mon_icemld', plotdir), '-dpng', '-r300');
    close;
    clear leg

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sice seasonality overlay
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(); clf; hold all; box on;
    % cmp = colCog(20);
    for i = 1:n_clims; clim=par.echam_clims{i};
        leg(i)=plot(1:12, circshift(sice_aa.(clim),shiftby), '-');
    end
    ylabel('Sea ice fraction (unitless)');
    title('ECHAM w/ ice')
    % legend(leg, '50 m', '45 m', '40 m', '40 m', '30 m', '25 m', 'location', 'southwest')
    % legend(leg, '50 m', '45 m', '40 m', '35 m', '30 m', '25 m', 'location', 'southwest')
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [0 1], 'yminortick', 'on', 'tickdir', 'out');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    print(sprintf('%s/alb/sice_mon_icemld', plotdir), '-dpng', '-r300');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % leg only
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    axis off;
    axis([100 101 100 101])
    title('');
    legend(flip(leg), '25 m', '40 m', '50 m', 'location', 'northwest', 'numcolumns', 3, 'orientation', 'horizontal');
    set(gcf, 'paperunits', 'inches', 'paperposition', [0 0 7 0.5])
    print(sprintf('%s/alb/icemld_legend', plotdir), '-dpng', '-r300');
    close;
    clear leg

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sice seasonality overlay (ALL)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(); clf; hold all; box on;
    % cmp = colCog(20);
    for i = 1:n_clims_all; clim=par.echam_clims_all{i};
        leg(i)=plot(1:12, circshift(sice_aa.(clim),shiftby), '-', 'color', colors{i});
    end
    ylabel('Sea ice fraction (unitless)');
    title('ECHAM w/ ice')
    % legend(leg, '50 m', '45 m', '40 m', '40 m', '30 m', '25 m', 'location', 'southwest')
    % legend(leg, '50 m', '45 m', '40 m', '35 m', '30 m', '25 m', 'location', 'southwest')
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [0 1], 'yminortick', 'on', 'tickdir', 'out');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    print(sprintf('%s/alb/sice_mon_icemld_all', plotdir), '-dpng', '-r300');
    close;
    clear leg

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % siced seasonality overlay
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(); clf; hold all; box on;
    % cmp = colCog(20);
    for i = 1:n_clims; clim=par.echam_clims{i};
        leg(i)=plot(1:12, circshift(siced_aa.(clim),shiftby), '-');
    end
    ylabel('Sea ice depth (m)');
    title('ECHAM w/ ice')
    % legend(leg, '50 m', '45 m', '40 m', '40 m', '30 m', '25 m', 'location', 'southwest')
    % legend(leg, '50 m', '45 m', '40 m','35 m', '30 m', '25 m', 'location', 'eastoutside')
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [0 8], 'yminortick', 'on', 'tickdir', 'out');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    print(sprintf('%s/alb/siced_mon_icemld', plotdir), '-dpng', '-r300');
    clear leg

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % siced seasonality overlay (ALL)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(); clf; hold all; box on;
    % cmp = colCog(20);
    for i = 1:n_clims_all; clim=par.echam_clims_all{i};
        % leg(i)=plot(1:12, circshift(siced_aa.(clim),shiftby), '-', 'color', 1/2*(1/2+i/n_clims)*[1 1 1]);
        leg(i)=plot(1:12, circshift(siced_aa.(clim),shiftby), '-', 'color', colors{i});
    end
    ylabel('Sea ice depth (m)');
    title('ECHAM w/ ice')
    % legend(leg, '50 m', '45 m', '40 m', '40 m', '30 m', '25 m', 'location', 'southwest')
    % legend(leg, '50 m', '45 m', '40 m','35 m', '30 m', '25 m', 'location', 'eastoutside')
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [0 8], 'yminortick', 'on', 'tickdir', 'out');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    print(sprintf('%s/alb/siced_mon_icemld_all', plotdir), '-dpng', '-r300');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % leg only
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    axis off;
    axis([100 101 100 101])
    title('');
    legend(flip(leg), '25 m', '30 m', '35 m', '40 m', '45 m', '50 m', 'location', 'northwest', 'numcolumns', 6, 'orientation', 'horizontal');
    set(gcf, 'paperunits', 'inches', 'paperposition', [0 0 7 0.5])
    print(sprintf('%s/alb/icemld_all_legend', plotdir), '-dpng', '-r300');
    close;
    clear leg

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % r1 seasonality overlay (ALL)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(); clf; hold all; box on;
    rcemax=0.1;
    raemin=0.9;
    ylim_lo =0.1;
    ylim_up =1.3;
    vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
    vertices = [1 ylim_lo; 12 ylim_lo ;12 rcemax; 1 rcemax];
    patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
    for i = 1:n_clims_all; clim=par.echam_clims_all{i};
        leg(i)=plot(1:12, circshift(r1_aa.(clim),shiftby), '-', 'color', colors{i});
    end
    ylabel('$R_1$ (unitless)');
    title('ECHAM w/ ice')
    % legend(leg, '50 m', '45 m', '40 m', '40 m', '30 m', '25 m', 'location', 'northeast')
    % legend(leg, '50 m', '40 m', '25 m', 'location', 'southwest')
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    print(sprintf('%s/alb/r1_mon_icemld_all', plotdir), '-dpng', '-r300');
    close;
    clear leg

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % r1 seasonality overlay
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(); clf; hold all; box on;
    rcemax=0.1;
    raemin=0.9;
    ylim_lo =0.1;
    ylim_up =1.3;
    vertices = [1 raemin; 12 raemin; 12 ylim_up; 1 ylim_up];
    patch(vertices(:,1), vertices(:,2), par.blue, 'edgecolor', 'none', 'facealpha', 0.5);
    vertices = [1 ylim_lo; 12 ylim_lo ;12 rcemax; 1 rcemax];
    patch(vertices(:,1), vertices(:,2), par.orange, 'edgecolor', 'none', 'facealpha', 0.5);
    for i = 1:n_clims; clim=par.echam_clims{i};
        % leg(i)=plot(1:12, circshift(r1_aa.(clim),shiftby), '-', 'color', 1/2*(1/2+i/n_clims)*[1 1 1]);
        leg(i)=plot(1:12, circshift(r1_aa.(clim),shiftby), '-');
    end
    ylabel('$R_1$ (unitless)');
    title('ECHAM w/ ice')
    % legend(leg, '50 m', '45 m', '40 m', '40 m', '30 m', '25 m', 'location', 'northeast')
    % legend(leg, '50 m', '40 m', '25 m', 'location', 'southwest')
    set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'ylim', [ylim_lo ylim_up], 'yminortick', 'on', 'tickdir', 'out');
    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
    print(sprintf('%s/alb/r1_mon_icemld', plotdir), '-dpng', '-r300');
    close;
    clear leg

end
