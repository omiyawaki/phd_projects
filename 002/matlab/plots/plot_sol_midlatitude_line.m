function plot_sol_midlatitude_line(type, par)
    make_dirs(type, par)

    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.(type).yr_span, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s', type, par.(type).yr_span);
    elseif strcmp(type, 'merra2')
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
    load(sprintf('%s/rad.mat', prefix)); % read grid data

    if any(strcmp(type, {'era5', 'era5c', 'erai'}))
        sol = rad.tsr; clear rad; % read tas
    elseif strcmp(type, 'gcm')
        sol = rad.rsdt-rad.rsut; clear rad; % read tas
    elseif strcmp(type, 'echam')
        sol = rad.srad0; clear rad; % read tas
    end

    sol = squeeze(nanmean(sol, 1)); % zonal mean

    lat_bound_list = [5 -5];

    for lb = 1:length(lat_bound_list); lat_bound = lat_bound_list(lb);
        dlat = 0.25; % step size for standard lat grid
        if lat_bound>0; lat_center=45; lat = [-lat_bound:dlat:lat_bound]+lat_center; shiftby=0; monlabel=par.monlabel;
        else; lat_center=-45; lat = [-lat_bound:-dlat:lat_bound]+lat_center; shiftby=6; monlabel=par.monlabelsh; end;
        clat = cosd(lat); % cosine of latitude for cosine weighting
        clat_mon = repmat(clat', [1 12]);

        [mesh_lat, mesh_mon] = meshgrid(1:12, lat);

        folder = sprintf('%s/sol/0_midlatitude_pm_lat_%g', par.plotdir, lat_center-lat_bound);
        if ~exist(folder, 'dir'); mkdir(folder); end;

        sol_lat = interp1(grid.dim2.lat, sol, lat);
        sol_lat = nansum(sol_lat.*clat_mon)/nansum(clat);

        dsol_lat = sol_lat - repmat(nanmean(sol_lat,2), [1 12]);

        % ALL lat x mon dependence of RCE and RAE
        figure(); clf; hold all; box on;
        % line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
        plot([1:12],  circshift(sol_lat ,shiftby,2), '-k');
        if any(strcmp(type, {'era5', 'era5c', 'erai'})); title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), lat_center-lat_bound, lat_center+lat_bound));
        elseif any(strcmp(type, 'merra2')); title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), lat_center-lat_bound, lat_center+lat_bound));
        elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), par.echam.(par.echam.clim), lat_center-lat_bound, lat_center+lat_bound));
        elseif strcmp(type, 'gcm');
            if contains(par.model, 'mmm')
                title(sprintf('CMIP5 %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.gcm.clim, lat_center-lat_bound, lat_center+lat_bound));
            else
                title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, lat_center-lat_bound, lat_center+lat_bound));
            end
        end;
        % xlabel('Month');
        ylabel(sprintf('$T_{2\\,\\mathrm{m}}$ (K)'));
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
        set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
        print(sprintf('%s/0_mon_sol', folder), '-dpng', '-r300');
        close;

        % ALL lat x mon dependence of RCE and RAE
        figure(); clf; hold all; box on;
        line([1 12], [0 0], 'linewidth', 0.5, 'color', 'k');
        plot([1:12],  circshift(dsol_lat ,shiftby,2), '-k');
        if any(strcmp(type, {'era5', 'era5c', 'erai'})); title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), lat_center-lat_bound, lat_center+lat_bound));
        elseif any(strcmp(type, 'merra2')); title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), lat_center-lat_bound, lat_center+lat_bound));
        elseif any(strcmp(type, 'echam')); title(sprintf('%s, %s, $\\phi=%g^\\circ$ to $%g^\\circ$', upper(type), par.echam.(par.echam.clim), lat_center-lat_bound, lat_center+lat_bound));
        elseif strcmp(type, 'gcm');
            if contains(par.model, 'mmm')
                title(sprintf('CMIP5 %s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.gcm.clim, lat_center-lat_bound, lat_center+lat_bound));
            else
                title(sprintf('%s, $\\phi=%g^\\circ$ to $%g^\\circ$', par.model, lat_center-lat_bound, lat_center+lat_bound));
            end
        end;
        % xlabel('Month');
        ylabel(sprintf('$\\Delta T_{2\\,\\mathrm{m}}$ (K)'));
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
        set(gca, 'xlim', [1 12], 'xtick', [1:12], 'xticklabels', monlabel, 'yminortick', 'on', 'tickdir', 'out');
        print(sprintf('%s/0_mon_dsol', folder), '-dpng', '-r300');
        close;

    end

end % for function
