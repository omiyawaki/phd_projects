function plot_thetaeq_zon_select(type, par)
% temperature profiles at selected locations and month
    make_dirs(type, par)

    % % load data
    % if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
    %     par.plotdir = sprintf('./figures/%s/%s/%s', type, par.(type).yr_span, par.lat_interp);
    %     prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
    %     load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
    %     load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/thetaeq_mon_lat.mat', type, par.lat_interp));
    %     plev = grid.dim3.plev/100;
    % elseif strcmp(type, 'gcm')
    %     par.plotdir = sprintf('./figures/%s/%s/%s/%s', type, par.model, par.gcm.clim, par.lat_interp);
    %     prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
    %     load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.model, par.gcm.clim));
    %     load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/thetaeq_mon_lat.mat', type, par.model, par.gcm.clim, par.lat_interp));
    %     plev = grid.dim3.plev/100;
    % elseif strcmp(type, 'echam')
    %     par.plotdir = sprintf('./figures/%s/%s/%s', type, par.echam.clim, par.lat_interp);
    %     prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
    %     load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/grid.mat', type, par.echam.clim));
    %     load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/thetaeq_mon_lat.mat', type, par.echam.clim, par.lat_interp));
    %     plev = grid.dim3.plev/100;
    % elseif strcmp(type, 'echam_ml')
    %     par.plotdir = sprintf('./figures/%s/%s/%s', type, par.(type).yr_span, par.lat_interp);
    %     prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
    %     load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
    %     load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/thetaeq_mon_lat.mat', type, par.lat_interp));
    %     plev = 1:47;
    % end

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);
    load(sprintf('%s/grid.mat', prefix));
    load(sprintf('%s/thetaeq_mon_lat.mat', prefix_proc));

    lat_pole = 85;
    lat_mid = 45;

    % for l = {'lo', 'l', 'o'}; land = l{1};
    for l = {'lo'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end
        for m = [1 7]; month = m(1);
            if month==1; mon_str = 'January';
            elseif month==7; mon_str = 'July'; end;

            thetaeqsi_mon.(land) = squeeze(thetaeqsi.(land)(:,month,:));

            thetaeqsi_sp = interp1(lat, thetaeqsi_mon.(land), -lat_pole); % sounding at -lat_pole S
            thetaeqsi_np = interp1(lat, thetaeqsi_mon.(land), lat_pole); % sounding at lat_pole N
            thetaeqsi_smid = interp1(lat, thetaeqsi_mon.(land), -lat_mid); % sounding at -lat_mid S
            thetaeqsi_nmid = interp1(lat, thetaeqsi_mon.(land), lat_mid); % sounding at lat_mid N
            thetaeqsi_eq = interp1(lat, thetaeqsi_mon.(land), 0); % sounding at equator

            % NH MID and HIGH ONLY
            figure(); clf; hold all; box on;
            if m == 1
                h_np = plot(thetaeqsi_np, grid.dim3.si, 'color', par.blue);
                h_nmid = plot(thetaeqsi_nmid, grid.dim3.si, 'color', 0.25*[1 1 1]);
            elseif m == 7
                h_np = plot(thetaeqsi_np, grid.dim3.si, 'color', 0.25*[1 1 1]);
                h_nmid = plot(thetaeqsi_nmid, grid.dim3.si, 'color', par.orange);
            end
            xlabel('$\theta_e$ (K)'); ylabel('$\sigma$ (unitless)');
            make_title_type_mon(type, mon_str, par);
            legend([h_np h_nmid], sprintf('%g N', lat_pole), sprintf('%g N', lat_mid), 'location', 'northeast');
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
            set(gca, 'fontsize', par.fs, 'xlim', [250 350], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[0 1000], 'xminortick', 'on')
            print(sprintf('%s/thetaeq_zon_sel/%s/%g/nh_only', plotdir, land, month), '-dpng', '-r300');
            close;

            % SH MID and HIGH ONLY
            figure(); clf; hold all; box on;
            h_sp = plot(thetaeqsi_sp, grid.dim3.si, 'color', par.blue);
            h_smid = plot(thetaeqsi_smid, grid.dim3.si, 'color', 0.25*[1 1 1]);
            xlabel('$\theta_e$ (K)'); ylabel('$\sigma$ (unitless)');
            make_title_type_mon(type, mon_str, par);
            legend([h_sp h_smid], sprintf('%g S', lat_pole), sprintf('%g S', lat_mid), 'location', 'northeast');
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
            set(gca, 'fontsize', par.fs, 'xlim', [250 350], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[0 1000], 'xminortick', 'on')
            print(sprintf('%s/thetaeq_zon_sel/%s/%g/sh_only', plotdir, land, month), '-dpng', '-r300');
            close;

        end
    end
end
