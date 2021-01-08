function plot_temp_zon_select(type, par)
% temperature profiles at selected locations and month
    make_dirs(type, par)

    % load data
    % if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
    %     plotdir = sprintf('./figures/%s/%s/%s', type, par.(type).yr_span, par.lat_interp);
    %     prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
    %     load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
    %     load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_mon_lat.mat', type, par.lat_interp));
    %     load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ma_mon_lat.mat', type, par.lat_interp));
    %     plev = grid.dim3.plev/100;
    % elseif strcmp(type, 'gcm')
    %     plotdir = sprintf('./figures/%s/%s/%s/%s', type, par.model, par.gcm.clim, par.lat_interp);
    %     prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
    %     load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.model, par.gcm.clim));
    %     load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/ta_mon_lat.mat', type, par.model, par.gcm.clim, par.lat_interp));
    %     load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/ma_mon_lat.mat', type, par.model, par.gcm.clim, par.lat_interp));
    %     plev = grid.dim3.plev/100;
    % elseif strcmp(type, 'echam')
    %     plotdir = sprintf('./figures/%s/%s/%s', type, par.echam.clim, par.lat_interp);
    %     prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
    %     load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/grid.mat', type, par.echam.clim));
    %     load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/ta_mon_lat.mat', type, par.echam.clim, par.lat_interp));
    %     load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/ma_mon_lat.mat', type, par.echam.clim, par.lat_interp));
    %     plev = grid.dim3.plev/100;
    % elseif strcmp(type, 'echam_ml')
    %     plotdir = sprintf('./figures/%s/%s/%s', type, par.(type).yr_span, par.lat_interp);
    %     prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
    %     load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
    %     load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_mon_lat.mat', type, par.lat_interp));
    %     load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ma_mon_lat.mat', type, par.lat_interp));
    %     plev = 1:47;
    % end
    %
    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);

    load(sprintf('%s/grid.mat', prefix));
    load(sprintf('%s/ta_mon_lat.mat', prefix_proc));
    load(sprintf('%s/ma_mon_lat.mat', prefix_proc));

    lat_pole = 85;
    lat_mid = 45;

    % for l = {'lo', 'l', 'o'}; land = l{1};
    for l = {'lo'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end
        for m = [1 6]; month = m(1);
            if month==1; mon_str = 'January';
            elseif month==6; mon_str = 'June';
            elseif month==7; mon_str = 'July'; end;

            tasi_mon.(land) = squeeze(tasi.(land)(:,month,:));
            masi_mon.(land) = squeeze(masi.(land)(:,month,:));

            % remove moist adiabat data below initialization level
            if ~strcmp(par.ma_init, 'surf')
                masi_mon.(land)(:, grid.dim3.si>par.ma_init) = nan;
            end

            tasi_sp(:,m) = interp1(lat, tasi_mon.(land), -lat_pole); % sounding at -lat_pole S
            tasi_np(:,m) = interp1(lat, tasi_mon.(land), lat_pole); % sounding at lat_pole N
            tasi_smid(:,m) = interp1(lat, tasi_mon.(land), -lat_mid); % sounding at -lat_mid S
            tasi_nmid(:,m) = interp1(lat, tasi_mon.(land), lat_mid); % sounding at lat_mid N
            tasi_eq(:,m) = interp1(lat, tasi_mon.(land), 0); % sounding at equator

            masi_sp(:,m) = interp1(lat, masi_mon.(land), -lat_pole); % sounding at -lat_pole S
            masi_np(:,m) = interp1(lat, masi_mon.(land), lat_pole); % sounding at lat_pole N
            masi_smid(:,m) = interp1(lat, masi_mon.(land), -lat_mid); % sounding at -lat_mid S
            masi_nmid(:,m) = interp1(lat, masi_mon.(land), lat_mid); % sounding at lat_mid N
            masi_eq(:,m) = interp1(lat, masi_mon.(land), 0); % sounding at equator

            % ALL
            figure(); clf; hold all; box on;
            if m == 1
                h_np = plot(tasi_np(:,m), grid.dim3.si, 'color', par.blue);
                h_nmid = plot(tasi_nmid(:,m), grid.dim3.si, 'color', 0.25*[1 1 1]);
                h_nmid_ma = plot(masi_nmid(:,m), grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
            elseif m==6 | m == 7
                h_np = plot(tasi_np(:,m), grid.dim3.si, 'color', 0.25*[1 1 1]);
                h_nmid = plot(tasi_nmid(:,m), grid.dim3.si, 'color', par.orange);
                h_nmid_ma = plot(masi_nmid(:,m), grid.dim3.si, ':', 'color', par.orange);
            end
            h_sp = plot(tasi_sp(:,m), grid.dim3.si, '--', 'color', par.blue);
            h_smid = plot(tasi_smid(:,m), grid.dim3.si, '--', 'color', 0.25*[1 1 1]);
            h_smid_ma = plot(masi_smid(:,m), grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
                title(sprintf('%s, %s', upper(type), mon_str));
            elseif strcmp(type, 'gcm')
                if contains(par.model, 'mmm')
                    title(sprintf('CMIP5 %s, %s', par.gcm.clim, mon_str));
                else
                    title(sprintf('%s, %s', par.model, mon_str));
                end
            elseif strcmp(type, 'echam')
                title(sprintf('%s, %s', upper(type), mon_str));
            end
            legend([h_np h_sp h_nmid h_smid], sprintf('%g N', lat_pole), sprintf('%g S', lat_pole), sprintf('%g N', lat_mid), sprintf('%g S', lat_mid), 'location', 'northeast');
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
            set(gca, 'fontsize', par.fs, 'xlim', [nanmin([tasi_np(:,m); tasi_nmid(:,m)]) inf], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
            print(sprintf('%s/temp_zon_sel/%s/%g/all', plotdir, land, month), '-dpng', '-r300');
            close;

            % NH HIGH ONLY
            figure(); clf; hold all; box on;
            if m == 1
                h_np = plot(tasi_np(:,m), grid.dim3.si, 'color', par.blue);
            elseif m==6 | m == 7
                h_np = plot(tasi_np(:,m), grid.dim3.si, 'color', 0.25*[1 1 1]);
            end
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
                title(sprintf('%s, %s', upper(type), mon_str));
            elseif strcmp(type, 'gcm')
                if contains(par.model, 'mmm')
                    title(sprintf('CMIP5 %s, %s', par.gcm.clim, mon_str));
                else
                    title(sprintf('%s, %s', par.model, mon_str));
                end
            elseif strcmp(type, 'echam')
                title(sprintf('%s, %s', upper(type), mon_str));
            end
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
            set(gca, 'fontsize', par.fs, 'xlim', [nanmin([tasi_np(:,m)]) inf], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
            print(sprintf('%s/temp_zon_sel/%s/%g/np', plotdir, land, month), '-dpng', '-r300');
            close;

            % NH MID ONLY
            figure(); clf; hold all; box on;
            if m == 1
                h_nmid = plot(tasi_nmid(:,m), grid.dim3.si, 'color', 0.25*[1 1 1]);
                h_nmid_ma = plot(masi_nmid(:,m), grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
            elseif m==6 | m == 7
                h_nmid = plot(tasi_nmid(:,m), grid.dim3.si, 'color', par.orange);
                h_nmid_ma = plot(masi_nmid(:,m), grid.dim3.si, ':', 'color', par.orange);
            end
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
                title(sprintf('%s, %s', upper(type), mon_str));
            elseif strcmp(type, 'gcm')
                if contains(par.model, 'mmm')
                    title(sprintf('CMIP5 %s, %s', par.gcm.clim, mon_str));
                else
                    title(sprintf('%s, %s', par.model, mon_str));
                end
            elseif strcmp(type, 'echam')
                title(sprintf('%s, %s', upper(type), mon_str));
            end
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
            set(gca, 'fontsize', par.fs, 'xlim', [nanmin([tasi_nmid(:,m)]) inf], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
            print(sprintf('%s/temp_zon_sel/%s/%g/nhmid', plotdir, land, month), '-dpng', '-r300');
            close;

            % NH MID and HIGH ONLY
            figure(); clf; hold all; box on;
            if m == 1
                h_np = plot(tasi_np(:,m), grid.dim3.si, 'color', par.blue);
                h_nmid = plot(tasi_nmid(:,m), grid.dim3.si, 'color', 0.25*[1 1 1]);
                h_nmid_ma = plot(masi_nmid(:,m), grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
            elseif m==6 | m == 7
                h_np = plot(tasi_np(:,m), grid.dim3.si, 'color', 0.25*[1 1 1]);
                h_nmid = plot(tasi_nmid(:,m), grid.dim3.si, 'color', par.orange);
                h_nmid_ma = plot(masi_nmid(:,m), grid.dim3.si, ':', 'color', par.orange);
            end
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
                title(sprintf('%s, %s', upper(type), mon_str));
            elseif strcmp(type, 'gcm')
                if contains(par.model, 'mmm')
                    title(sprintf('CMIP5 %s, %s', par.gcm.clim, mon_str));
                else
                    title(sprintf('%s, %s', par.model, mon_str));
                end
            elseif strcmp(type, 'echam')
                title(sprintf('%s, %s', upper(type), mon_str));
            end
            legend([h_np h_nmid], sprintf('%g N', lat_pole), sprintf('%g N', lat_mid), 'location', 'northeast');
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
            set(gca, 'fontsize', par.fs, 'xlim', [nanmin([tasi_np(:,m); tasi_nmid(:,m)]) inf], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
            print(sprintf('%s/temp_zon_sel/%s/%g/nh_only', plotdir, land, month), '-dpng', '-r300');
            close;

            % SH HIGH ONLY
            figure(); clf; hold all; box on;
            h_sp = plot(tasi_sp(:,m), grid.dim3.si, 'color', par.blue);
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
                title(sprintf('%s, %s', upper(type), mon_str));
            elseif strcmp(type, 'gcm')
                if contains(par.model, 'mmm')
                    title(sprintf('CMIP5 %s, %s', par.gcm.clim, mon_str));
                else
                    title(sprintf('%s, %s', par.model, mon_str));
                end
            elseif strcmp(type, 'echam')
                title(sprintf('%s, %s', upper(type), mon_str));
            end
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
            set(gca, 'fontsize', par.fs, 'xlim', [nanmin([tasi_sp(:,m)]) inf], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
            print(sprintf('%s/temp_zon_sel/%s/%g/sp', plotdir, land, month), '-dpng', '-r300');
            close;

            % SH MID ONLY
            figure(); clf; hold all; box on;
            h_smid = plot(tasi_smid(:,m), grid.dim3.si, 'color', 0.25*[1 1 1]);
            h_smid_ma = plot(masi_smid(:,m), grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
                title(sprintf('%s, %s', upper(type), mon_str));
            elseif strcmp(type, 'gcm')
                if contains(par.model, 'mmm')
                    title(sprintf('CMIP5 %s, %s', par.gcm.clim, mon_str));
                else
                    title(sprintf('%s, %s', par.model, mon_str));
                end
            elseif strcmp(type, 'echam')
                title(sprintf('%s, %s', upper(type), mon_str));
            end
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
            set(gca, 'fontsize', par.fs, 'xlim', [nanmin([tasi_smid(:,m)]) inf], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
            print(sprintf('%s/temp_zon_sel/%s/%g/shmid', plotdir, land, month), '-dpng', '-r300');
            close;

            % SH MID and HIGH ONLY
            figure(); clf; hold all; box on;
            h_sp = plot(tasi_sp(:,m), grid.dim3.si, 'color', par.blue);
            h_smid = plot(tasi_smid(:,m), grid.dim3.si, 'color', 0.25*[1 1 1]);
            h_smid_ma = plot(masi_smid(:,m), grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
                title(sprintf('%s, %s', upper(type), mon_str));
            elseif strcmp(type, 'gcm')
                if contains(par.model, 'mmm')
                    title(sprintf('CMIP5 %s, %s', par.gcm.clim, mon_str));
                else
                    title(sprintf('%s, %s', par.model, mon_str));
                end
            elseif strcmp(type, 'echam')
                title(sprintf('%s, %s', upper(type), mon_str));
            end
            legend([h_sp h_smid], sprintf('%g S', lat_pole), sprintf('%g S', lat_mid), 'location', 'northeast');
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
            set(gca, 'fontsize', par.fs, 'xlim', [nanmin([tasi_sp(:,m); tasi_smid(:,m)]) inf], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
            print(sprintf('%s/temp_zon_sel/%s/%g/sh_only', plotdir, land, month), '-dpng', '-r300');
            close;

            % NARROW NH MID and HIGH ONLY
            figure(); clf; hold all; box on;
            if m == 1
                h_np = plot(tasi_np(:,m), grid.dim3.si, 'color', par.blue);
                h_nmid = plot(tasi_nmid(:,m), grid.dim3.si, 'color', 0.25*[1 1 1]);
                h_nmid_ma = plot(masi_nmid(:,m), grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
            elseif m==6 | m == 7
                h_np = plot(tasi_np(:,m), grid.dim3.si, 'color', 0.25*[1 1 1]);
                h_nmid = plot(tasi_nmid(:,m), grid.dim3.si, 'color', par.orange);
                h_nmid_ma = plot(masi_nmid(:,m), grid.dim3.si, ':', 'color', par.orange);
            end
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
                title(sprintf('%s, %s', upper(type), mon_str));
            elseif strcmp(type, 'gcm')
                if contains(par.model, 'mmm')
                    title(sprintf('CMIP5 %s, %s', par.gcm.clim, mon_str));
                else
                    title(sprintf('%s, %s', par.model, mon_str));
                end
            elseif strcmp(type, 'echam')
                title(sprintf('%s, %s', upper(type), mon_str));
            end
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_vert)
            set(gca, 'fontsize', par.fs, 'xlim', [nanmin([tasi_np(:,m); tasi_nmid(:,m)]) inf], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
            print(sprintf('%s/temp_zon_sel/%s/%g/nh_only_vert', plotdir, land, month), '-dpng', '-r300');
            close;

            % NARROW SH MID and HIGH ONLY
            figure(); clf; hold all; box on;
            h_sp = plot(tasi_sp(:,m), grid.dim3.si, 'color', par.blue);
            h_smid = plot(tasi_smid(:,m), grid.dim3.si, 'color', 0.25*[1 1 1]);
            h_smid_ma = plot(masi_smid(:,m), grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
            xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
            if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
                title(sprintf('%s, %s', upper(type), mon_str));
            elseif strcmp(type, 'gcm')
                title(sprintf('%s, %s', par.model, mon_str));
                if contains(par.model, 'mmm')
                    title(sprintf('CMIP5 %s, %s', par.gcm.clim, mon_str));
                else
                    title(sprintf('%s, %s', par.model, mon_str));
                end
            elseif strcmp(type, 'echam')
                title(sprintf('%s, %s', upper(type), mon_str));
            end
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_vert)
            set(gca, 'fontsize', par.fs, 'xlim', [nanmin([tasi_sp(:,m); tasi_smid(:,m)]) inf], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
            print(sprintf('%s/temp_zon_sel/%s/%g/sh_only_vert', plotdir, land, month), '-dpng', '-r300');
            close;

        end

        % ALL NH
        figure(); clf; hold all; box on;
        h_nmid_jan = plot(tasi_nmid(:,1), grid.dim3.si, 'color', 0.25*[1 1 1]);
        h_nmid_ma_jan = plot(masi_nmid(:,1), grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
        h_nmid_jun = plot(tasi_nmid(:,6), grid.dim3.si, 'color', par.orange);
        h_nmid_ma_jun = plot(masi_nmid(:,6), grid.dim3.si, ':', 'color', par.orange);
        h_np_jan = plot(tasi_np(:,1), grid.dim3.si, 'color', par.blue);
        h_np_jun = plot(tasi_np(:,6), grid.dim3.si, 'color', 0.25*[1 1 1]);
        xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
        if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
            title(sprintf('%s, NH', upper(type)));
        elseif strcmp(type, 'gcm')
            if contains(par.model, 'mmm')
                title(sprintf('CMIP5 %s, NH', par.gcm.clim));
            else
                title(sprintf('%s, NH', par.model));
            end
        elseif strcmp(type, 'echam')
            title(sprintf('%s, NH', upper(type)));
        end
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
        set(gca, 'fontsize', par.fs, 'xlim', [200 290], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
        print(sprintf('%s/temp_zon_sel/%s/nh_all', plotdir, land), '-dpng', '-r300');
        close;

        % ALL SH
        figure(); clf; hold all; box on;
        h_smid_jan = plot(tasi_smid(:,1), grid.dim3.si, 'color', 0.25*[1 1 1]);
        h_smid_ma_jan = plot(masi_smid(:,1), grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
        h_smid_jun = plot(tasi_smid(:,6), grid.dim3.si, 'color', 0.25*[1 1 1]);
        h_smid_ma_jun = plot(masi_smid(:,6), grid.dim3.si, ':', 'color', 0.25*[1 1 1]);
        h_sp_jan = plot(tasi_sp(:,1), grid.dim3.si, 'color', par.blue);
        h_sp_jun = plot(tasi_sp(:,6), grid.dim3.si, 'color', par.blue);
        xlabel('T (K)'); ylabel('$\sigma$ (unitless)');
        if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
            title(sprintf('%s, SH', upper(type)));
        elseif strcmp(type, 'gcm')
            if contains(par.model, 'mmm')
                title(sprintf('CMIP5 %s, SH', par.gcm.clim));
            else
                title(sprintf('%s, SH', par.model));
            end
        elseif strcmp(type, 'echam')
            title(sprintf('%s, SH', upper(type)));
        end
        axis('tight');
        set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
        set(gca, 'fontsize', par.fs, 'xlim', [200 290], 'ydir', 'reverse', 'ytick', 1e-3*[0:100:1000], 'ylim', 1e-3*[200 1000], 'xminortick', 'on')
        print(sprintf('%s/temp_zon_sel/%s/sh_all', plotdir, land), '-dpng', '-r300');
        close;

    end
end
