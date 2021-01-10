function plot_temp_zon(type, par)
    make_dirs(type, par)

    % load data
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.(type).yr_span, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_lon_lat.mat', type, par.lat_interp));
        % load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ma_lon_lat.mat', type, par.lat_interp));
        plev = grid.dim3.plev/100;
    elseif strcmp(type, 'gcm')
        par.plotdir = sprintf('./figures/%s/%s/%s/%s', type, par.model, par.gcm.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s', type, par.model, par.gcm.clim);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/%s/grid.mat', type, par.model, par.gcm.clim));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/%s/ta_lon_lat.mat', type, par.model, par.gcm.clim, par.lat_interp));
        % load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/ma_lon_lat.mat', type, par.model, par.lat_interp));
        plev = grid.dim3.plev/100;
    elseif strcmp(type, 'echam')
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.echam.clim, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.echam.clim);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s/grid.mat', type, par.echam.clim));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/%s/ta_lon_lat.mat', type, par.echam.clim, par.lat_interp));
        plev = grid.dim3.plev/100;
    elseif strcmp(type, 'echam_ml')
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.(type).yr_span, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_lon_lat.mat', type, par.lat_interp));
        plev = 1:47;
    elseif strcmp(type, 'echam_pl')
        par.plotdir = sprintf('./figures/%s/%s/%s', type, par.(type).yr_span, par.lat_interp);
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/%s', type, par.(type).yr_span);
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/grid.mat', type));
        load(sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s/%s/ta_lon_lat.mat', type, par.lat_interp));
        plev = grid.dim3.plev/100;
    end

    for l = {'lo', 'l', 'o'}; land = l{1};
        if strcmp(land, 'lo'); land_text = 'Land + Ocean';
        elseif strcmp(land, 'l'); land_text = 'Land';
        elseif strcmp(land, 'o'); land_text = 'Ocean';
        end
        for t = {'ann', 'djf', 'jja', 'mam', 'son'}; time = t{1};
            ta_z.(land).(time) = squeeze(nanmean(ta_t.(land).(time), 1)); % zonal average
            tasi_z.(land).(time) = squeeze(nanmean(tasi_t.(land).(time), 1)); % zonal average
            if par.do_surf; tai_z.(land).(time) = squeeze(nanmean(tai_t.(land).(time), 1)); end % zonal average
            tsurf_z.(land).(time) = squeeze(nanmean(tsurf_t.(land).(time), 1)); % zonal average
            psurf_z.(land).(time) = squeeze(nanmean(psurf_t.(land).(time), 1)); % zonal average
            zsurf_z.(land).(time) = squeeze(nanmean(zsurf_t.(land).(time), 1)); % zonal average
            ta_sp = ta_z.(land).(time)(7,:); % sounding at -88.5 S
            ta_np = ta_z.(land).(time)(end-6,:); % sounding at 88.5 N
            ta_eq = interp1(lat, ta_z.(land).(time), 0); % sounding at equator

            z_int=par.z;
            taz_z.(land).(time) = squeeze(nanmean(taz_t.(land).(time), 1)); % zonal average
            taz_sp = taz_z.(land).(time)(7,:); % sounding at -88.5 S
            taz_np = taz_z.(land).(time)(end-6,:); % sounding at 88.5 N
            taz_eq = interp1(lat, taz_z.(land).(time), 0); % sounding at equator

            idx_f = min(find(~isnan(ta_z.(land).(time)(:,floor(size(ta_z.(land).(time),2)/2))))); % first non-nan latitude index
            [~, idx_fi] = min(abs(lat - ceil(lat(idx_f)))); % first non-nan integer latitude index
            step = 1; % step size of adjacent latitudes
            nums = 9; % number of latitudes to plot at once

            if par.do_surf; v_vec = {'si'};
            else v_vec = {'si'}; end
            for v = v_vec; vert = v{1};
                for h = {'sh', 'nh', 'eq', 'nmid', 'smid', 'nhall', 'shall'}; hemi = h{1};
                    figure(); clf; hold all; box on;
                    if strcmp(h, 'sh'); indices = idx_fi:step:step*(nums+1);
                    elseif strcmp(h, 'nh'); indices = length(lat)-idx_fi+1:-step:length(lat)-step*nums;
                    elseif strcmp(h, 'eq');
                        if mod(length(lat),2)==0; indices = length(lat)/2+step*floor(nums/2):-step:length(lat)/2-step*floor(nums/2);
                        elseif mod(length(lat),2)==1; indices = floor(length(lat)/2)+1+step*floor(nums/2):-step:floor(length(lat)/2)-step*floor(nums/2); end;
                    elseif strcmp(h, 'nmid');
                        if mod(length(lat),4)==0; indices = 3*length(lat)/4+step*floor(nums/2):-step:3*length(lat)/4-step*floor(nums/2);
                        elseif mod(length(lat),4)==1; indices = floor(3*length(lat)/4)+1+step*floor(nums/2):-step:floor(3*length(lat)/4)-step*floor(nums/2); end;
                    elseif strcmp(h, 'smid');
                        if mod(length(lat),4)==0; indices = length(lat)/4+step*floor(nums/2):-step:length(lat)/4-step*floor(nums/2);
                        elseif mod(length(lat),4)==1; indices = floor(length(lat)/4)+1+step*floor(nums/2):-step:floor(length(lat)/4)-step*floor(nums/2); end;
                    elseif strcmp(h, 'nhall');
                        if mod(length(lat),2)==0; indices = length(lat)-idx_fi+1:-step:length(lat)/2;
                        elseif mod(length(lat),2)==1; indices = length(lat)-idx_fi+1:-step:floor(length(lat)/2); end;
                    elseif strcmp(h, 'shall');
                        if mod(length(lat),2)==0; indices = idx_fi:step:length(lat)/2;
                        elseif mod(length(lat),2)==1; indices = idx_fi:step:floor(length(lat)/2); end;
                    end
                    c = 1;
                    for i = indices;
                        if any(strcmp(hemi, {'nhall', 'shall'}))
                            if strcmp(vert, 'p'); h{c} = plot(ta_z.(land).(time)(i,:), plev);
                            elseif strcmp(vert, 'z'); h{c} = plot(taz_z.(land).(time)(i,:), z_int*10^-3);
                            elseif strcmp(vert, 'si'); h{c} = plot(tasi_z.(land).(time)(i,:), grid.dim3.si);
                            elseif strcmp(vert, 'pi'); h{c} = plot(tai_z.(land).(time)(i,:), par.pa/1e2); end;
                        else
                            if strcmp(vert, 'p'); h{c} = plot(ta_z.(land).(time)(i,:), plev, 'color', (c/(1.3*nums))*[1 1 1]);
                            elseif strcmp(vert, 'z'); h{c} = plot(taz_z.(land).(time)(i,:), z_int*10^-3, 'color', (c/(1.3*nums))*[1 1 1]);
                            elseif strcmp(vert, 'si'); h{c} = plot(tasi_z.(land).(time)(i,:), grid.dim3.si, 'color', (c/(1.3*nums))*[1 1 1]);
                            elseif strcmp(vert, 'pi'); h{c} = plot(tai_z.(land).(time)(i,:), par.pa/1e2, 'color', (c/(1.3*nums))*[1 1 1]); end;
                        end
                        h_label(c) = "$"+string(lat(i))+"^\circ$";
                        c = c+1;
                    end
                    xlabel('T (K)');
                    if any(strcmp(vert, {'p', 'pi'})); ylabel('p (hPa)');
                    elseif any(strcmp(vert, {'si'})); ylabel('$\sigma$ (unitless)');
                    elseif strcmp(vert, 'z'); ylabel('z (km)'); end;
                    legend([h{:}], h_label,  'location', 'eastoutside');
                    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
                        title(sprintf('%s, %s, %s', upper(type), upper(time), land_text));
                    elseif strcmp(type, 'gcm')
                        title(sprintf('%s, %s, %s', par.model, upper(time), land_text));
                    end
                    axis('tight');
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_larger)
                    if any(strcmp(vert, {'p', 'pi'})); set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0:100:1000], 'ylim', [0 1000], 'xminortick', 'on', 'yminortick', 'on')
                    elseif any(strcmp(vert, {'si'})); set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0.1:0.1:1], 'ylim', [0.1 1], 'xminortick', 'on', 'yminortick', 'on')
                    elseif strcmp(vert, 'z'); set(gca, 'fontsize', par.fs, 'ytick', [0:2:20], 'ylim', [0 20], 'xminortick', 'on'); end;
                    print(sprintf('%s/temp_zon/%s/%s/%s/%s', par.plotdir, land, time, vert, hemi), '-dpng', '-r300');
                    close;
                    clear h h_label;
                end % hemi
            end % vert

            figure(); clf; hold all; box on;
            h_eq = plot(ta_eq, plev, 'color', par.orange);
            h_np = plot(ta_np, plev, 'color', par.blue);
            h_sp = plot(ta_sp, plev, '--', 'color', par.blue);
            xlabel('T (K)'); ylabel('p (hPa)');
            legend([h_eq h_np h_sp], 'Equator', '88.5', '$-88.5$', 'location', 'northeast');
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
            set(gca, 'fontsize', par.fs, 'ydir', 'reverse', 'ytick', [0:100:1000], 'ylim', [0 1000], 'xminortick', 'on')
            print(sprintf('%s/temp_zon/%s/%s/p/all', par.plotdir, land, time), '-dpng', '-r300');
            close;

            figure(); clf; hold all; box on;
            h_eq = plot(taz_eq, z_int*10^-3, 'color', par.orange);
            h_np = plot(taz_np, z_int*10^-3, 'color', par.blue);
            h_sp = plot(taz_sp, z_int*10^-3, '--', 'color', par.blue);
            xlabel('T (K)'); ylabel('z (km)');
            legend([h_eq h_np h_sp], 'Equator', '88.5', '-$88.5$', 'location', 'northeast');
            axis('tight');
            set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos_sq)
            set(gca, 'fontsize', par.fs, 'ytick', [0:2:20], 'ylim', [0 20], 'xminortick', 'on')
            print(sprintf('%s/temp_zon/%s/%s/z/all', par.plotdir, land, time), '-dpng', '-r300');
            close;
        end
    end
end
