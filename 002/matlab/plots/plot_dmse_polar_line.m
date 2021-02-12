function plot_dmse_polar_line(type, par)
    make_dirs(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    % load(sprintf('%s/sftlf.mat', prefix)); % read land fraction data
    % landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
    % par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

    % sftlf = nanmean(sftlf, 1); % zonal average
    % sftlf = repmat(sftlf', [1 12]); % expand land fraction data to time

    % lat_bound_list = [-85 -80 -70 70 80 85];
    lat_bound_list = [-80 80];

    for lb = 1:length(lat_bound_list); par.lat_bound = lat_bound_list(lb);

        [lat, clat, clat_mon, par] = make_polar_lat(par);

        tmp = load(sprintf('%s/dmse_poleward_of_lat_%g', prefix_proc, par.lat_bound)); dmse = tmp.dmse; clear tmp;
        
        % for l = {'lo', 'l', 'o'}; land = l{1};
        for l = {'lo'}; land = l{1};
            if strcmp(land, 'lo'); land_text = 'L+O';
            elseif strcmp(land, 'l'); land_text = 'L';
            elseif strcmp(land, 'o'); land_text = 'O';
            end
            if strcmp(type, 'echam');
                if strcmp(par.echam.clim, '20170908')
                    land_text = 'Snowball';
                else
                    land_text = par.echam.(par.echam.clim);
                end
            end;
            [mesh_lat, mesh_mon] = meshgrid(1:12, grid.dim2.lat);

            f_vec = assign_fw(type, par);
            for f = f_vec; fw = f{1};

                par.folder = sprintf('%s/dmse/%s/%s/0_poleward_of_lat_%g', plotdir, fw, land, par.lat_bound);
                if ~exist(par.folder, 'dir'); mkdir(par.folder); end;

                ymin = -170;
                ymax = 30;

                plot_dmse(dmse.ra_lat.(land).(fw), dmse.res_lat.(land).(fw), dmse.lh_lat.(land).(fw), dmse.sh_lat.(land).(fw), "", ymin, ymax, type, fw, par);
                plot_dmse(dmse.ra_lat.(land).(fw), dmse.res_lat.(land).(fw), dmse.lh_lat.(land).(fw), dmse.sh_lat.(land).(fw), "_noleg", ymin, ymax, type, fw, par);

            end % fw
        end % land
    end % lat bound

end % for function
