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

        tmp = load(sprintf('%s/dmse_poleward_of_lat_%g', prefix_proc, par.lat_bound));
        dmse = tmp.dmse;
        if strcmp(type, 'rea') | (strcmp(type, 'gcm') & strcmp(par.model, 'mmm'))
              dmse_std = tmp.dmse_std;
              dmse_min = tmp.dmse_min;
              dmse_max = tmp.dmse_max;
              dmse_25 = tmp.dmse_25;
              dmse_75 = tmp.dmse_75;
        end
        clear tmp;
        
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
                if strcmp(type, 'erai') & strcmp(fw, 'ceresrad')
                    ymax = 50;
                else
                    ymax = 30;
                end

                if strcmp(type, 'rea') 

                    % shade range
                    plot_dmse(dmse.ra_lat.(land).(fw), dmse.res_lat.(land).(fw), dmse.lh_lat.(land).(fw), dmse.sh_lat.(land).(fw), "", ymin, ymax, type, fw, par, dmse_min.ra_lat.(land).(fw), dmse_min.res_lat.(land).(fw), dmse_min.lh_lat.(land).(fw), dmse_min.sh_lat.(land).(fw), dmse_max.ra_lat.(land).(fw), dmse_max.res_lat.(land).(fw), dmse_max.lh_lat.(land).(fw), dmse_max.sh_lat.(land).(fw));
                    plot_dmse(dmse.ra_lat.(land).(fw), dmse.res_lat.(land).(fw), dmse.lh_lat.(land).(fw), dmse.sh_lat.(land).(fw), "_noleg", ymin, ymax, type, fw, par, dmse_min.ra_lat.(land).(fw), dmse_min.res_lat.(land).(fw), dmse_min.lh_lat.(land).(fw), dmse_min.sh_lat.(land).(fw), dmse_max.ra_lat.(land).(fw), dmse_max.res_lat.(land).(fw), dmse_max.lh_lat.(land).(fw), dmse_max.sh_lat.(land).(fw));

                elseif (strcmp(type, 'gcm') & strcmp(par.model, 'mmm'))

                    % shade IQR
                    plot_dmse(dmse.ra_lat.(land).(fw), dmse.res_lat.(land).(fw), dmse.lh_lat.(land).(fw), dmse.sh_lat.(land).(fw), "", ymin, ymax, type, fw, par, dmse_25.ra_lat.(land).(fw), dmse_25.res_lat.(land).(fw), dmse_25.lh_lat.(land).(fw), dmse_25.sh_lat.(land).(fw), dmse_75.ra_lat.(land).(fw), dmse_75.res_lat.(land).(fw), dmse_75.lh_lat.(land).(fw), dmse_75.sh_lat.(land).(fw));
                    plot_dmse(dmse.ra_lat.(land).(fw), dmse.res_lat.(land).(fw), dmse.lh_lat.(land).(fw), dmse.sh_lat.(land).(fw), "_noleg", ymin, ymax, type, fw, par, dmse_25.ra_lat.(land).(fw), dmse_25.res_lat.(land).(fw), dmse_25.lh_lat.(land).(fw), dmse_25.sh_lat.(land).(fw), dmse_75.ra_lat.(land).(fw), dmse_75.res_lat.(land).(fw), dmse_75.lh_lat.(land).(fw), dmse_75.sh_lat.(land).(fw));

                    % % shade std
                    % plot_dmse(dmse.ra_lat.(land).(fw), dmse.res_lat.(land).(fw), dmse.lh_lat.(land).(fw), dmse.sh_lat.(land).(fw), "", ymin, ymax, type, fw, par, dmse_std.ra_lat.(land).(fw), dmse_std.res_lat.(land).(fw), dmse_std.lh_lat.(land).(fw), dmse_std.sh_lat.(land).(fw));
                    % plot_dmse(dmse.ra_lat.(land).(fw), dmse.res_lat.(land).(fw), dmse.lh_lat.(land).(fw), dmse.sh_lat.(land).(fw), "_noleg", ymin, ymax, type, fw, par, dmse_std.ra_lat.(land).(fw), dmse_std.res_lat.(land).(fw), dmse_std.lh_lat.(land).(fw), dmse_std.sh_lat.(land).(fw));

                else
                    if strcmp(fw, 'ceresrad')
                        plot_dmse(dmse.ra_lat.(land).(fw), dmse.res_lat.(land).(fw), dmse.stf_lat.(land).(fw), dmse.stf_lat.(land).(fw), "", ymin, ymax, type, fw, par);
                        plot_dmse(dmse.ra_lat.(land).(fw), dmse.res_lat.(land).(fw), dmse.stf_lat.(land).(fw), dmse.stf_lat.(land).(fw), "_noleg", ymin, ymax, type, fw, par);
                    else
                        plot_dmse(dmse.ra_lat.(land).(fw), dmse.res_lat.(land).(fw), dmse.lh_lat.(land).(fw), dmse.sh_lat.(land).(fw), "", ymin, ymax, type, fw, par);
                        plot_dmse(dmse.ra_lat.(land).(fw), dmse.res_lat.(land).(fw), dmse.lh_lat.(land).(fw), dmse.sh_lat.(land).(fw), "_noleg", ymin, ymax, type, fw, par);
                    end
                end

            end % fw
        end % land
    end % lat bound

end % for function
