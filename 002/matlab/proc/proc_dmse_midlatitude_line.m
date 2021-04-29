function proc_dmse_midlatitude_line(type, par)
    make_dirs(type, par)

    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    % load(sprintf('%s/sftlf.mat', prefix)); % read land fraction data
    load(sprintf('%s/flux_z.mat', prefix_proc)); % load lat x mon RCAE data
    % load(sprintf('%s/%s/flux_t.mat', prefix_proc, par.lat_interp)); % load lat x lon RCAE data
    % landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
    % par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

    % sftlf = nanmean(sftlf, 1); % zonal average
    % sftlf = repmat(sftlf', [1 12]); % expand land fraction data to time
    
    lat_bound_list = [-10 10];
    lat_center = 50;

    for lb = 1:length(lat_bound_list); par.lat_bound = lat_bound_list(lb);
    
        [lat, clat, clat_mon, par] = make_midlatitude_lat(lat_center, par);

        savename = sprintf('%s/dmse_midlatitude_lat_%g_to_%g', prefix_proc, par.lat_center-par.lat_bound, par.lat_center+par.lat_bound);
                
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
            [mesh_lat, mesh_mon] = meshgrid(1:12, lat);

            f_vec = assign_fw(type, par);
            for f = f_vec; fw = f{1};

                dmse.ra = flux_z.(land).ra.(fw);
                dmse.ra_lat.(land).(fw) = interp1(grid.dim3.lat, dmse.ra, lat);
                dmse.ra_lat.(land).(fw) = nansum(dmse.ra_lat.(land).(fw).*clat_mon)/nansum(clat);
                dmse.tend = flux_z.(land).tend;
                dmse.tend_lat.(land).(fw) = interp1(grid.dim3.lat, dmse.tend, lat);
                dmse.tend_lat.(land).(fw) = nansum(dmse.tend_lat.(land).(fw).*clat_mon)/nansum(clat);
                dmse.res = flux_z.(land).res.(fw);
                dmse.res_lat.(land).(fw) = interp1(grid.dim3.lat, dmse.res, lat);
                dmse.res_lat.(land).(fw) = nansum(dmse.res_lat.(land).(fw).*clat_mon)/nansum(clat);

                if strcmp(fw, 'ceresrad')
                    dmse.stf = flux_z.(land).stf.(fw);
                    dmse.stf_lat.(land).(fw) = interp1(grid.dim3.lat, dmse.stf, lat);
                    dmse.stf_lat.(land).(fw) = nansum(dmse.stf_lat.(land).(fw).*clat_mon)/nansum(clat);
                else
                    [lh, sh] = rename_stf(type, flux_z, land);
                    dmse.lh_lat.(land).(fw) = interp1(grid.dim3.lat, lh, lat);
                    dmse.lh_lat.(land).(fw) = nansum(dmse.lh_lat.(land).(fw).*clat_mon)/nansum(clat);
                    dmse.sh_lat.(land).(fw) = interp1(grid.dim3.lat, sh, lat);
                    dmse.sh_lat.(land).(fw) = nansum(dmse.sh_lat.(land).(fw).*clat_mon)/nansum(clat);
                end

            end % for mse dse
        end % for land

        save(savename, 'dmse', '-v7.3')

    end % lat bound

end % function
