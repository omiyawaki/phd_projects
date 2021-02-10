function proc_dr1_midlatitude_line(type, par)
    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    % load(sprintf('%s/sftlf.mat', prefix)); % read land fraction data
    load(sprintf('%s/flux_z.mat', prefix_proc)); % load lat x mon RCAE data
    load(sprintf('%s/flux_t.mat', prefix_proc)); % load lat x lon RCAE data
    % landdata = load('/project2/tas1/miyawaki/matlab/landmask/land_mask.mat');
    % par.land = landdata.land_mask; par.landlat = landdata.landlat; par.landlon = landdata.landlon;

    % sftlf = nanmean(sftlf, 1); % zonal average
    % sftlf = repmat(sftlf', [1 12]); % expand land fraction data to time

    par.lat_bound_list = [-10 10];
    center = 50;

    for lb = 1:length(par.lat_bound_list); par.lat_bound = par.lat_bound_list(lb);

        savename = sprintf('%s/dr1_midlatitude_lat_%g_to_%g.mat', prefix_proc, par.lat_center-par.lat_bound, par.lat_center+par.lat_bound);
        if ~exist(par.folder, 'dir'); mkdir(par.folder); end;

        % for l = {'lo', 'l', 'o'}; land = l{1};
        for l = {'lo'}; land = l{1};
            if strcmp(land, 'lo'); land_text = 'L+O';
            elseif strcmp(land, 'l'); land_text = 'L';
            elseif strcmp(land, 'o'); land_text = 'O';
            end
            if strcmp(type, 'echam')
                if strcmp(par.echam.clim,'20170908'); echamtext='Snowball'; else; echamtext=par.echam.(par.echam.clim); end;
            end

            f_vec = assign_fw(type, par);

            for f = f_vec; fw = f{1};

                dlat = 0.25; % step size for standard lat grid
                if par.lat_bound>0; par.lat_center=center; lat = [-par.lat_bound:dlat:par.lat_bound]+par.lat_center; par.shiftby=0; par.monlabel=par.monlabelnh;
                else; par.lat_center=-center; lat = [-par.lat_bound:-dlat:par.lat_bound]+par.lat_center; par.shiftby=6; par.monlabel=par.monlabelsh; end;
                clat = cosd(lat); % cosine of latitude for cosine weighting
                clat_mon = repmat(clat', [1 12]);

                % R1 computed before zonal averaging
                r1_lat = interp1(grid.dim3.lat, flux_z.(land).r1.(fw), lat);
                r1_lat = nansum(r1_lat.*clat_mon)/nansum(clat);
                comp1r_lat = interp1(grid.dim3.lat, flux_z.(land).comp1.(fw), lat);
                comp1r_lat = nansum(comp1r_lat.*clat_mon)/nansum(clat);
                comp2r_lat = interp1(grid.dim3.lat, flux_z.(land).comp2.(fw), lat);
                comp2r_lat = nansum(comp2r_lat.*clat_mon)/nansum(clat);

                % R1 computed after zonal averaging
                r1z_lat = interp1(grid.dim3.lat, flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), lat);
                r1z_lat = nansum(r1z_lat.*clat_mon)/nansum(clat);

                % R1 computed at each lat x lon RES and RA
                r1_ann = repmat(nanmean(flux_z.(land).r1.(fw), 2), [1 12]);
                r1_ann_lat = interp1(grid.dim3.lat, r1_ann, lat);
                r1_ann_lat = nansum(r1_ann_lat.*clat_mon)/nansum(clat);

                % Alternate R1 (computed using zonally averaged RES and RA)
                r1z_ann = repmat(nanmean(flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), 2), [1 12]);
                r1z_ann_lat = interp1(grid.dim3.lat, r1z_ann, lat);
                r1z_ann_lat = nansum(r1z_ann_lat.*clat_mon)/nansum(clat);

                % annual mean values
                stf_ann = repmat(nanmean(flux_z.(land).stf.(fw),2), [1 12]);
                ra_ann = repmat(nanmean(flux_z.(land).ra.(fw),2), [1 12]);
                fm_ann = repmat(nanmean(flux_z.(land).res.(fw), 2), [1 12]);

                % compute deviation from annual mean
                dr1 = flux_z.(land).r1.(fw) - r1_ann;
                dr1_lat = interp1(grid.dim3.lat, dr1, lat);
                dr1_lat = nansum(dr1_lat.*clat_mon)/nansum(clat);

                dr1z = flux_z.(land).res.(fw)./flux_z.(land).ra.(fw) - r1z_ann;
                dr1z_lat = interp1(grid.dim3.lat, dr1z, lat);
                dr1z_lat = nansum(dr1z_lat.*clat_mon)/nansum(clat);

                % DIVFM and RA DECOMP
                ra_ann = repmat(nanmean(flux_z.(land).ra.(fw),2), [1 12]);
                comp1s = delta_fm./ra_ann;
                comp1s_lat = interp1(grid.dim3.lat, comp1s, lat);
                comp1s_lat = nansum(comp1s_lat.*clat_mon)/nansum(clat);

                delta_ra = flux_z.(land).ra.(fw) - repmat(nanmean(flux_z.(land).ra.(fw),2),[1 12]);
                ra_ann = repmat(nanmean(flux_z.(land).ra.(fw),2),[1 12]);
                comp2s = -fm_ann./(ra_ann).^2.*delta_ra;
                comp2s_lat = interp1(grid.dim3.lat, comp2s, lat);
                comp2s_lat = nansum(comp2s_lat.*clat_mon)/nansum(clat);

            end % for mse dse
        end % for land
    end % lat bounds

end % for function
