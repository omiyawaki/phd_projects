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

        dlat = 0.25; % step size for standard lat grid
        if par.lat_bound>0; par.lat_center=center; lat = [-par.lat_bound:dlat:par.lat_bound]+par.lat_center; par.shiftby=0;
        else; par.lat_center=-center; lat = [-par.lat_bound:-dlat:par.lat_bound]+par.lat_center; par.shiftby=6; end;
        clat = cosd(lat); % cosine of latitude for cosine weighting
        clat_mon = repmat(clat', [1 12]);

        savename = sprintf('%s/dr1_midlatitude_lat_%g_to_%g.mat', prefix_proc, par.lat_center-par.lat_bound, par.lat_center+par.lat_bound);

        for l = par.land_list; land = l{1};
            if strcmp(land, 'lo'); land_text = 'L+O';
            elseif strcmp(land, 'l'); land_text = 'L';
            elseif strcmp(land, 'o'); land_text = 'O';
            end
            if strcmp(type, 'echam')
                if strcmp(par.echam.clim,'20170908'); echamtext='Snowball'; else; echamtext=par.echam.(par.echam.clim); end;
            end

            f_vec = assign_fw(type, par);

            for f = f_vec; fw = f{1};

                % R1 computed before zonal averaging
                %dr1.r1_lat = interp1(grid.dim3.lat, flux_z.(land).r1.(fw), lat);
                %dr1.r1_lat = nansum(dr1.r1_lat.*clat_mon)/nansum(clat);
                %dr1.comp1r_lat = interp1(grid.dim3.lat, flux_z.(land).comp1.(fw), lat);
                %dr1.comp1r_lat = nansum(dr1.comp1r_lat.*clat_mon)/nansum(clat);
                %dr1.comp2r_lat = interp1(grid.dim3.lat, flux_z.(land).comp2.(fw), lat);
                %dr1.comp2r_lat = nansum(dr1.comp2r_lat.*clat_mon)/nansum(clat);

                % R1 computed after zonal averaging
                dr1.r1z_lat.(land).(fw) = interp1(grid.dim3.lat, flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), lat);
                dr1.r1z_lat.(land).(fw) = nansum(dr1.r1z_lat.(land).(fw).*clat_mon)/nansum(clat);

                % R1 computed at each lat x lon RES and RA
                %dr1.r1_ann = repmat(nanmean(flux_z.(land).r1.(fw), 2), [1 12]);
                %dr1.r1_ann_lat = interp1(grid.dim3.lat, dr1.r1_ann, lat);
                %dr1.r1_ann_lat = nansum(dr1.r1_ann_lat.*clat_mon)/nansum(clat);

                % Alternate R1 (computed using zonally averaged RES and RA)
                dr1.r1z_ann.(land).(fw) = repmat(nanmean(flux_z.(land).res.(fw)./flux_z.(land).ra.(fw), 2), [1 12]);
                dr1.r1z_ann_lat.(land).(fw) = interp1(grid.dim3.lat, dr1.r1z_ann.(land).(fw), lat);
                dr1.r1z_ann_lat.(land).(fw) = nansum(dr1.r1z_ann_lat.(land).(fw).*clat_mon)/nansum(clat);

                % annual mean values
                %stf_ann = repmat(nanmean(flux_z.(land).stf.(fw),2), [1 12]);
                %ra_ann = repmat(nanmean(flux_z.(land).ra.(fw),2), [1 12]);
                %fm_ann = repmat(nanmean(flux_z.(land).res.(fw), 2), [1 12]);

                % compute deviation from annual mean
                %dr1 = flux_z.(land).r1.(fw) - r1_ann;
                %ddr1.r1_lat = interp1(grid.dim3.lat, dr1, lat);
                %ddr1.r1_lat = nansum(ddr1.r1_lat.*clat_mon)/nansum(clat);

                dr1.dr1z.(land).(fw) = flux_z.(land).res.(fw)./flux_z.(land).ra.(fw) - dr1.r1z_ann.(land).(fw);
                dr1.dr1z_lat.(land).(fw) = interp1(grid.dim3.lat, dr1.dr1z.(land).(fw), lat);
                dr1.dr1z_lat.(land).(fw) = nansum(dr1.dr1z_lat.(land).(fw).*clat_mon)/nansum(clat);

                % DIVFM and RA DECOMP
                dr1.fm_ann.(land).(fw) = repmat(nanmean(flux_z.(land).res.(fw),2),[1 12]);
                dr1.ra_ann.(land).(fw) = repmat(nanmean(flux_z.(land).ra.(fw),2), [1 12]);
                
                dr1.delta_fm.(land).(fw) = flux_z.(land).res.(fw) - dr1.fm_ann.(land).(fw);
                dr1.comp1s.(land).(fw) = dr1.delta_fm.(land).(fw)./dr1.ra_ann.(land).(fw);
                dr1.comp1s_lat.(land).(fw) = interp1(grid.dim3.lat, dr1.comp1s.(land).(fw), lat);
                dr1.comp1s_lat.(land).(fw) = nansum(dr1.comp1s_lat.(land).(fw).*clat_mon)/nansum(clat);

                dr1.delta_ra.(land).(fw) = flux_z.(land).ra.(fw) - dr1.ra_ann.(land).(fw);
                dr1.ra_ann.(land).(fw) = repmat(nanmean(flux_z.(land).ra.(fw),2),[1 12]);
                dr1.comp2s.(land).(fw) = -dr1.fm_ann.(land).(fw)./(dr1.ra_ann.(land).(fw)).^2.*dr1.delta_ra.(land).(fw);
                dr1.comp2s_lat.(land).(fw) = interp1(grid.dim3.lat, dr1.comp2s.(land).(fw), lat);
                dr1.comp2s_lat.(land).(fw) = nansum(dr1.comp2s_lat.(land).(fw).*clat_mon)/nansum(clat);

                %%%%%%%%%% R2 %%%%%%%%%%%
                dr2.r2z_lat.(land).(fw) = interp1(grid.dim3.lat, flux_z.(land).stf.(fw)./flux_z.(land).ra.(fw), lat);
                dr2.r2z_lat.(land).(fw) = nansum(dr2.r2z_lat.(land).(fw).*clat_mon)/nansum(clat);

                dr2.r2z_ann.(land).(fw) = repmat(nanmean(flux_z.(land).stf.(fw)./flux_z.(land).ra.(fw), 2), [1 12]);
                dr2.r2z_ann_lat.(land).(fw) = interp1(grid.dim3.lat, dr2.r2z_ann.(land).(fw), lat);
                dr2.r2z_ann_lat.(land).(fw) = nansum(dr2.r2z_ann_lat.(land).(fw).*clat_mon)/nansum(clat);

                dr2.dr2z.(land).(fw) = flux_z.(land).stf.(fw)./flux_z.(land).ra.(fw) - dr2.r2z_ann.(land).(fw);
                dr2.dr2z_lat.(land).(fw) = interp1(grid.dim3.lat, dr2.dr2z.(land).(fw), lat);
                dr2.dr2z_lat.(land).(fw) = nansum(dr2.dr2z_lat.(land).(fw).*clat_mon)/nansum(clat);

                dr2.stf_ann.(land).(fw) = repmat(nanmean(flux_z.(land).stf.(fw),2),[1 12]);

                dr2.delta_stf.(land).(fw) = flux_z.(land).stf.(fw) - dr2.stf_ann.(land).(fw);
                dr2.comp1.(land).(fw) = dr2.delta_stf.(land).(fw)./dr1.ra_ann.(land).(fw);
                dr2.comp1_lat.(land).(fw) = interp1(grid.dim3.lat, dr2.comp1.(land).(fw), lat);
                dr2.comp1_lat.(land).(fw) = nansum(dr2.comp1_lat.(land).(fw).*clat_mon)/nansum(clat);

                dr2.comp2.(land).(fw) = -dr2.stf_ann.(land).(fw)./(dr1.ra_ann.(land).(fw)).^2.*dr1.delta_ra.(land).(fw);
                dr2.comp2_lat.(land).(fw) = interp1(grid.dim3.lat, dr2.comp2.(land).(fw), lat);
                dr2.comp2_lat.(land).(fw) = nansum(dr2.comp2_lat.(land).(fw).*clat_mon)/nansum(clat);

            end % for mse dse
        end % for land
        
        save(savename, 'dr1', 'dr2', '-v7.3');
        
    end % lat bounds

end % for function
