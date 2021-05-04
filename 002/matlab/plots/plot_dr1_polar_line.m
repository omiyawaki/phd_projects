function plot_dr1_polar_line(type, par)
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

    lat_bound_list = [-80 80];

    for lb = 1:length(lat_bound_list); par.lat_bound = lat_bound_list(lb);

        [lat, clat, clat_mon, par] =  make_polar_lat(par);

        load(sprintf('%s/dr1_poleward_of_lat_%g.mat', prefix_proc, par.lat_bound));
        if strcmp(type, 'rea') %| (strcmp(type, 'gcm') & strcmp(par.model, 'mmm'))
            load(sprintf('%s/dr2_poleward_of_lat_%g.mat', prefix_proc, par.lat_bound));
        end

        % for l = {'lo', 'l', 'o'}; land = l{1};
        for l = {'lo'}; land = l{1};
            if strcmp(land, 'lo'); land_text = 'L+O';
            elseif strcmp(land, 'l'); land_text = 'L';
            elseif strcmp(land, 'o'); land_text = 'O';
            end

            if strcmp(type, 'echam')
                if strcmp(par.echam.clim,'20170908'); echamtext='Snowball'; else; echamtext=par.echam.(par.echam.clim); end;
            end
            [mesh_lat, mesh_mon] = meshgrid(1:12, grid.dim2.lat);

            f_vec = assign_fw(type, par);
            for f = f_vec; fw = f{1};

                par.folder = sprintf('%s/dr1/%s/%s/0_poleward_of_lat_%g', plotdir, fw, land, par.lat_bound);
                if ~exist(par.folder, 'dir'); mkdir(par.folder); end;

                % Set y axis limits of plots
                ymin = 0.2;
                ymax = 2.5;

                ymin2 = -0.3;
                ymax2 = 0.8;

                % MAKE PLOTS
                if strcmp(type, 'rea') 

                    % shade range
                    plot_dr1(dr1.r1z_lat.(land).(fw), dr1.r1z_ann_lat.(land).(fw), dr1.dr1z_lat.(land).(fw), dr1.comp1s_lat.(land).(fw), dr1.comp2s_lat.(land).(fw), '', ymin, ymax, type, fw, par, dr1_min.dr1z_lat.(land).(fw), dr1_min.comp1s_lat.(land).(fw), dr1_min.comp2s_lat.(land).(fw), dr1_max.dr1z_lat.(land).(fw), dr1_max.comp1s_lat.(land).(fw), dr1_max.comp2s_lat.(land).(fw));
                    plot_dr1(dr1.r1z_lat.(land).(fw), dr1.r1z_ann_lat.(land).(fw), dr1.dr1z_lat.(land).(fw), dr1.comp1s_lat.(land).(fw), dr1.comp2s_lat.(land).(fw), '_noleg', ymin, ymax, type, fw, par, dr1_min.dr1z_lat.(land).(fw), dr1_min.comp1s_lat.(land).(fw), dr1_min.comp2s_lat.(land).(fw), dr1_max.dr1z_lat.(land).(fw), dr1_max.comp1s_lat.(land).(fw), dr1_max.comp2s_lat.(land).(fw));

                    plot_dr2(dr2.r2z_lat.(land).(fw), dr2.r2z_ann_lat.(land).(fw), dr2.dr2z_lat.(land).(fw), dr2.comp1_lat.(land).(fw), dr2.comp2_lat.(land).(fw), '', ymin2, ymax2, type, fw, par, dr2_min.dr2z_lat.(land).(fw), dr2_min.comp1_lat.(land).(fw), dr2_min.comp2_lat.(land).(fw), dr2_max.dr2z_lat.(land).(fw), dr2_max.comp1_lat.(land).(fw), dr2_max.comp2_lat.(land).(fw));
                    plot_dr2(dr2.r2z_lat.(land).(fw), dr2.r2z_ann_lat.(land).(fw), dr2.dr2z_lat.(land).(fw), dr2.comp1_lat.(land).(fw), dr2.comp2_lat.(land).(fw), '_noleg', ymin2, ymax2, type, fw, par, dr2_min.dr2z_lat.(land).(fw), dr2_min.comp1_lat.(land).(fw), dr2_min.comp2_lat.(land).(fw), dr2_max.dr2z_lat.(land).(fw), dr2_max.comp1_lat.(land).(fw), dr2_max.comp2_lat.(land).(fw));

                elseif (strcmp(type, 'gcm') & strcmp(par.model, 'mmm'))

                    % shade IQR
                    plot_dr1(dr1.r1z_lat.(land).(fw), dr1.r1z_ann_lat.(land).(fw), dr1.dr1z_lat.(land).(fw), dr1.comp1s_lat.(land).(fw), dr1.comp2s_lat.(land).(fw), '', ymin, ymax, type, fw, par, dr1_25.dr1z_lat.(land).(fw), dr1_25.comp1s_lat.(land).(fw), dr1_25.comp2s_lat.(land).(fw), dr1_75.dr1z_lat.(land).(fw), dr1_75.comp1s_lat.(land).(fw), dr1_75.comp2s_lat.(land).(fw));
                    plot_dr1(dr1.r1z_lat.(land).(fw), dr1.r1z_ann_lat.(land).(fw), dr1.dr1z_lat.(land).(fw), dr1.comp1s_lat.(land).(fw), dr1.comp2s_lat.(land).(fw), '_noleg', ymin, ymax, type, fw, par, dr1_25.dr1z_lat.(land).(fw), dr1_25.comp1s_lat.(land).(fw), dr1_25.comp2s_lat.(land).(fw), dr1_75.dr1z_lat.(land).(fw), dr1_75.comp1s_lat.(land).(fw), dr1_75.comp2s_lat.(land).(fw));

                    % plot_dr2(dr2.r2z_lat.(land).(fw), dr2.r2z_ann_lat.(land).(fw), dr2.dr2z_lat.(land).(fw), dr2.comp1_lat.(land).(fw), dr2.comp2_lat.(land).(fw), '', ymin2, ymax2, type, fw, par, dr2_25.dr2z_lat.(land).(fw), dr2_25.comp1_lat.(land).(fw), dr2_25.comp2_lat.(land).(fw), dr2_75.dr2z_lat.(land).(fw), dr2_75.comp1_lat.(land).(fw), dr2_75.comp2_lat.(land).(fw));
                    % plot_dr2(dr2.r2z_lat.(land).(fw), dr2.r2z_ann_lat.(land).(fw), dr2.dr2z_lat.(land).(fw), dr2.comp1_lat.(land).(fw), dr2.comp2_lat.(land).(fw), '_noleg', ymin2, ymax2, type, fw, par, dr2_25.dr2z_lat.(land).(fw), dr2_25.comp1_lat.(land).(fw), dr2_25.comp2_lat.(land).(fw), dr2_75.dr2z_lat.(land).(fw), dr2_75.comp1_lat.(land).(fw), dr2_75.comp2_lat.(land).(fw));

                    % % shade std
                    % % plot_dr1(dr1.r1z_lat.(land).(fw), dr1.r1z_ann_lat.(land).(fw), dr1.dr1z_lat.(land).(fw), dr1.comp1s_lat.(land).(fw), dr1.comp2s_lat.(land).(fw), '', ymin, ymax, type, fw, par, dr1_std.dr1z_lat.(land).(fw))
                    % % plot_dr1(dr1.r1z_lat.(land).(fw), dr1.r1z_ann_lat.(land).(fw), dr1.dr1z_lat.(land).(fw), dr1.comp1s_lat.(land).(fw), dr1.comp2s_lat.(land).(fw), '_noleg', ymin, ymax, type, fw, par, dr1_std.dr1z_lat.(land).(fw))
                    % plot_dr1(dr1.r1z_lat.(land).(fw), dr1.r1z_ann_lat.(land).(fw), dr1.dr1z_lat.(land).(fw), dr1.comp1s_lat.(land).(fw), dr1.comp2s_lat.(land).(fw), '', ymin, ymax, type, fw, par, dr1_std.dr1z_lat.(land).(fw), dr1_std.comp1s_lat.(land).(fw), dr1_std.comp2s_lat.(land).(fw));
                    % plot_dr1(dr1.r1z_lat.(land).(fw), dr1.r1z_ann_lat.(land).(fw), dr1.dr1z_lat.(land).(fw), dr1.comp1s_lat.(land).(fw), dr1.comp2s_lat.(land).(fw), '_noleg', ymin, ymax, type, fw, par, dr1_std.dr1z_lat.(land).(fw), dr1_std.comp1s_lat.(land).(fw), dr1_std.comp2s_lat.(land).(fw));

                    % plot_dr2(dr2.r2z_lat.(land).(fw), dr2.r2z_ann_lat.(land).(fw), dr2.dr2z_lat.(land).(fw), dr2.comp1_lat.(land).(fw), dr2.comp2_lat.(land).(fw), '', ymin-1, ymax-1, type, fw, par, dr2_std.dr2z_lat.(land).(fw), dr2_std.comp1_lat.(land).(fw), dr2_std.comp2_lat.(land).(fw));
                    % plot_dr2(dr2.r2z_lat.(land).(fw), dr2.r2z_ann_lat.(land).(fw), dr2.dr2z_lat.(land).(fw), dr2.comp1_lat.(land).(fw), dr2.comp2_lat.(land).(fw), '_noleg', ymin-1, ymax-1, type, fw, par, dr2_std.dr2z_lat.(land).(fw), dr2_std.comp1_lat.(land).(fw), dr2_std.comp2_lat.(land).(fw));

                else

                    % plot_dr1(dr1.r1z_lat.(land).(fw), dr1.r1z_ann_lat.(land).(fw), dr1.dr1z_lat.(land).(fw), dr1.comp1s_lat.(land).(fw), dr1.comp2s_lat.(land).(fw), '_legonly', ymin, ymax, type, fw, par)

                    plot_dr1(dr1.r1z_lat.(land).(fw), dr1.r1z_ann_lat.(land).(fw), dr1.dr1z_lat.(land).(fw), dr1.comp1s_lat.(land).(fw), dr1.comp2s_lat.(land).(fw), '', ymin, ymax, type, fw, par)
                    plot_dr1(dr1.r1z_lat.(land).(fw), dr1.r1z_ann_lat.(land).(fw), dr1.dr1z_lat.(land).(fw), dr1.comp1s_lat.(land).(fw), dr1.comp2s_lat.(land).(fw), '_noleg', ymin, ymax, type, fw, par)

                    plot_dr2(dr2.r2z_lat.(land).(fw), dr2.r2z_ann_lat.(land).(fw), dr2.dr2z_lat.(land).(fw), dr2.comp1_lat.(land).(fw), dr2.comp2_lat.(land).(fw), '', ymin-1, ymax-1, type, fw, par);
                    plot_dr2(dr2.r2z_lat.(land).(fw), dr2.r2z_ann_lat.(land).(fw), dr2.dr2z_lat.(land).(fw), dr2.comp1_lat.(land).(fw), dr2.comp2_lat.(land).(fw), '_noleg', ymin-1, ymax-1, type, fw, par);
                end

            end % for mse dse
        end % for land
    end % lat bound

end % for function
