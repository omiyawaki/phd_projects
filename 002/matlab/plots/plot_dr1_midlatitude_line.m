function plot_dr1_midlatitude_line(type, par)
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

    par.lat_bound_list = [-10 10];
    center = 50;

    for lb = 1:length(par.lat_bound_list); par.lat_bound = par.lat_bound_list(lb);
    
        [lat, clat, clat_mon, par] = make_midlatitude_lat(center, par);
        
        load(sprintf('%s/dr1_midlatitude_lat_%g_to_%g.mat', prefix_proc, par.lat_center-par.lat_bound, par.lat_center+par.lat_bound));
        
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
        
                par.folder = sprintf('%s/dr1/%s/%s/0_midlatitude_lat_%g_to_%g', plotdir, fw, land, par.lat_center-par.lat_bound, par.lat_center+par.lat_bound);
                if ~exist(par.folder, 'dir'); mkdir(par.folder); end;

                % Set y axis limits of plots
                if strcmp(type, 'echam')
                    ymin = -0.8;
                    ymax = 0.8;
                else
                    ymin = -0.6;
                    ymax = 0.8;
                end

                % MAKE PLOTS
                if strcmp(type, 'gcm') & strcmp(par.model, 'mmm')
                    plot_dr1(dr1.r1z_lat.(land).(fw), dr1.r1z_ann_lat.(land).(fw), dr1.dr1z_lat.(land).(fw), dr1.comp1s_lat.(land).(fw), dr1.comp2s_lat.(land).(fw), '', ymin, ymax, type, fw, par, dr1_std.dr1z_lat.(land).(fw))
                    plot_dr1(dr1.r1z_lat.(land).(fw), dr1.r1z_ann_lat.(land).(fw), dr1.dr1z_lat.(land).(fw), dr1.comp1s_lat.(land).(fw), dr1.comp2s_lat.(land).(fw), '_noleg', ymin, ymax, type, fw, par, dr1_std.dr1z_lat.(land).(fw))
                else
                    plot_dr1(dr1.r1z_lat.(land).(fw), dr1.r1z_ann_lat.(land).(fw), dr1.dr1z_lat.(land).(fw), dr1.comp1s_lat.(land).(fw), dr1.comp2s_lat.(land).(fw), '', ymin, ymax, type, fw, par)
                    plot_dr1(dr1.r1z_lat.(land).(fw), dr1.r1z_ann_lat.(land).(fw), dr1.dr1z_lat.(land).(fw), dr1.comp1s_lat.(land).(fw), dr1.comp2s_lat.(land).(fw), '_noleg', ymin, ymax, type, fw, par)
                end

            end % for mse dse
        end % for land
    end % lat bound

end % for function
