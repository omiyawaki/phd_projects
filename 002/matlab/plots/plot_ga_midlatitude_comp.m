function plot_ga_midlatitude_comp(type, par)
    prefix = make_prefix(type, par);
    prefix_proc = make_prefix_proc(type, par);
    plotdir = make_plotdir(type, par);

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/ga_frac_mon_lat.mat', prefix_proc));

    lat_list = [45 -45];
    si_up_list = [0.2 0.3 0.4];
    si_bl_list = [0.9 0.8 0.7];
    
    for f = {'mse_old'}; fw = f{1};
        for l = par.land_list; land = l{1};
            if strcmp(land, 'lo'); land_text = 'Land + Ocean';
            elseif strcmp(land, 'l'); land_text = 'Land';
            elseif strcmp(land, 'o'); land_text = 'Ocean';
            end

            for lb = 1:length(lat_list); lat_eval = lat_list(lb);
                if lat_eval > 0
                    leg_loc = 'north';
                    monlabel = par.monlabelnh;
                    shiftby = 0;
                else
                    leg_loc = 'southwest';
                    monlabel = par.monlabelsh;
                    shiftby = 6;
                end

                folder = sprintf('%s/ga_malr_diff/si_bl_comp/%s/', plotdir, land);
                if ~exist(folder, 'dir')
                    mkdir(folder);
                end

                % interpolate at select latitude
                % res = interp1(grid.dim2.lat, flux_z.(land).res.(fw), lat_eval);
                % ra = interp1(grid.dim2.lat, flux_z.(land).ra.(fw), lat_eval);
                % r1 = res./ra;

                for su = 1:length(si_up_list); par.si_up = si_up_list(su);
                    figure(); clf; hold all; box on;
                    for sb = 1:length(si_bl_list); par.si_bl = si_bl_list(sb);
                        ga = interp1(grid.dim3.lat, ga_frac.(land), lat_eval);
                        ga = permute(ga, [3 1 2]);
                        ga_v = interp1(grid.dim3.si, ga, linspace(par.si_bl, par.si_up, 101));
                        ga_v = squeeze(nanmean(ga_v, 1))';
                        plot(1:12, circshift(ga_v, shiftby));
                    end
                    ylabel(sprintf('$\\left[(\\Gamma_m - \\Gamma)/\\Gamma_m\\right]_{\\sigma_l\\rightarrow\\sigma_u}$ (\\%%)'));
                    axis('tight');
                    make_title_type_lat_pt_si_up(type, lat_eval, par);
                    legend('$\sigma_l=0.9$', '$\sigma_l=0.8$', '$\sigma_l=0.7$', 'location', leg_loc, 'fontsize', 7);
                    set(gcf, 'paperunits', 'inches', 'paperposition', par.ppos)
                    set(gca, 'fontsize', par.fs, 'xlim', [1 12], 'xtick', [1:12], 'xticklabel', monlabel) 
                    print(sprintf('%sga_mon_%g_%g.png', folder, par.si_up, lat_eval), '-dpng', '-r300');
                    close;
                end

            end % lat bounds

        end % land/ocean
    end % framework

end
