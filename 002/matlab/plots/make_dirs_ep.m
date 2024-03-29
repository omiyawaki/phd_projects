function make_dirs_ep(type, par)
    plotdir = make_plotdir(type, par);
    % make figure directories if it does not exist
    for i = 1:length(par.ep_swp); par.ep = par.ep_swp(i);
        f_vec = assign_fw(type, par);
        for f = f_vec; fw = f{1};
            for j = {'def', 'jak', 'jak30', 'jak10', 'pe', 'cp', 'w500', 'vh2', 'vh3', 'vh4'}; crit = j{1};
                for l = {'lo', 'l', 'o'}; land = l{1};
                    for k = {'ann', 'djf', 'jja', 'mam', 'son'}; time = k{1};
                        if ~exist(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp_%s', plotdir, par.ep, par.ga, fw, crit, land, time, num2str(par.ma_init)), 'dir')
                            mkdir(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/temp_%s', plotdir, par.ep, par.ga, fw, crit, land, time, num2str(par.ma_init)));
                        elseif ~exist(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/ga_fr', plotdir, par.ep, par.ga, fw, crit, land, time), 'dir')
                            mkdir(sprintf('%s/eps_%g_ga_%g/%s/%s/%s/%s/ga_fr', plotdir, par.ep, par.ga, fw, crit, land, time));
                        end
                    end
                end
            end
        end
    end
end
