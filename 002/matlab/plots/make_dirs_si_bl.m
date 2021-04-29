function make_dirs_si_bl(type, par)
    plotdir = make_plotdir(type, par);
    for f = par.(type).fw; fw = f{1};
        for l = par.land_list; land = l{1};
            for t = {'ann', 'djf', 'jja', 'mam', 'son', 'all'}; time = t{1};
                if ~exist(sprintf('%s/ga_malr_diff/si_bl_%g/%s/%s/%s', plotdir, par.si_bl, fw, land, time), 'dir')
                    mkdir(sprintf('%s/ga_malr_diff/si_bl_%g/%s/%s/%s', plotdir, par.si_bl, fw, land, time));
                elseif ~exist(sprintf('%s/dthedpa/si_bl_%g/%s/%s/%s', plotdir, par.si_bl, fw, land, time), 'dir')
                    mkdir(sprintf('%s/dthedpa/si_bl_%g/%s/%s/%s', plotdir, par.si_bl, fw, land, time));
                end
            end
        end
    end
end
