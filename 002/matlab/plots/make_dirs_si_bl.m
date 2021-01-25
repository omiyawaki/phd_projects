function make_dirs_si_bl(type, par)
    plotdir = make_plotdir(type, par);
    for l = {'lo', 'l', 'o'}; land = l{1};
        for t = {'ann', 'djf', 'jja', 'mam', 'son', 'all'}; time = t{1};
            if ~exist(sprintf('%s/ga_malr_diff/si_bl_%g/%s/%s', plotdir, par.si_bl, land, time), 'dir')
                mkdir(sprintf('%s/ga_malr_diff/si_bl_%g/%s/%s', plotdir, par.si_bl, land, time));
            end
        end
    end
end
