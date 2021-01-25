function disp_global_stf(par)
    root = '/project2/tas1/miyawaki/projects/002/data';

    for t = {'erai', 'era5'}; type = t{1};
        load(sprintf('%s/proc/%s/%s/flux_zt.mat', root, type, par.lat_interp));

        for l = {'lo', 'l', 'o'}; land = l{1}; disp(sprintf('---------- %s ----------',land));
            % output global averages
            for fn = {'slhf', 'sshf'}; fname = fn{1};
                disp( sprintf('%s %s is %g Wm^-2', upper(type), upper(fname), nansum(cosd(lat).*flux_zt.(land).ann.(fname)')/nansum(cosd(lat))) );
            end
        end
    end
end
