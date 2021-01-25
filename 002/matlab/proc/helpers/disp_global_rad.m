function disp_global_rad(par)
    root = '/project2/tas1/miyawaki/projects/002/data';
    load(sprintf('%s/proc/comp/comp_zt', root));

    % output global averages
    for fn = {'tsr', 'ssr', 'ttr', 'str'}; fname = fn{1};
        disp( sprintf('CERES %s is %g Wm^-2', upper(fname), nansum(cosd(lat).*ceres_zt.rad.(fname))/nansum(cosd(lat))) );
        disp( sprintf('ERA-I %s is %g Wm^-2', upper(fname), nansum(cosd(lat).*erai_zt.rad.(fname))/nansum(cosd(lat))) );
        disp( sprintf('ERA5 %s is %g Wm^-2', upper(fname), nansum(cosd(lat).*era5_zt.rad.(fname))/nansum(cosd(lat))) );
    end
end
