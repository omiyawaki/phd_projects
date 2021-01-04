function read_dondiv00(type, par)
% from 2000-03 to 2018-02
    if strcmp(type, 'erai')

        rootdir = "/project2/tas1/miyawaki/projects/002/data/raw/don/ERA_MHT"; % root directory of Donohoe MSE transport data
        means = load(sprintf('%s/means/1979_10means.mat', rootdir));
        lat = means.lat;
        latr = deg2rad(lat);
        dlat = latr(2)-latr(1);
        clat = cos(latr); clat(1)=nan; clat(end)=nan;

        startdate = datetime(2000,3,1);
        enddate = datetime(2018,2,1);
        dates = datetime(startdate:calmonths(1):enddate, 'format', 'yyyy_M');
        nfiles = length(dates);
        pb=CmdLineProgressBar("Reading Donohoe heat transport data..."); % track progress of this loop
        div_orig = nan([floor(nfiles/13)+1 length(lat) 12]);
        for ifile = 1:nfiles
            pb.print(ifile, nfiles);
            trans_orig = load(sprintf('%s/heat_transport/%s_heattrans.mat', rootdir, dates(ifile)));

            imonth = month(dates(ifile));
            iyear = year(dates(ifile)) - year(startdate) + 1;
            fmse = trans_orig.MME + trans_orig.SE + trans_orig.TE;
            % div_arg = fmse'.*clat;
            % div_orig(year,:,month) = 1./(2*pi*par.a^2*clat.^2).*gradient(div_arg, dlat);
            div_arg = fmse';
            div_orig(iyear,:,imonth) = 1./(2*pi*par.a^2*clat).*gradient(div_arg, dlat);
            % figure; clf; hold all;
            % plot(lat, fmse, 'k');
            % print(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/test_folder/fmse_test', type), '-dpng', '-r300')
            % figure; clf; hold all;
            % plot(lat, div_orig(iyear,:,imonth), 'k');
            % print(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/test_folder/test', type), '-dpng', '-r300')
            % clear fmse div_arg trans_orig
            % return
        end

        % div_orig = circshift(div_orig, -3, 3); % shift months so that January is the first entry (note that Donohoe ERA-I begins on Oct 1979)

        dondiv = squeeze(nanmean(div_orig, 1)); % take climatology
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/dondiv00.mat', type), 'dondiv', 'lat');

    else
        error('Donohoe MSE transport data are available only for ERA-I data.');
    end
end
