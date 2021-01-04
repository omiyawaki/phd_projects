function make_tempsi_from_ml(type, par)
% model level to sigma temperature
    if strcmp(type, 'echam_ml')
        prefix=sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s', type);
        prefix_proc=sprintf('/project2/tas1/miyawaki/projects/002/data/proc/%s', type);
        var = 't';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/ATM_*.ymonmean.nc'));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ta_orig = double(ncread(fullpath, var));
        var = 'aps';
        file=dir(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/echam_ml/ATM_*.ymonmean.nc'));
        fullpath=sprintf('%s/%s', file.folder, file.name);
        ps_orig = double(ncread(fullpath, var));
    else
        error('This code only works for data output in the model vertical grid.')
    end

    load(sprintf('%s/grid.mat', prefix)); % read grid data
    load(sprintf('%s/srfc.mat', prefix)); % load surface data

    % compute sigma from a and b
    ps_vert = repmat(ps_orig, [1 1 1 size(ta_orig, 3)]); % dims (lon x lat x time x plev)
    ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
    a = permute(repmat(grid.dim3.a, [1 size(ps_orig)]), [2 3 1 4]);
    b = permute(repmat(grid.dim3.b, [1 size(ps_orig)]), [2 3 1 4]);
    si = a./ps_vert + b;

    % interpolate to standard sigma levels
    si = permute(si, [3 1 2 4]); % bring si front
    ta_orig = permute(ta_orig, [3 1 2 4]);
    pb = CmdLineProgressBar("Interpolating temperature to new standard grid...");
    for lo=1:size(si,2)
        pb.print(lo, size(si,2));
        for la=1:size(si,3)
            for mo=1:size(si,4)
                tmp_ta = ta_orig(:,lo,la,mo);
                tmp_ta(end+1) = squeeze(srfc.temp2(lo,la,mo));

                tmp_si = si(:,lo,la,mo);
                tmp_si(end+1) = 1;

                ta_si.lin(:,lo,la,mo) = interp1(tmp_si, tmp_ta, grid.dim3.si, 'linear', nan);
                ta_si.cub(:,lo,la,mo) = interp1(tmp_si, tmp_ta, grid.dim3.si, 'pchip', nan);
                ta_si.spl(:,lo,la,mo) = interp1(tmp_si, tmp_ta, grid.dim3.si, 'spline', nan);
                ta_si.mak(:,lo,la,mo) = interp1(tmp_si, tmp_ta, grid.dim3.si, 'makima', nan);

                clear tmp_ta tmp_si

            end
        end
    end

    ta_si.lin = permute(ta_si.lin, [2 3 1 4]); % reorder to lon x lat x si x mon
    ta_si.cub = permute(ta_si.cub, [2 3 1 4]); % reorder to lon x lat x si x mon
    ta_si.spl = permute(ta_si.spl, [2 3 1 4]); % reorder to lon x lat x si x mon
    ta_si.mak = permute(ta_si.mak, [2 3 1 4]); % reorder to lon x lat x si x mon

    % ta_si = permute(ta_si, [2 3 1 4]); % reorder to lon x lat x si x mon

    if strcmp(type, 'echam_ml'); newdir=sprintf('/project2/tas1/miyawaki/projects/002/data/read/echam_ml'); end;
    if ~exist(newdir, 'dir'); mkdir(newdir); end
    filename='ta_si.mat';
    save(sprintf('%s/%s', newdir, filename), 'ta_si', '-v7.3');
end
