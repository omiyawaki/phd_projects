function make_dtempml(type, par)
    prefix = make_prefix(type, par);
    newdir = make_savedir(type, par);

    % tmp = load(sprintf('%s/grid.mat', prefix)); grid=tmp.grid; clear tmp; % read grid data
    tmp = load(sprintf('%s/ta_si.mat', prefix)); t0 = tmp.ta_si; clear tmp; % load control climate temperature
    if strcmp(par.gcm.clim, 'historical')
        par2 = par;
        par2.gcm.clim = 'rcp85';
        par2.gcm.yr_span = '207001-209912';
        prefix2 = make_prefix(type, par2);
        tmp = load(sprintf('%s/ta_si.mat', prefix2)); tf = tmp.ta_si; clear tmp; % load perturbed climate temperature
    end
    
    dta_si.spl = tf.spl - t0.spl;

    filename='dta_si.mat';
    save(sprintf('%s/%s', newdir, filename), 'dta_si', '-v7.3');
end
