function make_dmasi(type, par)
    prefix = make_prefix(type, par);
    newdir = make_savedir(type, par);

    % tmp = load(sprintf('%s/grid.mat', prefix)); grid=tmp.grid; clear tmp; % read grid data
    tmp = load(sprintf('%s/ma_si_%s_%s.mat', prefix, par.ma_init, par.ma_type)); t0 = tmp.ma_si; clear tmp; % load control climate temperature
    if strcmp(par.gcm.clim, 'historical')
        par2 = par;
        par2.gcm.clim = 'rcp85';
        par2.gcm.yr_span = '207001-209912';
        prefix2 = make_prefix(type, par2);
        tmp = load(sprintf('%s/ma_si_%s_%s.mat', prefix2, par.ma_init, par.ma_type)); tf = tmp.ma_si; clear tmp; % load perturbed climate temperature
    end
    
    dma_si = tf - t0;

    filename='dma_si.mat';
    save(sprintf('%s/%s', newdir, filename), 'dma_si', '-v7.3');
end
