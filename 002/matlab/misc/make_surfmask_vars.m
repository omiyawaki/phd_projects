function [ps_vert, pa] = make_surfmask_vars(grid, type, srfc, par) 

    if nargin < 4 % use grid pressure levels
        plev = grid.dim3.plev;
    else % use par.pa pressure levels (e.g. for moist adiabatic lapse rate calculation
        plev = par.pa';
    end

    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        ps_vert = repmat(srfc.sp, [1 1 1 length(plev)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(plev, [1 size(srfc.sp)]), [2 3 1 4]));
    elseif any(strcmp(type, {'merra2', 'merra2c'}))
        ps_vert = repmat(srfc.PS, [1 1 1 length(plev)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(plev, [1 size(srfc.PS)]), [2 3 1 4]));
    elseif strcmp(type, 'jra55')
        ps_vert = repmat(srfc.ps, [1 1 1 length(plev)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = double(permute(repmat(plev, [1 size(srfc.ps)]), [2 3 1 4]));
    elseif strcmp(type, 'gcm')
        ps_vert = repmat(srfc.ps, [1 1 1 length(plev)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(plev, [1 size(srfc.ps)]), [2 3 1 4]);
    elseif contains(type, 'echam')
        ps_vert = repmat(srfc.aps, [1 1 1 length(plev)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(plev, [1 size(srfc.aps)]), [2 3 1 4]);
    elseif contains(type, 'hahn')
        ps_vert = repmat(srfc.PS, [1 1 1 length(plev)]); % dims (lon x lat x time x plev)
        ps_vert = permute(ps_vert, [1 2 4 3]); % dims (lon x lat x plev x time)
        pa = permute(repmat(plev, [1 size(srfc.PS)]), [2 3 1 4]);
    end
        
end
