function [lat, clat, clat_mon, par] = make_midlatitude_lat(lat_center, par)

        dlat = 0.25; % step size for standard lat grid
        if par.lat_bound > 0 
            par.lat_center=lat_center;
            lat = [-par.lat_bound:dlat:par.lat_bound]+par.lat_center;
            par.shiftby=0;
            if isfield(par, 'monlabelnh')
                par.monlabel=par.monlabelnh;
            end
        else 
            par.lat_center=-lat_center;
            lat = [-par.lat_bound:-dlat:par.lat_bound]+par.lat_center;
            par.shiftby=6;
            if isfield(par, 'monlabelsh')
                par.monlabel=par.monlabelsh;
            end
        end
        clat = cosd(lat); % cosine of latitude for cosine weighting
        clat_mon = repmat(clat', [1 12]);

end