function [lat, clat, clat_mon, par] = make_polar_lat(par, pole_lat)

    if nargin == 1
        pole_lat = 90;
    end

        dlat = 0.25; % step size for standard lat grid
        if par.lat_bound > 0
            par.lat_pole = pole_lat;
            lat = par.lat_bound:dlat:par.lat_pole;
            par.shiftby=0;
            if isfield(par, 'monlabelnh')
                par.monlabel=par.monlabelnh;
            end
        else
            par.lat_pole = -pole_lat; 
            lat = par.lat_bound:-dlat:par.lat_pole;
            par.shiftby=6;
            if isfield(par, 'monlabelsh')
                par.monlabel=par.monlabelsh;
            end
        end
        clat = cosd(lat); % cosine of latitude for cosine weighting
        clat_mon = repmat(clat', [1 12]);

end