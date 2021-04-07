function [xmin, xmax] = make_lon_bounds(type)

    if strcmp(type, 'merra2')
        xmin = -180;
        xmax = 180;
    else
        xmin = 0;
        xmax = 360;
    end

end