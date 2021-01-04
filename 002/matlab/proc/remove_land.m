function land_mask = remove_land(lat, lon, nt)

    [lat2dgrid, lon2dgrid] = meshgrid(lat, lon);

    % land_mask =~ circshift(landmask(lat2dgrid, lon2dgrid, 100), length(lon)/2, 1);
    land_vec = land_or_ocean(lat2dgrid(:), lon2dgrid(:), 5, 0);
    land_mask = reshape(land_vec, [size(lat2dgrid)]);
    land_mask = double(land_mask);
    land_mask(land_mask==0) = nan;
    land_mask = repmat(land_mask, [1 1 nt]); % land mask in dims (lon x lat x time)

end % compute land mask (make land nan)
