function ocean_mask = remove_ocean(lat, lon, nt)

	[lat2dgrid, lon2dgrid] = meshgrid(lat, lon);

	% ocean_mask = circshift(landmask(lat2dgrid, lon2dgrid, 100), length(lon)/2, 1);
    ocean_vec = ~land_or_ocean(lat2dgrid(:), lon2dgrid(:), 5, 0);
    ocean_mask = reshape(ocean_vec, [size(lat2dgrid)]);
	ocean_mask = double(ocean_mask);
	ocean_mask(ocean_mask==0) = nan;
	ocean_mask = repmat(ocean_mask, [1 1 nt]); % ocean mask in dims (lon x lat x time)

end % compute ocean mask (make ocean nan)
