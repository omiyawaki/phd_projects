function read_div(type, par)
    if strcmp(type, 'era5') | strcmp(type, 'erai') | strcmp(type, 'era5c')
        div_vars=par.era.vars.div;
        div_vars_txt=par.era.vars.div_txt;
        for i=1:length(div_vars)
            % dimensions are (lon x lat x time)
            div.(div_vars_txt{i}) = double(ncread(sprintf('/project2/tas1/miyawaki/projects/002/data/raw/%s/div/%s_div_%s.ymonmean.nc', type, type, par.(type).yr_span), div_vars{i}));
            div.(div_vars_txt{i}) = div.(div_vars_txt{i});
        end
        save(sprintf('/project2/tas1/miyawaki/projects/002/data/read/%s/div.mat', type), 'div', 'div_vars_txt');

    else
        error('Divergence data are only available for ERA.');
    end
end
