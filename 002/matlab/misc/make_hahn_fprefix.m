function fprefix = make_hahn_fprefix(par) 

    if contains(par.hahn.clim, 'Control')
        if contains(par.hahn.clim, '1850')
            fprefix = 'ControlSOM1850';
        elseif contains(par.hahn.clim, '2xCO2')
            fprefix = 'ControlSOM2xCO2';
        end
    elseif contains(par.hahn.clim, 'Flat')
        if contains(par.hahn.clim, '1850')
            fprefix = 'FlatSOM1850';
        elseif contains(par.hahn.clim, '2xCO2')
            fprefix = 'FlatSOM2xCO2';
        end
    end
        
end
