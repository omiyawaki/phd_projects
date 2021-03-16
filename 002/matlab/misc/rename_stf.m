function [lh, sh] = rename_stf(type, flux, land, time)

    if nargin < 4

        if any(strcmp(type,{'era5', 'era5c', 'erai'}))
            lh = -flux.(land).slhf;
            sh = -flux.(land).sshf;
        elseif any(strcmp(type,{'hahn'}))
            lh = flux.(land).LHFLX;
            sh = flux.(land).SHFLX;
        elseif any(strcmp(type,{'merra2'}))
            lh = flux.(land).EFLUX;
            sh = flux.(land).HFLUX;
        elseif any(strcmp(type, {'gcm', 'jra55', 'rea'}))
            lh = flux.(land).hfls;
            sh = flux.(land).hfss;
        elseif strcmp(type, 'echam')
            lh = -flux.(land).ahfl;
            sh = -flux.(land).ahfs;
        end

    else

        if any(strcmp(type,{'era5', 'era5c', 'erai'}))
            lh = -flux.(land).(time).slhf;
            sh = -flux.(land).(time).sshf;
        elseif any(strcmp(type,{'hahn'}))
            lh = flux.(land).(time).LHFLX;
            sh = flux.(land).(time).SHFLX;
        elseif any(strcmp(type,{'merra2'}))
            lh = flux.(land).(time).EFLUX;
            sh = flux.(land).(time).HFLUX;
        elseif any(strcmp(type, {'gcm', 'jra55', 'rea'}))
            lh = flux.(land).(time).hfls;
            sh = flux.(land).(time).hfss;
        elseif strcmp(type, 'echam')
            lh = -flux.(land).(time).ahfl;
            sh = -flux.(land).(time).ahfs;
        end

    end

end
