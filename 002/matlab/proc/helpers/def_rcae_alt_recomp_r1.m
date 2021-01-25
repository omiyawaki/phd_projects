function rcae = def_rcae_alt_recomp_r1(type, flux, vh, par) % recalculate R1 at this step
    f_vec = assign_fw(type, par);
    for f = f_vec; fw = f{1};
        % calculate R1 again to test importance of order of operations
        if strcmp(fw, 'mse2'); r1 = flux.res.(fw)./flux.lw;
        else r1 = flux.res.(fw)./flux.ra.(fw); end

        % identify locations of RCE using threshold epsilon (ep)
        rcae.(fw).def = zeros(size(r1));
        rcae.(fw).def(r1 < par.ep) = 1;
        % identify locations of RAE using threshold gamma (ga)
        if strcmp(fw, 'dse'); rcae.(fw).def(flux.r1.mse > par.ga) = -1; % always use MSE for RAE
        else rcae.(fw).def(r1 > par.ga) = -1; end;

        % identify locations of RCE following Jakob et al. (2019) (abs(TEDIV) < 50 W/m^2)
        rcae.(fw).jak = zeros(size(r1));
        rcae.(fw).jak(flux.res.(fw) < 50) = 1;
        % identify locations of RAE using threshold gamma (ga)
        if strcmp(fw, 'dse'); rcae.(fw).jak(flux.r1.mse > par.ga) = -1; % always use MSE for RAE
        else rcae.(fw).jak(r1 > par.ga) = -1; end;

        % identify locations of rce following jakob et al. (2019) (abs(tediv) < 30 w/m^2)
        rcae.(fw).jak30 = zeros(size(r1));
        rcae.(fw).jak30(flux.res.(fw) < 30) = 1;
        % identify locations of rae using threshold gamma (ga)
        if strcmp(fw, 'dse'); rcae.(fw).jak30(flux.r1.mse > par.ga) = -1; % always use mse for rae
        else rcae.(fw).jak30(r1 > par.ga) = -1; end;

        % identify locations of rce following jakob et al. (2019) (abs(tediv) < 10 w/m^2)
        rcae.(fw).jak10 = zeros(size(r1));
        rcae.(fw).jak10(flux.res.(fw) < 10) = 1;
        % identify locations of rae using threshold gamma (ga)
        if strcmp(fw, 'dse'); rcae.(fw).jak10(flux.r1.mse > par.ga) = -1; % always use mse for rae
        else rcae.(fw).jak10(r1 > par.ga) = -1; end;

        % % add additional criteria for RCE that P-E>0
        % rcae.(fw).pe = zeros(size(r1));
        % if strcmp(type, 'era5') | strcmp(type, 'erai')
        %     rcae.(fw).pe(r1 < par.ep & (flux.cp+flux.lsp+flux.e>0)) = 1; % note that evap is defined negative into atmosphere
        % elseif strcmp(type, 'gcm')
        %     rcae.(fw).pe(r1 < par.ep & (flux.pr-flux.evspsbl>0)) = 1;
        % elseif strcmp(type, 'echam')
        %     rcae.(fw).pe(r1 < par.ep & (flux.aprc + flux.aprl -flux.evap>0)) = 1;
        % end
        % if strcmp(fw, 'dse'); rcae.(fw).pe(flux.r1.mse > par.ga) = -1;
        % else rcae.(fw).pe(r1 > par.ga) = -1; end;

        % % add additional criteria for RCE that (P_ls - E)<<1
        % rcae.(fw).cp = zeros(size(r1));
        % if strcmp(type, 'era5') | strcmp(type, 'erai')
        %     rcae.(fw).cp(r1 < par.ep & abs(flux.lsp./flux.cp<par.ep_cp)) = 1; % note that evap is defined negative into atmosphere
        % elseif strcmp(type, 'gcm')
        %     rcae.(fw).cp(r1 < par.ep & abs((flux.pr-flux.prc)./flux.prc<par.ep_cp)) = 1; % note that evap is defined negative into atmosphere
        % elseif strcmp(type, 'echam')
        %     rcae.(fw).cp(r1 < par.ep & abs(flux.aprl./flux.aprc<par.ep_cp)) = 1; % note that evap is defined negative into atmosphere
        % end
        % if strcmp(fw, 'dse') rcae.(fw).cp(flux.r1.mse > par.ga) = -1;
        % else rcae.(fw).cp(r1 > par.ga) = -1; end;

        % % add additional criteria for RCE that w500<0 (ascent)
        % rcae.(fw).w500 = zeros(size(r1));
        % if strcmp(type, 'era5') | strcmp(type, 'erai')
        %     rcae.(fw).w500(r1 < par.ep & flux.w500 < 0) = 1; % note that evap is defined negative into atmosphere
        % elseif strcmp(type, 'gcm')
        %     rcae.(fw).w500(r1 < par.ep & flux.w500 < 0) = 1; % note that evap is defined negative into atmosphere
        % end
        % if strcmp(fw, 'dse'); rcae.(fw).w500(flux.r1.mse > par.ga) = -1;
        % else rcae.(fw).w500(r1 > par.ga) = -1; end;

        % % add additional criteria for RCE that vh<max(vh)/2 (weak meridional velocity)
        % rcae.(fw).vas2 = zeros(size(r1));
        % if strcmp(type, 'era5') | strcmp(type, 'erai')
        %     rcae.(fw).vas2(r1 < par.ep & abs(flux.vas) < nanmax(nanmax(abs(flux.vas)))/2) = 1; % note that evap is defined negative into atmosphere
        % elseif strcmp(type, 'gcm')
        %     rcae.(fw).vas2(r1 < par.ep & flux.vas < nanmax(nanmax(abs(flux.vas)))/2) = 1; % note that evap is defined negative into atmosphere
        % end
        % if strcmp(fw, 'dse'); rcae.(fw).vas2(flux.r1.mse > par.ga) = -1;
        % else rcae.(fw).vas2(r1 > par.ga) = -1; end;

        % % add additional criteria for RCE that vh<max(vh)/4 (weak meridional velocity)
        % rcae.(fw).vas4 = zeros(size(r1));
        % if strcmp(type, 'era5') | strcmp(type, 'erai')
        %     rcae.(fw).vas4(r1 < par.ep & abs(flux.vas) < nanmax(nanmax(abs(flux.vas)))/4) = 1; % note that evap is defined negative into atmosphere
        % elseif strcmp(type, 'gcm')
        %     rcae.(fw).vas4(r1 < par.ep & flux.vas < nanmax(nanmax(abs(flux.vas)))/4) = 1; % note that evap is defined negative into atmosphere
        % end
        % if strcmp(fw, 'dse'); rcae.(fw).vas4(flux.r1.mse > par.ga) = -1;
        % else rcae.(fw).vas4(r1 > par.ga) = -1; end;

        if size(r1, 1)==length(par.lat_std) & size(r1, 2)==12 % this criteria only works for zonally and time averaged data because vh is only a function of lat
            % add additional criteria for RCE that horizontal transport is weak
            rcae.(fw).vh2 = zeros(size(r1));
            % if strcmp(type, 'era5') | strcmp(type, 'erai')
                rcae.(fw).vh2(r1 < par.ep & abs(vh.(fw)) < nanmax(nanmax(abs(vh.(fw))))/2 ) = 1;
            % elseif strcmp(type, 'gcm')
            %     rcae.(fw).vh2(r1 < par.ep & abs(vh.(fw)) < nanmax(nanmax(abs(vh.(fw))))/2 ) = 1;
            % end
            if strcmp(fw, 'dse'); rcae.(fw).vh2(flux.r1.mse > par.ga) = -1;
            else; rcae.(fw).vh2(r1 > par.ga) = -1; end;

            rcae.(fw).vh3 = zeros(size(r1));
            % if strcmp(type, 'era5') | strcmp(type, 'erai')
                rcae.(fw).vh3(r1 < par.ep & abs(vh.(fw)) < nanmax(nanmax(abs(vh.(fw))))/3 ) = 1;
            % elseif strcmp(type, 'gcm')
            %     rcae.(fw).vh3(r1 < par.ep & abs(vh.(fw)) < nanmax(nanmax(abs(vh.(fw))))/3 ) = 1;
            % end
            if strcmp(fw, 'dse') rcae.(fw).vh3(flux.r1.mse > par.ga) = -1;
            else rcae.(fw).vh3(r1 > par.ga) = -1; end;

            rcae.(fw).vh4 = zeros(size(r1));
            % if strcmp(type, 'era5') | strcmp(type, 'erai')
                rcae.(fw).vh4(r1 < par.ep & abs(vh.(fw)) < nanmax(nanmax(abs(vh.(fw))))/4 ) = 1;
            % elseif strcmp(type, 'gcm')
            %     rcae.(fw).vh4(r1 < par.ep & abs(vh.(fw)) < nanmax(nanmax(abs(vh.(fw))))/4 ) = 1;
            % end
            if strcmp(fw, 'dse'); rcae.(fw).vh4(flux.r1.mse > par.ga) = -1;
            else; rcae.(fw).vh4(r1 > par.ga) = -1; end;
        end

    end
end % define RCE and RAE
