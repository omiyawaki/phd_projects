% % PICONTROL
%par.gcm.clim = 'piControl'; % choose either piControl, historical, or abrupt4xCO2
%par.lat_interp = 'native'; % native: native model grid, don: Donohoe grid, ERA: native ERA grid, std: defined high resolution grid
%par.gcm_models = {'ACCESS1-0', 'ACCESS1-3',... % GCM model names
%                  'bcc-csm1-1', 'bcc-csm1-1-m',...
%                  'BNU-ESM', 'CanESM2', 'CCSM4',...
%                  'CNRM-CM5', 'CNRM-CM5-2',...
%                  'CSIRO-Mk3-6-0', 'FGOALS-g2', 'FGOALS-s2',...
%                  'GFDL-CM3', 'GFDL-ESM2G', 'GFDL-ESM2M',...
%                  'GISS-E2-H', 'GISS-E2-R',...
%                  'HadGEM2-ES', 'inmcm4',...
%                  'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'IPSL-CM5B-LR',...
%                  'MIROC5', 'MIROC-ESM',...
%                  'MPI-ESM-LR', 'MPI-ESM-MR', 'MPI-ESM-P',...
%                  'MRI-CGCM3', 'NorESM1-M'};

% HISTORICAL
par.gcm.clim = 'historical'; % choose either piControl, historical, or abrupt4xCO2
par.lat_interp = 'native'; % native: native model grid, don: Donohoe grid, ERA: native ERA grid, std: defined high resolution grid
%par.gcm_models = {'MPI-ESM-LR'};

% par.gcm_models = {'ACCESS1-0', 'ACCESS1-3',... % GCM model names
%                   'bcc-csm1-1', 'bcc-csm1-1-m',...
%                   'BNU-ESM', 'CanESM2', 'CCSM4',...
%                   'CESM1-BGC', 'CESM1-CAM5', 'CESM1-WACCM',...
%                   'CMCC-CESM', 'CMCC-CM',...
%                   'CNRM-CM5', 'CNRM-CM5-2',...
%                   'CSIRO-Mk3-6-0', 'FGOALS-g2', 'FGOALS-s2',...
%                   'GFDL-CM3', 'GFDL-ESM2G', 'GFDL-ESM2M',...
%                   'GISS-E2-H', 'GISS-E2-H-CC', 'GISS-E2-R', 'GISS-E2-R-CC',...
%                   'HadCM3', 'HadGEM2-CC', 'HadGEM2-ES', 'inmcm4',...
%                   'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'IPSL-CM5B-LR',...
%                   'MIROC5', 'MIROC-ESM', 'MIROC-ESM-CHEM',...
%                   'MPI-ESM-LR', 'MPI-ESM-MR', 'MPI-ESM-P',...
%                   'MRI-CGCM3', 'MRI-ESM1', 'NorESM1-M', 'NorESM1-ME'};

% subset with RCP85

par.gcm_models = {'ACCESS1-0', 'ACCESS1-3',... % GCM model names
                 'bcc-csm1-1', 'bcc-csm1-1-m',...
                 'BNU-ESM', 'CanESM2', 'CCSM4',...
                 'CESM1-BGC', 'CESM1-CAM5',...
                 'CMCC-CESM', 'CMCC-CM',...
                 'CNRM-CM5',...
                 'CSIRO-Mk3-6-0', 'FGOALS-g2',...
                 'GFDL-CM3', 'GFDL-ESM2G', 'GFDL-ESM2M',...
                 'GISS-E2-H', 'GISS-E2-H-CC', 'GISS-E2-R', 'GISS-E2-R-CC',...
                 'HadGEM2-CC', 'HadGEM2-ES', 'inmcm4',...
                 'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'IPSL-CM5B-LR',...
                 'MIROC5', 'MIROC-ESM', 'MIROC-ESM-CHEM',...
                 'MPI-ESM-LR', 'MPI-ESM-MR',...
                 'MRI-CGCM3', 'MRI-ESM1', 'NorESM1-M', 'NorESM1-ME'};


% % abrupt4xCO2
% par.gcm.clim = 'abrupt4xCO2'; % choose either piControl, historical, or abrupt4xCO2
% par.lat_interp = 'native'; % native: native model grid, don: Donohoe grid, ERA: native ERA grid, std: defined high resolution grid
% par.gcm_models = {'ACCESS1-0', 'ACCESS1-3',... % GCM model names
%                   'bcc-csm1-1', 'bcc-csm1-1-m',...
%                   'BNU-ESM', 'CanESM2', 'CCSM4',...
%                   'CNRM-CM5', 'CNRM-CM5-2',...
%                   'CSIRO-Mk3-6-0', 'FGOALS-g2', 'FGOALS-s2',...
%                   'GFDL-CM3', 'GFDL-ESM2G', 'GFDL-ESM2M',...
%                   'GISS-E2-H', 'GISS-E2-R',...
%                   'HadGEM2-ES', 'inmcm4',...
%                   'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'IPSL-CM5B-LR',...
%                   'MIROC5', 'MIROC-ESM',...
%                   'MPI-ESM-LR', 'MPI-ESM-MR', 'MPI-ESM-P',...
%                   'MRI-CGCM3', 'NorESM1-M'};

par.vars_gcm_2d = {'rsdt', 'rsut', 'rsus', 'rsds', 'rlus', 'rlds', 'rlut', 'hfls', 'hfss', 'ps', 'hurs'}; % 2D GCM variables to read
par.vars_gcm_3d = {'ta', 'hur', 'zg', 'ua', 'va'}; % 3d GCM variables
