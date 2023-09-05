%% Read in Vocals Rex Data


% By Andrew John Buggee
%%

function vocalsRex = readVocalsRex(filename)

info = ncinfo(filename);

% read the time vector
time = double(ncread(filename, 'Time'));                   % Measured in seconds

% make sure the time data is a row vector
time = reshape(time, 1, length(time));

% [hours, minutes]
startTime = [str2double(info.Variables(1).Attributes(3).Value(26:27)),...
    str2double(info.Variables(1).Attributes(3).Value(29:30))];      % Start time of the data log in UTC

% ----------------------------------------------------------------
% ----------- cloud droplet probe measurements (CDP) -------------
% ----------------------------------------------------------------

% What we measure IS the number concentration for each bin! we do NOT measure the differential dN/dr

num_concentration_CDP = ncread(filename, 'CCDP_RWO');          % #/cm^3 - number concentration for each bin

% The droplet number concentration is filtered into 30 bins from 2 to 52
% microns
% Below 20 microns, the bins are spaced by 1.18 microns
% above 22 microns, the bins are spaced by 2.28 microns
% find the variable information for the droplet concentration data
CCDP_RWO_info = ncinfo(filename, 'CCDP_RWO');

% Bins are defined as diameters. Divide by 2 to get the radius
drop_radius_bin_edges_CDP = double(CCDP_RWO_info.Attributes(10).Value./2);                        % microns

% The data tells us which bins to take
first_bin_CDP = CCDP_RWO_info.Attributes(8).Value;
last_bin_CDP = CCDP_RWO_info.Attributes(9).Value;

% There are 2 edges for each measurement
%
drop_radius_bin_edges_CDP = drop_radius_bin_edges_CDP(first_bin_CDP:last_bin_CDP+1);                         % microns
drop_radius_bin_center_CDP = drop_radius_bin_edges_CDP(1:end-1) + diff(drop_radius_bin_edges_CDP)/2;         % microns

%drop_size_dist_CDP2 = drop_size_dist_CDP;
num_concentration_CDP = num_concentration_CDP(first_bin_CDP:last_bin_CDP,1,:);          % #/cm^3



% ----------------------------------------------------------------
% ---------- Two-Dimensional Optical array probe (2DC) -----------
% ----------------------------------------------------------------

% This droplet probe filters droplets into 64 different size bins. The bins
% range from 25 microns to 1560 microns.

% ** IMPORTANT ** The first bin of the 2DC probe overlaps with the last bin
% of the CDP. So we ignore it. (See Painemal and Zuidema paragraph 15.)

% There have been some issues with naming the 2DC data different variables.
% Let's loop through all the variables and make sure we grab the correct
% one

for nn = 1:length(info.Variables)

    if strcmp('C1DCA_RPC', info.Variables(nn).Name)==true

        num_concentration_2DC = ncread(filename, 'C1DCA_RPC');        % #/L  - where I think L stands for liters. 1 liter = 1000 cubic cm - again we only measure the number concentration!

        % convert the 2DC data to inverse cubic centimeters
        num_concentration_2DC = num_concentration_2DC ./ 1000;         % # of droplets/cm^3 per bin

        % read in the 2DC variable info
        C1DCA_RPC_info = ncinfo(filename, 'C1DCA_RPC');

        % Grab the droplet size bins
        % These are the boundaries that define each bin size
        % Bins are defined as diameters. Divide by 2 to get the radius
        drop_radius_bin_edges_2DC = C1DCA_RPC_info.Attributes(12).Value./2;                        % microns

        % The data tells us which bins to take
        for jj = 1:length(C1DCA_RPC_info.Attributes)
            
            if strcmp(C1DCA_RPC_info.Attributes(jj).Name, 'FirstBin')==true
                first_bin_2DC = double(C1DCA_RPC_info.Attributes(jj).Value);
            end

            if strcmp(C1DCA_RPC_info.Attributes(jj).Name, 'LastBin')==true
                last_bin_2DC = double(C1DCA_RPC_info.Attributes(jj).Value);
            end


        end

        % Sometimes there isn't a frist or last bin. In that case, set them
        % to be the full breadth of the data
        if exist("first_bin_2DC", 'var')==false
            % set it to be 1
            first_bin_2DC = 1;
        end

        if exist("last_bin_2DC", 'var')==false
            % set it to be the length of drop_bin_edges
            last_bin_2DC = length(drop_radius_bin_edges_2DC)-1;
        end
        

        drop_radius_bin_edges_2DC = drop_radius_bin_edges_2DC(first_bin_2DC:last_bin_2DC+1);                          % microns
        drop_radius_bin_center_2DC = drop_radius_bin_edges_2DC(1:end-1) + diff(drop_radius_bin_edges_2DC)/2;          % microns

        num_concentration_2DC = num_concentration_2DC(first_bin_2DC:last_bin_2DC,1,:);



    end


end

% if after the entire for loop this variable is not found, we will ingore
% this data set. But we have to introduce some variables so the rest of the
% code doesn't have to change


if exist("num_concentration_2DC", 'var')==false


    for nn = 1:length(info.Variables)

        if strcmp('C1DC_RPC', info.Variables(nn).Name)==true


            % read in the 2DC variable info
            C1DC_RPC_info = ncinfo(filename, 'C1DC_RPC');

            % Grab the droplet size bins
            % These are the boundaries that define each bin size
            % Bins are defined as diameters. Divide by 2 to get the radius
            drop_radius_bin_edges_2DC = C1DC_RPC_info.Attributes(12).Value./2;                        % microns

            % The data tells us which bins to take
            first_bin_2DC = double(C1DC_RPC_info.Attributes(10).Value);
            last_bin_2DC = double(C1DC_RPC_info.Attributes(11).Value);

            drop_radius_bin_edges_2DC = drop_radius_bin_edges_2DC(first_bin_2DC:last_bin_2DC+1);                          % microns
            drop_radius_bin_center_2DC = drop_radius_bin_edges_2DC(1:end-1) + diff(drop_radius_bin_edges_2DC)/2;          % microns

            % ****** SET THE NUMBER CONCENTRATION TO ZERO FOR ALL BINS ******
            num_concentration_2DC = zeros((last_bin_2DC - first_bin_2DC +1),1, length(time));



        end


    end


end


% -------------------------------------------------------------------------
% ------- Grab Wind speed, direction and sensor lat,long, altititude ------
% -------------------------------------------------------------------------

% grab the horizontal wind speed and direction
horz_wind_speed = ncread(filename, 'WSC');                                    % m/s -  GPS corrected horizontal wind speed
horz_wind_direction = ncread(filename, 'WDC');                                 % degrees from north -  GPS corrected horizontal wind direction 'wind from'

% Grab position and altitude data
lat = ncread(filename, 'LAT');                                                        % Meausred in degrees North
long = ncread(filename, 'LON');                                                       % Meausred in degrees East
altitude = ncread(filename, 'ALTX');                                                  % Measured in meters above sea level


% Somestimes the altitude vector is not a vector at all but a matrix.
% Usually each row is roughly the same. Just take the first row
if numel(altitude)~=numel(time)
    altitude = altitude(1,:);
end

% Make sure that altitude is a row vector
altitude = reshape(altitude, 1, length(altitude));



% -------------------------------------------------------------------------
% ------- Grab the ambient air temperature and water vapor pressure ------
% -------------------------------------------------------------------------
water_vapor_pressure = ncread(filename, 'EDPUV');                                    % hPa - ambient water vapor pressure
ambient_air_temp = ncread(filename, 'ATX');                                    % C - ambient air tempertaure 



% -------------------------------------------------------------------------
% ------- Grab the short and longwave radiances looking up and down ------
% -------------------------------------------------------------------------

% lets grab the shortwave and longwave radiance looking down and looking up
shortwave_top = ncread(filename, 'SWT');
shortwave_bot = ncread(filename, 'SWB');

longwave_top = ncread(filename,'IRTC');              % These are corrected longwave irradiance values
longwave_bot = ncread(filename, 'IRBC');






% -------------------------------------------------------------------
% ----------------- ***** IMPORTANT STEP ***** ----------------------
% -------------------------------------------------------------------


% Combine data from both droplet probes to get the total droplet size
% distribution

% ***** There is a bin where no data exists. The CDP data ends at bin
% [24.0550 - 25.195] microns. The 2DC data starts with the bin
% [31.25 - 43.75] microns. Therefore we have to place a 0 inbetween these
% two data bins. That is, since we have no measurement, we define the
% number concentration to be 0 between [25.195 - 31.25] *****

% we also need to define the center point of the bin that has 0
% measurements
center_for_0 = drop_radius_bin_edges_CDP(end) + (drop_radius_bin_edges_2DC(1) - drop_radius_bin_edges_CDP(end))/2;                               % microns
drop_radius_bin_center = [drop_radius_bin_center_CDP, center_for_0, drop_radius_bin_center_2DC];                    % microns                                                   % microns
drop_radius_bin_edges = [drop_radius_bin_edges_CDP, drop_radius_bin_edges_2DC];                                                             % microns
%drop_radius_bin_edges2 = [drop_radius_bin_edges_CDP2, drop_radius_bin_edges_2DC2];                                                             % microns
%drop_radius_bin_center2 = drop_radius_bin_edges2(1:end-1) + diff(drop_radius_bin_edges2)./2;                                                  % microns
%drop_radius_bin_first_edge = [drop_radius_bin_edges_CDP(1:end-1), drop_radius_bin_edges_2DC(1:end-1)];                                      % microns

% don't forget the 0!
% Nc is measured in #/cm^3
Nc = [reshape(num_concentration_CDP,length(drop_radius_bin_center_CDP),[]); zeros(1,length(time)); reshape(num_concentration_2DC,length(drop_radius_bin_center_2DC),[])];          % Number of droplets in each bin
%Nc2 = [reshape(drop_size_dist_CDP2,length(drop_radius_bin_edges_CDP2),[]); reshape(drop_size_dist_2DC2,length(drop_radius_bin_edges_2DC2),[])];          % Number of droplets in each bin

droplet_matrix_center = repmat((drop_radius_bin_center)', 1, length(time))./1e4;            % cm                                                         % cm
%droplet_matrix_center2 = repmat((drop_radius_bin_center2)', 1, length(time));                                                                  % microns
%droplet_matrix_firstEdge = repmat((drop_radius_bin_first_edge)', 1, length(time))./1e4;                                                        % cm
%droplet_matrix_edges = repmat((drop_radius_bin_edges)', 1, length(time))./1e4;                                                                       % microns
%droplet_matrix_edges2 = repmat((drop_radius_bin_edges2)', 1, length(time))./1e4;                                                                       % microns


% -------------------------------------------------------------------
% ----------- To compute efffective radius we need dN/dr ------------
% -------------------------------------------------------------------

%dN_dr = diff(Nc2,1,1)./diff(droplet_matrix_edges2,1,1);                       % #/micron/cm^3     this is the droplet distribution

% Lets compute the effective radius, which is the 3rd moment to the 2nd
% moment

% Compute the ratio of the third moment to the second moment and convert
% back to microns
%re = 1e4 * trapz(drop_bins_cm, droplet_matrix.^3 .* nr, 1)./trapz(drop_bins_cm, droplet_matrix.^2 .* nr, 1);

%re = trapz(drop_bins, droplet_matrix.^3 .* Nc, 1)./trapz(drop_bins, droplet_matrix.^2 .* Nc, 1);

re = double(sum(droplet_matrix_center.^3 .* Nc, 1)./sum(droplet_matrix_center.^2 .* Nc,1) * 1e4);                 % microns

%re2 = sum(droplet_matrix_firstEdge.^3 .* Nc, 1)./sum(droplet_matrix_firstEdge.^2 .* Nc,1) * 1e4;             % microns

%re3 = trapz(drop_radius_bin_center2', droplet_matrix_center2.^3 .* dN_dr,1)./ trapz(drop_radius_bin_center2', droplet_matrix_center2.^2 .* dN_dr,1);              % microns

%re4 = sum(Nc./4 .* diff(droplet_matrix_edges.^4,1,1), 1)./sum(Nc./3 .*diff(droplet_matrix_edges.^3,1,1),1) *1e4;                                               % microns


% Lets compute the total number concetration at each time step by
% integrating over r
%total_Nc = trapz(drop_bins, Nc,1);       % cm^(-3)
total_Nc = sum(Nc,1);                     % cm^(-3)

% ------------------------------------------------------------------
% --------------- Compute liquid water content ---------------------
% ------------------------------------------------------------------

% Lets compute the liquid water content and liquid water path
rho_lw = 1e6;                                                   % g/m^3 - density of liquid water

% we have to convert re to cm in order to have the finals units be in grams
% per meter cubed

lwc = 4/3 * pi *  rho_lw * sum(Nc .* droplet_matrix_center.^3,1);                  % grams of liquid water/meter cubed of air

%lwc2 = 4/3 * pi *  rho_lw * sum(Nc .* droplet_matrix_firstEdge.^3,1);              % grams of liquid water/meter cubed of air





%% Collect all of the data
% We ignore the last bin in the data set. According to the nc info file,
% the data does not extend through all 31 rows.

vocalsRex.Nc = Nc;
vocalsRex.drop_radius_bin_edges = drop_radius_bin_edges;
vocalsRex.drop_radius_bin_center = drop_radius_bin_center;
vocalsRex.total_Nc = total_Nc;
vocalsRex.lwc = lwc;
vocalsRex.time = time;
vocalsRex.startTime = startTime;                                  % We have to assume that this is in UTC time as well
vocalsRex.time_utc = double((startTime(1) + startTime(2)/60) + vocalsRex.time./3600); % UTC time
vocalsRex.latitude = lat;
vocalsRex.longitude = long;
vocalsRex.altitude = altitude;
vocalsRex.re = re;
vocalsRex.horz_wind_speed = horz_wind_speed;
vocalsRex.horz_wind_direction = horz_wind_direction;
vocalsRex.ambient_air_temp = ambient_air_temp;
vocalsRex.water_vapor_presure = water_vapor_pressure;
vocalsRex.SWT = shortwave_top;
vocalsRex.SWB = shortwave_bot;
vocalsRex.LWT = longwave_top;
vocalsRex.LWB = longwave_bot;







end