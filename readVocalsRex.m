%% Read in Vocals Rex Data


% By Andrew John Buggee
%%

function vocalsRex = readVocalsRex(filename)

info = ncinfo(filename);

% read the time vector
time = ncread(filename, 'Time');                                                                                        % Measured in seconds

% [hours, minutes]
startTime = [str2double(info.Variables(1).Attributes(3).Value(26:27)),...
    str2double(info.Variables(1).Attributes(3).Value(29:30))];      % Start time of the data log in UTC

% ----------------------------------------------------------------
% ----------- cloud droplet probe measurements (CDP) -------------
% ----------------------------------------------------------------

drop_size_dist_CDP = ncread(filename, 'CCDP_RWO');                                                                       % Measured in #/micron/cm^3    

% The droplet number concentration is filtered into 30 bins from 2 to 52
% microns
% Below 20 microns, the bins are spaced by 1.18 microns
% above 22 microns, the bins are spaced by 2.28 microns
drop_bins_CDP = info.Variables(120).Attributes(10).Value;                        % microns

% The data tells us which bins to take
first_bin_CDP = info.Variables(120).Attributes(8).Value;
last_bin_CDP = info.Variables(120).Attributes(9).Value;

drop_bins_CDP = drop_bins_CDP(first_bin_CDP:last_bin_CDP);                                  % microns
drop_size_dist_CDP = drop_size_dist_CDP(first_bin_CDP:last_bin_CDP,1,:);



% ----------------------------------------------------------------
% ---------- Two-Dimensional Optical array probe (2DC) -----------
% ----------------------------------------------------------------

% This droplet probe filters droplets into 64 different size bins. The bins
% range from 25 microns to 1560 microns. 

% ** IMPORTANT ** The first bin of the 2DC probe overlaps with the last bin
% of the CDP. So we ignore it. (See Painemal and Zuidema paragraph 15.)

drop_size_dist_2DC = ncread(filename, 'C1DCA_RPC');                                             % Measured in #/micron/L  - where I think L stands for liters. 1 liter = 1000 cubic cm

% convert the 2DC data to inverse cubic centimeters
drop_size_dist_2DC = drop_size_dist_2DC ./ 1000;                                                % # of droplets/microns/cm^3

% Grab the droplet size bins
% These are the boundaries that define each bin size
drop_bins_2DC = info.Variables(119).Attributes(12).Value;                        % microns

% The data tells us which bins to take
first_bin_2DC = double(info.Variables(119).Attributes(10).Value);
last_bin_2DC = double(info.Variables(119).Attributes(11).Value);

drop_bins_2DC = drop_bins_2DC(first_bin_2DC:last_bin_2DC);                                  % microns
drop_size_dist_2DC = drop_size_dist_2DC(first_bin_2DC:last_bin_2DC,1,:);


% Grab position and timing data
lat = ncread(filename, 'LAT');                                                                                          % Meausred in degrees North
long = ncread(filename, 'LON');                                                                                         % Meausred in degrees East
altitude = ncread(filename, 'ALTX');                                                                                    % Measured in meters above sea level


% lets grab the shortwave and longwave radiance looking down and looking up
shortwave_top = ncread(filename, 'SWT');
shortwave_bot = ncread(filename, 'SWB');

longwave_top = ncread(filename,'IRTC');              % These are corrected longwave irradiance values
longwave_bot = ncread(filename, 'IRBC');



%     ***** IMPORTANT STEP *****

% Combine data from both droplet probes to get the total droplet size
% distribution


% Lets compute the effective radius, which is the 3rd moment to the 2nd
% moment
drop_bins = [drop_bins_CDP, drop_bins_2DC];

nr = [reshape(drop_size_dist_CDP,length(drop_bins_CDP),[]); reshape(drop_size_dist_2DC,length(drop_bins_2DC),[])];

droplet_matrix = repmat((drop_bins)', 1, length(time));            % microns


% Compute the ratio of the third moment to the second moment and convert
% back to microns
%re = 1e4 * trapz(drop_bins_cm, droplet_matrix.^3 .* nr, 1)./trapz(drop_bins_cm, droplet_matrix.^2 .* nr, 1);                      

re = trapz(drop_bins, droplet_matrix.^3 .* nr, 1)./trapz(drop_bins, droplet_matrix.^2 .* nr, 1);                      


% Lets compute the total number concetration at each time step by
% integrating over r
total_Nc = trapz(drop_bins, nr,1);       % cm^(-3)


% Lets compute the liquid water content and liquid water path
rho_lw = 1e6;                                                   % g/m^3 - density of liquid water

% we have to convert re to cm in order to have the finals units be in grams
% per meter cubed

lwc = 4/3 * pi *  rho_lw * sum(nr .* (droplet_matrix./1e4).^3,1);                    % grams of liquid water/meter cubed of air


%% Collect all of the data
% We ignore the last bin in the data set. According to the nc info file,
% the data does not extend through all 31 rows. 

vocalsRex.nr = nr;
vocalsRex.drop_bins = drop_bins;
vocalsRex.Nc = total_Nc;
vocalsRex.lwc = lwc;
vocalsRex.time = time;
vocalsRex.startTime = startTime;                                  % We have to assume that this is in UTC time as well
vocalsRex.latitude = lat;
vocalsRex.longitude = long;
vocalsRex.altitude = altitude;
vocalsRex.re = re;
vocalsRex.SWT = shortwave_top;
vocalsRex.SWB = shortwave_bot;
vocalsRex.LWT = longwave_top;
vocalsRex.LWB = longwave_bot;







end