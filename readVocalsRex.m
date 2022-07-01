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
drop_size_dist = ncread(filename, 'CCDP_RWO')./1e-4;                                                                       % Measured in #/cm^3    
lat = ncread(filename, 'LAT');                                                                                          % Meausred in degrees North
long = ncread(filename, 'LON');                                                                                         % Meausred in degrees East
altitude = ncread(filename, 'ALTX');                                                                                    % Measured in meters above sea level






% The droplet number concentration is filtered into 30 bins from 2 to 52
% microns
% Below 20 microns, the bins are spaced by 1.18 microns
% above 22 microns, the bins are spaced by 2.28 microns
drop_bins = info.Variables(120).Attributes(10).Value;                        % microns

% The data tells us which bins to take
first_bin = info.Variables(120).Attributes(8).Value;
last_bin = info.Variables(120).Attributes(9).Value;

drop_bins = drop_bins(first_bin:last_bin);
drop_size_dist = drop_size_dist(first_bin:last_bin,1,:);

% Lets compute the effective radius, which is the 3rd moment to the 2nd
% moment
nr = reshape(drop_size_dist,length(drop_bins),[]);

% Convert droplet size to cm
drop_bins_cm = drop_bins./1e4;                                          % cm
droplet_matrix = repmat((drop_bins_cm)', 1, length(time));            % cm


% Compute the ratio of the third moment to the second moment and convert
% back to microns
re = 1e4 * trapz(drop_bins_cm, droplet_matrix.^3 .* nr, 1)./trapz(drop_bins_cm, droplet_matrix.^2 .* nr, 1);                      


% Lets compute the total number concetration at each time step by
% integrating over r
total_Nc = trapz(drop_bins_cm, nr,1);       % cm^(-3)

% We ignore the last bin in the data set. According to the nc info file,
% the data does not extend through all 31 rows. 

vocalsRex.CDP_data = nr;
vocalsRex.Nc = total_Nc;
vocalsRex.time = time;
vocalsRex.startTime = startTime;                                  % We have to assume that this is in UTC time as well
vocalsRex.drop_bins = drop_bins;
vocalsRex.latitude = lat;
vocalsRex.longitude = long;
vocalsRex.altitude = altitude;
vocalsRex.re = re;





end