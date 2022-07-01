

function vocalsRex = cropVocalsRex2MODIS(vocalsRex, modisTime)

% Find the hours and minutes bewteen the vocalsRex start time and the MODIS
% start time

% Lets convert both times to seconds since midnight
seconds_modis = 3600*modisTime(1) + 60*modisTime(2);

seconds_Vocals = 3600*vocalsRex.startTime(1) + 60*vocalsRex.startTime(2);

% Vocals Rex has to have started flying before modis, or this calcualtion
% wouldn't make sense
if seconds_Vocals>seconds_modis
    error([newline, 'The VOCALS-REx data set started recording after MODIS recorded its data!', newline])
end

% Now we can assume MODIS recording time is greater (in seconds since
% midnight) than the Vocals Rex start time


diff_total_seconds = seconds_modis - seconds_Vocals;                     % sec

% Now that we have the number of seconds between the two data sets, lets
% find the vocalsRex data that is within 2.5 minutes before and after the
% MODIS recorded time
minutes_buffer = 5;               % Number of minutes to bracker around the MODIS recording time
time_window = minutes_buffer*60;               % seconds



index_timeWindow = vocalsRex.time>(diff_total_seconds-time_window) & vocalsRex.time<(diff_total_seconds+time_window);

% Now crop the data and only keep the values that fall within this
% timeWindow

% Record the seconds after VocalsRex started when the two data sets overlap
vocalsRex.overlapTime = diff_total_seconds;

vocalsRex.CDP_data = vocalsRex.CDP_data(index_timeWindow);
vocalsRex.Nc = vocalsRex.Nc(index_timeWindow);
vocalsRex.time = vocalsRex.time(index_timeWindow);
vocalsRex.latitude = vocalsRex.latitude(index_timeWindow);
vocalsRex.longitude = vocalsRex.longitude(index_timeWindow);
vocalsRex.altitude = vocalsRex.altitude(index_timeWindow);
vocalsRex.re = vocalsRex.re(index_timeWindow);





end