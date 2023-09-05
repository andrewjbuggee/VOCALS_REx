

function vocalsRex = cropVocalsRex2MODIS(vocalsRex, lwc_threshold, stop_at_max_lwc, Nc_threshold, modis)


% ----- Find all vertical profiles within VOCALS-REx data ------
vert_prof = find_verticalProfiles_VOCALS_REx(vocalsRex, lwc_threshold, stop_at_max_lwc, Nc_threshold);


min_per_hour = 60;                                                                      % min

modis_timeUTC = modis.time(1) + modis.time(2)/min_per_hour;



% Now combine all time vectors together so we can find which profile is
% closest in time to the MODIS data
time2modis = zeros(1, length(vert_prof.time_utc));
for nn = 1:length(vert_prof.time_utc)
    time2modis(nn) = min(abs(vert_prof.time_utc{nn} - modis_timeUTC)); % times are in hours
end

% Determine which profile is closest to the time MODIS data was recorded
[~, index_vertProf_closest2modis] = min(time2modis);


% Let's keep the vertical profile that is closest to MODIS
clear vocalsRex;

% to get the data we want, we need to convert the structure to a cell array
fields = fieldnames(vert_prof);
vert_prof_cell = struct2cell(vert_prof);

data2keep = cell(1, numel(vert_prof_cell));


% step through each field. If its a cell, only keep the index found above
for ii = 1:length(vert_prof_cell)
    if iscell(vert_prof_cell{ii})==true

        data2keep{ii} = vert_prof_cell{ii}{index_vertProf_closest2modis};

    else

        data2keep{ii} = vert_prof_cell{ii};

    end
end

% Convert back to a structure and keep the data closest in time to MODIS
vocalsRex = cell2struct(data2keep, fields, 2);


% Find the MODIS pixel closest to the vocalsRex data by using the median of
% the VOCALS-REx latitude and longitude
dist_btwn_VR_startTime_and_MODIS = sqrt((double(modis.geo.lat) - median(vocalsRex.latitude)).^2 + (double(modis.geo.long) - median(vocalsRex.longitude)).^2); 
[~, index_minDist] = min(dist_btwn_VR_startTime_and_MODIS, [], 'all');
% save this index so we know what MODIs pixel to use in comparison
vocalsRex.modisIndex_minDist_median = index_minDist;


% Find the MODIS pixel closest to the vocalsRex data by using the first
% latitude and longitude point of the vertical profile
dist_btwn_VR_startTime_and_MODIS = sqrt((double(modis.geo.lat) - (vocalsRex.latitude(1))).^2 + (double(modis.geo.long) - (vocalsRex.longitude(1))).^2); 
[~, index_minDist] = min(dist_btwn_VR_startTime_and_MODIS, [], 'all');
% save this index so we know what MODIs pixel to use in comparison
vocalsRex.modisIndex_minDist_first = index_minDist;


% Find the MODIS pixel closest to the vocalsRex data by using the last
% latitude and longitude point of the vertical profile
dist_btwn_VR_startTime_and_MODIS = sqrt((double(modis.geo.lat) - (vocalsRex.latitude(end))).^2 + (double(modis.geo.long) - (vocalsRex.longitude(end))).^2); 
[~, index_minDist] = min(dist_btwn_VR_startTime_and_MODIS, [], 'all');
% save this index so we know what MODIs pixel to use in comparison
vocalsRex.modisIndex_minDist_last = index_minDist;



% ---- Check to see if all 3 indices are unique. If not, delete the
% redundancy

[~, idx_unique] = unique([vocalsRex.modisIndex_minDist_first, vocalsRex.modisIndex_minDist_median, vocalsRex.modisIndex_minDist_last]);
idx_unique_logic = ismember([1,2,3], idx_unique);

% delete the indices that are not found above


if idx_unique_logic(1)==false
    
    % set an empty array
    vocalsRex.modisIndex_minDist_first = [];

elseif idx_unique_logic(2)==false
    
    % set an empty array
    vocalsRex.modisIndex_minDist_median = [];

elseif idx_unique_logic(3)==false
    
    % set an empty array
    vocalsRex.modisIndex_minDist_last = [];

end






% ---------------------------------------------------------------------
% -------------------------- OLD CODE --------------------------------
% ---------------------------------------------------------------------







% ---------------------------------------------------------------------
% ---- Below is the times associated with a MANUAL profile search -----
% ---------------------------------------------------------------------

% define the time within the vocals-rex profile where a nice cloud profile
% exists

% % convert the MODIS time to decimal time (fraction of 1 day)
% sec_per_hour = 3600;                                                                    % sec
% sec_per_min = 60;                                                                       % sec
% sec_per_day = 86400;                                                                    % sec

% % Convert start time of VOCALS-REx flight to same decimal time (fraction of
% % 1 day)
% startTime_dec = (vocalsRex.startTime(1)*sec_per_hour + vocalsRex.startTime(2)*sec_per_min)/sec_per_day;      % fraction of a day 
% 
% % Convert the rest of the Vocals-REx data record time into fractions of a
% % single day
% vocalsRex_time = startTime_dec + double(vocalsRex.time)./sec_per_day;



% if strcmp(modisFolder(end-16:end),'/2008_11_11_1430/')==true && strcmp(vocalsRexFolder(end-11:end),'/2008_11_11/')==true
% 
%     cloud_profile_time = 0.6069;
% 
%     index_delay = 4;
% 
% elseif strcmp(modisFolder(end-11:end),'/2008_11_09/')==true && strcmp(vocalsRexFolder(end-11:end),'/2008_11_09/')==true
% 
%     cloud_profile_time = 0.6120;
% 
%     index_delay = 0;
% 
% elseif strcmp(modisFolder(end-16:end),'/2008_11_11_1850/')==true && strcmp(vocalsRexFolder(end-11:end),'/2008_11_11/')==true
% 
%     cloud_profile_time = 0.7816;
% 
%     index_delay = -8;
% 
% else
% 
%     error([newline, 'What is the decimal time of the cloud profile you wish to look at?', newline])
% 
% end
% 
% 
% cloudProfile_secSinceStart = floor((cloud_profile_time - startTime_dec)*sec_per_day);                               % seconds since flight started up to the time indicated by painemal and zudema
% 
% 
% % Define the number of discrete points in time you wish to look at. The
% % center point will be at the start time calculated above
% windowLength = 100;
% 
% % Using the time defined above, find the location in space closest to
% % the airplanes location
% 
% % Using lat and long with can minimize the euclidean norm
% % ***** There are some special cases *****
% 
% % The values are closer on Novemeber 11th 2008 if we compare one time step
% % ahead. The re values are spatially similar, but the optical thickness
% % change quite a bit!
% dist_btwn_PZ_startTime_and_MODIS = sqrt((modis.geo.lat - vocalsRex.latitude(cloudProfile_secSinceStart+index_delay)).^2 + (modis.geo.long - vocalsRex.longitude(cloudProfile_secSinceStart+index_delay)).^2); 
% [~, index_minDist] = min(dist_btwn_PZ_startTime_and_MODIS, [], 'all');
% 
% 
% 
% 
% % ****** There are two ways we could define Cloud Top *****
% % Cloud top is tricky. The vocals rex data shows droplets being recorded
% % for several data points past the peak LWC value. But MODIS estimates for
% % optical depth better align with the cloud top being defined as the peak
% % LWC value. After this maximum, LWC tends to fall off dramatically
% 
% % Cloud top = location where re goes to 0 after LWC> 0.03 g/m^3
% % Cloud top = maximum value of LWC after LWC> 0.03 g/m^3
% 
% % First find the index cloud bottom using the definition above
% 
% % ---------------------------------------------------------------------
% % Cloud bottom = location where LWC = 0.03 g/m^3 - Painemal and Zuidema
% % (2011) defintion
% % ---------------------------------------------------------------------
% 
% lwc_lim = 0.03;                                                 % grams/m^3 - lower limit defining cloud base
% 
% LWC_window_index = cloudProfile_secSinceStart - windowLength/2 : cloudProfile_secSinceStart + windowLength/2;
% LWC_window = vocalsRex.lwc(LWC_window_index);
% index_minVal = find(LWC_window >= lwc_lim);
% index_cloudBottom = LWC_window_index(1) + index_minVal(1) -1;
% 
% 
% 
% % Now lets find the index for cloud top by finding the value of effective
% % radius goes to 0
% 
% window_data = vocalsRex.re(index_cloudBottom : index_cloudBottom + 100);
% index_nan = find(isnan(window_data));
% 
% index_cloudTop = index_cloudBottom + index_nan(1) - 2;
% 
% % Lets also find the index where LWC is a maxium after the index found for
% % cloud bottom
% 
% window_data = vocalsRex.lwc(index_cloudBottom : index_cloudBottom + 100);
% [~,index_maxLWC] = max(window_data);
% 
% index_cloudTop2 = index_cloudBottom + index_maxLWC - 1;
% 
% 
% 
% vocalsRex.total_Nc = vocalsRex.total_Nc(index_cloudBottom:index_cloudTop2);
% vocalsRex.time = vocalsRex.time(index_cloudBottom:index_cloudTop2);
% vocalsRex.latitude = vocalsRex.latitude(index_cloudBottom:index_cloudTop2);
% vocalsRex.longitude = vocalsRex.longitude(index_cloudBottom:index_cloudTop2);
% vocalsRex.altitude = vocalsRex.altitude(index_cloudBottom:index_cloudTop2);
% vocalsRex.re = vocalsRex.re(index_cloudBottom:index_cloudTop2);
% vocalsRex.lwc = vocalsRex.lwc(index_cloudBottom:index_cloudTop2);
% vocalsRex.startTime = vocalsRex.startTime;                               % We have to assume that this is in UTC time as well
% vocalsRex.modisIndex_minDist = index_minDist;




end