

function vocalsRex = cropVocalsRex2MODIS(vocalsRex, modis,vocalsRexFolder,modisFolder)



% convert the start time to decimal time
sec_per_hour = 3600;                                                                    % sec
sec_per_min = 60;                                                                       % sec
sec_per_day = 86400;                                                                    % sec

modis_dataTime = (modis.time(1)*sec_per_hour + modis.time(2)*sec_per_min)/sec_per_day;

startTime_dec = (vocalsRex.startTime(1)*sec_per_hour + vocalsRex.startTime(2)*sec_per_min)/sec_per_day;                      % fraction of a day 

% define the time within the vocals-rex profile where a nice cloud profile
% exists

if strcmp(modisFolder(96:end),'/2008_11_11_1430/')==true && strcmp(vocalsRexFolder(74:end),'/2008_11_11/')==true

    cloud_profile_time = 0.6069;

    index_delay = 4;

elseif strcmp(modisFolder(96:end),'/2008_11_09/')==true && strcmp(vocalsRexFolder(74:end),'/2008_11_09/')==true

    cloud_profile_time = 0.6120;

    index_delay = 0;

elseif strcmp(modisFolder(96:end),'/2008_11_11_1850/')==true && strcmp(vocalsRexFolder(74:end),'/2008_11_11/')==true

    cloud_profile_time = 0.7816;

    index_delay = -8;

else

    error([newline, 'What is the decimal time of the cloud profiel you wish to look at?', newline])

end


cloudProfile_secSinceStart = floor((cloud_profile_time - startTime_dec)*sec_per_day);                               % seconds since flight started up to the time indicated by painemal and zudema

modis_secSinceStart = floor((modis_dataTime - startTime_dec)*sec_per_day);


% Define the number of discrete points in time you wish to look at. The
% center point will be at the start time calculated above
windowLength = 100;

% Using the time defined above, find the location in space closest to
% the airplanes location

% Using lat and long with can minimize the euclidean norm
% ***** There are some special cases *****

% The values are closer on Novemeber 11th 2008 if we compare one time step
% ahead. The re values are spatially similar, but the optical thickness
% change quite a bit!
dist_btwn_PZ_startTime_and_MODIS = sqrt((modis.geo.lat - vocalsRex.latitude(cloudProfile_secSinceStart+index_delay)).^2 + (modis.geo.long - vocalsRex.longitude(cloudProfile_secSinceStart+index_delay)).^2); 
[~, index_minDist] = min(dist_btwn_PZ_startTime_and_MODIS, [], 'all');



% Cloud bottom = location where LWC = 0.03 g/m^3
% ****** There are two ways we could define LWC *****
% Cloud top is tricky. The vocals rex data shows droplets being recorded
% for several data points past the peak LWC value. But MODIS estimates for
% optical depth better align with the cloud top being defined as the peak
% LWC value. After this maximum, LWC tends to fall off dramatically

% Cloud top = location where re goes to 0 after LWC> 0.03 g/m^3
% Cloud top = maximum valude of LWC after LWC> 0.03 g/m^3

% First find the index cloud bottom using the definition above

lwc_lim = 0.03;                                                 % grams/m^3 - lower limit defining cloud base

LWC_window_index = cloudProfile_secSinceStart - windowLength/2 : cloudProfile_secSinceStart + windowLength/2;
LWC_window = vocalsRex.lwc(LWC_window_index);
index_minVal = find(LWC_window >= lwc_lim);
index_cloudBottom = LWC_window_index(1) + index_minVal(1) -1;



% Now lets find the index for cloud top using the definition above
% find the first value of re=0 after the above index

window_data = vocalsRex.re(index_cloudBottom : index_cloudBottom + 100);
index_nan = find(isnan(window_data));

index_cloudTop = index_cloudBottom + index_nan(1) - 2;

% Lets also find the index where LWC is a maxium after the index found for
% cloud bottom

window_data = vocalsRex.lwc(index_cloudBottom : index_cloudBottom + 100);
[~,index_maxLWC] = max(window_data);

index_cloudTop2 = index_cloudBottom + index_maxLWC - 1;



vocalsRex.total_Nc = vocalsRex.total_Nc(index_cloudBottom:index_cloudTop2);
vocalsRex.time = vocalsRex.time(index_cloudBottom:index_cloudTop2);
vocalsRex.latitude = vocalsRex.latitude(index_cloudBottom:index_cloudTop2);
vocalsRex.longitude = vocalsRex.longitude(index_cloudBottom:index_cloudTop2);
vocalsRex.altitude = vocalsRex.altitude(index_cloudBottom:index_cloudTop2);
vocalsRex.re = vocalsRex.re(index_cloudBottom:index_cloudTop2);
vocalsRex.lwc = vocalsRex.lwc(index_cloudBottom:index_cloudTop2);
vocalsRex.startTime = vocalsRex.startTime;                               % We have to assume that this is in UTC time as well
vocalsRex.modisIndex_minDist = index_minDist;




end