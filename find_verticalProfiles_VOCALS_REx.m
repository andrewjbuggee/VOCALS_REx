%% Finding VOCALS-REx data where the plane flies cleanly through a cloud layer
% Save these vertical profiles


function [vert_profs] = find_verticalProfiles_VOCALS_REx(vocalsRex)


% Let's also store the time vector in UTC format

UTC_starttime = vocalsRex.startTime(1) + vocalsRex.startTime(2)/60;   % hours.decimalhours in UTC format



% dz/dt must be non-zero. 
% Total Nc has to start at a value below 1
% Total Nc has to end at a value below 1

% First lets crop the data only to those portions where the plane is
% ascending or descending

dz_dt = diff(vocalsRex.altitude)./diff(double(vocalsRex.time))';

% Compute the mean with a sliding window for every 10 data points
% This will smooth out the data and make the horizontal flight segments
% easier to find

dz_dt_mean = movmean(dz_dt,20);


% The plane's vertical velocity exceeds 2 m/s on average when it ascends or
% descends. Use this the window the data. Find 

index_ascend_or_descend = find(abs(dz_dt_mean)>2);


% Find consecutive integers, which represent stand alone profiles where the
% plance is climbing or descending
index_consec = find(diff(index_ascend_or_descend)~=1);
% include a 1 to start
index_consec = [0, index_consec];



% -----------------------------------------------------------------------
% For each break in the data, create a profile
for ii = 1:length(index_consec)-1

    vert_profs.Nc{ii} = vocalsRex.total_Nc(index_ascend_or_descend(index_consec(ii)+1:(index_consec(ii+1))))';
    vert_profs.time{ii} = vocalsRex.time(index_ascend_or_descend(index_consec(ii)+1:(index_consec(ii+1))))';
    vert_profs.time_UTC{ii} = UTC_starttime + double(vert_profs.time{ii})./3600;
    vert_profs.lwc{ii} = vocalsRex.lwc(index_ascend_or_descend(index_consec(ii)+1:(index_consec(ii+1))))';
    vert_profs.latitude{ii} = vocalsRex.latitude(index_ascend_or_descend(index_consec(ii)+1:(index_consec(ii+1))))';
    vert_profs.longitude{ii} = vocalsRex.longitude(index_ascend_or_descend(index_consec(ii)+1:(index_consec(ii+1))))';
    vert_profs.altitude{ii} = vocalsRex.altitude(index_ascend_or_descend(index_consec(ii)+1:(index_consec(ii+1))))';
    vert_profs.re{ii} = vocalsRex.re(index_ascend_or_descend(index_consec(ii)+1:(index_consec(ii+1))))';
    vert_profs.SWT{ii} = vocalsRex.SWT(index_ascend_or_descend(index_consec(ii)+1:(index_consec(ii+1))))';
    vert_profs.SWB{ii} = vocalsRex.SWB(index_ascend_or_descend(index_consec(ii)+1:(index_consec(ii+1))))';
    vert_profs.LWT{ii} = vocalsRex.LWT(index_ascend_or_descend(index_consec(ii)+1:(index_consec(ii+1))))';
    vert_profs.LWB{ii} = vocalsRex.LWB(index_ascend_or_descend(index_consec(ii)+1:(index_consec(ii+1))))';

end

% Grab the non-time-stamped data as well
    vert_profs.drop_radius_bin_edges = vocalsRex.drop_radius_bin_edges;
    vert_profs.drop_radius_bin_center = vocalsRex.drop_radius_bin_center;
    vert_profs.startTime = vocalsRex.startTime;                                                                    % We have to assume that this is in UTC time as well



% -----------------------------------------------------------------------
% Now find which of these profiles start and end with Nc<1 AND have some values above 1 

index2delete = [];

for ii = 1:length(vert_profs.Nc)

    %[vert_profs.Nc{ii}(1), vert_profs.Nc{ii}(end), any(vert_profs.Nc{ii}>10^7), all(vert_profs.Nc{ii}<1)]


    if vert_profs.Nc{ii}(1)>1 || vert_profs.Nc{ii}(end)>1 || any(vert_profs.Nc{ii}>10^7) || all(vert_profs.Nc{ii}<1)

        % if this is true, mark the index for deletion
        index2delete = [index2delete, ii];

    end

end




% Let's delete all cells that met the above conditions
vert_profs.Nc(index2delete) = [];
vert_profs.time(index2delete) = [];
vert_profs.time_UTC(index2delete) = [];
vert_profs.lwc(index2delete) = [];
vert_profs.latitude(index2delete) = [];
vert_profs.longitude(index2delete) = [];
vert_profs.altitude(index2delete) = [];
vert_profs.re(index2delete) = [];
vert_profs.SWT(index2delete) = [];
vert_profs.SWB(index2delete) = [];
vert_profs.LWT(index2delete) = [];
vert_profs.LWB(index2delete) = [];




% -----------------------------------------------------------------------
% Usually at the end of each profile there are several zeros. Let's delete
% all 0 values at the end of the vector ONLY when all remaining values are
% 0

for ii = 1:length(vert_profs.Nc)



    % -----------------------------------------------------------------------
    % find the point where all remaining values are zero and truncate
    % -----------------------------------------------------------------------



    
end





% -----------------------------------------------------------------------
% Now sort through the vertical profiles and get rid of Data points where
% the LWC is below some threshold

LWC_threshold = 0;                              % g/m^3




for ii = 1:length(vert_profs.lwc)

    index2delete = [];

    % Find data points where the threshold is less than this and delete
    % those values
    
    index2delete = vert_profs.lwc{ii}<LWC_threshold;

    % delete all time-stamped data that meet the above logical statement

    vert_profs.Nc{ii}(index2delete) = [];
    vert_profs.time{ii}(index2delete) = [];
    vert_profs.time_UTC{ii}(index2delete) = [];
    vert_profs.lwc{ii}(index2delete) = [];
    vert_profs.latitude{ii}(index2delete) = [];
    vert_profs.longitude{ii}(index2delete) = [];
    vert_profs.altitude{ii}(index2delete) = [];
    vert_profs.re{ii}(index2delete) = [];
    vert_profs.SWT{ii}(index2delete) = [];
    vert_profs.SWB{ii}(index2delete) = [];
    vert_profs.LWT{ii}(index2delete) = [];
    vert_profs.LWB{ii}(index2delete) = [];
        

end







end