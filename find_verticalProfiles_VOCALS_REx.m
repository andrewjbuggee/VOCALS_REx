%% Finding VOCALS-REx data where the plane flies cleanly through a cloud layer
% Save these vertical profiles


% INPUTS:
% -------
%       (1) vocals-rex data set

%       (2) LWC_threshold - (g/m^3) this is a threshold that helps clean
%       the data. If set to 0, the vertical profiles will have values on
%       either end of the true vertical droplet profile that represent
%       measurements outside of what we want. the LWC_threshold can help
%       clean the data by limiting the data to have a certain liquid water
%       content threshold. For example. Painemal and Zuidema (2011) defined
%       the cloud top boundary as the level where the LWC was 0.03 g/m^3.
%       Beyond this, the excluded the data. This value will be used to
%       truncate data before and after the vertical profile.

%       (3) stop_at_max_lwc - this is a logical flag, whose values can either
%       be true or false. If true, each vertical profile found will be
%       truncated at the peak value of the liquid water content. Most
%       vertical droplet profiles will have a LWC that grows with altitude.
%       If this is true, the profile usually peters out shortly after the
%       peak LWC value. IF this is false, the data will be kept until the
%       the LWC is below the LWC threshold defined above.



function [vert_profs] = find_verticalProfiles_VOCALS_REx(vocalsRex, LWC_threshold, stop_at_max_lwc)


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
    
    % read in the total number of droplets per unit volume
    vert_profs.Nc{ii} = vocalsRex.total_Nc(index_ascend_or_descend(index_consec(ii)+1:(index_consec(ii+1))))';
    
    % read in the number of droplets for each size bin, and divide by the
    % width of the size bin to estimate the droplet distribution
    vert_profs.nr{ii} = vocalsRex.Nc(:,index_ascend_or_descend(index_consec(ii)+1:(index_consec(ii+1))))./...
        repmat(diff(double(vocalsRex.drop_radius_bin_edges))', 1, length(index_ascend_or_descend(index_consec(ii)+1:(index_consec(ii+1)))));

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
vert_profs.drop_radius_bin_edges = double(vocalsRex.drop_radius_bin_edges);
vert_profs.drop_radius_bin_center = double(vocalsRex.drop_radius_bin_center);
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
vert_profs.nr(index2delete) = [];
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

% for ii = 1:length(vert_profs.Nc)
%
%
%
%     % -----------------------------------------------------------------------
%     % find the point where all remaining values are zero and truncate
%     % -----------------------------------------------------------------------
%
%
%
%
% end





% -----------------------------------------------------------------------
% Now sort through the vertical profiles and get rid of Data points where
% the LWC is below some threshold

%LWC_threshold = 0;                              % g/m^3

% zero vector incase it is needed below
max_lwc = zeros(1, length(vert_profs.lwc));




for ii = 1:length(vert_profs.lwc)

    % Find the first data point and the last data point where the LWC
    % threshold is exceeded. That is, before the first index, every value
    % should be below the threshold. And all data points after the last
    % index should also be below the threshold

    indexes_above_threshold = vert_profs.lwc{ii}>=LWC_threshold;

    % find the first 0, the first data point that exceeds the LWC
    % threshold. The first value in the index below is where the
    % profile will start. The last index below will be the end of the
    % profile
    indexes2keep = find(indexes_above_threshold);


    % find the max LWC and the index associated with this value
    % The max value wil be used to remove profiles that don't have a max
    % LWC above the minimum threshold. The max_index will be used if the
    % profile needs to stop at the max LWC value
    [max_lwc(ii), index_max] = max(vert_profs.lwc{ii});



    if stop_at_max_lwc == true



        % The data will start at the first index above the minimum
        % threshold, found above, and end at the index of maximum LWC,
        % found in this if statement

        indexes2keep = [indexes2keep(1), index_max];


        %         index_max_2End = (index_max+1):numel(vert_profs.lwc{ii});
        %
        %
        %         index2delete = vert_profs.lwc{ii}<LWC_threshold;
        %
        %         % allow indices beyond the max LWC are set to true, and thus are
        %         % deleted
        %         index2delete(index_max_2End) = true;

    end


    % if none of the values within the vertical profile are above the LWC
    % threshold, delete this profile. But for now, just set it to be a
    % vector of 0's
    if sum(indexes2keep)==0

        vert_profs.Nc{ii} = 0;
        vert_profs.nr{ii} = 0;
        vert_profs.time{ii} = 0;
        vert_profs.time_UTC{ii} = 0;
        vert_profs.lwc{ii} = 0;
        vert_profs.latitude{ii} = 0;
        vert_profs.longitude{ii} = 0;
        vert_profs.altitude{ii} = 0;
        vert_profs.re{ii} = 0;
        vert_profs.SWT{ii} = 0;
        vert_profs.SWB{ii} = 0;
        vert_profs.LWT{ii} = 0;
        vert_profs.LWB{ii} = 0;

    else

        % delete all time-stamped data that meet the above logical statement

        vert_profs.Nc{ii} = vert_profs.Nc{ii}(indexes2keep(1):indexes2keep(end));
        vert_profs.nr{ii} = vert_profs.nr{ii}(:, indexes2keep(1):indexes2keep(end));
        vert_profs.time{ii} = vert_profs.time{ii}(indexes2keep(1):indexes2keep(end));
        vert_profs.time_UTC{ii} = vert_profs.time_UTC{ii}(indexes2keep(1):indexes2keep(end));
        vert_profs.lwc{ii} = vert_profs.lwc{ii}(indexes2keep(1):indexes2keep(end));
        vert_profs.latitude{ii} = vert_profs.latitude{ii}(indexes2keep(1):indexes2keep(end));
        vert_profs.longitude{ii} = vert_profs.longitude{ii}(indexes2keep(1):indexes2keep(end));
        vert_profs.altitude{ii} = vert_profs.altitude{ii}(indexes2keep(1):indexes2keep(end));
        vert_profs.re{ii} = vert_profs.re{ii}(indexes2keep(1):indexes2keep(end));
        vert_profs.SWT{ii} = vert_profs.SWT{ii}(indexes2keep(1):indexes2keep(end));
        vert_profs.SWB{ii} = vert_profs.SWB{ii}(indexes2keep(1):indexes2keep(end));
        vert_profs.LWT{ii} = vert_profs.LWT{ii}(indexes2keep(1):indexes2keep(end));
        vert_profs.LWB{ii} = vert_profs.LWB{ii}(indexes2keep(1):indexes2keep(end));


    end


end




% If the max LWC content found is below the LWC threshold, the vertical
% profile will be an empty vector, and we should delete it

find_max_less_than_threshold = max_lwc<LWC_threshold;

% if this is true, delete all entries for the index
if sum(find_max_less_than_threshold)>0

    % delete all time-stamped data that meet the above logical statement

    vert_profs.Nc(find_max_less_than_threshold) = [];
    vert_profs.nr(find_max_less_than_threshold) = [];
    vert_profs.time(find_max_less_than_threshold) = [];
    vert_profs.time_UTC(find_max_less_than_threshold) = [];
    vert_profs.lwc(find_max_less_than_threshold) = [];
    vert_profs.latitude(find_max_less_than_threshold) = [];
    vert_profs.longitude(find_max_less_than_threshold) = [];
    vert_profs.altitude(find_max_less_than_threshold) = [];
    vert_profs.re(find_max_less_than_threshold) = [];
    vert_profs.SWT(find_max_less_than_threshold) = [];
    vert_profs.SWB(find_max_less_than_threshold) = [];
    vert_profs.LWT(find_max_less_than_threshold) = [];
    vert_profs.LWB(find_max_less_than_threshold) = [];

end



% Store the LWC threshold used
vert_profs.lwc_threshold = LWC_threshold;            % g/m^3






end