%% Calculate ensemble statistics of cloud top and bottom radii measured by VOCALS-REx


% Andrew John Buggee

clear variables

%% Loop through all VOCALS-REx files and load all vertical profiles

% Macbook folder name
foldername = '/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/VOCALS_REx/';


% Oct-18-2008 Data
filename{1} = 'RF02.20081018.130300_213000.PNI.nc';

% Oct-25-2008 data
filename{2} = 'RF05.20081025.062900_152500.PNI.nc';

% Nov-02-2008 data
filename{3} = 'RF08.20081102.055700_152100.PNI.nc';

% Nov-09-2008 data 
filename{4} = 'RF11.20081109.125700_213600.PNI.nc';

% Nov-11-2008 data 
filename{5} = 'RF12.20081111.125000_214500.PNI.nc';

% Nov-13-2008 data 
filename{6} = 'RF13.20081113.125700_215700.PNI.nc';



% define the LWC threshold
LWC_threshold = 0.02;       % g/m^3

% do you want to cut off the profile at the maximum LWC value?
stop_at_max_LWC = false;



% Load data

for nn = 1:length(filename)

    vocalsRex = readVocalsRex([foldername,filename{nn}(6:9),'_',filename{nn}(10:11),...
                '_',filename{nn}(12:13),'/',filename{nn}]);

    % find the vertical profiles
    vert_profs = find_verticalProfiles_VOCALS_REx(vocalsRex, LWC_threshold, stop_at_max_LWC);

    % grab just the variables of interest for from each vertical profile
    if nn==1
        
        ensemble_profiles.altitude = vert_profs.altitude;
        ensemble_profiles.tau = vert_profs.tau;
        ensemble_profiles.re = vert_profs.re;
        ensemble_profiles.lwc = vert_profs.lwc;
        ensemble_profiles.Nc = vert_profs.Nc;
        ensemble_profiles.time = vert_profs.time;

    else
        
        ensemble_profiles.altitude = [ensemble_profiles.altitude, vert_profs.altitude];
        ensemble_profiles.tau = [ensemble_profiles.tau, vert_profs.tau];
        ensemble_profiles.re = [ensemble_profiles.re, vert_profs.re]; 
        ensemble_profiles.lwc = [ensemble_profiles.lwc, vert_profs.lwc];
        ensemble_profiles.Nc = [ensemble_profiles.Nc, vert_profs.Nc];
        ensemble_profiles.time = [ensemble_profiles.time, vert_profs.time];

    end


end

% store the LWC threshold

ensemble_profiles.lwc_threshold = LWC_threshold;



%% Plot a histogram of the optical depths and the geometric depths

cloud_optical_depth = zeros(1, length(ensemble_profiles.tau));
cloud_geometric_thickness = zeros(1, length(ensemble_profiles.altitude));


for nn = 1:length(ensemble_profiles.tau)

    cloud_optical_depth(nn) = ensemble_profiles.tau{nn}(end);
    cloud_geometric_thickness(nn) = abs(ensemble_profiles.altitude{nn}(1) - ensemble_profiles.altitude{nn}(end));

end


figure; 

subplot(1,2,1)
histogram(cloud_geometric_thickness, 'NumBins',20)
xlabel('Cloud Geometric Thickness (m)')
ylabel('Counts')
title([num2str(length(cloud_optical_depth)), ' vertical profiles'])

subplot(1,2,2)
histogram(cloud_optical_depth, 'NumBins',20)
xlabel('Cloud Optical Thickness')
ylabel('Counts')

set(gcf, 'Position', [0 0 1000 550])


%%  Segment re, LWC, Nc into N bins along optical depth

% In order to compute a median vertical profile, we have to first normalize
% the vertical extent so that all profiles lie between values [0,1]. Then
% we break up the vertical component in n discrete bins. Within each bin we
% can compute the mean, median and standard deviation

n_bins = 20; % number of segments the noramlized vertical component is broken up into

bin_edges = 0:1/n_bins:1;

% set up an empty cell array for all the values of each variable of interest
% within each segment boundaries. Let's do this for droplet size, total
% number concentration and liquid water content
vertically_segmented_attributes = cell(n_bins, 3);


normalized_tau = cell(1, length(ensemble_profiles.altitude));


for nn = 1:length(ensemble_profiles.altitude)

    % first we need to normalize the vertical component of all profiles
    normalized_tau{nn} = ensemble_profiles.tau{nn}./ensemble_profiles.tau{nn}(end);
    
    % the data is stored in altitude space. If we wish to have the data
    % oriented in optical depth space, we need to reverse the order,
    % because optical depth is meausred from cloud top to cloud bottom
    % so we need to check if the plane is ascending or descending so we
    % know whether the data starts at cloud top or cloud bottom
    % We want all data oriented in tau space. So look for ascending
    % profiles, and flip the vector. If the plane is descending, we don't
    % need to do anything
    if mean(diff(ensemble_profiles.altitude{nn})./diff(ensemble_profiles.time{nn}))>0
        % if true, then the plane is ascending, grab and flip the variables
        % of interest
        re = flipud(ensemble_profiles.re{nn});
        lwc = flipud(ensemble_profiles.lwc{nn});
        Nc = flipud(ensemble_profiles.Nc{nn});

    else
        % if false, then the plane is descending, just grab the variables
        % of interest
        re = ensemble_profiles.re{nn};
        lwc = ensemble_profiles.lwc{nn};
        Nc = ensemble_profiles.Nc{nn};

    end


    


    % for each profile, we need to segment the variables of interest into n
    % bins. 

    for bb = 1:length(bin_edges)-1

        % grab all re, LWC, and Nc values within each bin. Segment them
        % accordingly
        if bb==1
            index_segment = normalized_tau{nn}>=bin_edges(bb) & normalized_tau{nn}<=bin_edges(bb+1);

        else
            index_segment = normalized_tau{nn}>bin_edges(bb) & normalized_tau{nn}<=bin_edges(bb+1);
        end

        % store the effective radius values
        vertically_segmented_attributes{bb, 1} = [vertically_segmented_attributes{bb, 1}; re(index_segment)];

        % store the liquid water content values
        vertically_segmented_attributes{bb, 2} = [vertically_segmented_attributes{bb, 2}; lwc(index_segment)];
    
        % store the total droplet number concentration values
        vertically_segmented_attributes{bb, 3} = [vertically_segmented_attributes{bb, 3}; Nc(index_segment)];



    end



end



%% Plot the mean and standard deviation of each vertical bin for the effective radius


re_mean = zeros(n_bins, 1);
re_std = zeros(n_bins, 1);
bin_center = zeros(n_bins, 1);



for bb = 1:n_bins
    
    % compute the mean value for the current bin
    re_mean(bb) = mean(vertically_segmented_attributes{bb,1});       % microns - mean droplet radius value

    % compute the standard deviation of the current bin
    re_std(bb) = std(vertically_segmented_attributes{bb,1});         % microns - standard deviation
    
    % compute the bin center which is the tau location of the mean data
    bin_center(bb) = (bin_edges(bb+1) - bin_edges(bb))/2 + bin_edges(bb);

end

figure;
errorbar(re_mean, bin_center, re_std, re_std, 'horizontal', 'LineStyle','-',...
    'Marker','o', 'MarkerFaceColor', mySavedColors(1, 'fixed'), 'MarkerEdgeColor','k', ...
    'MarkerSize', 10, 'LineWidth', 2, 'Color', 'k')

set(gca, 'YDir', 'reverse')

grid on; grid minor
xlabel('Mean r_e')
ylabel('Normalized Optical Depth')

set(gcf, 'Position', [0 0 650 600])


%% Plot the median and standard deviation of each vertical bin for the effective radius


re_median = zeros(n_bins, 1);
re_std = zeros(n_bins, 1);
bin_center = zeros(n_bins, 1);



for bb = 1:n_bins
    
    % compute the mean value for the current bin
    re_median(bb) = median(vertically_segmented_attributes{bb,1});       % microns - mean droplet radius value

    % compute the standard deviation of the current bin
    re_std(bb) = std(vertically_segmented_attributes{bb,1});         % microns - standard deviation
    
    % compute the bin center which is the tau location of the mean data
    bin_center(bb) = (bin_edges(bb+1) - bin_edges(bb))/2 + bin_edges(bb);

end

figure;
errorbar(re_median, bin_center, re_std, re_std, 'horizontal', 'LineStyle','-',...
    'Marker','o', 'MarkerFaceColor', mySavedColors(1, 'fixed'), 'MarkerEdgeColor','k', ...
    'MarkerSize', 10, 'LineWidth', 2, 'Color', 'k')

set(gca, 'YDir', 'reverse')

grid on; grid minor
xlabel('Median r_e')
ylabel('Normalized Optical Depth')

set(gcf, 'Position', [0 0 650 600])


%% Plot the mean and median effective radius profiles

figure; 

plot(re_mean, bin_center, 'Color', mySavedColors(1, 'fixed'));
hold on
plot(re_median, bin_center, 'Color', mySavedColors(2, 'fixed'))

legend('Mean', 'Median', 'Location', 'best')

set(gca, 'YDir', 'reverse')

grid on; grid minor
xlabel('<r_e>')
ylabel('Normalized Optical Depth')

set(gcf, 'Position', [0 0 650 600]) 

