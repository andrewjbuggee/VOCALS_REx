%% Calculate ensemble statistics of cloud top and bottom radii measured by VOCALS-REx


% Andrew John Buggee

clear variables

%% Loop through all VOCALS-REx files and load all vertical profiles

% Macbook folder name
foldername = '/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data/';

% read all file names in the above listed folder
folder_contents = dir(foldername);

% grab just the vocals filenames

filename = cell(1, length(folder_contents));
backspace = 0;

for nn = 1:length(folder_contents)

    if length(folder_contents(nn).name)>1

        if strcmp(folder_contents(nn).name(1:2), 'RF')==true

            filename{nn - backspace} = folder_contents(nn).name;

        else

            % delete this cell array
            filename(nn - backspace) = [];
            % add one to backspace
            backspace = backspace+1;

        end

    else

        % delete this cell array
        filename(nn - backspace) = [];

        % add one to backspace
        backspace = backspace+1;

    end

end






% define the LWC threshold
LWC_threshold = 0.03;       % g/m^3

% do you want to cut off the profile at the maximum LWC value?
stop_at_max_LWC = false;

% define the total number concentration threshold
Nc_threshold = 1;       %  #-droplets/cm^3




% Load data

for nn = 1:length(filename)

    disp(['nn = ',num2str(nn)])

    vocalsRex = readVocalsRex([foldername,filename{nn}]);

    % find the vertical profiles
    vert_profs = find_verticalProfiles_VOCALS_REx(vocalsRex, LWC_threshold, stop_at_max_LWC, Nc_threshold);

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

% save the ensemble profiles
save([foldername,'ensemble_profiles_from_',num2str(length(filename)), '_files.mat'],...
    'ensemble_profiles', 'filename')




%% Ensemble profiles contains ALL vertical profiles whether they increase or decrease
% Let's split them up into increasing droplet size or decreasing droplet
% size

% Let's step through each profile found
for nn = 1:length(ensemble_profiles.re)

    % compute a moving average every 5 data points which will smooth out
    % the small scale variability
    mean_re_profile = movmean(ensemble_profiles.re{nn}, 5);

    % Now compute the average derivative of this smooth droplet profile with
    % respect to altitude
    dre_dz = mean(diff(mean_re_profile)./diff(ensemble_profiles.altitude{nn}));

    % sort out profiles between those where droplet size increase with
    % altitude, and those that decrease
    if dre_dz>0

        % grab just the variables of interest for from each vertical profile
        if exist("increasing_profiles", "var")==false

            increasing_profiles.altitude = ensemble_profiles.altitude(nn);
            increasing_profiles.tau = ensemble_profiles.tau(nn);
            increasing_profiles.re = ensemble_profiles.re(nn);
            increasing_profiles.lwc = ensemble_profiles.lwc(nn);
            increasing_profiles.Nc = ensemble_profiles.Nc(nn);
            increasing_profiles.time = ensemble_profiles.time(nn);

        else

            increasing_profiles.altitude = [increasing_profiles.altitude, ensemble_profiles.altitude(nn)];
            increasing_profiles.tau = [increasing_profiles.tau, ensemble_profiles.tau(nn)];
            increasing_profiles.re = [increasing_profiles.re, ensemble_profiles.re(nn)];
            increasing_profiles.lwc = [increasing_profiles.lwc, ensemble_profiles.lwc(nn)];
            increasing_profiles.Nc = [increasing_profiles.Nc, ensemble_profiles.Nc(nn)];
            increasing_profiles.time = [increasing_profiles.time, ensemble_profiles.time(nn)];

        end

    elseif dre_dz<0
        % then the profile decreases with increasing altitude, on average

        % grab just the variables of interest for from each vertical profile
        if exist("decreasing_profiles", "var")==false

            decreasing_profiles.altitude = ensemble_profiles.altitude(nn);
            decreasing_profiles.tau = ensemble_profiles.tau(nn);
            decreasing_profiles.re = ensemble_profiles.re(nn);
            decreasing_profiles.lwc = ensemble_profiles.lwc(nn);
            decreasing_profiles.Nc = ensemble_profiles.Nc(nn);
            decreasing_profiles.time = ensemble_profiles.time(nn);

        else

            decreasing_profiles.altitude = [decreasing_profiles.altitude, ensemble_profiles.altitude(nn)];
            decreasing_profiles.tau = [decreasing_profiles.tau, ensemble_profiles.tau(nn)];
            decreasing_profiles.re = [decreasing_profiles.re, ensemble_profiles.re(nn)];
            decreasing_profiles.lwc = [decreasing_profiles.lwc, ensemble_profiles.lwc(nn)];
            decreasing_profiles.Nc = [decreasing_profiles.Nc, ensemble_profiles.Nc(nn)];
            decreasing_profiles.time = [decreasing_profiles.time, ensemble_profiles.time(nn)];

        end

    end



end


% store the LWC threshold

increasing_profiles.lwc_threshold = LWC_threshold;
decreasing_profiles.lwc_threshold = LWC_threshold;




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

n_bins = 30; % number of segments the noramlized vertical component is broken up into

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



%% Make errorbar plot of the mean and standard deviation of effective radius at each vertical bin


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


%% Make errorbar plot of the median and standard deviation of effective radius at each vertical bin


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

% plot the standard deviation of the mean profile as an transparent area
% centered around the mean radius profile
x = [re_mean-re_std; flipud(re_mean + re_std)];
y = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(1,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the mean profile
plot(re_mean, bin_center, 'Color', mySavedColors(1, 'fixed'));
hold on

% plot the standard deviation of the mean profile as an transparent area
% centered around the mean radius profile
x = [re_median-re_std; flipud(re_median + re_std)];
y = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(2,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the median droplet profile
plot(re_median, bin_center, 'Color', mySavedColors(2, 'fixed'))

legend('','Mean', '', 'Median', 'Location', 'best')

set(gca, 'YDir', 'reverse')

grid on; grid minor
xlabel('$<r_e(z)>$', 'Interpreter','latex')
ylabel('Normalized Optical Depth')

set(gcf, 'Position', [0 0 650 600])



%% Plot the mean effective radius profile with a shaded area representing the standard deviation

figure;

% plot the standard deviation of the mean profile as an transparent area
% centered around the mean radius profile
x = [re_mean-re_std; flipud(re_mean + re_std)];
y = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(1,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the mean profile
plot(re_mean, bin_center, 'Color', mySavedColors(1, 'fixed'));

% plot a droplet profile fit over to see if it is adiabatic or not
profile_type = 'subadiabatic_aloft';
re_fit = create_droplet_profile2([re_mean(1), re_mean(end)], bin_center, 'optical_depth', profile_type);

hold on
plot(re_fit, bin_center, 'Color', 'k', 'LineWidth',1);

legend('','', 'subadiabatic', 'location', 'best')


set(gca, 'YDir', 'reverse')

grid on; grid minor
xlabel('<r_e>')
ylabel('Normalized Optical Depth')
title('Mean r_e')

set(gcf, 'Position', [0 0 650 600])


%% Make errorbar plot of the mean and standard deviation of liquid water content at each vertical bin


lwc_mean = zeros(n_bins, 1);
lwc_std = zeros(n_bins, 1);
bin_center = zeros(n_bins, 1);



for bb = 1:n_bins

    % compute the mean value for the current bin
    lwc_mean(bb) = mean(vertically_segmented_attributes{bb,2});       % microns - mean droplet radius value

    % compute the standard deviation of the current bin
    lwc_std(bb) = std(vertically_segmented_attributes{bb,2});         % microns - standard deviation

    % compute the bin center which is the tau location of the mean data
    bin_center(bb) = (bin_edges(bb+1) - bin_edges(bb))/2 + bin_edges(bb);

end

figure;
errorbar(lwc_mean, bin_center, lwc_std, lwc_std, 'horizontal', 'LineStyle','-',...
    'Marker','o', 'MarkerFaceColor', mySavedColors(1, 'fixed'), 'MarkerEdgeColor','k', ...
    'MarkerSize', 10, 'LineWidth', 2, 'Color', 'k')

set(gca, 'YDir', 'reverse')

grid on; grid minor
xlabel('Mean LWC')
ylabel('Normalized Optical Depth')

set(gcf, 'Position', [0 0 650 600])



%% Plot the mean liquid water content with a shaded area representing the standard deviation

figure;

% plot the standard deviation of the mean profile as an transparent area
% centered around the mean radius profile
x = [lwc_mean-lwc_std; flipud(lwc_mean + lwc_std)];
y = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(3,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the mean profile
plot(lwc_mean, bin_center, 'Color', mySavedColors(3, 'fixed'));



set(gca, 'YDir', 'reverse')

grid on; grid minor
xlabel('<LWC(z)>')
ylabel('Normalized Optical Depth')
title('Mean LWC')

set(gcf, 'Position', [0 0 650 600])


%% Make errorbar plot of the median and standard deviation of liquid water content at each vertical bin


lwc_median = zeros(n_bins, 1);
lwc_std = zeros(n_bins, 1);
bin_center = zeros(n_bins, 1);



for bb = 1:n_bins

    % compute the mean value for the current bin
    lwc_median(bb) = median(vertically_segmented_attributes{bb,2});       % microns - mean droplet radius value

    % compute the standard deviation of the current bin
    lwc_std(bb) = std(vertically_segmented_attributes{bb,2});         % microns - standard deviation

    % compute the bin center which is the tau location of the mean data
    bin_center(bb) = (bin_edges(bb+1) - bin_edges(bb))/2 + bin_edges(bb);

end

figure;
errorbar(lwc_median, bin_center, lwc_std, lwc_std, 'horizontal', 'LineStyle','-',...
    'Marker','o', 'MarkerFaceColor', mySavedColors(1, 'fixed'), 'MarkerEdgeColor','k', ...
    'MarkerSize', 10, 'LineWidth', 2, 'Color', 'k')

set(gca, 'YDir', 'reverse')

grid on; grid minor
xlabel('Median LWC')
ylabel('Normalized Optical Depth')

set(gcf, 'Position', [0 0 650 600])



%% Plot the median liquid water content with a shaded area representing the standard deviation

figure;

% plot the standard deviation of the mean profile as an transparent area
% centered around the mean radius profile
x = [lwc_median-lwc_std; flipud(lwc_median + lwc_std)];
y = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(4,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the median profile
plot(lwc_median, bin_center, 'Color', mySavedColors(4, 'fixed'));

% Plot a line between the first and last value of LWC
hold on
plot([lwc_median(2), lwc_median(end)], [bin_center(2), bin_center(end)], 'k', 'Linewidth', 1)

legend('','','Linear Fit', 'Location', 'best')



set(gca, 'YDir', 'reverse')

grid on; grid minor
xlabel('<LWC(z)>')
ylabel('Normalized Optical Depth')
title('Median LWC')

set(gcf, 'Position', [0 0 650 600])




%% Plot the mean and median Liquid Water content profiles

figure;

% plot the standard deviation of the mean profile as an transparent area
% centered around the mean radius profile
x = [lwc_mean-lwc_std; flipud(lwc_mean + lwc_std)];
y = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(6,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the mead profile
plot(lwc_mean, bin_center, 'Color', mySavedColors(6, 'fixed'));
hold on

% plot the standard deviation of the mean profile as an transparent area
% centered around the mean radius profile
x = [lwc_median-lwc_std; flipud(lwc_median + lwc_std)];
y = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(4,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the median droplet profile
plot(lwc_median, bin_center, 'Color', mySavedColors(4, 'fixed'))

legend('','Mean', '', 'Median', 'Location', 'best')

set(gca, 'YDir', 'reverse')

grid on; grid minor
xlabel('$<LWC(z)>$', 'Interpreter','latex')
ylabel('Normalized Optical Depth')

set(gcf, 'Position', [0 0 650 600])




%% Plot the mean and median Number Concentration profiles



Nc_median = zeros(n_bins, 1);
Nc_mean = zeros(n_bins, 1);
Nc_std = zeros(n_bins, 1);
bin_center = zeros(n_bins, 1);



for bb = 1:n_bins

    % compute the median value for the current bin
    Nc_median(bb) = median(vertically_segmented_attributes{bb,3});       % microns - mean droplet radius value

    % compute the mean value for the current bin
    Nc_mean(bb) = mean(vertically_segmented_attributes{bb,3});       % microns - mean droplet radius value

    % compute the standard deviation of the current bin
    Nc_std(bb) = std(vertically_segmented_attributes{bb,3});         % microns - standard deviation

    % compute the bin center which is the tau location of the mean data
    bin_center(bb) = (bin_edges(bb+1) - bin_edges(bb))/2 + bin_edges(bb);

end


figure;

% plot the standard deviation of the mean profile as an transparent area
% centered around the mean radius profile
x = [Nc_mean-Nc_std; flipud(Nc_mean + Nc_std)];
y = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(1,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the mean profile
plot(Nc_mean, bin_center, 'Color', mySavedColors(1, 'fixed'));
hold on

% plot the standard deviation of the mean profile as an transparent area
% centered around the mean radius profile
x = [Nc_median-Nc_std; flipud(Nc_median + Nc_std)];
y = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(2,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the median droplet profile
plot(Nc_median, bin_center, 'Color', mySavedColors(2, 'fixed'))

legend('','Mean', '', 'Median', 'Location', 'best')

set(gca, 'YDir', 'reverse')

grid on; grid minor
xlabel('$<N_c(z)>$ $(cm^{-3})$', 'Interpreter','latex')
ylabel('Normalized Optical Depth')

set(gcf, 'Position', [0 0 650 600])


%% Make a subplot showing the median profile of droplet size, liquid water content and number concentration

figure; 

% plot the median effective radius
subplot(1,3,1)

% plot the standard deviation of the median profile as an transparent area
% centered around the mean radius profile
x = [re_median-re_std; flipud(re_median + re_std)];
y = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(2,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the median droplet profile
plot(re_median, bin_center, 'Color', mySavedColors(2, 'fixed'))

set(gca, 'YDir', 'reverse')

grid on; grid minor
xlabel('$<r_e(z)>$ $(\mu m)$', 'Interpreter','latex')
ylabel('Normalized Optical Depth')


% plot the median liquid water content profile
subplot(1,3,2)

% plot the standard deviation of the median profile as an transparent area
% centered around the mean radius profile
x = [lwc_median-lwc_std; flipud(lwc_median + lwc_std)];
y = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(2,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the median droplet profile
plot(lwc_median, bin_center, 'Color', mySavedColors(2, 'fixed'))

set(gca, 'YDir', 'reverse')

grid on; grid minor
xlabel('$<LWC(z)>$ $(g/m^{3})$', 'Interpreter','latex')
ylabel('Normalized Optical Depth')

% set the figure title
title('Median Profiles')



% plot the median droplet number concentration
subplot(1,3,3)

% plot the standard deviation of the median profile as an transparent area
% centered around the mean radius profile
x = [Nc_median-Nc_std; flipud(Nc_median + Nc_std)];
y = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(2,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the median droplet profile
plot(Nc_median, bin_center, 'Color', mySavedColors(2, 'fixed'))

set(gca, 'YDir', 'reverse')

grid on; grid minor
xlabel('$<N_c(z)>$ $(cm^{-3})$', 'Interpreter','latex')
ylabel('Normalized Optical Depth')

% set the size of the figure
set(gcf, 'Position', [0 0 1200 625])


%% What is the average difference between cloud top and cloud bottom radii?

% I've broken up all of the data into vertical segments. 
% What is the average radius value for the bottom most segment?
re_top_median = mean(vertically_segmented_attributes{1,1});
re_bottom_median = median(vertically_segmented_attributes{end,1});

for nn = 1:length(ensemble_profiles.re)

    diff_top_bottom_1(nn) = ensemble_profiles.re{nn}(1) - ensemble_profiles.re{nn}(end);
    diff_top_bottom_2(nn) = ensemble_profiles.re{nn}(2) - ensemble_profiles.re{nn}(end);

end



