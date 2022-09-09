%% Scratch work to figure out why my calcualtions of VOCALS Rex data differ from Painemal and Zudema

% By Andrew John Buggee

%%  Grab VOCALS-REx Data

clear variables

scriptPlotting_wht;

% filename to open
filename = 'RF11.20081109.125700_213600.PNI.nc';
% load vocals rex data
vr = readVocalsRex(filename);

%% Grab MODIS data

% This is the modis data set taken on day 314 at decimal time 0.611
modis_folderName = './MODIS_data/2008_11_09/';


[modis,L1B_500m_fileName] = retrieveMODIS_data(modis_folderName);

modis_dataTime = (modis.time(1)*3600 + modis.time(2)*60)/86400;


%% Painemal and Zudema use decimal time
% This is used to identify when along a flights path they assumed a cloud
% and made their calculations

% Lets use their decimal indicator to find where in the data to compare

pz_time = 0.6120;

% convert the start time to decimal time
sec_per_hour = 3600;                                                                    % sec
sec_per_min = 60;                                                                       % sec
sec_per_day = 86400;                                                                    % sec

startTime_dec = (vr.startTime(1)*3600 + vr.startTime(2)*60)/86400;                      % fraction of a day 

pz_secSinceStart = floor((pz_time - startTime_dec)*sec_per_day);                               % seconds since flight started up to the time indicated by painemal and zudema

modis_secSinceStart = floor((modis_dataTime - startTime_dec)*sec_per_day);


% Define the number of discrete points in time you wish to look at. The
% center point will be at the start time calculated above
windowLength = 100;

% Using the time defined above, find the location in space closest to
% the airplanes location

% Using lat and long with can minimize the euclidean norm
[minDist, index_minDist] = min(sqrt((modis.geo.lat - vr.latitude(pz_secSinceStart)).^2 + (modis.geo.long - vr.longitude(pz_secSinceStart)).^2), [], 'all');

% grab the value for droplet size at the location closest to the planes
% locaiton at the time specified above
re_16_modis = modis.cloud.effRadius16(index_minDist);
re_17_modis = modis.cloud.effRadius17(index_minDist);

%% Make a plot of the number concentration and effective radius aat the time defined above

% I want a subplot with the number concentration and altitude, and the
% effective radius with altitude
nice_blue = [0 0.4470 0.741];

figure; subplot(1,2,1)
semilogy(vr.time, vr.re); xlabel('Time (sec)'); ylabel('r_{e} (\mum)')
grid on; grid minor; hold on; 

% Plot the modis droplet estimate as a constant vertical line
xl0 = yline(re_16_modis,'--',['MODIS $$r_{e} = $$',num2str(re_16_modis), '$$\mu m$$'], 'Fontsize',18, 'Interpreter','latex','LineWidth',2,'Color',nice_blue);
xl0.LabelVerticalAlignment = 'top';

% now plot the altitude of the plane on the right hand side
yyaxis right; plot(vr.time, vr.altitude); ylabel('Altitude (m agl)')
% zoom in on a particular window

xlim([pz_secSinceStart - windowLength/2, pz_secSinceStart + windowLength/2]);



% plot a vertical line marking the bottom of the cloud
% This occurs when the LWC is atleast 0.03 c/m^3
lwc_lim = 0.03;                                                 % grams/m^3 - lower limit defining cloud base

LWC_window_index = pz_secSinceStart - windowLength/4 : pz_secSinceStart + windowLength/4;
LWC_window = vr.lwc(LWC_window_index);
[minVal, Ind_minVal] = min(abs(LWC_window - lwc_lim));

% Plot a vertical line!
xl = xline(vr.time(LWC_window_index(1) + Ind_minVal -1),'--','$$LWC = 0.03\; g/m^{3}$$', 'Fontsize',18, 'Interpreter','latex','LineWidth',2);
xl.LabelVerticalAlignment = 'top';


subplot(1,2,2)
semilogy(vr.time, vr.Nc); xlabel('Time (sec)'); ylabel('N_{c} (cm^{-3})')
grid on; grid minor; hold on; yyaxis right; plot(vr.time, vr.altitude); ylabel('Altitude (m agl)')
% zoom in on a particular window
xlim([pz_secSinceStart - windowLength/2, pz_secSinceStart + windowLength/2]);

% Plot a vertical line!
xl = xline(vr.time(LWC_window_index(1) + Ind_minVal -1),'--','$$LWC = 0.03\; g/m^{3}$$', 'Fontsize',18, 'Interpreter','latex','LineWidth',2);
xl.LabelVerticalAlignment = 'top';
xl.LabelHorizontalAlignment = "left";

set(gcf, 'Position', [0 0 1300 400])

%% Using Painemal and Zudema definition for cloud top and bottom

lwc_lim = 0.03;                                                 % grams/m^3 - lower limit defining cloud base

% I want a subplot with the number concentration and liquid water content,
% and the effective radius with liquid water content

figure; subplot(1,2,1)
semilogy(vr.time, vr.re); xlabel('Time (sec)'); ylabel('r_{e} (\mum)')
grid on; grid minor; hold on; yyaxis right; plot(vr.time, vr.lwc); ylabel('LWC (g/m^{3})')
% zoom in on a particular window

xlim([pz_secSinceStart - windowLength/2, pz_secSinceStart + windowLength/2]);

subplot(1,2,2)
semilogy(vr.time, vr.Nc); xlabel('Time (sec)'); ylabel('N_{c} (cm^{-3})')
grid on; grid minor; hold on; yyaxis right; plot(vr.time, vr.lwc); ylabel('LWC (g/m^{3})')
% zoom in on a particular window
xlim([pz_secSinceStart - windowLength/2, pz_secSinceStart + windowLength/2]);

set(gcf, 'Position', [0 0 1300 400])

%% Compute the Liquid water path using the threshold defined by Painemal and Zudema

% Let's find all times when lwc is greater than the threshold above
% This will define each individual cloud in our time series

index_lwc = vr.lwc>lwc_lim;

index_minLength = 10;                            % each segment has to have atleast 3 data points to be counted as a cloud

% We want to find where these logical ones start and end. A 1 next to a 0
% will have a different of 1. A 1 next to a 1 or a 0 next to a 0 will have
% a difference of 0
% We add zeros to either end of our vector so that the index is the same as
% the original vector. The difference vector shrinks the length of the
% input by 1

% So the odd entries are the start of the cloud, the even entires are the
% end of the cloud

b = find(diff([0 index_lwc 0]));
c = b(2:2:end) - b(1:2:end);

d = find(c>=index_minLength);

cloud_index = cell(length(d),1);

for ii = 1:length(d)
    cloud_index{ii} = b(2*d(ii) - 1):1:b(2*d(ii))-1;
end



%% Lets plot the time series of all clouds found in the above section

figure;

colororder('default')

for ii = 1:length(cloud_index)
    hold on
    plot(vr.time(cloud_index{ii}), vr.re(cloud_index{ii}),'Color',[0, 0.4470, 0.7410]); 
   
end
xlabel('Time (sec)'); ylabel('r_{e} (\mum)')

yyaxis right

for ii = 1:length(cloud_index)
    hold on; 
    plot(vr.time(cloud_index{ii}), vr.altitude(cloud_index{ii}), 'Color',[0.8500, 0.3250, 0.0980],'LineStyle','-'); 
end
ylabel('Altitude (m agl)')
set(gca,'ycolor',[0.8500, 0.3250, 0.0980]) 

grid on; grid minor; 

set(gcf, 'Position', [0 0 1300 400])




%% Make plot of the size distribution at time = 6256 when the number concentration is constant


figure;
histogram(vr.nr(:,pz_secSinceStart), vr.drop_bins); xlabel('Droplet Size (\mum)'); ylabel('n(r) (cm^{-3})')
grid on; grid minor; hold on;
set(gca,'YScale','log')
