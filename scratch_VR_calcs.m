%% Scratch work to figure out why my calcualtions of VOCALS Rex data differ from Painemal and Zudema

% By Andrew John Buggee
%%

clear variables

scriptPlotting_wht;

% filename to open
filename = 'RF11.20081109.125700_213600.PNI.nc';
% load vocals rex data
vr = readVocalsRex(filename);

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

% Define the number of discrete points in time you wish to look at. The
% center point will be at the start time calculated above
windowLength = 100;

%% Make a plot of the number concentration and effective radius around time = 6256 

% I want a subplot with the number concentration and altitude, and the
% effective radius with altitude

figure; subplot(1,2,1)
semilogy(vr.time, vr.re); xlabel('Time (sec)'); ylabel('r_{e} (\mum)')
grid on; grid minor; hold on; yyaxis right; plot(vr.time, vr.altitude); ylabel('Altitude (m agl)')
% zoom in on a particular window

xlim([pz_secSinceStart - windowLength/2, pz_secSinceStart + windowLength/2]);

subplot(1,2,2)
semilogy(vr.time, vr.Nc); xlabel('Time (sec)'); ylabel('N_{c} (cm^{-3})')
grid on; grid minor; hold on; yyaxis right; plot(vr.time, vr.altitude); ylabel('Altitude (m agl)')
% zoom in on a particular window
xlim([pz_secSinceStart - windowLength/2, pz_secSinceStart + windowLength/2]);

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








%% Make plot of the size distribution at time = 6256 when the number concentration is constant


figure;
histogram(vr.nr(:,pz_secSinceStart), vr.drop_bins); xlabel('Droplet Size (\mum)'); ylabel('n(r) (cm^{-3})')
grid on; grid minor; hold on;
