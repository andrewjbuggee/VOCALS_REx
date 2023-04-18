%% Load some VOCALS-REx Data

% Oct-18-2008 Data
filename = 'RF02.20081018.130300_213000.PNI.nc';




% Load data
vocalsRex = readVocalsRex(filename);

%% Plot the VOCALS-REx flight path

% Plot MODIS measured relfectance at 650nm with VOCALS flight path
figure; geoscatter(vocalsRex.latitude, vocalsRex.longitude, 10,'.')
set(gca, 'FontSize',25)
set(gca, 'FontWeight', 'bold')
set(gcf, 'Position', [0 0 800 800])
title('VOCALS-REx flight path','Interpreter','latex', 'FontSize', 40)
geolimits([min(vocalsRex.latitude)-5 max(vocalsRex.latitude)+5],[min(vocalsRex.longitude)-5 max(vocalsRex.longitude)+5])



%% Plot total number concentration versus UTC time and altitude versus time for entire data set

UTC_starttime = vocalsRex.startTime(1) + vocalsRex.startTime(2)/60;   % hours.decimalhours in UTC format

figure; 
semilogy(UTC_starttime + double(vocalsRex.time)./3600, vocalsRex.total_Nc); 
grid on; grid minor; 
xlabel('UTC Time','Interpreter','latex')
ylabel('Total $N_c$  $(cm^{3})$','Interpreter','latex')
title('VOCALS-REx flight Data', 'Interpreter','latex')

hold on;
yyaxis right; 
plot(UTC_starttime + double(vocalsRex.time)./3600, vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

set(gcf, 'Position',[0 0 1300, 450])



%% Plot effective radius versus UTC time and altitude versus time for entire data set

UTC_starttime = vocalsRex.startTime(1) + vocalsRex.startTime(2)/60;   % hours.decimalhours in UTC format



figure; 
semilogy(UTC_starttime + double(vocalsRex.time)./3600, vocalsRex.re); 
grid on; grid minor; 
xlabel('UTC Time','Interpreter','latex')
ylabel('$r_e$  $(\mu m)$','Interpreter','latex')
title('VOCALS-REx flight Data', 'Interpreter','latex')

hold on;
yyaxis right; 
plot(UTC_starttime + double(vocalsRex.time)./3600, vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

set(gcf, 'Position',[0 0 1300, 450])



%% Plot total number concentration versus time and altitude versus time for entire data set


figure; 
semilogy((vocalsRex.time), vocalsRex.total_Nc); 
grid on; grid minor; 
xlabel('Time (sec since takeoff)','Interpreter','latex')
ylabel('Total $N_c$  $(cm^{3})$','Interpreter','latex')
title('VOCALS-REx flight Data', 'Interpreter','latex')

hold on;
yyaxis right; 
plot((vocalsRex.time), vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

set(gcf, 'Position',[0 0 1300, 450])


%% Plot total number concentration versus time and altitude versus time for entire data set
%  Below this, create a subplot with effective radius versus time and
%  altitude versus time
% ----- VOCALS TIME -------



figure; subplot(2,1,1)
semilogy(double(vocalsRex.time), vocalsRex.total_Nc); 
grid on; grid minor; 
ylabel('Total $N_c$  $(cm^{3})$','Interpreter','latex')
title('VOCALS-REx flight Data', 'Interpreter','latex')

hold on;
yyaxis right; 
plot(double(vocalsRex.time), vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

% grab current axes
ax1 = gca;

subplot(2,1,2)
semilogy(double(vocalsRex.time), vocalsRex.re); 
grid on; grid minor; 
xlabel('Time (hours)','Interpreter','latex')
ylabel('$r_e$ $(\mu m)$','Interpreter','latex')

hold on;
yyaxis right; 
plot(double(vocalsRex.time), vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

% grab current axes
ax2 = gca;

linkaxes([ax1 ax2],'x');


set(gcf, 'Position',[0 0 1300, 600])



%% Plot total number concentration versus time and altitude versus time for entire data set
%  Below this, create a subplot with effective radius versus time and
%  altitude versus time

UTC_starttime = vocalsRex.startTime(1) + vocalsRex.startTime(2)/60;   % hours.decimalhours in UTC format



figure; subplot(2,1,1)
semilogy(UTC_starttime + double(vocalsRex.time)./3600, vocalsRex.total_Nc); 
grid on; grid minor; 
ylabel('Total $N_c$  $(cm^{3})$','Interpreter','latex')
title('VOCALS-REx flight Data', 'Interpreter','latex')

hold on;
yyaxis right; 
plot(UTC_starttime + double(vocalsRex.time)./3600, vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

% grab current axes
ax1 = gca;

subplot(2,1,2)
semilogy(UTC_starttime + double(vocalsRex.time)./3600, vocalsRex.re); 
grid on; grid minor; 
xlabel('Time (hours)','Interpreter','latex')
ylabel('$r_e$ $(\mu m)$','Interpreter','latex')

hold on;
yyaxis right; 
plot(UTC_starttime + double(vocalsRex.time)./3600, vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

% grab current axes
ax2 = gca;

linkaxes([ax1 ax2],'x');


set(gcf, 'Position',[0 0 1300, 600])



%% Plot total number concentration versus time and LWC versus time on two subplots
% Overlay the airplane altitude on both


UTC_starttime = vocalsRex.startTime(1) + vocalsRex.startTime(2)/60;   % hours.decimalhours in UTC format


figure; subplot(2,1,1)
semilogy(UTC_starttime + double(vocalsRex.time)./3600, vocalsRex.total_Nc); 
grid on; grid minor; 
ylabel('Total $N_c$  $(cm^{3})$','Interpreter','latex')
title('VOCALS-REx flight Data', 'Interpreter','latex')

hold on;
yyaxis right; 
plot(UTC_starttime + double(vocalsRex.time)./3600, vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

% grab current axes
ax1 = gca;

subplot(2,1,2)
semilogy(UTC_starttime + double(vocalsRex.time)./3600, vocalsRex.lwc); 
grid on; grid minor; 
xlabel('Time (hours)','Interpreter','latex')
ylabel('$LWC$ $(g/m^{3})$','Interpreter','latex')

hold on;
yyaxis right; 
plot(UTC_starttime + double(vocalsRex.time)./3600, vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

% grab current axes
ax2 = gca;

linkaxes([ax1 ax2],'x');


set(gcf, 'Position',[0 0 1275, 600])



%% Find data for when the plane ascends or descends cleanly through a cloud deck

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

% Plot the derivative 
UTC_starttime = vocalsRex.startTime(1) + vocalsRex.startTime(2)/60;   % hours.decimalhours in UTC format

time_UTC = UTC_starttime + double(vocalsRex.time)./3600;            % Time in UTC


figure; 
plot(time_UTC(2:end), dz_dt)
hold on;
plot(time_UTC(2:end), dz_dt_mean)
grid on; grid minor; 
xlabel('UTC Time','Interpreter','latex')
ylabel('$dz/dt$ $(m/s)$','Interpreter','latex')
title('VOCALS-REx flight Data', 'Interpreter','latex')

hold on;
yyaxis right; 
plot(time_UTC, vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

set(gcf, 'Position',[0 0 1300, 450])




% Link axes with total number concentration
figure; subplot(2,1,1)
plot(time_UTC(2:end), dz_dt_mean)
grid on; grid minor; 
xlabel('UTC Time','Interpreter','latex')
ylabel('$dz/dt$ $(m/s)$','Interpreter','latex')
title('VOCALS-REx flight Data', 'Interpreter','latex')

hold on;
yyaxis right; 
plot(time_UTC, vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

% grab current axes
ax1 = gca;

% Plot total number concentration versus time and altitude
subplot(2,1,2)
semilogy(time_UTC, vocalsRex.total_Nc); 
grid on; grid minor; 
xlabel('Time (hours)','Interpreter','latex')
ylabel('Total $N_c$ $(cm^{-3})$','Interpreter','latex')

hold on;
yyaxis right; 
plot(time_UTC, vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

% grab current axes
ax2 = gca;

linkaxes([ax1 ax2],'x');


set(gcf, 'Position',[0 0 1300, 450])


%% Plot the individual vertical profiles versus their original time

[vert_profs] = find_verticalProfiles_VOCALS_REx(vocalsRex);

% Let's see what they look like
figure;
for ii=1:length(vert_profs.Nc)
    semilogy(vert_profs.time_UTC{ii}, vert_profs.Nc{ii}, '.')
    hold on;
end
grid on; grid minor
xlabel('Time (UTC)','Interpreter','latex')
ylabel('Total $N_c$  $(cm^{3})$','Interpreter','latex')
title('Vertical profiles only', 'Interpreter','latex')
set(gcf, 'Position',[0 0 1300, 450])




%% Plot the individual vertical profiles normalized by cloud top height
% Place them all on the same panel


% Let's see what they look like
figure;


% First plot the LWC
subplot(1,3,1)

for ii=1:length(vert_profs.lwc)

    % Normalize the data by diving the altitude by the cloud top height
    z_norm = (vert_profs.altitude{ii} - min(vert_profs.altitude{ii}))./(max(vert_profs.altitude{ii})- min(vert_profs.altitude{ii}));

    l = semilogx(vert_profs.lwc{ii}, z_norm, 'LineWidth',3);
    % Set the transparency to 50%
    l.Color(4) = 0.5;

    hold on;

end

grid on; grid minor
ylabel('$Z/Z_{top}$','Interpreter','latex')
xlabel('$LWC$  $(g/m^{3})$','Interpreter','latex')


% Plot the effective radius next
subplot(1,3,2)

for ii=1:length(vert_profs.re)

    % Normalize the data by diving the altitude by the cloud top height
    z_norm = (vert_profs.altitude{ii} - min(vert_profs.altitude{ii}))./max(vert_profs.altitude{ii});

    semilogx(vert_profs.re{ii}, z_norm, 'LineWidth',2)
    hold on;

end

grid on; grid minor
xlabel('$r_e$  $(\mu m)$','Interpreter','latex')


% Plot the total droplet number concentration last
subplot(1,3,3)

for ii=1:length(vert_profs.Nc)

    % Normalize the data by diving the altitude by the cloud top height
    z_norm = (vert_profs.altitude{ii} - min(vert_profs.altitude{ii}))./max(vert_profs.altitude{ii});

    semilogx(vert_profs.Nc{ii}, z_norm, 'LineWidth',2)
    hold on;

end

grid on; grid minor
xlabel('$N_c$  $(cm^{-3})$','Interpreter','latex')



set(gcf, 'Position',[0 0 1300, 450])


%% Get all indices for brushed data

% brushedData is the data grabbed from the plot

index_v = [];

for ii=1:size(brushedData,1)

    is_true = brushedData(ii,1)==vocalsRex.time;

    index_v = [index_v, is_true];
    
end

% sum all columns up into one logical vector!
index_v = logical(sum(index_v,2));


% Save all new data

% create structure with just vertical profiles
vert_profs.time = vocalsRex.time(index_v);
vert_profs.total_Nc = vocalsRex.total_Nc(index_v)';
vert_profs.lwc = vocalsRex.lwc(index_v)';
vert_profs.latitude = vocalsRex.latitude(index_v);
vert_profs.longitude = vocalsRex.longitude(index_v);
vert_profs.altitude = vocalsRex.altitude(index_v)';
vert_profs.re = vocalsRex.re(index_v)';


%% Plot just the vertical profiles taken from above
% Plot their LWC versus altitude


figure; 
plot(vert_profs.lwc, vert_profs.altitude/1e3, '.'); 
grid on; grid minor; 
xlabel('LWC $(g/m^3)$','Interpreter','latex')
ylabel('Altitude $(km)$','Interpreter','latex')
title('VOCALS-REx flight Data', 'Interpreter','latex')
set(gcf, 'Position',[0 0 1300, 450])


%% Plot droplet distribution at the base of some vert profile

% Plot one distribuiton at the bottom of the cloud
vocals_time1 = 5916-1;
index_time1 = vocalsRex.time==vocals_time1;

figure;
h = histogram('BinEdges',vocalsRex.drop_radius_bin_edges ,'BinCounts',vocalsRex.Nc(:,index_time1)); 
h.FaceColor = mySavedColors(11, 'fixed');
h.FaceAlpha = 0.5;
h.EdgeAlpha = 1;
xlabel('Droplet Radius ($\mu m$)', 'Interpreter','latex');
ylabel('$N_c$ ($cm^{-3}$)', 'Interpreter','latex');
grid on; grid minor; hold on;
xlim([0,30])
ylim([10^(-2) 10^2])
set(gca, 'YScale', 'log')
title('VOCALS-REx flight Data', 'Interpreter','latex')
set(gcf, 'Position',[0 0 1300, 450])

% Compute the effective radius and plot it as a solid vertical line
re1 = 1e4 * sum((vocalsRex.drop_radius_bin_center*1e-4).^3 .* vocalsRex.Nc(:,index_time1)')./...
      sum((vocalsRex.drop_radius_bin_center*1e-4).^2 .* vocalsRex.Nc(:,index_time1)');
xline(re1,'k--', 'LineWidth',2)


% plot another distribution at the top of the cloud
vocals_time2 = 5964-1;
index_time2 = vocalsRex.time==vocals_time2;

figure;
h = histogram('BinEdges',vocalsRex.drop_radius_bin_edges ,'BinCounts',vocalsRex.Nc(:,index_time2)); 
h.FaceColor = mySavedColors(14, 'fixed');
h.FaceAlpha = 0.5;
h.EdgeAlpha = 1;
xlabel('Droplet Radius ($\mu m$)', 'Interpreter','latex');
ylabel('$N_c$ ($cm^{-3}$)', 'Interpreter','latex');
grid on; grid minor; hold on;
xlim([0,30])
ylim([10^(-2) 10^2])
set(gca, 'YScale', 'log')
title('VOCALS-REx flight Data', 'Interpreter','latex')
set(gcf, 'Position',[0 0 1300, 450])

% Compute the effective radius and plot it as a solid vertical line
re2 = 1e4 * sum((vocalsRex.drop_radius_bin_center*1e-4).^3 .* vocalsRex.Nc(:,index_time2)')./...
      sum((vocalsRex.drop_radius_bin_center*1e-4).^2 .* vocalsRex.Nc(:,index_time2)');
xline(re2,'k--', 'LineWidth',2)



% Plot both of them on the same figure!


f1 = figure;
h1 = histogram('BinEdges',vocalsRex.drop_radius_bin_edges ,'BinCounts',vocalsRex.Nc(:,index_time1)); 
h1.FaceColor = mySavedColors(11, 'fixed');
h1.FaceAlpha = 0.7;
h1.EdgeAlpha = 1;
hold on
xline(re1,'k--', 'LineWidth',4)

h2 = histogram('BinEdges',vocalsRex.drop_radius_bin_edges ,'BinCounts',vocalsRex.Nc(:,index_time2)); 
h2.FaceColor = mySavedColors(14, 'fixed');
h2.FaceAlpha = 0.4;
h2.EdgeAlpha = 1;
xline(re2,'k--', 'LineWidth',4)

xlabel('Droplet Radius ($\mu m$)', 'Interpreter','latex', 'FontSize',32);
ylabel('$n(r)$ ($\mu m^{-1} \, cm^{-3}$)', 'Interpreter','latex', 'FontSize',32);
grid on; grid minor; hold on;
xlim([0,25])
ylim([10^(-2) 10^2])
set(gca, 'YScale', 'log')
set(gcf, 'Position',[0 0 1000, 600])

% Label first xline
annotation(f1,'textbox',[0.229125 0.917777777777778 0.0911875000000002 0.0711111111111111],...
    'String','$r_e = 4.7 \mu m$',...
    'Interpreter','latex',...
    'FontSize',18,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Label second x line
annotation(f1,'textbox',[0.317406250000002 0.917777777777778 0.0911875000000001 0.0711111111111111],...
    'String','$r_e = 7.2 \mu m$',...
    'Interpreter','latex',...
    'FontSize',18,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');


