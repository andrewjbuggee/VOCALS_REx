%% Plot total number concentration versus time and altitude versus time for entire data set


figure; 
semilogy(double(vocalsRex.time)./3600, vocalsRex.total_Nc); 
grid on; grid minor; 
xlabel('Time (hours)','Interpreter','latex')
ylabel('Total $N_c$  $(cm^{3})$','Interpreter','latex')
title('VOCALS-REx flight Data', 'Interpreter','latex')

hold on;
yyaxis right; 
plot(double(vocalsRex.time)./3600, vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

set(gcf, 'Position',[0 0 1300, 450])



%% Plot effective radius versus time and altitude versus time for entire data set


figure; 
semilogy(double(vocalsRex.time)./3600, vocalsRex.re); 
grid on; grid minor; 
xlabel('Time (hours)','Interpreter','latex')
ylabel('$r_e$  $(\mu m)$','Interpreter','latex')
title('VOCALS-REx flight Data', 'Interpreter','latex')

hold on;
yyaxis right; 
plot(double(vocalsRex.time)./3600, vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

set(gcf, 'Position',[0 0 1300, 450])


%% Plot total number concentration versus time and altitude versus time for entire data set
%  Below this, create a subplot with effective radius versus time and
%  altitude versus time


figure; subplot(2,1,1)
semilogy(double(vocalsRex.time)./3600, vocalsRex.total_Nc); 
grid on; grid minor; 
ylabel('Total $N_c$  $(cm^{3})$','Interpreter','latex')
title('VOCALS-REx flight Data', 'Interpreter','latex')

hold on;
yyaxis right; 
plot(double(vocalsRex.time)./3600, vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

% grab current axes
ax1 = gca;

subplot(2,1,2)
semilogy(double(vocalsRex.time)./3600, vocalsRex.re); 
grid on; grid minor; 
xlabel('Time (hours)','Interpreter','latex')
ylabel('$r_e$ $(\mu m)$','Interpreter','latex')

hold on;
yyaxis right; 
plot(double(vocalsRex.time)./3600, vocalsRex.altitude)
ylabel('Altitude ($m$)','Interpreter','latex')

% grab current axes
ax2 = gca;

linkaxes([ax1 ax2],'x');


set(gcf, 'Position',[0 0 1300, 600])
