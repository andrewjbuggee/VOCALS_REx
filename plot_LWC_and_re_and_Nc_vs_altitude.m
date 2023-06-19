%% Plot the Liquid Water Content and the Effective Radius versus Altitude 


% By Andrew John Buggee
%%

function [] = plot_LWC_and_re_and_Nc_vs_altitude(vert_profiles, indices)


% Plot the number of curves within the vert_profiles structure

N_cuvres = length(indices);

legend_str = cell(1, N_cuvres);

figure;
for nn = 1:N_cuvres
    
    % First plot the LWC
    subplot(1,3,1); plot(vert_profiles.lwc{indices(nn)}, vert_profiles.altitude{indices(nn)})
    hold on

    % next plot the effective radius
    subplot(1,3,2); plot(vert_profiles.re{indices(nn)}, vert_profiles.altitude{indices(nn)}); 
    hold on

    % Lastly, plot the total droplet number concentration
    subplot(1,3,3); plot(vert_profiles.Nc{indices(nn)}, vert_profiles.altitude{indices(nn)}); 
    hold on

    legend_str{nn} = ['idx = ', num2str(indices(nn))];


end

% Make each subplot pretty
subplot(1,3,1)
grid on; grid minor; 
xlabel('LWC ($g/m^3$)', 'Interpreter','latex');
ylabel('Altitude ($m$)', 'Interpreter','latex');

% in the first subplot, define the indices being plotted
legend(legend_str, 'Interpreter','latex', 'Location','best')



subplot(1,3,2)
grid on; grid minor; 
xlabel('$r_e$ ($\mu m$)', 'Interpreter','latex')

subplot(1,3,3)
grid on; grid minor; 
xlabel('$N_c$ ($m^{-3}$)', 'Interpreter','latex')

% set plot size
set(gcf, 'Position', [0 0 1000 550])




end