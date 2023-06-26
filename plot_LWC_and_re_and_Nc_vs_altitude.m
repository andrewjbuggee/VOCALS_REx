%% Plot the Liquid Water Content and the Effective Radius versus Altitude


% By Andrew John Buggee
%%

function [] = plot_LWC_and_re_and_Nc_vs_altitude(vert_profiles, indices, normalize_altitude)

% Check to make sure there are 3 inputs


if nargin~=3
    error([newline,'Wrong number of inputs. Need 3: vertical profiles, indices to plot, normalization flag', newline])
end



% Plot the number of curves within the vert_profiles structure

N_cuvres = length(indices);

legend_str = cell(1, N_cuvres);

figure;
for nn = 1:N_cuvres

    % if normalize altitude is true, all altitude vectors will be
    % normalized between 0 and 1

    if normalize_altitude==true

        norm_alt = (vert_profiles.altitude{indices(nn)} - vert_profiles.altitude{indices(nn)}(1))./...
            (vert_profiles.altitude{indices(nn)}(end) - vert_profiles.altitude{indices(nn)}(1));

        % First plot the LWC
        ax1 = subplot(1,3,1); plot(vert_profiles.lwc{indices(nn)}, norm_alt);
        hold on

        % next plot the effective radius
        ax2 = subplot(1,3,2); plot(vert_profiles.re{indices(nn)}, norm_alt);
        hold on

        % Lastly, plot the total droplet number concentration
        ax3 = subplot(1,3,3); plot(vert_profiles.Nc{indices(nn)}, norm_alt);
        hold on

    else

        % First plot the LWC
        ax1 = subplot(1,3,1); plot(vert_profiles.lwc{indices(nn)}, vert_profiles.altitude{indices(nn)});
        hold on

        % next plot the effective radius
        ax2 = subplot(1,3,2); plot(vert_profiles.re{indices(nn)}, vert_profiles.altitude{indices(nn)});
        hold on

        % Lastly, plot the total droplet number concentration
        ax3 = subplot(1,3,3); plot(vert_profiles.Nc{indices(nn)}, vert_profiles.altitude{indices(nn)});
        hold on

    end

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
% include a title in the middle plot
title(['LWC $\geq$ ', num2str(vert_profiles.lwc_threshold),' $g/m^{3}$'], 'interpreter', 'latex')


subplot(1,3,3)
grid on; grid minor;
xlabel('$N_c$ ($m^{-3}$)', 'Interpreter','latex')

% set plot size
set(gcf, 'Position', [0 0 1000 550])

% link the yaxes so that they all have the same bounds
linkaxes([ax1 ax2 ax3],'y')




end