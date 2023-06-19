%% Plot the droplet distribution as a histogram for multiple any number of heights within a cloud




% By Andrew John Buggee
%%

function [] = plot_dropDistribution_topAndBottom_hist(vocalsRex, vert_profs, indexs2plot)



% Define the min and max radius values to plot
r_min = 0;      % microns
r_max = 30;     % microns



for nn = 1:length(indexs2plot)


    f1 = figure;

    % define the index for cloud top and cloud bottom
    times2plot = [1, numel(vert_profs.Nc{indexs2plot(nn)})];
    
    % set the effective radius vector
    re = zeros(1, 2);

    for tt = 1:length(times2plot)

        % plot some droplet distribution within the cloud
        index_time = vert_profs.time{indexs2plot(nn)}(times2plot(tt))+1;


        % Compute the effective radius for the two distributions and plot it as a solid vertical line
        re(tt) = vert_profs.re{indexs2plot(nn)}(times2plot(tt));

        % Plot the distribution at cloud bottom first
        h1 = histogram('BinEdges',vocalsRex.drop_radius_bin_edges ,'BinCounts',vocalsRex.Nc(:,index_time));
        h1.FaceColor = mySavedColors(tt, 'fixed');
        h1.FaceAlpha = 0.7;
        h1.EdgeAlpha = 1;
        hold on
        xline(re(tt),':', 'LineWidth',4, 'Color',mySavedColors(tt, 'fixed'))

    end
    

    % what profile are we plotting?
    title(['Vertical Profile ', num2str(indexs2plot(nn))], 'Interpreter','latex')

    % set axes limits and labels
    xlabel('Droplet Radius ($\mu m$)', 'Interpreter','latex', 'FontSize',32);
    ylabel('$n(r)$ ($cm^{-3}$)', 'Interpreter','latex', 'FontSize',32);
    grid on; grid minor; hold on;
    xlim([r_min, r_max])
    ylim([10^(-2) 10^2])
    set(gca, 'YScale', 'log')
    set(gcf, 'Position',[0 0 1000, 600])
    
    
    legend('Cloud Bottom', ['$r_e$ = ',num2str(round(re(1),2)), ' $\mu m$'], 'Cloud Top',...
            ['$r_e$ = ',num2str(round(re(2),2)), ' $\mu m$'], 'Location','best', 'Interpreter','latex')



end











end