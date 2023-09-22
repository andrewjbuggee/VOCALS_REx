%% Estimate the uncertainty of the effective radius based on the measurement uncertainty of 
% Droplet Measurement Technologies Cloud Droplet Probe (Lance et al.
% (2010)).


% By Andrew John Buggee

%%


function [droplet_radius_uncertainty] = cloud_droplet_probe_uncertainty_estimate(vocalsRex)



% -------------------------------------------------------------------
% ------------ Simplest 1st order uncertainty estimate --------------
% -------------------------------------------------------------------

% radius uncertainty is a function of the droplet number concentration
% and the radius itself
% (figure 12 of Lance et al. (2010))

number_concentration = [10, 100, 200, 300, 400];            % #-droplets/cm^3
re = [2.5, 3.75, 5, 10];                                % microns
re_uncertainty = [5, 20, 40, 50, 60;...
                  -1, 10, 20, 30, 38;...
                  -5, 5, 10, 20, 22;...
                  -5, -1, 5, 10, 12];                   % percent uncertainty

% 2D interpolation to get the CDP uncertainty
[Nc, Re] = meshgrid(number_concentration, re);

droplet_radius_uncertainty = interp2(Nc, Re, re_uncertainty, vocalsRex.Nc, vocalsRex.re, "spline");



end

