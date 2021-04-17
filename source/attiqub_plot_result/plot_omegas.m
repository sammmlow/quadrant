%%#######################################################################
% #######################################################################
% ###                                                                 ###
% ###  @@@@@@@@ @@@@@@@ @@@@@@@ @@@@@ @@@@@@@@@ @@@    @@@ @@@@@@@@   ###
% ###  @@@@@@@@ @@@@@@@ @@@@@@@ @@@@@ @@    (@@ @@@    @@@ @@@  @@@   ###
% ###  @@    @@   @@@     @@@    @@@  @@    (@@ @@@    @@@ @@@@@@@@@  ###
% ###  @@@@@@@@   @@@     @@@    @@@  @@@@@@@@@ @@@    @@@ @@@@@@@@@  ###
% ###  @@    @@   @@@     @@@   @@@@@ @@@@@@@@@ @@@@@@@@@@ @@@   @@@  ###
% ###  @@    @@   @@@     @@@   @@@@@     (@@@  @@@@@@@@@@ @@@@@@@@@  ###
% ###                                                                 ###
% ### Project ATTIQUB: Plotting Results Toolbox                       ###
% ### Plots the body angular velocity (rad/s) as seen in the          ###
% ### reference frame and the inertial frame.                         ###
% ###                                                                 ###
% ### By Samuel Low (13-01-2021), DSO National Laboratories           ###
% ###                                                                 ###
% #######################################################################
% #######################################################################

function plot_omegas( t, OHM_BR_ARRAY, OHM_BN_ARRAY )

%% Inputs:
% ---> t: time axis
% ---> OHM_BR_ARRAY: Matrix of body-reference angular velocities (tx3)
% ---> OHM_BN_ARRAY: Matrix of body-inertial angular velocities (tx3)

figure('units','normalized','outerposition',[0 0 1 1])

%% Plot the body-reference angular velocities.
subplot(1,2,1); grid on; hold on;
OHM_BR_ARRAY = OHM_BR_ARRAY';
plot(t, OHM_BR_ARRAY(1,:), 'LineWidth', 1);
plot(t, OHM_BR_ARRAY(2,:), 'LineWidth', 1);
plot(t, OHM_BR_ARRAY(3,:), 'LineWidth', 1);
plot(t, OHM_BR_ARRAY(4,:), 'LineWidth', 1);
title('Angular Velocity Body-to-Reference (rad/s)');
legend('\omega_{BR-1}','\omega_{BR-2}','\omega_{BR-3}','\omega_{BR-Norm}');

%% Plot the body-inertial angular velocities.
subplot(1,2,2); grid on; hold on;
OHM_BN_ARRAY = OHM_BN_ARRAY';
plot(t, OHM_BN_ARRAY(1,:), 'LineWidth', 1);
plot(t, OHM_BN_ARRAY(2,:), 'LineWidth', 1);
plot(t, OHM_BN_ARRAY(3,:), 'LineWidth', 1);
plot(t, OHM_BN_ARRAY(4,:), 'LineWidth', 1);
title('Angular Velocity Body-to-Inertial (rad/s)');
legend('\omega_{BN-1}','\omega_{BN-2}','\omega_{BN-3}','\omega_{BN-Norm}');

end

