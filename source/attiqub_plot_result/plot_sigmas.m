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
% ### Plots the body attitude as seen in the reference frame and      ###
% ### the inertial frame.                                             ###
% ###                                                                 ###
% ### By Samuel Low (13-01-2021), DSO National Laboratories           ###
% ###                                                                 ###
% #######################################################################
% #######################################################################

function plot_sigmas( satAttCoord, t, ATT_BR_ARRAY, ATT_BN_ARRAY )

%% Inputs:
% ---> t: time axis
% ---> ATT_BR_ARRAY: Matrix of body-reference satellite attitudes (txN)
% ---> ATT_BN_ARRAY: Matrix of body-inertial satellite attitudes (txN)
% ---> satAttCoord: Integer 1, 2 or 3
%       |-> 1 --> Quarternions, or Euler Parameters (1x4 set)
%       |-> 2 --> Classical Rodrigues Parameters (1x3 set)
%       |-> 3 --> Modified Rodrigues Parameters (1x3 set)

figure('units','normalized','outerposition',[0 0 1 1])

%% Get the attitude identified for the plot title.
if satAttCoord == 1
    ATT_TITLE = '(QRT)';
elseif satAttCoord == 2
    ATT_TITLE = '(CRP)';
elseif satAttCoord == 3
    ATT_TITLE = '(MRP)';
end

%% Plot the body-reference attitudes.
subplot(1,2,1); grid on; hold on;
ATT_BR_ARRAY = ATT_BR_ARRAY';
plot(t, ATT_BR_ARRAY(1,:), 'LineWidth', 1);
plot(t, ATT_BR_ARRAY(2,:), 'LineWidth', 1);
plot(t, ATT_BR_ARRAY(3,:), 'LineWidth', 1);

if satAttCoord == 1
    plot(t, ATT_BR_ARRAY(4,:), 'LineWidth', 1);
    legend('\beta_{BR-0}','\beta_{BR-1}','\beta_{BR-2}','\beta_{BR-3}');
else
    legend('\sigma_{BR-1}','\sigma_{BR-2}','\sigma_{BR-3}');
end
title(['Attitude Error Body-to-Reference ' ATT_TITLE]);

%% Plot the body-inertial attitudes.
subplot(1,2,2); grid on; hold on;
ATT_BN_ARRAY = ATT_BN_ARRAY';
plot(t, ATT_BN_ARRAY(1,:), 'LineWidth', 1);
plot(t, ATT_BN_ARRAY(2,:), 'LineWidth', 1);
plot(t, ATT_BN_ARRAY(3,:), 'LineWidth', 1);

if satAttCoord == 1
    plot(t, ATT_BN_ARRAY(4,:), 'LineWidth', 1);
    legend('\beta_{BN-0}','\beta_{BN-1}','\beta_{BN-2}','\beta_{BN-3}');
else
    legend('\sigma_{BN-1}','\sigma_{BN-2}','\sigma_{BN-3}');
end
title(['Attitude Error Body-to-Inertial ' ATT_TITLE]);

end

