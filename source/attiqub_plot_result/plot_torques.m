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
% ### Plots the control feedback body torques, individual reaction    ###
% ### wheel spin torques (the momentum envelope), and individual      ###
% ### reaction wheel spin rates.                                      ###
% ###                                                                 ###
% ### By Samuel Low (13-01-2021), DSO National Laboratories           ###
% ###                                                                 ###
% #######################################################################
% #######################################################################

function plot_torques( t, RW_Num, TORQUE_ARRAY, MOTTRQ_ARRAY )

%% Inputs:
% ---> t: time axis
% ---> RW_Num: The number of reaction wheels (integer)
% ---> TORQUE_ARRAY: Matrix of body-frame control torques (t x 4).
% ---> MOTTRQ_ARRAY: Matrix of individual RW motor torques (t x RW_Num).

figure('units','normalized','outerposition',[0 0 1 1])

%% Plot the desired attitude control body torques.
subplot(1,2,1); grid on; hold on;
TORQUE_ARRAY = TORQUE_ARRAY';
plot(t, TORQUE_ARRAY(1,:), 'LineWidth', 1);
plot(t, TORQUE_ARRAY(2,:), 'LineWidth', 1);
plot(t, TORQUE_ARRAY(3,:), 'LineWidth', 1);
plot(t, TORQUE_ARRAY(4,:), 'LineWidth', 1);
title('Attitude Control Body Torque (N m)');
legend('\tau_{1}','\tau_{2}','\tau_{3}','\tau_{Norm}');

%% Plot the individual reaction wheel motor torques.
subplot(1,2,2); grid on; hold on;
MOTTRQ_ARRAY = MOTTRQ_ARRAY';
rwLegend = cell( RW_Num, 1 );
for Num = 1:RW_Num
    rwLegend{Num} = strcat( 'RW', num2str(Num) );
    plot(t, MOTTRQ_ARRAY(Num,:), 'LineWidth', 1);
end
legend(rwLegend);
title('Individual Reaction Wheel Motor Torques (N m)');

end

