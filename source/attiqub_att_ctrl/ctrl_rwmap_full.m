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
% ### Project ATTIQUB: Attitude Control and Torque Toolbox            ###
% ### Maps the desired body attitude control torque into individual   ###
% ### reaction wheel spin torques using a minimum norm solution       ###
% ###                                                                 ###
% ### Adapted from J Martin & H Schaub,                               ###
% ### University of Colorado Boulder, AVS Lab                         ###
% ### By Samuel Low (21-12-2020), DSO National Laboratories           ###
% ###                                                                 ###
% #######################################################################
% #######################################################################

function [motorTorqueArray] = ctrl_rwmap_full(Gsx, Trq)

% By Euler's Rotational Equations, [RW_Gsx] * RW_Trq' = -L,
% Where L is the desired attitude control torque, as seen in the B-frame.
% Num = Total number of reaction wheels (Num)
% MxT = Maximum spin torque per RW (Nm)
% MxR = Maximum spin rate per RW (rad/s)
% Ohm = Initial spin rates on individual reaction wheels (1 x Num)
% Gsx = Projection matrix of spin axes vectors per RW [g1...gN] (3 x Num)
% Trq = Computed desired attitude control body torque (Nm)

motorTorqueArray = Gsx' * ( ( Gsx * Gsx' ) \ ( -1 * Trq' ) );

% TO-DO: Create a reduced motor torque array.

end

