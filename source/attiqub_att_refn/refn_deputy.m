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
% ### Project ATTIQUB: Attitude Reference Coordinates Toolbox         ###
% ### Extracts the attitude (DCM) and angular velocity (rad/s) based  ###
% ### on pointing to a deputy spacecraft in orbit.                    ###
% ###                                                                 ###
% ### By Samuel Low (06-01-2021), DSO National Laboratories           ###
% ###                                                                 ###
% #######################################################################
% #######################################################################

function [DCM_RN, OHM_RN] = refn_deputy(p1, p2, dcm_prev, k, ts)

% For the inertial to inter-satellite communications frame, denoted as
% CN, we need to first resolve the three axes of the pointing frame by
% determining the relative position vectors from the main to the
% secondary spacecraft. The antenna of the main spacecraft is assumed
% to be pointing in the negative b1 axis (body-frame).

% Inputs:
% - Position vectors (inertial) of the chief and deputy (p1 and p2)
% - DCM of deputy pointing reference [CN] from the previous loop
% - Flag to indicate first epoch or not, for DCM_CN_DOT estimation
% - Time step used in the dynamics loop

% Output: Reference attitude in 3x3 DCM and reference angular velocity
% All outputs taken with respect to the inertial frame (i.e., RN frame).

% Relative position vector, main pointing to secondary.
pos_relative = p2' - p1';

% Cross the relative position with the inertial Z-axis.
pos_cross = cross( pos_relative, [ 0 0 1 ] );

% Get the three axes of the comms pointing frame.
DCM_CN_X = -1 * pos_relative / norm( pos_relative );
DCM_CN_Y = pos_cross / norm( pos_cross );
DCM_CN_Z = cross( DCM_CN_X, DCM_CN_Y );

% Form the actual DCM describing the formation of this frame.
DCM_CN = [ DCM_CN_X ; DCM_CN_Y ; DCM_CN_Z ];

% We need to estimate the rate of change of how the main-to-
% secondary spacecraft rotation frame, relative to the inertial,
% changes with time. The only numerical way to do this is to do a
% linearised interpolation with the DCM_CN of the previous step.
if k == 1
    DCM_CN_PREV = DCM_CN;
else
    DCM_CN_PREV = dcm_prev;
end

% Update the current DCM rate.
DCM_CN_DOT = ( DCM_CN - DCM_CN_PREV ) / ts;

% Finally, we need to compute the angular velocity of relative
% pointing frame, as seen in the inertial frame.
OHM_CN_TILDE = -1 * DCM_CN' * DCM_CN_DOT;
OHM_CN = [ OHM_CN_TILDE(3,2) OHM_CN_TILDE(1,3) OHM_CN_TILDE(2,1) ];

% Now, let us set the correct reference DCM and angular velocity.
DCM_RN = DCM_CN; % R => Reference, N => Inertial
OHM_RN = OHM_CN; % R => Reference, N => Inertial

end

