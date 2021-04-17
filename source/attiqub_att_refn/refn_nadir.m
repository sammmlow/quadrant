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
% ### Extracts the nadir attitude (DCM) and angular velocity (rad/s)  ###
% ###                                                                 ###
% ### By Samuel Low (06-01-2021), DSO National Laboratories           ###
% ###                                                                 ###
% #######################################################################
% #######################################################################

function [DCM_RN, OHM_RN] = refn_nadir(pos1, vel1, trueAnom, i, R, w )

% For the nadir-pointing mode, the nadir is always fixed in definition
% with respect to the hill frame of the main spacecraft. Thus, we need
% to compute the hill frame first, and then rotate hill to nadir. The hill
% frame requires the S/C position and velocity in the inertial frame.

% Inputs: 04x Angular Keplerian elements (deg), position and velocity (m).
%         - pos1 -> Inertial frame position vector 1x3 (m)
%         - vel1 -> Inertial frame velocity vector 1x3 (m)
%         - i    -> Inclination (degrees)
%         - R    -> Right Angle of Asc Node (degrees)
%         - w    -> Argument of Perigee (degrees)
%         - v    -> True Anomaly (degrees)
% Output: Reference attitude in 3x3 DCM and reference angular velocity
% All outputs taken with respect to the inertial frame (i.e., RN frame).

% Initialise and construct the hill frame DCM.
DCM_HN = eye(3);
DCM_HN = DCM_HN * dcm_rotate_z( w + trueAnom );
DCM_HN = DCM_HN * dcm_rotate_x( i );
DCM_HN = DCM_HN * dcm_rotate_z( R );

% Nadir-pointing from hill frame is fixed by definition.
DCM_RH = [ [ -1  0  0 ] ;  ...
           [  0  1  0 ] ;  ...
           [  0  0 -1 ] ];

% Since the nadir-pointing frame is a fixed rotation away from the
% hill frame, we will obtain the inertial-to-nadir frame through
% the matrix products of the DCMs RH and HN below.
DCM_RN = DCM_RH * DCM_HN; 

% Now, we also need to compute the angular velocity of the main
% satellite, as seen in the N-frame (inertial). By the definition
% of angular velocity:
OHM_HN = cross( pos1, vel1 ) / dot( pos1, pos1 );

% Note, the angular velocity OHM_RN = OHM_HN, because the nadir
% pointing frame is actually a fixed rotation from the standard
% hill frame. If you want to prove this by DCM multiplication,
% then DCM_NH * d(trueAnom)/dt = DCM_NR * * d(trueAnom)/dt. The
% same answer would have been obtained.
OHM_RN = OHM_HN';

% Now, let us set the correct reference DCM and angular velocity.
% DCM_RN; % R => Reference, N => Inertial
% OHM_RN; % R => Reference, N => Inertial

end

