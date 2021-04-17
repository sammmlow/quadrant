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
% ### Extracts the sun attitude (DCM) and angular velocity (rad/s)    ###
% ###                                                                 ###
% ### By Samuel Low (06-01-2021), DSO National Laboratories           ###
% ###                                                                 ###
% #######################################################################
% #######################################################################

function [DCM_RN, OHM_RN] = refn_sunpt()

% For the sun-pointing mode, we have defined the sun to always be
% located when the main spacecraft points in the n2 unit direction
% vector (and where n3 is the direction vector pointing out of the
% ecliptic, and n1 x n2 = n3). The DCM_SN is defined below.

% Inputs: None. Since the rotation from inertial to sun-pointing is a
% fixed rotation (not time-varying), external inputs are not needed.

% Output: Reference attitude in 3x3 DCM and reference angular velocity
% All outputs taken with respect to the inertial frame (i.e., RN frame).

% Sun-pointing from inertial frame is fixed by definition.
DCM_SN = [ [ -1  0  0 ] ;  ...
           [  0  0  1 ] ;  ...
           [  0  1  0 ] ];

% Now, let us set the correct reference DCM and angular velocity.
        DCM_RN = DCM_SN;          % R => Reference, N => Inertial
        OHM_RN = [ 0.0 0.0 0.0 ]; % R => Reference, N => Inertial

end

