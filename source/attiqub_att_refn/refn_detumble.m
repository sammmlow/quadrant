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
% ### Extracts the zero attitude (DCM) and angular velocity (rad/s)   ###
% ###                                                                 ###
% ### By Samuel Low (06-01-2021), DSO National Laboratories           ###
% ###                                                                 ###
% #######################################################################
% #######################################################################

function [DCM_RN, OHM_RN] = refn_detumble()

% In the regulation problem, we are not interested in having a pointing
% reference. We simply wish to stop de-tumbling. Thus, the goal is to aim
% towards an angular velocity (in the B/N frame, observed in the N frame).

DCM_RN = eye( 3 );          % R => Reference, N => Inertial
OHM_RN = [ 0.0 0.0 0.0 ];   % R => Reference, N => Inertial

end

