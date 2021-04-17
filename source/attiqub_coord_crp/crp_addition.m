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
% ### Project ATTIQUB: Classical Rodrigues Parameters Toolbox         ###
% ### Addition of two CRPs (rotations) together into one              ###
% ###                                                                 ###
% ### By Samuel Low (15-12-2020), DSO National Laboratories           ###
% ###                                                                 ###
% #######################################################################
% #######################################################################

function [crp] = crp_addition(c2, c1)

% Input: c1 and c2 are two arrays each holding 3 elements
% Rotational addition in this order: output = c2 * c1
% If there was a rotation FN, done by sequentially [FB]*[BN]
% Then c1 represents the CRP rotation as the [BN] DCM
% And c2 represents the CRP rotation as the [FB] DCM

crp = ( c2 + c1 - ( cross(c2,c1) ) ) / ( 1 - dot(c2,c1) );

end

