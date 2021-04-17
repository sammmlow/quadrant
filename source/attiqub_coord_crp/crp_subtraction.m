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
% ### Finding the CRP attitude difference between two rotations       ###
% ###                                                                 ###
% ### By Samuel Low (15-12-2020), DSO National Laboratories           ###
% ###                                                                 ###
% #######################################################################
% #######################################################################

function [c2] = crp_subtraction(crp, c1)

% Input: crp and c1 are two arrays each holding 3 elements
% If there was a rotation DCM [FN], done by sequentially [FB]*[BN]
% Then c1 represents quarternions in the [BN] DCM
% And crp represents final quarternion in the [FN] DCM
% Thus, the inputs are defined for [FN(crp)] = [FB(c2)] * [BN(c1)]
% And so, the output of this function is to find [FB] in CRP format

c2 = ( crp - c1 + cross(crp,c1) ) / ( 1 + dot(crp,c1) );

end

