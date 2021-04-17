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
% ### Project ATTIQUB: Modified Rodrigues Parameters Toolbox          ###
% ### Finding the MRP attitude difference between two rotations       ###
% ###                                                                 ###
% ### By Samuel Low (15-12-2020), DSO National Laboratories           ###
% ###                                                                 ###
% #######################################################################
% #######################################################################

function [m2] = mrp_subtraction(mrp, m1)

% Input: mrp and m1 are two arrays each holding 3 elements
% If there was a rotation DCM [FN], done by sequentially [FB]*[BN]
% Then m1 represents the MRP rotation in the [BN] DCM
% And mrp represents final MRP in the [FN] DCM
% Thus, the inputs are defined for [FN(mrp)] = [FB(m2)] * [BN(m1)]
% And so, the output of this function is to find [FB] in MRP format

m1x = 1 - dot(m1,m1);
mrx = 1 - dot(mrp,mrp);
m2n = ( m1x * mrp ) - ( mrx * m1 ) + 2 * ( cross(mrp,m1) );
m2 = m2n / ( 1 + ( dot(m1,m1) * dot(mrp,mrp) ) + ( 2 * dot(m1,mrp) ) );

% Warning, if the denominator reaches near zero, the MRP subtraction goes
% undefined. Thus, you should check if it goes near zero and then switch
% one of the MRPs to its shadow set.

% Another warning, at times you will find that the MRP subtraction formula
% gives you the shadow set instead of the original. Do not be alarmed if
% the output of MRP addition does not match the MRP obtained manually by
% DCM multiplication. Simply use the "mrp_set_shadow" function to verify.

end

