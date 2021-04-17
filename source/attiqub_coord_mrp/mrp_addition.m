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
% ### Addition of two MRPs (rotations) together into one              ###
% ###                                                                 ###
% ### By Samuel Low (15-12-2020), DSO National Laboratories           ###
% ###                                                                 ###
% #######################################################################
% #######################################################################

function [mrp] = mrp_addition(m2, m1)

% We will go through DCMs for the addition of MRPs.
% Input: m1 and m2 are two arrays each holding 3 elements
% Rotational addition in this order: output = m2 * m1
% If there was a rotation FN, done by sequentially [FB]*[BN]
% Then m1 represents the MRP rotation as the [BN] DCM
% And m2 represents the MRP rotation as the [FB] DCM

m1x = 1 - dot(m1,m1);
m2x = 1 - dot(m2,m2);
mrpi = ( m1x * m2 ) + ( m2x * m1 ) - 2 * ( cross(m2,m1) );
mrp = mrpi / ( 1 + ( dot(m1,m1) * dot(m2,m2) ) - ( 2 * dot(m1,m2) ) );

% Warning, if the denominator reaches near zero, the MRP addition goes
% undefined. Thus, you should check if it goes near zero and then switch
% one of the MRPs to its shadow set.

% Another warning, at times you will find that the MRP addition formula
% gives you the shadow set instead of the original. Do not be alarmed if
% the output of MRP addition does not match the MRP obtained manually by
% DCM multiplication. Simply use the "mrp_set_shadow" function to verify.

end

