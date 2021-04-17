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
% ### Project ATTIQUB: Directional Cosine Matrices Toolbox            ###
% ### Conversion of a direction cosine matrix into CRP form           ###
% ###                                                                 ###
% ### By Samuel Low (15-12-2020), DSO National Laboratories           ###
% ###                                                                 ###
% #######################################################################
% #######################################################################

function [crp] = dcm_conv2crp(dcm)

% Input DCM is a 3x3 matrix.
% Output is a 1x3 CRP vector.

% Calculate zeta = 2 * beta_0 from the Euler parameters
zeta = sqrt( trace(dcm) + 1 );

% We can directly extract CRP coordinates from the DCM elements.
crp = [ (dcm(2,3)-dcm(3,2)) (dcm(3,1)-dcm(1,3)) (dcm(1,2)-dcm(2,1)) ];
crp = ( 1 / (zeta^2) ) * crp;

end

