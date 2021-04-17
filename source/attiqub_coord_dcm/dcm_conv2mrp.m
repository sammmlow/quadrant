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
% ### Conversion of a direction cosine matrix into MRP form           ###
% ###                                                                 ###
% ### By Samuel Low (15-12-2020), DSO National Laboratories           ###
% ###                                                                 ###
% #######################################################################
% #######################################################################

function [mrp] = dcm_conv2mrp(dcm)

% Input DCM is a 3x3 matrix.
% Output is a 1x3 MRP vector.

% This function will always provide the inner MRP set where |mrp| <= 1.
% Realise also that this transformation introduces a 0/0 singularity,
% at principal rotation of 180 degrees, because then zeta becomes 0.
zeta = sqrt( trace(dcm) + 1 );

% Check if a singularity occurs for 180 degree rotation.
if zeta ~= 0.0
    mt = ( transpose(dcm) - dcm ) / ( zeta * ( zeta + 2 ) );
    mrp = [ mt(3,2) mt(1,3) mt(2,1) ];

% If there was a singularity, convert it to a DCM as an intermediate step.
else
    q = dcm_conv2qrt(dcm);
    mrp = qrt_conv2mrp(q);
end

end

