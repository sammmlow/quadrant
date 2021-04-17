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
% ### Conversion of a direction cosine matrix into quarternion form   ###
% ###                                                                 ###
% ### By Samuel Low (13-12-2020), DSO National Laboratories           ###
% ###                                                                 ###
% #######################################################################
% #######################################################################

function [quarternion] = dcm_conv2qrt(dcm)

% Input DCM is a 3x3 matrix.
% Output is a 1x4 quarternion vector.

% First, calculate the individual betas of the quarternion vector.
b0sq = 0.25 * ( 1 + trace(dcm) );
b1sq = 0.25 * ( 1 + 2*dcm(1,1) - trace(dcm) );
b2sq = 0.25 * ( 1 + 2*dcm(2,2) - trace(dcm) );
b3sq = 0.25 * ( 1 + 2*dcm(3,3) - trace(dcm) );

% Next, calculate the intermediate betas used in Sheppard's method.
b0b1 = ( dcm(2,3) - dcm(3,2) ) / 4;
b0b2 = ( dcm(3,1) - dcm(1,3) ) / 4;
b0b3 = ( dcm(1,2) - dcm(2,1) ) / 4;
b1b2 = ( dcm(1,2) + dcm(2,1) ) / 4;
b3b1 = ( dcm(3,1) + dcm(1,3) ) / 4;
b2b3 = ( dcm(2,3) + dcm(3,2) ) / 4;

% Check to see which beta is largest.
bmax = max( [ b0sq b1sq b2sq b3sq ] );

% Divide by the largest beta.
if bmax == b0sq
    b0 = sqrt(b0sq);
    b1 = b0b1 / b0;
    b2 = b0b2 / b0;
    b3 = b0b3 / b0;
    
elseif bmax == b1sq
    b1 = sqrt(b1sq);
    b0 = b0b1 / b1;
    b2 = b1b2 / b1;
    b3 = b3b1 / b1;
    
elseif bmax == b2sq
    b2 = sqrt(b2sq);
    b0 = b0b2 / b2;
    b1 = b1b2 / b2;
    b3 = b2b3 / b2;
    
elseif bmax == b3sq
    b3 = sqrt(b3sq);
    b0 = b0b3 / b3;
    b1 = b3b1 / b3;
    b2 = b2b3 / b3;
end

quarternion = [b0 b1 b2 b3];

end
