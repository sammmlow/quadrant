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
% ### Returns the differential equation in CRP vector form            ###
% ### Input in the CRP coordinates and the angular velocity vectors   ###
% ###                                                                 ###
% ### By Samuel Low (15-12-2020), DSO National Laboratories           ###
% ###                                                                 ###
% #######################################################################
% #######################################################################

function [crpdot] = crp_derivative(crp, omega)

% Input: 1x3 CRP coordinates and 1x3 angular velocity vector
% Returns the 1x3 CRP dot (derivative) vector

crp_matrix = zeros(3,3);
crp_matrix(1,1) = 1 + crp(1)^2;
crp_matrix(1,2) = crp(1) * crp(2) - crp(3);
crp_matrix(1,3) = crp(1) * crp(3) + crp(2);
crp_matrix(2,1) = crp(2) * crp(1) + crp(3);
crp_matrix(2,2) = 1 + crp(2)^2;
crp_matrix(2,3) = crp(2) * crp(3) - crp(1);
crp_matrix(3,1) = crp(3) * crp(1) - crp(2);
crp_matrix(3,2) = crp(3) * crp(2) + crp(1);
crp_matrix(3,3) = 1 + crp(3)^2;

% Compute the derivative of the CRPs
crpdot = 0.5 * transpose( crp_matrix * transpose(omega) );
    
end

