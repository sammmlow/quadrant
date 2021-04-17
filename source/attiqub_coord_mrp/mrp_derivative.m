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
% ### Returns the differential equation in MRP vector form            ###
% ### This definition is used when |sigma| < 1                        ###
% ### Input in the MRP coordinates and the angular velocity vectors   ###
% ###                                                                 ###
% ### By Samuel Low (15-12-2020), DSO National Laboratories           ###
% ###                                                                 ###
% #######################################################################
% #######################################################################

function [mrpdot] = mrp_derivative(mrp, omega)

% Input: 1x3 MRP coordinates and 1x3 angular velocity vector
% Returns the 1x3 CRP dot (derivative) vector

% Express the MRP with cross-product operator in tilde form.
mt = [ 0.0        (-1*mrp(3))   mrp(2)      ;
       mrp(3)      0.0         (-1*mrp(1))  ;
      (-1*mrp(2))  mrp(1)       0.0        ];

% Compute the B matrix used in the MRP kinematic ODE. Note, unlike NumPy, 
% there isn't an explicit method to compute the outer product. Thus, 
% we use X' * X for the outer (whereas X * X' refers to inner)
B = ( ( 1 - dot(mrp,mrp) ) * eye(3) ) + (2 * mt) + 2 * (mrp' * mrp);

mrpdot = 0.25 * transpose( B * transpose(omega) );

end

