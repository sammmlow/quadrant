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
% ### This definition is used when switching to the shadow set        ###
% ### Input in the MRP coordinates and the angular velocity vectors   ###
% ###                                                                 ###
% ### By Samuel Low (15-12-2020), DSO National Laboratories           ###
% ###                                                                 ###
% #######################################################################
% #######################################################################

function [mrpdot_s] = mrp_deriv_shadow(mrp, omega)

% Input: 1x3 MRP coordinates and 1x3 angular velocity vector
% Returns the 1x3 CRP dot (derivative) vector

% In this segment below, we will first express the non-shadow MRP-dot.
% Express the MRP with cross-product operator in tilde form.
mt = [ 0.0        (-1*mrp(3))   mrp(2)      ;
       mrp(3)      0.0         (-1*mrp(1))  ;
      (-1*mrp(2))  mrp(1)       0.0        ];

% Compute the B matrix used in the MRP kinematic ODE. Note, unlike NumPy, 
% there isn't an explicit method to compute the outer product. Thus, 
% we use X' * X for the outer (whereas X * X' refers to inner)
B = ( ( 1 - dot(mrp,mrp) ) * eye(3) ) + (2 * mt) + 2 * (mrp' * mrp);

% Now, we compute the actual non-shadow MRP-dot.
mrpdot = 0.25 * transpose( B * transpose(omega) );

% Finlly, we perform the shadow set conversion.
mn = dot(mrp,mrp);
ms1 = -1 * ( mrpdot ) / mn;
ms2 = ( ( ( 0.5 * ( (1 + mn) / mn^2 ) * (mrp' * mrp) ) * omega' ))' ;
mrpdot_s = ms1 + ms2;

end

