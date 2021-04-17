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
% ### Calculate the 3x1 angular velocity vector, omega, given the     ###
% ### input vectors theta (3-1-3 Euler angle set) and its theta-dot,  ###
% ### which is the time derivative vector of the 3-1-3 theta.         ###
% ###                                                                 ###
% ### By Samuel Low (10-12-2020), DSO National Laboratories           ###
% ###                                                                 ###
% #######################################################################
% #######################################################################

function [omega] = get_omega_313(theta_vect, theta_vect_dot)

% Input vector is a 3-1-3 Euler angle set, and its time derivative.
%ran = theta_vect(1);
inc = theta_vect(2);
lat = theta_vect(3);

% Compute the transitional B-matrix that relates angular velocities to
% the Euler angle rate components.
B = [ (sind(lat)*sind(inc))   cosd(lat)        0.0 ;
      (cosd(lat)*sind(inc))   (-1*sind(lat))   0.0 ;
      cosd(inc)               0.0              1.0 ];

% Now, we can compute the angular velocity vector.
omega = transpose( B * transpose(theta_vect_dot) );

% Note that the Euler angle kinematic differential equations encounter a
% singularity at theta[2] = +/-90 for the 3-2-1 Euler angle set, or any
% asymmetric Euler angle set actually, and a singularity at theta[2] = 0
% or 180 for the 3-1-3 Euler angle set, or any symmetric Euler angle set.
    
end