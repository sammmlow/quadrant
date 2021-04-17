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
% ### Project ATTIQUB: Quarternions Toolbox                           ###
% ### Returns the differential equation in quarternion vector form    ###
% ### Input in the quarternion and the angular velocity vectors       ###
% ###                                                                 ###
% ### By Samuel Low (13-12-2020), DSO National Laboratories           ###
% ###                                                                 ###
% #######################################################################
% #######################################################################

function [quartdot] = qrt_derivative(quart, omega)

% Input: 1x4 quarternion coordinates and 1x3 angular velocity vector
% Returns the 1x4 quarternion dot (derivative) vector

quartmat = [ -1*quart(2) -1*quart(3) -1*quart(4) ;
                quart(1) -1*quart(4)    quart(3) ;
                quart(4)    quart(1) -1*quart(2) ;
             -1*quart(3)    quart(2)    quart(1) ];

% Note that the transformation matrix quartmat that relates d(beta)/dt
% and the angular velocity vector is orthogonal and singularity free.
% This means the inverse transformation can also be defined.

quartdot = 0.5 * transpose( quartmat * transpose(omega) );

end

