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
% ### Project ATTIQUB: Two Body Keplerian Toolbox                     ###
% ### Obtaining the eccentric anomaly by solving Kepler's Equation    ###
% ###                                                                 ###
% ### By Samuel Low (21-12-2020), DSO National Laboratories           ###
% ###                                                                 ###
% #######################################################################
% #######################################################################

function [E2] = twobody_kepler(M, e)

% Given some mean anomaly, M, we wish to find the eccentric anomaly E
% from the relation: "M = E - e*sin(E)", where M is input in degrees,
% and e is the unit-less eccentricity of the orbit.

% We convert to radians in order to improve the precision.
meanAnom = deg2rad(M);
E1 = meanAnom;

% Initialise convergence residual
residual = 1.0; 

% Now we solve for eccentric anomaly E2 via Newton's method.
while residual >= 0.000001
    fn = E1 - ( e * sin(E1) ) - meanAnom;
    fd = 1 - ( e * cos(E1) );
    E2 = E1 - ( fn / fd );
    residual = abs( E2 - E1 );
    E1 = E2;
end

E2 = rad2deg(E2);

end

