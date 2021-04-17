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
% ### Addition of two quarternions (rotations) together into one      ###
% ###                                                                 ###
% ### By Samuel Low (13-12-2020), DSO National Laboratories           ###
% ###                                                                 ###
% #######################################################################
% #######################################################################

function [quarternion] = qrt_addition(quart2, quart1)

% Input: quart1 and quart2 are two arrays each holding 4 elements
% Rotational addition in this order: output = quart2 * quart1
% If there was a rotation FN, done by sequentially [FB]*[BN]
% Then quart1 represents quarternions in the [BN] DCM
% And quart2 represents quarternions in the [FB] DCM

quart1vec = [ quart1(1) quart1(2) quart1(3) quart1(4) ];

% We will now generate the rotation matrix equivalent for quart2.
% Realise that quart2mat is orthogonal => the inverse is the transpose.

quart2mat = [ quart2(1) -1*quart2(2) -1*quart2(3) -1*quart2(4) ;
              quart2(2)    quart2(1)    quart2(4) -1*quart2(3) ;
              quart2(3) -1*quart2(4)    quart2(1)    quart2(2) ;
              quart2(4)    quart2(3) -1*quart2(2)    quart2(1) ];

% Now we can compute the quarternion [FB].
quarternion = transpose( quart2mat * transpose(quart1vec) );

end

