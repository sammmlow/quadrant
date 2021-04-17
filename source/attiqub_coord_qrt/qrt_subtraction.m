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
% ### Finding the quarternion difference between two rotations        ###
% ###                                                                 ###
% ### By Samuel Low (13-12-2020), DSO National Laboratories           ###
% ###                                                                 ###
% #######################################################################
% #######################################################################

function [quarternion] = qrt_subtraction(quart2, quart1)

% Input: quart1 and quart2 are two arrays each holding 4 elements
% If there was a rotation DCM [FN], done by sequentially [FB]*[BN]
% Then quart1 represents quarternions in the [BN] DCM
% And quart2 represents final quarternion in the [FN] DCM
% Thus, the inputs are defined for [FN(q2)] = [FB(quarternion)] * [BN(q1)]
% We would need to find the 4x4 matrix in the middle, [FB]
% We start by inverting the quarternion from [BN] to [NB]

quart1vec = [ quart1(1) -1*quart1(2) -1*quart1(3) -1*quart1(4) ];

% Realise that quart2mat is orthogonal => inverse is the transpose
% We find the 4x4 matrix for [FN]
    
quart2mat = [ quart2(1) -1*quart2(2) -1*quart2(3) -1*quart2(4) ;
              quart2(2)    quart2(1)    quart2(4) -1*quart2(3) ;
              quart2(3) -1*quart2(4)    quart2(1)    quart2(2) ;
              quart2(4)    quart2(3) -1*quart2(2)    quart2(1) ];
    
% Now we have [FB] = [FN]*[NB]
quarternion = transpose( quart2mat * transpose(quart1vec) );

end

