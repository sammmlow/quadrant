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
% ### Project ATTIQUB: Plotting Planetary Bodies Toolbox              ###
% ### Plots the orbit into the current figure given the ephemeris     ###
% ###                                                                 ###
% ### Adapted by Samuel Low (06-01-2021), DSO National Laboratories   ###
% ###                                                                 ###
% #######################################################################
% #######################################################################

function plot_orbit( pos_array, colorStr )

% This function takes in an ephemeris matrix of size N x 3, where N refers
% to the number of epochs observed, and 3 refers to the XYZ components.
% This function also inputs the hex colour (string) for the orbit plot.
% The figure number should be called outside of this function (in case
% multiple orbit plots are necessary, together with the central body).

pos = pos_array';

orbPlot = plot3( pos(1,:),  ...
                 pos(2,:),  ...
                 pos(3,:));

orbPlot.LineWidth = 1;
orbPlot.Color = colorStr;

end
