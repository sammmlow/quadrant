% #######################################################################
% #######################################################################
% ###                                                                 ###
% ###  @@@@@@@@ @@@@@@@ @@@@@@@ @@@@@ @@@@@@@@@ @@@    @@@ @@@@@@@@   ###
% ###  @@@@@@@@ @@@@@@@ @@@@@@@ @@@@@ @@    (@@ @@@    @@@ @@@  @@@   ###
% ###  @@    @@   @@@     @@@    @@@  @@    (@@ @@@    @@@ @@@@@@@@@  ###
% ###  @@@@@@@@   @@@     @@@    @@@  @@@@@@@@@ @@@    @@@ @@@@@@@@@  ###
% ###  @@    @@   @@@     @@@   @@@@@ @@@@@@@@@ @@@@@@@@@@ @@@   @@@  ###
% ###  @@    @@   @@@     @@@   @@@@@     (@@@  @@@@@@@@@@ @@@@@@@@@  ###
% ###                                                                 ###
% ###  Version 0.5 (Alpha)                                            ###
% ###  By Samuel Low, DSO National Laboratories                       ###
% ###  First Build: 10-12-2020                                        ###
% ###  Last Modified: 13-01-2021                                      ###
% ###                                                                 ###
% #######################################################################
% #######################################################################

function attiqub_engine()

clc; clear; close all;
load('parameters.mat');

% RW_Gsx is defined as RW_G = [ g1 g2 g3 ... gm ... gN ] where [gm] refers
% to the spin axis unit direction vector, as seen in the B-frame.
% TO-DO: Account for saturated RWs based on max torque and max speed.
% TO-DO: Perform reduced torque mapping [CB] matrix
% TO-DO: Enable/disable option for secondary spacecraft
% TO-DO: Allow for computation of the reaction wheel speeds.
% TO-DO: enable a mode for automatic feedback gain selection.
% TO-DO: For ISC antenna, make it adjustable to any axis?
% TO-DO: Should we make the ISC antenna location customisable?

%% #######################################################################
% ########################################################################
% ###                                                                  ###
% ###    FILE PATH RESOLUTION IN ORDER TO ACCESS ATTIQUB LIBRARIES     ###
% ###                                                                  ###
% ########################################################################
% ########################################################################

% ATTIQUB must resolve all file paths in order to access libraries.
attqubPath = mfilename('fullpath');
[directory, ~, ~]  = fileparts(attqubPath);
srcPaths = {[ directory '\graphics\'                   ];  ...
            [ directory '\source\'                     ];  ...
            [ directory '\source\attiqub_att_ctrl\'    ];  ...
            [ directory '\source\attiqub_att_refn\'    ];  ...
            [ directory '\source\attiqub_coord_dcm\'   ];  ...
            [ directory '\source\attiqub_coord_qrt\'   ];  ...
            [ directory '\source\attiqub_coord_crp\'   ];  ...
            [ directory '\source\attiqub_coord_mrp\'   ];  ...
            [ directory '\source\attiqub_orbit_prop\'  ];  ...
            [ directory '\source\attiqub_plot_body\'   ];  ...
            [ directory '\source\attiqub_plot_result\' ]};
        
% Add all the paths in the srcPaths cell array.
for n = 1 : length(srcPaths)
    path = string(srcPaths(n));
    addpath(path);
end

%% #######################################################################
% ########################################################################
% ###                                                                  ###
% ###    THE INITIALISATION OF ATTIQUB BEGINS AFTER THIS CODE BLOCK.   ###
% ###                                                                  ###
% ########################################################################
% ########################################################################

% Initialize the start time as zero.
tStart = 0;

% Initialize the gravitational constant based on the primary body.
gravBodies = [ 3.9860e+14 ... % 1 = Earth
               4.9041e+12 ... % 2 = Moon
               4.2828e+13 ... % 3 = Mars
               3.2487e+14 ... % 4 = Venus
               2.1925e+13 ];  % 5 = Mercury

% Initialize the radius of the primary body.
radBodies = [ 6371140 ... % 1 = Earth
              1737100 ... % 2 = Moon
              3389500 ... % 3 = Mars
              6051800 ... % 4 = Venus
              2439700 ];  % 5 = Mercury

% Select the correct gravitational constant.
mu = gravBodies( body );
rad = radBodies( body );

% Compute the scenario's time axis array.
t = linspace( tStart, tFinal, ( (1 / tStepD) * ( tFinal - tStart ) + 1 ) );

% Get the length of the number of epochs.
tN = length(t);

% Compute the Keplerian mean motion of the main and secondary spacecraft
meanMotion1 = sqrt( mu / ( orb1A^3 ) );
meanMotion2 = sqrt( mu / ( orb2A^3 ) );

% Initialize a position and vector array across time.
pos1_array = zeros(tN,3); % Position of main spacecraft
vel1_array = zeros(tN,3); % Velocity of main spacecraft
pos2_array = zeros(tN,3); % Position of secondary spacecraft
vel2_array = zeros(tN,3); % Velocity of secondary spacecraft

% Initialize the array of attitude coordinates for plotting.
if satAttCoord == 1
    ATT_CRD_NUM = 4;
elseif satAttCoord == 2
    ATT_CRD_NUM = 3;
elseif satAttCoord == 3
    ATT_CRD_NUM = 3;
else
    fprintf('Invalid attitude coordinate chosen! \n');
end

% Initialize the array of attitude coordinates for plotting. 
ATT_BR_ARRAY = zeros( tN, ATT_CRD_NUM ); % Body-to-reference
ATT_BN_ARRAY = zeros( tN, ATT_CRD_NUM ); % Body-to-inertial

% Initialize the array of final angular velocity vectors for plotting.
OHM_BR_ARRAY = zeros(tN,4);
OHM_BN_ARRAY = zeros(tN,4);

% Initialize the array of body torque vectors for plotting.
TORQUE_ARRAY = zeros(tN,4);

% Initialize the individual reaction wheel motor torques for plotting.
MOTTRQ_ARRAY = zeros(tN,RW_Num);

% Initialize the reference-to-inertial DCM.
DCM_RN = eye(3);

% Initialize the attitude integration error.
INTG_ERR = zeros(1,ATT_CRD_NUM);

% Initialize the control torque.
CTRL_TORQ = zeros(1,3);

% Check if the chief or deputy orbits below the planet surface.
if ( 1 - orb1E ) * orb1A <= rad
    fprintf('Warning! \n');
    fprintf('Chief spacecraft periapsis below planet surface! \n');
end
if ( 1 - orb2E ) * orb2A <= rad
    fprintf('Warning! \n');
    fprintf('Deputy spacecraft periapsis below planet surface! \n');
end

%% #######################################################################
% ########################################################################
% ###                                                                  ###
% ###   THE ATTQUB ENGINE BEGINS LOOPING AFTER THIS CODE BLOCK, FOR    ###
% ###   THE NUMERICAL INTEGRATION OF ATTITUDES AND ORBIT PROPAGATION.  ###
% ###   TO THE USER, PLEASE DO NOT CHANGE ANYTHING BELOW.              ###
% ###                                                                  ###
% ########################################################################
% ########################################################################

for k = 1 : tN
    
    % Get the current time in seconds.
    time = t(k);
    
    %% ###################################################################
    % ####################################################################
    % ###                                                              ###
    % ###         ORBIT PROPAGATION BLOCK IN THE DYNAMICS LOOP         ###
    % ###                                                              ###
    % ####################################################################
    % ####################################################################
    
    % Mean anomaly 1: Computation of Chief's M (-180 : +180 deg).
    orb1M0 = mod( orb1M + 180 + rad2deg( meanMotion1 * time), 360 ) - 180;
    
    % Orbit propagation 1: Chief's position, velocity and true anomaly.
    [ p1, v1, nn1 ] = twobody_posvel( orb1A, orb1E, orb1I, ...
                                      orb1R, orb1W, orb1M0, mu );
    
    % Mean anomaly 2: Computation of Deputy's M (-180 : +180 deg).
    orb2M0 = mod( orb2M + 180 + rad2deg( meanMotion2 * time), 360 ) - 180;
    
    % Orbit propagation 2: Deputy's position, velocity and true anomaly.
    [ p2, v2, nn2 ] = twobody_posvel( orb2A, orb2E, orb2I, ...
                                      orb2R, orb2W, orb2M0, mu );
    
    % Store the position-velocities of the main and deputy spacecraft.
    for m = [ 1 2 3 ] % [ 1 2 3 ] ==> [ X Y Z ]
        pos1_array(k,m) = p1(m);
        pos2_array(k,m) = p2(m);
        vel1_array(k,m) = v1(m);
        vel2_array(k,m) = v2(m);
    end
    
    %% ###################################################################
    % ####################################################################
    % ###                                                              ###
    % ###  SELECTION OF REFERENCE ATTITUDE AND ANGULAR VELOCITY RATES  ###
    % ###                                                              ###
    % ####################################################################
    % ####################################################################
    
    % Satellite Pointing Mode 1: Nadir Pointing
    if satPointing == 1
        [DCM_RN, OHM_RN] = refn_nadir( p1, v1, nn1, orb1I, orb1R, orb1W );
    
    % Satellite Pointing Mode 2: Sun Pointing
    elseif satPointing == 2
        [DCM_RN, OHM_RN] = refn_sunpt();
    
    % Satellite Pointing Mode 3: Inter-satellite Pointing (Deputy)
    elseif satPointing == 3
        [DCM_RN, OHM_RN] = refn_deputy( p1, p2, DCM_RN, k, tStepD );
        
    % Satellite Pointing Mode 4: De-tumbling or regulation (default)
    else
        [DCM_RN, OHM_RN] = refn_detumble();
    
    end % End of satellite pointing reference attitude computations.
    
    %% ###################################################################
    % ####################################################################
    % ###                                                              ###
    % ###    COMPUTATION OF BODY TO REFERENCE / INERTIAL ATTITUDES,    ###
    % ###    ANGULAR RATES, AND APPLICATION OF ATTITUDE CONTROL LAW    ###
    % ###                                                              ###
    % ####################################################################
    % ####################################################################
    
    % If attitude coordinates are in quarternions ...
    if satAttCoord == 1
        
        [ ATT_BR,   ...
          OHM_BR,   ...
          ATT_BN,   ...
          OHM_BN,   ...
          INTG_ERR, ...
          CTRL_TORQ ] = ctrl_feedback_qrt( time,       ... % Time
                                           tStepD,     ... % Dynamics step
                                           tStepC,     ... % Controls step
                                           Kp,         ... % P-Gain
                                           Ki,         ... % I-Gain
                                           Kd,         ... % D-Gain
                                           ctrlLaw,    ... % Control law
                                           satInertia, ... % Inertia matrix
                                           ATT_BN,     ... % Attitude [BN]
                                           OHM_BN,     ... % Omega [BN]
                                           OHM_RN,     ... % Omega [RN]
                                           DCM_RN,     ... % DCM [RN]
                                           INTG_ERR,   ... % Integral Error
                                           CTRL_TORQ );    % Control torque
    
    % If attitude coordinates are in CRPs ...
    elseif satAttCoord == 2
        
        [ ATT_BR,   ...
          OHM_BR,   ...
          ATT_BN,   ...
          OHM_BN,   ...
          INTG_ERR, ...
          CTRL_TORQ ] = ctrl_feedback_crp( time,       ... % Time
                                           tStepD,     ... % Dynamics step
                                           tStepC,     ... % Controls step
                                           Kp,         ... % P-Gain
                                           Ki,         ... % I-Gain
                                           Kd,         ... % D-Gain
                                           ctrlLaw,    ... % Control law
                                           satInertia, ... % Inertia matrix
                                           ATT_BN,     ... % Attitude [BN]
                                           OHM_BN,     ... % Omega [BN]
                                           OHM_RN,     ... % Omega [RN]
                                           DCM_RN,     ... % DCM [RN]
                                           INTG_ERR,   ... % Integral Error
                                           CTRL_TORQ );    % Control torque
    
    % ####################################################################
    % ####################################################################
    
    % If attitude coordinates are in MRPs.
    elseif satAttCoord == 3 
        
        [ ATT_BR,   ...
          OHM_BR,   ...
          ATT_BN,   ...
          OHM_BN,   ...
          INTG_ERR, ...
          CTRL_TORQ ] = ctrl_feedback_mrp( time,       ... % Time
                                           tStepD,     ... % Dynamics step
                                           tStepC,     ... % Controls step
                                           Kp,         ... % P-Gain
                                           Ki,         ... % I-Gain
                                           Kd,         ... % D-Gain
                                           ctrlLaw,    ... % Control law
                                           satInertia, ... % Inertia matrix
                                           ATT_BN,     ... % Attitude [BN]
                                           OHM_BN,     ... % Omega [BN]
                                           OHM_RN,     ... % Omega [RN]
                                           DCM_RN,     ... % DCM [RN]
                                           INTG_ERR,   ... % Integral Error
                                           CTRL_TORQ );    % Control torque
    end
    
    %% ###################################################################
    % ####################################################################
    % ###                                                              ###
    % ###   APPEND COMPUTED ATTITUDES WITH RESPECT TO THE REFERENCE    ###
    % ###  AND INERTIAL, THE ANGULAR VELOCITIES (RAD/S), THE CONTROL   ###
    % ###  BODY TORQUES, AND TORQUES MAPPED ONTO THE REACTION WHEELS   ###
    % ###                                                              ###
    % ####################################################################
    % ####################################################################
    
    % Map the computed control torques onto individual reaction wheels.
    if mod( time, tStepC ) == 0.0
        CTRL_TORQ_RW = ctrl_rwmap_full( RW_Gsx, CTRL_TORQ );
    end
    
    % Append all attitude errors, angular velocities, and torques.
    for m = 1 : ATT_CRD_NUM
        ATT_BR_ARRAY(k,m) = ATT_BR(m); % Attitude Error [BR]
        ATT_BN_ARRAY(k,m) = ATT_BN(m); % Attitude Error [BN]
        if m < 4
            OHM_BR_ARRAY(k,m) = OHM_BR(m); % Angular Velocity Error [BR]
            OHM_BN_ARRAY(k,m) = OHM_BN(m); % Angular Velocity Error [BN]
            TORQUE_ARRAY(k,m) = CTRL_TORQ(m); % Control Torque (Nm)
        end
    end
    
    % Append norms of attitude errors, angular velocities, and torques.
    OHM_BR_ARRAY(k,4) = norm( OHM_BR );
    OHM_BN_ARRAY(k,4) = norm( OHM_BN );
    TORQUE_ARRAY(k,4) = norm( CTRL_TORQ );
    
    % Append all individual RW motor torques.
    for n = 1 : RW_Num
        MOTTRQ_ARRAY(k,n) = CTRL_TORQ_RW(n);
    end
    
    %% ###################################################################
    % ####################################################################
    % ###                                                              ###
    % ###  SANITY CHECK, IN CASE SPACECRAFT IS SPINNING OUT OF CONTROL ###
    % ###                                                              ###
    % ####################################################################
    % ####################################################################
    
    if ( OHM_BR * OHM_BR' ) > 6.2832 % Two pi radians per second
        fprintf('Spacecraft spinning faster than 1 RPS! \n');
        fprintf('Tune down your control gains! \n');
        fprintf('Ending program at t = %.2f seconds. \n', time);
        break
    end
end

%% #######################################################################
% ########################################################################
% ###                                                                  ###
% ###   PLOTTING OF ORBITS, ATTITUDE ERRORS, ANGLE RATES AND TORQUES   ###
% ###                                                                  ###
% ########################################################################
% ########################################################################

% Plot the control feedback torque, and the mapped RW spin torques.
plot_torques( t, RW_Num, TORQUE_ARRAY, MOTTRQ_ARRAY )

% Plot the body-reference and body-inertial angular velocities.
plot_omegas( t, OHM_BR_ARRAY, OHM_BN_ARRAY );

% Plot the body-reference and body-inertial attitudes.
plot_sigmas( satAttCoord, t, ATT_BR_ARRAY, ATT_BN_ARRAY );

% Plot the central body.
plot_body( body );

% Plot the orbits and set the orbit colors (in hex).
% Orbits share the same figure object as central body.
plot_orbit( pos1_array, '#809bff' );
plot_orbit( pos2_array, '#52ff9a' );

end
