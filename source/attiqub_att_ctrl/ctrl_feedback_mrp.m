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
% ### Project ATTIQUB: Attitude Control and Torque Toolbox            ###
% ### Feedback control computation and update (using MRPs)            ###
% ###                                                                 ###
% ### By Samuel Low (07-01-2020), DSO National Laboratories           ###
% ###                                                                 ###
% #######################################################################
% #######################################################################

function [ ATT_BR,   ...
           OHM_BR,   ...
           ATT_BN,   ...
           OHM_BN,   ...
           INTG_ERR, ...
           CTRL_TORQ ] = ctrl_feedback_mrp( TIME,     ... 
                                            STEPD,    ...
                                            STEPC,    ...
                                            Kp,       ...
                                            Ki,       ...
                                            Kd,       ...
                                            CTRL,     ...
                                            INERT,    ... 
                                            ATT_BN,   ... 
                                            OHM_BN,   ... 
                                            OHM_RN,   ... 
                                            DCM_RN,   ... 
                                            INTG_ERR, ...
                                            PREV_CTRL )

% ------------------------------------------------------------------------
% MEMO ON SIMULATING SPACECRAFT ATTITUDES AND EXECUTION OF CONTROL LAW
% Written by Samuel Low (07-01-2020), DSO National Laboratories.
% ------------------------------------------------------------------------
%
% The inputs of this function are:
% --> TIME   ... ... Current time in the dynamics loop
% --> STEPD  ... ... Time step in the dynamics loop
% --> STEPC  ... ... Time step in the controls software loop
% --> Kp     ... ... Feedback control proportional gain
% --> Ki     ... ... Feedback control integral gain
% --> Kd     ... ... Feedback control derivative gain
% --> CTRL   ... ... Control law choice (integer)
% --> INERT  ... ... Spacecraft 3x3 inertia matrix
% --> ATT_BN ... ... Attitude error in inertial-to-body frame [BN]
% --> OHM_BN ... ... Angular velocity vector (rad/s) inertial-to-reference
% --> OHM_RN ... ... Angular velocity vector (rad/s) inertial-to-reference
% --> DCM_RN ... ... Direction cosine matrix inertial-to-reference [RN]
% --> INTG_ERR   ... Accumulated ATT_BR error integrated over time
% --> PREV_CTRL  ... Control torque vector from previous loop (Nm)
%
% The outputs of this function are the attitude errors and angular rates:
% --> ATT_BR ... ... Reference-to-Body attitude
% --> OHM_BR ... ... Reference-to-Body angular velocity (rad/s)
% --> ATT_BN ... ... Inertial-to-Body attitude
% --> OHM_BN ... ... Inertial-to-Body angular velocity (rad/s)
% --> ATT_BN_DOT ... Inertial-to-Body attitude rate
% --> INTG_ERR   ... Reference-to-Body attitude integration error
% --> CTRL_TORQ  ... Control torque response (Nm) in current feedback
%
% ------------------------------------------------------------------------
%
% To understand the feedback control logic in this function, which is to
% be embedded in the numerical integrator loop, let us work backwards.
%
% In general for tracking, the reference-to-body outputs ATT_BR and OHM_BR
% are to be updated in the feedback. This means ATT_BR_DOT is required,
% which can be computed using ATT_BR and OHM_BR using the differential
% kinematic equation at that instant. OHM_BR can be obtained via vector
% addition (as seen in the B-frame) of OHM_BN and -OHM_RN. OHM_RN should
% be provided via the "attiqub_att_refn" libraries before calling this
% function. OHM_BN must be updated in each loop, and thus requires the
% computation of OHM_BN_DOT for the update. OHM_BN_DOT can be computed by
% providing inputs OHM_BN, the satellite inertia I, and the control law U,
% into Euler's rotational equation of motion. If OHM_RN_DOT is needed, the
% numerical first order gradient will suffice as an estimate.
%
% Take note that OHM_RN and OHM_RN_DOTs must be seen in the B-frame!
%
% Keep in mind also that if we are working with MRPs, we should switch to
% the appropriate shadow and non-shadow set values by checking if the MRP
% norm exceeds 1. This is not necessary with CRPs, since CRP attitude
% coordinates are unique (i.e., the shadow set = non-shadow set). In the
% case of quarternions, it is necessary to check that the Beta-0 term is
% greater than or equal to 0, in order to ensure that the short rotation
% is called, instead of the long rotation (think: the -359 degree case).
%
% Now, we understand what are the ingredients we need to solve for ATT_BR
% in each loop, we can update the respective derivatives, retrieve the
% appropriate reference attitudes, repeat again in the next loop.
%
% --------------------------------------------------------------------

% Initialise the reference-to-body and reference-to-inertial attitudes.
DCM_BN = mrp_conv2dcm( ATT_BN );
DCM_BR = DCM_BN * DCM_RN';
ATT_BR = dcm_conv2mrp( DCM_BR );

% Switch body-to-reference attitudes to shadow set if necessary.
if norm( ATT_BR ) >= 1
    ATT_BR = mrp_set_shadow( ATT_BR );
end

% The original OHM_RN is seen in the N-frame. Convert to B-frame.
OHM_RN = transpose(DCM_BN * OHM_RN');

% Update omega of body-to-reference (B-frame)
OHM_BR = OHM_BN - OHM_RN;

% Trigger the control response at the appropriate step.
if mod( TIME, STEPC ) == 0.0
    
    % Mode 1: Trigger the Lyapunov PD Control Response
    if CTRL == 1
        CTRL_TORQ = ctrl_response_PD( Kp, Ki, Kd, ATT_BR, OHM_BR );
        
    % Mode 2: Trigger the Lyapunov PID Control Response
    elseif CTRL == 2
        CTRL_TORQ = ctrl_response_PID( Kp, Ki, Kd, ATT_BR, OHM_BR, ...
                                       INTG_ERR );
        
    % Mode 3: Trigger the Model Predictive Control Response
    elseif CTRL == 3
        CTRL_TORQ = ctrl_response_MPC( Kp, Ki, Kd, ATT_BR, OHM_BR );
        
    % Default Mode: Lyapunov PD Control Response
    else
        CTRL_TORQ = (-Kd * ATT_BR) - ( Kp * OHM_BR ); 
    end
    
% Else, return the control torque from the previous loop.
else
    CTRL_TORQ = PREV_CTRL;
end

% Solve for OHM_BN_DOT using Euler's equation of motion.
I_OHM = transpose( INERT * OHM_BN' );
OHM_BN_DOT = ( INERT \ ( CTRL_TORQ - cross( OHM_BN, I_OHM ) )' )';

% Solve for ATT_BN_DOT using the differential kinematic equation.
if norm( ATT_BN ) >= 1
    ATT_BN = mrp_set_shadow( ATT_BN );
    ATT_BN_DOT = mrp_deriv_shadow( ATT_BN, OHM_BN );
else
    ATT_BN_DOT = mrp_derivative( ATT_BN, OHM_BN );
end

% Actual numerical integration via first order derivative update.
OHM_BN = OHM_BN + ( STEPD * OHM_BN_DOT );
ATT_BN = ATT_BN + ( STEPD * ATT_BN_DOT );
INTG_ERR = INTG_ERR + ( ATT_BR * STEPD );

end

