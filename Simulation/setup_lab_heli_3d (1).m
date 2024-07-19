%% SETUP_LAB_HELI_3D
%
% This script sets the model parameters and designs a PID position
% controller using LQR for the Quanser 3-DOF Helicopter plant.
%
% Copyright (C) 2007 Quanser Consulting Inc.
% Quanser Consulting Inc.
%
clear all;
%
%% ############### USER-DEFINED 3-DOF HELI CONFIGURATION ###############
% Amplifier Gain used for yaw and pitch axes: set VoltPAQ to 3.
K_AMP = 3;
% Amplifier Maximum Output Voltage (V)
VMAX_AMP = 24;
% Digital-to-Analog Maximum Voltage (V): set to 10 for Q4/Q8 cards
VMAX_DAC = 10;
% Initial elvation angle (rad)
elev_0 = -27.5*pi/180;
%
%% ############### USER-DEFINED CONTROLLER/FILTER DESIGN ###############
% Type of Controller: set it to 'LQR_AUTO' or 'MANUAL'  
CONTROLLER_TYPE = 'LQR_AUTO';    % LQR controller design: automatic mode
%CONTROLLER_TYPE = 'MANUAL';    % controller design: manual mode
% Is Active Disturbance Mass used: 'YES' or 'NO'
% note: larger gain controller calculated for ADS system
WITH_ADS = 'YES';
% Specifications of a second-order low-pass filter
wcf = 2 * pi * 20; % filter cutting frequency
zetaf = 0.9;        % filter damping ratio
% Anti-windup: integrator saturation (V)
SAT_INT_ERR_ELEV = 7.5;
SAT_INT_ERR_TRAVEL = 7.5;
%
%% ############### USER-DEFINED DESIRED COMMAND SETTINGS ###############
% Note: These limits are imposed on both the program and joystick commands.
% Elevation position command limit (deg)
CMD_ELEV_POS_LIMIT_LOWER = elev_0*180/pi;
CMD_ELEV_POS_LIMIT_UPPER = -CMD_ELEV_POS_LIMIT_LOWER;
% Maximum Rate of Desired Position (rad/s)
CMD_RATE_LIMIT = 45.0 * pi / 180;
%
%% ############### USER-DEFINED JOYSTICK SETTINGS ###############
% Joystick input X sensitivity used for travel (deg/s/V)
K_JOYSTICK_X = 40.0;
% Joystick input Y sensitivity used for elevation (deg/s/V)
K_JOYSTICK_Y = 45.0;
%
%% ############### ACTIVE DISTURBANCE SYSTEM ###############
% Amplifier Gain for ADS: set VoltPAQ to 3.
K_AMP_X = 3;
% Maximum Output Voltage (V)
VMAX_AMP_X = 15;
% Load ADS parameters
[ K_EC_X, POS_HOME, CAL_VEL, CAL_TIME, wcf_x, kp, X_MAX, V_MAX  ] = setup_ads_configuration( );
%
%% ############### MODELING ###############
% These parameters are used for model representation and controller design.
[ Kf, m_h, m_w, m_f, m_b, Lh, La, Lw, g, K_EC_T, K_EC_P, K_EC_E ] = setup_heli_3d_configuration();
%
% For the following state vector: X = [ elevation; pitch; travel, elev_dot, pitch_dot, travel_dot]
% Initialization the state-Space representation of the open-loop System
HELI3D_ABCD_eqns;
% Augment state: Xi = [ elevation; pitch; travel, elev_dot, pitch_dot, travel_dot, elev_int, travel_int]
Ai = A;
Ai(7,1) = 1; % elevation integrator 
Ai(8,3) = 1; % travel integrator 
Ai(8,8) = 0;
Bi = B;
Bi(8,2) = 0;

Ts = 0.1;
sys_cont =  ss(A, B, C, D);
sys_disc = c2d(sys_cont, Ts, 'zoh');
Ad = sys_disc.A;
Bd = sys_disc.B;
Cd = sys_disc.C;
Dd = sys_disc.D;


% Speichere diese Matrizen im Workspace
assignin('base', 'Ad', Ad);
assignin('base', 'Bd', Bd);
assignin('base', 'Cd', Cd);
assignin('base', 'Dd', Dd);


% state and controller input constraints
uf_min = 23;
uf_max = 23;

ub_min = 23;
ub_max = 23;

elevation_min = deg2rad(45); % min constraint of x1
elevation_max = deg2rad(45); %
pitch_min = deg2rad(30);
pitch_max =  deg2rad(30);
%travel_min = -Inf; % deg2rad(-180);
%travel_max = Inf; % deg2rad(180);

%mpc controller 
kf = 15;              % prediction horizon

q1 = 100;            % penalize the elevation angle 
q2 = 200;             % penalize the pitch angle
q3 = 120;             % penalize the travel angle
q4 = 10;             % penalize the rate of change of the elevation angle
q5 = 10;             % penalize the rate of change of the pitch angle
q6 = 10;            % penalize the rate of change of the travel angle

r1 = 0.1;            % penalize the power of the front motor, 
r2 = 0.1;             % penalize the power of the back motor

Q = diag([q1 q2 q3 q4 q5 q6]); % Weighting matrix for the states
R = diag([r1 r2]);             % Weighting matrix for the inputs


% State and Input Constraints
% state constraints

wx = [
    1 0 0 0 0 0;    % constrain epsilon 
    -1 0 0 0 0 0;   % constrain -epsilon
    0 1 0 0 0 0;    % constrain rho
    0 -1 0 0 0 0;   % constrain -rho
    ];

omx =  [elevation_max;elevation_min;pitch_max;pitch_min];
% input constraints
wu = [1 0;
       -1 0;
        0 1; 
        0 -1];

omu = [uf_max;uf_min;ub_max;ub_min];
omxf = omx;
wxf = wx;

S = Q;


% Calculate the new equilibrium points
[nx, nu] = size(Bd);
[ny, ~] = size(Cd);

% Calculate the new equilibrium points
%Ad_total = [Ad - eye(nx) Bd; Cd zeros(ny, nu)];
%rhs = [zeros(nx, 1);xref];


%xe_ue = Ad_total \ rhs;
%xe = xe_ue(1:nx);
%ue = xe_ue(nx+1:end);

% shitf states and inputs
%del_omx = omx -wx*xe;
%del_omxf = omx -wx*xe;
%del_omu = omu - wu*ue;
%delta_x0 = x0 - xe;
%delta_xref = xref - Cd*xe;


% condensed matrix computation

Cneu = [C;zeros(3,6)];
       
[Qkf, Rkf] = Compute_Qkf_Rkf(kf, Q, S, R);
[Wx, Omx, Wu, Omu] = Compute_Wx_Omx_Wu_Omu(kf, wx, wxf, wu, omx, omxf, omu);
[Akf, Bkf] = Compute_Akf_Bkf(kf, A, B);
[D_cal, E_cal, Ckf] = Compute_D_cal_E_cal_Ckf(kf, B, Cneu);

%[Akf, Bkf] = Compute_Akf_Bkf(kf,Ad,Bd);
%[Qkf, Rkf] = Compute_Qkf_Rkf(kf, Q, S, R);
%[Wx, Omx, Wu, Omu] = Compute_Wx_Omx_Wu_Omu(kf,wx,wxf,wu,omx,omxf,omu);


[n,m] = size(Bd);

%Akf = [eye(n);Ad;Ad^2];
%Bkf = [zeros(n,m), zeros(n,m);Bd, zeros(n,m);Ad*Bd,Bd];
%Qkf = blkdiag(Q,Q,S);
%Rkf = blkdiag(R,R);
%Wx = blkdiag(wx,wx,wxf);
%Omx = [omx; omx; omxf];
%Wu = blkdiag(wu,wu);
%delOmu = [del_omu; del_omu];


%
%% ############### DISPLAY ###############
if strcmp ( CONTROLLER_TYPE, 'LQR_AUTO' )
    % LQR Controller Design Specifications
    if strcmp(WITH_ADS, 'NO')
        Q = diag([200 1 20 0 0 2 10 0.1]);
        R = 0.05*diag([1 1]);        
    elseif strcmp(WITH_ADS, 'YES')
        Q = diag([200 1 20 0 0 2 10 0.1]);
        R = 0.025*diag([1 1]);    
    end
    % Automatically calculate the LQR controller gain
    K = lqr( Ai, Bi, Q, R );    
    % Display the calculated gains
    disp( ' ' )
    disp( 'Calculated LQR controller gain elements: ' )
    K    
elseif strcmp ( CONTROLLER_TYPE, 'MANUAL' )
    % Automatically calculate the LQR controller gain
    K = zeros(2,6);
    % Display the calculated gains
    disp( ' ' )
    disp( 'Calculated LQR controller gain elements: ' )
    K
    disp( ' ' )
    disp( 'STATUS: manual mode' ) 
    disp( 'The model parameters of your 3 DOF Hover system have been set.' )
    disp( 'You can now design your state-feedback position controller.' )
    disp( ' ' )
else
    error( 'Error: Please set the type of controller that you wish to implement.' )
end



