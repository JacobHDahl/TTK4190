% Project in TTK4190 Guidance and Control of Vehicles 
%
% Author:           Jacubi Dahl
% Study program:    Kyb :|

clear;
clc;
addpath(genpath("C:\Users\jacob\Documents\Student\Fart�y\TTK4190\MSS"));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h  = 0.1;    % sampling time [s]
Ns = 10000;  % no. of samples

psi_ref = 10 * pi/180;  % desired yaw angle (rad)
U_d = 9;                % desired cruise speed (m/s)
               
% ship parameters 
m = 17.0677e6;          % mass (kg)
Iz = 2.1732e10;         % yaw moment of inertia about CO (kg m^3)
xg = -3.7;              % CG x-ccordinate (m)
L = 161;                % length (m)
B = 21.8;               % beam (m)
T = 8.9;                % draft (m)
%KT = 0.7;               % propeller coefficient (-)
Dia = 3.3;              % propeller diameter (m)
rho = 1025;             % density of water (kg/m^3)
visc = 1e-6;            % kinematic viscousity at 20 degrees (m/s^2)
eps = 0.001;            % a small number added to ensure that the denominator of Cf is well defined at u=0
k = 0.1;                % form factor giving a viscous correction
t_thr = 0.05;           % thrust deduction number

% rudder limitations
delta_max  = 40 * pi/180;        % max rudder angle      (rad)
Ddelta_max = 5  * pi/180;        % max rudder derivative (rad/s)

% added mass matrix about CO
Xudot = -8.9830e5;
Yvdot = -5.1996e6;
Yrdot =  9.3677e5;
Nvdot =  Yrdot;
Nrdot = -2.4283e10;
MA = -[ Xudot 0    0 
        0 Yvdot Yrdot
        0 Nvdot Nrdot ];

% rigid-body mass matrix
MRB = [ m 0    0 
        0 m    m*xg
        0 m*xg Iz ];
    
Minv = inv(MRB + MA); % Added mass is included to give the total inertia

% ocean current in NED
beta_Vc = deg2rad(45);  % current direction (rad)
Vc = 1;              

% wind expressed in NED
Vw = 10;                   % wind speed (m/s)
betaVw = deg2rad(135);     % wind direction (rad)
rho_a = 1.247;             % air density at 10 deg celsius
cy = 0.95;                 % wind coefficient in sway
cn = 0.15;                 % wind coefficient in yaw
A_Lw = 10 * L;             % projected lateral area

% linear damping matrix (only valid for zero speed)
T1 = 20; T2 = 20; T6 = 10;

Xu = -(m - Xudot) / T1;
Yv = -(m - Yvdot) / T2;
Nr = -(Iz - Nrdot)/ T6;
D = diag([-Xu -Yv -Nr]);         % zero speed linear damping

% rudder coefficients (Section 9.5)
b = 2;
AR = 8;
CB = 0.8;

lambda = b^2 / AR;
tR = 0.45 - 0.28*CB;
CN = 6.13*lambda / (lambda + 2.25);
aH = 0.75;
xH = -0.4 * L;
xR = -0.5 * L;

X_delta2 = 0.5 * (1 - tR) * rho * AR * CN;
Y_delta = 0.25 * (1 + aH) * rho * AR * CN; 
N_delta = 0.25 * (xR + aH*xH) * rho * AR * CN;   

% input matrix
Bu = @(u_r,delta) [ (1-t_thr)  -u_r^2 * X_delta2 * delta
                        0      -u_r^2 * Y_delta
                        0      -u_r^2 * N_delta            ];
                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
% Heading Controller
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rudder control law
wb = 0.06;
zeta = 1;
wn = 1 / sqrt( 1 - 2*zeta^2 + sqrt( 4*zeta^4 - 4*zeta^2 + 2) ) * wb;
T_reg = 175;
K_reg = 1/130;
m_reg = T_reg/K_reg;

%Controller gains:
Kp = m_reg * wn^2;
Kd = (2*zeta*wn*T_reg - 1)/K_reg;
Ki = wn*Kp/10;

%Kd seems very large but what I got

% linearized sway-yaw model (see (7.15)-(7.19) in Fossen (2021)) used
% for controller design. The code below should be modified.
C_RB_lin = [0,0,0;
            0,0,m*U_d;
            0,0,m*xg*U_d];
C_A_lin = [0,0,0;
            0,0,-Xudot*U_d;
            0,(Xudot-Yvdot)*U_d,-Yrdot*U_d];
N_lin = C_RB_lin + C_A_lin + D;
M_lin = MRB + MA;
b_lin = [-2*U_d*Y_delta -2*U_d*N_delta]';

%mystuff
syms s;
H = s*M_lin(2:3,2:3) + N_lin(2:3,2:3);
Hinv = inv(H);
TF_delta_r = (Hinv(2,1)*b_lin(1) + Hinv(2,2)*b_lin(2)); %didnt work

% nu_dot =  - Minv * N_lin * nu + Minv * b_lin * delta
% x_dot = Ax + Bu 
% A = -Minv * N_lin 
% B = Minv*b_lin
%C = [0,1];
%D = 0;

%[b,a] = ss2tf(-Minv(2:3,2:3)*N_lin(2:3,2:3),Minv(2:3,2:3)*b_lin,C,D);




% initial states
eta = [0 0 0]';
nu  = [0.1 0 0]';
delta = 0;
n = 0;
Q_m = 0;

%reference generation
xd = [0,0,0]';
w_ref = 0.03;
zeta_ref = 1;
wn_ref  = w_ref;%1 / sqrt( 1 - 2*zeta_ref^2 + sqrt( 4*zeta_ref^4 - 4*zeta_ref^2 + 2) ) * w_ref;


Ad = [ 0 1 0;
       0 0 1;
      -wn_ref^3,  -(2*zeta_ref +1)*wn_ref^2,  -(2*zeta_ref+1)*wn_ref ];
  
Bd = [0, 0, wn_ref^3]';

int_error=0;
e=0;

Ja = 0;
PD = 1.5;
AEAO = 0.65;
blades = 4;
[KT, KQ] = wageningen(Ja,PD,AEAO,blades);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simdata = zeros(Ns+1,14);                % table of simulation data
nu_r_Data = zeros(Ns+1,3);
ref_data = zeros(Ns+1,3);
for i=1:Ns+1

    if(i/10 > 500)
        psi_ref = -20*pi/180;
    else
        psi_ref = 10*pi/180;
    end
    
    
    
    t = (i-1) * h;                      % time (s)
    R = Rzyx(0,0,eta(3));
    
    % current (should be added here)
    v_c_n = [Vc*cos(beta_Vc), Vc*sin(beta_Vc), 0]';
    v_c = R' * v_c_n;
    nu_r = nu - v_c; %comment out -vc to have no current
    nu_r_Data(i,:) = nu_r;
    
    % wind (should be added here)
    if t > 200
        gamma_w = eta(3) - 3.14 - betaVw;
        
        C_Y = cy*sin(gamma_w);
        C_N = cn*sin(2*gamma_w);
        Ywind = 0.5*rho_a*Vw^2 * C_Y * A_Lw;
        Nwind = 0.5*rho_a*Vw^2 * C_N*A_Lw*L;
    else
        Ywind = 0;
        Nwind = 0;
    end
    %tau_env = [0 Ywind Nwind]';
    tau_env = [0 0 0]'; %switch to add wind
    
    % state-dependent time-varying matrices
    CRB = m * nu(3) * [ 0 -1 -xg 
                        1  0  0 
                        xg 0  0  ];
                    
    % coriolis due to added mass
    CA = [  0   0   Yvdot * nu_r(2) + Yrdot * nu_r(3)
            0   0   -Xudot * nu_r(1) 
          -Yvdot * nu_r(2) - Yrdot * nu_r(3)    Xudot * nu_r(1)   0];
    N = CRB + CA + D;
    
    % nonlinear surge damping
    Rn = L/visc * abs(nu_r(1));
    Cf = 0.075 / ( (log(Rn) - 2)^2 + eps);
    Xns = -0.5 * rho * (B*L) * (1 + k) * Cf * abs(nu_r(1)) * nu_r(1);
    
    % cross-flow drag
    Ycf = 0;
    Ncf = 0;
    dx = L/10;
    Cd_2D = Hoerner(B,T);
    for xL = -L/2:dx:L/2
        vr = nu_r(2);
        r = nu_r(3);
        Ucf = abs(vr + xL * r) * (vr + xL * r);
        Ycf = Ycf - 0.5 * rho * T * Cd_2D * Ucf * dx;
        Ncf = Ncf - 0.5 * rho * T * Cd_2D * xL * Ucf * dx;
    end
    d = -[Xns Ycf Ncf]';
    
    % reference models
    xd_dot = Ad * xd + Bd * psi_ref; 
    psi_d = xd(1);
    r_d = xd(2);
    u_d = U_d;
    
    
    
    
    % thrust 
    
    
    thr = rho * Dia^4 * KT * abs(n) * n;    % thrust command (N)
    torq = rho*Dia^5 * KQ * abs(n)*n;
    
    
        
    % control law
    e_psi   = ssa(eta(3)-psi_d);          % Difference between desired and actual psi
    e_r     = nu(3)-r_d;
    delta_c = -(Kp*e_psi + Ki*int_error + Kd*e_r);
    
    
    % ship dynamics
    u = [ thr delta ]';
    tau = Bu(nu_r(1),delta) * u;
    nu_dot = Minv * (tau_env + tau - N * nu_r - d); 
    eta_dot = R * nu;    
    
    % Rudder saturation and dynamics (Sections 9.5.2)
    if abs(delta_c) >= delta_max
        delta_c = sign(delta_c)*delta_max;
    end
    
    delta_dot = delta_c - delta;
    if abs(delta_dot) >= Ddelta_max
        delta_dot = sign(delta_dot)*Ddelta_max;
    end    
    
    % propeller dynamics
    t_dumb = 0.05; %change?
    T_d  = ((U_d) * Xu)/(t_dumb -1);
    abs_n_n = T_d / (rho * Dia^4 * KT);
    n_d = sign(abs_n_n)*sqrt(abs(abs_n_n));
    
    Q_d = rho * Dia^5 * KQ * abs(n_d)*n_d;  
    Im = 100000; Tm = 10; Km = 0.6;         % propulsion parameters
    Y = Q_d * (1/Km);
    [Am, Bm, Cm, Dm] = tf2ss(Km, [Tm, 1]);
    
    e_dot = Am*e + Bm*Y;
    Q_m = Cm*e + Dm*Y;
    
    Q_f = 0;    %?????
    Q = torq;
    n_c = 10;                               % propeller speed (rps)
    n_dot = (Q_m - Q - Q_f)/(Im);           
    
    % store simulation data in a table (for testing)
    simdata(i,:) = [t n_c delta_c n delta eta' nu' u_d psi_d r_d];       
     
    % Euler integration
    eta = euler2(eta_dot,eta,h);
    nu  = euler2(nu_dot,nu,h);
    delta = euler2(delta_dot,delta,h);   
    n  = euler2(n_dot,n,h);    
    
    %my stuff
    xd = euler2(xd_dot,xd,h);
    int_error = euler2(e_psi,int_error,h);
    ref_data(i,:) = xd';
    n_dot = euler2(n_dot,n,h);
    e_dot = euler2(e_dot,e,h);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t       = simdata(:,1);                 % s
n_c     = 60 * simdata(:,2);            % rpm
delta_c = (180/pi) * simdata(:,3);      % deg
n       = 60 * simdata(:,4);            % rpm
delta   = (180/pi) * simdata(:,5);      % deg
x       = simdata(:,6);                 % m
y       = simdata(:,7);                 % m
psi     = (180/pi) * simdata(:,8);      % deg
u       = simdata(:,9);                 % m/s
v       = simdata(:,10);                % m/s
r       = (180/pi) * simdata(:,11);     % deg/s
u_d     = simdata(:,12);                % m/s
psi_d   = (180/pi) * simdata(:,13);     % deg
r_d     = (180/pi) * simdata(:,14);     % deg/s

figure(1)
figure(gcf)
subplot(311)
plot(y,x,'linewidth',2); axis('equal')
title('North-East positions (m)'); xlabel('time (s)'); 
subplot(312)
plot(t,psi,t,psi_d,'linewidth',2);
title('Actual and desired yaw angles (deg)'); xlabel('time (s)');
subplot(313)
plot(t,r,t,r_d,'linewidth',2);
title('Actual and desired yaw rates (deg/s)'); xlabel('time (s)');

figure(2)
figure(gcf)
subplot(311)
plot(t,u,t,u_d,'linewidth',2);
title('Actual and desired surge velocities (m/s)'); xlabel('time (s)');
subplot(312)
plot(t,n,t,n_c,'linewidth',2);
title('Actual and commanded propeller speed (rpm)'); xlabel('time (s)');
subplot(313)
plot(t,delta,t,delta_c,'linewidth',2);
title('Actual and commanded rudder angles (deg)'); xlabel('time (s)');

figure(3) 
figure(gcf)
subplot(211)
plot(t,u,'linewidth',2);
title('Actual surge velocity (m/s)'); xlabel('time (s)');
subplot(212)
plot(t,v,'linewidth',2);
title('Actual sway velocity (m/s)'); xlabel('time (s)');

beta_c_arg = v./u;
beta_c = atan(beta_c_arg);
v_relative = nu_r_Data(:,2);
root_arg = nu_r_Data(:,1).^2 + nu_r_Data(:,2).^2;
U_relative = sqrt(root_arg);
beta_arg= v_relative./U_relative;
beta = asin(beta_arg);

figure(4)
plot(beta_c)
hold on
plot(beta)
legend("beta_c","beta")
hold off

