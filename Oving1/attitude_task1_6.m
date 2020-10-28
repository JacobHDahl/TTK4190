% M-script for numerical integration of the attitude dynamics of a rigid 
% body represented by unit quaternions. The MSS m-files must be on your
% Matlab path in order to run the script.
%
% System:                      .
%                              q = T(q)w
%                              .
%                            I w - S(Iw)w = tau
% Control law:
%                            tau = constant
% 
% Definitions:             
%                            I = inertia matrix (3x3)
%                            S(w) = skew-symmetric matrix (3x3)
%                            T(q) = transformation matrix (4x3)
%                            tau = control input (3x1)
%                            w = angular velocity vector (3x1)
%                            q = unit quaternion vector (4x1)
%
% Author:                   2018-08-15 Thor I. Fossen and Håkon H. Helgesen

%% USER INPUTS
addpath(genpath("C:\Users\jacob\Documents\Student\Fartøy\MSS"));
h = 0.1;                     % sample time (s)
N  = 4000;                    % number of samples. Should be adjusted

% model parameters
m = 180;
r = 2;
I = m*r^2*eye(3);            % inertia matrix
I_inv = inv(I);

% constants
deg2rad = pi/180;   
rad2deg = 180/pi;

phi = -5*deg2rad;            % initial Euler angles
theta = 10*deg2rad;
psi = -20*deg2rad;

q = euler2q(phi,theta,psi);   % transform initial Euler angles to q

w = [0 0 0]';                 % initial angular rates

table = zeros(N+1,14);        % memory allocation

kp = 20;                       %proportional gain
kd = 400;                      %derivative gain


table_d = zeros(N+1,3);




%% FOR-END LOOP
for i = 1:N+1
   t = (i-1)*h;                  % time
   
   phi_d = 0*deg2rad;                       %desired phi
   theta_d = 15*cos(0.1*t)*deg2rad;        %desired theta
   psi_d = 10*sin(0.05*t)*deg2rad;          %desired psi
   
   phi_d_dot = 0*deg2rad;                   %desired phi_dot
   theta_d_dot = -1.5*sin(0.1*t)*deg2rad;   %desired theta_dot
   psi_d_dot = 0.5*cos(0.05*t)*deg2rad;     %desired psi_dot
   
   Theta_D = [phi_d_dot, theta_d_dot, psi_d_dot]';

   q_d = euler2q(phi_d,theta_d,psi_d); %desired quaternion
   
   inv_matr = eye(4)*-1;                %matrix for inverting the quaternion
   inv_matr(1,1) = 1;
   
   q_d_inv = inv_matr*q_d;              %inverting quaternion
   
   T_inv = [1, 0, -sin(theta);             %T_inv from 2.41 in Fossen
    0, cos(phi), cos(theta)*sin(phi);
    0, -sin(phi), cos(theta)*cos(phi)];

    
    w_d = T_inv * Theta_D;                  %desired omega
    
    w_err = w-w_d;                          %error in omega
   
   
   q_err = quatprod(q_d_inv,q);         %calculating the error-quaternion
   x = [q_err(2:4)' w_err']';
   Kp = -kp*eye(3);
   Kd = -kd*eye(3);
   
   K = [Kp, Kd];
   
   tau = K*x;           % control law

   [phi,theta,psi] = q2euler(q); % transform q to Euler angles
   [J,J1,J2] = quatern(q);       % kinematic transformation matrices
   
   q_dot = J2*w;                        % quaternion kinematics
   w_dot = I_inv*(Smtrx(I*w)*w + tau);  % rigid-body kinetics
   
   table(i,:) = [t q' phi theta psi w' tau'];  % store data in table
   table_d(i,:) = [phi_d,theta_d,psi_d];
   
   q = q + h*q_dot;	             % Euler integration
   w = w + h*w_dot;
   
   q  = q/norm(q);               % unit quaternion normalization
end 

%% PLOT FIGURES
t       = table(:,1);  
q       = table(:,2:5); 
phi     = rad2deg*table(:,6);
theta   = rad2deg*table(:,7);
psi     = rad2deg*table(:,8);
w       = rad2deg*table(:,9:11);  
tau     = table(:,12:14);

phi_d = rad2deg*table_d(:,1);
theta_d = rad2deg*table_d(:,2);
psi_d = rad2deg*table_d(:,3);

phi_err = phi_d - phi;
theta_err = theta_d - theta;
psi_err = psi_d - psi;


figure (1); clf;
hold on;
plot(t, phi, 'b');
plot(t, theta, 'r');
plot(t, psi, 'g');
hold off;
grid on;
legend('\phi', '\theta', '\psi');
title('Euler angles');
xlabel('time [s]'); 
ylabel('angle [deg]');

figure (2); clf;
hold on;
plot(t, w(:,1), 'b');
plot(t, w(:,2), 'r');
plot(t, w(:,3), 'g');
hold off;
grid on;
legend('x', 'y', 'z');
title('Angular velocities');
xlabel('time [s]'); 
ylabel('angular rate [deg/s]');

figure (3); clf;
hold on;
plot(t, tau(:,1), 'b');
plot(t, tau(:,2), 'r');
plot(t, tau(:,3), 'g');
hold off;
grid on;
legend('x', 'y', 'z');
title('Control input');
xlabel('time [s]'); 
ylabel('input [Nm]');

figure (4); clf;
hold on;
plot(t, phi_err, 'b');
plot(t, theta_err, 'r');
plot(t, psi_err, 'g');
hold off;
grid on;
legend('\phi_{err}', '\theta_{err}', '\psi_{err}');
title('Euler angles error');
xlabel('time [s]'); 
ylabel('angle [deg]');