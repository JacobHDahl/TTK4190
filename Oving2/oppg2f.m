addpath(genpath("C:\Users\jacob\Documents\Student\Fartøy\TTK4190\MSS"));
Kp_phi = -2;
Kd_phi = 1.935;
Ki_phi = 0;

a2 = -0.65;
a1 = 2.87;
Vg = 161.1;
g = 9.81;

omega_n_phi = 1.14;
Wx = 10;
omega_n_x = omega_n_phi/Wx;
zeta_x = 2;

Kp_x = 2*zeta_x*omega_n_x*Vg/g;
Ki_x = omega_n_x^2 *Vg/g;

dist = 1.5;
e_phi_sat = 15;
act_sat = pi/12;

A = [-0.322, 0.052, 0.028, -1.12,  0.002;
       0,       0,    1,  -0.001,      0;
       -10.6,   0,  -2.87,  0.46,  -0.65;
       6.87,   0,  -0.04,  -0.32,  -0.02;
       0,      0,      0,      0,   -7.5];
   
B = [0,0,0,0,7.5]';
C = [1,0,0,0,0;
    0,1,0,0,0;
    0,0,1,0,0;
    0,0,0,1,0];
    
sim("succ_loop_whole_model")
