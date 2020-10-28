addpath(genpath("C:\Users\jacob\Documents\Student\Fartøy\MSS"));
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
sim("succ_loop_sim.slx")
