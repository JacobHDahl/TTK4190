addpath(genpath("C:\Users\jacob\Documents\Student\Fartøy\TTK4190\MSS"));

Kp = -2;
Kd = 1.935;
Ki = -1;
a2 = -0.65;
a1 = 2.87;

Ki_range = [0:0.1:10];

sys1 = tf([Kp*a2, Ki*a2],[1,a1+a2*Kd, a2*Kp, a2*Ki]);
sysEvans = tf([a2],[1,a1+a2*Kd,a2*Kp,0]);
P = pole(sys1);
P_real = real(P);
P_imag = imag(P);
[r,k] = rlocus(sysEvans);

figure(1)
rlocus(sysEvans,Ki_range)

figure(2)
plot(P_real,P_imag,'*')
grid()


