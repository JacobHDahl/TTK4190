m = 180;
R_33 = 2.0;
kp = 2;
kd = 40;
MRSQ = m*R_33*R_33;

A = [0, 0, 0, 0.5, 0, 0;
    0, 0, 0, 0, 0.5, 0,;
    0, 0, 0, 0, 0, 0.5;
    -kp/MRSQ, 0, 0, -kd/MRSQ, 0, 0;
    0, -kp/MRSQ, 0, 0, -kd/MRSQ, 0;
    0, 0, -kp/MRSQ, 0, 0, -kd/MRSQ];

e = eig(A)