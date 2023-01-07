function X = CWHdyn(t,states)

global n

vx = states(1);
vy = states(2);
vz = states(3);
x = states(4);
y = states(5);
z = states(6);

% Hill's equations of motion 
x_ddot = 3*n^2*x + 2*n*vy;
y_ddot = -2*n*vx;
z_ddot = -n^2*z;
x_dot  = vx;
y_dot  = vy;
z_dot  = vz;

X = [x_ddot y_ddot z_ddot x_dot y_dot z_dot]';

end