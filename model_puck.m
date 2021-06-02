clear;
g = 9.80665;
d_t = 0.005;
t_sim = 10; %seconds
sim_steps = round(t_sim / d_t);
noise = randi([-10, 10], n, 1);

%Initial conditions, ideal System
t = 0;
x_ideal = 0.05; % meters
y_ideal = 0.05; % meters