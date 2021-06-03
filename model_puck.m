clear;


d_t = 0.05;
t_sim = 10; % [s]
sim_steps = round(t_sim / d_t);
noise = randi([-10, 10], sim_steps, 1);

%physical parameters
width = 1.5; %[m]
length = 2; %[m]
m = 0.05; % mass of the puck [kg]
g = 9.80665;
% friction_coef = 0.01034409; % for aerohockey table
friction_coef = 0.3;
a = friction_coef * g;


%Initial conditions, ideal System
t = 0;
x_0 = length/10; % [m]
y_0 = width/2; % [m]
v0 = 1; % m/s
alpha = (pi/180) * 37.5; % [rad] angle to x-axis

v_x0 = v0 * cos(alpha);
v_y0 = v0 * sin(alpha);

ax = a * cos(alpha);
ay = a * sin(alpha);

X_vec = zeros(4, sim_steps);
X_vec(:, 1) = [x_0, y_0, v_x0, v_y0].';
A_mat = [1 0 d_t 0; 0 1 0 d_t; 0 0 1 0; 0 0 0 1]; 
U_vec = [-ax; -ay];
B_mat = [(0.5 * d_t * d_t) 0; 0 (0.5 * d_t * d_t); d_t 0; 0 d_t];

for step = 1:sim_steps
    X_vec(:, step + 1) = A_mat * X_vec(:, step) + B_mat * U_vec;
	
	%bounce off the X walls
	if X_vec(1, step + 1) > length || X_vec(1, step + 1) < 0
		A_mat(1, 3) = A_mat(1,3) * (-1);
		B_mat(1, 1) = B_mat(1, 1) * (-1);
		X_vec(:, step + 1) = A_mat * X_vec(:, step) + B_mat * U_vec;
	end
	%bounce off the Y walls
	if X_vec(2, step + 1) > width || X_vec(2, step + 1) < 0
		A_mat(2, 4) = A_mat(2, 4) * (-1);
		B_mat(1, 2) = B_mat(1, 2) * (-1);
		X_vec(:, step + 1) = A_mat * X_vec(:, step) + B_mat * U_vec;
	end
	
	%stop if v < 0
	if X_vec(3, step + 1) < 0
		X_vec(3, step + 1 ) = 0;
		U_vec = [0;0];
	end
	if X_vec(4, step + 1) < 0
		X_vec(4, step + 1 ) = 0;
		U_vec = [0;0];
	end
end

%plot ideal system
% subplot(2,1,1);
plot(X_vec(1,:), X_vec(2,:), 'LineWidth', 2);
axis([0, length, 0, width])
title('Bewegung des Pucks', 'FontSize', 14)
xlabel('x/m', 'FontSize', 14)
ylabel('y/m', 'FontSize', 14)
% str1=['v_{x0} = ' num2str(vx_0) 'm/s'];
% str2=['v_{y0} = ' num2str(vy_0) 'm/s'];
% str3=['t_{end} = ' num2str(t(k)) 's'];
% text(0.03, 0.7, str1, 'FontSize', 12)
% text(0.03, 0.6, str2, 'FontSize', 12)
% text(0.03, 0.5, str3, 'FontSize',12)