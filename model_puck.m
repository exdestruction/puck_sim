clear;


d_t = 0.01;
t_sim = 10; % [s]
sim_steps = round(t_sim / d_t);
noise = randi([-5, 5], 2, sim_steps);
noise = noise / 120;

%physical parameters
width = 1.5; %[m]
length = 2; %[m]
m = 0.05; % mass of the puck [kg]
g = 9.80665;
friction_coef = 0.01034409; % for aerohockey table
a = friction_coef * g;


%Initial conditions, ideal System
t = 0;
x_0 = length/10; % [m]
y_0 = width/2; % [m]
v0 = 1; % m/s
alpha = (pi/180) * (-38.9); % [rad] angle to x-axis

v_x0 = v0 * cos(alpha);
v_y0 = v0 * sin(alpha);

a_x = a * cos(alpha);
a_y = a * sin(alpha);

X_vec = zeros(4, sim_steps);
X_vec(:, 1) = [x_0, y_0, v_x0, v_y0].';
A_mat = [1 0 d_t 0; 0 1 0 d_t; 0 0 1 0; 0 0 0 1]; 
U_vec = [-a_x; -a_y];
B_mat = [(0.5 * d_t * d_t) 0; 0 (0.5 * d_t * d_t); d_t 0; 0 d_t];

%simulation ideal system trajectory
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
	if abs(X_vec(3, step + 1)) < 0.0005
		X_vec(3, step + 1 ) = 0;
		U_vec = [0;0];
	end
	if abs(X_vec(4, step + 1)) < 0.0005
		X_vec(4, step + 1 ) = 0;
		U_vec = [0;0];
	end
end

H_mat = [1 0 0 0; 0 1 0 0]; % Measurement matrix (Messmatrix)
Z_meas = zeros(2, sim_steps); % Z - measurements

%initial estimation
X_predicted = [0 0 0 0].'; % state of the system after prediction step, mu_hat
X_corrected = [0 0 0 0].'; % state of the system after correction step, mu
K_mat = zeros(4,4); %Kalman matrix
Q_mat = [0.01 0; 0 0.01]; % variance matrix of measurenemts error
R_mat = [0.00001 0 0 0;0 0.00001 0 0;0 0 0.00001 0; 0 0 0 0.00001]; %covariance matrix 
P_predicted = diag([1 1 1 1]); % covariance matrix (not measurement's error), sigma_hat matrix
I_mat = diag([1 1 1 1]);

%simulation estimated trajectory
for step = 1:sim_steps
	
	% Measurements (Messung)
	Z_meas(:, step) = H_mat * X_vec(:, step) + noise(:, step);
	
    % Correction with the measurements (Korrektur mit der Messung)
    K_mat = P_predicted * H_mat' / (H_mat * P_predicted * H_mat' + Q_mat);
    X_corrected(:, step) = X_predicted(:, step) + K_mat * (Z_meas(:,step) - H_mat * X_predicted(:, step)); 
    P_corrected = (I_mat - K_mat * H_mat) * P_predicted; %sigma matrix
    
    % Prediction (Prädiktion)
    X_predicted(:, step + 1) = A_mat * X_corrected(:, step) + B_mat * U_vec;
	P_predicted = A_mat * P_corrected * A_mat' + R_mat;
end

%plot ideal system
subplot(2,1,1);
plot(X_vec(1,:), X_vec(2,:), 'LineWidth', 2);
axis([0, length, 0, width])
title('Trajectory', 'FontSize', 14)
xlabel('x, [m]', 'FontSize', 14)
ylabel('y, [m]', 'FontSize', 14)


% Kalman Filter
subplot(2,1,2);
plot(X_vec(1,:), X_vec(2,:), '--', ...
	Z_meas(1,:), Z_meas(2,:), 'g+' , ...
	X_predicted(1, 2:end), X_predicted(2, 2:end), 'r', 'LineWidth', 1);
axis([0, length, 0, width])
title('Estimated trajectory', 'FontSize', 14)
xlabel('x, [m]', 'FontSize', 14)
ylabel('y, [m]', 'FontSize', 14)
