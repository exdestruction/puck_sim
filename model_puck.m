clear;
d_t = 0.01;
t_sim = 10; % [s]
sim_steps = round(t_sim / d_t);
noise = randi([-5, 5], 2, sim_steps);
noise = noise / 120;
% noise = randi([0, 0], 2, sim_steps)
% physical parameters
width_table = 1.5; %[m]
length_table = 2; %[m]
m = 0.05; % mass of the puck [kg]
g = 9.80665;
f = 0.01034409; % friction coef for aerohockey table

%puck geometry
h = 1;
r = 5;

% components of motion laws when type D trajectory
delta = 2 * f * h / r;
phi_0 = 0.975 + (0.992 - 0.975) * delta;
alpha = cos(phi);
beta = sin(phi);
p_a = @(phi) real(integral2(@(ro, nu) ...
	(cos(phi) - sin(phi) .* (ro/r) .* sin(nu)).^2 .* ro ./(pi .* ((sin(phi)).^2 * (ro/r).^2 - 2 .* (cos(phi)) .* (sin(phi)) .* (ro/r) .* sin(nu) + (cos(phi)).^2).^(3/2)),...
	0, r, 0, 2*pi));
p_b = @(phi) real(integral2(@(ro, nu) ...
	((ro/r) .* cos(nu)).^2 .* ro ./ (pi .* ((sin(phi)).^2 * (ro/r).^2 - 2 .* (cos(phi)) .* (sin(phi)) .* (ro/r) .* sin(nu) + (cos(phi)).^2).^(1/2)),...
	0, r, 0, 2*pi));
p_a_hat = @(phi) real(integral2(@(ro, nu) ...
	( (sin(phi) .* (ro/r) - cos(phi) .* sin(nu)).^2) .*(ro/r).^2 .* ro ./ (pi .* ((sin(phi)).^2 * (ro/r).^2 - 2 .* (cos(phi)) .* (sin(phi)) .* (ro/r) .* sin(nu) + (cos(phi)).^2).^(3/2)),...
	0, r, 0, 2*pi));
p_c_hat = @(phi) real(integral2(@(ro, nu) ...
	(sin(phi).^2 .* (ro/r).^2 .* cos(2*nu) + 2 .* cos(phi) .* sin(phi) .* (ro/r) .* sin(nu).^3 - cos(phi).^2 .* sin(nu).^2) .*(ro/r).^2 .* ro ./ (pi .* ((sin(phi)).^2 * (ro/r).^2 - 2 .* (cos(phi)) .* (sin(phi)) .* (ro/r) .* sin(nu) + (cos(phi)).^2).^(3/2)),...
	0, r, 0, 2*pi));
p_c =@(phi) p_b(r, phi) - p_a_hat(r, phi);

p_d = @(phi) 1 - delta.^2 .* sin(phi).^2 .* p_b(r, phi) .* p_c(r, phi);

p_parall =@(phi) p_a(phi) / p_d(phi);
p_omega =@(phi) 2.*(p_a_hat(phi) + (delta.^2 * cos(phi).^2 .* p_a(phi) .* p_b(phi) .* p_c_hat(phi) ./ p_d(phi)));
p_orth = @(phi) delta .* p_a(phi) .* p_b(phi) ./ p_d(phi);
J = f * N / m;


% r * w / v = tan(phi);

% motion laws
v_x = Q * cos(phi) * cos(psi);
v_y = Q * cos(phi) * sin(psi);
omega = Q * sin(phi) / r;

Q_dot = -((cos(phi))^2 * p_parall  + (sin(phi))^2 * p_omega)* J;
phi_dot = sin(phi) * cos(phi) * (p_parall - p_omega) * J / Q;
psi_dot = -sin(phi) * p_orth * J / Q;


%Initial conditions, ideal System
t = 0;
x_0 = length_table/10; % [m]
y_0 = width_table/2; % [m]
v_0 = 1; % m/s
omega_0 = 10; % rad/s
psi_0 = (pi/180) * (41.4); % [rad] angle to x-axis


X_vec = zeros(2, sim_steps);
X_vec(:, 1) = [x_0, y_0].';
A_mat = [1 0; 0 1]; 
U_vec = [];
B_mat = [];

%simulation ideal system trajectory
for step = 1:sim_steps
    X_vec(:, step + 1) = A_mat * X_vec(:, step) + B_mat * U_vec;
	
	%bounce off the X walls
	if X_vec(1, step + 1) > length_table || X_vec(1, step + 1) < 0
		A_mat(1, 3) = A_mat(1,3) * (-1);
		B_mat(1, 1) = B_mat(1, 1) * (-1);
		X_vec(:, step + 1) = A_mat * X_vec(:, step) + B_mat * U_vec;
	end
	
	%bounce off the Y walls
	if X_vec(2, step + 1) > width_table || X_vec(2, step + 1) < 0
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
R_mat = [0.0025 0 0 0; 0 0.0025 0 0; 0 0 0 0; 0 0 0 0]; %covariance matrix 
P_predicted = diag([1 1 1 1]); % covariance matrix (not measurement's error), sigma_hat matrix
I_mat = diag([1 1 1 1]);

%simulation estimated trajectory
for step = 1:sim_steps
	% Measurements (Messung)
	Z_meas(:, step) = H_mat * X_vec(:, step) + noise(:, step);
	
    % Correction with the measurements (Korrektur mit der Messung)
    K_mat = P_predicted * H_mat' / (H_mat * P_predicted * H_mat' + Q_mat);
    X_corrected(:, step) = X_predicted(:, step) + ...
		K_mat * (Z_meas(:,step) - H_mat * X_predicted(:, step)); 
    P_corrected = (I_mat - K_mat * H_mat) * P_predicted; %sigma matrix
    
    % Prediction (Prädiktion)
    X_predicted(:, step + 1) = A_mat * X_corrected(:, step) + B_mat * U_vec;
	P_predicted = A_mat * P_corrected * A_mat' + R_mat;
end

%plot ideal system
subplot(2,1,1);
plot(X_vec(1,:), X_vec(2,:), 'LineWidth', 2);
axis([0, length_table, 0, width_table])
title('Trajectory', 'FontSize', 14)
xlabel('x, [m]', 'FontSize', 14)
ylabel('y, [m]', 'FontSize', 14)


% Kalman Filter
subplot(2,1,2);
plot(X_vec(1,:), X_vec(2,:), '--', ...
	Z_meas(1,:), Z_meas(2,:), 'g+' , ...
	X_predicted(1, 2:end), X_predicted(2, 2:end), 'r', 'LineWidth', 1);
axis([0, length_table, 0, width_table])
title('Estimated trajectory', 'FontSize', 14)
xlabel('x, [m]', 'FontSize', 14)
ylabel('y, [m]', 'FontSize', 14)

