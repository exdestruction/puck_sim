clear;
d_t = 0.01;
t_sim = 10; % [s]
sim_steps = round(t_sim / d_t);
noise = normrnd(0, 1, 2, sim_steps);
noise = noise / 100;


% physical parameters
width_table = 1.5; %[m]
length_table = 2; %[m]
m = 0.05; % mass of the puck [kg]
g = 9.80665;

f = 1; % friction coef for aerohockey table

% Acceleration and also a parameter in integral functions
J = f * g; 

%puck geometry
h = 0.5; % m
r = 1; % m
delta = 2 * f * h / r;


%Calculating integrals on [0; pi/2] to reduce computing and increase
%stability
[p_a_res, p_b_res, p_a_hat_res, p_c_hat_res, p_c_res] = ...
	calculate_integrals(r);

% functions for calculating integral functions of phi
p_d = @(phi) (1 - delta.^2 .* sin(phi).^2 .* ...
	p_b_res(fix(phi/0.01)+1) .* p_c_res(fix(phi/0.01)+1));


p_d_res = [];
p_parall_res = [];
p_omega_res = [];
p_orth_res = [];
for phi = 0:0.01:pi/2
	p_d_res = [p_d_res; p_d(phi)];
end

pause(1);
p_parall =@(phi) p_a_res(fix(phi/0.01)+1) ./ p_d_res(fix(phi/0.01)+1);
p_omega =@(phi) 2.*(p_a_hat_res(fix(phi/0.01)+1) + ...
	(delta.^2 * cos(phi).^2 .* p_a_res(fix(phi/0.01)+1) ...
	.* p_b_res(fix(phi/0.01)+1) .* p_c_hat_res(fix(phi/0.01)+1) ...
	./ p_d_res(fix(phi/0.01)+1)));
p_orth = @(phi) delta .* p_a_res(fix(phi/0.01)+1) .* p_b_res(fix(phi/0.01)+1) ...
	./ p_d_res(fix(phi/0.01)+1);

for phi = 0:0.01:pi/2
	p_parall_res = [p_parall_res; p_parall(phi)];
	p_omega_res = [p_omega_res; p_omega(phi)];
	p_orth_res = [p_orth_res; p_orth(phi)];
end


Q_dot =@(phi) -((cos(phi)).^2 .* p_parall_res(fix(phi/0.01)+1) ...
	+ (sin(phi)).^2 .* p_omega_res(fix(phi/0.01)+1))* J;
phi_dot =@(phi) sin(phi) .* cos(phi) .* (p_parall_res(fix(phi/0.01)+1) ...
	- p_omega_res(fix(phi/0.01)+1)) .* J;
psi_dot =@(phi) -sin(phi) .* p_orth_res(fix(phi/0.01)+1) .* J;


%Initial conditions, ideal System
t = 0;
x_0 = length_table/10; % [m]
y_0 = width_table/2; % [m]

% components of motion laws when type D trajectory
phi_0 = 0.975 + (0.992 - 0.975) * delta;

%initial data
v_0 = 6; % m/s
w_0 = 20; % rad/s

psi_0 = (pi/180) * (42); % [rad] angle to x-axis
psi = psi_0;

phi_start = atan(r * w_0 / v_0);
phi = phi_start;

Q_0 = v_0 / cos(phi_start);
Q = Q_0;

v_x0 = v_0 * cos(psi_0);
v_y0 = v_0 * sin(psi_0);

X_vec(:, 1) = [x_0, y_0, v_x0, v_y0].';
for step = 1:sim_steps-1
	phi_next = phi(end) + phi_dot(phi(end))/Q(end) * d_t;
	
	Q_next = Q(end) + Q_dot(phi_next) * d_t;
	
	%stop if v < 0 (Q < 0)
	if Q_next < 0
		break
	end
	
	psi_next = psi(end) + psi_dot(phi_next)/Q_next * d_t;
	
	v_x = Q_next .* cos(phi_next) .* cos(psi_next);
	v_y = Q_next .* cos(phi_next) .* sin(psi_next);
	
	x = X_vec(1, step) + v_x * d_t;
	y = X_vec(2, step) + v_y * d_t;
	
	X_vec(:, step + 1) = [x; y; v_x; v_y];
	
	% bounce off the X walls
	if x > length_table || x < 0
		psi_next = pi - psi_next;
	end
	% bounce off the Y walls
	if y > width_table || y < 0
		psi_next = -psi_next;
	end
	psi = [psi; psi_next];
	Q = [Q; Q_next];
	phi = [phi; phi_next];

end

H_mat = [1 0 0 0; 0 1 0 0]; % Measurement matrix (Messmatrix)
Z_meas = zeros(2, sim_steps); % Z - measurements

%Model matrices
% a = f*g;
% ax = a;
% ay = a;
A_mat = [1 0 d_t 0; 0 1 0 d_t; 0 0 1 0; 0 0 0 1];
% B_mat = [(1/2 .* d_t.^2) 0; 0 (1/2 .* d_t.^2); d_t 0; 0 d_t];
% U_vec = [-ax; -ay];

%initial estimation
X_predicted = [0 0 0 0].'; % state of the system after prediction step, mu_hat
X_corrected = [0 0 0 0].'; % state of the system after correction step, mu
K_mat = zeros(4,4); %Kalman matrix
Q_mat = [0.1 0; 0 0.1]; % variance matrix of measurenemts error
R_mat = [0.04 0 0 0; 0 0.04 0 0; 0 0 0.01 0; 0 0 0 0.01]; %covariance matrix 
P_predicted = diag([1 1 1 1]); % covariance matrix (not measurement's error), sigma_hat matrix
I_mat = diag([1 1 1 1]);

%simulation estimated trajectory
for step = 1:length(X_vec(1,:))
	% Measurements (Messung)
	Z_meas(:, step) = H_mat * X_vec(:, step) + noise(:, step);
	
    % Correction with the measurements (Korrektur mit der Messung)
    K_mat = P_predicted * H_mat' / (H_mat * P_predicted * H_mat' + Q_mat);
    X_corrected(:, step) = X_predicted(:, step) + ...
		K_mat * (Z_meas(:,step) - H_mat * X_predicted(:, step)); 
    P_corrected = (I_mat - K_mat * H_mat) * P_predicted; %sigma matrix
    
    % Prediction (Prädiktion)
	X_predicted(:, step + 1) = A_mat * X_corrected(:, step);
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

function list = clean_data(data)
	list = data;
	prev_idx = 0;
	next_idx = 0;
	for ii = 1:length(list)
		if isinf(list(ii)) || isnan(list(ii))
			prev_idx = ii-1;
			for jj = ii:length(list)
				if ~(isinf(list(jj)) || isnan(list(jj)))
					next_idx = jj;
					break
				end
			end
		end
		for kk = 1:(next_idx - prev_idx)
			derivative = (list(next_idx) - list(prev_idx)) / (next_idx - prev_idx);
			list(prev_idx+kk) = list(prev_idx) + kk * derivative;
		end
		ii = next_idx;
	end
end

function [p_a_res, p_b_res, p_a_hat_res, p_c_hat_res, p_c_res] = ...
	calculate_integrals(r)
	% calculates values for p_a, p_b, p_a_hat, p_c_hat, p_c 
	% for phi on [0; pi/2]
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
	p_c =@(phi) p_b(phi) - p_a_hat(phi);


	
	p_a_res = [];
	p_b_res = [];
	p_a_hat_res = [];
	p_c_hat_res = [];
	p_c_res = [];
	p_d_res = [];
	for phi = 0:0.01:pi/2
		p_a_data = p_a(phi);
		p_b_data = p_b(phi);
		p_a_hat_data = p_a_hat(phi);
		p_c_hat_data = p_c_hat(phi);
		p_c_data = p_c(phi);
		p_a_res = [p_a_res; real(p_a_data)];
		p_b_res = [p_b_res; real(p_b_data)];
		p_a_hat_res = [p_a_hat_res; real(p_a_hat_data)];
		p_c_hat_res = [p_c_hat_res; real(p_c_hat_data)];
		p_c_res = [p_c_res; real(p_c_data)];	

	end
	p_a_res = clean_data(p_a_res);
	p_b_res = clean_data(p_b_res);
	p_a_hat_res = clean_data(p_a_hat_res);
	p_c_hat_res = clean_data(p_c_hat_res);
	p_c_res = clean_data(p_c_res);
end



