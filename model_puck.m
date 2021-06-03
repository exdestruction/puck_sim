clear;


d_t = 0.05;
t_sim = 10; % [s]
sim_steps = round(t_sim / d_t);
noise = randi([-10, 10], 2, sim_steps);
noise = noise / 100;

%physical parameters
width = 1.5; %[m]
length = 2; %[m]
m = 0.05; % mass of the puck [kg]
g = 9.80665;
friction_coef = 0.01034409; % for aerohockey table
% friction_coef = 0.3;
a = friction_coef * g;


%Initial conditions, ideal System
t = 0;
x_0 = length/10; % [m]
y_0 = width/2; % [m]
v0 = 1; % m/s
alpha = (pi/180) * 37.5; % [rad] angle to x-axis

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
	if X_vec(3, step + 1) < 0
		X_vec(3, step + 1 ) = 0;
		U_vec = [0;0];
	end
	if X_vec(4, step + 1) < 0
		X_vec(4, step + 1 ) = 0;
		U_vec = [0;0];
	end
end

H = [1 0 0 0; 0 1 0 0]; % Measurement matrix (Messmatrix)
Y_meas = zeros(2, sim_steps); % preallocation memory 

%initial estimation
K_mat = zeros(4,4); %Kalman matrix
R_mat = [0.01 0; 0 0.01]; % covariance matrix
Pm_mat = []; 

%simulation estimated trajectory
for step = 1:sim_steps
	
	% Measurement (Messung)
	Y_meas(:, step) = H * X_vec(:, step) + noise(:, step);
    
%     %Wurfparabel mit Stoerung
%     if (0.7 < x_ideal(k) && x_ideal(k)<0.8)
%         Vec_Y(2) = 0.1 + rand;
%         Mat_R(1) = 1;
%         Mat_R(4) = 1;
%     else
%         Mat_R(1)=0.01;
%         Mat_R(4)=0.01;
%     end
    
%     %Ende Wurfparabel mit Stoerung
%     if(Vec_Y(2)<0)
%         Vec_Y(2) = 0;
%     end
%     x_mes(k) = Vec_Y(1);
%     y_mes(k) = Vec_Y(2);
%     
%     %−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
%     %Korrektur mit der Messung
%     Inverse = inv(Mat_R + Mat_H * Mat_Pm * Mat_H');
%     Mat_K = Mat_Pm * Mat_H' * Inverse;
%     Vec_Xp = Vec_Xm + Mat_K * (Vec_Y - Mat_H * Vec_Xm);
%     Mat_Pp = (Mat_I - Mat_K * Mat_H) * Mat_Pm;
%     
%     %−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
%     %Prädiktion
%     Vec_Xm = Mat_A * Vec_Xp + Mat_B * Vec_u;
%     Array_Xm((k*4) + 1) = Vec_Xm(1);
%     Array_Xm((k*4) + 2) = Vec_Xm(2);
%     Array_Xm((k*4) + 3) = Vec_Xm(3);
%     Array_Xm((k*4) + 4) = Vec_Xm(4);
%     x_kal(k) = Vec_Xm(1);
%     y_kal(k) = Vec_Xm(2);
%     Mat_Pm = Mat_A * Mat_Pp * Mat_A' + Mat_Q;
end

%plot ideal system
subplot(2,1,1);
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

%KalmanFilter
subplot(2,1,2);
plot(X_vec(1,:), X_vec(2,:), '--', ...
	Y_meas(1,:), Y_meas(2,:), 'g+' ,'LineWidth', 1);
axis([0, length, 0, width])
title('Bewegung des Pucks', 'FontSize', 14)
xlabel('x/m', 'FontSize', 14)
ylabel('y/m', 'FontSize', 14)
% str1=['v_{x0} = ' num2str(vx_0) 'm/s'];
% str2=['v_{y0} = ' num2str(vy_0) 'm/s'];
% str3=['t_{end} = ' num2str(t(k)) 's'];
% text(0.03, 0.7, str1, 'FontSize', 12)
% text(0.03,0.6, str2, 'FontSize', 12)
% text(0.03,0.5,str3, 'FontSize', 12)