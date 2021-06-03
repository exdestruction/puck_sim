clear;
g = 9.80665;
del_t = 0.005;
k = 0;
n = 134;
random = randi([-5, 5], n, 1);
%−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
%Fuer ideales System
t = 0;
x_ideal = 0;
y_ideal = 0;

v_0 = 3.8;
vx_0 = 0.5 * v_0;
vy_0 = 0.86 * v_0;


vx_ideal = vx_0;

vy_ideal = vy_0;
%−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
%Fuer Kalman Filter
Mat_K = zeros(4,4);

%Startschaetzung
Mat_Pm = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
Mat_R = [0.01 0; 0 0.01];
Mat_H = [1 0 0 0; 0 1 0 0];

%Startschaetzung
Vec_Xm = [0 0 1 1].' ;
Array_Xm = zeros(4, n+1);
Array_Xm(1) = Vec_Xm(1);
Array_Xm(2) = Vec_Xm(2);
Array_Xm(3) = Vec_Xm(3);
Array_Xm(4) = Vec_Xm(4);
x_kal = 0;
y_kal = 0;
x_mes = 0;
y_mes = 0;
Vec_Xp = [0 0 0 0].';
Vec_Y = [0 0].';
Mat_Pp = zeros(4,4);
Mat_I = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
Mat_A = [1 0 del_t 0; 0 1 0 del_t; 0 0 1 0; 0 0 0 1];
Mat_B = [0 0; 0 (0.5 * del_t * del_t); 0 0; 0 del_t];
Vec_u = [0 -g].';
Mat_Q = [0.00001 0 0 0;0 0.00001 0 0;0 0 0.00001 0; 0 0 0 0.00001];
for k = 1:n
    %−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
    %ideales System
    t(k+1) = t(k)+ del_t;
    x_ideal(k+1) = x_ideal(k)+ del_t * vx_ideal(k);
    y_ideal(k+1) = y_ideal(k)+ del_t * vy_ideal(k) - 0.5 * g * del_t * del_t;
    vx_ideal(k+1) = vx_ideal(k);
    vy_ideal(k+1) = vy_ideal(k) - g * del_t;
    
    %−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
    %Messung
    Vec_Y(1) = x_ideal(k);
    rand = (random(k)/100);
    Vec_Y(2) = y_ideal(k) + rand;
    
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
%ideal
subplot(2,1,1);
plot(x_ideal, y_ideal, 'LineWidth', 2);
axis([0, 1.4, 0, 0.8])
title('IdealeWurfparabel', 'FontSize', 14)
xlabel('x/m', 'FontSize', 14)
ylabel('y/m', 'FontSize', 14)
str1=['v_{x0} = ' num2str(vx_0) 'm/s'];
str2=['v_{y0} = ' num2str(vy_0) 'm/s'];
str3=['t_{end} = ' num2str(t(k)) 's'];
text(0.03, 0.7, str1, 'FontSize', 12)
text(0.03, 0.6, str2, 'FontSize', 12)
text(0.03, 0.5, str3, 'FontSize',12)

%KalmanFilter
subplot(2,1,2);
plot(x_ideal, y_ideal, '--', ...
     x_mes, y_mes,'g+', ...
     x_kal, y_kal, 'r', 'LineWidth', 2);
axis([0, 1.4, 0, 0.8])
xlabel('x/m', 'FontSize', 14)
ylabel('y/m', 'FontSize', 14)
title('Wurfparabel rekonstruiert durch Kalman Filter mit Suchfenster', 'FontSize', 14)
str1=['v_{x0} = ' num2str(vx_0) 'm/s'];
str2=['v_{y0} = ' num2str(vy_0) 'm/s'];
str3=['t_{end} = ' num2str(t(k)) 's'];
text(0.03, 0.7, str1, 'FontSize', 12)
text(0.03,0.6, str2, 'FontSize', 12)
text(0.03,0.5,str3, 'FontSize', 12)
