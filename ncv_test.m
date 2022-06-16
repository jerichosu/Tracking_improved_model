% Set up model 1 for mode 1: NCV model, one dimensional, T=1
% Start with continuous-time; convert to discrete-time
% Ac = [0 1; 0 0]; 
% Bc = [0; 1]; 
% Swc = 1;
% Z = [-Ac Bc*Swc*Bc'; zeros(size(Ac)) Ac'];
% C = expm(Z*1); 
% 
% 
% Ad = C(3:4,3:4)'; 
% Sw1 = Ad*C(1:2,3:4);
% Bd = [1^2/2; 1]; 
% Cd = [1 0]; 
% Sv1 = 0.2;
% 
% 
% %   Generate the true system data
% xtrue = zeros([2,201]);
% z = zeros(1,200);
% for k = 1:200,
% xtrue(:,k+1) = Ad*xtrue(:,k) + chol(Sw1,'lower')*randn([2 1]);
% z(k) = Cd*xtrue(:,k) + sqrt(Sv1)*randn(1);
% end
% 
% figure(1); clf; plot(xtrue(1,:)); hold on;


clear all
close all
clc


maxT=200; % max sim time

x = zeros(4,maxT); % storage

% initial position, velocity
x(:,1) = [0; 
          5; 
          0; 
          7];


T = 1; % sampling time


A = [1 T 0 0; %x -> position information
     0 1 0 0; %x_dot -> speed information
     0 0 1 T; %y -> position information
     0 0 0 1];%y_dot -> speed information
 
 
F = [T^2/2 0; 
     T     0; 
     0     T^2/2; 
     0     T];
 
C = [1 0 0 0;
     0 0 1 0];
 
G = [1 0;
     0 1];

% state noise covariance matrix
V = [0.85 0;
     0    0.3];
% measurement noise covariance matrix
W = [0.28 0;
     0    0.65];
           

for k = 1:maxT, % simulate model
    x(:,k+1) = A*x(:,k) + F*V*randn(2,1);
    y(:,k) = C*x(:,k) + G*W*randn(2,1);
end

% figure,
% plot(x(1,:),x(3,:));
% title('Traj. of Discrete-time NCV Model');
% xlabel('x position'); ylabel('y position');

% figure,
% plot(y(1,:),y(2,:));
% title('Traj. of Discrete-time NCV Model');
% xlabel('x position'); ylabel('y position');


% set up Kalman filters to estiamte the position

x_hat(:,1) = [0.5 20 7 4]';
P(:,:,1) = 10*eye(4);

kmax = 200;
for k = 1:kmax
    
    K = A*P(:,:,k)*C' * inv(C*P(:,:,k)*C' + G*W*G');
    x_hat(:,k+1) = A*x_hat(:,k) + K*(y(:,k) - C*x_hat(:,k));
    P(:,:,k+1) = A*P(:,:,k)*A' + K*C*P(:,:,k)*A' + F*V*F';
    
end


figure,
plot(x(1,:),x(3,:));
title('Traj. of Discrete-time NCV Model');
xlabel('x position'); ylabel('y position');


figure,
plot(x_hat(1,:),x_hat(3,:));
title('Traj. of Discrete-time NCV Model');
xlabel('x position estimate'); ylabel('y position estimate');















