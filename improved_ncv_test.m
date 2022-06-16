clear all
close all
clc

% load measurements.mat % this is our sensor measurments (position of x and y)

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

y_raw = y;

y(:, 50:80) = NaN;
y(:, 130:160) = NaN;

figure,
plot(y_raw(1,:),y_raw(2,:),'LineWidth',3);
hold on
plot(y(1,:),y(2,:),'LineWidth',2);
legend('Raw position traj','Detected position traj (have miss detections)')
title('Traj. of Discrete-time NCV Model');
xlabel('traj of x'); ylabel('traj of y');
%%

% set up Kalman filters to estiamte the position
T = 1; % sampling time
A = [1 T 0 0; %x -> position information
     0 1 0 0; %x_dot -> speed information
     0 0 1 T; %y -> position information
     0 0 0 1];%y_dot -> speed information
 
% F = [T^2/2 0; 
%      T     0; 
%      0     T^2/2; 
%      0     T];
 
F = [T^2/2; 
     T; 
     T^2/2; 
     T];
 
 
C = [1 0 0 0;
     0 0 1 0];
 
% G = [1 0;
%      0 1];
 
G = [1;
     1];

% % state noise covariance matrix
% V = [0.85 0;
%      0    0.3];
% % measurement noise covariance matrix
% W = [0.28 0;
%      0    0.65];
 
% state noise covariance matrix
V = 0.5;
% measurement noise covariance matrix
W = 0.2;

x_hat(:,1) = [0.5 20 7 4]';
P(:,:,1) = 10*eye(4);

kmax = 200;
for k = 1:kmax
    
    if isnan(y(1,k))
%         K = A*P(:,:,k)*C' * inv(C*P(:,:,k)*C' + G*W*G');
        x_hat(:,k+1) = A*x_hat(:,k);
        P(:,:,k+1) = A*P(:,:,k)*A' + K*C*P(:,:,k)*A' + F*V*F';
    else
        K = A*P(:,:,k)*C' * inv(C*P(:,:,k)*C' + G*W*G');
        x_hat(:,k+1) = A*x_hat(:,k) + K*(y(:,k) - C*x_hat(:,k));
        P(:,:,k+1) = A*P(:,:,k)*A' + K*C*P(:,:,k)*A' + F*V*F';
    end
    
end


figure,
plot(y_raw(1,:),y_raw(2,:),'LineWidth',3);
hold on
plot(x_hat(1,:),x_hat(3,:), 'g','LineWidth',2);
legend('Raw position traj','estiamted traj')


% title('Traj. of Discrete-time NCV Model');
% xlabel('traj of x'); ylabel('traj of y');
% 
% title('Traj. of Discrete-time NCV Model');
% xlabel('x position estimate'); ylabel('y position estimate');


%%
%%
% set up Kalman filters WITH THE IMPROVED MODEL to estiamte the position
T = 1; % sampling time

% state noise covariance matrix
V = 0.5;
% measurement noise covariance matrix
W = 0.2;
%tunning factor S
S = 1;

Am = [1 T                     0 0; %x -> position information
      0 sqrt((S^2 - V^2)/S^2) 0 0; %x_dot -> speed information
      0 0                     1 T; %y -> position information
      0 0                     0 sqrt((S^2 - V^2)/S^2)];%y_dot -> speed information
 
F = [T^2/2; 
     T; 
     T^2/2; 
     T];
 
 
C = [1 0 0 0;
     0 0 1 0];
 
G = [1;
     1];

x_hat_m(:,1) = [0.5 20 7 4]';
P(:,:,1) = 10*eye(4);

kmax = 200;
for k = 1:kmax
    
    if isnan(y(1,k))
        K = Am*P(:,:,k)*C' * inv(C*P(:,:,k)*C' + G*W*G');
        x_hat_m(:,k+1) = Am*x_hat_m(:,k);
        P(:,:,k+1) = Am*P(:,:,k)*Am' + K*C*P(:,:,k)*Am' + F*V*F';
    else
        K = Am*P(:,:,k)*C' * inv(C*P(:,:,k)*C' + G*W*G');
        x_hat_m(:,k+1) = Am*x_hat_m(:,k) + K*(y(:,k) - C*x_hat_m(:,k));
        P(:,:,k+1) = Am*P(:,:,k)*Am' + K*C*P(:,:,k)*Am' + F*V*F';
    end
    
end


figure,
plot(y_raw(1,:),y_raw(2,:),'LineWidth',3);
hold on
plot(x_hat(1,:),x_hat(3,:), 'g','LineWidth',2);
hold on
plot(x_hat_m(1,:),x_hat_m(3,:), '--','LineWidth',2);
legend('Raw position traj','estiamted traj using NCV model', 'estiamted traj using improved NCV model')
