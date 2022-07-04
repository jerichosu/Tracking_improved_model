% clear all
% close all
% clc
% 
% alpha_0 = 0.99;
% alpha(:,1) = 0.95;
% 
% kmax = 200;
% 
% for k = 1:kmax
%     alpha(:,k+1) = alpha_0 * alpha(:,k) + (1 - alpha_0);
% end
% 
% figure,
% plot(1:kmax+1, alpha)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc


kmax=200; % max sim time

x = zeros(4,kmax); % storage

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
           

for k = 1:kmax, % simulate model
    x(:,k+1) = A*x(:,k) + F*V*randn(2,1);
    y(:,k) = C*x(:,k) + G*W*randn(2,1);
end

y_raw = y;

y(:, 50:80) = NaN;
y(:, 130:160) = NaN;

figure,
plot(y_raw(1,:),y_raw(2,:),'LineWidth',3);
hold on
plot(y(1,:),y(2,:),'r','LineWidth',2);
legend('Raw position traj','Detected position traj (have miss detections)')
title('Traj. of Discrete-time NCV Model');
xlabel('traj of x'); ylabel('traj of y');

%%


% F = [1 1 1 1]';
% 
% 
% G = [1;
%      1];


% state noise covariance matrix
V = 0.5;
% measurement noise covariance matrix
W = 0.2;


x_hat(:,1) = [0.5 2 7 4]';
P(:,:,1) = 10*eye(4);

Am = [1 T    0 0; 
     0 0.87 0 0;
     0 0    1 T; 
     0 0    0 0.87];



for k = 1:kmax
    
    if isnan(y(1,k))
%         K = A*P(:,:,k)*C' * inv(C*P(:,:,k)*C' + G*W*G');
        x_hat(:,k+1) = Am*x_hat(:,k);
        P(:,:,k+1) = Am*P(:,:,k)*Am' + F*V*F';
    else
        K = Am*P(:,:,k)*C' / (C*P(:,:,k)*C' + G*W*G');
        x_hat(:,k+1) = Am*x_hat(:,k) + K*(y(:,k) - C*x_hat(:,k));
        P(:,:,k+1) = Am*P(:,:,k)*Am' - K*C*P(:,:,k)*Am' + F*V*F';
        p_trace(:,k) = trace(P(:,:,k));
    end
    
end


figure,
plot(y_raw(1,:),y_raw(2,:),'LineWidth',3);
hold on
plot(x_hat(1,:),x_hat(3,:), 'g','LineWidth',2);
legend('Raw position traj','estiamted traj')

figure,
plot(1:kmax, p_trace)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%
% clear all
% close all
% clc
% 
% maxT=200; % max sim time
% 
% x = zeros(4,maxT); % storage
% 
% % initial position, velocity
% x(:,1) = [0; 
%           0.1; 
%           0; 
%           0.1];
% 
% T = 1; % sample period
% 
% A = [1 T 0 0; 
%      0 1 0 0;
%      0 0 1 T; 
%      0 0 0 1];
% 
% F = [T^2/2 0; 
%      T     0; 
%      0     T^2/2; 
%      0     T];
% 
% G = [1 0;
%      0 1];
% 
% C = [1 0 0 0;
%      0 0 1 0];
% 
% % state noise covariance matrix
% V = [0.85 0;
%      0    0.3];
% % measurement noise covariance matrix
% W = [0.28 0;
%      0    0.65];
% 
% for k = 2:maxT % simulate model
%     x(:,k) = A*x(:,k-1) + F*V*randn(2,1);
%     y(:,k-1) = C*x(:,k-1) + G*W*randn(2,1);
% end
% 
% figure,
% plot(x(1,:),x(3,:));
% title('Discrete-time NCV sim.');
% xlabel('x'); ylabel('y');
% 
% figure,
% plot(y(1,:),y(2,:));
% title('Discrete-time NCV measurements');
% xlabel('x'); ylabel('y');
% 
% 
% x_hat(:,1) = [0.5 20 7 4]';
% P(:,:,1) = 100*eye(4);
% W = 1;
% V = 1;
% 
% 
% % for k = 2:maxT
% % 
% %         K = A*P(:,:,k-1)*C' / (C*P(:,:,k-1)*C' + G*W*G');
% %         x_hat(:,k) = A*x_hat(:,k-1) + K*(y(:,k-1) - C*x_hat(:,k-1));
% %         P(:,:,k) = A*P(:,:,k-1)*A' - K*C*P(:,:,k-1)*A' + F*V*F'; 
% % 
% % end
% 
% Am = [1 T    0 0; 
%      0 0.87 0 0;
%      0 0    1 T; 
%      0 0    0 0.87];
% 
% 
% for k = 1:maxT-1
% 
%         K = Am*P(:,:,k)*C' / (C*P(:,:,k)*C' + G*W*G');
%         x_hat(:,k+1) = Am*x_hat(:,k) + K*(y(:,k) - C*x_hat(:,k));
%         P(:,:,k+1) = Am*P(:,:,k)*Am' - K*C*P(:,:,k)*Am' + F*V*F'; 
% 
% end
% 
% 
% 
% figure,
% plot(x(1,:),x(3,:), 'LineWidth', 2.5)
% hold on
% plot(x_hat(1,:), x_hat(3,:), 'LineWidth',2)
% title('estimated results');
% xlabel('x1'); ylabel('x2');
% legend('states', 'state estimates')

