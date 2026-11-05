clear;
clc;
close all;
global m R M g L Jr J x0 vars A B C D k1 k2 klqr Tau
%% Defaults for this blog post
width = 8;     % Width in inches
height = 4;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz =16 ;      % Fontsize
lw = 1;      % LineWidth
msz = 8;       % MarkerSize
%%
m=0.11 ; % 	mass of the ball   0.11 kg 	 
R=0.015; 	%radius of the ball 	0.015 m 
M=1;	%Mass of the beam	1 kg
g=9.8; 	%gravitational acceleration 	9.8 m/s^2 
L=1.0; 	%length of the beam 	1.0 m 
Jr=1e-5 ;	%ball's moment of inertia 	9.99e-6 kgm^2 
J=2e-3; 	%beam's moment of inertia 	
r=0; 	%ball position coordinate 	
rdot=0;
theta=0; 	%beam angle coordinate [0.5 0.2 1 0.3]
thetadot=0;
% x0=[r;rdot;theta;thetadot];
x0=[0.5 0.2 1 0.3];
% k1=0;
% k2=10;
% vars=[m R M g L Jr J x0];
%%
A=[0 1 0 0
  (-(Jr/R)*m*g)/(m*(J+Jr)+(J*Jr)/(R^2)) 0 -((J+Jr)*m*g)/(m*(J+Jr)+((J*Jr)/(R^2))) 0
   0 0 0 1
   -((m+(Jr/(R^2)))*m*g)/(m*(J+Jr)+(J*Jr)/(R^2)) 0 -((Jr/R)*m*g)/(m*(J+Jr)+((J*Jr)/(R^2))) 0];
B=[0
   (Jr/R)/(m*(J+Jr)+((J*Jr)/(R^2)))
   0
   (m+((Jr)/(R^2)))/(m*(J+Jr)+((J*Jr)/(R^2)))];
C=eye(4);
D=zeros(4,1);
%%
% sys=ss(A,B,C,D);
% t=0:0.1:5;
% u=1.*ones(1,numel(t));
% [y,~,~]=lsim(sys,u,t,x0);
% plot(t,y(:,3),'linewidth',2)
% C=eye(4);
% D=zeros(4,1);
% C=[0 1 0 1];
% C=[1 0 0 0];
% C=[0 0 1 0];
% C=[1 0 0 1];
% C=[0 1 1 0];
%% calculation of transition matrix
for t=0:0.1:1
M1=expm(A*t);
i=int8(t*10+1);
y(1:4,i)=M1*[0.2 0.05 1 0.01]';
end
% figure(1)
% plot((0:.1:1),y(1,:),'linewidth',2)
% 
% grid on
%% ctrb
Co = ctrb(A,B)
Rank=rank(obsv(A,C))
k = acker(A,B,[-2+2j,-2-2j,-7,-7])
Eigen=eig(A-B*k)
%% LQR
sys1=ss(A,B,C,D);
Q=diag([1 0 80/180*pi 0]);
R=1;N=0;
Tau=1;
[T1,X1] = ode45(@nonlinear_model1,[0 5],[0.5 0.2 1 0.3]);
klqr=lqr(sys1,Q,R,N)
[T2,X2] = ode45(@nonlinear_model,[0 5],[0.5 0.2 1 0.3]);
figure(2)
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*300]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties

subplot(4,1,1)
plot(T1,X1(:,1),'linewidth',2)
hold on
plot(T2,X2(:,1),'linewidth',2)
legend('r','r_c_t_r_l');
grid on
ylabel('r - m')
xlabel('Time- sec')
subplot(4,1,2)
plot(T2,X2(:,1),'linewidth',2)
legend('r_c_t_r_l');
grid on
ylabel('r - m')
xlabel('Time- sec')
subplot(4,1,3)
plot(T1,X1(:,3),'linewidth',2)
hold on
plot(T2,X2(:,3),'linewidth',2)
grid on
ylabel('\theta - rad')
xlabel('Time- sec')
legend('teta','teta_c_t_r_l');
subplot(4,1,4)
plot(T2,X2(:,3),'linewidth',2)
ylabel('\theta - rad')
xlabel('Time- sec')
grid on
legend('teta_c_t_r_l');
figure(3)
plot(T2,-klqr*[X2(:,1) X2(:,2) X2(:,3) X2(:,4)]','linewidth',2)
grid on
legend('u');
ylabel('Torque - N.m')
xlabel('Time- sec')
%%
sim('BallBeam1')