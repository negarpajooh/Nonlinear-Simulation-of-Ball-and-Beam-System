function dx = nonlinear_model1(t,x)
global m R M g L Jr J x0 vars A B C D k1 k2 Tau
m=0.11 ; % 	mass of the ball   0.11 kg 	 
R=0.015; 	%radius of the ball 	0.015 m 
M=1;	%Mass of the beam	1 kg
g=9.8; 	%gravitational acceleration 	9.8 m/s^2 
L=1.0; 	%length of the beam 	1.0 m 
Jr=1e-5 ;	%ball's moment of inertia 	9.99e-6 kgm^2 
J=2e-3; 	%beam's moment of inertia 
%-------------------------------------------------------------
dx = zeros(4,1);    % a column vector
%u=-Kx
dx(1) = x(2) ;
dx(2) =[(J+m*x(1)^2+Jr)*(m*x(1)*x(4)^2-m*g*sin(x(3)))+Jr/R*...
    (Tau-2*m*x(1)...
    *x(2)*x(4)-m*g*x(1)*cos(x(3)))]/(m*(J+Jr)+m*(m+Jr/R^2)*x(1)^2+J*Jr/R^2);
dx(3) = x(4);
dx(4) = [Jr/R^2*(m*x(1)*x(4)^2-m*g*sin(x(3)))+(m+Jr/R^2)*...
    (Tau-2*m*x(1)...
    *x(2)*x(4)-m*g*x(1)*cos(x(3)))]/(m*(J+Jr)+m*(m+Jr/R^2)*x(1)^2+J*Jr/R^2);
   