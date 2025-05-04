clc;close all;clear;
f = @(x) sin(x); % Continuous time system xdot = f(x)
F = @(x,delt) x + delt*f(x); % Discrete time system with step delt, euler integration
g = @(x) x.^2; % Lyapunov-like observable
x = (linspace(-1,1,101))'; % Legendre basis on L2[-1,1]
delt = 0.01; % Time-step 
goF = g(F(x,delt));
phi_1 = @(x) 1/sqrt(2);
phi_2 = @(x) sqrt(3/2)*x;
phi_3 = @(x) sqrt(5/8)*(3*x.^2-1);
phi_4 = @(x) sqrt(7/8)*(5*x.^3-3*x);

K_prime = 15/4*x.*f(x)*(3*x.^2-1)'*(1/(length(x)-1));
Kog = (eye(length(x)) + delt*K_prime)*g(x);
figure;
set(gcf,"WindowState","maximized")
set(gca,"FontSize",18)
hold on;
grid on;
plot(x,goF,"LineWidth",3,"Color",[0 0 1])
plot(x,Kog,"-o","LineWidth",1.5)
legend("$g \circ F(x)$","$K_F \cdot g(x)$","interpreter","latex","FontSize",20,"location","best")
title("$g(x) = x^2$,  system : $\dot{x} = sin(x)$, $\Delta t = 0.01$","interpreter","latex","FontSize",25)
 
% function K = KoopK(x,f,phi_2)
% % Function to find the Koopman operator for a non-linear system with F(x) governing evolution
% % Inputs : x - (n x 1) state 
% %          f - xdot = f(x)
% %          phi - (1 x Nphi) cell of basis functions, sufficient to represent
% %          g(x) exactly
% % Output : K - Koopman operator for f(x), size (n x n)
%     n = length(x)-1;
%     K = 
%     K = K*(1/n); % delta xi in kernel matrix is 1/n here
%     K = real(K);
% end


