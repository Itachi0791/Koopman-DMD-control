clc;%clear;%close all;
f = @(x) -tanh(x); % Continuous time system xdot = f(x)
f_prime = @(x) -3*x.^2;
F = @(x,delt) x + delt*f(x)+(delt)^2/2*f_prime(x); % Euler integration till 2nd term
g = @(x) x.^2; % Lyapunov-like observable
% % Legendre polynomials
num_basis = 10; phi = cell(1,num_basis);
% Define basis functions as a cell array
for i = 1:num_basis
    phi{i} = @(x)  norm_legendre(i-1, x);
end
x = (linspace(-1,1,100))'; % Legendre basis on L2[-1,1]
delt = 1e-2; % Time-step 
goF = g(F(x,delt));
Kog = g(x) + 2*delt*x.*f(x);
% figure;
% set(gcf,"WindowState","maximized")
% set(gca,"FontSize",18)
% hold on;
% grid on;
% plot(x,goF,"LineWidth",3,"Color",[0 0 1])
% plot(x,Kog,"-o","LineWidth",1.5)
% legend("$g \circ F(x)$","$K_F \cdot g(x)$","interpreter","latex","FontSize",20,"location","best")
% title("$g(x) = x^2$,  system : $\dot{x} = sin(x)$, $\Delta t = 0.01$","interpreter","latex","FontSize",25)

%% Kf formulation
x = (linspace(-1,1,10001))';
Kf_prime = ( 3/2*f(x)*(x)' + 15/4*(x.*f(x))*(3*x.^2-1)' +...
           7/8*((15*x.^2-3).*f(x))*(5*x.^3-3*x)' )/(length(x)-1);
Kf = eye(length(x)) + delt*Kf_prime;
% Kf2 = (phi_1(F(x,delt))*ones(size(x))*(phi_1(x))' + phi_2(F(x,delt))*(phi_2(x))'...
%       + phi_3(F(x,delt))*(phi_3(x))' + phi_4(F(x,delt))*(phi_4(x))')/(length(x)-1);
X = [phi{1}(x) phi{2}(x) phi{3}(x) phi{4}(x)];
A_bar_prime = (X')*(Kf')*pinv(X');
%% A bar formulation
f_coeffs = zeros(num_basis,1);
for i=1:num_basis
    f_coeffs(i) = integral(@(x) f(x).*norm_legendre(i-1,x),-1,1);
end
% Initialize A_bar
A_bar = zeros(num_basis, num_basis); A_bar_est = A_bar;
% Compute A_bar elements using numerical integration
for i = 0:num_basis-1
    for j = 0:num_basis-1
        % % Define the integrand
        % integrand = @(x) phi{i+1}(F(x,delt)) .* phi{j+1}(x);
        % % Compute the integral using MATLAB's integral function
        % A_bar(i+1, j+1) = integral(integrand, -1, 1, 'ArrayValued', true);
        sum1 = 0; sum2 = 0;
        for k = 0:num_basis-1
            for q = 1 : floor((i+1)/2)
                sum1 = sum1 + f_coeffs(k+1)*sqrt(2*k+1)*(2*i-4*q+3)*Wigner3j([i-2*q+1,k,j],[0,0,0]).^2;
                for r = 1 : floor((k+1)/2)
                    sum2 = sum2 + f_coeffs(k+1)*sqrt(2*k+1)*(2*i-4*q+3)*(2*k-4*r+3)*Wigner3j([i-2*q+1,k-2*r+1,j],[0,0,0]).^2;
                end
            end
        end
        A_bar_est(i+1,j+1) = (i==j) + sqrt((2*i+1)*(2*j+1)/2)*(delt*sum1 + delt^2/2*sum2);
    end
end
Ns = 200;
%% Prediction from Simulation
%rng(1)
x0 = 2*rand-1; phi_x = zeros(num_basis,Ns);
for i=1:num_basis
    phi_x(i,1) = phi{i}(x0);
end
phi_x_prime = phi_x;
[t, x] = ode45(@(t,x) f(x), 0:delt:(Ns-1)*delt, x0);
for i=1:Ns-1
    phi_x(:,i+1) = A_bar_est*phi_x(:,i);
   % phi_x_prime(:,i+1) = A_bar_prime*phi_x_prime(:,i);
end
C_x = zeros(1,num_basis); C_x(2) = sqrt(2/3); % Mapping from phi_1,phi_2,phi_3 to x
%C_g = [sqrt(2)/3, 0, 2/3*sqrt(2/5), 0]; % Mapping from phi_1,phi_2,phi_3 to g
x_pred = C_x*phi_x; %x_pred_prime = C_x*phi_x_prime;
%g_x = C_g*phi_x; %g_x_prime = C_g*phi_x_prime;
figure;  
set(gcf,"WindowState","maximized")
set(gca,"FontSize",18)
hold on;
grid on;
plot(0:delt:(Ns-1)*delt,x,"LineWidth",3.5,"Color",[0 0 1])
plot(0:delt:(Ns-1)*delt,x_pred,"-o","LineWidth",1.5,"Color",[1 0 0])
%plot(0:delt:(Ns-1)*delt,x_pred_prime,"-o","LineWidth",1.5,"Color",[0 1 0])
xlabel("Time(t)","FontSize",20,"interpreter","latex", "FontWeight", "bold");
ylabel("$x(t)$","FontSize",20,"interpreter","latex", "FontWeight", "bold")
title("Time evolution prediction","interpreter","latex","FontSize",25)
legend("Actual","Predicted (using $\bar{A}$)","interpreter","latex")
%legend("Actual","Predicted (using $\bar{A}$)","Predicted (using $\bar{A}'$)","interpreter","latex")
% figure;
% set(gcf,"WindowState","maximized")
% set(gca,"FontSize",18)
% hold on;
% grid on;
% plot(0:delt:(Ns-1)*delt,g(x),"LineWidth",2.5,"Color",[0 0 1])
% plot(0:delt:(Ns-1)*delt,g_x,"-o","LineWidth",1.5,"Color",[1 0 0])
% plot(0:delt:(Ns-1)*delt,g_x_prime,"-o","LineWidth",1.5,"Color",[0 1 0])
% xlabel("Time(t)","FontSize",20,"interpreter","latex", "FontWeight", "bold");
% ylabel("$g(x(t))$","FontSize",20,"interpreter","latex", "FontWeight", "bold")
% title("Time evolution of Lyapunov-like observable, $g(x)=x^2$","interpreter","latex","FontSize",25)
% legend("Actual","Predicted (using $\bar{A}$)","Predicted (using $\bar{A}'$)","interpreter","latex")

%%
% clc; clear; close all;
% 
% % Define functions
% f = @(x) sin(x); % Continuous time system xdot = f(x)
% F = @(x,delt) x + delt*f(x); % Discrete time system with step delt, Euler integration
% g = @(x) x.^2; % Lyapunov-like observable
% x = linspace(-1, 1, 100); % Legendre basis on L2[-1,1]
% 
% % Array of delt values
% delt_values = [0.01, 0.05, 0.1, 0.2]; % Specify step sizes
% 
% % Initialize figure
% figure;
% set(gcf, "WindowState", "maximized");
% set(gca, "FontSize", 18);
% hold on;
% grid on;
% 
% % Plot g∘F(x) once, as it doesn't depend on delt
% goF = g(F(x, delt_values(1))); % Just using the first delt for F(x)
% plot(x, goF, "LineWidth", 4, "Color", [0, 0, 1], 'DisplayName', "$g \circ F(x)$");
% 
% % Compute and plot Kog for each delt
% %colors = lines(length(delt_values)); % Generate distinct colors
% for i = 1:length(delt_values)
%     delt = delt_values(i);
%     Kog = g(x) + 2 * delt * x .* f(x); % Compute Kog
%     plot(x, Kog, "-o", "LineWidth", 1, ...
%         'DisplayName', sprintf("$K_{F} \\cdot g(x)$   $(\\Delta t = %.2f)$", delt));
% end
% 
% % Add legend and title
% legend("interpreter", "latex", "FontSize", 18, "Location", "best");
% title("$g(x) = x^2$,  system : $\dot{x} = sin(x)$", ...
%     "interpreter", "latex", "FontSize", 25);

%%
function val = in_prod(k,l,j)
    val = sqrt((2*j+1)*(2*k+1)*(2*l+1)/2)*(Wigner3j([k,l,j],[0 0 0]))^2;
end



