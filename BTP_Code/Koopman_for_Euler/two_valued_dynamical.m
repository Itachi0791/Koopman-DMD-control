clc;clear;close all;
f = @(x) [x(2,:);-sin(x(1,:))]; % Continuous time system xdot = f(x)
F = @(x,delt) x + delt*f(x); % Discrete time system with step delt, euler integration
g = @(x) x(1,:).^2 + x(2,:).^2; % Lyapunov-like observable
x1 = linspace(-1, 1, 100); % Linearly spaced values for x1
x2 = linspace(-1, 1, 100); % Linearly spaced values for x2
[X1, X2] = meshgrid(x1, x2); % Create a 2D grid
x = [X1(:)'; X2(:)']; % Legendre basis on L2[-1,1]
delt = 0.01; % Time-step 
goF = g(F(x,delt));
Kog = g(x) + 4*delt*dot(x,f(x));
% Reshape results back to 2D grid for surface plotting
goF_grid = reshape(goF, size(X1));
Kog_grid = reshape(Kog, size(X1));

% Overlay the surfaces
figure;
set(gcf,"WindowState","maximized")
set(gca,"FontSize",18)
surf(X1, X2, goF_grid, 'FaceAlpha', 0.6, 'EdgeColor', 'none'); % goF surface
hold on;
plot3(X1(:), X2(:), Kog_grid(:), 'ro', 'MarkerSize', 4); % Kog surface

% Customize the plot
title("$g(x) = x_1^2+x_2^2$,  Pendulum system, $\Delta t = 0.01$","interpreter","latex","FontSize",25)
xlabel('$x_1$',"interpreter","latex","FontSize",20); ylabel('$x_2$',"interpreter","latex","FontSize",20);
legend({'$g \circ F(x)$', '$K_F \cdot g(x)$'},"interpreter","latex","FontSize",20, 'Location', 'northwest');
view(3); % Ensure a 3D view
axis equal;

%% A_bar formulation
% Basis functions in 2D, x1[-pi to pi], x2[-1 to 1]
% phi_1 = @(x1, x2) 1 / sqrt(4 * pi)*ones(size(x1));
% phi_2 = @(x1, x2) sin(x1) / sqrt(pi);
% phi_3 = @(x1, x2) sqrt(3 / (4 * pi)) * x2;
% phi_4 = @(x1, x2) cos(x1) / sqrt(pi);
% phi_5 = @(x1, x2) sqrt(3 / (2 * pi)) * x2 .* sin(x1);
% phi_6 = @(x1, x2) sqrt(3 / (2 * pi)) * x2 .* cos(x1);

% Define basis functions as a cell array
%phi = {phi_1, phi_2, phi_3, phi_4, phi_5, phi_6};

% % Basis functions in 2D, x1[-pi to pi], x2[-1 to 1]
phi_1 = @(x1, x2) 1/2*ones(size(x1));
phi_2 = @(x1, x2) sqrt(3)/2*x1;
phi_3 = @(x1, x2) sqrt(3)/2*x2;
phi_4 = @(x1, x2) sqrt(5)/2*(3*x1.^2-1);
phi_5 = @(x1, x2) sqrt(5)/2*(3*x2.^2-1);

% Define basis functions as a cell array
phi = {phi_1, phi_2, phi_3, phi_4, phi_5};

% Discrete-time system
F1 = @(x1, x2, delt) x1 + delt * x2;
F2 = @(x1, x2, delt) x2 - delt * sin(x1); 

% Initialize A_bar
num_basis = length(phi);
A_bar = zeros(num_basis, num_basis);

% Compute A_bar elements using numerical integration
for i = 1:num_basis
    for j = 1:num_basis
        % Define the integrand
        integrand = @(x1, x2) ...
            phi{i}(F1(x1, x2, delt), F2(x1, x2, delt)) .* phi{j}(x1, x2);

        % Compute the integral over the 2D domain
        A_bar(i, j) = integral2(integrand, -pi, pi, -1, 1);
    end
end

Ns = 200;
x0 = rand(2,1)-0.5; phi_x = zeros(num_basis,Ns);
for i=1:num_basis
    phi_x(i,1) = phi{i}(x0(1),x0(2));
end
[t, x] = ode45(@(t,x) f(x), 0:delt:(Ns-1)*delt, x0);
x = x';
for i=1:Ns-1
    phi_x(:,i+1) = A_bar*phi_x(:,i);
end
C_g = [sqrt(2)/3, 0, 2/3*sqrt(2/5), 0]; % Mapping from phi_1,phi_2,phi_3 to g
x2_pred = [0,0,2*sqrt(pi/3),0,0,0]*phi_x;
x1_pred = atan2([0,1,0,0,0,0]*phi_x,[0,0,0,1,0,0]*phi_x);
% x1_pred = [0,2*sqrt(3),0,0,0]*phi_x;
% x2_pred = [0,0,2*sqrt(3),0,0]*phi_x;
figure;
set(gcf,"WindowState","maximized")
sgtitle("Pendulum Dynamical Prediction using $\bar{A}$","FontSize",25,"interpreter","latex")
subplot(2,1,1)
hold on
grid on
set(gca,"FontSize",20)
plot(0:delt:(Ns-1)*delt,x(1,:),'LineWidth',2)
plot(0:delt:(Ns-1)*delt,x1_pred,'-o','LineWidth',1)
xlabel("Time (t)","FontSize",20);
ylabel("Angle ($\theta$)","FontSize",20,"Interpreter","latex")
legend("Actual","Predicted","FontSize",18)
subplot(2,1,2)
hold on
grid on
set(gca,"FontSize",20)
plot(0:delt:(Ns-1)*delt,x(2,:),'LineWidth',2)
plot(0:delt:(Ns-1)*delt,x2_pred,'-o','LineWidth',1)
xlabel("Time (t)","FontSize",20);
ylabel("Anglular Velocity ($\dot{\theta}$)","FontSize",20,"Interpreter","latex")
legend("Actual","Predicted","FontSize",18)





