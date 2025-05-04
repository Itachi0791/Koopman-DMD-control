clc; clear;
close all
addpath('./Resources')
addpath('./Resources/qpOASES-3.1.0/interfaces/matlab') 

%% Dynamics
% Two Link Manipulator dynamics: xdot = alpha*x - beta*xy + ux ,
%                                ydot = -gamma*y + delta*xy + uy
% State x = [th1; th2; th1dot; th2dot]
% Parameters for the Manipulator 
m1 = 5; m2 = 4;                 % Masses of the links
l1 = 1; l2 = 0.8;               % Lengths of the links
I1 = m1*l1^2/3; I2 = m2*l2^2/3; % Inertia of links
g = 9.81;                       % Gravity

% Inertia matrix
H = @(theta1, theta2) [I1 + m1*l1^2/4 + I2 + m2*(l2^2/4 + l1^2 + l1*l2*cos(theta2)), I2 + m2*(l2^2/4 + l1*l2*cos(theta2)/2); 
                       I2 + m2*(l2^2/4 + l1*l2*cos(theta2)/2), I2 + m2*l2^2/4];

% Coriolis matrix
C = @(theta1, theta2, theta1_dot, theta2_dot) [-m2*l1*l2*sin(theta2)*theta2_dot, -m2*l1*l2*sin(theta2)*theta2_dot/2;
                                               m2*l1*l2*sin(theta2)*theta1_dot/2, 0];

% Gravity vector
G = @(theta1, theta2) [m1*g*l1*cos(theta1)/2 + m2*g*(l1*cos(theta1) + l2/2*cos(theta1 + theta2));
                       m2*g*l2/2*cos(theta1 + theta2)];

% Dynamics function for a single trajectory
f_u_single = @(t, x, u) ([x(3); x(4); ...
                          H(x(1), x(2)) \ (-C(x(1), x(2), x(3), x(4)) * x(3:4) - G(x(1), x(2))) + u]);

n = 4;  % Number of states
m = 2;  % Number of control inputs

%% Discretization
deltaT = 0.01;

% Runge-Kutta 4 for a single trajectory
rk4_single = @(t, x, u) (x + (deltaT / 6) * ...
    (f_u_single(t, x, u) + ...
     2 * f_u_single(t, x + deltaT/2 * f_u_single(t, x, u), u) + ...
     2 * f_u_single(t, x + deltaT/2 * f_u_single(t, x + deltaT/2 * f_u_single(t, x, u), u), u) + ...
     f_u_single(t, x + deltaT * f_u_single(t, x, u), u)));

% Vectorized RK4 for all trajectories
f_ud = @(t, X, U) cell2mat(arrayfun(@(traj) ...
    rk4_single(t, X(:, traj), U(:, traj)), ...
    1:size(X, 2), 'UniformOutput', false));

rng(115123)

%% Data Collection
disp('Starting data collection')
tic;
Ntraj = 1000; Nsim = 4000;

Cy = eye(n); % Output matrix: y = Cy*x
nD = 0; % Number of delays
ny = size(Cy, 1); % Number of outputs

% Random control input forcing between 0 and 1
Ubig = 10*rand(Nsim, m, Ntraj)-5;

% Random initial condition
theta_max = pi; theta_min = -pi; theta_dot_max = 4; theta_dot_min = -4;
Xcurrent = [theta_min + (theta_max - theta_min) * rand(2, Ntraj); ...
            theta_dot_min + (theta_dot_max - theta_dot_min) * rand(2, Ntraj)];

n_zeta = (nD + 1) * ny + nD * m; % Dimension of delay-embedded "state"

% Preallocate memory for performance
X = zeros(n_zeta, Nsim * Ntraj); 
Y = zeros(n_zeta, Nsim * Ntraj); 
U = zeros(m, Nsim * Ntraj); 

index = 1; % Index for storing data
if nD == 0
    % Simplified "state" without delay embedding
    zeta_current = Cy * Xcurrent;
    for i = 1:Nsim
        Xnext = f_ud(0, Xcurrent, squeeze(Ubig(i, :, :)));
        zeta_prev = zeta_current;
        zeta_current = Cy * Xnext;
        X(:, index:index + Ntraj - 1) = zeta_prev;
        Y(:, index:index + Ntraj - 1) = zeta_current;
        U(:, index:index + Ntraj - 1) = squeeze(Ubig(i, :, :));
        index = index + Ntraj;
        Xcurrent = Xnext;
    end
else
    % Delay-embedded "state" zeta_k = [y_{k}; u_{k-1}; y_{k-1}; ... ; u_{k-nd}; y_{k-nd}]
    zeta_current = [Cy * Xcurrent; NaN(nD * (ny + m), Ntraj)];

    for i = 1:Nsim
        Xnext = f_ud(0, Xcurrent, squeeze(Ubig(i, :, :)));
        zeta_prev = zeta_current;
        zeta_current = [[Cy * Xnext; squeeze(Ubig(i, :, :))]; zeta_current(1:end - ny - m, :)];

        % Store data after delay embedding is complete
        if i > nD
            X(:, index:index + Ntraj - 1) = zeta_prev;
            Y(:, index:index + Ntraj - 1) = zeta_current;
            U(:, index:index + Ntraj - 1) = squeeze(Ubig(i, :, :));
            index = index + Ntraj;
        end
        Xcurrent = Xnext;
    end
end

% Trim unused preallocated space
X = X(:, 1:index - 1);
Y = Y(:, 1:index - 1);
U = U(:, 1:index - 1);
save("Resources\data.mat","X","Y","U")
fprintf('Data collection DONE. Time taken - %f s\n', toc);

%% Lift
load('Resources\data.mat')
params = load('Resources\observable_nn_params.mat');
liftFun = @(X) [X;observables(X,params)];
disp('Starting LIFTING')
tic
Xlift = liftFun(X);
Ylift = liftFun(Y);
Nlift = size(Xlift,1);
fprintf('Lifting DONE. Time Taken - %f s\n',toc);

%% Regression

disp('Starting REGRESSION for A,B,C')
tic

W = [Ylift ; X];
V = [Xlift ; U];
VVt = V*V';
WVt = W*V';
ABC = WVt * pinv(VVt);
Alift = ABC(1:Nlift,1:Nlift);
Blift = ABC(1:Nlift,Nlift+1:end);
Clift = ABC(Nlift+1:end,1:Nlift);
fprintf('Regression for A, B, C DONE. Time taken - %f s\n',toc);

% Residual
fprintf( 'Regression residual : %f \n', norm(Ylift - Alift*Xlift - Blift*U,'fro') / norm(Ylift,'fro') );

%% Predictor comparison
close all;
Tmax = 5;
Nsim = Tmax / deltaT;
uprbs = zeros(Nsim,m);%14*rand(Nsim, m)-7; 
u_dt = @(i)(uprbs(i + 1,:));
f_cont_d = @(t, xx)(f_ud(t, xx, u_dt(t)));

x0 = [theta_min + (theta_max - theta_min) * rand(2, 1); ...
      theta_dot_min + (theta_dot_max - theta_dot_min) * rand(2, 1)];
x = x0;

if nD == 0
    % No delay: Initial condition is just Cy*x
    xstart = Cy * x;
    xp = x; % For consistency in defining xloc and x_true
    urand = zeros(m, 1); % Assign a default value to urand for local linearization
else
    % Delayed initial condition (assume random control input in the past)
    xstart = [Cy * x; NaN(nD * (ny + m), 1)];
    for i = 1:nD
        urand = 2 * rand(m, 1) - 1;
        xp = f_ud(0, x, urand);
        xstart = [Cy * xp; urand; xstart(1:end-ny-m)];
        x = xp;
    end
end

% Initial conditions
x_true = xp; % Start true dynamics simulation from xp
xlift = liftFun(xstart);

% Simulation
for i = 0:Nsim-1
    % True dynamics
    x_true = [x_true, f_ud(0, x_true(:, end), u_dt(i))];
    
    % Koopman predictor
    xlift = [xlift, Alift * xlift(:, end) + Blift * u_dt(i)'];
end

th1_true = x_true(1, :); th2_true = x_true(2, :);
th1dot_true = x_true(3, :); th2dot_true = x_true(4, :);
th1_pred = Clift(1, :) * xlift; th2_pred = Clift(2, :) * xlift;
th1dot_pred = Clift(3, :) * xlift; th2dot_pred = Clift(4, :) * xlift;

% figure
% stairs((0:Nsim-1)*deltaT,u_dt(0:Nsim-1),'linewidth',2); hold on
% title('Control input'); xlabel('time [s]')

figure
lw_koop = 2;
subplot(2,2,1);
plot((0:Nsim)*deltaT, th1_true,'-b','linewidth', lw_koop); hold on
plot((0:Nsim)*deltaT, th1_pred, '--r','linewidth',lw_koop)
ylabel('$\theta_1$','interpreter','latex','fontsize',20);
xlabel('time (s)','interpreter','latex','fontsize',14)
title('Koopman Prediction','interpreter','latex','fontsize',20)
LEG = legend('True','Koopman');
set(LEG,'Interpreter','latex','location','northeast','fontsize',14)
set(gca,'FontSize',16);
subplot(2,2,2);
plot((0:Nsim)*deltaT,th2_true,'-b','linewidth', lw_koop); hold on
plot((0:Nsim)*deltaT,th2_pred, '--r','linewidth',lw_koop)
title('Koopman Prediction','interpreter','latex','fontsize',20)
ylabel('$\theta_2$','interpreter','latex','fontsize',20);
xlabel('time (s)','interpreter','latex','fontsize',14)
LEG = legend('True','Koopman');
set(LEG,'Interpreter','latex','location','northeast','fontsize',14)
set(gca,'FontSize',16);
subplot(2,2,3);
plot((0:Nsim)*deltaT, th1dot_true,'-b','linewidth', lw_koop); hold on
plot((0:Nsim)*deltaT, th1dot_pred, '--r','linewidth',lw_koop)
ylabel('$\dot{\theta_1}$','interpreter','latex','fontsize',20);
xlabel('time (s)','interpreter','latex','fontsize',14)
title('Koopman Prediction','interpreter','latex','fontsize',20)
LEG = legend('True','Koopman');
set(LEG,'Interpreter','latex','location','northeast','fontsize',14)
set(gca,'FontSize',16);
subplot(2,2,4);
plot((0:Nsim)*deltaT,th2dot_true,'-b','linewidth', lw_koop); hold on
plot((0:Nsim)*deltaT,th2dot_pred, '--r','linewidth',lw_koop)
ylabel('$\dot{\theta_2}$','interpreter','latex','fontsize',20);
xlabel('time (s)','interpreter','latex','fontsize',14)
title('Koopman Prediction','interpreter','latex','fontsize',20)
LEG = legend('True','Koopman');
set(LEG,'Interpreter','latex','location','northeast','fontsize',14)
set(gca,'FontSize',16);

%% Feedback control 
n_setpoints = 1;       % Number of setpoints
Tmax = 5 * n_setpoints; % Total simulation time
% Initial condition (start near bottom position)
x0 = [0;0;0;0]; 
Nsim = Tmax / deltaT;   % Total number of simulation steps

% Generate random setpoints for theta_ref
theta_refs = [0;0];%2*pi * (rand(2, n_setpoints) - 0.5); % Random values in range [-pi, pi]

% Calculate duration of each setpoint
duration = Nsim / n_setpoints;

% Create yrr
tspan = linspace(0,Tmax,Nsim);
% theta_ref_trajectory = pi*sin(tspan);
% thetadot_ref_trajectory = pi*cos(tspan);
% x0 = [theta_ref_trajectory(1);0];%thetadot_ref_trajectory(1)];
theta_ref_trajectory = repelem(theta_refs, 1, duration); % Repeat each setpoint for 'duration' steps
thetadot_ref_trajectory = zeros(2, Nsim);            % Zero velocity for thetadot

% Combine into yrr
yrr = [theta_ref_trajectory; thetadot_ref_trajectory];

% Define Koopman controller
C = zeros(n,Nlift); 
C(1:n,1:n) = eye(n); % Extract both states from lifted state

% Weight matrices (adjusted for two-state tracking)
Q = diag([0000, 0000, 0, 0]); % Higher weight on angle tracking
QN = diag([1000, 1000, 1, 1]); % Higher terminal weight
R = diag([150,150]);

% Prediction horizon
Tpred = 0.2;
Np = round(Tpred / deltaT);

% Constraints
u_min = [-10;-10]; u_max = [10;10];
theta_min = -pi; theta_max = pi;
thetadot_min = -5; thetadot_max = 5;

xlift_min = [theta_min;theta_min;thetadot_min; thetadot_min; -inf*ones(Nlift-n,1)];
xlift_max = [theta_max;theta_max;thetadot_max; thetadot_max; inf*ones(Nlift-n,1)];

% Build Koopman MPC controller
koopmanMPC = getMPC(Alift,Blift,C,zeros(n,1),Q,R,Q,Np,u_min,u_max,xlift_min,xlift_max,'qpoases');

% Simulation initialization
x_koop = x0;
XX_koop = x0; 
UU_koop = [];
zeta0 = Cy*x_koop;
% Initial condition for the delay-embedded state (assuming zero control in the past)
if nD ~= 0
    zeta0 = [Cy*x_koop ; NaN(nD*(ny+m),1)];
    for i = 1:nD
        upast = zeros(m,1);
        xp = f_ud(0,x_koop,upast);
        zeta0 = [Cy*xp ; upast ; zeta0(1:n_zeta-ny-m)]; 
        x_koop = xp;
    end
    x0 = x_koop;
end

x_koop = x0; x_loc = x0;
zeta = zeta0; % Delay-embedded "state"

% Simulation loop
for i = 0:Nsim-1
    if(mod(i,round(Nsim/5)) == 0)
        fprintf('Closed-loop simulation: iterate %i out of %i \n', i+round(Nsim/5), Nsim)
    end
    
    % Current value of the reference signal
    yr = yrr(:,i+1);

    % Koopman MPC
    xlift = liftFun(zeta); % Lift
    u_koop = koopmanMPC(xlift,yr); % Get control input
    x_koop = f_ud(0,x_koop,u_koop); % Update true state
    
    % Lift current state
    if nD == 0
        zeta = x_koop; % No delay embedding
    else
        zeta = [Cy*x_koop; u_koop; zeta(1:n_zeta-ny-m)]; % With delay embedding
    end
    xlift = liftFun(zeta);

    % Store values
    XX_koop = [XX_koop x_koop];
    UU_koop = [UU_koop u_koop];
end

%% Plot results
figure('Position', [100, 100, 800, 600])
time = (0:Nsim)*deltaT;
th1 = XX_koop(1,:); th2 =  XX_koop(2,:);
thdot1 = XX_koop(1,:); thdot2 =  XX_koop(2,:);
% Plot theta
subplot(3,2,1)
plot(time, th1, 'LineWidth', 3)
hold on
stairs(time(1:end-1), yrr(1,:), 'LineWidth', 2)
plot(time, theta_max*ones(1,Nsim+1), '--k', 'LineWidth', 2)
plot(time, theta_min*ones(1,Nsim+1), '--k', 'LineWidth', 2)
ylabel('$\theta_1$','Interpreter','latex','FontSize',17,'FontWeight','bold')
legend('Actual', 'Reference', 'Constraints','FontSize',12)
grid on
set(gca,'fontsize',14)
%title('Pendulum KMPC','FontSize',17)

subplot(3,2,2)
plot(time, th2, 'LineWidth', 3)
hold on
stairs(time(1:end-1), yrr(2,:), 'LineWidth', 2)
plot(time, theta_max*ones(1,Nsim+1), '--k', 'LineWidth', 2)
plot(time, theta_min*ones(1,Nsim+1), '--k', 'LineWidth', 2)
ylabel('$\theta_2$','Interpreter','latex','FontSize',17,'FontWeight','bold')
legend('Actual', 'Reference', 'Constraints','FontSize',12)
set(gca,'fontsize',14)
grid on
subplot(3,2,3)
plot(time, thdot1, 'LineWidth', 3)
hold on
stairs(time(1:end-1), yrr(3,:), 'LineWidth', 2)
plot(time, thetadot_max*ones(1,Nsim+1), '--k', 'LineWidth', 2)
plot(time, thetadot_min*ones(1,Nsim+1), '--k', 'LineWidth', 2)
ylabel('$\dot{\theta_1}$','Interpreter','latex','FontSize',17,'FontWeight','bold')
legend('Actual', 'Reference', 'Constraints','FontSize',12)
grid on
set(gca,'fontsize',14)
%title('Pendulum KMPC','FontSize',17)

subplot(3,2,4)
plot(time, thdot2, 'LineWidth', 3)
hold on
stairs(time(1:end-1), yrr(4,:), 'LineWidth', 2)
plot(time, thetadot_max*ones(1,Nsim+1), '--k', 'LineWidth', 2)
plot(time, thetadot_min*ones(1,Nsim+1), '--k', 'LineWidth', 2)
ylabel('$\dot{\theta_2}$','Interpreter','latex','FontSize',17,'FontWeight','bold')
legend('Actual', 'Reference', 'Constraints','FontSize',12)
set(gca,'fontsize',14)
grid on
% Plot control input
subplot(3,2,5)
plot(time(1:end-1), u_max*ones(1,Nsim), '--k', 'LineWidth', 2)
hold on
plot(time(1:end-1), u_min*ones(1,Nsim), '--k', 'LineWidth', 2)
plot(time(1:end-1), UU_koop(1,:), 'LineWidth', 3,'Color','r')
ylabel('$u_1$','Interpreter','latex','FontSize',17,'FontWeight','bold')
xlabel('Time (s)','FontSize',17,'FontWeight','bold')
legend('Constraints','FontSize',12)
set(gca,'fontsize',14)
grid on
subplot(3,2,6)
plot(time(1:end-1), u_max*ones(1,Nsim), '--k', 'LineWidth', 2)
hold on
plot(time(1:end-1), u_min*ones(1,Nsim), '--k', 'LineWidth', 2)
plot(time(1:end-1), UU_koop(2,:), 'LineWidth', 3,'Color','r')
ylabel('$u_2$','Interpreter','latex','FontSize',17,'FontWeight','bold')
xlabel('Time (s)','FontSize',17,'FontWeight','bold')
legend('Constraints','FontSize',12)
set(gca,'fontsize',14)
grid on

%% Functions
params = load('Resources\observable_nn_params.mat');

% Define your input X (n x N matrix)
X = randn(4, 10);  % Example input matrix with 4 features and 10 samples

% Compute the observables (i.e., the feedforward output)
y = observables(X, params);

% Display the output
disp(y);
function y = observables(X, weights_biases)
    % This function performs a feedforward pass through a neural network
    % X is the input matrix (n x N)
    % weights_biases is the structure containing the weights and biases
    
    % Get the field names in the weights_biases structure
    weight_fields = fieldnames(weights_biases);
    
    % Start the feedforward pass
    layer_input = X;
    % Iterate through each layer based on the field names
    for i = 0:2:(length(weight_fields))-1
        % Extract the weights and biases for the current layer
        W = weights_biases.(sprintf('layer%d_weights', i));  % Access layeri_weights
        b = weights_biases.(sprintf('layer%d_biases', i));   % Access layeri_biases
        
        % Apply the linear transformation: W * input + b
        layer_input = W * layer_input + b';
        
        % Apply the activation function (Tanh)
        layer_input = selu(layer_input);
    end
        
    % The output is the result of the last layer transformation
    y = layer_input;
end

function y = selu(x)
    % SELU activation function
    % x is the input (matrix, vector, or scalar)
    
    % Define the constants for SELU
    lambda = 1.0507;   % Scaling factor
    alpha = 1.6733;    % Stretching factor
    
    % Apply the SELU activation function element-wise
    y = lambda * (x .* (x > 0) + alpha * (exp(x) - 1) .* (x <= 0));
end

