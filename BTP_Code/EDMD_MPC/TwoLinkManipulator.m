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
Ntraj = 1; Nsim = 10;

Cy = eye(n); % Output matrix: y = Cy*x
nD = 0; % Number of delays
ny = size(Cy, 1); % Number of outputs

% Random control input forcing between 0 and 1
Ubig = zeros(Nsim,m,Ntraj);%10*rand(Nsim, m, Ntraj)-5;

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

fprintf('Data collection DONE. Time taken - %f s\n', toc);

%% Basis functions
basisFunction = 'rbf';
Nrbf = 100;
%cent = kmeansplusplus(X, Nrbf);
cent = rand(n_zeta,Nrbf)*2 - 1; % RBF centers
rbf_type = 'thinplate';
liftFun = @(xx)( [xx;rbf(xx,cent,rbf_type)] );
Nlift = Nrbf + n_zeta;
%type = 'hermite';
%liftFun = @(xx)( [xx;lifting_functions(xx,type,4,true)] );
%% Modified Basis functions for Two-Link Manipulator
% Function that computes the custom observables for a single state
% State vector structure: [theta1; theta2; theta1_dot; theta2_dot]
% Define the custom observables function
custom_observables = @(theta1, theta2, theta1_dot, theta2_dot)([
    % Single terms
    sin(theta1);
    cos(theta1);
    sin(theta2);
    cos(theta2);

    % Degree 2 combinations
    sin(theta1).^2;
    cos(theta1).^2;
    sin(theta2).^2;
    cos(theta2).^2;
    theta1_dot.^2;
    theta2_dot.^2;
    sin(theta1).*cos(theta1);
    sin(theta2).*cos(theta2);
    sin(theta1).*sin(theta2);
    sin(theta1).*cos(theta2);
    cos(theta1).*sin(theta2);
    cos(theta1).*cos(theta2);
    sin(theta1).*theta1_dot;
    cos(theta1).*theta1_dot;
    sin(theta2).*theta1_dot;
    cos(theta2).*theta1_dot;
    sin(theta1).*theta2_dot;
    cos(theta1).*theta2_dot;
    sin(theta2).*theta2_dot;
    cos(theta2).*theta2_dot;
    theta1_dot.*theta2_dot;

    % Degree 3 combinations
    sin(theta1).^3;
    cos(theta1).^3;
    sin(theta2).^3;
    cos(theta2).^3;
    theta1_dot.^3;
    theta2_dot.^3;
    sin(theta1).^2.*cos(theta1);
    sin(theta1).*cos(theta1).^2;
    sin(theta2).^2.*cos(theta2);
    sin(theta2).*cos(theta2).^2;
    sin(theta1).*sin(theta2).*cos(theta2);
    cos(theta1).*sin(theta2).*cos(theta2);
    sin(theta1).*sin(theta2).*cos(theta1);
    cos(theta1).*cos(theta2).*sin(theta2);
    sin(theta1).*theta1_dot.^2;
    cos(theta1).*theta1_dot.^2;
    sin(theta2).*theta1_dot.^2;
    cos(theta2).*theta1_dot.^2;
    sin(theta1).*theta2_dot.^2;
    cos(theta1).*theta2_dot.^2;
    sin(theta2).*theta2_dot.^2;
    cos(theta2).*theta2_dot.^2;
    theta1_dot.^2.*theta2_dot;
    theta1_dot.*theta2_dot.^2;
    sin(theta1).*theta1_dot.*theta2_dot;
    cos(theta1).*theta1_dot.*theta2_dot;
    sin(theta2).*theta1_dot.*theta2_dot;
    cos(theta2).*theta1_dot.*theta2_dot;

    % Degree 4 combinations
    sin(theta1).^4;
    cos(theta1).^4;
    sin(theta2).^4;
    cos(theta2).^4;
    theta1_dot.^4;
    theta2_dot.^4;
    sin(theta1).^3.*cos(theta1);
    sin(theta1).*cos(theta1).^3;
    sin(theta2).^3.*cos(theta2);
    sin(theta2).*cos(theta2).^3;
    sin(theta1).^2.*sin(theta2).^2;
    cos(theta1).^2.*cos(theta2).^2;
    sin(theta1).^2.*cos(theta2).^2;
    cos(theta1).^2.*sin(theta2).^2;
    sin(theta1).*cos(theta1).*sin(theta2).*cos(theta2);
    sin(theta1).*sin(theta2).*theta1_dot.^2;
    cos(theta1).*cos(theta2).*theta1_dot.^2;
    sin(theta1).*cos(theta2).*theta1_dot.^2;
    cos(theta1).*sin(theta2).*theta1_dot.^2;
    sin(theta1).*sin(theta2).*theta2_dot.^2;
    cos(theta1).*cos(theta2).*theta2_dot.^2;
    sin(theta1).*cos(theta2).*theta2_dot.^2;
    cos(theta1).*sin(theta2).*theta2_dot.^2;
    sin(theta1).*theta1_dot.*theta2_dot.^2;
    cos(theta1).*theta1_dot.*theta2_dot.^2;
    sin(theta2).*theta1_dot.*theta2_dot.^2;
    cos(theta2).*theta1_dot.*theta2_dot.^2;
    sin(theta1).^2.*theta1_dot.^2;
    cos(theta1).^2.*theta1_dot.^2;
    sin(theta2).^2.*theta2_dot.^2;
    cos(theta2).^2.*theta2_dot.^2;
    theta1_dot.^2.*theta2_dot.^2;
    sin(theta1).*cos(theta1).*theta1_dot.*theta2_dot;
    sin(theta2).*cos(theta2).*theta1_dot.*theta2_dot;
    sin(theta1).*sin(theta2).*theta1_dot.*theta2_dot;
    cos(theta1).*cos(theta2).*theta1_dot.*theta2_dot;
    sin(theta1).*theta1_dot.^3;
    cos(theta1).*theta1_dot.^3;
    sin(theta2).*theta2_dot.^3;
    cos(theta2).*theta2_dot.^3;
]);



% Create the lifting function for the delay-embedded state
liftFun = @(xx)( lift_with_delays_and_observables(xx, nD, ny, m, custom_observables) );

% Calculate new lifting dimension
% Original dimension + 16 observables for each delay
%Nlift = n_zeta + 16*(nD + 1); 

function phi = lift_with_delays_and_observables(zeta, nD, ny, m, obs_fun)
    % zeta structure: [y_k; u_{k-1}; y_{k-1}; ...; u_{k-nD}; y_{k-nD}]
    % where y_k = [theta1; theta2; theta1_dot; theta2_dot]

    % Initialize the lifted state vector with the original delay-embedded data
    phi = zeta;

    % Get the number of samples
    n_samples = size(zeta, 2);

    % Extract and process each delayed state component
    for i = 0:nD
        % Calculate position in zeta for current delay
        if i == 0
            idx = 1:ny; % Current state (no delay)
        else
            idx = (ny + m)*i + 1:(ny + m)*i + ny; % Delayed states
        end

        % Extract states for current delay
        theta1 = zeta(idx(1), :);     % First joint angle
        theta2 = zeta(idx(2), :);     % Second joint angle
        theta1_dot = zeta(idx(3), :); % First joint velocity
        theta2_dot = zeta(idx(4), :); % Second joint velocity

        % Compute observables for current delay
        obs = obs_fun(theta1, theta2, theta1_dot, theta2_dot);

        % Append observables to phi
        phi = [phi; obs];
    end
end

%% Lift
disp('Starting LIFTING')
tic
Xlift = liftFun(X);
Ylift = liftFun(Y);
Nlift = size(Xlift,1);
save('koopman_data.mat', 'Xlift', 'Ylift', 'U');
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
Y_pred = Alift*Xlift + Blift*U;
P = Ylift*pinv(Y_pred);
Alift = P*Alift; Blift = P*Blift;
fprintf('Regression for A, B, C DONE. Time taken - %f s\n',toc);

% Residual
fprintf( 'Regression residual : %f \n', norm(Ylift - Alift*Xlift - Blift*U,'fro') / norm(Ylift,'fro') );

%% Multi Step Koopman
% Tp = 10; % No. of steps 
% [A_opt, B_opt] = Multi_Step_Koopman(Xlift, U, Ylift, Ntraj, Nsim, Tp, Nlift, m);
%load("DataFiles\koopman_matrices.mat")
%Alift = double(A_lift); Blift = double(B_lift); Clift = [eye(n),zeros(n,Nlift-n)];
%% Predictor comparison
close all;
Tmax = 0.05;
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
Tpred = 0.1;
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

%% Helper Functions

function [loss] = multi_step_loss(A_lift, B_lift, X_lift, U, Y_lift, N_sim, N_traj, Tp)    
    % Compute multi-step loss
    loss = 0;
    k = 1; % Dummy Counter
    while k <= N_sim * N_traj - Tp
         for i = k : (k + N_traj - Tp - 1)
            Z = X_lift(:,i);
            loss = loss + norm(Y_lift(:,i) - A_lift*Z-B_lift*U(:,i), 'fro')^2;
            for j = 0 : Tp - 1
                Z = A_lift*Z + B_lift*U(:,i+j);
                loss = loss + norm(Y_lift(:,i+j+1) - Z, 'fro')^2;      
            end
         end
         k = k + N_traj;
    end
end

function [A_opt, B_opt] = Multi_Step_Koopman(X_lift, U, Y_lift, N_sim, N_traj, Tp, N_lift, m)
    % Define initial guesses for A_lift and B_lift
    A_init = randn(N_lift, N_lift);  % Random initialization
    B_init = randn(N_lift, m);

    % Flatten A_lift and B_lift into a single vector
    theta_init = [A_init(:); B_init(:)];

    % Define loss function wrapper
    loss_func = @(theta) loss_wrapper(theta, X_lift, U, Y_lift, N_sim, N_traj, Tp, N_lift, m);

    % Optimization options
     % Optimization options
    options = optimoptions('fmincon', ...
        'Algorithm', 'sqp', ...  % Faster than 'interior-point'
        'Display', 'iter', ...  % Show iterations
        'MaxIterations', 500, ...
        'OptimalityTolerance', 1e-6, ...
        'StepTolerance', 1e-6);
    % Optimize using fminunc
    1
    theta_opt = fmincon(loss_func, theta_init, options);

    % Reshape back to matrices
    A_opt = reshape(theta_opt(1:N_lift*N_lift), N_lift, N_lift);
    B_opt = reshape(theta_opt(N_lift*N_lift+1:end), N_lift, m);
end

function loss = loss_wrapper(theta, X_lift, U, Y_lift, N_sim, N_traj, Tp, n_lift, n_u)
    % Extract A_lift and B_lift from theta
    A_lift = reshape(theta(1:n_lift*n_lift), n_lift, n_lift);
    B_lift = reshape(theta(n_lift*n_lift+1:end), n_lift, n_u);

    % Compute the multi-step loss
    loss = multi_step_loss(A_lift, B_lift, X_lift, U, Y_lift, N_sim, N_traj, Tp);
end

