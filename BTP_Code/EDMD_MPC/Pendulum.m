clc; clear;
close all
addpath('./Resources')
addpath('./Resources/qpOASES-3.1.0/interfaces/matlab') 
%% Dynamics
% Pendulum dynamics: θ'' = -(g/l)*sin(θ) - (b/ml^2)*θ' + (1/ml^2)*u
% State x = [θ; θ']
% Parameters
g = 9.81;  % gravity
l = 1.0;   % length
m = 0.1;   % mass
b = 0;   % damping coefficient

f_u = @(t,x,u)([x(2,:); 
                -(g/l)*sin(x(1,:)) - (b/(m*l^2))*x(2,:) + (1/(m*l^2))*u]);
n = 2;  % number of states
m = 1;  % number of control inputs

%% Discretization
tic
deltaT = 0.01;
% Runge-Kutta 4
k1 = @(t,x,u) (f_u(t,x,u));
k2 = @(t,x,u) (f_u(t,x + k1(t,x,u)*deltaT/2,u));
k3 = @(t,x,u) (f_u(t,x + k2(t,x,u)*deltaT/2,u));
k4 = @(t,x,u) (f_u(t,x + k1(t,x,u)*deltaT,u));
f_ud = @(t,x,u) (x + (deltaT/6) * (k1(t,x,u) + 2*k2(t,x,u) + 2*k3(t,x,u) + k4(t,x,u)));

rng(115123)
disp('Starting data collection')

Ntraj = 200; Nsim = 1000;

Cy = eye(n); % Output matrix: y = Cy*x
nD = 0; % Number of delays
ny = size(Cy,1); % Number of outputs

% Random control input forcing
Ubig = 2 * rand(Nsim, Ntraj) - 1;

% Random initial condition
Xcurrent = [2*rand(1, Ntraj) - 1;8*rand(1,Ntraj)-4];%(rand(n, Ntraj) * 2 - 1);
X = []; Y = []; U = [];

n_zeta = (nD + 1) * ny + nD * m; % Dimension of delay-embedded "state"
if nD == 0
    % Simplified "state" without delay embedding
    zeta_current = Cy * Xcurrent;
    for i = 1:Nsim
        Xnext = f_ud(0, Xcurrent, Ubig(i,:));
        zeta_prev = zeta_current;
        zeta_current = Cy * Xnext;
        X = [X zeta_prev];
        Y = [Y zeta_current];
        U = [U Ubig(i,:)];
        Xcurrent = Xnext;
    end
else
    % Delay-embedded "state" zeta_k = [y_{k}; u_{k-1}; y_{k-1}; ... ; u_{k-nd}; y_{k-nd}]
    zeta_current = [Cy * Xcurrent; NaN(nD * (ny + m), Ntraj)];

    for i = 1:Nsim
        Xnext = f_ud(0, Xcurrent, Ubig(i,:));
        zeta_prev = zeta_current;
        zeta_current = [[Cy * Xnext; Ubig(i,:)]; zeta_current(1:end-ny-m, :)];

        % Store data after delay embedding is complete
        if i > nD
            X = [X zeta_prev];
            Y = [Y zeta_current];
            U = [U Ubig(i,:)];
        end
        Xcurrent = Xnext;
    end
end

fprintf('Data collection DONE. Time taken : %f s \n', toc);

%% Basis functions
% basisFunction = 'rbf';
% Nrbf = 50;
% cent = rand(n_zeta,Nrbf)*2 - 1; % RBF centers
% rbf_type = 'thinplate';
% liftFun = @(xx)( [xx;rbf(xx,cent,rbf_type)] );
% Nlift = Nrbf + n_zeta;

%% Modified Basis functions 
% First, create a function that computes the custom observables for a single state
custom_observables = @(theta, theta_dot)([sin(theta); 
                                        cos(theta);
                                        theta_dot.*sin(theta);
                                        theta_dot.*cos(theta)]);

% Create the lifting function for the delay-embedded state
liftFun = @(xx)( lift_with_delays_and_observables(xx, nD, ny, m, custom_observables) );
Nlift = n_zeta + 4*(nD + 1); % Original dimension + 4 observables for each delay

% Function to compute lifting with delays and custom observables
function phi = lift_with_delays_and_observables(zeta, nD, ny, m, obs_fun)
    % zeta structure: [y_k; u_{k-1}; y_{k-1}; ...; u_{k-nD}; y_{k-nD}]
    % where y_k = [theta; theta_dot]
    
    % Initialize the lifted state vector with the original delay-embedded data
    phi = zeta;
    
    % Get the number of samples
    n_samples = size(zeta, 2);
    
    % Extract and process each delayed state component
    for i = 0:nD
        % Calculate position in zeta for current delay
        if i == 0
            idx = 1:ny;  % Current state (no delay)
        else
            idx = (ny + m)*i + 1:(ny + m)*i + ny;  % Delayed states
        end
        
        % Extract theta and theta_dot for current delay
        theta = zeta(idx(1), :);  % First row of current state component
        theta_dot = zeta(idx(2), :);  % Second row of current state component
        
        % Compute observables for current delay
        obs = obs_fun(theta, theta_dot);
        
        % Append observables to phi
        phi = [phi; obs];
    end
end

%% Lift
disp('Starting LIFTING')
tic
Xlift = liftFun(X);
Ylift = liftFun(Y);

fprintf('Lifting DONE. Time taken : %f s \n', toc);
%% Regression

disp('Starting REGRESSION for A,B,C')
tic
% W = [Ylift ; X];
% V = [Xlift ; U];
% VVt = V*V';
% WVt = W*V';
% ABC = WVt * pinv(VVt);
% Alift = ABC(1:Nlift,1:Nlift);
% Blift = ABC(1:Nlift,Nlift+1:end);
% Clift = ABC(Nlift+1:end,1:Nlift);
K = [Ylift;U]*pinv([Xlift;U]);
Alift = K(1:Nlift,1:Nlift); Blift = K(1:Nlift,Nlift+1:end); Clift = [eye(n),zeros(n,Nlift-n)];
fprintf('Regression for A, B, C DONE. Time taken : %f s \n',toc);

% Residual
fprintf( 'Regression residual : %f \n', norm(Ylift - Alift*Xlift - Blift*U,'fro') / norm(Ylift,'fro') );

%% Rank 1 Update
XXT_inv1 = pinv([Xlift;U]*[Xlift;U]'); 
XXT_inv = XXT_inv1;
K_bar = K; n_iter = 10000;
Norms = zeros(n_iter,1);
for iter = 1:n_iter
    XXT_inv_prev = XXT_inv;
    K_prev = K_bar;
    x = [2*pi*rand-pi;8*rand-4]; u = 2*rand-1; y=f_ud(0,x,u);
    cx = [liftFun(x);u] ; cy = [liftFun(y);u];
    scalar = 1 + cx'*XXT_inv_prev*cx;
    rank1_upd = eye(Nlift+m) - (cx*cx')*XXT_inv_prev/scalar; 
    XXT_inv = XXT_inv_prev*rank1_upd;
    K_bar = (K_prev + cy*cx'*XXT_inv_prev)*rank1_upd;
    Norms(iter) = rms(K_bar(:)-K_prev(:));
end
figure;
plot(1:length(Norms),(Norms),"LineWidth",2)
set(gca,"FontSize",20)
xlabel('No. of iterations',"FontWeight","bold")
ylabel('$||K_{i+1} - K_i ||$',"FontWeight","bold","Interpreter","latex",'Rotation', 0)
Alift = K(1:Nlift,1:Nlift); Blift = K(1:Nlift,Nlift+1:end);
%% Predictor comparison
close all;
Tmax = 5;
Nsim = Tmax / deltaT;

uprbs = 2 * rand(Nsim, 1) - 1; % (2*myprbs(Nsim,0.5) - 1);
u_dt = @(i)(uprbs(i + 1));
f_cont_d = @(t, xx)(f_ud(t, xx, u_dt(t)));

x0 = [2*rand - 1;8*rand-4];
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

% Local linearization
xloc = xp; % Use the last state (xp) as the linearization point
x = sym('x', [2 1]); syms u;
Aloc = double(subs(jacobian(f_ud(0, x, u), x), [x; u], [xloc; urand]));
Bloc = double(subs(jacobian(f_ud(0, x, u), u), [x; u], [xloc; urand]));
cloc = double(subs(f_ud(0, x, u), [x; u], [xloc; urand])) - Aloc * xloc - Bloc * urand;

% Initial conditions
x_true = xp; % Start true dynamics simulation from xp
xlift = liftFun(xstart);

% Simulation
for i = 0:Nsim-1
    % True dynamics
    x_true = [x_true, f_ud(0, x_true(:, end), u_dt(i))];
    
    % Koopman predictor
    xlift = [xlift, Alift * xlift(:, end) + Blift * u_dt(i)];
    
    % Local linearization predictor
    xloc = [xloc, Aloc * xloc(:, end) + Bloc * u_dt(i) + cloc];
end

theta_true = x_true(1, :);
theta_dot_true = x_true(2, :);
theta_pred = Clift(1, :) * xlift;
theta_dot_pred = Clift(2, :) * xlift;

% figure
% stairs((0:Nsim-1)*deltaT,u_dt(0:Nsim-1),'linewidth',2); hold on
% title('Control input'); xlabel('time [s]')

figure
lw_koop = 2;
subplot(2,1,1);
plot((0:Nsim)*deltaT,theta_true,'-b','linewidth', lw_koop); hold on
plot((0:Nsim)*deltaT,theta_pred, '--r','linewidth',lw_koop)
title('$\theta$ (Angle)','interpreter','latex','fontsize',20);
xlabel('time (s)','interpreter','latex','fontsize',14)
LEG = legend('True','Koopman');
set(LEG,'Interpreter','latex','location','northeast','fontsize',14)
set(gca,'FontSize',14);
subplot(2,1,2);
plot((0:Nsim)*deltaT,theta_dot_true,'-b','linewidth', lw_koop); hold on
plot((0:Nsim)*deltaT,theta_dot_pred, '--r','linewidth',lw_koop)
title('$\dot{\theta}$ (Angular Velocity)','interpreter','latex','fontsize',20);
xlabel('time (s)','interpreter','latex','fontsize',14)
LEG = legend('True','Koopman');
set(LEG,'Interpreter','latex','location','northeast','fontsize',14)
set(gca,'FontSize',16);

%% Feedback control 
close all;
n_setpoints = 4;       % Number of setpoints
Tmax = 2 * n_setpoints; % Total simulation time
% Initial condition (start near bottom position)
x0 = [0; 0]; % Slightly perturbed from bottom position
Nsim = Tmax / deltaT;   % Total number of simulation steps

% Generate random setpoints for theta_ref
theta_refs = 2*pi * (rand(1, n_setpoints) - 0.5); % Random values in range [-pi, pi]

% Calculate duration of each setpoint
duration = Nsim / n_setpoints;

% Create yrr
tspan = linspace(0,Tmax,Nsim);
% theta_ref_trajectory = pi*sin(tspan);
% thetadot_ref_trajectory = pi*cos(tspan);
%x0 = [theta_ref_trajectory(1);0];%thetadot_ref_trajectory(1)];
theta_ref_trajectory = repelem(theta_refs, duration); % Repeat each setpoint for 'duration' steps
thetadot_ref_trajectory = zeros(1, Nsim);            % Zero velocity for thetadot

% Combine into yrr
yrr = [theta_ref_trajectory; thetadot_ref_trajectory];

% % Define step reference for both states
% theta_ref = pi*0.75 ; % Target angle
% thetadot_ref = 0; % Target angular velocity (zero for setpoint)
% 
% % Create reference trajectories
% yrr = [theta_ref*ones(1,Nsim);  % Reference for theta
%        thetadot_ref*ones(1,Nsim)]; % Reference for thetadot

% Define Koopman controller
C = zeros(n,Nlift); 
C(1:n,1:n) = eye(n); % Extract both states from lifted state

% Weight matrices (adjusted for two-state tracking)
Q = diag([1000, 1]); % Higher weight on angle tracking
QN = diag([10000, 1]); % Higher terminal weight
R = 5;

% Prediction horizon
Tpred = 0.2;
Np = round(Tpred / deltaT);

% Constraints
u_min = -2; u_max = 2;
theta_min = -pi; theta_max = pi;
thetadot_min = -4; thetadot_max = 4;

xlift_min = [theta_min; thetadot_min; -inf*ones(Nlift-2,1)];
xlift_max = [theta_max; thetadot_max; inf*ones(Nlift-2,1)];

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
theta = XX_koop(1,:);
% Plot theta
subplot(3,1,1)
plot(time, theta, 'LineWidth', 3)
hold on
stairs(time(1:end-1), yrr(1,:), 'LineWidth', 2)
plot(time, theta_max*ones(1,Nsim+1), '--k', 'LineWidth', 2)
plot(time, theta_min*ones(1,Nsim+1), '--k', 'LineWidth', 2)
ylabel('$\theta$ (Angle)','Interpreter','latex','FontSize',17,'FontWeight','bold')
legend('Actual', 'Reference', 'Constraints','FontSize',12)
grid on
set(gca,'fontsize',14)
title('Pendulum KMPC','FontSize',17)

theta_dot = XX_koop(2,:);
% Plot thetadot
subplot(3,1,2)
plot(time, theta_dot, 'LineWidth', 3)
hold on
stairs(time(1:end-1), yrr(2,:), 'LineWidth', 2)
plot(time, thetadot_max*ones(1,Nsim+1), '--k', 'LineWidth', 2)
plot(time, thetadot_min*ones(1,Nsim+1), '--k', 'LineWidth', 2)
ylabel('$\dot{\theta}$ (Angular velocity)','Interpreter','latex','FontSize',17,'FontWeight','bold')
legend('Actual', 'Reference', 'Constraints','FontSize',12)
set(gca,'fontsize',14)
grid on

% Plot control input
subplot(3,1,3)
plot(time(1:end-1), u_max*ones(1,Nsim), '--k', 'LineWidth', 2)
hold on
plot(time(1:end-1), u_min*ones(1,Nsim), '--k', 'LineWidth', 2)
plot(time(1:end-1), UU_koop, 'LineWidth', 3,'Color','r')
ylabel('Control input u','Interpreter','latex','FontSize',17,'FontWeight','bold')
xlabel('Time (s)','FontSize',17,'FontWeight','bold')
legend('Constraints','FontSize',12)
set(gca,'fontsize',14)
grid on
%% Animate the pendulum

% Visualize and animate pendulum motion
figure('Position', [100, 100, 700, 525]);
hold on;

% Initialize pendulum graphics
bob = plot(0, 0, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b'); % Pendulum bob
rod = line([0, 0], [0, 0], 'Color', 'k', 'LineWidth', 2); % Pendulum rod
time_text = text(-l * 1.1, l * 1.1, 'Time: 0.00 s', 'FontSize', 18, 'Color', 'r'); % Time text

% Set axis limits and aspect ratio
axis equal;
set(gca,'FontSize',20)
grid on;
axis([-l, l, -l, l] * 1.2);
xlabel('X axis (m)');
ylabel('Y axis (m)');
title('Koopman Control for Pendulum');

% % Create a video writer object
% video = VideoWriter('pendulum_animation.avi'); % Specify the file name
% video.FrameRate = 1 / deltaT; % Set the frame rate
% open(video); % Open the video writer

for k = 1:length(theta)
    % Calculate pendulum position
    xpos = l * sin(theta(k));
    ypos = -l * cos(theta(k));

    % Update graphics
    set(bob, 'XData', xpos, 'YData', ypos); % Update bob position
    set(rod, 'XData', [0, xpos], 'YData', [0, ypos]); % Update rod position
    set(time_text, 'String', sprintf('Time : %.2f s', (k-1)*deltaT)); % Update time text

    % % Capture the current frame
    % frame = getframe(gcf); % Capture the figure as a frame
    % writeVideo(video, frame); % Write the frame to the video

    % Pause to control animation speed
    pause(deltaT);
end

hold off;
% Close the video writer
%close(video);

%disp('Animation saved as pendulum_animation.avi');