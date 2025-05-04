clc; clear;
close all
addpath('./Resources')
addpath('./Resources/qpOASES-3.1.0/interfaces/matlab') 
%% Dynamics
% CartPole dynamics 
% State x = [x; theta; xdot; theta_dot]
% Parameters
m_p = 1;   % Mass of the pole (kg)
m_c = 4.0;   % Mass of the cart (kg)
l = 1;     % Length of the pole (m)
g = -9.81;    % Acceleration due to gravity (m/s^2)

f_u = @(t, x, u) ([...
    x(3,:); ... % x_dot
    x(4,:); ... % theta_dot
    (u + m_p * l * (x(4,:).^2 .* sin(x(2,:)) - ... % x_ddot
        ((g * sin(x(2,:)) + cos(x(2,:)) .* ...
        ((-u - m_p * l * x(4,:).^2 .* sin(x(2,:))) / (m_c + m_p))) ...
        ./ (l * (4/3 - (m_p * cos(x(2,:)).^2) / (m_c + m_p)))) .* cos(x(2,:)))) ...
    ./ (m_c + m_p); ...
    (g * sin(x(2,:)) + cos(x(2,:)) .* ... % theta_ddot
        ((-u - m_p * l * x(4,:).^2 .* sin(x(2,:))) / (m_c + m_p))) ...
    ./ (l * (4/3 - (m_p * cos(x(2,:)).^2) / (m_c + m_p)))]);

n = 4;  % number of states
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

Ntraj = 2000; Nsim = 2000;

Cy = eye(n); % Output matrix: y = Cy*x
nD = 0; % Number of delays
ny = size(Cy,1); % Number of outputs

% Random control input forcing
Ubig = 20 * rand(Nsim, Ntraj) - 10;

% Random initial condition
Xcurrent = [2*rand(1, Ntraj) - 1;2*pi*rand(1,Ntraj)-pi;4*rand(1, Ntraj) - 2;4*rand(1, Ntraj) - 2];%(rand(n, Ntraj) * 2 - 1);
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
% First, create a function that computes the custom observables for a single state
custom_observables = @(x,theta,xdot, theta_dot)([sin(theta); 
                                        cos(theta);
                                        theta_dot.*sin(theta);
                                        theta_dot.*cos(theta); ...
                                        x.^2;xdot.^2;x.*xdot;x.^2.*xdot;x.*xdot.^2]);

% Create the lifting function for the delay-embedded state
liftFun = @(xx)( lift_with_delays_and_observables(xx, nD, ny, m, custom_observables) );
%Nlift = n_zeta + 7*(nD + 1); % Original dimension + 4 observables for each delay

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
        x = zeta(idx(1),:); theta = zeta(idx(2), :); 
        xdot = zeta(idx(3),:);theta_dot = zeta(idx(4), :); 

        % Compute observables for current delay
        obs = obs_fun(x, theta, xdot, theta_dot);

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
fprintf('Lifting DONE. Time taken : %f s \n', toc);
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
fprintf('Regression for A, B, C DONE. Time taken : %f s \n',toc);

% Residual
fprintf( 'Regression residual : %f \n', norm(Ylift - Alift*Xlift - Blift*U,'fro') / norm(Ylift,'fro') );

%% Predictor comparison
close all;
Tmax = 5;
Nsim = Tmax / deltaT;

uprbs = 20 * rand(Nsim, 1) - 10; % (2*myprbs(Nsim,0.5) - 1);
u_dt = @(i)(uprbs(i + 1));
f_cont_d = @(t, xx)(f_ud(t, xx, u_dt(t)));

x0 = 2*rand(n,1)-1;
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
    xlift = [xlift, Alift * xlift(:, end) + Blift * u_dt(i)];
    
end

pos_true = x_true(1, :);theta_true = x_true(2, :);
vel_true = x_true(3, :);omega_true = x_true(4, :);
pos_pred = Clift(1, :) * xlift;theta_pred = Clift(2, :) * xlift;
vel_pred = Clift(3, :) * xlift;omega_pred = Clift(4, :) * xlift;

% figure
% stairs((0:Nsim-1)*deltaT,u_dt(0:Nsim-1),'linewidth',2); hold on
% title('Control input'); xlabel('time [s]')

figure
lw_koop = 2;
subplot(2,2,1);
plot((0:Nsim)*deltaT,pos_true,'-b','linewidth', lw_koop); hold on
plot((0:Nsim)*deltaT,pos_pred, '--r','linewidth',lw_koop)
title('$x$ (Position)','interpreter','latex','fontsize',20);
xlabel('time (s)','interpreter','latex','fontsize',14)
LEG = legend('True','Koopman');
set(LEG,'Interpreter','latex','location','northeast','fontsize',14)
set(gca,'FontSize',14);
subplot(2,2,2);
plot((0:Nsim)*deltaT,theta_true,'-b','linewidth', lw_koop); hold on
plot((0:Nsim)*deltaT,theta_pred, '--r','linewidth',lw_koop)
title('$\theta$ (Angle)','interpreter','latex','fontsize',20);
xlabel('time (s)','interpreter','latex','fontsize',14)
LEG = legend('True','Koopman');
set(LEG,'Interpreter','latex','location','northeast','fontsize',14)
set(gca,'FontSize',16);
subplot(2,2,3);
plot((0:Nsim)*deltaT,vel_true,'-b','linewidth', lw_koop); hold on
plot((0:Nsim)*deltaT,vel_pred, '--r','linewidth',lw_koop)
title('$\dot{x}$ (Velocity)','interpreter','latex','fontsize',20);
xlabel('time (s)','interpreter','latex','fontsize',14)
LEG = legend('True','Koopman');
set(LEG,'Interpreter','latex','location','northeast','fontsize',14)
set(gca,'FontSize',14);
subplot(2,2,4);
plot((0:Nsim)*deltaT,omega_true,'-b','linewidth', lw_koop); hold on
plot((0:Nsim)*deltaT,omega_pred, '--r','linewidth',lw_koop)
title('$\dot{\theta}$ (Angular Velocity)','interpreter','latex','fontsize',20);
xlabel('time (s)','interpreter','latex','fontsize',14)
LEG = legend('True','Koopman');
set(LEG,'Interpreter','latex','location','northeast','fontsize',14)
set(gca,'FontSize',16);

%% Feedback control 
close all;
n_setpoints = 4;       % Number of setpoints
Tmax = 7 * n_setpoints; % Total simulation time
% Initial
x0 = [0;0;0;0];
Nsim = Tmax / deltaT;   % Total number of simulation steps

% Generate random setpoints for x_ref
x_refs = [-5,5,-3,2];%9*rand(1,n_setpoints)-4.5;
theta_refs = zeros(1, n_setpoints); % Pole up control

% Calculate duration of each setpoint
duration = Nsim / n_setpoints;

% Create yrr
tspan = linspace(0,Tmax,Nsim);
x_ref_trajectory = repelem(x_refs,duration); 
theta_ref_trajectory = repelem(theta_refs, duration); % Repeat each setpoint for 'duration' steps
xdot_ref_trajectory = zeros(1,Nsim);
thetadot_ref_trajectory = zeros(1, Nsim);  % Zero velocity for thetadot

% Combine into yrr
yrr = [x_ref_trajectory; theta_ref_trajectory; xdot_ref_trajectory; thetadot_ref_trajectory];

% Define Koopman controller
C = zeros(n,Nlift); 
C(1:n,1:n) = eye(n); % Extract states from lifted state

% Weight matrices
Q = diag([10, 20, 1, 1]); 
QN = diag([1000, 10000, 1, 1]); % Higher terminal weight
R = 0.01;

% Prediction horizon
Tpred = 3;
Np = round(Tpred / deltaT);

% Constraints
u_min = -20; u_max = 20;
x_min = -6; x_max = 6;
theta_min = -pi; theta_max = pi;
xdot_min = -5; xdot_max = 5;
thetadot_min = -8; thetadot_max = 8;

xlift_min = [x_min; theta_min; xdot_min; thetadot_min; -inf*ones(Nlift-n,1)];
xlift_max = [x_max; theta_max; xdot_max; thetadot_max; inf*ones(Nlift-n,1)];

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
x = XX_koop(1,:);theta = XX_koop(2,:);xdot = XX_koop(3,:);theta_dot = XX_koop(4,:);
%% Plot results
figure('Position', [100, 100, 800, 600])
time = (0:Nsim)*deltaT;
subplot(2,2,1)
plot(time, x, 'LineWidth', 3)
hold on
stairs(time(1:end-1), yrr(1,:), 'LineWidth', 2)
plot(time, x_max*ones(1,Nsim+1), '--k', 'LineWidth', 2)
plot(time, x_min*ones(1,Nsim+1), '--k', 'LineWidth', 2)
ylabel('$x$ (Position)','Interpreter','latex','FontSize',17,'FontWeight','bold')
legend('Actual', 'Reference', 'Constraints','FontSize',12)
grid on
set(gca,'fontsize',14)
title('Pendulum KMPC','FontSize',17)
subplot(2,2,2)
plot(time, theta, 'LineWidth', 3)
hold on
stairs(time(1:end-1), yrr(2,:), 'LineWidth', 2)
plot(time, theta_max*ones(1,Nsim+1), '--k', 'LineWidth', 2)
plot(time, theta_min*ones(1,Nsim+1), '--k', 'LineWidth', 2)
ylabel('$\theta$ (Angle)','Interpreter','latex','FontSize',17,'FontWeight','bold')
legend('Actual', 'Reference', 'Constraints','FontSize',12)
grid on
set(gca,'fontsize',14)
title('Pendulum KMPC','FontSize',17)
subplot(2,2,3)
plot(time, xdot, 'LineWidth', 3)
hold on
stairs(time(1:end-1), yrr(3,:), 'LineWidth', 2)
plot(time, xdot_max*ones(1,Nsim+1), '--k', 'LineWidth', 2)
plot(time, xdot_min*ones(1,Nsim+1), '--k', 'LineWidth', 2)
ylabel('$\dot{x}$ (Velocity)','Interpreter','latex','FontSize',17,'FontWeight','bold')
legend('Actual', 'Reference', 'Constraints','FontSize',12)
set(gca,'fontsize',14)
grid on
subplot(2,2,4)
plot(time, theta_dot, 'LineWidth', 3)
hold on
stairs(time(1:end-1), yrr(4,:), 'LineWidth', 2)
plot(time, thetadot_max*ones(1,Nsim+1), '--k', 'LineWidth', 2)
plot(time, thetadot_min*ones(1,Nsim+1), '--k', 'LineWidth', 2)
ylabel('$\dot{\theta}$ (Angular velocity)','Interpreter','latex','FontSize',17,'FontWeight','bold')
legend('Actual', 'Reference', 'Constraints','FontSize',12)
set(gca,'fontsize',14)
grid on

% Plot control input
%subplot(3,1,3)
figure;
plot(time(1:end-1), u_max*ones(1,Nsim), '--k', 'LineWidth', 2)
hold on
plot(time(1:end-1), u_min*ones(1,Nsim), '--k', 'LineWidth', 2)
plot(time(1:end-1), UU_koop, 'LineWidth', 3,'Color','r')
ylabel('Control input u','Interpreter','latex','FontSize',17,'FontWeight','bold')
xlabel('Time (s)','FontSize',17,'FontWeight','bold')
legend('Constraints','FontSize',12)
set(gca,'fontsize',14)
grid on

%% Animation
figure('Position', [100, 100, 700, 525]);
hold on;
axis equal;
title('Koopman Control for Cart-Pole')
xlabel("X-axis (m)")
ylabel("Y-axis (m)")
set(gca,'fontsize',20);
%xlim([-2, 2]);
ylim([-1.5, 1.5]);
grid on;

cartWidth = 0.3; % Cart width
cartHeight = 0.2; % Cart height
poleWidth = 0.02; % Pole width
wheelRadius = 0.05; % Radius of the wheels

cart = rectangle('Position', [0, 0, cartWidth, cartHeight], 'FaceColor', 'blue');
pole = line([0, 0], [0, 0], 'LineWidth', 2, 'Color', 'red');
time_text = text(-l * 1.5, l * 1.3, 'Time: 0.00 s', 'FontSize', 20, 'Color', 'b'); % Time text

% Wheels (represented as circles)
theta_circle = linspace(0, 2*pi, 50);
wheel1 = fill(wheelRadius * cos(theta_circle), wheelRadius * sin(theta_circle), 'black');
wheel2 = fill(wheelRadius * cos(theta_circle), wheelRadius * sin(theta_circle), 'black');

% % Create a video writer object
% video = VideoWriter('cart_pole_animation.avi'); % Specify the file name
% video.FrameRate = 1 / deltaT; % Set the frame rate
% open(video); % Open the video writer

for i = 1:length(time)
    % Update cart position
    cart.Position = [x(i) - cartWidth/2, -cartHeight/2, cartWidth, cartHeight];

    % Update pole position
    poleX = x(i) + l * sin(theta(i));
    poleY = l * cos(theta(i));
    pole.XData = [x(i), poleX];
    pole.YData = [0, poleY];

    % Update wheel positions
    set(wheel1, 'XData', x(i) - cartWidth/4 + wheelRadius * cos(theta_circle), ...
                'YData', -cartHeight/2 - wheelRadius + wheelRadius * sin(theta_circle));
    set(wheel2, 'XData', x(i) + cartWidth/4 + wheelRadius * cos(theta_circle), ...
                'YData', -cartHeight/2 - wheelRadius + wheelRadius * sin(theta_circle));
    set(time_text, 'Position', [x(i)-l * 1.5, l * 1.3]); % Dynamically move with the cart
    set(time_text, 'String', sprintf('Time : %.2f s', (i-1)*deltaT)); % Update time text

    % % Capture the current frame
    % frame = getframe(gcf); % Capture the figure as a frame
    % writeVideo(video, frame); % Write the frame to the video

    % Pause to control animation speed
    pause(deltaT);
end

% % Close the video writer
% close(video);
% 
% disp('Animation saved as cart_pole_animation.avi');
