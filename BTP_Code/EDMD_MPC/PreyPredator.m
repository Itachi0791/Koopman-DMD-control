clc; clear;
close all
addpath('./Resources')
addpath('./Resources/qpOASES-3.1.0/interfaces/matlab') 

%% Dynamics
% Lotka-Volterra model dynamics: xdot =  alpha*x - beta*xy  + ux ,
%                                ydot = -gamma*y + delta*xy + uy
% State x = [x; y]
% Parameters for the Lotka-Volterra model
alpha = 2/3;    % Prey birth rate
beta = 4/3;     % Predator-prey interaction rate
gamma = 1;      % Predator death rate
delta = 1;      % Predator reproduction rate

f_u = @(t,x,u)([alpha * x(1) - beta * x(1) * x(2) + u(1);   % Prey equation 
               -gamma * x(2) + delta * x(1) * x(2) + u(2);  % Predator equation 
                ]);
n = 2;  % number of states
m = 2;  % number of control inputs

%% Discretization
tic
deltaT = 0.02;
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
nD = 2; % Number of delays
ny = size(Cy,1); % Number of outputs

% Random control input forcing between 0 and 1
Ubig = rand(Nsim, m, Ntraj);

% Random initial condition
x_max = 3.2; x_min = 0.1; y_max = 2; y_min = 0.1;
Xcurrent = [x_min + (x_max-x_min)*rand(1, Ntraj);y_min + (y_max-y_min)*rand(1, Ntraj)];
X = []; Y = []; U = [];

n_zeta = (nD + 1) * ny + nD * m; % Dimension of delay-embedded "state"
if nD == 0
    % Simplified "state" without delay embedding
    zeta_current = Cy * Xcurrent;
    for i = 1:Nsim
        Xnext = f_ud(0, Xcurrent, squeeze(Ubig(i,:,:)));
        zeta_prev = zeta_current;
        zeta_current = Cy * Xnext;
        X = [X zeta_prev];
        Y = [Y zeta_current];
        U = [U squeeze(Ubig(i,:,:))];
        Xcurrent = Xnext;
    end
else
    % Delay-embedded "state" zeta_k = [y_{k}; u_{k-1}; y_{k-1}; ... ; u_{k-nd}; y_{k-nd}]
    zeta_current = [Cy * Xcurrent; NaN(nD * (ny + m), Ntraj)];

    for i = 1:Nsim
        Xnext = f_ud(0, Xcurrent, squeeze(Ubig(i,:,:)));
        zeta_prev = zeta_current;
        zeta_current = [[Cy * Xnext; squeeze(Ubig(i,:,:))]; zeta_current(1:end-ny-m, :)];

        % Store data after delay embedding is complete
        if i > nD
            X = [X zeta_prev];
            Y = [Y zeta_current];
            U = [U squeeze(Ubig(i,:,:))];
        end
        Xcurrent = Xnext;
    end
end

fprintf('Data collection DONE. Time taken - %f s\n', toc);

%% Basis functions
basisFunction = 'rbf';
Nrbf = 100;
cent = rand(n_zeta,Nrbf)*2 - 1; % RBF centers
rbf_type = 'invmultquad';
liftFun = @(xx)( [xx;rbf(xx,cent,rbf_type)] );
Nlift = Nrbf + n_zeta;

%% Lift
disp('Starting LIFTING')
tic
Xlift = liftFun(X);
Ylift = liftFun(Y);

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
Tmax = 20;
Nsim = Tmax / deltaT;

uprbs = zeros(Nsim,2);%0.2*rand(Nsim, 2); 
u_dt = @(i)(uprbs(i + 1,:));
f_cont_d = @(t, xx)(f_ud(t, xx, u_dt(t)));

x0 = [1;0.5];%[x_min + (x_max-x_min)*rand;y_min + (y_max-y_min)*rand];
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

% % Local linearization
% xloc = xp; % Use the last state (xp) as the linearization point
% x = sym('x', [2 1]); syms u;
% Aloc = double(subs(jacobian(f_ud(0, x, u), x), [x; u], [xloc; urand]));
% Bloc = double(subs(jacobian(f_ud(0, x, u), u), [x; u], [xloc; urand]));
% cloc = double(subs(f_ud(0, x, u), [x; u], [xloc; urand])) - Aloc * xloc - Bloc * urand;

% Initial conditions
x_true = xp; % Start true dynamics simulation from xp
xlift = liftFun(xstart);

% Simulation
for i = 0:Nsim-1
    % True dynamics
    x_true = [x_true, f_ud(0, x_true(:, end), u_dt(i))];
    
    % Koopman predictor
    xlift = [xlift, Alift * xlift(:, end) + Blift * u_dt(i)'];
    
    % % Local linearization predictor
    % xloc = [xloc, Aloc * xloc(:, end) + Bloc * u_dt(i) + cloc];
end

p_true = x_true(1, :);
q_true = x_true(2, :);
p_pred = Clift(1, :) * xlift;
q_pred = Clift(2, :) * xlift;

% figure
% stairs((0:Nsim-1)*deltaT,u_dt(0:Nsim-1),'linewidth',2); hold on
% title('Control input'); xlabel('time [s]')

figure
lw_koop = 2;
subplot(2,1,1);
plot((0:Nsim)*deltaT, p_true,'-b','linewidth', lw_koop); hold on
plot((0:Nsim)*deltaT, p_pred, '--r','linewidth',lw_koop)
ylabel('x (Prey)','interpreter','latex','fontsize',20);
title('Koopman Prediction','interpreter','latex','fontsize',20)
LEG = legend('True','Koopman');
set(LEG,'Interpreter','latex','location','northeast','fontsize',14)
set(gca,'FontSize',16);
subplot(2,1,2);
plot((0:Nsim)*deltaT,q_true,'-b','linewidth', lw_koop); hold on
plot((0:Nsim)*deltaT,q_pred, '--r','linewidth',lw_koop)
ylabel('y (Predators)','interpreter','latex','fontsize',20);
xlabel('time (s)','interpreter','latex','fontsize',14)
LEG = legend('True','Koopman');
set(LEG,'Interpreter','latex','location','northeast','fontsize',14)
set(gca,'FontSize',16);

