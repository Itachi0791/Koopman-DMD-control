clc; clear;
close all
addpath('./Resources')
% addpath('./Resources/qpOASES-3.1.0/interfaces/matlab') 
g = 1;  % gravity
l = 1.0;   % length
m = 0.1;   % mass
b = 0;     % damping coefficient
dt = 0.1; % Time step
N_sim = 20; % No. of simulations
N_traj = 200; % No. of trajectory points in each sim

traj_data = generate_pendulum_trajectories(N_sim, N_traj, m, l, b, g, dt);
%save("DataFiles\pendulum_data2.mat","traj_data")
%load("DataFiles\pendulum_data2.mat")

%% Data conversion with Delay Embedding
tic
disp('Starting Data Extraction.')

n = 2;         % State dimension
m = 1;         % Input dimension
delt = 0.01;   % Time step size (if needed)
nD = 0;        % Delay depth
%y_bias = min(min(st))

% Compute delay-embedded state size
n_zeta = (nD + 1) * n + nD * m;

% Total samples after delay embedding
total_samples = N_sim * (N_traj - nD);  

% Preallocate matrices for efficiency
X = zeros(n_zeta, total_samples);
Y = zeros(n_zeta, total_samples);
U = zeros(m, total_samples);

sample_idx = 1;  % Index for storing data efficiently

for s = 1:N_sim
    % Extract trajectory data
    traj_x = traj_data{s}.x;
    traj_y = traj_data{s}.y;
    traj_u = traj_data{s}.u;

    if nD == 0
        % No delay embedding, just store raw data
        idx_range = (s-1)*N_traj+1 : s*N_traj;
        X(:, idx_range) = traj_x;
        Y(:, idx_range) = traj_y;
        U(:, idx_range) = traj_u;
    else
        traj_x = [zeros(n,nD) traj_x];
        traj_y = [zeros(n,nD) traj_y];
        traj_u = [zeros(m,nD) traj_u];
        % Initialize delay-embedded vector
        zeta_prev = zeros(n_zeta, 1);
    
        for i = nD+1:N_traj
            % Construct delay-embedded state
            zeta_prev(1:n) = traj_x(:, i);  % Current state (y_k)
            
            % Fill in past states and inputs
            for d = 1:nD
                zeta_prev(n + (d-1)*(n + m) + (1:m)) = traj_u(:, i-d);  % Past input (u_{k-d})
                zeta_prev(n + (d-1)*(n + m) + m + (1:n)) = traj_y(:, i-d);  % Past state (y_{k-d})
            end
            
            % Update current delay-embedded state
            zeta_current = [traj_y(:, i); traj_u(:, i); zeta_prev(1:end-n-m)];
            
            % Store in preallocated matrices
            X(:, sample_idx) = zeta_prev(:);
            Y(:, sample_idx) = zeta_current(:);
            U(:, sample_idx) = traj_u(:, i);
            
            sample_idx = sample_idx + 1;  % Increment index
        end
    end
end

fprintf('Data extraction DONE. Time taken: %f s \n', toc);

%% Lifting to Koopman Hilbert Space
degree=3; 
x_bar = 5*max(abs(X(1,:)));y_bar=5*max(abs(X(2,:)));
custom_observables = generate_observables(15);
mn_vals = generateMNArray(degree);

if nD == 0
   % liftFun = @(xx) [xx;custom_observables(xx(1,:),xx(2,:))];
    liftFun = @(xx) [xx;legendre_observables(xx(1,:),xx(2,:),mn_vals,x_bar,y_bar)];
else
    % Create the lifting function for the delay-embedded state
    liftFun = @(xx)( lift_with_delays_and_observables(xx, nD, n, m, custom_observables) );
end

% Function to compute lifting with delays and custom observables
function phi = lift_with_delays_and_observables(zeta, nD, n, m, custom_observables)
    % zeta structure: [y_k; u_{k-1}; y_{k-1}; ...; u_{k-nD}; y_{k-nD}]
    % where y_k = [x; y] (n-dimensional state)
    
    % Initialize the lifted state vector with the original delay-embedded data
    phi = zeta;
    
    % Get the number of samples
    n_samples = size(zeta, 2);
    
    % Extract current state (y_k)
    y_k = zeta(1:n, :);  % Current state (n x n_samples)
    
    % Compute observables for the current state
    obs = custom_observables(y_k(1, :), y_k(2, :));  % Assuming x and y are the first two states
    
    % Append observables to phi
    phi = [phi; obs];
    
    % Process delayed states and inputs
    for i = 1:nD
        % Extract delayed input (u_{k-i})
        u_idx = n + (i-1)*(n + m) + (1:m);
        u_k_i = zeta(u_idx, :);  % Past input (m x n_samples)
        
        % Extract delayed state (y_{k-i})
        y_idx = n + (i-1)*(n + m) + m + (1:n);
        y_k_i = zeta(y_idx, :);  % Past state (n x n_samples)
        
        % Compute observables for the delayed state
        obs_delayed = custom_observables(y_k_i(1, :), y_k_i(2, :));  % Assuming x and y are the first two states
        
        % Append observables to phi
        phi = [phi; obs_delayed];
    end
end

disp('Starting LIFTING')
tic
Xlift = liftFun(X);
Ylift = liftFun(Y);
Nlift = size(Xlift,1);
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
Clift = [eye(n),zeros(n,Nlift-n)];%Clift = ABC(Nlift+1:end,1:Nlift);
K = [Ylift;U]*pinv([Xlift;U]);
Alift = K(1:Nlift,1:Nlift); Blift = K(1:Nlift,Nlift+1:end); Clift = [eye(n),zeros(n,Nlift-n)];
fprintf('Regression for A, B, C DONE. Time taken : %f s \n',toc);
Blift(1) = 0;Blift(2)=1;
% Residual
fprintf( 'Regression residual : %f \n', norm(Ylift - Alift*Xlift - Blift*U,'fro')/ norm(Ylift,'fro') );

%% Predictor comparison
close all;
traj_id = randi(N_sim); % Trajectory to test
n_traj = N_traj-1; % Number of time steps to simulate
u_traj = traj_data{traj_id}.u(:, 1:n_traj); % Input trajectory
x_true = traj_data{traj_id}.x(:, 1:n_traj+1); % True state trajectory

% Initialize state and lifted state
x = zeros(n, n_traj+1); % State trajectory (n x n_traj+1)
x(:, 1) = x_true(:, 1); % Initial state
x_start = x_true(:, 1); % Initial state for lifting

% Initialize lifted state
if nD == 0
    % No delay embedding: lift the initial state directly
    x_lift = zeros(Nlift, n_traj+1); % Lifted state trajectory (Nlift x n_traj+1)
    x_lift(:, 1) = liftFun(x_start); % Lift the initial state
else
    % With delay embedding: construct the initial lifted state
    % The initial lifted state includes the current state and past states/inputs
    x_lift = zeros(Nlift, n_traj+1); % Lifted state trajectory (Nlift x n_traj+1)
    
    % Construct the initial zeta vector for delay embedding
    zeta_init = zeros(n_zeta, 1); % Initial zeta vector (n_zeta x 1)
    zeta_init(1:n) = x_start; % Current state (y_k)
    
    % Fill in past states and inputs (if available)
    for d = 1:nD
        if d <= size(traj_data{traj_id}.x, 2) - 1
            zeta_init(n + (d-1)*(n + m) + (1:m)) = traj_data{traj_id}.u(:, d); % Past input (u_{k-d})
            zeta_init(n + (d-1)*(n + m) + m + (1:n)) = traj_data{traj_id}.x(:, d); % Past state (y_{k-d})
        else
            % If past data is not available, assume zeros (or other initial conditions)
            zeta_init(n + (d-1)*(n + m) + (1:m)) = zeros(m, 1); % Zero past input
            zeta_init(n + (d-1)*(n + m) + m + (1:n)) = zeros(n, 1); % Zero past state
        end
    end
    
    % Lift the initial zeta vector
    x_lift(:, 1) = liftFun(zeta_init);
end

% Simulation loop
for i = 1:n_traj
    % Koopman prediction
    x_lift(:, i+1) = Alift * x_lift(:, i) + Blift * u_traj(:, i); % Predict next lifted state
    
    % Extract the predicted state from the lifted state
    x_pred = Clift(1:n, :) * x_lift(:, i+1); % Predicted state (n x 1)
    
    % If nD > 0, construct the new zeta vector for the next step
    if nD > 0
        % Construct the new zeta vector for delay embedding
        zeta_new = zeros(n_zeta, 1); % New zeta vector (n_zeta x 1)
        zeta_new(1:n) = x_pred; % Current state (y_{k+1})
        
        % Fill in past states and inputs
        for d = 1:nD
            if i - d + 1 >= 1
                zeta_new(n + (d-1)*(n + m) + (1:m)) = u_traj(:, i - d + 1); % Past input (u_{k-d+1})
                zeta_new(n + (d-1)*(n + m) + m + (1:n)) = x(:, i - d + 1); % Past state (y_{k-d+1})
            % else
            %     % If past data is not available, assume zeros (or other initial conditions)
            %     zeta_new(n + (d-1)*(n + m) + (1:m)) = zeros(m, 1); % Zero past input
            %     zeta_new(n + (d-1)*(n + m) + m + (1:n)) = zeros(n, 1); % Zero past state
            end
        end
        
        % Lift the new zeta vector
        x_lift(:, i+1) = liftFun(zeta_new);
    else
        % If nD == 0, lift the predicted state directly
        x_lift(:, i+1) = liftFun(x_pred);
    end
end

% Extract true and predicted states for plotting
x_val_true = x_true(1, :); % True x trajectory
y_val_true = x_true(2, :); % True y trajectory
x_val_pred = Clift(1, :) * x_lift; % Predicted x trajectory
y_val_pred = Clift(2, :) * x_lift; % Predicted y trajectory

% Plot results
dt = 0.1; % Time step size (in seconds)
figure
lw_koop = 2; % Line width for plotting
subplot(1, 2, 1);
plot((0:n_traj) * dt, x_val_true, '-b', 'linewidth', lw_koop); hold on
plot((0:n_traj) * dt, x_val_pred, '--r', 'linewidth', lw_koop)
title('$x$', 'interpreter', 'latex', 'fontsize', 20);
xlabel('time (s)', 'interpreter', 'latex', 'fontsize', 14)
LEG = legend('True', 'Koopman');
set(LEG, 'Interpreter', 'latex', 'location', 'northeast', 'fontsize', 14)
set(gca, 'FontSize', 14);

subplot(1, 2, 2);
plot((0:n_traj) * dt, y_val_true, '-b', 'linewidth', lw_koop); hold on
plot((0:n_traj) * dt, y_val_pred, '--r', 'linewidth', lw_koop)
title('$y$', 'interpreter', 'latex', 'fontsize', 20);
xlabel('time (s)', 'interpreter', 'latex', 'fontsize', 14)
LEG = legend('True', 'Koopman');
set(LEG, 'Interpreter', 'latex', 'location', 'northeast', 'fontsize', 14)
set(gca, 'FontSize', 16);

%%
% controb = [Blift,Alift*Blift,Alift^2*Blift];
% %x_test_3 = Alift^3*Xlift(:,1) + controb*[U(:,1);U(:,2);U(:,3)];
% %x_test_2 = Alift^2*Xlift(:,1) + controb*[U(:,1);U(:,2);zeros(5,1)];
% x_test_1 = Alift*Xlift(:,1) + controb*[U(:,1);zeros(5,1);zeros(5,1)];
%%

function custom_observables = generate_observables(maxDegree)
    % Generate all monomials up to maxDegree, excluding degree 1 terms (x, y)
    
    terms = {'ones(size(x))'};  % Include constant term (degree 0)
    
    % Loop through degrees from 2 to maxDegree
    for total_degree = 2:maxDegree
        for x_power = 0:total_degree
            y_power = total_degree - x_power;
            
            % Exclude linear terms (x^1, y^1)
            if (x_power == 1 && y_power == 0) || (x_power == 0 && y_power == 1)
                continue;
            end
            
            % Generate the monomial term
            if x_power == 0
                term = sprintf('y.^%d', y_power);
            elseif y_power == 0
                term = sprintf('x.^%d', x_power);
            else
                term = sprintf('x.^%d .* y.^%d', x_power, y_power);
            end
            
            % Add the term to the list
            terms{end+1} = term;
        end
    end
    
    % Combine all terms into a single anonymous function
    expr = strjoin(terms, ';');
    custom_observables = eval(['@(x, y) [', expr, ']']);
end

function cellArray = generateMNArray(degree)
    % Initialize empty cell array
    total_elements = 0;
    
    % Calculate total number of elements needed
    for d = 0:degree
        total_elements = total_elements + (d + 1);
    end
    
    cellArray = cell(1, total_elements);
    idx = 1;
    
    % Generate elements for each degree level
    for d = 0:degree
        % Get all combinations of m and n that sum to current degree
        combinations = [];
        for m = 0:d
            n = d - m;
            combinations = [combinations; m n];
        end
        
        % Sort combinations based on the specified pattern
        % For degree d, order should be: (d,0), (d-1,1), ..., (0,d)
        [~, order] = sort(combinations(:,1), 'descend');
        combinations = combinations(order,:);
        
        % Create structs and store in cell array
        for i = 1:size(combinations, 1)
            s = struct('m', combinations(i,1), 'n', combinations(i,2));
            cellArray{idx} = s;
            idx = idx + 1;
        end
    end
end

function vals = legendre_observables(x, y, mn_vals, x_bar, y_bar)
    % Precompute normalization factor
    norm_factor = 1 / sqrt(x_bar * y_bar);
    
    % Precompute scaled x and y
    x_scaled = x / x_bar;
    y_scaled = y / y_bar;
    
    % Initialize the output array
    num_basis = length(mn_vals);
    vals = zeros(num_basis, numel(x));
    
    % Vectorized computation of Legendre polynomials
    for i = 1:num_basis
        m = mn_vals{i}.m;
        n = mn_vals{i}.n;
        
        % Compute normalized Legendre polynomials for m and n
        P_m = norm_legendre(m, x_scaled);
        P_n = norm_legendre(n, y_scaled);
        
        % Multiply and scale by the normalization factor
        vals(i, :) = norm_factor * P_m .* P_n;
        %vals(i, :) = 1 * P_m .* P_n;
    end
end

function P_m = norm_legendre(m, x)
    % Vectorized computation of normalized Legendre polynomial
    P = legendre(m, x);
    P_m = sqrt((2 * m + 1) / 2) * P(1, :); % Extract the first row (corresponding to m-th order)
end

function traj_data = generate_pendulum_trajectories(Nsim, Ntraj, m, l, b, g, dt)
    traj_data = cell(1, Nsim);
    
    for i = 1:Nsim
        traj = struct('x', [], 'u', [], 'y', []);
        
        % Initial state: random small perturbation around stable point
        x0 = [pi/4 * (rand); 0.4*rand-0.2]; 
        %x0 = [0; 0]; 
        
        for j = 1:Ntraj
            % Generate random control input u in range [-1, 1]
           % u = 0.1*rand-0.05;
            u = 0.5*rand-0.25;
            
            % Solve dynamics using ode45 over a single time step
            [~, x_next] = ode45(@(t, x) pendulum_dynamics(t, x, u, m, l, b, g), [0, dt], x0);
            x_next = x_next(end, :)'; % Take final state after integration
            
            % Convert theta to Cartesian coordinates
            % x_cart = l * [sin(x0(1)); -cos(x0(1));cos(x0(1)).*x0(2); sin(x0(1)).*x0(2)];
            % y_cart = l * [sin(x_next(1)); -cos(x_next(1));cos(x_next(1)).*x_next(2); sin(x_next(1)).*x_next(2)];

            x_cart = l * [sin(x0(1));cos(x0(1)).*x0(2)];
            y_cart = l * [sin(x_next(1));cos(x_next(1)).*x_next(2)];

           % x_cart = l * [sin(x0(1));-cos(x0(1))];
            %y_cart = l * [sin(x_next(1));-cos(x_next(1))];

            traj.x(:, j) = x0;
            traj.u(j) = u;
            traj.y(:, j) = x_next;
            
            % % Store data
            % traj.x(:, j) = x_cart;
            % traj.u(j) = u;
            % traj.y(:, j) = y_cart;
            
            % Update state for next iteration
            x0 = x_next;
        end
        
        traj_data{i} = traj;
    end
end

function dxdt = pendulum_dynamics(~, x, u, m, l, b, g)
    dxdt = [x(2);
            -(g/l) * sin(x(1)) - (b/(m*l^2)) * x(2) + (1/(m*l^2)) * u];
end
