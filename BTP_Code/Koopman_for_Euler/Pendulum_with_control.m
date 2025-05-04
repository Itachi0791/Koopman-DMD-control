clc;clear;close all;
gamma=0;g=5;l=1;
% Pendulum Equation : xdot = y, ydot = -gamma*y -g/l*sin(x)
f1 = @(x,y) y;
f2 = @(x,y) -gamma*y-g/l*sin(x);
f = @(x,y) [f1(x,y);f2(x,y)];
F1 = @(x,y,delt) x + delt*f1(x,y);
F2 = @(x,y,delt) y + delt*f2(x,y);
F = @(x,y,delt) [x;y] + delt*f(x,y);
delt = 0.01; degree=5;
mn_vals = generateMNArray(degree);
num_basis = length(mn_vals);
phi = cell(1,num_basis);
x_bar = pi; y_bar = pi;
for i = 1:num_basis
    phi{i} = @(x,y) 1/sqrt(x_bar*y_bar)* norm_legendre(mn_vals{i}.m,x/x_bar).*norm_legendre(mn_vals{i}.n,y/y_bar);
end
f1_coeffs = zeros(num_basis,1); f2_coeffs = zeros(num_basis,1);
for i=1:num_basis
    f1_coeffs(i) = integral2(@(x, y) f1(x, y).*phi{i}(x, y), -x_bar, x_bar, -y_bar, y_bar);
    f2_coeffs(i) = integral2(@(x, y) f2(x, y).*phi{i}(x, y), -x_bar, x_bar, -y_bar, y_bar);
end

%% A_bar computation

Ax = zeros(1,num_basis);Ay = zeros(1,num_basis);
for i = 1:num_basis
    Ax(i) = (i==2)/phi{2}(1,1) + delt*f1_coeffs(i);
    Ay(i) = (i==3)/phi{3}(1,1) + delt*f2_coeffs(i);
end
A = [Ax;Ay]; % A gives [x,y] from Phi_x

% %% Prediction from Simulation
% %rng(1)
% close all;
% Ns = 500;
% x0 = [pi/4;0];%[x_bar*rand-x_bar/2;2*y_bar*rand-y_bar];
% phi_x = zeros(num_basis,Ns);
% for i=1:num_basis
%     phi_x(i,1) = phi{i}(x0(1),x0(2));
% end
% x_euler = zeros(2,Ns); x_euler(:,1)=x0;x_pred=zeros(2,Ns); x_pred(:,1) = x0;
% [~, x] = ode45(@(t,x)[x(2); -gamma*x(2) - (g/l) * sin(x(1))], 0:delt:(Ns-1)*delt, x0);
% x=x';
% C = zeros(2,num_basis); C(1,2) = 1/phi{2}(1,1); C(2,3) = 1/phi{3}(1,1); % Mapping from phi to [x,y]
% for i=1:Ns-1
%     x_euler(:,i+1) = F(x_euler(1,i),x_euler(2,i),delt);
%     %x_pred(:,i+1) = C*phi_x(:,i+1);
%     x_pred(:,i+1) = A*phi_x(:,i);
%     for j=1:num_basis
%         phi_x(j,i+1) = phi{j}(x_pred(1,i+1),x_pred(2,i+1));
%         phi_x(j,i+1) = phi{j}(x_pred(1,i+1),x_pred(2,i+1));
%     end
% 
% end
% t = 0:delt:(Ns-1)*delt;
% 
% figure;
% sgtitle("Koopman Prediction(Two Rows)", "FontSize", 22, "interpreter", "latex", "FontWeight", "bold");
% set(gcf,"WindowState","maximized")
% subplot(2,1,1) 
% plot(t,x(1,:),"-g","LineWidth",3)
% hold on
% grid on
% plot(t,x_euler(1,:),"-b","LineWidth",3)
% plot(t,x_pred(1,:),"-ro","LineWidth",0.5)
% xlabel("Time(t)","FontSize",20,"interpreter","latex", "FontWeight", "bold");
% ylabel("$\theta(t)$","FontSize",25,"interpreter","latex", "FontWeight", "bold","Rotation",0)
% legend("Actual","Euler","Predicted (using $\bar{A}$)","interpreter","latex","FontSize",15)
% set(gca,"FontSize",20)
% subplot(2,1,2)
% plot(t,x(2,:),"-g","LineWidth",3)
% hold on
% plot(t,x_euler(2,:),"-b","LineWidth",3)
% grid on
% plot(t,x_pred(2,:),"-ro","LineWidth",0.5)
% xlabel("Time(t)","FontSize",20,"interpreter","latex", "FontWeight", "bold");
% ylabel("$\dot{\theta}(t)$","FontSize",25,"interpreter","latex", "FontWeight", "bold","Rotation",0)
% legend("Actual","Euler","Predicted (using $\bar{A}$)","interpreter","latex","FontSize",15)
% set(gca,"FontSize",18)

%% Prediction with Control Input
B = [0;1];
B_lift = delt*B;%[1/phi{2}(1,1),0;0,1/phi{3}(1,1)]*B;
f_u = @(x,y,u) [f1(x,y);f2(x,y)] + B*u;
Ns = 200; x0 = [pi/2*rand-pi/4;2*rand-1];
x = zeros(2,Ns+1); x(:,1) = x0; x_pred = x; u = 30*rand(1,Ns)-15;
phi_x = zeros(num_basis,Ns+1);
for i=1:num_basis
    phi_x(i,1) = phi{i}(x0(1),x0(2));
end
for i = 1:Ns 
    f_sys = @(t,x) [x(2); -gamma*x(2) - (g/l) * sin(x(1))] + B*u(i);
    [~,y] = ode45(f_sys,[0,delt],x(:,i));
    x(:,i+1) = (y(end,:))';
    x_pred(:,i+1) = A*phi_x(:,i) + B_lift*u(i);
    for j=1:num_basis
        phi_x(j,i+1) = phi{j}(x_pred(1,i+1),x_pred(2,i+1));
        phi_x(j,i+1) = phi{j}(x_pred(1,i+1),x_pred(2,i+1));
    end
end
t = 0:delt:Ns*delt;
figure;
sgtitle("Koopman Prediction($\Delta t =0.01$, 21 observables)", "FontSize", 22, "interpreter", "latex", "FontWeight", "bold");
set(gcf,"WindowState","maximized")
subplot(2,1,1) 
plot(t,x(1,:),"-g","LineWidth",3)
hold on
grid on
%plot(t,x_euler(1,:),"-b","LineWidth",3)
plot(t,x_pred(1,:),"-ro","LineWidth",0.5)
xlabel("Time(t)","FontSize",25,"interpreter","latex", "FontWeight", "bold");
ylabel("$\theta(t)$","FontSize",25,"interpreter","latex", "FontWeight", "bold","Rotation",0)
legend("ODE45","Predicted (using $\bar{A}$)","interpreter","latex","FontSize",15)
set(gca,"FontSize",20)
subplot(2,1,2)
plot(t,x(2,:),"-g","LineWidth",3)
hold on
%plot(t,x_euler(2,:),"-b","LineWidth",3)
grid on
plot(t,x_pred(2,:),"-ro","LineWidth",0.5)
xlabel("Time(t)","FontSize",25,"interpreter","latex", "FontWeight", "bold");
ylabel("$\dot{\theta}(t)$","FontSize",25,"interpreter","latex", "FontWeight", "bold","Rotation",0)
legend("ODE45","Predicted (using $\bar{A}$)","interpreter","latex","FontSize",15)
set(gca,"FontSize",18)

%% Functions
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

function val = in_prod(k,l,j)
    % computes < phi_k(x)*phi_l(x),phi_j(x) >, where phi is norm_legendre
    val = sqrt((2*j+1)*(2*k+1)*(2*l+1)/2)*(Wigner3j([k,l,j],[0 0 0]))^2;
end

function val = deriv_in_prod(i,k,j)
    val = 0;
     % computes < phi_i_prime(x)*phi_k(x),phi_j(x) >, where phi is norm_legendre
    for q = 1:floor((i+1)/2)
        val = val + sqrt(2*i-4*q+3)*in_prod(i-2*q+1,k,j);
    end
    val = sqrt(2*i+1)*val;
end

function P_m = norm_legendre(m, x)
    % This function returns the normalized m-th Legendre polynomial at x.
    % Inputs:
    %   m - The order of the Legendre polynomial (integer)
    %   x - The value(s) at which to evaluate the polynomial (real numbers between -1 and 1)
    % Output:
    %   P_m - The value(s) of the normalized m-th Legendre polynomial at x (same size as x)
    
    % Check if x contains values outside the valid range
    % if any(x(:) < -1) || any(x(:) > 1)
    %     error('All elements of x must be between -1 and 1.');
    % end
    
    % Compute the associated Legendre function for each element of x
    % and extract the value corresponding to m
    size_x = size(x); % Store the original size of x
    x_flat = x(:);    % Flatten x for easier iteration
    P_m_flat = zeros(size(x_flat)); % Preallocate the output
    
    for i = 1:numel(x_flat)
        P = legendre(m, x_flat(i));
        P_m_flat(i) = sqrt((2*m+1)/2)*P(1);
    end
    
    % Reshape the output to match the input size
    P_m = reshape(P_m_flat, size_x);
end