clc;clear;close all;
mu=1;
% Vanderpol Equation : xdot = y, ydot = -x +mu*(1-x^2)*y
f1 = @(x,y) y;
f2 = @(x,y) -x + mu*(1-x.^2).*y;
f = @(x,y) [f1(x,y);f2(x,y)];
F1 = @(x,y,delt) x + delt*f1(x,y);
F2 = @(x,y,delt) y + delt*f2(x,y);
F = @(x,y,delt) [x;y] + delt*f(x,y);
delt = 0.02; degree=5;
mn_vals = generateMNArray(degree);
num_basis = length(mn_vals);
phi = cell(1,num_basis);
x_bar = 6; y_bar = 6;
for i = 1:num_basis
    phi{i} = @(x,y) 1/sqrt(x_bar*y_bar)* norm_legendre(mn_vals{i}.m,x/x_bar).*norm_legendre(mn_vals{i}.n,y/y_bar);
end
f1_coeffs = zeros(num_basis,1); f2_coeffs = zeros(num_basis,1);
for i=1:num_basis
    f1_coeffs(i) = integral2(@(x, y) f1(x, y).*phi{i}(x, y), -x_bar, x_bar, -y_bar, y_bar);
    f2_coeffs(i) = integral2(@(x, y) f2(x, y).*phi{i}(x, y), -x_bar, x_bar, -y_bar, y_bar);
end

%% A_bar computation

A_bar = zeros(num_basis, num_basis); A_bar_est = A_bar;
Ax = zeros(1,num_basis);Ay = zeros(1,num_basis);
for i = 1:num_basis
    for j = 1:num_basis
        %integrand = @(x,y) phi{i}(F1(x,y,delt),F2(x,y,delt)).*phi{j}(x,y);
        %A_bar(i, j) = integral2(integrand, -x_bar, x_bar, -y_bar, y_bar);
        sum = 0;
        mi = mn_vals{i}.m;ni = mn_vals{i}.n;mj = mn_vals{j}.m;nj = mn_vals{j}.n;
        for k = 1:num_basis
            mk = mn_vals{k}.m;nk = mn_vals{k}.n;
            sum = sum + y_bar*f1_coeffs(k)*deriv_in_prod(mi,mk,mj)*in_prod(ni,nk,nj)...
                      + x_bar*f2_coeffs(k)*deriv_in_prod(ni,nk,nj)*in_prod(mi,mk,mj);
        end
        A_bar_est(i,j) = (i==j) + delt/((x_bar*y_bar)^1.5)*sum;
    end
    Ax(i) = (i==2) + delt*phi{2}(1,1)*f1_coeffs(i);
    Ay(i) = (i==3) + delt*phi{3}(1,1)*f2_coeffs(i);
end

%% Prediction from Simulation
%rng(1)
Ns = 1000;
x0 = [0.2;-0.1];%[2*x_bar/2*rand-x_bar/2;2*y_bar/2*rand-y_bar/2];
phi_x = zeros(num_basis,Ns);
for i=1:num_basis
    phi_x(i,1) = phi{i}(x0(1),x0(2));
end
phi_x_2 = phi_x;
x_euler = zeros(2,Ns); x_euler(:,1)=x0;x_pred1=zeros(2,Ns); x_pred1(:,1) = x0;x_pred2 = x_pred1;
[~, x] = ode45(@(t,x)[x(2);-x(1) + mu*(1-x(1).^2).*x(2)], 0:delt:(Ns-1)*delt, x0);
x=x';
C = zeros(2,num_basis); C(1,2) = 1/phi{2}(1,1); C(2,3) = 1/phi{3}(1,1); % Mapping from phi_x to [x,y]
for i=1:Ns-1
    x_euler(:,i+1) = F(x_euler(1,i),x_euler(2,i),delt);
    phi_x(:,i+1) = A_bar_est*phi_x(:,i);
    x_pred1(:,i+1) = C*phi_x(:,i+1);
    x_pred2(:,i+1) = [C(1,2)*Ax;C(2,3)*Ay]*phi_x_2(:,i);
    for j=1:num_basis
        phi_x(j,i+1) = phi{j}(x_pred1(1,i+1),x_pred1(2,i+1));
        phi_x_2(j,i+1) = phi{j}(x_pred2(1,i+1),x_pred2(2,i+1));
    end

end
t = 0:delt:(Ns-1)*delt;
% figure;
% sgtitle("Koopman Prediction(Full Matrix)", "FontSize", 22, "interpreter", "latex", "FontWeight", "bold");
% set(gcf,"WindowState","maximized")
% subplot(2,1,1)
% plot(t,x(1,:),"-b","LineWidth",3)
% hold on
% grid on
% plot(t,x_pred(1,:),"-ro","LineWidth",1)
% xlabel("Time(t)","FontSize",20,"interpreter","latex", "FontWeight", "bold");
% ylabel("$\theta(t)$","FontSize",25,"interpreter","latex", "FontWeight", "bold","Rotation",0)
% legend("Actual","Predicted (using $\bar{A}$)","interpreter","latex","FontSize",15)
% set(gca,"FontSize",20)
% subplot(2,1,2)
% plot(t,x(2,:),"-b","LineWidth",3)
% hold on
% grid on
% plot(t,x_pred(2,:),"-ro","LineWidth",1)
% xlabel("Time(t)","FontSize",20,"interpreter","latex", "FontWeight", "bold");
% ylabel("$\dot{\theta}(t)$","FontSize",25,"interpreter","latex", "FontWeight", "bold","Rotation",0)
% legend("Actual","Predicted (using $\bar{A}$)","interpreter","latex","FontSize",15)
% set(gca,"FontSize",18)

figure;
sgtitle("Koopman Prediction(Wide Matrix, $\Delta t = 0.02$, 21 observables)", "FontSize", 22, "interpreter", "latex", "FontWeight", "bold");
set(gcf,"WindowState","maximized")
plot(x(1,:),x(2,:),"-g","LineWidth",3)
hold on
plot(x_euler(1,:),x_euler(2,:),"-b","LineWidth",3)
grid on
plot(x_pred2(1,:),x_pred2(2,:),"-ro","LineWidth",1)
xlabel("x(t)","FontSize",20,"interpreter","latex", "FontWeight", "bold");
ylabel("$\dot{x}(t)$","FontSize",25,"interpreter","latex", "FontWeight", "bold","Rotation",0)
legend("ODE45","Euler","Predicted (using $\bar{A}$)","interpreter","latex","FontSize",20)
set(gca,"FontSize",20)


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

% integrand = @(x, y) phi{3}(x, y) .* phi{1}(x, y);
% 
% result = integral2(integrand, -x_bar, x_bar, -y_bar, y_bar);