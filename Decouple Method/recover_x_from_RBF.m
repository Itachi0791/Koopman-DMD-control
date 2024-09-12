function x_sol = recover_x_from_RBF(rbf_values, centers)
    % Inputs:
    %   rbf_values: an mx1 vector of known RBF values
    %   centers: an mx2 matrix where each row is a center c_i for RBF i
    sigma = 1; % Gaussian width (standard deviation)
    
    % Initial guess for x
    x0 = [0; 0];  % You can modify the initial guess if needed

    % Define the system of equations as a function handle
    rbf_system = @(x) exp(-sum((x' - centers).^2, 2) / (2 * sigma^2)) - rbf_values;

    % Specify options to use the Levenberg-Marquardt algorithm explicitly
    options = optimoptions('fsolve', 'Algorithm', 'levenberg-marquardt', 'Display', 'off'); 

    % Solve the system using fsolve
    x_sol = fsolve(rbf_system, x0, options);
end

