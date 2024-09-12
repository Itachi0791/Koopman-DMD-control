function x = recover_x(Xi)
% Recover_x : Retrives the state space vector from the lifted space of
% observables for different time instants
% Input: Xi - (m x tsteps), Row i ha observables g_i, column j is for time
% instant j
% Output : x - (tsteps x 1), state space for each time instant
    m = size(Xi,1);
    n_max = (m - 1) / 2;
    % Generate the indices
    positive_indices = (1:n_max)';          % [1, 2, 3, ..., n_max]'
    negative_indices = -(1:n_max)';         % [-1, -2, -3, ..., -n_max]' 
    % Interleave the positive and negative indices
    basis_indices = zeros(m, 1);
    basis_indices(1:2:end) = [0; positive_indices]; % [0, 1, 2, ...]'
    basis_indices(2:2:end) = negative_indices;      % [-1, -2, ...]'
    coeffs = zeros(m,1);Cm_inv = zeros(m,m);
    coeffs(1) = 0.5;
    for i = 2:m
        k = basis_indices(i);
        coeffs(i) = exp(-2*pi*1i*k)*(-2*pi*1i*k + exp(2*pi*1i*k)-1)/(4*pi^2*k^2);
    end
    Cm_inv(1,1) = 1;
    for i = 2:2:m-1
        Cm_inv(i, i) = 1;Cm_inv(i, i+1) = 1i;
        Cm_inv(i+1, i) = 1;Cm_inv(i+1, i+1) = -1i;
    end
    x = real((coeffs'*Cm_inv*Xi)');
end