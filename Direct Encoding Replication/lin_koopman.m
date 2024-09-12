function A = lin_koopman(m)
% lin_koopman : Computes state evolution matrix A for the lifted space
% input : m(odd) - No. of observables
% Output : A (m x m) - Lifted state evolution matrix
    Am_bar = zeros(m,m); Cm = zeros(m,m); Cm_inv = zeros(m,m);
    n_max = (m - 1) / 2;
    positive_indices = (1:n_max)';          % [1, 2, 3, ..., n_max]'
    negative_indices = -(1:n_max)';         % [-1, -2, -3, ..., -n_max]' 
    basis_indices = zeros(m, 1);
    basis_indices(1:2:end) = [0; positive_indices]; % [0, 1, 2, ...]'
    basis_indices(2:2:end) = negative_indices;      % [-1, -2, ...]'

    for i = 1:m
        for j = 1:m
            k = basis_indices(i); n = basis_indices(j);
            
            if n+3*k == 0 || 33*k - 35*n == 0
                Am_bar(i,j) = integral(@(x) exp(1i*2*pi*(k*F_x(x)-n*x)), 0, 0.3, 'RelTol', 1e-6, 'AbsTol', 1e-10) + ...
                              integral(@(x) exp(1i*2*pi*(k*F_x(x)-n*x)), 0.3, 1, 'RelTol', 1e-6, 'AbsTol', 1e-10);
            else
                % Use logarithmic formulation to reduce overflow risk
                term1 = log(exp((pi*(k*1i - n*6i))/10) * 1i - exp((pi*k*19i)/10) * 1i) - log(2*pi*(3*k + n));
                term2 = log(35) + log(exp((pi*(k*29i - n*100i))/50) * (exp((pi*k*33i)/25) * 1i - exp((pi*n*7i)/5) * 1i)) - log(2*pi*(33*k - 35*n));
                Am_bar(i,j) = exp(term1) - exp(term2);
            end           
        end
    end

    % Build Cm and Cm_inv matrices
    Cm(1,1) = 1;
    for i = 2:2:m-1
        Cm(i, i) = 0.5; Cm(i, i+1) = 0.5;
        Cm(i+1, i) = -0.5i; Cm(i+1, i+1) = 0.5i;
    end

    Cm_inv(1,1) = 1;
    for i = 2:2:m-1
        Cm_inv(i, i) = 1; Cm_inv(i, i+1) = 1i;
        Cm_inv(i+1, i) = 1; Cm_inv(i+1, i+1) = -1i;
    end

    % Calculate A matrix
    A = Cm * Am_bar * Cm_inv;
    A = real(A);
end
