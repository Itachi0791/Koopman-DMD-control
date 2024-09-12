function A = A_negx(m,delt)
% A_negx : Computes state evolution matrix A for the lifted space
% input : m(odd) - No. of observables, delt - Time step for discrete map
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
            val = (n - k + delt*k);
            if k == n
                Am_bar(i,j) = 1 - pi*1i*k*delt;
            else
                p = 2*pi*1i*(k-n);
                Am_bar(i,j) = (1/p^2 + (exp(p)*(p - 1))/p^2)*2*pi*1i*delt*k;
            end           
            % if val == 0
            %     Am_bar(i,j) = 1 ;
            % else
            %     t1 = log(exp(-2*pi*(- k*1i))); t2 = log(exp(-2*pi*(n*1i))); t3 = log(exp(-2*pi*(delt*k*1i)));
            %     Am_bar(i,j) = (exp(t1)*exp(t2)*exp(t3)*1i - 1i)/(2*pi*(n - k + delt*k));
            %     %Am_bar(i,j) = (exp(-2*pi*(- k*1i + n*1i + delt*k*1i))*1i - 1i)/(2*pi*(n - k + delt*k));
            % end
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