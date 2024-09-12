function A = lin_system(delt, centers)
% Input : delt is time step, centers of gaussian RBFs
% Output : A is m x m , where m is no. of Gaussian RBFs (no. of centers)
    m = size(centers,1);
    Q = zeros(m,m); R = zeros(m,m);
    for i = 1:m
        for j = 1:m
            c1 = (centers(i,:))';  c2 = (centers(j,:))'; 
            [Q(i,j),R(i,j)] =  rbf_inner_product(delt, c1, c2);
        end
    end
    A = Q*pinv(R);
end

function [Q, R] = rbf_inner_product(delt, c1, c2)
    % Function to compute the inner product of two RBF kernels numerically
    % Inputs:
    %   delt: Time step
    %   c1: 2x1 vector for center of the first RBF kernel [c1; c2]
    %   c2: 2x1 vector for center of the second RBF kernel [c3; c4]
    % Output:
    %   Q: Result of first RBF kernel inner product
    %   R: Result of second RBF kernel inner product
    
    % Define numerical values for the parameters
    sigma = 1;  % Example value for sigma, adjust if needed
    
    % Define the transformation F(x) = (I + A*delta)*x
    A = [0 1; -2 -3]; % Given matrix A
    I = eye(2);       % Identity matrix
    
    % Define the product of RBF kernels using element-wise operations
    product_fn1 = @(x1, x2) exp(-((([x1; x2]' * (I + A * delt)') - c1') * ...
                           ((I + A * delt) * [x1; x2] - c1)) / (2 * sigma^2)) ...
                         .* exp(-(([x1; x2] - c2)' * ([x1; x2] - c2)) / (2 * sigma^2));
    
    product_fn2 = @(x1, x2) exp(-(([x1; x2] - c1)' * ([x1; x2] - c1)) / (2 * sigma^2)) ...
                         .* exp(-(([x1; x2] - c2)' * ([x1; x2] - c2)) / (2 * sigma^2));

    % Compute the double integral using integral2 for numerical evaluation
    Q = integral2(@(x1, x2) arrayfun(product_fn1, x1, x2), -10, 10, -10, 10); 
    R = integral2(@(x1, x2) arrayfun(product_fn2, x1, x2), -10, 10, -10, 10);
end

