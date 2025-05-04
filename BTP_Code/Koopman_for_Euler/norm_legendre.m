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