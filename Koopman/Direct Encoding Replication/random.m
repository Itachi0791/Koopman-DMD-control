% Define symbolic variables
syms x k n

% Define the piecewise function F_x(x)
F = piecewise(0 <= x & x <= 0.3, 0.95 - 3*x, x > 0.3, (33/35)*x + 1/140);

% Define the function y(x)
y = exp(1i*2*pi*(k*F - n*x));

% Compute the integral of y(x) over the interval [0, 1]
integral_result = int(y, x, 0, 1);

% Simplify the result (optional)
integral_result = simplify(integral_result);

% Display the result
disp('The integral of y(x) from 0 to 1 is:');
disp(integral_result);
