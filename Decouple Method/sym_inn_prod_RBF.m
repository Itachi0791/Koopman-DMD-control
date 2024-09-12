syms x1 x2 delta c1 c2 c3 c4 sigma real

% Define the symbolic variables
x = [x1; x2];
A = [0 1; -2 -3];

% Define the transformation F(x) = (I + A*delta)*x
I = eye(2);
F_x = (I + A*delta) * x;

% Define the RBF kernel function gk(x)
gk = exp(-norm(x - [c1; c2])^2 / (2 * sigma^2));

% Define the RBF kernel function gn(x)
gn = exp(-norm(x - [c3; c4])^2 / (2 * sigma^2));

% Compute the inner product of gk(F(x)) and gn(x) over the domain X
inner_product = int(int(gk * gn, x1, -10, 10), x2, -10, 10);

% Simplify the result
inner_product_simplified = simplify(inner_product);

% Display the result
disp(inner_product_simplified);
