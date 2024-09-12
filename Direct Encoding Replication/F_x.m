function F = F_x(x)
    F = zeros(size(x));
    idx1 = (0 <= x) & (x <= 0.3);
    F(idx1) = 0.95 - 3*x(idx1);
    idx2 = x > 0.3;
    F(idx2) = 33/35*x(idx2) + 1/140;
end