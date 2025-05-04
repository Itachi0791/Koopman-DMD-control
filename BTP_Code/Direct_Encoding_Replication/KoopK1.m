 close all;clear;clc;
m=101;
n=200;
x=(linspace(0,1,n))';
Fx = F_x(x);
figure;
plot(x,Fx,"LineWidth",2.5)
xlabel("x","FontSize",20);
ylabel("F(x)","FontSize",20);
grid on;
set(gca,"FontSize",15)
%Fx = x;
%g = @(x) cos(2*pi*x);
g = @(x) x.^2;
K = KoopK(x,Fx,m);
goFx = g(Fx);
Kfog = K*g(x);
figure;
plot(x,goFx,"LineWidth",2.5);hold on;
scatter(x,Kfog,"magenta","filled")
xlabel("x","FontSize",20)
title("g(x) = x^2","FontSize",20)
legend("g(F(x))","K_F g(x)","FontSize",20)
grid on
set(gca,"FontSize",15)
hold off
% G = zeros(m,n);
% G(1,:) = 1; 
% for i = 2:m
%     for j =1:n
%         if mod(i,2) == 0 
%             G(i,j) = cos(pi*i*x(j));
%         else 
%             G(i,j) = sin(pi*(i-1)*x(j));
%         end
%     end
% end
% %G(1,:) = x';
% A = G*K'*pinv(G);
 
function K = KoopK(x,Fx,m)
% Function to find the Koopman operator for a non-linear system with F(x) governing evolution
% Inputs : x - (n x 1) state 
%          Fx - F(x) (n x 1) values
%          m - set of values k for the orthonormal trigonometric basis functions
% Output : K - Koopman operator for F(x), size (n x n)
    n = length(x)-1;
    band = (m-1)/2;
    K = zeros(n+1,n+1);
    for p = -band:band
        
        hp = exp(2*pi*1i*p*x);
        hFp = exp(2*pi*1i*p*Fx);
        K = K + hFp*hp' ;
        %hp'*cos(2*pi*x)
        %norm(hFp*hp')
    end
    K = K*(1/n); % delta xi in kernel matrix is 1/n here
    K = real(K);
end

