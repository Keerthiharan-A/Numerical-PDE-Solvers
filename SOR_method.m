clc
clear
L = 24;  % Length of the plate in cm
% Taking del x = del y = del = 0.5cm
del = 0.5;
N = L/del; % Nx = Ny = N
init = zeros(N+1,N+1);
u_exact = zeros(N+1,N+1);
init(N+1,:) = 20;  % Top boundary
u_exact(N+1,:) = 20;
for i = 2:N
    for j = 2:N
        x = (j-1)*del;
        y = (i-1)*del;
        u_exact(i,j) = (80/pi)*inf_sum(x,y);
    end
end
w = linspace(1,1.95,20);  % Range of omega values

iterations = zeros(size(w));

for k = 1:length(w)
    omega = w(k);
    u_new = init;
    u = init;
    % Take 1e-8 to be the maximum error tolerance 
    tol = 1e-8;
    error = tol + 1;  % Initialize error to enter the loop
    iter = 0;
    while error >= tol
        for i = 2:N
            for j = 2:N
                u_new(i,j) = (1 - omega)*u(i,j) + (omega/4)*(u_new(i-1,j) + u(i+1,j) + u_new(i,j-1) + u(i,j+1));
            end
        end
        error = norm(abs(u_new - u));
        u = u_new;
        iter = iter + 1;
    end
    
    iterations(k) = iter;
end
% Create a plot of omega versus number of iterations
figure;
plot(w, iterations, '-o');
title('Number of Iterations for convergence vs. Overrelaxation Parameter (\omega)');
xlabel('\omega');
ylabel('Number of Iterations');
grid on;

% Error between successive iteration is 9.551076420699293e-09
absolute_error = abs(u_exact - u);

% Display the maximum absolute error
absolute_error_norm = norm(absolute_error);
fprintf('Absolute error Norm: %.6f\n', absolute_error_norm);
% Absolute error Norm: 0.349822
% Plot the surface of the numerical solution
figure;
surf((0:N)*del, (0:N)*del, u);
title('Surface Plot of Numerical Solution using SOR method (u)');
xlabel('X');
ylabel('Y');
zlabel('Temperature');
xlim([0, L]);
ylim([0, L]);

% Plot the surface of the absolute error
figure;
surf((0:N)*del, (0:N)*del, absolute_error);
title('Surface Plot of Absolute Error');
xlabel('X');
ylabel('Y');
zlabel('Absolute Error');
xlim([0, L]);
ylim([0, L]);

function val = inf_sum(x,y)
    val = 0;
    n = 1;
    while n < 100
        term = (sin(n*pi*x/24)*sinh(n*pi*y/24)/(n*sinh(n*pi)));
        val = val + term;
        n = n + 2;
    end
end
