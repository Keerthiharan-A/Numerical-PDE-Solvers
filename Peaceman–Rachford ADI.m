clc
clear
% Defining range and step size
n = 10;
x = linspace(0,1,n+1);
y = linspace(0,1,n+1);
t = linspace(0,0.2,n*n+1);
del_t = 0.002;
del_x = 0.1;
del_y = 0.1;

% Exact solution (Analytical soln)
exact_soln = zeros(n+1,n+1,n*n+1);
for k = 1:n*n+1
    for j = 1:n+1
        for i = 1:n+1
            exact_soln(i,j,k) = exp(-2*pi^2*t(k))*sin(pi*x(i))*sin(pi*y(j));
        end
    end
end

% Peaceman-Rachford ADI method
% Equation is Ut = (Uxx + Uyy)
rx = del_t / (2*(del_x)^2);
ry = del_t / (2*(del_y)^2);
u = zeros(n+1,n+1,n+1);

% Finding the boundary values
for j = 2:n
    for i = 2:n
        u(i,j,1) = sin(pi*x(i))*sin(pi*y(j));
    end
end

for k = 1:n*n
    % Step 1
    a(1) = 0;
    b(1) = 1+2*rx;
    c(1) = -rx;
    for i = 2:n-2
        a(i) = -rx;
        b(i) = (1+2*rx);
        c(i) = -rx;
    end
    a(n-1) = -rx;
    c(n-1) = 0;
    b(n-1) = 1+2*rx;
    u_half = zeros(n+1,n+1);
    for i = 2:n
        for j = 2:n
            d(j-1) = ry*u(i,j-1,k) + (1-2*ry)*u(i,j,k) + ry*u(i,j+1,k);
        end
        u_half(i,2:n) = Thomas_algorithm(a,b,c,d);
    end
    % Step 2
    b(1) = 1+2*ry;
    c(1) = -ry;
    for i = 2:n-2
        a(i) = -ry;
        b(i) = (1+2*ry);
        c(i) = -ry;
    end
    a(n-1) = -ry;
    b(n-1) = 1+2*ry;
    for j = 2:n
        for i = 2:n
            d(i-1) = rx*u_half(i-1,j) + (1-2*rx)*u_half(i,j) + rx*u_half(i+1,j);
        end
        u(2:n,j,k+1) = Thomas_algorithm(a,b,c,d);
    end
end

%%% Comparing the analytical and numerically obtained solution


% Compute error along each column for every time step excluding the first and last rows
error_columns = zeros(n - 1, length(t) - 1);
for k = 2:length(t)
    for i = 2:n
        error_column_sum = 0;
        for j = 2:n
            error_column_sum = error_column_sum + (exact_soln(i, j, k) - u(i, j, k))^2;
        end
        error_columns(i - 1, k - 1) = sqrt(error_column_sum / (n - 1));
    end
end

% As expected the solution and error are symmetric about x and y axis also about x = 0.5 

% Plot error along each column for every time step (it will be same for along rows too)
figure;
for i = 2:n 
    subplot(ceil((n - 1) / 4), 4, i - 1);
    plot(t(2:end), error_columns(i - 1, :), 'b-', 'LineWidth', 1.5); % Adjusted time array
    title(['RMSE along x =  ', num2str((i-1)*del_x)]);
    xlabel('Time');
    ylabel('Error');
    grid on;
end
% Plot the differrence between exact and numerical solutions
figure;
for i = 2:2:10
    exact_solution = exact_soln(:, :, 1+i*10);
    numerical_solution = u(:, :, 1+i*10);
    % Plot side by side
    error = abs(exact_soln - u);
    subplot(2, 3, i/2);
    surf(x, y, error(:,:,1+i*10));
    title(['Absolute error at t = ', num2str(t(1+i*10))]);
    xlabel('x');
    ylabel('y');
    zlabel('u');
    axis([0 1 0 1 -1 1]);
    colormap('jet');
    shading interp;
    colorbar;
end

% Using Thomas algorithm for solving the tri-diagonal system
function soln = Thomas_algorithm(a,b,c,d)
    l = length(a);
    for i = 2:l
        b(i) = b(i) - a(i)*c(i-1)/b(i-1);
        d(i) = d(i) - a(i)*d(i-1)/b(i-1);
    end
    soln = zeros(1,l);
    soln(l) = d(l)/b(l);
    for i = l-1:-1:1
        soln(i) = (d(i)-c(i)*soln(i+1))/b(i);
    end
end
