function ex1_2()
    %h([1 2 3], [4 5 6])
    %grad_descent(0.001)
    ex2_4()
end

function ex_2_1()
% Surface plot of h
    x = -2:0.05:2;
    y = -1:0.05:3;
    h = h(x, y');
    figure()
    surf(x, y, h);
end

function ex2_4()
    % Contour plot of gradient descent on h(x, y)
    x = -2:0.005:2;
    y = -1:0.005:3;
    H = h(x, y');
    figure()
    contour(x, y, H, 50)
    hold on 
    % Gradient descent plots for different eta's
    [x, y, h] = grad_descent(0.001);
    plot(x, y)
    xlim([-2 2])
    ylim([-1 3])
    hold off
end

function [x, y, H] = grad_descent(eta)

    % Convergence threshold
    precision = 0.00001;
    previous_step_size = 1;
    max_iterations = 1/precision; 

    % Starting point
    x = -1.75; 
    y = -0.5; 
    H(1) = h(x, y);
    i = 2;
    while previous_step_size > precision && i < max_iterations
       x = grad_x(x, y, eta);
       y = grad_y(x, y, eta);
       %x(i) = x(i-1) - eta * (-400 * x(i-1) * y(i-1) + 400 * x(i-1).^3 + 2 * x(i-1) - 2);
       %y(i) = y(i-1) - eta * (200 * y(i-1) - 200 * x(i-1).^2);
       H(i) = h(x(i), y(i));
       previous_step_size = abs(x(i) - x(i-1)) + abs(y(i) - y(i-1)); 
       i = i + 1;
    end
end

function x = grad_x(x, y, eta)
    n = length(x);
    x(n+1) = x(n) - eta * (-400 * x(n) * y(n) + 400 * x(n).^3 + 2 * x(n) - 2);
end

function y = grad_y(x, y, eta)
    n = length(y);
    y(n+1) = y(n) - eta * (200 * y(n) - 200 * x(n).^2);
end


function h = h(x, y)
    h = 100 * (y - x.^2).^2 + (1 - x).^2;
end

    
    