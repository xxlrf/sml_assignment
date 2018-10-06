function ex1_1()
    % Create training and test set
    training_set_10 = create_dataset(10);
    training_set_40 = create_dataset(10);
    test_set = create_dataset(100);

    %ex1_3(training_set_10, test_set)
    %ex1_3(training_set_40, test_set)
    ex1_5(training_set_10, test_set)
end

function ex1_3(training_data, test_data)
%     D_10 = create_dataset(10)
%     w = PolCurFit(D_10, 9)
%     fplot(@(x) y(x, w), [0,1], 'r')
%     hold on
%     scatter(linspace(0, 1, 10), D_10(2,:), 'b')
%     fplot(@(x) f(x), [0,1], 'g')
%     hold off
    W = zeros(9, 9);
    E = zeros(2, 9);
    for M = 0:9
        % Calculate w
        w = PolCurFit(training_data, M);
        W(1:length(w), M+1) = w;
        
        % RMS error
        E(1, M+1) = RMS(training_data, w);
        E(2, M+1) = RMS(test_data, w);
                
        % Plot y(x, w) if M in [0 1 3 9]
        if ismember(M,[0, 1, 3, 9])
            %Plot function f(x) and training data
            figure()
            fplot(@(x) f(x), [0,1], 'g') % Plot f(x) over interval [0 1]
            hold on
            scatter(linspace(0, 1, length(training_data)), training_data(2,:), 'b')
            fplot(@(x) y(x, w), [0,1], 'r')
            xlabel('x')
            ylabel('t')
            ylim([-0.5 2.5])
            num2str(M)
            title(sprintf('Exercise 1.3, M = %1i', M))
          
            hold off
        end
    end
    
    % Plot the RMS for M = 0:9
    figure()
    plot(0:9, E(1,:), '-ob')
    hold on
    plot(0:9, E(2,:), '-or')
    xlabel('M')
    ylabel('E_R_M_S')
    ylim([0 1])
    xlim([-0.5 9.5])
    legend('Training', 'Test')
    title("Exercise 1.4")
    hold off
end

function ex1_5(training_data, test_data)
    
    j = 1;
    for i = 40:-0.5:20
        lambda = 1/exp(i);
        w = PolCurFitLambda(training_data, 9, lambda)
        
%         % Plot curve
%         figure()
%         fplot(@(x) f(x), [0,1], 'g') % Plot f(x) over interval [0 1]
%         hold on
%         scatter(linspace(0, 1, length(training_data)), training_data(2,:), 'b')
%         fplot(@(x) y(x, w), [0,1], 'r')
%         xlabel('x')
%         ylabel('t')
%         ylim([-0.5 2.5])
% 
%         hold off
        
        % RMS error
        E(1, j) = RMS(training_data, w);
        E(2, j) = RMS(test_data, w);
        j = j + 1;
    end
        
    % Plot the RMS for M = 0:9
    figure()
    plot(-40:0.5:-20, E(1,:), '-b')
    hold on
    plot(-40:0.5:-20, E(2,:), '-r')
    xlabel('ln\lambda')
    ylabel('E_R_M_S')
    ylim([0 1])
    xlim([-39, -19])
    legend('Training', 'Test')
    title("Exercise 1.5")
    hold off
end

function y = f(x)
% F  Function f(x) = 1 + sin(6(x-2))
    y = 1 + sin(6*(x - 2));
end

function data = create_dataset(N)
% CREATE_DATASET  Creates a dataset of N input/output pairs
% 
%   data = set of (x; y) values where
%       x = input value x (N evenly spaced values between [0, 1]
%       y = f(x) + gaussian noise
    input = linspace(0, 1, N); % N input values of f spaced uniformly in range [0, 1]
    D_N = f(input);
    noise = normrnd(0, 0.3, [1,N]); % Noise taken from a Gaussian distribution with m = 1, s = 0.3
    data = [input; D_N + noise];
end

% % Manually loop over n for verification
% function w = PolCurFit(D_N, M)
% 
%     A = zeros(M+1, M+1)
%     T = zeros(1, M+1)
% 
%     for i = 0:M % Loop over w_i
%         
%         % Loop over n to get T_i
%         for n = 1:length(D_N)
%             T(i + 1) = T(i+1) + D_N(2,n) * D_N(1,n).^i            
%         end
%         
%         for j = 0:M % Loop over A_i_j
%             
%             % Loop over n to get A_i_j
%             for n = 1:length(D_N)
%                 A(i+1, j+1) = A(i+1, j+1) + D_N(1,n) .^ (i+j)                
%             end
%             
%         end
%     end
%     
%     w = A\T';
%     
% end

function w = PolCurFitLambda(D_N, M, lambda)
% POLCURFIT Fits a polynomial curve of dimension M to data given by D_N
%
%   Output w: column vector of length M with parameters
    A = zeros(M+1, M+1);
    T = zeros(1, M+1);

    for i = 0:M % Loop over w_i where i from [0, M]
        
        % Loop over n to get T_i
        for n = 1:length(D_N)
            T(i + 1) = T(i+1) + D_N(2,n) * D_N(1,n).^i;            
        end
        
        for j = 0:M % Loop over A_i_j
            
            % Loop over n to get A_i_j
            for n = 1:length(D_N)
                A(i+1, j+1) = A(i+1, j+1) + D_N(1,n) .^ (i+j);               
            end
            
        end
    end
    
    for i = 0:M
       A(i+1, i+1) = A(i+1, i+1) + lambda; 
    end
        
    
    w = A\T';
    
end


function w = PolCurFit(D_N, M)
% POLCURFIT Fits a polynomial curve of dimension M to data given by D_N
%
%   Output w: column vector of length M with parameters
    x = D_N(1,:);
    t = D_N(2,:);
    
    for i = 0:M
        T(i+1) = sum(t .* x.^i);
        for j = 0:M
            A(i+1,j+1) = sum(x.^(i+j));            
        end
    end
    
    w = A\T';
end

% function w = PolCurFitLambda(D_N, M, lambda)
% % POLCURFIT Fits a polynomial curve of dimension M to data given by D_N
% %
% %   Output w: column vector of length M with parameters
%     x = D_N(1,:);
%     t = D_N(2,:);
%     
%     for i = 0:M
%         T(i+1) = sum(t .* x.^i);
%         for j = 0:M
%             A(i+1,j+1) = sum(x.^(i+j));            
%         end
%         A(i+1, i+1) = A(i+1, i+1) + lambda;
%     end
%     
%     w = A\T';
% end

function RMS = RMS(data, w)
    RMS = sqrt(2*E(data, w)/length(data));
end

function error = E(data, w)
% E  Calculates the sum of squared error for datapairs (x, t) and
% parameters y
% 
%   data = set of (x; t) input and true values
%   w = parameters to be used in y(x, w)
    error = sum(0.5 * (y(data(1,:), w) - data(2,:)).^2);
end

function y = y(x, w)
% Y  M-th order polynomial of vector x and vector w.
%   Y = Y(x, w) predicted value y using input vector x and parameter vector w

    y = 0;
    for j = 0:length(w)-1 % Loop over the dimensions (M equals the number of weights)
        y = y + w(j+1) * x.^j;
    end
end