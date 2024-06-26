clear all
success_rate = 0;
success = 0;

number_of_iteration_list = [];

time_elapsed_list = [];

number_of_execution = 1;

i = 0;

while(i < number_of_execution)
    close all
    clc
    
    X1=-50:0.1:50;
    X2=-50:0.1:50;
    [x1,x2]=meshgrid(X1,X2);
    
    F = x1.^2 + 2*x2.^2 - 0.3*cos(3*pi*x1) - 0.4*cos(4*pi*x2) + 0.7;
    realFMin = min(min(F))
    
    max_iteration = 100;
    
    % mesh(x1,x2,F)
    
    figure
    contourf(x1,x2,F)
    hold on
    
    fprintf('Newton-Raphson Algorithm\n');
    
    % Initial guesses
    % [0, 1] -> [-2, 2]
    x01 = -10 + 20 * rand(2, 1);
    x02 = -5 + 10 * rand(2, 1);
    x03 = -2 + 4 * rand(2, 1);
    
    x0 = [x01, x02, x03];
    colors = {'r', 'g', 'k'};
    
    epsilon=10^(-4);
    
    for j = 1:3
        x = x0(:, j);
        k = 1;

        path = x; % Store the path for visualization
    
        tic
        fprintf('k=1, x1=%f, x2=%f, f(x)=%f\n',x(1),x(2),func(x))
        plot(x(1), x(2), '.', 'Color', colors{j})
    
        x_next = x - inv(hessianfunc(x)) * gradfunc(x);
        
        path = [path, x_next];

        fprintf('k=2, x1=%f, x2=%f, f(x)=%f, error=%f\n',x_next(1),x_next(2), ...
        func(x_next),norm(gradfunc(x_next)))
    
        % plot(x_next(1),x_next(2),'*', 'Color', colors{j})
        
        k=3;
    
        while(norm(gradfunc(x_next))>epsilon) % Main step of our algorithm
            x = x_next;
            x_next = x - inv(hessianfunc(x)) * gradfunc(x);
        
            path = [path, x_next];

            fprintf('k=%d, x1=%f, x2=%f, f(x)=%f, error=%f\n',k,x_next(1), ...
                x_next(2),func(x_next),norm(gradfunc(x_next)))
        
            % plot(x_next(1),x_next(2),'*', 'Color', colors{j})
            k = k + 1;
        
            if(k > max_iteration)
                break;
            end
        end
        time_elapsed = toc;

        plot(path(1, :), path(2, :), 'o-', 'Color', colors{j}, 'MarkerFaceColor', colors{j}, 'LineWidth', 1, 'MarkerSize', 2.5);
        plot(path(1, end), path(2, end), 'x', 'Color', colors{j}, 'MarkerSize', 7, 'LineWidth', 2); % final point


        if(k <= max_iteration)
            success = success + 1;
        end
        
        time_elapsed_list = [time_elapsed_list, time_elapsed];
        number_of_iteration_list = [number_of_iteration_list, k];
    end
    
    % title('Newton-Raphson Algorithm')
    set(gca,'fontsize',24)
    set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);
    
    i = i+1;
end

success_rate = success / (number_of_execution * 3)
avarage_of_number_of_iteration = mean(number_of_iteration_list)
avarage_of_time_elapsed = mean(time_elapsed_list)
