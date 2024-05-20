clear all
close all
clc

X1=-50:0.1:50;
X2=-50:0.1:50;
[x1,x2]=meshgrid(X1,X2);

F = x1.^2 + 2*x2.^2 - 0.3*cos(3*pi*x1) - 0.4*cos(4*pi*x2) + 0.7;
realFMin = min(min(F))

max_iteration = 4;

mesh(x1,x2,F)

figure
contourf(x1,x2,F)
hold on

%% Newton-Raphson
fprintf('Newton-Raphson Algorithm\n');

% Initial guesses
% [0, 1] -> [-2, 2]
x01 = -2 + 4 * rand(2, 1);
x02 = -2 + 4 * rand(2, 1);
x03 = -2 + 4 * rand(2, 1);

x = x03;
epsilon=10^(-4);

tic

fprintf('k=1, x1=%f, x2=%f, f(x)=%f\n',x01(1),x01(2),func(x01))
plot(x01(1),x01(2),'r.')

fprintf('k=2, x1=%f, x2=%f, f(x)=%f\n',x02(1),x02(2),func(x02))
plot(x02(1),x02(2),'g.')

fprintf('k=3, x1=%f, x2=%f, f(x)=%f\n',x03(1),x03(2),func(x03))
plot(x03(1),x03(2),'k.')

x_next = x - inv(hessianfunc(x)) * gradfunc(x);

fprintf('k=4, x1=%f, x2=%f, f(x)=%f, error=%f\n',x_next(1),x_next(2), ...
    func(x_next),norm(gradfunc(x_next)))

plot(x_next(1),x_next(2),'m*')

k=6;

while(norm(gradfunc(x_next))>epsilon) % Main step of our algorithm
    x = x_next;
    x_next = x - inv(hessianfunc(x)) * gradfunc(x);

    fprintf('k=%d, x1=%f, x2=%f, f(x)=%f, error=%f\n',k,x_next(1), ...
        x_next(2),func(x_next),norm(gradfunc(x_next)))

    plot(x_next(1),x_next(2),'m*')
    k = k + 1;

    if(k > max_iteration)
        break;
    end
end

toc

title('Newton-Raphson Algorithm')
set(gca,'fontsize',35)
set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);

%% Hestenes-Stiefel Algorithm

figure
contourf(x1,x2,F)
hold on

fprintf('Hestenes-Stiefel Algorithm\n');

% Initial guesses
% [0, 1] -> [-50, 50]
x01 = -50 + 100 * rand(2, 1);
x02 = -50 + 100 * rand(2, 1);
x03 = -50 + 100 * rand(2, 1);

x = x03;
epsilon=10^(-4);

tic

fprintf('k=1, x1=%f, x2=%f, f(x)=%f\n',x01(1),x01(2),func(x01))
plot(x01(1),x01(2),'r.')

fprintf('k=2, x1=%f, x2=%f, f(x)=%f\n',x02(1),x02(2),func(x02))
plot(x02(1),x02(2),'g.')

fprintf('k=3, x1=%f, x2=%f, f(x)=%f\n',x03(1),x03(2),func(x03))
plot(x03(1),x03(2),'k.')

g = gradfunc(x);
d = -g;

% alpha argmin procedure
alpha = 0:0.01:1;
funcalpha = zeros(length(alpha), 1);
for i=1:length(alpha)
    funcalpha(i) = func(x + alpha(i)*d);
end
[val, ind] = min(funcalpha);
alpha = alpha(ind);
% end of alpha argmin procedure

x_next = x + alpha * d;
g_next = gradfunc(x_next);

beta = (g_next' * (g_next - g)) / (d' * (g_next - g)); 

d_next = -g_next + beta * d;

fprintf('k=4, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n',x_next(1),x_next(2), ...
    func(x_next),norm(gradfunc(x_next)))

plot(x_next(1),x_next(2),'m*')
k=5;

while(norm(gradfunc(x_next))>epsilon)
    x = x_next;
    g = g_next;
    d = d_next;

    % alpha argmin procedure
    alpha = 0:0.01:1;
    funcalpha = zeros(length(alpha), 1);
    for i=1:length(alpha)
        funcalpha(i) = func(x + alpha(i)*d);
    end
    [val, ind] = min(funcalpha);
    alpha = alpha(ind);
    % end of alpha argmin procedure

    x_next = x + alpha * d;
    g_next = gradfunc(x_next);
    beta = (g_next' * (g_next - g)) / (d' * (g_next - g));
    d_next = -g_next + beta * d;

    fprintf('k=%d, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n',k,x_next(1), ...
        x_next(2),func(x_next),norm(gradfunc(x_next)))

    plot(x_next(1),x_next(2),'m*')
    k=k+1;
end
toc

title('Hestenes-Stiefel Algorithm')
set(gca,'fontsize',35)
set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);

%% Polak-Ribiere Algorithm

figure
contourf(x1,x2,F)
hold on

fprintf('Polak-Ribiere Algorithm\n');

% Initial guesses
% [0, 1] -> [-50, 50]
x01 = -50 + 100 * rand(2, 1);
x02 = -50 + 100 * rand(2, 1);
x03 = -50 + 100 * rand(2, 1);

x = x03;
epsilon=10^(-4);

tic

fprintf('k=1, x1=%f, x2=%f, f(x)=%f\n',x01(1),x01(2),func(x01))
plot(x01(1),x01(2),'r.')

fprintf('k=2, x1=%f, x2=%f, f(x)=%f\n',x02(1),x02(2),func(x02))
plot(x02(1),x02(2),'g.')

fprintf('k=3, x1=%f, x2=%f, f(x)=%f\n',x03(1),x03(2),func(x03))
plot(x03(1),x03(2),'k.')

g = gradfunc(x);
d = -g;

% alpha argmin procedure
alpha = 0:0.01:1;
funcalpha = zeros(length(alpha), 1);
for i=1:length(alpha)
    funcalpha(i) = func(x + alpha(i)*d);
end
[val, ind] = min(funcalpha);
alpha = alpha(ind);
% end of alpha argmin procedure

x_next = x + alpha * d;
g_next = gradfunc(x_next);

beta = (g_next' * (g_next - g)) / (g' * g); 

d_next = -g_next + beta * d;

fprintf('k=4, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n',x_next(1),x_next(2), ...
    func(x_next),norm(gradfunc(x_next)))

plot(x_next(1),x_next(2),'m*')
k=5;

while(norm(gradfunc(x_next))>epsilon)
    x = x_next;
    g = g_next;
    d = d_next;

    % alpha argmin procedure
    alpha = 0:0.01:1;
    funcalpha = zeros(length(alpha), 1);
    for i=1:length(alpha)
        funcalpha(i) = func(x + alpha(i)*d);
    end
    [val, ind] = min(funcalpha);
    alpha = alpha(ind);
    % end of alpha argmin procedure

    x_next = x + alpha * d;
    g_next = gradfunc(x_next);
    beta = (g_next' * (g_next - g)) / (g' * g);
    d_next = -g_next + beta * d;

    fprintf('k=%d, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n',k,x_next(1), ...
        x_next(2),func(x_next),norm(gradfunc(x_next)))

    plot(x_next(1),x_next(2),'m*')
    k=k+1;
end
toc

title('Polak-Ribiere Algorithm')
set(gca,'fontsize',35)
set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);

%% Fletcher-Reeves Algorithm

figure
contourf(x1,x2,F)
hold on

fprintf('Fletcher-Reeves Algorithm\n');

% Initial guesses
% [0, 1] -> [-50, 50]
x01 = -50 + 100 * rand(2, 1);
x02 = -50 + 100 * rand(2, 1);
x03 = -50 + 100 * rand(2, 1);

x = x03;
epsilon=10^(-4);

tic

fprintf('k=1, x1=%f, x2=%f, f(x)=%f\n',x01(1),x01(2),func(x01))
plot(x01(1),x01(2),'r.')

fprintf('k=2, x1=%f, x2=%f, f(x)=%f\n',x02(1),x02(2),func(x02))
plot(x02(1),x02(2),'g.')

fprintf('k=3, x1=%f, x2=%f, f(x)=%f\n',x03(1),x03(2),func(x03))
plot(x03(1),x03(2),'k.')

g = gradfunc(x);
d = -g;

% alpha argmin procedure
alpha = 0:0.01:1;
funcalpha = zeros(length(alpha), 1);
for i=1:length(alpha)
    funcalpha(i) = func(x + alpha(i)*d);
end
[val, ind] = min(funcalpha);
alpha = alpha(ind);
% end of alpha argmin procedure

x_next = x + alpha * d;
g_next = gradfunc(x_next);

beta = (g_next' * g_next) / (g' * g); 

d_next = -g_next + beta * d;

fprintf('k=4, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n',x_next(1),x_next(2), ...
    func(x_next),norm(gradfunc(x_next)))

plot(x_next(1),x_next(2),'m*')
k=5;

while(norm(gradfunc(x_next))>epsilon)
    x = x_next;
    g = g_next;
    d = d_next;

    % alpha argmin procedure
    alpha = 0:0.01:1;
    funcalpha = zeros(length(alpha), 1);
    for i=1:length(alpha)
        funcalpha(i) = func(x + alpha(i)*d);
    end
    [val, ind] = min(funcalpha);
    alpha = alpha(ind);
    % end of alpha argmin procedure

    x_next = x + alpha * d;
    g_next = gradfunc(x_next);
    beta = (g_next' * g_next) / (g' * g);
    d_next = -g_next + beta * d;

    fprintf('k=%d, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n',k,x_next(1), ...
        x_next(2),func(x_next),norm(gradfunc(x_next)))

    plot(x_next(1),x_next(2),'m*')
    k=k+1;
end
toc

title('Fletcher-Reeves Algorithm')
set(gca,'fontsize',35)
set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);

%% Gradient Descent Algorithm

figure
contourf(x1,x2,F)
hold on

fprintf('Gradient Descent Algorithm\n');

% Initial guesses
% [0, 1] -> [-50, 50]
x01 = -50 + 100 * rand(2, 1);
x02 = -50 + 100 * rand(2, 1);
x03 = -50 + 100 * rand(2, 1);

x = x03;
epsilon=10^(-4);

tic

fprintf('k=1, x1=%f, x2=%f, f(x)=%f\n',x01(1),x01(2),func(x01))
plot(x01(1),x01(2),'r.')

fprintf('k=2, x1=%f, x2=%f, f(x)=%f\n',x02(1),x02(2),func(x02))
plot(x02(1),x02(2),'g.')

fprintf('k=3, x1=%f, x2=%f, f(x)=%f\n',x03(1),x03(2),func(x03))
plot(x03(1),x03(2),'k.')

g = gradfunc(x);

alpha = 0.01;

x_next = x - alpha * g;
g_next = gradfunc(x_next);

fprintf('k=4, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n',x_next(1),x_next(2), ...
    func(x_next),norm(gradfunc(x_next)))

plot(x_next(1),x_next(2),'m*')
k=5;

while(norm(gradfunc(x_next))>epsilon)
    x = x_next;
    g = g_next;

    x_next = x - alpha * g;
    g_next = gradfunc(x_next);

    fprintf('k=%d, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n',k,x_next(1), ...
        x_next(2),func(x_next),norm(gradfunc(x_next)))

    plot(x_next(1),x_next(2),'m*')
    k=k+1;
end

toc

title('Gradient Descent Algorithm')
set(gca,'fontsize',35)
set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);
