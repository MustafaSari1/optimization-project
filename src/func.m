function y = func(x)
    x1 = x(1);
    x2 = x(2);
    y = x1.^2 + 2*x2.^2 - 0.3*cos(3*pi*x1) - 0.4*cos(4*pi*x2) + 0.7;
end