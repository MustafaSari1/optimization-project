function ypr = gradfunc(x)
    syms x1 x2;
    F = x1.^2 + 2*x2.^2 - 0.3*cos(3*pi*x1) - 0.4*cos(4*pi*x2) + 0.7;
    g = gradient(F);
    ypr = double(subs(g, [x1 x2], [x(1), x(2)]));
end