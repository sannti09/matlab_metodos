

% Ejemplo de uso:
F = @(t, w) w - t^2 + 1;
a = 0;
b = 2;
w0 = 0.5;
h = 0.2;
n = 4;

[t, w] = metodoTaylor(F, a, b, w0, h, n);






function [t, w] = metodoTaylor(F, a, b, w0, h, n)
    N = (b - a) / h;

    t = zeros(1, N+1);
    w = zeros(1, N+1);

    t(1) = a;
    w(1) = w0;

    for i = 1:N
        t(i+1) = a + i*h;
        dw = zeros(1, n);
        dw(1) = F(t(i), w(i));

        for j = 2:n
            dw(j) = 0;
            for k = 1:(j-1)
                dw(j) = dw(j) + ((t(i)^k) / factorial(k)) * dw(j-k);
            end
            dw(j) = dw(j) / factorial(j-1);
        end
        w(i+1) = w(i) + h * sum(dw);
    end

    disp([t' w']);
end

