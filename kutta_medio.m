f = @(x, y) x*(y^(1/2));

x0 = 1;
y0 = 4;
h = 0.2;
n = 1;

approx = kutta_med(f, x0, y0, n, h);

disp('Aproximación final:');
disp(approx);

function approx = kutta_med(f, x0, y0, n, h)
    for i = 1:n
        x_eval = [x0, y0];
        
        k1 = f(x_eval(1), x_eval(2));
        x1_val = [x0 + (1/2)*h, y0 + (1/2)*k1*h];
        k2 = f(x1_val(1), x1_val(2));
        
        x1 = x0 + h;
        y1 = y0 + h*k2;
        
        y0 = y1;
        x0 = x1;
    end
    
    approx = [x0, y0];
end


