f = @(y, x, dy, d2y, d3y) d3y + 3*d2y + 2*dy + y - x^3;

x0 = 0;
y0 = 1;
dy0 = 0;
d2y0 = 0;
d3y0 = 0;
h = 0.1;
n = 10;

approx = runge_kutta_4(f, x0, y0, dy0, d2y0, d3y0, n, h);

disp('Aproximación final:');
disp(approx);

function approx = runge_kutta_4(f, x0, y0, dy0, d2y0, d3y0, n, h)
    for i = 1:n
        x = x0 + (i-1)*h;
        x_eval = [y0, x, dy0, d2y0, d3y0];
        
        k1 = h*f(x_eval(1), x_eval(2), x_eval(3), x_eval(4), x_eval(5));
        k2 = h*f(x_eval(1) + 0.5*k1, x_eval(2) + 0.5*h, x_eval(3) + 0.5*k1, x_eval(4) + 0.5*k1, x_eval(5) + 0.5*k1);
        k3 = h*f(x_eval(1) + 0.5*k2, x_eval(2) + 0.5*h, x_eval(3) + 0.5*k2, x_eval(4) + 0.5*k2, x_eval(5) + 0.5*k2);
        k4 = h*f(x_eval(1) + k3, x_eval(2) + h, x_eval(3) + k3, x_eval(4) + k3, x_eval(5) + k3);

        y = y0 + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
        
        dy = dy0 + (1/6)*(f(x_eval(1), x_eval(2), x_eval(3), x_eval(4), x_eval(5)) + 2*f(x_eval(1) + 0.5*k1, x_eval(2) + 0.5*h, x_eval(3) + 0.5*k1, x_eval(4) + 0.5*k1, x_eval(5) + 0.5*k1) + 2*f(x_eval(1) + 0.5*k2, x_eval(2) + 0.5*h, x_eval(3) + 0.5*k2, x_eval(4) + 0.5*k2, x_eval(5) + 0.5*k2) + f(x_eval(1) + k3, x_eval(2) + h, x_eval(3) + k3, x_eval(4) + k3, x_eval(5) + k3))*h/6;
        
        d2y = d2y0 + (1/6)*(f(x_eval(1), x_eval(2), x_eval(3), x_eval(4), x_eval(5)) + 2*f(x_eval(1) + 0.5*k1, x_eval(2) + 0.5*h, x_eval(3) + 0.5*k1, x_eval(4) + 0.5*k1, x_eval(5) + 0.5*k1) + 2*f(x_eval(1) + 0.5*k2, x_eval(2) + 0.5*h, x_eval(3) + 0.5*k2, x_eval(4) + 0.5*k2, x_eval(5) + 0.5*k2) + f(x_eval(1) + k3, x_eval(2) + h, x_eval(3) + k3, x_eval(4) + k3, x_eval(5) + k3))*h/6;
        
        d3y = d3y0 + (1/6)*(f(x_eval(1), x_eval(2), x_eval(3), x_eval(4), x_eval(5)) + 2*f(x_eval(1) + 0.5*k1, x_eval(2) + 0.5*h, x_eval(3) + 0.5*k1, x_eval(4) + 0.5*k1, x_eval(5) + 0.5*k1) + 2*f(x_eval(1) + 0.5*k2, x_eval(2) + 0.5*h, x_eval(3) + 0.5*k2, x_eval(4) + 0.5*k2, x_eval(5) + 0.5*k2) + f(x_eval(1) + k3, x_eval(2) + h, x_eval(3) + k3, x_eval(4) + k3, x_eval(5) + k3))*h/6;

        y0 = y;
        dy0 = dy;
        d2y0 = d2y;
        d3y0 = d3y;
    end
    
    approx = [y0, dy0, d2y0, d3y0];
end



