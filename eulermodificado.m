f = @(x, y) x - y;

x0 = 0;
xf = 1;
y0 = 2;

n = 5;  % Aumenta el número de iteraciones para mayor precisión

h = (xf - x0) / n;
approx = euler_modificado(f, x0, y0, h, n);

disp('Aproximación final:');
disp(approx(end, :));

function approx = euler_modificado(f, x0, y0, h, n)
    approx = zeros(n+1, 2);
    approx(1,:) = [x0, y0];
    
    for i = 1:n
        x = approx(i, 1);
        y = approx(i, 2);
        
        y_bar = y + h * f(x, y);  % Calcular y-barra
        
        y_new = y + (h/2) * (f(x, y) + f(x + h, y_bar));  % Actualizar y
        
        approx(i+1, :) = [x + h, y_new];
    end
end



