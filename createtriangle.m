function createtriangle(x_a, y_a, dy_dx, dx, color)
    x_b = x_a * 10^dx;
    y_b = y_a * 10^-(dy_dx*dx);
    plot(gca, [x_a;x_b], [y_a;y_b],'LineWidth',3, 'Color', color);
    plot(gca, [x_a;x_b], [y_b;y_b],'LineWidth',1, 'Color', color);
    plot(gca, [x_a;x_a], [y_a;y_b],'LineWidth',1, 'Color', color); 
    text(x_a*0.95,0.5*y_a,num2str(dy_dx));
    text(x_a*1.2,0.5*y_b,num2str(1));
end


