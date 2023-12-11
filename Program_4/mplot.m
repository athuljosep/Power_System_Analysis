function mplot(x1,y1,x2,y2)
x_label = 'Iteration Count'; % x axis label
y_label = 'Error'; % y axis label
legend_name = {'Newton Raphson','Fast Decoupled'}; % legend names

figure('Renderer', 'painters', 'Position', [10 10 1000 400])
plot(x1,y1,'-xb','LineWidth',1.5)
hold on
plot(x2,y2,'-xr','LineWidth',1.5)
xlabel(x_label,'FontSize',18,'FontName','Times New Roman')
ylabel(y_label,'FontSize',18,'FontName','Times New Roman')
legend (legend_name,'Location','northeast')
set(gca,'fontsize',16,'Fontname','Times New Roman','GridAlpha',0.5)
ax = gca

ax.XRuler.Axle.LineWidth = 1.5;
ax.YRuler.Axle.LineWidth = 1.5;
grid
grid minor
% legend (legend_name,'Location','southeast')
saveas(gca,'plot.png')
end