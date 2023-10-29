clc
clear all; close all;

% Initializing 14 bus and importing data
n_bus = 14;
bus_data = importdata('ieee14bus.txt');
bus_data = bus_data.data;
branch_data = importdata('ieee14branch.txt');
branch_data = branch_data.data;
   
% Ybus formation
t = 1; % 0 for without tap, 1 for with tap
Y = y_bus_calc(n_bus,bus_data,branch_data,t);

% Scheduled power calculation
base_MVA = 100;
P_inj = (bus_data(:,8) - bus_data(:,6)) / base_MVA;
Q_inj = (bus_data(:,9) - bus_data(:,7)) / base_MVA;

% Finding bus types
pv_i = find(bus_data(:,3) == 2);
pq_i = find(bus_data(:,3) == 0);
n_pv = length(pv_i);
n_pq = length(pq_i);

% Initializing Voltage magnitude and angles
V = bus_data(:,11);
V(find(V(:)==0)) = 1;
T = zeros(n_bus,1);

[V1_data,T_data,T1] = NR(bus_data,V,T,P_inj,Q_inj,n_bus,Y,n_pq,pq_i);
[V2_data,T_data,T2] = FD(bus_data,V,T,P_inj,Q_inj,n_bus,Y,n_pq,pq_i);

% P,Q calculation after convergence
% V = V_data(:,size(V_data,2));
% T = T_data(:,size(T_data,2));
% P = zeros(n_bus,1);
% Q = zeros(n_bus,1);
% for i = 1:n_bus
%     for j = 1:n_bus
%         P(i) = P(i) + V(i)*V(j)*abs(Y(i,j))*cos(T(i)-T(j)-angle(Y(i,j)));
%         Q(i) = Q(i) + V(i)*V(j)*abs(Y(i,j))*sin(T(i)-T(j)-angle(Y(i,j)));
%     end
% end
% P
% Q

x_label = 'Iteration Count'; % x axis label
y_label = 'Error'; % y axis label
legend_name = {'Newton Raphson','Fast Decoupled'}; % legend names

figure('Renderer', 'painters', 'Position', [10 10 1000 400])
plot([1:size(V1_data,2)],T1,'-xb','LineWidth',1.5)
hold on
plot([1:size(V2_data,2)],T2,'-xr','LineWidth',1.5)
xlabel(x_label,'FontSize',18,'FontName','Times New Roman')
ylabel(y_label,'FontSize',18,'FontName','Times New Roman')
legend (legend_name,'Location','northeast')
set(gca,'fontsize',16,'Fontname','Times New Roman','GridAlpha',0.5)
ax = gca;

ax.XRuler.Axle.LineWidth = 1.5;
ax.YRuler.Axle.LineWidth = 1.5;
grid
grid minor
% legend (legend_name,'Location','southeast')
saveas(gca,'plot.png')
% loglog([1:size(V1_data,2)],T1)
% loglog([1:size(V2_data,2)],T2)