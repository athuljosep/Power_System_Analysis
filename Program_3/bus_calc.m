clc
clear all; close all;

% Initializing 14 bus and importing data
% n_bus = 14;
% bus_data = importdata('ieee14bus.txt').data;
% branch_data = importdata('ieee14branch.txt').data;

n_bus = 3;
bus_data = importdata('3bus.txt').data;
branch_data = importdata('3branch.txt').data;

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
fr = branch_data(:,1);
to = branch_data(:,2);
n_br = length(fr);

B = branch_data(:,9);

% Initializing Voltage magnitude and angles
V = bus_data(:,11);
V(find(V(:)==0)) = 1;
T = zeros(n_bus,1);

[V_data,T_data,T1] = NR(bus_data,V,T,P_inj,Q_inj,n_bus,Y,n_pq,pq_i);

% P,Q calculation after convergence
[P,Q] = PQ_calc(V_data(:,size(V_data,2)),T_data(:,size(T_data,2)),Y);


for j = 2:size(V_data,2)
    temp1 = PQ_calc(V_data(:,j-1),T_data(:,j-1),Y)
    temp2 = PQ_calc(V_data(:,j),T_data(:,j),Y)
    temp2-temp1
end

% calculating line power flows
[Pij Pji Qij Qji] = PQline_calc(V_data(:,size(V_data,2)),T_data(:,size(T_data,2)),Y,n_br,fr,to,B);

[V_data,T_data,T1] = FD(bus_data,V,T,P_inj,Q_inj,n_bus,Y,n_pq,pq_i);
