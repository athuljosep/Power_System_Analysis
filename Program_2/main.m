clc
clear all; close all;

% Initializing 14 bus and importing data
n_bus = 14;
bus_data = importdata('ieee14bus.txt').data;
branch_data = importdata('ieee14branch.txt').data;
   
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

[V_data,T_data,T1] = NR(bus_data,V,T,P_inj,Q_inj,n_bus,Y,n_pq,pq_i);

V = V_data(:,size(V_data,2));
T = T_data(:,size(T_data,2));

n = 4;
P_inj(n) = 0;

theta = 0;
lambda = 0;
vec = [theta; V(n); lambda];
K = [0;1];

sigma = 0.1;
for k = 1:20
    % Predictor
    J = J_calc(bus_data,V,T,Y,n_bus,n_pq,pq_i);
    a = n-1;
    b = n_bus-1 + find(pq_i == n);
    J_p = [J(a,a) J(a,b); J(b,a) J(b,b)];
    pre(:,k) = vec(:,k) + sigma*inv([J_p K; 0 0 1])*[0;0;1];
    

    % Corrector
    T(n) = pre(1,k);
    V(n) = pre(2,k);
    P_inj(n) = pre(3,k);
    [V_data,T_data,T1] = NR(bus_data,V,T,P_inj,Q_inj,n_bus,Y,n_pq,pq_i);
    V = V_data(:,size(V_data,2));
    T = T_data(:,size(T_data,2));
    vec(:,k+1) = [T(n);V(n);P_inj(n)];
end
pre
vec
plot(vec(3,:),vec(2,:))
% P,Q calculation after convergence
% [P,Q] = PQ_calc(V_data(:,size(V_data,2)),T_data(:,size(T_data,2)),Y)

% plotting convergence curves
% mplot([1:size(V1_data,2)],T1,[1:size(V2_data,2)],T2)