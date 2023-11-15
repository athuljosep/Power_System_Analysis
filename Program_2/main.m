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
sl_i = find(bus_data(:,3) == 3);
pv_i = find(bus_data(:,3) == 2);
pq_i = find(bus_data(:,3) == 0);
n_pv = length(pv_i);
n_pq = length(pq_i);

% Initializing Voltage magnitude and angles
V = bus_data(:,11);
V(find(V(:)==0)) = 1;
T = zeros(n_bus,1);

% [V_data,T_data,T1] = NR(bus_data,V,T,P_inj,Q_inj,n_bus,Y,n_pq,pq_i);

% V = V_data(:,size(V_data,2));
% T = T_data(:,size(T_data,2));

n = 4;

P_k = P_inj;
P_k(sl_i) = [];
K = [P_k;Q_inj(pq_i)];
ek1 = [zeros(1,n_bus-1+n_pq) 1];
ek2 = zeros(1,n_bus+n_pq);
ek2_i = find(pq_i == n);
ek2(n_bus-1+ek2_i) = -1;
sigma1 = 0.1;
sigma2 = 0.005;
lambda = 0;

i = 1;
dvrg = 0;
while(dvrg == 0 & i < 100)
    % step 1
    % Predictor
    theta = T(2:n_bus);
    voltage = V(pq_i);
    vec = [theta; voltage; lambda];
    J = J_calc(bus_data,V,T,Y,n_bus,n_pq,pq_i);
    pre = vec + sigma1*inv([J -K; ek1])*ek1';
    T = [0; pre(1:n_bus-1)];
    for j = 1:n_pq
        V(pq_i(j)) = pre(n_bus+j-1);
    end
    lambda = pre(end);
    % Corrector
    [V_data,T_data,T1,dvrg] = NR(bus_data,V,T,P_inj*lambda,Q_inj*lambda,n_bus,Y,n_pq,pq_i);
    if dvrg == 1
        V = prev_V;
        T = prev_T;
        lambda = prev_lambda;
        i = i-1;
    else
        V = V_data(:,size(V_data,2));
        T = T_data(:,size(T_data,2));
        prev_V = V;
        prev_T = T;
        prev_lambda = lambda;
        y(i) = V(n);
        x(i) = lambda;
        i = i+1;
    end
end

% step 2
for k = i:100
    k
    % Predictor
    theta = T(2:n_bus);
    voltage = V(pq_i);
    vec = [theta; voltage; lambda]
    J = J_calc(bus_data,V,T,Y,n_bus,n_pq,pq_i);
    pre = vec + sigma2*inv([J K; ek2])*ek1';
    T = [0; pre(1:n_bus-1)];
    for j = 1:n_pq
        V(pq_i(j)) = pre(n_bus+j-1);
    end
    lambda = pre(end);
    T;
    V;
    lambda;
    % Corrector
    J = J_calc(bus_data,V,T,Y,n_bus,n_pq,pq_i);
    [del_P, del_Q] = dpdq_calc(bus_data,V,T,P_inj*lambda,Q_inj*lambda,n_bus,Y);
    [del_P del_Q 0]';
    corr = inv([J lambda*K; ek2])*[del_P del_Q 0]'
    if lambda <= 0.75*prev_lambda
        V = prev_V;
        T = prev_T;
        lambda = prev_lambda;
    else
        T = [0; corr(1:n_bus-1)];
        for j = 1:n_pq
            V(pq_i(j)) = corr(n_bus+j-1);
        end
        lambda = lambda+corr(end);
        T;
        V;
        lambda;
        prev_V = V;
        prev_T = T;
        prev_lambda = lambda;
        y(k) = V(n);
        x(k) = lambda;
    end
end
plot(x,y)
