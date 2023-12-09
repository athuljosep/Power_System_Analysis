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

V = V_data(:,size(V_data,2));
T = T_data(:,size(T_data,2));



% calculating line power flows
[Pij Pji Qij Qji] = PQline_calc(V,T,Y,n_br,fr,to,B);
z_actual = [P;Q;Pij;Pji;Qij;Qji];

V_tol = 2/100;
P_tol = 5/100;


noise = P_tol*randn(length(z_actual), 1)
variance=noise.^2
z_meas = z_actual+noise
W=zeros(length(variance),length(variance))
for i=1:length(variance)
    for j=1:length(variance)
        if i==j
            W(i,j)=1/variance(i);
        end
    end
end

% Initializing SE parameters
error = 1;
iter = 0;

% z_meas(2) = 10


while error >= 1e-3
        % Pbus and Qbus
        [P,Q] = PQ_calc(V,T,Y);
        
        % calculating line power flows
        [Pij Pji Qij Qji] = PQline_calc(V,T,Y,n_br,fr,to,B);

        % forming the error vector
%         Err_vec = [V_m - V ; P_inj_m - P_inj ; Q_inj_m - Q_inj ; Pij_m - Pij ; Pji_m - Pji ; Qij_m - Qij ; Qji_m - Qji];
%          Err_vec = [P_inj_m - P_inj ; Q_inj_m - Q_inj ; Pij_m - Pij ; Pji_m - Pji ; Qij_m - Qij ; Qji_m - Qji];
        
        % Forming the matrix of partial derivatives
        
        Hx = Hx_calc2(V, T, n_bus,bus_data,branch_data,Y,fr,to,pq_i);

%         Hx = Hx_calc(V, T, P, Q, n_bus,bus_data,branch_data,T,Y,n_pq,pq_i);

        % Weight matrix
%         W = (diag([ones(1, length(V)) * V_tol, ones(1,length(Hx) - length(V)) * P_tol]))^2;

        % minimizing the error
        x = [T(2:end); V(pq_i)]; 
        x_prev = x;
        h = [P;Q;Pij;Pji;Qij;Qji];
        x = x + inv(Hx' * W * Hx) * Hx' * W * (z_meas-h);
        T(2:end) = x(1:n_bus-1);
        V(pq_i) = x(n_bus:end);
        iter = iter + 1;
        error = max(abs(x - x_prev));
        
end

disp('-----------------------');
disp('    Voltages & angles:');
V_T = [V, T];

[P,Q] = PQ_calc(V,T,Y);

% calculating line power flows
[Pij Pji Qij Qji] = PQline_calc(V,T,Y,n_br,fr,to,B);



chi_squared_value = 0;
for i=1:length(variance)
    chi_squared_value = chi_squared_value + (z_meas(i) - h(i))^2 / variance(i);
end
chi_squared_value



k = length(z_meas) - length(x);
chi_squared_limit = chi2inv(0.99, k)

if  chi_squared_value >= chi_squared_limit
   disp('Measurement contains Bad Data')
else
   disp('No Bad Data')
end
%    
% errors = (V - bus_data(:, 5))./bus_data(:, 5) * 100;
% plot(errors)



