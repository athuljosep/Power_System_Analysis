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

% Creating indexes for functions
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

% Executing NR for Power Flow results
[V_data,T_data,T1] = NR(bus_data,V,T,P_inj,Q_inj,n_bus,Y,n_pq,pq_i);
V = V_data(:,size(V_data,2));
T = T_data(:,size(T_data,2));

% P,Q calculation after convergence
[P,Q] = PQ_calc(V,T,Y);

% calculating line power flows
[Pij Pji Qij Qji] = PQline_calc(V,T,Y,n_br,fr,to,B)
% for_graph = [[Pij Pji Qij Qji]]
z = [V;P;Q;Pij;Pji;Qij;Qji];

% Voltage & Power measurement error
V_tol = 2/100;
P_tol = 10/100;

% Introducing noise into measurement
n = [V_tol*randn(length([V]), 1);P_tol*randn(length([P;Q;Pij;Pji;Qij;Qji]), 1)];
var = n.^2;
z_m = z + n;
W = zeros(length(var),length(var));
for i=1:length(var)
    for j=1:length(var)
        if i == j
            W(i,j) = 1/var(i);
        end
    end
end

% Initializing SE parameters
error = 1;
iter = 0;
SE_tol = 0.001;

% Adding Bad Data
z_m(2) = 1.5;
z_m(26) = 10;
z_m(54) = 10;

while error >= SE_tol
        % calculating bus powers
        [P,Q] = PQ_calc(V,T,Y);
        
        % calculating line power flows
        [Pij Pji Qij Qji] = PQline_calc(V,T,Y,n_br,fr,to,B);

        % Generating Hx matric
        Hx = Hx_calc(V, T, Y,B,P, Q, n_bus,n_br,fr,to);

        x = [T(2:end); V]; 
        x_prev = x;
        h = [V;P;Q;Pij;Pji;Qij;Qji];
        x = x + inv(Hx' * W * Hx) * Hx' * W * (z_m-h);
        T(2:end) = x(1:n_bus-1);
        V = x(n_bus:end);
        iter = iter + 1;
        error = max(abs(x - x_prev));
end

% calculating chi squared value
chi_squared_value = 0;
for i=1:length(var)
    chi_squared_value = chi_squared_value + (z_m(i) - h(i))^2 / var(i);
end
chi_squared_value

% Finding chi squared limit
chi_squared_limit = chi2inv(0.99, length(z_m) - length(x))

% Checking bad data present or not
if  chi_squared_value >= chi_squared_limit
   disp('Measurement contains Bad Data')
else
   disp('No Bad Data')
end