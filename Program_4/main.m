clc
clear all; close all;
sympref('FloatingPointOutput',true);

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
V = V_data(:,size(V_data,2));
T = T_data(:,size(T_data,2));

% P,Q calculation after convergence
[P,Q] = PQ_calc(V,T,Y);

% Case 1: Total P = 259 MW
P_tot = 259;
P_eq = [8.0 0.004; 6.4 0.0048]; 
[P_1, P_2, cost, x] = ED(P_eq,P_tot);
disp('Economic Dispatch with total P = 259MW')
disp(['Share of Generator-1 is ' num2str(P_1) 'MW']) 
disp(['Share of Generator-2 is ' num2str(P_2) 'MW'])
disp(['Cost of Generation is ' num2str(x) '/MWhr'])
fprintf('Total Cost is $%f/hr\n', cost)
disp('------------------------------------------------------------')
 
% Case 2: Total P = P1 + P2
P_tot = (P(1) + P(2) + bus_data(2,6)/base_MVA) * base_MVA;
[P_1, P_2, cost, x] = ED(P_eq,P_tot);
disp('Economic Dispatch with total P = P1 + P2')
disp(['Share of Generator-1 is ' num2str(P_1) 'MW']) 
disp(['Share of Generator-2 is ' num2str(P_2) 'MW'])
disp(['Cost of Generation is ' num2str(x) '/MWhr'])
fprintf('Total Cost is $%f/hr\n', cost)
disp('------------------------------------------------------------')

syms T_2
PG0 = [P_1; P_2];
QG0 = [Q(1); Q(2)]*base_MVA;   
bus_no = [1 2];
eq = P_tot == bus_data(2,6) + base_MVA*(real(Y(2,2))*V(2)^2 + abs(Y(2,1)*V(1)*V(2))*cos(angle(Y(2,1)) + T_2));
T_2_soln= solve(eq, T_2);
T_2_Val = T_2_soln(1);

% Case 3: OPF without line limits
constraint = 0;
[P_1,P_2, cost] = OPF(V,T,Y,PG0,QG0,bus_no,T_2_Val,base_MVA,constraint,bus_data,n_bus,n_pq,pq_i,P_inj,Q_inj,P_eq);
disp('Optimal Power Flow without line limits')
disp(['Share of Generator-1 is ' num2str(P_1) 'MW']) 
disp(['Share of Generator-2 is ' num2str(P_2) 'MW'])
fprintf('Total Cost is $%f/hr\n', cost)
disp('------------------------------------------------------------')
    
% Case 4: OPF with line limits
constraint = 1;
[P_1,P_2, cost] = OPF(V,T,Y,PG0,QG0,bus_no,T_2_Val,base_MVA,constraint,bus_data,n_bus,n_pq,pq_i,P_inj,Q_inj,P_eq);
disp('Optimal Power Flow with line limits')
disp(['Share of Generator-1 is ' num2str(P_1) 'MW']) 
disp(['Share of Generator-2 is ' num2str(P_2) 'MW'])
fprintf('Total Cost is $%f/hr\n', cost)
disp('------------------------------------------------------------')  