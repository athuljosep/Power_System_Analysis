% Test code

clc
clear all; close all;

% Initializing 14 bus and importing data
n_bus = 14;
bus_data = importdata('ieee14bus.txt').data;
branch_data = importdata('ieee14branch.txt').data;
   
% Ybus formation
t = 1; % 0 for without tap, 1 for with tap
Y = y_bus_calc(n_bus,bus_data,branch_data,t);

base_MVA = 100;
P_inj = (bus_data(:,8) - bus_data(:,6)) / base_MVA;
Q_inj = (bus_data(:,9) - bus_data(:,7)) / base_MVA;

sl_i = find(bus_data(:,3) == 3);
pv_i = find(bus_data(:,3) == 2);
pq_i = find(bus_data(:,3) == 0);
n_pv = length(pv_i);
n_pq = length(pq_i);

V = [    1.0600
    1.0450
    1.0100
    0.4491
    0.3824
    1.0700
    0.5587
    1.0900
    0.4623
    0.5238
    0.7793
    0.9803
    0.9218
    0.5583];

T = [         0
   -0.5324
   -1.2975
   -1.1671
   -0.8669
   -3.1734
   -2.2285
   -2.2285
   -2.6728
   -2.8606
   -3.0670
   -3.2105
   -3.1895
   -3.1053];


lambda = 2.8622;
[V_data,T_data,T1,dvrg] = NR(bus_data,V,T,P_inj*lambda,Q_inj*lambda,n_bus,Y,n_pq,pq_i);
% J = J_calc(bus_data,V,T,Y,n_bus,n_pq,pq_i);

dvrg
% T1;
