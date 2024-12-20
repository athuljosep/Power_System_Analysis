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

% n=4
% CPF(n,n_bus,bus_data,P_inj,Q_inj,Y,sl_i,pq_i,n_pq);

% Continuation Power Flow for each load bus
for i = 1:n_pq 
    n = pq_i(i);
    CPF(n,n_bus,bus_data,P_inj,Q_inj,Y,sl_i,pq_i,n_pq);
end