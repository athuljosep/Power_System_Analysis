clc
clear all; close all;

% Initializing 14 bus and importing data
n_bus = 14;
bus_data = importdata('ieee14bus.txt').data;
branch_data = importdata('ieee14branch.txt').data;

base_MVA = 100;
P_inj = (bus_data(:,8) - bus_data(:,6)) / base_MVA;
Q_inj = (bus_data(:,9) - bus_data(:,7)) / base_MVA;
   
% Ybus formation
t = 1; % 0 for without tap, 1 for with tap
Y = y_bus_calc(n_bus,bus_data,branch_data,t);

sl_i = find(bus_data(:,3) == 3);
pv_i = find(bus_data(:,3) == 2);
pq_i = find(bus_data(:,3) == 0);
n_pv = length(pv_i);
n_pq = length(pq_i);

V=[1.0600
1.0450
1.0100
0.7906
0.7770
1.0700
0.8674
1.0900
0.7906
0.8014
0.9169
0.9845
0.9432
0.7545
 
 
];

T=[0
-0.6265
-1.3882
-1.1291
-0.9697
-1.5896
-1.4433
-1.4433
-1.6105
-1.6351
-1.6191
-1.6575
-1.6598
-1.7482
];


V = [1.0600
    1.0450
    1.0100
    0.7590
    0.7437
    1.0700
    0.8429
    1.0900
    0.7604
    0.7751
    0.9036
    0.9816
    0.9372
    0.7283]

T = [0
   -0.5742
   -1.3883
   -1.1304
   -0.9643
   -1.6524
   -1.4792
   -1.4792
   -1.6638
   -1.6928
   -1.6793
   -1.7208
   -1.7226
   -1.8165]

lambda=4.0082;

[del_P, del_Q] = dpdq_calc(bus_data,V,T,P_inj*lambda,Q_inj*lambda,n_bus,Y)