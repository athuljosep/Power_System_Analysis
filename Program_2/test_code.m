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

sl_i = find(bus_data(:,3) == 3);
pv_i = find(bus_data(:,3) == 2);
pq_i = find(bus_data(:,3) == 0);
n_pv = length(pv_i);
n_pq = length(pq_i);

V = [    1.0600
    1.0450
    1.0100
    0.4470
    0.3780
    1.0700
    0.6012
    1.0900
    0.5128
    0.5675
    0.7999
    0.9787
    0.9212
    0.5834];

T = [         0
   -0.6247
   -1.5170
   -1.3735
   -1.0696
   -3.0504
   -2.3104
   -2.3104
   -2.6929
   -2.8348
   -2.9776
   -3.0957
   -3.0802
   -3.0507
    ];


P_inj = [    6.3711
    0.5017
   -2.5824
   -1.3104
   -0.2083
   -0.3070
         0
         0
   -0.8087
   -0.2467
   -0.0959
   -0.1672
   -0.3701
   -0.4085];

Q_inj = [   -0.4633
    0.8142
    0.1206
    0.1069
   -0.0439
    0.1288
         0
    0.4770
   -0.4551
   -0.1590
   -0.0493
   -0.0439
   -0.1590
   -0.1371
    ];

[V_data,T_data,T1,dvrg] = NR(bus_data,V,T,P_inj,Q_inj,n_bus,Y,n_pq,pq_i);


dvrg
T1
