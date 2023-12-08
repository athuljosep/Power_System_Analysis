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
for i = 1:n_br
    a = fr(i);
    b = to(i);
    
    Pij(i,1) = -(V(a)^2)*real(Y(a,b)) + abs(V(a)*V(b)*Y(a,b))*cos(angle(Y(a,b)) + T(b) - T(a));
    Pji(i,1) = -(V(b)^2)*real(Y(b,a)) + abs(V(b)*V(a)*Y(b,a))*cos(angle(Y(b,a)) + T(a) - T(b));
   
    Qij(i,1) = -((V(a)^2)*(abs(B(b))/2 - imag(Y(a,b))) + abs(V(a)*V(b)*Y(a,b))*sin(angle(Y(a,b)) + T(a) - T(b)));
    Qji(i,1) = -((V(b)^2)*(abs(B(a))/2 - imag(Y(b,a))) + abs(V(b)*V(a)*Y(b,a))*sin(angle(Y(b,a)) + T(b) - T(a)));
end

V_tol = 2/100;
P_tol = 5/100;

V_m = V + V*(2*rand - 1)*V_tol;
P_inj_m = P_inj + P_inj*(2*rand - 1)*P_tol;
Q_inj_m = Q_inj + Q_inj*(2*rand - 1)*P_tol;
Pij_m = Pij + Pij*(2*rand - 1)*P_tol;
Pji_m = Pji + Pji*(2*rand - 1)*P_tol;
Qij_m = Qij + Qij*(2*rand - 1)*P_tol;
Qji_m = Qji + Qji*(2*rand - 1)*P_tol;

V_m = [1.05993063580404;1.04493161737285;1.00993390770008;1.01760426332505;1.01944714898401;1.06992998142484;1.06145007406503;1.08992867266642;1.05586262839580;1.05091585629349;1.05683736076918;1.05511951488222;1.05031298056000;1.03546218777418];
P_inj_m = [2.38377333324533;0.187706764192726;-0.966228261582230;-0.490294170951492;-0.0779547217412415;-0.114880642566040;0;0;-0.302587406758766;-0.0923148020619965;-0.0359002008018875;-0.0625689213975754;-0.138472203092995;-0.152832283413750];
Q_inj_m = [-0.173288021693307;0.304535754100071;0.0451164080148254;0.0399895434676861;-0.0164059665508456;0.0481925267431089;0;0.178414886240446;-0.170211902965023;-0.0594716287468153;-0.0184567123697013;-0.0164059665508456;-0.0594716287468153;-0.0512686454713925];
Pij_m = [1.64095482447598,0.789819239173311,0.766046314921584,0.587121061176332,0.434248963723849,-0.243562361858336,-0.639699438089355,0.293648697563946,0.168190154751257,0.461142066005819,0.0769133939511572,0.0814402659935928,0.185639292229293,3.47184809818362e-15,0.293648719579583,0.0546788406192374,0.0985975244078439,-0.0395935226504773,0.0168847224311412,0.0590332442952940]';
Pji_m = [-1.50149807445794,-0.715863533350507,-0.697824152366375,-0.535856598964357,-0.399642274211587,0.232815020125079,0.606882671500350,-0.276260725274226,-0.158231024080903,-0.433836222214063,-0.0718141866195875,-0.0759112707591972,-0.172559934795839,-2.47468775344497e-15,-0.276260745986239,-0.0513144206686193,-0.0916162172179976,0.0373728555951265,-0.0158229413452751,-0.0550055306106929]';
Qij_m = [0.766362692400290,0.388256906401231,0.361438522721456,0.320342027811948,0.253529105736288,-0.139086197766226,-0.242998971930279,-0.212518019760864,-0.0663138882064762,-0.206562308179906,0.0802475953937794,0.0770822338427663,0.187459754975283,-0.173093122219019,0.0582796374420552,0.0667547264690764,0.0964341691121495,-0.0388439288406370,0.0154405317434403,0.0558175818866599]';
Qji_m = [-0.666218500184412,-0.309022238228266,-0.300964220883450,-0.302217280261940,-0.260844797522527,0.104452412457026,0.215464219372961,0.209935561961953,0.0607769695752364,0.224769579656482,-0.0831956415928852,-0.0797954326446385,-0.183546042093953,0.171367095506841,-0.0483917298098923,-0.0640296770150756,-0.0905754565937821,0.0377381430067003,-0.0148317329744392,-0.0527464432027977]';

% Initializing SE parameters
eps = 0.001;
err = 1;
iter = 0;

V_m(3) = 4;
Pij_m(9) = 20;
Qij_m(12) = 20;

while err >= eps
        % Pbus and Qbus
        [P,Q] = PQ_calc(V,T,Y);

        % Pline and Qline
        for i = 1:n_br
            a = fr(i);
            b = to(i);
    
            Pij(i,1) = -(V(a)^2)*real(Y(a,b)) + abs(V(a)*V(b)*Y(a,b))*cos(angle(Y(a,b)) + T(b) - T(a));
            Pji(i,1) = -(V(b)^2)*real(Y(b,a)) + abs(V(b)*V(a)*Y(b,a))*cos(angle(Y(b,a)) + T(a) - T(b));
   
            Qij(i,1) = -((V(a)^2)*(abs(B(b))/2 - imag(Y(a,b))) + abs(V(a)*V(b)*Y(a,b))*sin(angle(Y(a,b)) + T(a) - T(b)));
            Qji(i,1) = -((V(b)^2)*(abs(B(a))/2 - imag(Y(b,a))) + abs(V(b)*V(a)*Y(b,a))*sin(angle(Y(b,a)) + T(b) - T(a)));
        end

        % forming the error vector
        
        Err_vec = [V_m - V ; P_inj_m - P_inj ; Q_inj_m - Q_inj ; Pij_m - Pij ; Pji_m - Pji ; Qij_m - Qij ; Qji_m - Qji];
        % Forming the matrix of partial derivatives
        Hx = Hx_calc(V, T, P, Q, n_bus,bus_data,branch_data,T,Y,n_pq,pq_i);

        % Weight matrix
        R = (diag([ones(1, length(V)) * V_tol, ones(1,length(Hx) - length(V)) * P_tol]))^2;

        % minimizing the error
        x = [T(2:end); V]; 
        x_prev = x;
        x = x + inv(Hx' * inv(R) * Hx) * Hx' * inv(R) * Err_vec;
        T(2:end) = x(1:n_bus-1);
        V = x(n_bus:end);
        iter = iter + 1;
        err = max(abs(x - x_prev));
end

disp('-----------------------');
disp('    Voltages & angles:');
V_T = [V, T];

[P,Q] = PQ_calc(V,T,Y);

for i = 1:n_br
            a = fr(i);
            b = to(i);
    
            Pij(i,1) = -(V(a)^2)*real(Y(a,b)) + abs(V(a)*V(b)*Y(a,b))*cos(angle(Y(a,b)) + T(b) - T(a));
            Pji(i,1) = -(V(b)^2)*real(Y(b,a)) + abs(V(b)*V(a)*Y(b,a))*cos(angle(Y(b,a)) + T(a) - T(b));
   
            Qij(i,1) = -((V(a)^2)*(abs(B(b))/2 - imag(Y(a,b))) + abs(V(a)*V(b)*Y(a,b))*sin(angle(Y(a,b)) + T(a) - T(b)));
            Qji(i,1) = -((V(b)^2)*(abs(B(a))/2 - imag(Y(b,a))) + abs(V(b)*V(a)*Y(b,a))*sin(angle(Y(b,a)) + T(b) - T(a)));
end

Err_vec = [V_m - V ; P_inj_m - P_inj ; Q_inj_m - Q_inj ; Pij_m - Pij ; Pji_m - Pji ; Qij_m - Qij ; Qji_m - Qji];

f = 0;
for i = 1:length(Err_vec)
    f = f + R(i, i) * Err_vec(i).^2;
end

disp('-----------------------');
disp('Error is:');
disp(f)

k = length(Err_vec) - length(x)

Err_thresh = chi2inv(0.99, k)

if f*100 >= Err_thresh
   disp('Bad data detected')
else
   disp('Measurements are good')
end
   
errors = (V - bus_data(:, 5))./bus_data(:, 5) * 100;
% plot(errors)



