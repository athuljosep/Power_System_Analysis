function CPF(n,n_bus,bus_data,P_inj,Q_inj,Y,sl_i,pq_i,n_pq)
% Initializing Voltage magnitude and angles
V = bus_data(:,11);
V(find(V(:)==0)) = 1;
T = zeros(n_bus,1);



P_k = P_inj;
P_k(sl_i) = [];
K = [P_k;Q_inj(pq_i)];
ek1 = [zeros(1,n_bus-1+n_pq) 1];
ek2 = zeros(1,n_bus+n_pq);
ek2_i = find(pq_i == n);
ek2(n_bus-1+ek2_i) = -1;
ek3 = -ek1;
sigma1 = 0.1;
sigma2 = 0.005;
lambda = 0;
i = 1;

% step 1
[V_data,T_data,T1,dvrg] = NR(bus_data,V,T,P_inj*lambda,Q_inj*lambda,n_bus,Y,n_pq,pq_i);
V = V_data(:,size(V_data,2));
T = T_data(:,size(T_data,2));
y(i) = V(n);
x(i) = lambda;
i = i+1;
dvrg = 0;
while(dvrg == 0 & i < 100)
    
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
        ph1 = i-1;
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
f_stop = 0.75;
if (n == 12 | n ==13)
    f_stop = 0.98
    sigma2 = 0.0005
end
while(lambda > f_stop*prev_lambda & i < 200)
    
    % Predictor
    theta = T(2:n_bus);
    voltage = V(pq_i);
    vec = [theta; voltage; lambda];
    J = J_calc(bus_data,V,T,Y,n_bus,n_pq,pq_i);
    pre = vec + (sigma2*inv([J -K; ek2])*ek1');
    T = [0; pre(1:n_bus-1)];
    for j = 1:n_pq
        V(pq_i(j)) = pre(n_bus+j-1);
    end
    lambda = pre(end);

    % Corrector
    J = J_calc(bus_data,V,T,Y,n_bus,n_pq,pq_i);
    [del_P, del_Q] = dpdq_calc(bus_data,V,T,P_inj*lambda,Q_inj*lambda,n_bus,Y);
    corr = inv([J -lambda*K; ek2])*[del_P del_Q 0]';
    lambda = lambda + corr(end);
    if lambda <= f_stop*prev_lambda
        V = prev_V;
        T = prev_T;
        ph2 = i-1;
        break;
    else
        T = T + [0; corr(1:n_bus-1)];
        for j = 1:n_pq
            V(pq_i(j)) = V(pq_i(j)) + corr(n_bus+j-1);
        end
        prev_V = V;
        prev_T = T;
        y(i) = V(n);
        x(i) = lambda;
        i = i+1;
    end
end

% step 3
dvrg = 0;
while(lambda >= 0 & i < 250)
   
    % Predictor
    i;
    theta = T(2:n_bus);
    voltage = V(pq_i);
    vec = [theta; voltage; lambda];
    J = J_calc(bus_data,V,T,Y,n_bus,n_pq,pq_i);
    pre = vec + sigma1*inv([J -K; ek3])*ek1';
    T = [0; pre(1:n_bus-1)];
    for j = 1:n_pq
        V(pq_i(j)) = pre(n_bus+j-1);
    end
    lambda = pre(end);
    
    % Corrector
    [V_data,T_data,T1,dvrg] = NR(bus_data,V,T,P_inj*lambda,Q_inj*lambda,n_bus,Y,n_pq,pq_i);
    V = V_data(:,size(V_data,2));
    T = T_data(:,size(T_data,2));
    y(i) = V(n);
    x(i) = lambda;
    i = i+1;
end

% plotting code
x_label = '\lambda'; % x axis label
y_label = 'Voltage (p.u.)'; % y axis label
legend_name = {'CPF Phase 1','CPF Phase 2','CPF Phase 3'}; % legend names
title_name = ['PV Curve for bus number ' num2str(n)];
figure('Renderer', 'painters', 'Position', [10 10 800 600])
plot(x(1:ph1),y(1:ph1),'-ob','LineWidth',1.5)
hold on
plot(x(ph1:ph2),y(ph1:ph2),'-ok','LineWidth',1.5)
plot(x(ph2:end),y(ph2:end),'-or','LineWidth',1.5)
xlabel(x_label,'FontSize',18,'FontName','Times New Roman')
ylabel(y_label,'FontSize',18,'FontName','Times New Roman')
legend (legend_name,'Location','northeast')
set(gca,'fontsize',16,'Fontname','Times New Roman','GridAlpha',0.5)
ax = gca
xlim([0 4.5])
title(title_name)
ax.XRuler.Axle.LineWidth = 1.5;
ax.YRuler.Axle.LineWidth = 1.5;
grid
grid minor
saveas(gca,[title_name '.png'])

end