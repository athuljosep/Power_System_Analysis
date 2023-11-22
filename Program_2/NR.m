function [V_data,T_data,Tol_data,dvrg] = NR(bus_data,V,T,P_inj,Q_inj,n_bus,Y,n_pq,pq_i)
% Initializing index
i = 0;
Tol = 1;
del_T = zeros(n_bus,1);
del_V = zeros(n_bus,1);
dvrg = 0;
iter_limit = 50;

% Iteration loop
while(Tol > 1e-3 & i < 50)
    i = i+1;

    V = V+del_V;
    T = T+del_T;
    T_data(:,i) = T;
    V_data(:,i) = V;
    [del_P, del_Q] = dpdq_calc(bus_data,V,T,P_inj,Q_inj,n_bus,Y);
    if(i==2)

        [e, f] = dpdq_calc(bus_data,V,T,P_inj,Q_inj,n_bus,Y);
    end
    dpdq = [del_P, del_Q]; % mismatch calculation
    c= dpdq';
    J = J_calc(bus_data,V,T,Y,n_bus,n_pq,pq_i); % Jacobian calculation
    if(i<5)
        J;
    end
    delta = fwd_bwd(J,dpdq); % finding errors
    
    a=delta';
    b = inv(J)*dpdq';
    del_T = [0 delta(1:n_bus-1)]';
    for j = 1:n_pq
        del_V(pq_i(j)) = delta(n_bus+j-1);
    end
    Tol = max(abs(delta)); % updating error for convergence
    Tol_data(i) = Tol;
end
i;
if i >= iter_limit
    dvrg = 1;
end
end
