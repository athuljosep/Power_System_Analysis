function [V_data,T_data,Tol_data] = FD(bus_data,V,T,P_inj,Q_inj,n_bus,Y,n_pq,pq_i)
% Initializing index
B = imag(Y);
B_T =  - B(2:n_bus,2:n_bus);
B_V = - B(pq_i,pq_i);
i = 0;
Tol = 1;
del_T = zeros(n_bus,1);
del_V = zeros(n_bus,1);

% Iteration loop
while(Tol > 1e-5 & i < 100)
    i = i+1;
    V = V+del_V;
    T = T+del_T;
    T_data(:,i) = T;
    V_data(:,i) = V;
    [del_P, del_Q] = dpdq_calc(bus_data,V,T,P_inj,Q_inj,n_bus,Y);
    P_T = del_P'./V(2:n_bus);
    d_T = fwd_bwd(B_T,P_T);
    Q_V = del_Q'./V(pq_i);
    d_V = fwd_bwd(B_V,Q_V);
    del_T = [0 d_T]'; % angle calculation
    for j = 1:n_pq
        del_V(pq_i(j)) = d_V(j); % magnitude calculation
    end
    Tol = max(abs([del_T; del_V]));
    Tol_data(i) = Tol;
end
end