function [P1,P2, cost] = OPF(V,T,Y,PG0,QG0,bus_no,T_2_Val,base_MVA,limit,bus_data,n_bus,n_pq,pq_i,P_inj,Q_inj,P_eq)
syms PG_i QG_i P12
u = sym('u', [1 1]);
u_r = num2cell(u);
x0 = sym('x', [n_bus+n_pq 1]);
x = x0(2:end);

r = 0.06;
p_coef = 1;
corr = 100;
tol = 0.01;
P12_x = 100.0;
P12_Limit = 5;
iter = 1;

i = bus_no(1);
j = bus_no(2); 
u_x = PG0(j);
xVals = [T_2_Val; T(3:end); V(pq_i)];
P1 = PG0(i);
Q1_Val = QG0(i)/base_MVA;
J = J_calc(bus_data,V,T,Y,n_bus,n_pq,pq_i);
V = V(2:end);
C_i = P_eq(1, 2)*PG_i^2 + P_eq(1, 1)*PG_i; 
C_j(u) = P_eq(2, 2)*u^2 + P_eq(2, 1)*u; 
f(u_r{:}) = C_i + C_j(u_r{:});

if limit == 1     
    f(u_r{:}) = f(u_r{:}) + p_coef*(P12 - P12_Limit)^2;
end
display(f(u_r{:}));

while corr > tol || P12_x*limit > P12_Limit
    if iter > 1
        bus_data(2,8) = u_x;
        P_inj = (bus_data(:,8) - bus_data(:,6)) / base_MVA;
        V = bus_data(:,11);
        V(find(V(:)==0)) = 1;
        T = zeros(n_bus,1);
        
        [V_data,T_data,T1] = NR(bus_data,V,T,P_inj,Q_inj,n_bus,Y,n_pq,pq_i);
        V = V_data(:,size(V_data,2));
        T = T_data(:,size(T_data,2));
        [P,Q] = PQ_calc(V,T,Y);

        xVals = [T(2:n_bus); V(pq_i)];
        P1 = P(1)*base_MVA;
        Q1_Val = Q(1)*base_MVA;
        
        J = J_calc(bus_data,V,T,Y,n_bus,n_pq,pq_i);
        P12_x = abs(abs(V(2)*V(1)*Y(1, 2))*cos(-T(2) - angle(Y(1, 2))) - real(Y(1, 2))*(V(1))^2);
        P12_x = P12_x*base_MVA;
        dfdxP12 = p_coef*[2*(P12_x-5)*abs(V(1)*V(2)*Y(1,2))*sin(T(1)-T(2)-angle(Y(1,2))); zeros(n_bus+n_pq-2, 1)];
    end
    dgdu = [1; zeros(n_bus + n_pq - 2, 1)];
    dgdx = - J;
            
    df_du_temp = jacobian(f(u_r{:}), u);
    dfdu = subs(df_du_temp, u, u_x);
    
    P_inj_1_temp = Pinj_calc(x, V, Y);
    P1_temp = P_inj_1_temp + bus_data(2,6)/base_MVA;
    
    Q_inj_1_temp = Qinj_calc(x, V, Y);
    Q1_temp = Q_inj_1_temp + bus_data(1,7)/base_MVA;
    
    df_dP1_temp = transpose(jacobian(f(u_r{:}), PG_i));
    dfdP1 = subs(df_dP1_temp, PG_i, P1);
    dP1_dx_temp = jacobian(P1_temp, x);
    dP1dx = subs(dP1_dx_temp, x, xVals);
    
    df_dQ1_temp = transpose(jacobian(f(u_r{:}), QG_i));
    dfdQ1 = subs(df_dQ1_temp, QG_i, Q1_Val);
    dQ1_dx_temp = jacobian(Q1_temp, x(14:end));
    dQ1dx = subs(dQ1_dx_temp, x, xVals);

    dfdx = dfdP1*dP1dx;
    
    if limit == 1
        if iter > 1
            dfdx = dfdx + transpose(dfdxP12);
        end
    end
    
    del_C = dfdu - transpose(dgdu)*inv(transpose(dgdx))*transpose(dfdx);
    corr = r*del_C;
    u_x = u_x - corr;
    cost = subs(f(u_x), PG_i, P1);

    iter = iter + 1
    if iter > 50
        break;
    end         
end

P2 = double(u_x);
P12_x = abs(abs(V(2)*V(1)*Y(1, 2))*cos(-T(2) - angle(Y(1, 2))) - real(Y(1, 2))*(V(1))^2);
P12_x = P12_x*base_MVA
if limit == 1
    cost = subs(f(P2), [PG_i P12], [P1 P12_x]);
else
    cost = subs(f(P2), PG_i, P1);
end
end