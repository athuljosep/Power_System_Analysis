function J = J_full_calc(V, T, P, Q,n_bus,Y)
for m = 1:n_bus 
    for n = 1:n_bus
        if m ~= n
            J11(m, n) = -abs(V(m) * V(n) * Y(m, n)) * sin(angle(Y(m, n)) + T(n) - T(m));
            J21(m, n) = -abs(V(m) * V(n) * Y(m, n)) * cos(angle(Y(m, n)) + T(n) - T(m));
            J12(m, n) = abs(V(m) * V(n) * Y(m, n)) * cos(angle(Y(m, n)) + T(n) - T(m));
            J22(m, n) = -abs(V(m) * V(n) * Y(m, n)) * sin(angle(Y(m, n)) + T(n) - T(m));
        else
            J11(m, m) = -Q(m) - (abs(V(m))^2) * imag(Y(m,m));
            J21(m, m) = P(m) - (abs(V(m))^2) * real(Y(m,m));
            J12(m, m) = P(m) + (abs(V(m))^2) * real(Y(m,m));
            J22(m, m) = Q(m) - (abs(V(m))^2) * imag(Y(m,m));
        end
    end
end
J=[J11 J12; J21 J22];
end