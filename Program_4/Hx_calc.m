function H = Hx_calc(V, T, Y,B,P, Q, n_bus,n_br,fr,to)

J = J_full_calc(V, T, P, Q, n_bus, Y);

H1 = zeros(n_bus, n_bus-1);
H2 = eye(n_bus, n_bus);
H3 = J(1:n_bus, 2:n_bus);
H4 = J(1:n_bus, n_bus+1:end);
H5 = J(n_bus+1:end, 2:n_bus);
H6 = J(n_bus+1:end, n_bus+1:end);

for i = 1:n_br
    a = fr(i);
    b = to(i);
    if a~=1
        H7(i, a-1) = abs(V(a) * V(b) * Y(a,b)) * sin(angle(Y(a,b)) + T(b) - T(a));
        H9(i, a-1) = -abs(V(b) * V(a) * Y(b,a)) * sin(angle(Y(b,a)) + T(a) - T(b));  
        H11(i, a-1) = abs(V(a) * V(b) * Y(a,b)) * cos(angle(Y(a,b)) + T(b) - T(a));
        H13(i, a-1) = -abs(V(b) * V(a) * Y(b,a)) * cos(angle(Y(b,a)) + T(a) - T(b));
    end
    if b~=1
        H7(i, b-1) = -abs(V(a) * V(b) * Y(a,b)) * sin(angle(Y(a,b)) + T(b) - T(a));
        H9(i, b-1) = abs(V(b) * V(a) * Y(b,a)) * sin(angle(Y(b,a)) + T(a) - T(b));
        H11(i, b-1) = -abs(V(a) * V(b) * Y(a,b)) * cos(angle(Y(a,b)) + T(b) - T(a));  
        H13(i, b-1) = abs(V(b) * V(a) * Y(b,a)) * cos(angle(Y(b,a)) + T(a) - T(b));
    end
    
    H8(i, b) = abs(V(a) * Y(a,b)) * cos(angle(Y(a,b)) + T(b) - T(a));
    H8(i, a) = -2*V(a) * real(Y(a,b)) + abs(V(b) * Y(a,b)) * cos(angle(Y(a,b)) + T(b) - T(a));
    
    H10(i,a) = abs(V(b) * Y(b,a)) * cos(angle(Y(b,a)) + T(a) - T(b));
    H10(i,b) = -2*V(b) * real(Y(b,a)) + abs(V(a) * Y(b,a))*cos(angle(Y(b,a)) + T(a) - T(b));
    
    H12(i,b) = -abs(V(a) * Y(a,b)) * sin(angle(Y(a,b)) + T(b) - T(a));
    H12(i,a) = -2 * V(a) * (0.5*abs(B(i)) - imag(Y(a,b))) - abs(V(b) * Y(a,b)) * sin(angle(Y(a,b)) + T(b) - T(a));
    
    H14(i,a) = -abs(V(b) * Y(b,a)) * sin(angle(Y(b,a)) + T(a) - T(b));
    H14(i,b) = -2 * V(b) * (0.5 * abs(B(i)) - imag(Y(b,a))) - abs(V(a) * Y(b,a)) * sin(angle(Y(b,a)) + T(a) - T(b));

end

H=[ H1, H2; H3, H4; H5, H6; H7, H8; H9, H10; H11, H12; H13, H14];

end