function P_inj_temp = Pinj_calc(x, V, Y)
P_inj_temp = real(Y(1, 1))*abs(V(1)^2);
listOfNeighbours = [2,5];
for k = listOfNeighbours
    if k == 5
        Vk = x(15);
    else
        Vk = V(k);
    end
    P_inj_temp = P_inj_temp + abs( Y(1, k)*Vk* V(1) )*cos( angle(Y(1, k)) + x(k-1));
end
P_inj_temp;
end
