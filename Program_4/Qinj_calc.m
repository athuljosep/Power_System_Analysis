function Q_inj_1_Formula = Qinj_calc(x, V, Y)
    Q_inj_1_Formula = -imag(Y(1, 1))*abs(V(1)^2);
    listOfNeighbours = [2,5];
    for k = listOfNeighbours
        if k == 5
            Vk = x(15);
        else
            Vk = V(k);
        end
        Q_inj_1_Formula = Q_inj_1_Formula - abs( Y(1, k) * Vk * V(1) )*sin( angle(Y(1, k)) + x(k-1));
    end
end