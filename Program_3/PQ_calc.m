function [P,Q] = PQ_calc(V,T,Y)
n_bus = size(V,1);
P = zeros(n_bus,1);
Q = zeros(n_bus,1);
for i = 1:n_bus
    for j = 1:n_bus
        P(i) = P(i) + V(i)*V(j)*abs(Y(i,j))*cos(T(i)-T(j)-angle(Y(i,j)));
        Q(i) = Q(i) + V(i)*V(j)*abs(Y(i,j))*sin(T(i)-T(j)-angle(Y(i,j)));
    end
end
end