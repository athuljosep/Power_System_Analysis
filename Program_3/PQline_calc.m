function [Pij Pji Qij Qji] = PQline_calc(V,T,Y,n_br,fr,to,B)
for i = 1:n_br
    a = fr(i);
    b = to(i);
    
    Pij(i,1) = -(V(a)^2)*real(Y(a,b)) + abs(V(a)*V(b)*Y(a,b))*cos(angle(Y(a,b)) + T(b) - T(a));
    Pji(i,1) = -(V(b)^2)*real(Y(b,a)) + abs(V(b)*V(a)*Y(b,a))*cos(angle(Y(b,a)) + T(a) - T(b));
   
%     Qij(i,1) = -((V(a)^2)*(abs(B(b))/2 - imag(Y(a,b))) + abs(V(a)*V(b)*Y(a,b))*sin(angle(Y(a,b)) + T(a) - T(b)));
%     Qji(i,1) = -((V(b)^2)*(abs(B(a))/2 - imag(Y(b,a))) + abs(V(b)*V(a)*Y(b,a))*sin(angle(Y(b,a)) + T(b) - T(a)));
    

    Qij(i,1) = - V(a)*V(b)*abs(Y(a,b))*sin(angle(Y(a,b)) + T(b) - T(a)) + V(a).^2*abs(Y(a,b))*sin(angle(Y(a,b)));
    Qji(i,1) = - V(b)*V(a)*abs(Y(b,a))*sin(angle(Y(b,a)) + T(a) - T(b)) + V(b).^2*abs(Y(b,a))*sin(angle(Y(b,a)));
    
  
end
end