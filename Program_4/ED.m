function [P_1, P_2, cost, x] = ED(P_eq,P_total)
y = sym('y', [2 1]);
syms P_tot x

y_x = num2cell(y);
C_1(y) = P_eq(1, 2)*y(1)^2 + P_eq(1, 1)*y(1); 
C_2(y) = P_eq(2, 2)*y(2)^2 + P_eq(2, 1)*y(2); 

yfull = [y; x];
f(y) = C_1(y_x{:}) + C_2(y_x{:});
c(y) = sum([y_x{:}].*[1 1]) - P_tot; 
c(y) = subs(c(y_x{:}), P_tot, P_total); 

y_x_all = num2cell(yfull);
C(yfull) = C_1(y_x{:}) + C_2(y_x{:}) - x*c(y_x{:}); 
Delta_C(yfull) = jacobian(subs(C(y_x_all{:}), P_tot, P_total), yfull).'; 

y_x_soln = solve(Delta_C(y_x_all{:}), y_x_all{:});
P_1 = double(y_x_soln.y1);
P_2 = double(y_x_soln.y2);
x = double(y_x_soln.x);
cost = f(P_1, P_2);       
end