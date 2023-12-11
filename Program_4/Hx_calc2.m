function H = Hx_calc2(voltage, delta, nbus, bus_data, branch_data, Y,fb,tb,pq)

nsb = 13;
npq = 9;
npv = 4;
nonslackbus = [2 3 4 5 6 7 8 9 10 11 12 13 14];
nbranch = 20;

theta = angle(Y);
H1=zeros(nbus,nsb);
H3=zeros(nbus,nsb);
H5=zeros(nbranch,nsb);
H7=zeros(nbranch,nsb);
H9=zeros(nbranch,nsb);
H11=zeros(nbranch,nsb);
H2=zeros(nbus,npq);
H4=zeros(nbus,npq);
H6=zeros(nbranch,npq);
H8=zeros(nbranch,npq);
H10=zeros(nbranch,npq);
H12=zeros(nbranch,npq);
for a=1:nbus
        for b=1:nsb
            k=nonslackbus(b);
            if k==a
                H1(a,b)=H1(a,b)-voltage(a)*voltage(k)*abs(Y(a,k))*sin(theta(a,k)+delta(k)-delta(a));
                H3(a,b)=H3(a,b)-voltage(a)*voltage(k)*abs(Y(a,k))*cos(theta(a,k)+delta(k)-delta(a));
                for m=1:nbus
                H1(a,b)=H1(a,b)+voltage(a)*voltage(m)*abs(Y(a,m))*sin(theta(a,m)+delta(m)-delta(a));
                H3(a,b)=H3(a,b)+voltage(a)*voltage(m)*abs(Y(a,m))*cos(theta(a,m)+delta(m)-delta(a));
                end
            else
            H1(a,b)=-voltage(a)*voltage(k)*abs(Y(a,k))*sin(theta(a,k)+delta(k)-delta(a));
            H3(a,b)=-voltage(a)*voltage(k)*abs(Y(a,k))*cos(theta(a,k)+delta(k)-delta(a));
            end
        end
end
for a=1:nbranch
    i=fb(a);
    j=tb(a);
    for b=1:nsb
        k=nonslackbus(b);
        if i==k
            H5(a,b)=voltage(i)*voltage(j)*abs(Y(i,j))*sin(theta(i,j)+delta(j)-delta(i));
            H7(a,b)=-voltage(j)*voltage(i)*abs(Y(j,i))*sin(theta(j,i)+delta(i)-delta(j));
            H9(a,b)=voltage(i)*voltage(j)*abs(Y(i,j))*cos(theta(i,j)+delta(j)-delta(i));
            H11(a,b)=-voltage(j)*voltage(i)*abs(Y(j,i))*cos(theta(j,i)+delta(i)-delta(j));
        end
        if j==k
            H5(a,b)=-voltage(i)*voltage(j)*abs(Y(i,j))*sin(theta(i,j)+delta(j)-delta(i));
            H7(a,b)=voltage(j)*voltage(i)*abs(Y(j,i))*sin(theta(j,i)+delta(i)-delta(j));
            H9(a,b)=-voltage(i)*voltage(j)*abs(Y(i,j))*cos(theta(i,j)+delta(j)-delta(i));
            H11(a,b)=voltage(j)*voltage(i)*abs(Y(j,i))*cos(theta(j,i)+delta(i)-delta(j));
        end
    end
end
    for a=1:nbus
        for b=1:npq
            k=pq(b);
            if k==i 
                H2(a,b)=H2(a,b)+voltage(k)*abs(Y(a,k))*cos(theta(a,k)+delta(k)-delta(a));
                H4(a,b)=H4(a,b)-voltage(k)*abs(Y(a,k))*sin(theta(a,k)+delta(k)-delta(a));
                for m=1:nbus
                    H2(a,b)=H2(a,b)+voltage(m)*abs(Y(a,m))*cos(theta(a,m)+delta(m)-delta(a));
                    H4(a,b)=H4(a,b)-voltage(m)*abs(Y(a,m))*sin(theta(a,m)+delta(m)-delta(a));
                end
            else
            H2(a,b)=voltage(a)*abs(Y(a,k))*cos(theta(a,k)+delta(k)-delta(a));
            H4(a,b)=-voltage(a)*abs(Y(a,k))*sin(theta(a,k)+delta(k)-delta(a));
            end
        end
    end
    for a=1:nbranch
    i=fb(a);
    j=tb(a);
    for b=1:npq
        k=pq(b);
        if i==k
            H6(a,b)=voltage(j)*abs(Y(i,j))*cos(theta(i,j)+delta(j)-delta(i))-2*voltage(i)*abs(Y(i,j))*cos(theta(i,j));
            H8(a,b)=voltage(j)*abs(Y(j,i))*cos(theta(j,i)+delta(i)-delta(j));
            H10(a,b)=2*voltage(i)*abs(Y(i,j))*sin(theta(i,j))-voltage(j)*abs(Y(i,j))*sin(theta(i,j)+delta(j)-delta(i));
            H12(a,b)=-voltage(j)*abs(Y(j,i))*sin(theta(j,i)+delta(i)-delta(j));
        end
        if j==k
            H6(a,b)=voltage(i)*abs(Y(i,j))*cos(theta(i,j)+delta(j)-delta(i));
            H8(a,b)=voltage(i)*abs(Y(j,i))*cos(theta(j,i)+delta(i)-delta(j))-2*voltage(j)*abs(Y(j,i))*cos(theta(j,i));
            H10(a,b)=-voltage(i)*abs(Y(i,j))*sin(theta(i,j)+delta(j)-delta(i));
            H12(a,b)=2*voltage(j)*abs(Y(j,i))*sin(theta(j,i))-voltage(i)*abs(Y(j,i))*sin(theta(j,i)+delta(i)-delta(j));
        end
    end
    end
    H=[H1 H2;H3 H4;H5 H6;H7 H8;H9 H10;H11 H12];

end