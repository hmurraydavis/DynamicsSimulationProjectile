function states = A4Q2sphpend_fun(T, ZZ) %ZZ=[t1,t1d,t1dd,t2,t2d,t2dd]
    
    g = 9.81; % gravitational acceleration in m/s^2
    m1 = 3.0; m2=4; %kg
    l1=5; l2=3;

    % unpack vectors:
    t1=ZZ(1); t2=ZZ(4);
    t1d=ZZ(2); t2d=ZZ(5);
    t1dd=ZZ(3); t2dd=ZZ(6);
    t1dd= ( -(m1+m2)*g*sin(t1)-m2*l2*t2dd*cos(t2-t1)+m2*l2*t2d^2*sin(t2-t1) )/(m1+m2)/l1;
    t2dd=( m2*l1*t1dd*cos(t2-t1)+m2*l1*t1d^2*sin(t2-t1)+m2*g*sin(t2) )/-m2/l2;
    
    n1=-1*(m1+m2)*g*sin(t1);
    n2=-m2*l2*t2d^2*sin(t2-t1);
    n3=m2*g*sin(t2)*cos(t2-t1);
    n4=-m2*l1*t1d^2*sin(t2-t1)*cos(t2-t1);
    n5=(m1+m2)*l1-m2*l1*cos(t2-t1)^2;
    
    t1ddC=(n1+n2+n3+n4)/n5;
    
    n6= m2*g*sin(t1)*cos(t2-t1);
    n7= -m2^2*l2*t2d^2*sin(t2-t1)*cos(t2-t1);
    n8= -m2*g*sin(t2);
    n9= -m2*l1*t1d^2*sin(t2-t1);
    n10= (-1*m2^2*cos(t2-t1)^2*l2/(m1+m2))+(m2*l2); %denom
    
    t2ddC= (n6+n7+n8+n9)/n10;
    at1=[at1,t1ddC]; at2=[at2.t2ddC];

    Ftens2=(m1*l1*t1ddC+m1*g*sin(t1))/(sin(t2-t1));
    Ftens1=(Ftens2*cos(t2-t1))+(m1*g*cos(t1))+(m1*l1*t1d^2);
    tens2=[tens2,Ftens2]; tens1=[tens1,Ftens1]
    
   
    states = [t1d;t1ddC;t1dd;t2d;t2ddC;t2dd];
    %
end
