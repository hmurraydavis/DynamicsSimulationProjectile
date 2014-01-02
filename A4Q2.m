    clear all
    close all

    global at1=[]; global at2=[];
    global tens1=[]; global tens2=[];

    g = 9.81; % gravitational acceleration in m/s^2
    m1 = 3.0; m2=4; %kg
    l1=5; l2=3;
    theta1 = 30 * pi/180; %initial phi value = 30 degrees, convert to rad
    theta2 = 45 * pi/180;
    theta2_dot = 1/ (sin(theta1)*l2); %velocity component in theta direction = 0.5 m/s
    
    % Define state variables: 
    %  z1 = theta1, z2 = d(theta1)/dt, z3 = theta2, z4 = d(theta2)/dt

    % Specify initial conditions. 
    z1_0 = theta1;  
    z2_0 = 0; 
    z3_0 = theta2;    %
    z4_0 = theta2_dot;

    Z_0 = [z1_0, z2_0, 0,z3_0, z4_0,0]; %[t1,t1d,t1dd,t2,t2d,t2dd]
    % Define simulation parameters
    t_span = [0:0.01:15];  % max time span for simulation 

    [t, zout] = ode45(@A4Q2_sphpend_fun, t_span, Z_0);


    % x-y position in cartesian coordinates
    t11=zout(:,1);t22=zout(:,3);
    x1 =l1.*sin(t11); x2=x1-(l2.*sin(t22));
    y1=l1.*cos(t11); y2=y1+(l2*cos(t22));
    
    % x-y velocity in cartesian coordinates
    vt11=zout(:,2);vt22=zout(:,4);
    vx1 =l1.*sin(vt11); vx2=vx1-(l2.*sin(vt22));
    vy1=l1.*cos(vt11); vy2=vy1+(l2*cos(vt22));
    
    % x-y acceleration in cartesian coordinates
    ax1 =l1.*sin(at1); ax2=ax1-(l2.*sin(at2));
    ay1=l1.*cos(at1); ay2=ay1+(l2*cos(at2));

    %acceleration at pivots:
    apx1 =sin(at1); apx2=apx1-(sin(at2));
    apy1=cos(at1); apy2=apy1+(cos(at2));

     hold on
    figure;
    %energy graph
         plot(GPE+KE)
         plot(KE)
         plot(GPE)
         xlabel('Time (s)', 'FontSize', 16)
         ylabel('Energy (J)', 'FontSize', 16)
         title('Total Energy of a 2D pendulum, 30 degrees start', 'FontSize', 20)
         legend('Gravitational and Kinetic Energy','Kinetic Energy','Gravitational Energy')
    
    figure;
    %position graph
        plot3(x1,-y1,t);
        plot3(x2,-y2,t);
        legend('Upper Pendulum', 'Lower Pendulum')
        xlabel('x position (m)', 'FontSize', 16)
        ylabel('y position (m)', 'FontSize', 16)
        title('Double Pendulum, 30, 0 degrees start', 'FontSize', 20)

    figure;
    %velocity graph
        plot(t,vx1);
        plot(t,vy1);
        plot(t,vx2);
        plot(t,vy2);
        legend('X velocity--Upper Pendulum','Y velocity--Upper Pendulum', 'X velocity--Lower Pendulum','Y velocity--Lower Pendulum')
        xlabel('Time (s)', 'FontSize', 16)
        ylabel('Velocity (m/s)', 'FontSize', 16)
        title('Velocity Double Pendulum, 30, 0 degrees start', 'FontSize', 20)

    figure;
    %acceleration graph
        plot(linspace(0,length(ax1)),ax1);
        plot(linspace(0,length(ay1)),ay1);
        plot(linspace(0,length(ax2)),ax2);
        plot(linspace(0,length(ay2)),ay2);
        legend('X acceleration--Upper Pendulum','Y acceleration--Upper Pendulum', 'X acceleration--Lower Pendulum','Y acceleration--Lower Pendulum')
        xlabel('Time (s)', 'FontSize', 16)
        ylabel('Acceleration (m/s^2)', 'FontSize', 16)
        title('Acceleration Double Pendulum, 30, 0 degrees start', 'FontSize', 20)

    figure;
    %Rod tension graph:
        plot(linspace(0,length(tens1)), tens1);
        plot(linspace(0,length(tens2)), tens2);
        legend('Tension Rod 1','Tension Rod 2')
        xlabel('Time (s)', 'FontSize', 16)
        ylabel('Tension (N)', 'FontSize', 16)
        title('Rod Tension Double Pendulum, 30, 0 degrees start', 'FontSize', 20)

    figure;
    %reaction force at pivots graph
        plot(linspace(0,length(apx1)),m1*apx1);
        plot(linspace(0,length(apy1)),m1*apy1);
        plot(linspace(0,length(apx2)),m2*apx2);
        plot(linspace(0,length(apy2)),m2*apy2);
        legend('X Reaction--Upper Pivot','Y Reaction--Upper Pivot', 'X Reaction--Lower Pivot','Y Reaction--Lower Pivot')
        xlabel('Time (s)', 'FontSize', 16)
        ylabel('Force (N)', 'FontSize', 16)
        title('Pivotal Reaction Forces Double Pendulum, 30, 0 degrees start', 'FontSize', 20)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function states = A4Q2_sphpend_fun(T, ZZ) %ZZ=[t1,t1d,t1dd,t2,t2d,t2dd]
    
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
