from math import pi, sin
import numpy as np
from scipy import integrate
from matplotlib.pylab import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

at1 = []
at2 = []
tens1=[]
tens2=[]

g = 9.81 # gravitational acceleration in m/s**2
m1 = 3.0; m2=4 #kg
l1=5; l2=3
theta1 = 30 * pi/180 #initial phi value = 30 degrees, convert to rad
theta2 = 45 * pi/180
theta2_dot = 1/ (sin(theta1)*l2) #velocity component in theta direction = 0.5 m/s

def A4Q2_sphpend_fun(T, ZZ):
    global at1
    global at2
    global tens1
    global tens2
    
    g = 9.81; # gravitational acceleration in m/s**2
    m1 = 3.0; m2=4; #kg
    l1=5; l2=3;

    # unpack vectors:
    t1=ZZ[0]; t2=ZZ[3];
    t1d=ZZ[1]; t2d=ZZ[4];
    t1dd=ZZ[2]; t2dd=ZZ[5];
    #t1dd= ( -(m1+m2)*g*sin(t1)-m2*l2*t2dd*cos(t2-t1)+m2*l2*t2d**2*sin(t2-t1) )/(m1+m2)/l1;
    #t2dd=( m2*l1*t1dd*cos(t2-t1)+m2*l1*t1d**2*sin(t2-t1)+m2*g*sin(t2) )/-m2/l2;
    
    mm1 = np.matrix([[m2*l1*cos(t2-t1),  m2*l1],
                     [(m1+m2)*l1,        m2*l2*cos(t2-t1)]]) #my matrix 1
    
    bbm = np.matrix([[(-m2*g*sin(t2))-(m2*l1*(t1d**2)*sin(t2-t1))],
                     [(-(m1+m2)*g*sin(t1))+(m2*l2*(t2d**2)*sin(t2-t1))]]) #big, bad matrix
    ans = mm1.I * bbm
    t1dd = ans[0]
    t2dd = ans[1]

    n1=-1*(m1+m2)*g*sin(t1);
    n2=-m2*l2*t2d**2*sin(t2-t1);
    n3=m2*g*sin(t2)*cos(t2-t1);
    n4=-m2*l1*t1d**2*sin(t2-t1)*cos(t2-t1);
    n5=(m1+m2)*l1-m2*l1*cos(t2-t1)**2;
    
    t1ddC=(n1+n2+n3+n4)/n5;
    
    n6= m2*g*sin(t1)*cos(t2-t1);
    n7= -m2**2*l2*t2d**2*sin(t2-t1)*cos(t2-t1);
    n8= -m2*g*sin(t2);
    n9= -m2*l1*t1d**2*sin(t2-t1);
    n10= (-1*m2**2*cos(t2-t1)**2*l2/(m1+m2))+(m2*l2); #denom
    
    t2ddC= (n6+n7+n8+n9)/n10;
    at1 += [t1ddC]
    at2 += [t2ddC];

    Ftens2=(m1*l1*t1ddC+m1*g*sin(t1))/(sin(t2-t1));
    Ftens1=(Ftens2*cos(t2-t1))+(m1*g*cos(t1))+(m1*l1*t1d**2);
    tens2 += [Ftens2]
    tens1 += [Ftens1]
    
    n = len(ZZ)
    dydt = np.zeros((n,1))
    dydt[0] = t1d
    dydt[1] = t1ddC
    dydt[2] = t1dd
    dydt[3] = t2d
    dydt[4] = t2ddC
    dydt[5] = t2dd
    
    return dydt
    #return mat[[t1d],[t1ddC],[t1dd],[t2d],[t2ddC],[t2dd]]

r = integrate.ode(A4Q2_sphpend_fun).set_integrator('dopri5')
#r = integrate.ode(A4Q2_sphpend_fun).set_integrator('vode', method='bdf')

t_start = 0.0
t_final = 15.0
delta_t = 0.01

num_steps = np.floor((t_final - t_start)/delta_t) + 1

z1_0 = theta1
z2_0 = 0
z3_0 = theta2
z4_0 = theta2_dot

r.set_initial_value([z1_0, z2_0, 0,z3_0, z4_0, 0], t_start)

t = np.zeros((num_steps, 1))
zout = np.zeros((num_steps, 1))
t11 = np.zeros((num_steps, 1))
t22 = np.zeros((num_steps, 1))
vt11  = np.zeros((num_steps, 1))
vt22 = np.zeros((num_steps, 1))

t[0] = t_start
zout[0] = theta1

k = 1
while r.successful() and k < num_steps:
    r.integrate(r.t + delta_t)
    t[k] = r.t
    zout[k] = r.y[0]
    t11[k] = r.y[0]
    t22[k] = r.y[2]
    k += 1
# All done!  Plot the trajectories:
#plot(t, zout)
'''grid('on')
xlabel('Time [minutes]')
ylabel('Concentration [mol/L]')'''

# x-y position in cartesian coordinates
x1 = l1 * sin(t11);
x2 = x1 - (l2 * sin(t22))
y1 = l1 * cos(t11)
y2 = y1+(dot(l2, cos(t22)))

# x-y velocity in cartesian coordinates

vx1 =l1 * sin(vt11)
vx2=vx1-(l2 * sin(vt22))
vy1=l1 * cos(vt11)
vy2=vy1+(dot(l2, cos(vt22)))

# x-y acceleration in cartesian coordinates
ax1 =l1*sin(at1); ax2=ax1-(l2*sin(at2));
ay1=l1*cos(at1); ay2=ay1+(dot(l2, cos(at2)));

#acceleration at pivots:
apx1 =sin(at1); apx2=apx1-(sin(at2));
apy1=cos(at1); apy2=apy1+(cos(at2));

plot(t, x1)
xlabel("time")
ylabel("x1 position")
show()

plot(t, y1)
xlabel("time")
ylabel("y1 position")
show()

plot(t, x2)
xlabel("time")
ylabel("x2 position")
show()

plot(t, y2)
xlabel("time")
ylabel("y2 position")
show()

plot(t,vx1)
plot(t,vx2)
ylabel('Velocity (m/s)')
title('Velocity Double Pendulum, 30, 0 degrees start')
show()

plot(ax1)
plot(ax2)
xlabel('Time (s)')
ylabel('Acceleration (m/s**2)')
show()
