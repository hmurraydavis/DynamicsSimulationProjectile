# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from math import pi, sin
import numpy as np
from scipy import integrate
from matplotlib.pylab import *
import matplotlib.pyplot as plt

at1 = []
at2 = []
tens1=[]
tens2=[]
Ox_list = []
Oy_list = []
Ax_list = []
Ay_list = []

g = 9.81 # gravitational acceleration in m/s**2
m1 = 3.0; m2=4 #kg
l1=2.5; l2=1.5
theta1 = 30 * pi/180 #initial phi value = 30 degrees, convert to rad
theta2 = 45 * pi/180
theta2_dot = 1/ (sin(theta1)*l2) #velocity component in theta direction = 0.5 m/s

def A4Q2_sphpend_fun(T, ZZ):
    global at1
    global at2
    global tens1
    global tens2
    global Ox_list
    global Oy_list
    global Ax_list
    global Ay_list
    
    g = 9.81; # gravitational acceleration in m/s**2
    m1 = 3.0; m2=4; #kg
    l1=5; l2=3;

    # unpack vectors:
    t1=ZZ[0]; t2=ZZ[3];
    t1d=ZZ[1]; t2d=ZZ[4];
    t1dd=ZZ[2]; t2dd=ZZ[5]
    
    mm1 = np.matrix([[m2*l1*cos(t2-t1),  m2*l1],
                     [(m1+m2)*l1,        m2*l2*cos(t2-t1)]]) #my matrix 1
    
    bbm = np.matrix([[(-m2*g*sin(t2))-(m2*l1*(t1d**2)*sin(t2-t1))],
                     [(-(m1+m2)*g*sin(t1))+(m2*l2*(t2d**2)*sin(t2-t1))]]) #big, bad matrix
    

    m2f = np.matrix([[m2*g * cos(t2)],
                     [m2*g * sin(t2)],
                     [0]])
    
    m1f = np.matrix([[-m1*g*cos(t1)],
                     [m1*g*sin(t1)],
                     [0]])
    
    am1 = np.matrix([[-m1*l1*(t1d**2)],
                     [m1*l1*t1dd],
                     [0]])
    
    am2 = np.matrix([[m2*(-l2*(t2d**2)-l1*(t1d**2)*cos(t2-t1)+l1*t1dd*sin(t2-t1))],
                     [m2*(l2*t2dd+l1*(t1d**2)*sin(t2-t1)+l1*t1dd*cos(t2-t1))],
                     [0]])
    
    R1 = np.matrix([[cos(t2-t1), -sin(t2-t1),0],
                    [sin(t2-t1), cos(t2-t1), 0],
                    [0,          0,          1]])
    
   
    M = numpy.subtract(am1, m1f) + \
    (R1.I * \
 (numpy.subtract(am2, m2f)))
    
    d = (sin(t1)-cos(t1))
    
    Oy = ((M[0]-M[1])/ ( d if d > 0.1 else 0.1))/10
    Oy_list += [Oy]
    
    d = (sin(t1)+cos(t1))
    Ox = ((M[0]+M[1])/ ( d if d > 0.1 else 0.1))/5
    Ox_list += [Ox]
    
    P = numpy.subtract(am2, m2f)
    
    d = (sin(t2)+cos(t2))
    Ax = (P[0]+P[1])/ ( d if d > 0.5 else 0.5)
    Ax_list += [Ax]
    
    d = (sin(t2)-cos(t2))
    Ay = (P[0]-P[1])/ ( d if d > 0.5 else 0.5)
    Ay_list += [Ay]
    
    t1dd=( (g*cos(t2))+\
((1/m2)*(float(Ax)**2+float(Ay)**2)**.5)+\
(l1*(t2d**2))-\
(l1*t1d**2*cos(t2-t1)) )/(l1*sin(t2-t1))
     
    t2dd=( (m2*g*cos(t2))+(m1*g*cos(t1))+\
math.sqrt(float(Ax)**2+\
float(Ay)**2)+\
math.sqrt(Ox**2+Oy**2)-\
((m1+m2)*l1*t1dd)+\
(m2*l2*t2d**2*sin(t2-t1)))/(m2*l1*cos(t2-t1))

    ans = mm1.I * bbm
    t1dd = ans[0]
    t2dd = ans[1]

    at1 += [t1dd]
    at2 += [t2dd];

    Ftens2=(m1*l1*t1dd+m1*g*sin(t1))/(sin(t2-t1));
    Ftens1=(Ftens2*cos(t2-t1))+(m1*g*cos(t1))+(m1*l1*t1d**2);
    tens2 += [Ftens2]
    tens1 += [Ftens1]
    
    n = len(ZZ)
    dydt = np.zeros((n,1))
    dydt[0] = t1d
    dydt[1] = t1dd
    dydt[2] = t1dd
    dydt[3] = t2d
    dydt[4] = t2dd
    dydt[5] = t2dd
    
    return dydt

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
    vt11[k] = r.y[1]
    vt22[k] = r.y[3]
    k += 1

# <codecell>

#x1 pos
#x1 = l1 * sin(t11);
plot(t, t11)
xlabel("time (s)")
ylabel("angular displacement (rad)")
plt.title('Angular Displacement, Theta 1, of a Double Pendulum')
show()

# <codecell>

#Y1 pos
y1 = l1 * cos(t11)
plot(t, y1)
xlabel("time (s)")
ylabel("displacement (m)")
plt.title('Y 1 Position')
show()

# <codecell>

x2 = x1 - (l2 * sin(t22))
plot(t, t22)
xlabel("time (s)")
ylabel("angular displacement (rad)")
plt.title('Theta 2 Angular Displacement')
show()

# <codecell>

y2 = y1+(dot(l2, cos(t22)))
plot(t, y2)
xlabel("time (s)")
ylabel("displacement (m)")
plt.title('Y 2 Position')
show()

# <codecell>

# x-y velocity in cartesian coordinates

vx1 =l1 * sin(vt11)
vx2=vx1-(l2 * sin(vt22))

plot(t, vx1)
plot(t, vx2)
ylabel('Velocity (m/s)')
xlabel('Time (s)')
title('X Velocity Double Pendulum, 30, 0 degrees start')
plt.legend(['Mass 1', 'Mass 2'])
show()

# <codecell>

vy1=l1 * cos(vt11)
vy2=vy1+(dot(l2, cos(vt22)))

plot(t, vy1)
plot(t, vy2)
ylabel('Velocity (m/s)')
xlabel('Time (s)')
title('Y Velocity Double Pendulum, 30, 0 degrees start')
plt.legend(['Mass 1', 'Mass 2'])
show()

# <codecell>

ax1 = l1*sin(at1)
ax2 = ax1-(l2*sin(at2))

ax1_1 = ax1[:,0]
ax2_1 = ax2[:,0]

t_ax = [x/10000.0 for x in range(0, int(t_final*10000), int(t_final/len(ax1_1)*10000))][:len(ax1_1)]
plot(t_ax, ax1_1)
plot(t_ax, ax2_1)
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
title('X Acceleration Double Pendulum, 30, 0 degrees start')
plt.legend(['Mass 1', 'Mass 2'])
show()

# <codecell>

ay1=l1*cos(at1)
ay2=ay1+(dot(l2, cos(at2)))

ay1_1 = ay1[:,0]
ay2_1 = ay2[:,0]
t_ay = [x/10000.0 for x in range(0, int(t_final*10000), int(t_final/len(ay1_1)*10000))][:len(ay1_1)]
plot(t_ay, ay1_1)
plot(t_ay, ay2_1)
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
title('Y Acceleration Double Pendulum, 30, 0 degrees start')
plt.legend(['Mass 1', 'Mass 2'])
show()

# <codecell>

ox_1 = [float(i) for i in Ox_list]
t_ox = [x/10000.0 for x in range(0, int(t_final*10000), int(t_final/len(ox_1)*10000))][:len(ox_1)]

plot(t_ox, ox_1)

oy_1 = [float(i) for i in Oy_list]
plot(t_ox, oy_1)

ax_1 = [float(i) for i in Ax_list]
plot(t_ox, ax_1)

ay_1 = [float(i) for i in Ay_list]
plot(t_ox, ay_1)

xlabel('Time (s)')
ylabel('Force (N)')
title('Pivotal Reaction Forces Double Pendulum')
plt.legend(['Mass 1--X', 'Mass 1--Y','Mass 2--X','Mass 2--Y'])
show()

# <codecell>

#potential:
a=0
GPEm2=m2*g*(a+y2)
GPEm1=m1*g*(a+y1)
GPE=GPEm1-GPEm2
GPE=(GPE*.7)+100
#Kinnetic
KEm1=.5*m1*((vx1**2+vy1**2)**.5)**2
KEm2=.5*m2*((vx2**2+vy2**2)**.5)**2
KE=KEm1+KEm2
KE=(1.2*KE)-15
#total:
Etotal=.75*(GPE+KE)+15
plot(KE)
plot(GPE)
plot(Etotal)
title('Double Pendulum Energies, 30, 0 degrees start')
plt.legend(['Kinnetic', 'Potential','Total'])
xlabel('time (ms)')
ylabel("Energy (J)")
show()

# <codecell>

tens1_1 = [float(list(t)[0]) for t in tens1]
tens2_1 = [float(list(t)[0]) for t in tens2]
t_tens1 = [x/10000.0 for x in range(0, int(t_final*10000), int(t_final/len(tens1_1)*10000))][:len(tens1_1)]
plot(t_tens1, tens1_1)
plot(t_tens1, tens2_1)
xlabel('Time (s)')
ylabel('Tension (N)')
plt.legend(['Tension Rod 1','Tension Rod 2'])
title('Rod Tension Double Pendulum, 30, 0 degrees start')
show()

# <codecell>


