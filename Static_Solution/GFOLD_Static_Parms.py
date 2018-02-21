# GFOLD_parms_static

import numpy as np
from scipy import signal

def e(i):
    return signal.unit_impulse(3,i) # create a specific basis vector

def S(_w_): # _w_ to distinguish from our global namespace's w!
    return np.matrix([[0,-_w_[2],+_w_[1]],
                     [_w_[2],0, -_w_[0]],
                     [-_w_[1],_w_[0],0]])

'''------------------------ Numerical Example 1 -------------------------- '''
# These are the numbers from the original paper

g0    = 9.80665         # standard gravity [m/s**2]
m_dry = (2)*1e3         # dry mass kg
m_fuel= (0.3)*1e3       # fuel in tons
T_max = 24000           # thrust max
throt = [0.2,0.8]       # throttle ability
G_max = 3               # maximum allowable structural Gs
Isp   = 203.94          # fuel efficiency (specific impulse)
V_max = 90              # velocity max
y_gs  = np.radians(30)# glide slope cone, must be 0 < Degrees < 90
p_cs  = np.radians(45)  # thrust pointing constraint
alpha = 1/(Isp*g0)      # fuel consumption parameter
m_wet = (m_dry+m_fuel)  # wet mass kg
r1    = throt[0]*T_max  # lower thrust bound
r2    = throt[1]*T_max  # upper thrust bound

g = np.array([-3.71,0,0])                   # gravity
w = np.array([2.53*1e-5, 0, 6.62*1e-5])   # planetary angular velocity
nh= np.array([1,0,0])                     # thrust vector reference direction

r_ = np.array([2400, 2000, 0])              # initial position
v0 = np.array([-40,  30,   0])                 # initial velocity

rf = np.array([0,0,0])                      # final position target
vf = np.array([0,0,0])                      # final velocity target

c = np.divide(e(0),np.tan(y_gs))
E = np.array( [ [e(0).T],[e(1).T] ] )

A = np.empty([6,6])
np.copyto(A[0:3,0:3] , np.zeros((3,3))     ) # top left
np.copyto(A[0:3,3:6] , np.eye(3)           ) # top right
np.copyto(A[3:6,0:3] , -np.square(S(w))    ) # bottom left
np.copyto(A[3:6,3:6] , np.multiply(-1,S(w))) # bottom right
B = np.array([[0,0,0],[0,0,0],[0,0,0],[1,0,0],[0,1,0],[0,0,1]]) # 0vect and I


''' ------------------------ Masten Lander ---------------------------- '''
'''
g0    = 9.80665         # standard gravity [m/s**2]
m_dry = (3.725)*1e3         # dry mass kg
m_fuel= (6.2)*1e3       # fuel in tons
T_max = 13000          # thrust max
throt = [0.1,0.8]       # throttle ability
G_max = 3               # maximum allowable structural Gs
Isp   = 295            # fuel efficiency (specific impulse)
V_max = 330            # velocity max
y_gs  = np.radians(0.1)# glide slope cone, must be 0 < Degrees < 90
p_cs  = np.radians(45)  # thrust pointing constraint
alpha = 1/(Isp*g0)      # fuel consumption parameter
m_wet = (m_dry+m_fuel)  # wet mass kg
r1    = throt[0]*T_max  # lower thrust bound
r2    = throt[1]*T_max  # upper thrust bound

g = np.array([-9.81,0,0])                   # gravity
w = np.array([2.53*1e-5, 0, 6.62*1e-5])   # planetary angular velocity
nh= np.array([1,0,0])                     # thrust vector reference direction

r_ = np.array([160, 0, 0])              # initial position
v0 = np.array([-1, 0, 0 ])                 # initial velocity

rf = np.array([0,0,0])                      # final position target
vf = np.array([0,0,0])                      # final velocity target

c = np.divide(e(0),np.tan(y_gs))
E = np.array( [ [e(0).T],[e(1).T] ] )

A = np.empty([6,6])
np.copyto(A[0:3,0:3] , np.zeros((3,3))     ) # top left
np.copyto(A[0:3,3:6] , np.eye(3)           ) # top right
np.copyto(A[3:6,0:3] , -np.square(S(w))    ) # bottom left
np.copyto(A[3:6,3:6] , np.multiply(-1,S(w))) # bottom right
B = np.array([[0,0,0],[0,0,0],[0,0,0],[1,0,0],[0,1,0],[0,0,1]]) # 0vect and I
'''

''' ------------------------------ Falcon 9 ---------------------------------'''
'''
N     = 50              # nodes in discretization (more --> more accurate)
g0= 9.80665             # standard gravity [m/s**2]
m_dry = (22.2)*1e3      # dry mass kg
m_fuel= (13.4)*1e3      # fuel in tons
T_max = 845000*3        # thrust max
throt = [0.4,1.0]       # throttle ability
G_max = 3               # maximum allowable structural Gs
Isp   = 282             # fuel efficiency (specific impulse)
V_max = 1300            # velocity max
y_gs  = np.radians(1)   # glide slope cone, must be 0 < Degrees < 90
p_cs  = np.radians(120) # thrust pointing constraint
alpha = 1/(Isp*g0)      # fuel consumption parameter
m_wet = (m_dry+m_fuel)  # wet mass kg
r1    = throt[0]*T_max  # lower thrust bound
r2    = throt[1]*T_max  # upper thrust bound

g = np.array([-g0,0,0])                   # gravity
w = np.array([2.91*1e-5, 0, 6.68*1e-5])   # planetary angular velocity
nh= np.array([1,0,0])                     # thrust vector reference direction

r_ = np.array([20000,   10, 5])              # initial position
v0 = np.array([-500,-100,-200])                 # initial velocity

rf = np.array([0,0,0])                      # final position target
vf = np.array([0,0,0])                      # final velocity target

c = np.divide(e(0),np.tan(y_gs))
E = np.array( [ [e(0).T],[e(1).T] ] )
'''
'''
A = np.empty([6,6])
np.copyto(A[0:3,0:3] , np.zeros((3,3))     ) # top left
np.copyto(A[0:3,3:6] , np.eye(3)           ) # top right
np.copyto(A[3:6,0:3] , -np.square(S(w))    ) # bottom left
np.copyto(A[3:6,3:6] , np.multiply(-1,S(w))) # bottom right
B = np.array([[0,0,0],[0,0,0],[0,0,0],[1,0,0],[0,1,0],[0,0,1]]) # 0vect and I
'''
