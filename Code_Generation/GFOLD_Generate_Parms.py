# GFOLD_parms_static_gen

import numpy as np

def e(i): # create a specific basis vector
    if i==0:
        return [1,0,0]
    if i==1:
        return [0,1,0]
    if i==2:
        return [0,0,1]

def S_mat(_w_): # _w_ to distinguish from our global namespace's w!
    return np.matrix([[0,-_w_[2],+_w_[1]],
                     [_w_[2],0, -_w_[0]],
                     [-_w_[1],_w_[0],0]])

def A(w):
    _A_ = np.empty([6,6])
    np.copyto(_A_[0:3,0:3] , np.zeros((3,3))     ) # top left
    np.copyto(_A_[0:3,3:6] , np.eye(3)           ) # top right
    np.copyto(_A_[3:6,0:3] , -np.square(S(w))    ) # bottom left
    np.copyto(_A_[3:6,3:6] , np.multiply(-1,S(w))) # bottom right
    return _A_

''' Numerical Example 1 '''
s = [ # scalars
    #'N'     : 100,             # Deprecated, replaced by N_tf, static
    ['tf'    , 90],
    ['g0'    , 9.80665],         # standard gravity [m/s**2]
    ['m_dry' , (2)*1e3],         # dry mass kg
    ['m_fuel', (0.3)*1e3],       # fuel in tons
    ['T_max' , 24000],           # thrust max
    ['tmin' ,  0.2],             # throttle ability
    ['tmax' ,  0.8],
    ['G_max' , 3],              # maximum allowable structural Gs
    ['Isp'   , 203.94 ],        # fuel efficiency (specific impulse)
    ['V_max' , 90 ] ,           # velocity max
    ['y_gs'  , np.radians(30)],  # glide slope cone, must be 0 < Degrees < 90
    ['p_cs'  , np.cos(np.radians(45))],  # thrust pointing constraint
]
v = [ # vectors

    ['g' , np.array([-3.71,0,0])],                 # gravity
    ['w' , np.array([2.53*1e-5, 0, 6.62*1e-5])] ,  # planetary angular velocity
    ['nh', np.array([1,0,0])   ],                  # thrust vector reference direction

    ['r0' , np.array([2400, 2000, 0]) ],             # initial position
    ['v0' , np.array([-40,  30,   0]) ],             # initial velocity
    #['v0' , np.array([0,  0,   0]) ],             # initial velocity
    #['r0' , np.array([2400, 0, 0]) ],             # initial position
    #['v0' , np.array([-40,  0,   0]) ],             # initial velocity

    ['rf3', np.array([0,0,0])   ]    ,               # final position target for p4
    ['rf' , np.array([0,0,0])   ]    ,               # final position target
    ['vf' , np.array([0,0,0])   ]                    # final velocity target
]

sk = [k[0] for k in s]
sv = [n[1] for n in s]
# derived values:
s += [
        ['alpha' , 1/(sv[sk.index('Isp')]*sv[sk.index('g0')])    ],     # fuel consumption parameter
        ['m_wet' , (sv[sk.index('m_dry')]+sv[sk.index('m_fuel')])],  # wet mass kg
        ['r1'    , sv[sk.index('tmin')]*sv[sk.index('T_max')] ], # lower thrust bound
        ['r2'    , sv[sk.index('tmax')]*sv[sk.index('T_max')] ],   # upper thrust bound
        #['z0' , np.log(sv[sk.index('m_dry')]+sv[sk.index('m_fuel')])] # initial log(mass) constraint
        ['z0' , np.log(sv[sk.index('m_dry')]+sv[sk.index('m_fuel')])] # initial log(mass) constraint
]
v += [
        ['c' , np.divide(e(0),np.tan(sv[sk.index('y_gs')]))],
]
S,Sk,n=[],{},0
for loople in (s): # 'loople' = a list who wants to be a tuple, but wants assignment too :)
    key = loople[0]
    value=loople[1]
    Sk[key] = n
    S.append( value)
    n+=1
S=np.matrix(S)

V,Vk,n=[],{},0
for loople in (v):
    key = loople[0]
    value=loople[1]
    Vk[key] = n
    V.append( value)
    n+=1

V = np.matrix(V).transpose() # form into shape (width,height) not (height,width)

print('MAKE S HAVE SHAPE',S.shape)
print('MAKE V HAVE SHAPE',V.shape)
