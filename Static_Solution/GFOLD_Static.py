# GFOLD_static_p3p4

from cvxpy import *
from time import time
import numpy as np
from GFOLD_parms_static import *
from EvilPlotting import *

''' As defined in the paper...

 PROBLEM 3: Minimum Landing Error (tf roughly solved)
 MINIMIZE : norm of landing error vector
 SUBJ TO  :
            0) initial conditions satisfied (position, velocity)
            1) final conditions satisfied (altitude, velocity)
            2) dynamics always satisfied
            3) x stays in cone at all times
            4) relaxed convexified mass and thrust constraints
            5) thrust pointing constraint
            6) sub-surface flight constraint

 PROBLEM 4: Minimum Fuel Use
 MAXIMIZE : landing mass, opt variables are dynamical and
 SUBJ TO  :
            0) same constraints as p1, plus:
            1) landing point must be equal or better than that found by p1

'''

def GFOLD(inputs): # PRIMARY GFOLD SOLVER

    #dt = 0.24 #1e0 # dynamics precision ----> BEWARE OF MEMORY OVERFLOW!

    if inputs[-1]=='p3':
        program = 3
        tf_,r0,prog_flag = inputs
    elif inputs[-1]=='p4':
        program = 4
        tf_,r0,rf_,prog_flag=inputs
        #N =int(tf_/dt)

    N  = 250  # Need not be fixed in STATIC runs, MUST be fixed in code-gen
    dt = 4.5  # Integration dt

    print('N = ',N)
    print('dt= ',dt)

    x0=Parameter()
    x0=np.array([r0[0],r0[1],r0[2],v0[0],v0[1],v0[2]])

    x =Variable(6,N) # state vector (3position,3velocity)
    u =Variable(3,N) # u = Tc/mass because Tc[:,n]/m[n] is not allowed by DCP
    z= Variable(1,N)  # z = ln(mass)
    s= Variable(1,N) # thrust slack parameter

    con = []  # CONSTRAINTS LIST

    con += [x[0:3:1,0]  == x0[0:3:1]]
    con += [x[3:6,0]  == x0[3:6]]
    con += [x[3:6,N-1]== vf] # don't forget to slow down, buddy!

    con += [s[N-1] == 0] # thrust at the end must be zero
    con += [u[:,0] == s[0]*np.array([1,0,0]) ] # thrust direction starts straight
    con += [u[:,N-1] == s[N-1]*np.array([1,0,0]) ] # and ends straight
    con += [z[0] == log(m_wet)] # convexified (7)

    if program==3:
        #con += [np.multiply(e(0).T,r[:,N-1]) == rf] # end altitude (fully general)
        con += [x[0,N-1] == 0]

    elif program==4:
        con += [x[0:3,N-1] == rf_] # force landing point equal to found p1 pt
        #con += [norm(E*(x[0:3,N-1]-rf))<=norm(rf_-rf)] # CONVEX <= CONVEX (?)

    for n in range(0,N-1): # any t in [0,tf] maps to any n in [0,N-1]

        # Leapfrog Integration Method
        #    accurate +/- sqrt( (dt*df/dr)**2 + 1)
        #    https://goo.gl/jssWkB
        #    https://en.wikipedia.org/wiki/Leapfrog_integration

        # Dynamics --> v = A(w)*x + B*(g + u)

        con += [x[3:6,n+1] == x[3:6,n] + (dt/2)*((u[:,n]+g) + (u[:,n+1]+g))]
        con += [x[0:3,n+1] == x[0:3,n] + (dt/2)*(x[3:6,n+1]+x[3:6,n])]

        #con += [ norm(E*(x[0:3,n]-rf)) - c.T*(x[0:3,n]-rf) <= 0 ] # glideslope, full generality # (5)
        con += [ norm( (x[0:3,n]-rf)[0:2] ) - c.T[0]*(x[0,n]-rf[0])  <= 0 ] # specific, but faster

        con += [ norm(x[3:6,n]) <= V_max ] # velocity
        con += [z[n+1] == z[n] - (alpha*dt/2)*(s[n] + s[n+1])] # mass decreases
        con += [norm(u[:,n]) <= s[n]] # limit thrust magnitude & also therefore, mass

        # Thrust pointing constraint
        #con += [ nh.T*u[:,n] >= np.cos(p_cs)*G[n]  ] # full generality
        con += [ u[0,n] >= np.cos(p_cs)*s[n]  ]

        if n > 0:
            z0_term = m_wet - alpha * r2 * (n) * dt  # see ref [2], eq 34,35,36
            z1_term = m_wet - alpha * r1 * (n) * dt
            z0 = log( z0_term )
            z1 = log( z1_term )
            mu_1 = r1/(z1_term)
            mu_2 = r2/(z0_term)

            # lower thrust bound
            con += [s[n] <= mu_2 * (1 - (z[:,n] - z0))] # upper thrust bound
            con += [z[n] >= z0] # Ensures physical bounds on z are never violated
            con += [z[n] <= z1]


    con += [x[0,0:N-1] >= 0] # no, this is not the Boring Company!

    if program == 3:
        print('-----------------------------')
        objective=Minimize(norm(x[0:3,N-1]-rf))
        problem=Problem(objective,con)
        obj_opt=problem.solve(solver=ECOS,verbose=True,feastol=5e-20)#solver=SCS,max_iters=5000,verbose=True,use_indirect=False)
        print('-----------------------------')
    elif program == 4:
        print('-----------------------------')
        objective=Maximize(z[N-1])
        problem=Problem(objective,con)
        obj_opt=problem.solve(solver=ECOS,verbose=True)#solver=SCS,max_iters=5000,verbose=True,use_indirect=False,warm_start=True) # OK to warm start b/c p1 gave us a decent answer probably
        print('-----------------------------')

    if program == 3:
        #return obj_opt,(N/dt),x[0:3,N-1]
        if z.value is not None:
            m     = map(np.exp,z.value.tolist()[0]) # make a mass iterable fm z
            return obj_opt,x,u,m,(N/dt),s,z # N/dt is tf
        else:
            return obj_opt,None,None,None,(N/dt),None,None #
    elif program == 4:
        if z.value is not None:
            m     = map(np.exp,z.value.tolist()[0]) # make a mass iterable fm z
            return obj_opt,x,u,m,(N/dt),s,z # N/dt is tf
        else:
            return obj_opt,None,None,None,None,(N/dt),None,None

def P3_P4(r0=r_):

    start_ = time()

    tf_min = (m_dry*np.linalg.norm(v0/r2))
    tf_max = (m_wet-m_dry)/(alpha*r1)
    print 'min tf :%f max tf: %f' % (tf_min,tf_max)

    t = [60,60]

    obj,x,u,m,tf,s,z = GFOLD((t,r0,'p3')) # EXECUTE PROBLEM 3
    print 'p3 object :%f after %f sec' % (obj,time()-start_)
    print 'tf : %f' % (tf)
    print 'rf :'
    if x is None: # Better luck next time! :)
        print('      '+str(None))
        return None
    for r in x[0:3,-1].value:
        print('      '+str(r))
    print

    '''
    obj,x,u,m,tf,s,z = GFOLD((tf,r0,rf,'p4')) # EXECUTE PROBLEM 4
    print 'p4 object :%f after %f sec' % (obj,time()-start_)
    print 'tf : %f' % (tf)
    print 'rf :'
    for r in x[0:3,-1].value:
        print('      '+str(r))
    print

    print 'gfold took: %f sec'%(time()-start_)
    '''

    # Debugging stuff:

    #obj,r,v,u,m = yielded[[vector[0] for vector in yielded].index(min(tf_yield))]
    #tf_opt = tf_array[[vector[0] for vector in yielded].index(min(tf_yield))]

    #for var in (tf,r,v,u,m):
    #    print('var =',var)
    #    print(' ')
    #    print('varval =',var.value)

    x=x.value
    u=u.value
    s=s.value
    z=z.value
    plot_run3D(tf,x,u,m,s,z)
    return obj,x,u,m,tf

if __name__ == '__main__':
    P3_P4(r_)
