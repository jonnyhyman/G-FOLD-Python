# GFOLD_static_p3p4_gen_precalcz

import numpy as np
import GFOLD_parms_static_gen as parameters
from EvilPlotting import *

'''

    This code can do both static runs (tests) AND initialize code gen.
    If doing code generation, you must still compile the generated code.

 PROBLEM 1: Minimum Landing Error (tf roughly solved)
 MINIMIZE : norm of landing error vector
 SUBJ TO  :
            0) initial conditions satisfied (position, velocity)
            1) final conditions satisfied (altitude, velocity)
            2) dynamics always satisfied
            3) x stays in cone at all times
            4) relaxed convexified mass and thrust constraints
            5) thrust pointing constraint
            6) sub-surface flight constraint

 PROBLEM 2: Minimum Fuel Use
 MAXIMIZE : landing mass, opt variables are dynamical and
 SUBJ TO  :
            0) same constraints as p1, plus:
            1) landing point must be equal or better than that found by p1
'''

VERSION = 1.0

test = 1  # are we doing a static run or a generation run?

if test:
    from cvxpy import *
else:
    from cvxpy_codegen import *

def GFOLD_C_GEN(prog_flag,_s_,_v_): # PRIMARY GFOLD SOLVER

    N_tf=250  # MUST BE FIXED FOR CODE GEN TO WORK

    sk,vk=parameters.Sk,parameters.Vk
    if not test:
        dt=Parameter(1,1,name='dt') # determines tf implicitly dt = tf/N, tf = dt*N(const)
        S=Parameter(1,17,name='S') # contains all parms_static scalar variables
        V=Parameter(3,9,name='V') # contains all parms_static vect variables
        z0=Parameter(N_tf,name='z0')
        z1=Parameter(N_tf,name='z1')
        mu_1=Parameter(N_tf,name='mu_1')
        mu_2=Parameter(N_tf,name='mu_2')
        z0_term=Parameter(N_tf,name='z0_term')
        z1_term=Parameter(N_tf,name='z1_term')

    else:
        V=_v_ # for cvxpy testing
        S=_s_ # for cvxpy testing

        dt=Parameter(1,1,name='dt') # determines tf implicitly dt = tf/N,
                                    # tf = dt*N(const)

        dt.value = float(S[0,sk['tf']])/(N_tf)

        # Precalculate Z limits, then pass in as a PARAMETER

        z0=Parameter(N_tf,name='z0')
        z1=Parameter(N_tf,name='z1')
        mu_1=Parameter(N_tf,name='mu_1')
        mu_2=Parameter(N_tf,name='mu_2')
        z0_term=Parameter(N_tf,name='z0_term')
        z1_term=Parameter(N_tf,name='z1_term')

        z0_term_, z1_term_ = np.zeros(N_tf),np.zeros(N_tf)
        z0_, z1_           = np.zeros(N_tf),np.zeros(N_tf)
        mu_1_, mu_2_       = np.zeros(N_tf),np.zeros(N_tf)

        for n in range(0,N_tf-1):
            z0_term_[n] = S[0,sk['m_wet']] - S[0,sk['alpha']] * S[0,sk['r2']] * (n) * dt.value  # see ref [2], eq 34,35,36
            z1_term_[n] = S[0,sk['m_wet']] - S[0,sk['alpha']] * S[0,sk['r1']] * (n) * dt.value
            z0_[n] = np.log( z0_term_[n] )
            z1_[n] = np.log( z1_term_[n] )
            mu_1_[n] = S[0,sk['r1']]/(z1_term_[n])
            mu_2_[n] = S[0,sk['r2']]/(z0_term_[n])

        z0_term.value = z0_term_
        z1_term.value = z1_term_
        z0.value = z0_
        z1.value = z1_
        mu_1.value = mu_1_
        mu_2.value = mu_2_

    # new variables here for brevity in the dynamics equations
    c=vk['c']
    g=vk['g']
    rf=vk['rf']
    alpha=sk['alpha']
    #print(c,g,rf,alpha)

    if prog_flag=='p3':
        program = 3
    elif prog_flag=='p4':
        program = 4

    x = Variable(6,N_tf,name='x') # state vector (3position,3velocity)
    u = Variable(3,N_tf,name='u') # u = Tc/mass because Tc[:,n]/m[n] is not allowed by DCP
    z = Variable(1,N_tf,name='z')  # z = ln(mass)
    s = Variable(1,N_tf,name='s') # thrust slack parameter

    con = []

    con += [x[0:3,0]  ==   V[:,vk['r0']]]
    con += [x[3:6,0]  ==   V[:,vk['v0']]]
    con += [x[3:6,N_tf-1]==V[:,vk['vf']]] # don't forget to slow down, buddy!

    con += [s[N_tf-1] == 0] # thrust at the end must be zero
    con += [u[:,0] == s[0]*np.array([1,0,0]) ] # thrust direction starts straight
    con += [u[:,N_tf-1] == s[N_tf-1]*np.array([1,0,0]) ] # and ends straight
    con += [z[0] == S[0,sk['z0']]] # convexified (7)

    if program==3:
        con += [x[0,N_tf-1] == 0] # end altitude

    elif program==4:

        # force landing point equal to found program 3 point
        con += [x[0:3,N_tf-1] == V[:,vk['rf3']]]

    for n in range(0,N_tf-1): # any t in [0,tf] maps to any n in [0,N-1]

        # Leapfrog Integration Method
        #    accurate +/- sqrt( (dt*df/dr)**2 + 1)
        #    https://goo.gl/jssWkB
        #    https://en.wikipedia.org/wiki/Leapfrog_integration

        con += [x[3:6,n+1] == x[3:6,n] + (dt*0.5)*((u[:,n]+V[:,g]) + (u[:,n+1]+V[:,g]))]
        con += [x[0:3,n+1] == x[0:3,n] + (dt*0.5)*(x[3:6,n+1]+x[3:6,n])]

        con += [ norm((x[0:3,n]-V[:,rf])[0:2] ) - V[0,c]*(x[0,n]-V[0,rf])  <= 0 ] # glideslope constraint
        con += [ norm(x[3:6,n]) <= S[0,sk['V_max']] ] # velocity

        con += [z[n+1] == z[n] - (S[0,alpha]*dt*0.5)*(s[n] + s[n+1])] # mass decreases
        con += [norm(u[:,n]) <= s[n]] # limit thrust magnitude & also therefore, mass

        # Thrust pointing constraint
        con += [  u[0,n] >= S[0,sk['p_cs']]*s[n]  ]
        if n > 0:
            # lower thrust bound
            #con += [s[n] >= mu_1 * (1 - (z[:,n] - z0) + (1/2)*square(z[:,n] - z0))]
            con += [s[n] <= mu_2[n]* (1 - (z[:,n] - z0[n]))] # upper thrust bound
            con += [z[n] >= z0[n]] # Ensures physical bounds on z are never violated
            con += [z[n] <= z1[n]]

    con += [x[0,0:N_tf-1] >= 0] # no, this is not the Boring Company!

    if program == 3:
        print('-----------------------------')
        if test:

            objective=Minimize(norm(x[0:3,N_tf-1]-V[:,rf]))
            problem=Problem(objective,con)
            obj_opt=problem.solve(solver=ECOS,verbose=True)
            print(obj_opt)

        else:

            objective=Minimize(norm(x[0:3,N_tf-1]-V[:,rf]))
            problem=Problem(objective,con)
            obj_opt=problem.codegen('GFOLD_'+prog_flag+'_')

        print('-----------------------------')


    elif program == 4:
        print('-----------------------------')
        if test:

            objective=Minimize(z[N_tf-1])
            problem=Problem(objective,con)
            obj_opt=problem.solve(solver=ECOS,verbose=True)
            print(obj_opt)

        else:

            objective=Maximize(z[N_tf-1])
            problem=Problem(objective,con)
            obj_opt=problem.codegen('GFOLD_'+prog_flag)
            print('-----------------------------')

if __name__ == '__main__':

    PROGRAM_TO_COMPILE = 'p3'
    GFOLD_C_GEN(PROGRAM_TO_COMPILE, parameters.S, parameters.V)
