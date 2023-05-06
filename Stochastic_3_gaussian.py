import matplotlib.pyplot as plt
import numpy as np
import scipy.constants
import time
import scipy.integrate as integrate
import scipy as sp
from scipy.integrate import odeint
from matplotlib.pyplot import cm
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.transforms import Bbox
from mpl_toolkits.mplot3d import axes3d
from scipy.optimize import curve_fit

def my_max_by_weight(sequence):
    if not sequence:
        raise ValueError('empty sequence')

    maximum = sequence[0]

    for item in sequence:
        # Compare elements by their weight stored
        # in their second element.
        if item[1] > maximum[1]:
            maximum = item
    return maximum

N = 10 # no of eggs per brooder

def F(N,r_o,alpha):
    #f=  (1/(1+ np.exp(-r_o) ) )
    #print(r_o)
    f = np.power(r_o,alpha)
    #print (f)
    return(f)

def S(M,b,gamma,r_i):  ### brooding period
    s = np.exp(-(M*b)/(1+(gamma*r_i)))
    #print('r_i is', r_i)
    #print ('s is', s)
    return (s)

def S_future(M,gamma,r_i):
    s_fut = np.exp(-(M)/(1+(gamma*r_i)))
    return (s_fut)

def W_prime(r_prime):
    #print (i_prime)
    w=Fofr[ r_range.index(r_prime) ]
    return(w)

def g_func(time):
    g = (0.3/(2+ 4*np.exp(-1/2*time) ) ) #1 logistic
    g = g*100
    g = np.floor(g)
    g = int(g)
    #print ('g is',g)
    return g


pop = 200


#end time
TT= 21 #use odd numbers for proper plotting
T = 15
#start time
#t=1

#the i_critical/min and max values
r_min =0

C=20 # no. of categories

r_max = 1

#########

M = 3.5# 1 or 3.5  Factors causing mortality
gamma = 4 #2 or 4 amount of decreasing/increases mortality. Controld mortality or its rate


#g = 10#int(C/4) ### gain in reserves in the non-brooding period until the next brooding season


alpha = 0.5

phi_max = 100 # maximum number of offspring produced (on average for a species)


##sigma and g### #prob_gaussian##


g = np.array([6,7,8,9,10,11,12,13,14])
#g = np.array([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17])
l = len(g)
g_mean= np.mean(g)
#g_func = g/20 ### normalisation to keep gE[0,1.0]
#g_func_mean = np.mean(g_func)

sigma = 8

prob= np.arange(l,dtype=float)

for p in range (0,l):
    prob[p] = np.exp( (-1 * np.power(g[p]-g_mean,2)  ) / ( 2* sigma * sigma ) )


prob_sum = np.sum(prob) ###
prob = prob/prob_sum ### ensure all probs all to 1
#print(prob)


brooding_period = np.array([0.1,0.5,0.9])  ### brooding period array [0.0,1.0]; ## brooding period; months or days. Proportion of brooding period at a time t. 1-b would be non-brooding period in a time period t.
lim =len(brooding_period)

opt_immunity = np.zeros((lim,pop,T-1))
opt_offspring = np.zeros((lim,pop,T-1))

sur = np.zeros((lim,T-1))### no. of alive individuals at a certain time for different brooding periods


for B in range (0,lim):
    b = brooding_period[B]
    print ('b is', b)

    np.random.seed(200)

    ### initializing arrays for forward dynamics ####
    opt_ro = np.zeros((pop,T-1),dtype=float)
    opt_ri = np.zeros((pop,T-1),dtype=float)
    state = np.zeros((pop,T-1),dtype=float)

    ### initializing arrays for rerserves and functions forward dynamics ####
    roo = np.zeros((pop,T-1),dtype=float)   # gives fecundity
    rii = np.zeros((pop,T-1),dtype=float)  # immune system , proxy for immunocompetence
    rr = np.zeros((pop,T-1),dtype=float)   # total reserves

    fertility_current = np.zeros((pop,T-1),dtype=float)
    Cum_fertility = np.zeros((pop,T-1),dtype=float)
    survivability_current = np.zeros((pop,T-1),dtype=float)
    survivability_future = np.zeros((pop,T-1),dtype=float)
    Cum_survivability = np.zeros((pop,T-1),dtype=float)


    s = np.ones((pop,T-1),dtype=float)### checks if the individual is alive or not; s=1 means alive by default; updates can cause s =0 meaning dead
    #print ('s is', s)


    #the range of values of the state variable inclusively between the critical value and the capacity value
    r_range=range(r_min,C+1)
    #range function acts like interval notation [m,n) w/ defalt step size 1 (an optional 3rd argument can be added to specify the step size)

    #initializing the fittness function at each state x
    Fofr=[]
    #print ('Fofr is' , Fofr)
    #filling Fofx with the End condition (t=T) fittness values
    for r in r_range:
        Fofr.append(0)

    #printing t=T fitnesses
    #print('')
    #print('________ t=%2d ________________'%T)
    #print('State   F(x,t)      r_i*  r_o*')
    counter=0
    for endfitness in Fofr: #I use the counter to index into the the x_rangre list at the appropriate place
        #print ('%1.4f     %1.4f      %s  %s'%(r_range[counter]*( r_max /C), endfitness,'N/A', 'N/A') )
        counter=counter+1

    subPltIndex=1
    #Fofr is a list of fitnesses based on a particular state the fitness list is in the same order as the state list

    R_I = np.zeros((l,TT-1,C),dtype=float)
    R_O = np.zeros((l,TT-1,C),dtype=float)

    #### BACKWARD
    #iterating from t=T-1 to t=1 in steps of -1
    for t in range(TT-1, 0, -1):
        #print ('t',t)
        time = t

        #print('g is', g)

        maxW=[]

        for r in r_range[1:]:
            #print('r',r)

            MAX = np.zeros((l,C+1,C+1),dtype=float)
            for ri in range (0,r+1):
                #print('ri',ri)
                R = r - ri
                #print('R',R)

                for ro in range(0,R+1):
                    #print('ro',ro)
                    rt = r- ro -ri
                    #print ('rt is', rt)

                    r_prime = np.arange(l,dtype=float) # r' for each g in gaussian distribution
                    r_prime = rt + g

                    #print(r_prime)
                    for p in range(0,l):
                        if r_prime[p] >C: r_prime[p] = C
                        if r_prime[p] <0: r_prime[p] = 0

                    #print ('r_prime is ', r_prime)
                    #print('r',r)
                    #print('ri',ri)
                    #print('ro',ro)
                    RO = ro * r_max /C
                    RI = ri * r_max /C
                    RT = rt * r_max /C
                    R_PRIME =  r_prime * r_max /C

                    Max = np.arange(l,dtype=float) # Max value for each g in gaussian distribution
                    for p in range (0,l):
                        Dri=[]
                        Dro=[]
                        MAX[p,ri,ro] =  ( ( phi_max* F(N,RO,alpha) * S(M,b,gamma,RI) )  +  S_future(M,gamma,RI) * W_prime(r_prime[p]-1)   )
                        #print('MAX is',MAX)
                        ind = np.unravel_index(np.argmax(MAX[p,:,:], axis=None), MAX.shape)
                        #print('ind is',ind)
                        dro = ind[2]
                        dri = ind[1]
                        DRO = dro * r_max /C
                        DRI = dri * r_max /C
                        R_I[p][t-1][r-1] = DRI
                        R_O[p][t-1][r-1] = DRO
                        #print (t-1)
                        #print(r-1)
                        Max[p] =  np.max(MAX[p,:,:])


            W = 0
            Max = prob * Max
            W = np.sum(Max)
            #print('W is', W)
            maxW.append(W)
            #print('maxW',maxW)

            #print(t)
            #print(r)
        #housekeeping to keep the indicies of the updated version of Fofx the same as the original version
        #print('Fofr0',Fofr)
        Fofr=maxW
        #print('Fofr is',Fofr)

    ##### FORWARD ##########
    for i in range (0,pop):  ### for every nth individual
        #print('i or ind is ',i)

        time = np.arange(T-1)+1  ### time range;  reproductive lifetime
        ### assigning inital reserves STATES to initial time.
        init_res = g_func(0) ### initial reserves
        #print ('g_func is', init_res)
        state[:,0] = init_res
        #print ('state[:,0]', state[:,0])
        #### start of forward routine #####
        for t in range(0,T-1):
            #print('t is', t)
            #print('t-1 is', t-1)
            #print ('s[i,t] is', s[i,t])
            if(s[i,t]==1.0):
                #print ('g is', g)
                prob_random = np.random.random()  ### generates a random number in range [0,1] which gives the prob to go to g_h
                #print ('prob_random is', prob_random)

                prob_gen = np.cumsum(prob)

                for p in range(1,l):
                    #print (p-1, p)
                    if((prob_random>prob_gen[p-1])and(prob_random<prob_gen[p])):
                        g_gain= g[p]
                        p_g = p



                if(t>0): state[i,t] = state[i,t-1] - opt_ro[i,t-1] - opt_ri[i,t-1]  + g_gain

                if (state[i,t]>C):
                    state[i,t] = C

                if (state[i,t]<0):
                    state[i,t] = 0

                #print('state[i,t] is', state[i,t])

                #print('opt_ro[i,t] is', opt_ro[i,t])
                #print('opt_ri[i,t] is', opt_ri[i,t])

                rr[i,t] = state[i,t] *r_max/C
                #print('r[i,t] is',rr[i,t])

                R = int(rr[i,t] *C/r_max)
                #print ('R is', R)

                roo[i,t] = R_O[p_g][t][R-1]
                rii[i,t] = R_I[p_g][t][R-1]
                #print('ro[i,t] is',roo[i,t])
                #print('ri[i,t] is ',rii[i,t])
                opt_immunity[B,i,t] = R_I[p_g][t][R-1]
                opt_offspring [B,i,t] = R_O[p_g][t][R-1]


                survivability_current[i,t] = S(M,b,gamma,rii[i,t])
                survivability_future[i,t] = S_future(M,gamma,rii[i,t])
                Cum_survivability[i,t] = np.power(survivability_future[i,t],t)

                fertility_current[i,t] = F(N,roo[i,t],alpha) * S(M,b,gamma,rii[i,t])
                Cum_fertility[i,0] = Cum_survivability[i,0] * F(N,roo[i,0],alpha)
                Cum_fertility[i,t] = Cum_fertility[i,t-1] + ( Cum_survivability[i,t] * F(N,roo[i,t],alpha) )



                opt_ro[i,t] = roo[i,t] * C/r_max
                opt_ri[i,t] = rii[i,t] * C/r_max
                #print('opt_ro[i,t] updated is', opt_ro[i,t])
                #print('opt_ri[i,t] updated is', opt_ri[i,t])

                #print ('survivability_current[i,t] is ', survivability_current[i,t])
                x = np.random.random()  ### generates a random number in range [0,1] which gives the survival prob
                #print ('x is', x)
                #if (survivability_future[i,t]<x):  ### individual is dead; will not perform any funtions like resource allocation etc
                    #s[i,t+1:] = 0
        #print ('s is',s)



    for t in range(0,T-1):
        #alive[B][t] =  np.sum(s[:,t])
        sur[B][t] =  np.mean(Cum_survivability[:,t])



#######################################

#########  PLOTS ########
print("M is",M)
print("gamma is",gamma)

for B in range (0,lim):
    immunity1 = opt_immunity[B,:,:]
    m1 = immunity1.mean(0)
    v1 = immunity1.var(0)
    stdv1 = immunity1.std(0)
    upper1 = m1+v1
    lower1 = m1-v1


    offspring1 = opt_offspring[B,:,:]
    m1_o = offspring1.mean(0)
    v1_o = offspring1.var(0)
    stdv1_o = offspring1.std(0)
    upper1_o = m1_o+v1_o
    lower1_o = m1_o-v1_o

    #m_o = np.array([m1_o,m2_o,m3_o])
    #m_o = m_o.mean(0)

    #plt.figure()

    fig, ax1 = plt.subplots()

    ax2 = ax1.twinx()

    #ax1.plot(time,upper1_o,lw = 6,color='lightgrey',zorder=1)
    #ax1.plot(time,lower1_o,lw = 6,color='lightgrey',zorder=2)
    #ax1.fill_between(time,upper1_o,lower1_o,color='lightgrey',zorder=3)
    #ax1.plot(time,m1_o,lw = 2,color='darkgreen',zorder=7)
    plt.errorbar(time,m1_o,color='darkgreen',yerr=stdv1_o, fmt='o',ecolor='darkgreen',markersize=5.5, capsize=7)#uplims=True, lolims=True)
    ####
    #ax1.plot(time,upper1,lw = 6,color='lightgrey',zorder=4)
    #ax1.plot(time,lower1,lw = 6,color='lightgrey',zorder=5)
    #ax1.fill_between(time,upper1,lower1,color='lightgrey',zorder=6)
    #ax1.plot(time,m1,lw = 2,color='darkorange',zorder=8)
    plt.errorbar(time,m1,color='darkorange',yerr=stdv1, fmt='o',ecolor='darkorange',markersize=5.5, capsize=7)#uplims=True, lolims=True)

    ax2.scatter(time,sur[B,:],lw = 1,color='black',label='$f=%s$'%brooding_period[B],zorder=9)

    print("f is",brooding_period[B])
    print("rii is",round(m1[2],2))
    print("roo is",round(m1_o[2],2))
    print("survival is",sur[B,2])
    ax1.tick_params(axis="y", labelsize=20)
    ax2.tick_params(axis="y", labelsize=20)

    ax1.tick_params(axis="x", labelsize=20)
    ax2.tick_params(axis="x", labelsize=20)

    #plt.xticks(fontsize = 50)

    #plt.legend()
    plt.title('M = %s, gamma=%s,b=%0.1f'%(M,gamma,brooding_period[B]))
    ax1.set_ylim([-0.05,1.05])
    ax1.set_xlim([0.05,15.05])


    ax2.set_ylim([-0.05,1.05])
    ax2.set_xlim([0.05,15.05])

    #plt.xlabel('time',fontsize=18)
    #plt.ylabel('$r_i$',fontsize=18)
    #save = '../vandana/Dimorphism/figures/Stochastic_ri_b_%s_gamma_%s_M_%s_AllAlive.pdf' %(brooding_period[B],gamma,M)
    #save = '../vandana/Dimorphism/figures/Stochastic_ri_b_%s_gamma_%s_M_%s.pdf' %(brooding_period[B],gamma,M)
    save = '/Users/vandana/Dimorphism/figures/Stochastic_ri_f_%s_gamma_%s_M_%s_AllAlive.pdf' %(brooding_period[B],gamma,M)
    #save = 'Stochastic_ri_b_%s_gamma_%s_M_%s.pdf' %(brooding_period[B],gamma,M)
    plt.savefig(save,format='pdf')
    plt.show()
