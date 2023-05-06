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

M = 4.0 # Factors causing mortality
gamma = 2 #amount of decreasing/increases mortality. Controld mortality or its rate


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

Res = np.zeros((pop,lim),dtype=float)
ROO = np.zeros((pop,lim),dtype=float)
RII = np.zeros((pop,lim),dtype=float)
alive = np.zeros((lim,T-1),dtype=float)### no. of alive individuals at a certain time for different brooding periods

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
        state[:,0] = init_res
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
                if (survivability_future[i,t]<x):  ### individual is dead; will not perform any funtions like resource allocation etc
                    s[i,t+1:] = 0

        #print ('s is',s)
        Res[i][B] = rr[i][1]
        ROO[i][B] = roo[i][1]
        RII[i][B] = rii[i][1]

    #print(s)

    for t in range(0,T-1):
        alive[B][t] =  np.sum(s[:,t])

    m = rii.mean(0)
    v = rii.var(0)
    stdv= rii.std(0)
    upper = m+v
    lower = m-v
    plt.figure()
    #for nn in range (0,pop):
        #plt.plot(time,rii[nn,:],color='darkorange')

    #plt.plot(time,upper,lw = 6,color='lightgrey',zorder=1)
    #plt.plot(time,lower,lw = 6,color='lightgrey',zorder=2)
    #plt.fill_between(time,upper,lower,color='lightgrey',zorder=3)
    #plt.scatter(time,m,lw = 0.2,color='darkorange',zorder=4)
    plt.errorbar(time,m,color='darkorange',yerr=stdv, fmt='o',ecolor='black',markersize=5.5, capsize=7)#uplims=True, lolims=True)
    plt.xticks([2,4,6,8,10,12,14] , fontsize = 20)
    plt.yticks([0.0,0.2,0.4,0.6,0.8,1.0] ,fontsize = 20)
    #plt.legend()
    plt.title('M = %s, gamma=%s,b=%0.1f'%(M,gamma,b))
    plt.ylim([-0.05,0.65])
    plt.xlabel('time',fontsize=18)
    plt.ylabel('$r_i$',fontsize=18)
    save = '/Users/vandana/Dimorphism/figures/Stochastic_ri_b_%s_AllAlive.pdf' %(b)
    plt.savefig(save,format='pdf')
    #plt.show()



    plt.figure()
    for nn in range (0,pop):
        plt.scatter(time,rii[nn,:])
    plt.xticks([2,4,6,8,10,12,14] , fontsize = 20)
    plt.yticks([0.0,0.2,0.4,0.6,0.8,1.0] ,fontsize = 20)
    #plt.legend()
    plt.title('M = %s, gamma=%s,b=%0.1f'%(M,gamma,b))
    plt.ylim([-0.05,0.65])
    plt.xlabel('time',fontsize=18)
    plt.ylabel('$r_i$',fontsize=18)
    #save = '/Users/vandana/Dimorphism/figures/Stochastic_ri_b_%s_individuals.pdf' %(b)
    #plt.savefig(save,format='pdf')
    #plt.show()


####### plots#########


for i in range (0,pop):
    #alphas = i/(10) +0.5
    #print(alphas)
    #plt.plot(brooding_period,Res[i,:],color='blue',alpha = alphas,label='r,ind=%0.1f'%i)
    #plt.plot(brooding_period,Res[i,:],label='r,ind=%0.1f'%i)
    #plt.plot(brooding_period,RII[i,:],color='darkorange',alpha = alphas,label='ri,ind=%0.1f'%i)
    plt.plot(brooding_period,RII[i,:],label='ri,ind=%0.1f'%i)
    #plt.plot(brooding_period,ROO[i,:],color='green',alpha = alphas,label='ro,ind=%0.1f'%i)
#plt.legend()
plt.title('M = %s, gamma=%s'%(M,gamma))
plt.xlabel('brooding period',fontsize=18)
plt.ylabel('reserves for immunity $r_i(r(t=1),t=1$)',fontsize=18)
save = '/Users/vandana/Dimorphism/figures/Stochastic_b_ri_ro_Mortality_%s.pdf' %M
plt.savefig(save,format='pdf')
plt.show()



#########  PLOTS ########

hue = []
hue = ['deepskyblue','royalblue', 'blue']
markers = ['.', '+', 'x']
period = np.array(['short brooding period', 'intermediate brooding period', 'long brooding period'])
##########
def func2(x, a, c):
    return (a  * x) + c
expo =[]

for i in range (0,lim):
    b = brooding_period[i]
    alphas = b+0.1#((b*10-1)/(2*(brooding_period[i]-brooding_period[0]))) + 0.5
    plt.scatter(time,np.log(alive[i,:]),color=hue[i],marker=markers[i],s=30,label=period[i])#label='$f $ = %0.1f'%b)
plt.legend()
plt.yscale('log')
plt.title('M = %s, gamma=%s'%(M,gamma))
plt.xlabel('time',fontsize=18)
plt.ylabel('# alive individuals',fontsize=18)
plt.xticks( fontsize = 20)
plt.yticks(fontsize = 20)
#save = '../vandana/Dimorphism/figures/Stochastic_alive_Mortality_%s_gamma_%s_AllAlive.pdf' %(M,gamma)
save = '/Users/vandana/Dimorphism/figures/Stochastic_alive_Mortality_%s_gamma_%s.pdf' %(M,gamma)
plt.savefig(save,format='pdf')
plt.show()


#########   ########

def func1(x, a, b, c):
    return a * np.exp(-b * x) + c

expo =[]

for i in range (0,lim):
    b = brooding_period[i]
    #alphas = b+0.1#((b*10-1)/(2*(brooding_period[i]-brooding_period[0]))) + 0.5
    #plt.scatter(time,np.log(alive[i,:]),color=hue[i],alpha = alphas,label='$f $ = %0.1f'%b)
    #plt.scatter(time,np.log(alive[i,:]),color=hue[i],alpha = alphas,label='$f $ = %0.1f'%b)
    #plt.scatter(time,np.log(alive[i,:]),color='black',marker=markers[i],label='$f $ = %0.1f'%b)
    xdata= time
    ydata = alive[i,:]/pop
    popt, pcov = curve_fit(func1, xdata, ydata)
    plt.scatter(time,((alive[i,:])/pop),color=hue[i],marker=markers[i],s=30,label='$M_{emergent}$ for '+""+period[i]+'=%0.2f'%popt[1])#label='$f $ = %0.1f'%b)
    #plt.plot(xdata, func(xdata, *popt), color=hue[i],lw=0.8,label='%0.2f*exp(-%0.2f*time)'%(popt[0],popt[1]) )
plt.legend()
plt.title('M = %s, gamma=%s'%(M,gamma))
plt.xlabel('time',fontsize=18)
plt.ylabel('Fraction of individuals alive',fontsize=18)
plt.xticks( fontsize = 20)
plt.yticks(fontsize = 20)
#save = '/Users/vandana/Dimorphism/figures/Stochastic_alive_Mortality_%s_gamma_%s_AllAlive.pdf' %(M,gamma)
save = '/Users/vandana/Dimorphism/figures/Stochastic_alive_Mortality_Scatter_%s_gamma_%s.pdf' %(M,gamma)
plt.savefig(save,format='pdf')
plt.show()
