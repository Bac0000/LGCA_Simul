import numpy as np                                             
from math import sqrt
import random as rd
from random import randint
from random import choice
from copy import deepcopy as dp

# Creates a list of indexes for a table of dimensions M x M
def indexes(M:int):
    v=[]
    for i in range(0,M):
        for j in range(0,M):
            coppia=[i,j]
            v.append(dp(coppia))
    return v

# Find for which indexes i, element l[i] of list l is equal to el
def find_indexes(l:list,el):
   v=[]
   for i in range(0,len(l)):
      if (l[i]==el):
         v.append(i)
   return v

   
# Craetes a random initial configuration for the I0 initial infected inside the cells
# of the grid. it asupres that the configuration is admissible
def Init_LGCA(N:int,S0:int,I0:int):
   i_vec=[]
   s_vec=[]
   r_vec=[]
   aus=indexes(int(sqrt(N)))
   for i in range(0,I0):
      scelta=choice(aus)
      i_vec.append(scelta)
      aus.remove(scelta)
   s_vec=aus
   return [s_vec,i_vec,r_vec]


# Defines the Neighboors of an hexagon cells, in which individual of such cell can move to
# This are the six hexagons whixh share one side with the original cell
def Neighboors(cella:list[int,int],dim:int):
   n=int(cella[1])
   m=int(cella[0])
   #print("n","m",n,m)
   #print("n",n,"m",m)
   up= (n+1) % dim
   down=(n-1) % dim
   right= (m+1) % dim
   left= (m-1) % dim
   #print(up,down,right,left,dim)
   remainder=m % 2 
   if (remainder==0):
    v2=[right,down]
    v4=[left,down]
   elif(remainder==1):
    v2=[right,up]
    v4=[left,up]
   v1=[right,n]
   v3=[left,n]
   v5=[m,up]
   v6=[m,down]
   #print("io")
   return [v1,v2,v3,v4,v5,v6]

# Run an LGCA simulation for 100 days, or until there are no more Infected
# Individuals move to a neighbooring cell at each time, all individuals in the same cell choose a different path for the next
# time step. With a proper initial distribution, this ensures that no more then six individuals share the same cell at each time step
# Infections can only happen inside the same cell, the number of infected inside a cell determin the probability
# of susceptible individuals in the same cell to get infected. Recovery may happen in all cells
# At the border of the grid, Pacman effect are taken into acount to keep the restriction of the number of individuals
# in the same cell.
def LGCA(N:int,alpha_LG:float,beta:float,S0:int,I0:int):
   Init=Init_LGCA(N,S0,I0)
   aus=indexes(int(sqrt(N)))
   s_vec=Init[0]
   i_vec=Init[1]
   r_vec=Init[2]
   i_new=[]
   r_new=[]
   s_new=[]
   stato=[[dp(s_vec),dp(i_vec),dp(r_vec)]]
   k=0
   count=0
   while(k<100 and len(i_vec)>0):
      i_new=[]
      r_new=[]
      s_new=[]
      for coppia in aus:
          I_coppia=find_indexes(i_vec,coppia)
          S_coppia=find_indexes(s_vec,coppia)
          R_coppia=find_indexes(r_vec,coppia)
          ni=len(I_coppia)
          ns=len(S_coppia)
          nr=len(R_coppia)
          numer=ni+ns+nr
          count=count+numer
          neig=Neighboors(coppia,sqrt(N))
          out=rd.sample(neig,k=numer)
          for s in range(0,ns):
            part=out[s]
            phi=(1-alpha_LG)**(ni)               # susceptibles can either get infected or stay susceptibles. 
            psi=1-phi                            # infection probability depend on beta and number of infected
            bool=np.random.binomial(n=1,p=psi)   # in the cell
            if (bool==0):
               s_new.append(dp(part))
            else:
               i_new.append(dp(part))
          for i in range(ns,ns+ni):
            parti=out[i]                         # Infected people either pass to recovered with rate gamma, 
            bool=np.random.binomial(n=1,p=beta)  # or remain infected
            if (bool==0):
               i_new.append(dp(parti))
            else:
               r_new.append(dp(parti))           # Recovered cannot pass to other states
          for r in range(ns+ni,ns+ni+nr):
            partr=out[r]
            r_new.append(dp(partr))
      stato.append([dp(s_new),dp(i_new),dp(r_new)])
      i_vec=dp(i_new)
      r_vec=dp(r_new)
      s_vec=dp(s_new)
      print(len(i_vec),len(s_vec),len(r_vec))
      k=k+1
      count=0
   return stato

# Constructs neighboors for a square grid. In general the nine neighboors of a square cells include the eight
# Surrounding cell plus itself. This means that it is posible that individuals remain fixed inside the cell
# at different times
def Square_Neighs(couple:list[int,int],dim:int):
   DIM=dim-1
   j=couple[0]
   i=couple[1]
   n1=[j+1,i+1]
   n2=[j,i+1]
   n3=[j-1,i+1]
   n4=[j+1,i]
   n5=[j,i]
   n6=[j-1,i]
   n7=[j+1,i-1]
   n8=[j,i-1]
   n9=[j-1,i-1]
   neig=[n1,n2,n3,n4,n5,n6,n7,n8,n9]
   if (j==DIM):
      neig=[n2,n3,n5,n6,n8,n9]
   if (j==0):
      neig=[n2,n1,n5,n4,n8,n7]
   if (i==0):
      neig=[n1,n2,n3,n4,n5,n6]
   if (i==DIM):
      neig=[n6,n7,n8,n4,n5,n9]
   if (j==0 and i==0):
      neig=[n1,n2,n4,n5]
   if (j==DIM and i==DIM):
      neig=[n5,n6,n8,n9]
   if (j==0 and i==DIM):
      neig=[n4,n5,n7,n8]
   if (j==DIM and i==0):
      neig=[n3,n5,n2,n6]
   return neig



# Creates a squared grid
def gen_couples(n:int,M:int):
      vec=[]
      if(n==0):
         return vec
      else:
        for ind in range(0,n):
            i=randint(0,M-1)
            j=randint(0,M-1)
            couple=[j,i]
            vec.append(couple)
        return vec
   

# Initializes a random distribution for the individuals on the grid
# Theoretically we could choose which ever distribution of individuals on the grid we prefer,
# Which was not possible in the LGCA standard model
def Init_new(M,S0,I0,R0):
   i_vec=gen_couples(I0,M)
   s_vec=gen_couples(S0,M)
   r_vec=gen_couples(R0,M)

   config=[s_vec,i_vec,r_vec]
   return config
      

# This is the modified LGCA version, movement between neighbooring cell is randomized without any restrictions
# i.e each individual a time t+1 moves to randomly choosen cell from the neighborhood, which includes the cell itself
# There is no PacMan effect here, hence in the borders of the grid the number of neighboor is reduced.
# For the epidemic diffusion, the same rule applies: only individuals in the same cell interact with each other
# Note that the beta coefficient has to be modified according to the dimension of the grid.
# Before, as the number of cells corresponded with the number of individuals, we simply had
#                         beta_LGCA = beta_SIR * N = (beta/N)*N= beta
# here instead, if we have a grid of dimension M x M, we will have
#                         beta_mod = beta_SIR * M^2 = (beta/N)*M^2
# As input _beta_ we will give beta_SIR.

def Modified_LGCA(M:int,_beta_:float,gamma:float,S0,I0,R0):
   Init=Init_new(M,S0,I0,R0)
   beta=_beta_*M*M
   grid=indexes(M)
   s_vec=Init[0]
   i_vec=Init[1]
   r_vec=Init[2]
   print(len(i_vec),len(s_vec),len(r_vec))
   i_new=[]
   r_new=[]
   s_new=[]
   stato=[[dp(s_vec),dp(i_vec),dp(r_vec)]]
   k=0
   while(k<100 and len(i_vec)>0): # 100 days estimation, if infected goes to zero, stop earlier
      i_new=[]
      r_new=[]
      s_new=[]
      for coppia in grid:
          I_coppia=find_indexes(i_vec,coppia)
          S_coppia=find_indexes(s_vec,coppia)
          R_coppia=find_indexes(r_vec,coppia)
          ni=len(I_coppia)
          ns=len(S_coppia)
          nr=len(R_coppia)
          tot=ni+ns+nr
          neig=Square_Neighs(coppia,M)
          out=rd.choices(neig,k=tot)
          for s in range(0,ns):                    # susceptibles can either get infected or stay susceptibles. 
            part=out[s]                            # infection probability depend on beta and number of infected
            phi=(1-beta)**(ni)                     # in the cell
            psi=1-phi
            bool=np.random.binomial(n=1,p=psi)
            if (bool==0):
               s_new.append(dp(part))
            else:
               i_new.append(dp(part))
          for i in range(ns,ns+ni):               # Infected people either pass to recovered with rate gamma, 
            parti=out[i]                          # or remain infected
            bool=np.random.binomial(n=1,p=gamma)
            if (bool==0):
               i_new.append(dp(parti))
            else:
               r_new.append(dp(parti))
          for r in range(ns+ni,ns+ni+nr):
            partr=out[r]
            r_new.append(dp(partr))

      stato.append([dp(s_new),dp(i_new),dp(r_new)])
      i_vec=dp(i_new)
      r_vec=dp(r_new)
      s_vec=dp(s_new)
      # print(len(i_vec),len(s_vec),len(r_vec)) # If you want to see the interactions, un-comment this line

      k=k+1
   return stato

# Number of Infected for each iteration of the LGCA model
def count_I(Stati):
    v=[]
    for el in Stati:
       I=el[1]
       v.append(len(I))
    return v

# Number of Recovered for LGCA model in each iteration
def count_R(Stati):
    v=[]
    for el in Stati:
       R=el[2]
       v.append(len(R))
    return v


# Finds number of new positives per day, to compare with real data
def find_incid(Stati):
    I=count_I(Stati)
    R=count_R(Stati)
    n=len(I)
    incid=[15]
    for i in range(1,n):
        I_today=I[i]
        R_new=R[i]-R[i-1]
        I_yest=I[i-1]-R_new
        I_new=I_today-I_yest
        incid.append(I_new)
    return incid       

