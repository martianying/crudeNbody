#AIM

#---try to plot the histogram of total energy vs particle number at a target time.---#


#---to import everything before getting started---#

import numpy as np
import matplotlib.pyplot as plt
import random
from itertools import repeat


#partical number
pn=60

#time steps
N=70

#intervals
t=0.01

#softening re
re=0.1


#---define functions---#


#FUNCTION ONE#

#create a function in order to get the new x position array of each particle later 
def add(array1,array2):
    return np.array(array1.tolist()+array2.tolist())


T=np.array([])
while len(T)<pn*3:
    t1=np.random.random_sample(3)-1/2
    if (t1[0])**2+(t1[1])**2+(t1[2])**2>1/2:
        T=T
    else:
        T=add(T,t1)
NtriPosition=T.reshape(pn,3)
#-------------------------------------------------------------copy and paste the code above--------------------------------------------------------#

#---parameters---#


#bins
bins1=18

#target times choices
target_time1=1
target_time2=13
target_time3=40



#---define functions---#

#FUNCTION TWO#

#to define a function from the x-y-z position information(3*pn_elments) to the potential energy under this position.(pn_elements)

def PtoP(A):
    b=0
    N=[]
    for i in range(pn):
        for j in range(pn):
            if np.array_equal(A[i],A[j])==True:
                b=b
            else:
                b=b-1/(((A[i][0]-A[j][0])**2+(A[i][1]-A[j][1])**2+(A[i][2]-A[j][2])**2+re)**1/2)
        N=N+list([b])
        b=0
    return np.array(N)
#both types of A and N are arrays.



#FUNCTION THREE#

#define a function in order to get the 3N-acceleration matrix from the 3N-position matrix
def PtoA(A):
    b=0
    N=[]
    n=0
    for i in range(pn):
        for n in range(3):
            for j in range(pn):
                if np.array_equal(A[i],A[j])==True:
                    b=b
                else:
                    #b=b+(A[j][n]-A[i][n])/((A[j][0]-A[i][0])**2+(A[j][1]-A[i][1])**2+(A[j][2]-A[i][2])**2)**(3/2)
                    b=b-(A[i][n]-A[j][n])/((A[i][0]-A[j][0])**2+(A[i][1]-A[j][1])**2+(A[i][2]-A[j][2])**2+re)**(3/2)
            
            N=N+list([b])
            b=0
    N=np.array(N).reshape(pn,3)
    return N


#---main body---#


#the PLOT-LOOP

for time in [target_time1,target_time2,target_time3]:
    NtriVelocity=np.zeros((pn,3)) # initial state of velocity of n particles in 3 directions
    C=[]


#the initial x position of each particle
    Xpool=NtriPosition[:,0]

#to get the initial accelerations of three directions of each particle
    NtriA=PtoA(NtriPosition)
    Vpool=np.repeat(0,3*pn).tolist()
    full=T.tolist()
    for n in range(N):
        for i in range(pn):
            for j in range(3):
                NtriPosition[i][j]=NtriPosition[i][j]+NtriVelocity[i][j]*t
                NtriVelocity[i][j]=NtriVelocity[i][j]+1/2 *(NtriA[i][j]*t) 
                C=C+list([NtriPosition[i][j]])
                Vpool=Vpool+list([NtriVelocity[i][j]])
# the form of C after finishing n=0(when t=0*0.01) is like [no.1x,no.1y,no.1z,no.2x,no.2y,no2.z......no40.z]
#the new position matrix of each particle 
        NtriPosition=np.array(C).reshape(pn,3)
#the new x-component of position array of each particle
        Xpool=add(Xpool,NtriPosition[:,0])
    

# reshape the listC in order to get the new 3N-acceleration matrix.
        NtriA=PtoA(NtriPosition)
        full=full+C
        C=[]  # clearing the data in matrix C to start the next round
    
    
    Xpool=np.array(Xpool).reshape(N+1,pn)

    Vpool=np.array(Vpool).reshape(N+1,3*pn)

#we must initialize the value of NtriPosition
    NtriPosition=T.reshape(pn,3)


#---to get the potential energy---#


#to get the potential energy from the position imformation
    us=np.array(full).reshape(N+1,3*pn)
    transit=us[time].reshape(pn,3)
    potential_energy=PtoP(transit)

#----to get the kinetic energy---#


    kinetic_energy_xyz=Vpool[time]

    kinetic_energy_x=np.array([])
    kinetic_energy_y=np.array([])
    kinetic_energy_z=np.array([])

    for i in range(0,3*pn,3):
        kinetic_energy_x=add(kinetic_energy_x,np.array([kinetic_energy_xyz[i]]))
    for i in range(1,3*pn,3):
        kinetic_energy_y=add(kinetic_energy_y,np.array([kinetic_energy_xyz[i]]))
    for i in range(2,3*pn,3):
        kinetic_energy_z=add(kinetic_energy_z,np.array([kinetic_energy_xyz[i]]))

#the array of each particle's kinetic energy
    kinetic_energy=(1/2)*(kinetic_energy_x**2+kinetic_energy_y**2+kinetic_energy_z**2)


#---to get the total energy---#


    values=kinetic_energy+potential_energy   #to get the plot of kinetic energy or potential energy, CHANGE the defanation of values here!


#plot the histogram of the total energy and particle number.
    if time==target_time1:
        plt.hist(values,bins=bins1,color='purple',histtype='step')
        plt.xlabel("total energy")
        plt.ylabel("particle number")
    elif time==target_time2:
        plt.hist(values,bins=bins1,color='red',histtype='step')
        plt.xlabel("total energy")
        plt.ylabel("particle number")
    else:
        plt.hist(values,bins=bins1,color='pink',histtype='step')
        plt.xlabel("total energy")
        plt.ylabel("particle number")
    NtriPosition=T.reshape(pn,3)
plt.show()

