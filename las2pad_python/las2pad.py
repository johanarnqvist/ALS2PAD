#   This software may be used, copied, or redistributed as long as it is not sold and this note is reproduced on each copy made. 
#   This routine is provided as is without any express or implied warranties whatsoever.
#   For more details See Arnqvist et al. Robust processing of airborne laser scans to plant area density profiles, 2020, biogeosciences.
# Please send a notification of any bugs to johan.arnqvist[at]geo.uu.se

# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 14:39:50 2019

@author: johar477
"""


import laspy as lp
import numpy as np

#Define a function that works like matlabs find for 1 d arrays. This is done to speed up calculations.
def f1d(x):
    return np.array(np.nonzero(x)[0])


#Use a timer to show the progress
import time

timer1=time.time()

dx=20 #Set the grid size in x [m]
dy=20 #Set the grid size in y [m]
dz=5 #Set the grid size in z [m]

hmax=40 #The highest trees to search for in m
k=0.5 #Extinction coefficient

#Load the file. Specify the path to the file you want analyse.

folderget=('Enter_folder_path')
fname='file_name.las'
A=lp.file.File(folderget+fname)

Int=np.array(A.intensity)
x=np.array(A.x)
y=np.array(A.y)
z=np.array(A.z)
C=np.array(A.classification)
Rnum=np.array(A.return_num)
Num=np.array(A.num_returns)
phi=np.array(A.scan_angle_rank)

Isum=np.array(Int) #Assign vector for the sum of the intensities

####------------Display how many points have a certain return number and create scaled intensities
rmax=np.max(Num)
for i in range(rmax):
    print(['A total of ' + np.str(np.round(np.sum(Num==(i+1))/np.size(x)*100, decimals=2)) + ' % of the pulses have number of returns = ' + str(i+1)])

#Rescale the intensities according to the total return intensity in one pulse
for i in np.arange(1,rmax):
    p=f1d((Rnum==(i+1)) & (Num==(i+1)))
    maxnum=p.size
    for ii in np.arange(i):
        #Take ii steps back to make sure that a consecutive number of returns lead up to Num
        pprev=f1d((Rnum==(i+1-ii)) & (Num==(i+1)))+ii 
        p=np.intersect1d(p,pprev)
    print([str(np.round(p.size/maxnum*100,decimals=2)) + ' % of ' + str(i+1) + ':th order returns where in order'])
    if Num[0]!=1: #Check that the start of the firs pulse is in the file
        p=p[1:]
    #Sum up the intensity over all shots
    sum1=np.zeros(p.size)
    for ii in range(i+1):
        sum1=sum1+Int[p-ii]
    #Make a vector for the total intensity of each shot (equivalent to Num)
    for ii in range(i+1):
        Isum[p-ii]=sum1

#Remove the values which have sero intensity for the entire shot
p=np.nonzero(Isum!=0)
print([str(x.size-len(p[0])) + ' values where removes because those shots had a total backscatter intensity = 0'])
x=x[p]
y=y[p]
z=z[p]
C=C[p]
Rnum=Rnum[p]
Num=Num[p]
phi=phi[p]
phi=phi/180*np.pi
Int=Int[p]
Isum=Isum[p]
Iscale=Int/Isum #This scales according to share of the total intensity. For points that do not have the order of returns sorted, the scale will be 1.

########----------------------------------------Grid the data and run the PAD routines

ll=np.floor(np.array([np.min(x), np.min(y)])) #Coordinates of lower left corner (in whole metres)
ur=np.ceil(np.array([np.max(x), np.max(y)])) #Coordinates of upper right corner (in whole metres)

sx=np.int(np.ceil((ur[0]-ll[0])/dx)) #Size in x (in number of grid points)
sy=np.int(np.ceil((ur[1]-ll[1])/dy)) #Size in y (in number of gridpoints)

gh=np.zeros(sx*sy) #Allocate variables to be computed
th=np.zeros(sx*sy)
PAI=np.zeros(sx*sy)
wtrflag=np.zeros(sx*sy)

sz=np.int(np.ceil(hmax/dz)) #Number of vertical levels
rz=range(sz) #Range in z (except ground value)

PAD=np.zeros((sx*sy,sz)) #Allocate space for the gridded average intensities (scaled if ISRR is used)
Igp=np.zeros(sz+1) #Intensities within the gridbox

grp=f1d(C==2) #Ground hits
wtp=f1d(C==9) #Water hits

timer2=time.time() #Start a timer for each completed row
for ix in range(sx):
    p1=f1d(((x-ll[0])//dx)==ix)
    p1g=f1d(((x[grp]-ll[0])//dx)==ix)
    p1g=grp[p1g]
    for iy in range(sy):
        p2=f1d(((y[p1]-ll[1])//dy)==iy)
        p=p1[p2] #Point within the gridbox
        p2=f1d(((y[p1g]-ll[1])//dy)==iy)
        pgr=p1g[p2] #Points within hte gridbox that are ground hits
        if pgr.size!=0: #There must be at least one ground hit to determine the ground height
            zg=np.median(z[pgr]) #The ground height is determined as the median of the ground level in the gridbox
            gh[sy*ix+iy]=zg
            zgb=z[p]-zg #The heights of the returns falling wihtin the gridbox
            th[sy*ix+iy]=np.max(zgb)
            Igp[0]=np.sum(Iscale[pgr])
            for iz in rz:
                Igp[iz+1]=np.sum(Iscale[p[zgb<((iz+1)*dz)]])
            P=Igp/Igp[-1]
            PAI[sy*ix+iy]=-np.mean(np.abs(np.cos(phi[p])))*np.log(P[0])/k
            PAD[sy*ix+iy,:]=-np.diff(-np.mean(np.abs(np.cos(phi[p])))*np.log(P)/k/dz)
        elif np.size(np.intersect1d(wtp,p))!=0: #If the points are water reflections, take the water height
            zg=np.min(z[p]) #The ground height is determined as the median of the ground level in the gridbox
            gh[sy*ix+iy]=zg
            PAI[sy*ix+iy]=0
            PAD[sy*ix+iy,:]=0
            wtrflag[sy*ix+iy]=1
        else:
            gh[sy*ix+iy]=np.nan
            PAI[sy*ix+iy]=np.nan
            PAD[sy*ix+iy,:]=np.nan
            wtrflag[sy*ix+iy]=np.nan
    #Print the progress
    tstamp=time.time()-timer2
    timer2=time.time()
    print(['Time to process row ' + str(ix+1) + ' is ' + str(tstamp) + ' seconds'])

#Turn the data into North as vertical and East as horizontal.

gh=np.reshape(gh,(sx,sy)).transpose()
wtrflag=np.reshape(wtrflag,(sx,sy)).transpose()
th=np.reshape(th,(sx,sy)).transpose()
PAI=np.reshape(PAI,(sx,sy)).transpose()
PAD=np.reshape(PAD,(sx,sy,sz)).transpose()

#Print the final time 
print(['Time to process file ' + fname + ' is ' + str(time.time()-timer1) + ' seconds'])

            
            
        
    
