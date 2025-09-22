###############################################
# 3D RANDOM WALK
# Adding Bias to the system
# Running Comparison on the basis of walkers
###############################################
import math
import numpy as np
from random import random
import matplotlib.pyplot as plt

def func1(x,m):
    return m*x
def makehist(p,b):
    plt.hist(p, bins=b, edgecolor='black')
    plt.xlabel('Value in R')
    plt.ylabel('Frequency')
    plt.title('Histogram of R')
    plt.show()
def _sum(arr, index):
    return sum(arr[:index+1])
# def gauss(prob,std=0.01): # THIS IS THE **NEW** SELECTIVE GAUSSIAN DISTRIBUTION
#     p = 0.0
#     prob2 = prob - ((std**4/(prob**3))/2)
#     while True:
#         p = np.random.normal(prob2,std)
#         if (p > 0):
#             break 
#     return p
def gauss(desired_mean,desired_std_dev): # THIS IS ACTUALLY THE LOGNORMAL DISTRIBUTION
    mu = np.log(desired_mean) - 0.5 * np.log((desired_std_dev/desired_mean)**2 + 1)
    sigma = np.sqrt(np.log((desired_std_dev/desired_mean)**2 + 1))
    value = np.random.lognormal(mean=mu, sigma=sigma)
    return value
# def gauss(prob,std=0.01): # THIS IS THE SELECTIVE GAUSSIAN DISTRIBUTION
#     p = 0.0
#     while True:
#         p = np.random.normal(prob,std)
#         if (p > 0):
#             break 
#     return p
def covariance(x, y):
    x = np.array(x)
    y = np.array(y)
    mean_x = np.mean(x)
    mean_y = np.mean(y)
    sub_x = x - mean_x
    sub_y = y - mean_y
    numerator = np.dot(sub_x, sub_y)
    denominator = len(x) - 1 if len(x) > 1 else 1
    cov = numerator / denominator
    return cov
def trimarr(matrix,sa,wa):
    matrixtemp = np.array(matrix,dtype=float)
    matrixtemp = np.delete(matrixtemp, slice(sa,-1), 1)
    matrixtemp = np.delete(matrixtemp, slice(wa,-1), 0)
    if (sa != steps):
        matrixtemp = np.delete(matrixtemp,-1,1)
    if (wa != walkers):
        matrixtemp = np.delete(matrixtemp,-1,0)
    return matrixtemp
def neighbor(q,tz,bz,ma=2):
    narr = np.zeros([6],dtype=int)
    narr[0] = 1 if (Lattice[q[0]-1][q[1]][q[2]] == 0) else 0
    narr[1] = 1 if (Lattice[q[0]+1][q[1]][q[2]] == 0) else 0
    narr[2] = 1 if (Lattice[q[0]][q[1]-1][q[2]] == 0) else 0
    narr[3] = 1 if (Lattice[q[0]][q[1]+1][q[2]] == 0) else 0
    narr[4] = 1 if ((Lattice[q[0]][q[1]][q[2]-1] == 0) and (q[2] > bz)) else 0
    narr[5] = 1 if ((Lattice[q[0]][q[1]][q[2]+1] == 0) and (q[2] < tz)) else 0
    return narr
def move(q,rarrold,tz,bz):
    rsum = np.sum(rarrold)
    rarr = rarrold/rsum
    mov = random()
    movement = 0
    if (q[2] == bz):
        q = initpoint-1
        movement = 4
    else :
        if (0 <= mov < _sum(rarr,0)): #-x
            q[0] -= 1
            movement = 0
        elif (_sum(rarr,0) <= mov < _sum(rarr,1)): #+x
            q[0] += 1
            movement = 1
        elif (_sum(rarr,1) <= mov < _sum(rarr,2)): #-y
            q[1] -= 1
            movement = 2
        elif (_sum(rarr,2) <= mov < _sum(rarr,3)) : #+y
            q[1] += 1
            movement = 3
        elif (_sum(rarr,3) <= mov < _sum(rarr,4)): #-z
            q[2] -= 1
            movement = 4
        elif (_sum(rarr,4) <= mov < _sum(rarr,5)) : #+z
            q[2] += 1
            movement = 5
    return q,movement
def cuboidpore(poresidex,poresidey,poresidez):
    pcube = np.zeros([(poresidex+2),(poresidey+2),(poresidez+2)],dtype=int)
    pcube[:,:,0] = 1
    pcube[:,:,-1] = 1
    pcube[:,0,:] = 1
    pcube[:,-1,:] = 1
    pcube[0,:,:] = 1
    pcube[-1,:,:] = 1
    return pcube # len(pcube) = poresidex, len(pcube)[0] = poresidey, len(pcube)[0][0] = poresidez
def makepore(L,v):
    f = open('test{0}/PoreCoordinates.xyz'.format(v),'w')
    porecoord = np.zeros([3])
    for li1 in range (len(L)):
        for li2 in range (len(L[0])):
            for li3 in range (len(L[0][0])):
                if (Lattice[li1][li2][li3] == 0):
                    f.write("H   {0}   {1}   {2} \n".format(li1,li2,li3))
                    newrow = [li1, li2, li3]
                    porecoord = np.vstack([porecoord, newrow])
    f.close()
    porecoord = np.delete(porecoord, 0, 0) #removing the first zeros line
    porecoord = np.unique(porecoord, axis=0) #getting unique values only
    return porecoord
def makeRpore(L,v):
    f = open('test{0}/RPoreCoordinates.xyz'.format(v),'w')
    porecoord = np.zeros([4])
    for li1 in range (len(L)-1):
        for li2 in range (len(L[0])-1):
            for li3 in range (len(L[0][0])-1):
                if (Lattice[li1][li2][li3] == 0):
                    f.write("H   {0}   {1}   {2}   {3}\n".format(li1,li2,li3,L[li1][li2][li3]))
                    newrow = [li1, li2, li3, L[li1][li2][li3]]
                    porecoord = np.vstack([porecoord, newrow])
    f.close()
    porecoord = np.delete(porecoord, 0, 0) #removing the first zeros line
    porecoord = np.unique(porecoord, axis=0) #getting unique values only
    return porecoord
def grampore():
    f = open('GramPore.xyz','r')
    Pore = np.zeros([3])
    lines = f.readlines()
    for i, line in enumerate(lines):
        vline = line.rstrip()
        words = vline.split()
        if (words[0] == 'H'):
            ycoord = int(float(words[1]))
            zcoord = int(float(words[2]))
            xcoord = int(float(words[3]))
            # Lattice[xcoord][ycoord][zcoord] = 0
            newrow = [xcoord, ycoord, zcoord]
            Pore = np.vstack([Pore, newrow])
    f.close()
    Pore = np.delete(Pore, 0, 0) #removing the first zeros line
    Pore = np.unique(Pore, axis=0) #getting unique values only
    xmax = int(np.amax(Pore[:,0]))
    ymax = int(np.amax(Pore[:,1]))
    zmax = int(np.amax(Pore[:,2]))
    xmin = int(np.amin(Pore[:,0]))
    ymin = int(np.amin(Pore[:,1]))
    zmin = int(np.amin(Pore[:,2]))
    LX = xmax-xmin+2
    LY = ymax-ymin+2
    LZ = zmax-zmin+2
    Lattice = np.full((int(LX),int(LY),int(LZ)), 2, dtype=int)
    for i in range (len(Pore)):
        Lattice[int(Pore[i][0]-xmin-1)][int(Pore[i][1]-ymin-1)][int(Pore[i][2]-zmin-1)] = 0
    return Lattice,LX,LY,LZ
def navpore():
    # #get coords of point
    filename = "ionpore3.pdb"
    LX = 160
    LY = 165
    LZ = 122
    Lattice = np.full((int(LX),int(LY),int(LZ)), 2, dtype=int)
    NAtoms = 0
    XYZFile = np.zeros([3])
    f = open(filename,'r')
    lines = f.readlines()
    for i, line in enumerate(lines):
        vline = line.rstrip()
        words = vline.split()
        if (words[0] == 'ATOM'):
            xcoord = math.floor(float(vline[30:38]))
            ycoord = math.floor(float(vline[38:46]))
            zcoord = math.floor(float(vline[46:54]))
            Lattice[xcoord][ycoord][zcoord] = 1
            newrow = [xcoord, ycoord, zcoord]
            XYZFile = np.vstack([XYZFile, newrow])
            NAtoms += 1
    f.close()
    XYZFile = np.delete(XYZFile, 0, 0) #removing the first zeros line
    XYZFile = np.unique(XYZFile, axis=0) #getting unique values only
    # Read the Pore Details
    f = open('ILPore.xyz','r')
    Pore = np.zeros([3])
    lines = f.readlines()
    for i, line in enumerate(lines):
        vline = line.rstrip()
        words = vline.split()
        if (words[0] == 'H'):
            xcoord = int(float(words[1]))
            ycoord = int(float(words[2]))
            zcoord = int(float(words[3]))
            Lattice[xcoord][ycoord][zcoord] = 0
            newrow = [xcoord, ycoord, zcoord]
            Pore = np.vstack([Pore, newrow])
    f.close()
    Pore = np.delete(Pore, 0, 0) #removing the first zeros line
    Pore = np.unique(Pore, axis=0) #getting unique values only
    return Lattice,LX,LY,LZ
def prob(scenario,sr):
    parr = np.zeros([6])
    for pir in range (6):
        if (pir == 4):
            if (scenario == 1): parr[pir] = sr
            elif (scenario == 2): parr[pir] = sr*2
            elif (scenario == 3): parr[pir] = gauss(sr,(0.01*sr))
            elif (scenario == 4): parr[pir] = gauss(sr,(0.1*sr))
            elif (scenario == 5): parr[pir] = gauss(sr*2,(0.01*sr))
            elif (scenario == 6): parr[pir] = gauss(sr*2,(0.1*sr))
            elif (scenario == 7): parr[pir] = gauss(sr,(0.25*sr))
            elif (scenario == 8): parr[pir] = gauss(sr*2,(0.25*sr))
            elif (scenario == 9): parr[pir] = gauss(sr,(0.5*sr))
            elif (scenario == 10): parr[pir] = gauss(sr*2,(0.5*sr))
            elif (scenario == 11): parr[pir] = gauss(sr,(0.75*sr))
            elif (scenario == 12): parr[pir] = gauss(sr*2,(0.75*sr))
            elif (scenario == 13): parr[pir] = gauss(sr,sr)
            elif (scenario == 14): parr[pir] = gauss(sr*2,sr)
            elif (scenario == 15): parr[pir] = sr*2
            elif (scenario == 16): parr[pir] = gauss(sr*2,(0.01*sr))
            elif (scenario == 17): parr[pir] = gauss(sr*2,(0.1*sr))
            elif (scenario == 18): parr[pir] = gauss(sr*2,(0.25*sr))
            elif (scenario == 19): parr[pir] = gauss(sr*2,(0.5*sr))
            elif (scenario == 20): parr[pir] = gauss(sr*2,(0.75*sr))
            elif (scenario == 21): parr[pir] = gauss(sr*2,sr)
            elif (scenario == 22): parr[pir] = sr*3
            elif (scenario == 23): parr[pir] = gauss(sr*3,(0.01*sr))
            elif (scenario == 24): parr[pir] = gauss(sr*3,(0.1*sr))
            elif (scenario == 25): parr[pir] = gauss(sr*3,(0.25*sr))
            elif (scenario == 26): parr[pir] = gauss(sr*3,(0.5*sr))
            elif (scenario == 27): parr[pir] = gauss(sr*3,(0.75*sr))
            elif (scenario == 28): parr[pir] = gauss(sr*3,sr)
            elif (scenario == 29): parr[pir] = sr*3
            elif (scenario == 30): parr[pir] = gauss(sr*3,(0.01*sr))
            elif (scenario == 31): parr[pir] = gauss(sr*3,(0.1*sr))
            elif (scenario == 32): parr[pir] = gauss(sr*3,(0.25*sr))
            elif (scenario == 33): parr[pir] = gauss(sr*3,(0.5*sr))
            elif (scenario == 34): parr[pir] = gauss(sr*3,(0.75*sr))
            elif (scenario == 35): parr[pir] = gauss(sr*3,sr)
            elif (scenario == 36): parr[pir] = sr*4
            elif (scenario == 37): parr[pir] = gauss(sr*4,(0.01*sr))
            elif (scenario == 38): parr[pir] = gauss(sr*4,(0.1*sr))
            elif (scenario == 39): parr[pir] = gauss(sr*4,(0.25*sr))
            elif (scenario == 40): parr[pir] = gauss(sr*4,(0.5*sr))
            elif (scenario == 41): parr[pir] = gauss(sr*4,(0.75*sr))
            elif (scenario == 42): parr[pir] = gauss(sr*4,sr)
        elif (pir == 5):
            if (scenario == 1): parr[pir] = sr
            elif (scenario == 2): parr[pir] = sr/2
            elif (scenario == 3): parr[pir] = gauss(sr,(0.01*sr))
            elif (scenario == 4): parr[pir] = gauss(sr,(0.1*sr))
            elif (scenario == 5): parr[pir] = gauss(sr/2,(0.01*sr))
            elif (scenario == 6): parr[pir] = gauss(sr/2,(0.1*sr))
            elif (scenario == 7): parr[pir] = gauss(sr,(0.25*sr))
            elif (scenario == 8): parr[pir] = gauss(sr/2,(0.25*sr))
            elif (scenario == 9): parr[pir] = gauss(sr,(0.5*sr))
            elif (scenario == 10): parr[pir] = gauss(sr/2,(0.5*sr))
            elif (scenario == 11): parr[pir] = gauss(sr,(0.75*sr))
            elif (scenario == 12): parr[pir] = gauss(sr/2,(0.75*sr))
            elif (scenario == 13): parr[pir] = gauss(sr,sr)
            elif (scenario == 14): parr[pir] = gauss(sr/2,sr)
            elif (scenario == 15): parr[pir] = sr
            elif (scenario == 16): parr[pir] = gauss(sr,(0.01*sr))
            elif (scenario == 17): parr[pir] = gauss(sr,(0.1*sr))
            elif (scenario == 18): parr[pir] = gauss(sr,(0.25*sr))
            elif (scenario == 19): parr[pir] = gauss(sr,(0.5*sr))
            elif (scenario == 20): parr[pir] = gauss(sr,(0.75*sr))
            elif (scenario == 21): parr[pir] = gauss(sr,sr)
            elif (scenario == 22): parr[pir] = sr
            elif (scenario == 23): parr[pir] = gauss(sr,(0.01*sr))
            elif (scenario == 24): parr[pir] = gauss(sr,(0.1*sr))
            elif (scenario == 25): parr[pir] = gauss(sr,(0.25*sr))
            elif (scenario == 26): parr[pir] = gauss(sr,(0.5*sr))
            elif (scenario == 27): parr[pir] = gauss(sr,(0.75*sr))
            elif (scenario == 28): parr[pir] = gauss(sr,sr)
            elif (scenario == 29): parr[pir] = sr/3
            elif (scenario == 30): parr[pir] = gauss(sr/3,(0.01*sr))
            elif (scenario == 31): parr[pir] = gauss(sr/3,(0.1*sr))
            elif (scenario == 32): parr[pir] = gauss(sr/3,(0.25*sr))
            elif (scenario == 33): parr[pir] = gauss(sr/3,(0.5*sr))
            elif (scenario == 34): parr[pir] = gauss(sr/3,(0.75*sr))
            elif (scenario == 35): parr[pir] = gauss(sr/3,sr)
            elif (scenario == 36): parr[pir] = sr
            elif (scenario == 37): parr[pir] = gauss(sr,(0.01*sr))
            elif (scenario == 38): parr[pir] = gauss(sr,(0.1*sr))
            elif (scenario == 39): parr[pir] = gauss(sr,(0.25*sr))
            elif (scenario == 40): parr[pir] = gauss(sr,(0.5*sr))
            elif (scenario == 41): parr[pir] = gauss(sr,(0.75*sr))
            elif (scenario == 42): parr[pir] = gauss(sr,sr)
        else:
            if (scenario == 1): parr[pir] = sr
            elif (scenario == 2): parr[pir] = sr
            elif (scenario == 3): parr[pir] = gauss(sr,(0.01*sr))
            elif (scenario == 4): parr[pir] = gauss(sr,(0.1*sr))
            elif (scenario == 5): parr[pir] = gauss(sr,(0.01*sr))
            elif (scenario == 6): parr[pir] = gauss(sr,(0.1*sr))
            elif (scenario == 7): parr[pir] = gauss(sr,(0.25*sr))
            elif (scenario == 8): parr[pir] = gauss(sr,(0.25*sr))
            elif (scenario == 9): parr[pir] = gauss(sr,(0.5*sr))
            elif (scenario == 10): parr[pir] = gauss(sr,(0.5*sr))
            elif (scenario == 11): parr[pir] = gauss(sr,(0.75*sr))
            elif (scenario == 12): parr[pir] = gauss(sr,(0.75*sr))
            elif (scenario == 13): parr[pir] = gauss(sr,sr)
            elif (scenario == 14): parr[pir] = gauss(sr,sr)
            elif (scenario == 15): parr[pir] = sr
            elif (scenario == 16): parr[pir] = gauss(sr,(0.01*sr))
            elif (scenario == 17): parr[pir] = gauss(sr,(0.1*sr))
            elif (scenario == 18): parr[pir] = gauss(sr,(0.25*sr))
            elif (scenario == 19): parr[pir] = gauss(sr,(0.5*sr))
            elif (scenario == 20): parr[pir] = gauss(sr,(0.75*sr))
            elif (scenario == 21): parr[pir] = gauss(sr,sr)
            elif (scenario == 22): parr[pir] = sr
            elif (scenario == 23): parr[pir] = gauss(sr,(0.01*sr))
            elif (scenario == 24): parr[pir] = gauss(sr,(0.1*sr))
            elif (scenario == 25): parr[pir] = gauss(sr,(0.25*sr))
            elif (scenario == 26): parr[pir] = gauss(sr,(0.5*sr))
            elif (scenario == 27): parr[pir] = gauss(sr,(0.75*sr))
            elif (scenario == 28): parr[pir] = gauss(sr,sr)
            elif (scenario == 29): parr[pir] = sr
            elif (scenario == 30): parr[pir] = gauss(sr,(0.01*sr))
            elif (scenario == 31): parr[pir] = gauss(sr,(0.1*sr))
            elif (scenario == 32): parr[pir] = gauss(sr,(0.25*sr))
            elif (scenario == 33): parr[pir] = gauss(sr,(0.5*sr))
            elif (scenario == 34): parr[pir] = gauss(sr,(0.75*sr))
            elif (scenario == 35): parr[pir] = gauss(sr,sr)
            elif (scenario == 36): parr[pir] = sr
            elif (scenario == 37): parr[pir] = gauss(sr,(0.01*sr))
            elif (scenario == 38): parr[pir] = gauss(sr,(0.1*sr))
            elif (scenario == 39): parr[pir] = gauss(sr,(0.25*sr))
            elif (scenario == 40): parr[pir] = gauss(sr,(0.5*sr))
            elif (scenario == 41): parr[pir] = gauss(sr,(0.75*sr))
            elif (scenario == 42): parr[pir] = gauss(sr,sr)
    return parr
def makelattice(s):
    L = []
    sts = []
    # Create a finite 3D lattice and write pore file
    if (s == 1):
        sts = np.array([21,21,21],dtype=int)
        L = cuboidpore(sts[0],sts[1],sts[2])
    elif (s == 2):
        sts = np.array([67,67,67],dtype=int)
        L = cuboidpore(sts[0],sts[1],sts[2])
    elif (s == 3):
        sts = np.array([6,21,6],dtype=int)
        L = cuboidpore(sts[0],sts[1],sts[2])
    elif (s == 4):
        sts = np.zeros([3],dtype=int)
        L,sts[0],sts[1],sts[2] = grampore()
    elif (s == 5):
        sts = np.zeros([3],dtype=int)
        L,sts[0],sts[1],sts[2] = navpore()
    L = np.array(L,dtype=int)
    sts = np.array(sts,dtype=int)
    return L,sts
def initials(sys):
    tsl,bsl,sfs,sfe,ccs,cce,igs,ige = 0,0,0,0,0,0,0,0
    initp = [0,0,0]
    if (sys == 2): # C67
        tsl,bsl,sfs,sfe,ccs,cce,igs,ige = (47,3,41,37,37,21,21,9)
        initp = [33,33,50]
    elif (sys == 4): # Gram
        tsl,bsl,sfs,sfe,ccs,cce,igs,ige = (16,2,14,12,12,6,6,4)
        initp = [3,3,18]
    elif (sys == 5): # IL
        tsl,bsl,sfs,sfe,ccs,cce,igs,ige = (75,31,69,65,65,49,49,37)
        initp = [83,76,78]
    return tsl,bsl,sfs,sfe,ccs,cce,igs,ige,initp
def zones(sys,z,tw,tsl,bsl,sfs,ccs,igs,ige):
    # # tsl > sfs > sfe = ccs > cce = igs > ige > bsl
    crs = np.zeros([6],dtype=int) # tsl sfs ccs igs ige bsl
    ts,ss,te,se = np.zeros([4]),np.zeros([4],dtype=int),np.zeros([4]),np.zeros([4],dtype=int) # sf cc ig il
    if (sys == 4):
        round1 = 0
        for i0 in range(len(z)-1):
            if ((z[i0] == bsl) and (z[i0+1] >= tsl)):
                round1 = i0
                break
        zr = [(z[i1]) for i1 in range (0,round1)]
    newz = zr if (sys == 4) else z # changing for Gramicidin
    newstepend = round1 if ((sys == 4) and (round1 > 0)) else steps
    newtend = round1 if ((sys == 4) and (round1 > 0)) else -1
    for i in range(len(newz)):
        if (z[i] == (tsl)):
            if (crs[0] == 0): # start il
                ts[3] = tw[i]
                ss[3] = i
            crs[0] += 1
        elif (z[i] == sfs):
            if (crs[1] == 0): # start sf
                ts[0] = tw[i]
                ss[0] = i
            crs[1] += 1
        elif (z[i] == ccs):
            if (crs[2] == 0): # start cc, end sf
                ts[1] = tw[i]
                ss[1] = i
                te[0] = tw[i] 
                se[0] = i
            crs[2] += 1
        elif (z[i] == igs):
            if (crs[3] == 0): # start ig, end cc
                ts[2] = tw[i]
                ss[2] = i
                te[1] = tw[i] 
                se[1] = i
            crs[3] += 1
        elif (z[i] == ige): # end ig
            if (crs[4] == 0):
                te[2] = tw[i]
                se[2] = i
            crs[4] += 1            
        elif (z[i] == bsl): # end il
            if (crs[5] == 0):
                te[3] = tw[i]
                se[3] = i
            crs[5] += 1
    se = [(newstepend if ((int(se[k]) == 0) and (int(ss[k] != 0))) else se[k]) for k in range (4)]
    te = [(tw[newtend] if ((int(te[k]) == 0) and (int(ts[k] != 0))) else te[k]) for k in range (4)]
    wt = te - ts
    ps =  se - ss
    return ts,te,wt,ss,se,ps,crs
def velocity_autocorrelation2(vs,t0=1):
    # Only does v(t0).v(t0+t)
    # Default t0 = 1
    t0 = int(t0)
    vacfvs = [((vs[ivs]*vs[ivs+t0]) if ((ivs+t0) < len(vs)) else 0) for ivs in range (len(vs))]
    return vacfvs
def velocity_autocorrelation3(vs):
    num_steps = len(vs)
    vacf_values = np.zeros(num_steps)
    for t in range(num_steps):
        for dt in range(num_steps - t):
            vacf_values[t] += np.dot(vs[dt], vs[dt + t])
        vacf_values[t] /= num_steps - t
    return vacf_values
def zonevacfs(vacfversion, vacfwindow, walkers, stepmax, sfvels, ccvels, igvels, ilvels):
    sftvacf = np.zeros([4,walkers,stepmax[0]+1])# if stepmax[0] > 0 else np.zeros([4,walkers,1])
    cctvacf = np.zeros([4,walkers,stepmax[1]+1])# if stepmax[1] > 0 else np.zeros([4,walkers,1])
    igtvacf = np.zeros([4,walkers,stepmax[2]+1])# if stepmax[2] > 0 else np.zeros([4,walkers,1])
    iltvacf = np.zeros([4,walkers,stepmax[3]+1])# if stepmax[3] > 0 else np.zeros([4,walkers,1])
    for i5 in range (walkers):
        for j5 in range (3):
            if (vacfversion == 2):
                sftvacf[j5][i5] = velocity_autocorrelation2(sfvels[j5][i5],(vacfwindow*len(sfvels[j5][i5]))) if (stepmax[0] != 0) else [0]
                cctvacf[j5][i5] = velocity_autocorrelation2(ccvels[j5][i5],(vacfwindow*len(ccvels[j5][i5]))) if (stepmax[1] != 0) else [0]
                igtvacf[j5][i5] = velocity_autocorrelation2(igvels[j5][i5],(vacfwindow*len(igvels[j5][i5]))) if (stepmax[2] != 0) else [0]
                iltvacf[j5][i5] = velocity_autocorrelation2(ilvels[j5][i5],(vacfwindow*len(ilvels[j5][i5]))) if (stepmax[3] != 0) else [0]
            elif (vacfversion == 3):
                sftvacf[j5][i5] = velocity_autocorrelation3(sfvels[j5][i5]) if (stepmax[0] != 0) else [0]
                cctvacf[j5][i5] = velocity_autocorrelation3(ccvels[j5][i5]) if (stepmax[1] != 0) else [0]
                igtvacf[j5][i5] = velocity_autocorrelation3(igvels[j5][i5]) if (stepmax[2] != 0) else [0]
                iltvacf[j5][i5] = velocity_autocorrelation3(ilvels[j5][i5]) if (stepmax[3] != 0) else [0]
        sftvacf[3][i5] = sftvacf[0][i5] + sftvacf[1][i5] + sftvacf[2][i5]
        cctvacf[3][i5] = cctvacf[0][i5] + cctvacf[1][i5] + cctvacf[2][i5]
        igtvacf[3][i5] = igtvacf[0][i5] + igtvacf[1][i5] + igtvacf[2][i5]
        iltvacf[3][i5] = iltvacf[0][i5] + iltvacf[1][i5] + iltvacf[2][i5]
    return sftvacf,cctvacf,igtvacf,iltvacf
def zonemvacfs(walkers, steps, stepStart, stepEnd, stepmax, sftvacf, cctvacf, igtvacf, iltvacf):
    msfvacf = np.zeros([4,stepmax[0]+1])# if stepmax[0] > 0 else [0]
    mccvacf = np.zeros([4,stepmax[1]+1])# if stepmax[1] > 0 else [0]
    migvacf = np.zeros([4,stepmax[2]+1])# if stepmax[2] > 0 else [0]
    milvacf = np.zeros([4,stepmax[3]+1])# if stepmax[3] > 0 else [0]
    sfD,ccD,igD,ilD = np.zeros([4]),np.zeros([4]),np.zeros([4]),np.zeros([4])
    sfc,ccc,igc,ilc = 0,0,0,0
    # l1sf = open("test{0}/p7sf_{1}_{2}.dat".format(version,walkers,steps),'w')
    # l1cc = open("test{0}/p7cc_{1}_{2}.dat".format(version,walkers,steps),'w')
    # l1ig = open("test{0}/p7ig_{1}_{2}.dat".format(version,walkers,steps),'w')
    # l1il = open("test{0}/p7il_{1}_{2}.dat".format(version,walkers,steps),'w')
    for i6 in range (walkers):
        sff,ccf,igf,ilf = 0,0,0,0
        for j6 in range (steps):
            if (((stepEnd[i6][0]-stepStart[i6][0]) > 0) and (stepStart[i6][0] <= j6 <= stepEnd[i6][0])):
                msfvacf[0][sff],sfc = ((msfvacf[0][sff] + sftvacf[0][i6][sff]),(sfc+1)) if (sftvacf[0][i6][sff] != 0) else (msfvacf[0][sff],sfc)
                msfvacf[1][sff],sfc = ((msfvacf[1][sff] + sftvacf[1][i6][sff]),(sfc+1)) if (sftvacf[1][i6][sff] != 0) else (msfvacf[1][sff],sfc)
                msfvacf[2][sff],sfc = ((msfvacf[2][sff] + sftvacf[2][i6][sff]),(sfc+1)) if (sftvacf[2][i6][sff] != 0) else (msfvacf[2][sff],sfc)
                msfvacf[3][sff],sfc = ((msfvacf[3][sff] + sftvacf[3][i6][sff]),(sfc+1)) if (sftvacf[3][i6][sff] != 0) else (msfvacf[3][sff],sfc)
                # l1sf.write("{0} {1} {2} {3} {4}\n".format(j6,msfvacf[0][sff],msfvacf[1][sff],msfvacf[2][sff],msfvacf[3][sff]))
                sff += 1
            if (((stepEnd[i6][1]-stepStart[i6][1]) > 0) and (stepStart[i6][1] <= j6 <= stepEnd[i6][1])):
                mccvacf[0][ccf],ccc = ((mccvacf[0][ccf] + cctvacf[0][i6][ccf]),(ccc+1)) if (cctvacf[0][i6][ccf] != 0) else (mccvacf[0][ccf],ccc)
                mccvacf[1][ccf],ccc = ((mccvacf[1][ccf] + cctvacf[1][i6][ccf]),(ccc+1)) if (cctvacf[1][i6][ccf] != 0) else (mccvacf[1][ccf],ccc)
                mccvacf[2][ccf],ccc = ((mccvacf[2][ccf] + cctvacf[2][i6][ccf]),(ccc+1)) if (cctvacf[2][i6][ccf] != 0) else (mccvacf[2][ccf],ccc)
                mccvacf[3][ccf],ccc = ((mccvacf[3][ccf] + cctvacf[3][i6][ccf]),(ccc+1)) if (cctvacf[3][i6][ccf] != 0) else (mccvacf[3][ccf],ccc)
                # l1cc.write("{0} {1} {2} {3} {4}\n".format(j6,mccvacf[0][ccf],mccvacf[1][ccf],mccvacf[2][ccf],mccvacf[3][ccf]))
                ccf += 1
            if (((stepEnd[i6][2]-stepStart[i6][2]) > 0) and (stepStart[i6][2] <= j6 <= stepEnd[i6][2])):
                migvacf[0][igf],igc = ((migvacf[0][igf] + igtvacf[0][i6][igf]),(igc+1)) if (igtvacf[0][i6][igf] != 0) else (migvacf[0][igf],igc)
                migvacf[1][igf],igc = ((migvacf[1][igf] + igtvacf[1][i6][igf]),(igc+1)) if (igtvacf[1][i6][igf] != 0) else (migvacf[1][igf],igc)
                migvacf[2][igf],igc = ((migvacf[2][igf] + igtvacf[2][i6][igf]),(igc+1)) if (igtvacf[2][i6][igf] != 0) else (migvacf[2][igf],igc)
                migvacf[3][igf],igc = ((migvacf[3][igf] + igtvacf[3][i6][igf]),(igc+1)) if (igtvacf[3][i6][igf] != 0) else (migvacf[3][igf],igc)
                # l1ig.write("{0} {1} {2} {3} {4}\n".format(j6,migvacf[0][igf],migvacf[1][igf],migvacf[2][igf],migvacf[3][igf]))
                igf += 1
            if (((stepEnd[i6][3]-stepStart[i6][3]) > 0) and (stepStart[i6][3] <= j6 <= stepEnd[i6][3])):
                milvacf[0][ilf],ilc = ((milvacf[0][ilf] + iltvacf[0][i6][ilf]),(ilc+1)) if (iltvacf[0][i6][ilf] != 0) else (milvacf[0][ilf],ilc)
                milvacf[1][ilf],ilc = ((milvacf[1][ilf] + iltvacf[1][i6][ilf]),(ilc+1)) if (iltvacf[1][i6][ilf] != 0) else (milvacf[1][ilf],ilc)
                milvacf[2][ilf],ilc = ((milvacf[2][ilf] + iltvacf[2][i6][ilf]),(ilc+1)) if (iltvacf[2][i6][ilf] != 0) else (milvacf[2][ilf],ilc)
                milvacf[3][ilf],ilc = ((milvacf[3][ilf] + iltvacf[3][i6][ilf]),(ilc+1)) if (iltvacf[3][i6][ilf] != 0) else (milvacf[3][ilf],ilc)
                # l1il.write("{0} {1} {2} {3} {4}\n".format(j6,milvacf[0][ilf],milvacf[1][ilf],milvacf[2][ilf],milvacf[3][ilf]))
                ilf += 1
    sfD = [(np.sum((msfvacf[d1]/sfc)) if (sfc != 0) else 0) for d1 in range (4)]
    ccD = [(np.sum((mccvacf[d2]/ccc)) if (ccc != 0) else 0) for d2 in range (4)]
    igD = [(np.sum((migvacf[d3]/igc)) if (igc != 0) else 0) for d3 in range (4)]
    ilD = [(np.sum((milvacf[d4]/ilc)) if (ilc != 0) else 0) for d4 in range (4)]
    # l1sf.close()
    # l1cc.close()
    # l1ig.close()
    # l1il.close()
    return msfvacf,mccvacf,migvacf,milvacf,sfD,ccD,igD,ilD
def zonevels(walkers, steps, stepStart, stepEnd, velx, vely, velz, stepmax):
    sfvels = np.zeros([3,walkers,stepmax[0]+1]) if stepmax[0] > 0 else [0]
    ccvels = np.zeros([3,walkers,stepmax[1]+1]) if stepmax[1] > 0 else [0]
    igvels = np.zeros([3,walkers,stepmax[2]+1]) if stepmax[2] > 0 else [0]
    ilvels = np.zeros([3,walkers,stepmax[3]+1]) if stepmax[3] > 0 else [0]
    for i4 in range (walkers):
        sff,ccf,igf,ilf = 0,0,0,0
        for j4 in range (steps):
            if (((stepEnd[i4][0]-stepStart[i4][0]) > 0) and (stepStart[i4][0] <= j4 <= stepEnd[i4][0])):
                sfvels[0][i4][sff] = velx[i4][j4]
                sfvels[1][i4][sff] = vely[i4][j4]
                sfvels[2][i4][sff] = velz[i4][j4]
                sff += 1
            if (((stepEnd[i4][1]-stepStart[i4][1]) > 0) and (stepStart[i4][1] <= j4 <= stepEnd[i4][1])):
                ccvels[0][i4][ccf] = velx[i4][j4]
                ccvels[1][i4][ccf] = vely[i4][j4]
                ccvels[2][i4][ccf] = velz[i4][j4]
                ccf += 1
            if (((stepEnd[i4][2]-stepStart[i4][2]) > 0) and (stepStart[i4][2] <= j4 <= stepEnd[i4][2])):
                igvels[0][i4][igf] = velx[i4][j4]
                igvels[1][i4][igf] = vely[i4][j4]
                igvels[2][i4][igf] = velz[i4][j4]
                igf += 1
            if (((stepEnd[i4][3]-stepStart[i4][3]) > 0) and (stepStart[i4][3] <= j4 <= stepEnd[i4][3])):
                ilvels[0][i4][ilf] = velx[i4][j4]
                ilvels[1][i4][ilf] = vely[i4][j4]
                ilvels[2][i4][ilf] = velz[i4][j4]
                ilf += 1
    return sfvels,ccvels,igvels,ilvels
def makefileA(version, steparr, walkerarr, t, sd, tvacfx, tvacfy, tvacfz, tvacfr):
    Dvacfx = np.zeros([len(walkerarr),len(steparr)])
    Dvacfy = np.zeros([len(walkerarr),len(steparr)])
    Dvacfz = np.zeros([len(walkerarr),len(steparr)])
    Dvacfr = np.zeros([len(walkerarr),len(steparr)])
    for i in range (len(walkerarr)):
        for j in range(len(steparr)):
            ttemp = trimarr(t,steparr[j],walkerarr[i])
            sdtemp = trimarr(sd,steparr[j],walkerarr[i])
            vacfxtemp = trimarr(tvacfx,steparr[j],walkerarr[i])
            vacfytemp = trimarr(tvacfy,steparr[j],walkerarr[i])
            vacfztemp = trimarr(tvacfz,steparr[j],walkerarr[i])
            vacfrtemp = trimarr(tvacfr,steparr[j],walkerarr[i])
            msdt = [(sum(sdtemp[:,k1])/walkerarr[i]) for k1 in range (0,steparr[j])]
            tt = [(sum(ttemp[:,k8])/walkerarr[i]) for k8 in range (0,steparr[j])]
            vacfxt = [(sum(vacfxtemp[:,k13])/walkerarr[i]) for k13 in range (0,steparr[j])]
            vacfyt = [(sum(vacfytemp[:,k14])/walkerarr[i]) for k14 in range (0,steparr[j])]
            vacfzt = [(sum(vacfztemp[:,k15])/walkerarr[i]) for k15 in range (0,steparr[j])]
            vacfrt = [(sum(vacfrtemp[:,k16])/walkerarr[i]) for k16 in range (0,steparr[j])]
            l2 = open("test{0}/p35_A_{1}_{2}.dat".format(version,walkerarr[i],steparr[j]),'w')
            for k17 in range(steparr[j]):
                l2.write("{0:4f} {1:4f} {2:4f} {3:4f} {4:4f} {5:4f} \n".format(msdt[k17],tt[k17],vacfxt[k17],vacfyt[k17],vacfzt[k17],vacfrt[k17]))
            l2.close()
            Dvacfy[i][j] = np.sum(vacfyt)
            Dvacfz[i][j] = np.sum(vacfzt)
            Dvacfx[i][j] = np.sum(vacfxt)
            Dvacfr[i][j] = np.sum(vacfrt)
def gettrajs(version, walkers, steps, sites, t, path, frameskip):
    for tr in range (walkers):
        if ((tr % frameskip) == 0):
            ltrarc = open("test{0}/trajs/pTraj_{0}_{1}.arc".format(version,tr),'w')
            for tri in range (steps):
                ltrarc.write("1 \n    {0:5f}  {1:5f}  {2:5f}   90.000000   90.000000   90.000000 \n".format(sites[0],sites[1],sites[2]))
                ltrarc.write("1 Na+ {0} {1} {2} 352 \n".format(path[tr][tri][0],path[tr][tri][1],path[tr][tri][2],t[tr][tri]))
            ltrarc.close()
def fullvacfs(vacfversion, vacfwindow, walkers, steps, velx, vely, velz):
    tvacfx = np.zeros([walkers,steps])
    tvacfy = np.zeros([walkers,steps])
    tvacfz = np.zeros([walkers,steps])
    tvacfr = np.zeros([walkers,steps])
    dvarx = np.zeros([walkers,steps])
    dvary = np.zeros([walkers,steps])
    dvarz = np.zeros([walkers,steps])
    dvarxyz = np.zeros([walkers,steps])
    for i3 in range (walkers):
        if (vacfversion == 2):
            tvacfx[i3] = velocity_autocorrelation2(velx[i3],(vacfwindow*len(velx[i3])))
            tvacfy[i3] = velocity_autocorrelation2(vely[i3],(vacfwindow*len(vely[i3])))
            tvacfz[i3] = velocity_autocorrelation2(velz[i3],(vacfwindow*len(velz[i3])))
        elif (vacfversion == 3):
            tvacfx[i3] = velocity_autocorrelation3(velx[i3])
            tvacfy[i3] = velocity_autocorrelation3(vely[i3])
            tvacfz[i3] = velocity_autocorrelation3(velz[i3])
    return tvacfx,tvacfy,tvacfz,tvacfr
def statprob(ver,sr):
    pr = 0.0
    if (ver in (1,2,15,22,29,36)): pr = sr
    elif (ver in (3,5,16,23,30,37)): pr = gauss(sr,(0.01*sr))
    elif (ver in (4,6,17,24,31,38)):  pr = gauss(sr,(0.1*sr))
    elif (ver in (7,8,18,25,32,39)) :  pr = gauss(sr,(0.25*sr))
    elif (ver in (9,10,19,26,33,40)) :  pr = gauss(sr,(0.5*sr))
    elif (ver in (11,12,20,27,34,41)) :  pr = gauss(sr,(0.75*sr))
    elif (ver in (13,14,21,28,35,42)) :  pr = gauss(sr,(1*sr))
    return pr
def ratelattice(L,s,v,sr):
    latt = np.zeros((int(s[0])+1,int(s[1])+1,int(s[2])+1),dtype=float)
    for s1 in range (s[0]):
        for s2 in range (s[1]):
            for s3 in range (s[2]):
                if (L[s1][s2][s3] == 0):
                    latt[s1][s2][s3] = statprob(v,sr)
    return latt
def getrarr(q,RL,ver):
    tr1 = np.zeros([6])
    tr1[0],tr1[1],tr1[2],tr1[3] = (RL[int(q[0])-1][int(q[1])][int(q[2])],RL[int(q[0])+1][int(q[1])][int(q[2])],RL[int(q[0])][int(q[1])-1][int(q[2])],RL[int(q[0])][int(q[1])+1][int(q[2])])
    if (ver in (1,3,4,7,9,11,13)) : tr1[4],tr1[5] = ((RL[int(q[0])][int(q[1])][int(q[2])-1]),(RL[int(q[0])][int(q[1])][int(q[2])+1]))
    elif (ver in (2,5,6,8,10,12,14)) : tr1[4],tr1[5] = ((2*(RL[int(q[0])][int(q[1])][int(q[2])-1])),((RL[int(q[0])][int(q[1])][int(q[2])+1])/2))
    elif (ver in (15,16,17,18,19,20,21)) : tr1[4],tr1[5] = ((2*(RL[int(q[0])][int(q[1])][int(q[2])-1])),(RL[int(q[0])][int(q[1])][int(q[2])+1]))
    elif (ver in (22,23,24,25,26,27,28)) : tr1[4],tr1[5] = ((3*(RL[int(q[0])][int(q[1])][int(q[2])-1])),(RL[int(q[0])][int(q[1])][int(q[2])+1]))
    elif (ver in (29,30,31,32,33,34,35)) : tr1[4],tr1[5] = ((3*(RL[int(q[0])][int(q[1])][int(q[2])-1])),((RL[int(q[0])][int(q[1])][int(q[2])+1])/3))
    elif (ver in (36,37,38,39,40,41,42)) : tr1[4],tr1[5] = ((4*(RL[int(q[0])][int(q[1])][int(q[2])-1])),(RL[int(q[0])][int(q[1])][int(q[2])+1]))
    return tr1

vacfversion = 2 # = 0 for TSx
vacfwindow = 0.1
version = 1
system = 5
walkers = 5000
steps = 500
sitesize = 1
timestep = 1
frameskip = 100 # print trajs after skipping every frameskip steps
setrate = 30000 # 3*10^11 - this is divided separately for timestep calculation
fracarr1 = np.array([1])
fracarr2 = np.array([1])
steparr = steps*fracarr1
walkerarr = walkers*fracarr2
steparr = np.array(steparr,dtype=int)
walkerarr = np.array(walkerarr,dtype=int)

# Get Lattice and Pore
Lattice, sites = makelattice(system) # Get the topography of the Lattice, here Lattice[x][y][z] = 0 means pore
RLatt = ratelattice(Lattice,sites,version,setrate)
mainaxis = 2
Pore = makepore(Lattice,version)
# RPore = makeRpore(RLatt,version)
# Set some Global Pore Properties
topslice,botslice,SFstart,SFend,CCstart,CCend,IGstart,IGend,initpoint = initials(system)
initpoint = np.array(initpoint,dtype=int)
# Time and Steps arrays here
t = np.zeros([walkers,steps])
waittime = np.zeros([walkers,4]) # [SF,CC,IG,Total]
tStart = np.zeros([walkers,4]) # [SF,CC,IG,Total]
tEnd = np.zeros([walkers,4]) # [SF,CC,IG,Total]
permsteps = np.zeros([walkers,4],dtype=int) # [SF,CC,IG,Total]
stepStart = np.zeros([walkers,4],dtype=int) # [SF,CC,IG,Total]
stepEnd = np.zeros([walkers,4],dtype=int) # [SF,CC,IG,Total]
movs = np.zeros([6],dtype=int) # Total movements in each direction
crosses = np.zeros([walkers,6],dtype=int) # Crosses at tsl,sfs,ccs(sfe),igs(cce),ige,bsl

# listofR = [x for x in RPore[:,3] if x > 0]
# makehist(listofR,50)

# # Run the RW
path = [] #Absolute path
traj = [] #Lattice trajectory
for w in range(walkers):
    pathw = initpoint
    pos = np.array(initpoint)
    for i in range(steps):
        # r1arr = prob(version,setrate) # these are the set rates, not probabilities
        r1arr = getrarr(pos,RLatt,version)
        neighbourarr = neighbor(pos,topslice,botslice,mainaxis) # these check the neighbours
        ratearr = np.multiply(r1arr,neighbourarr) # these will give the rates to give the probabilites of movement
        ratesum = np.sum(ratearr)
        ratesum = 1 if (ratesum == 0) else ratesum
        rho = random()
        timestep = ((-1)*math.log(rho))/(ratesum/10000) ############## set rate changed to 30000 from 3, therefore time multiplied by 10000 (rate/10000)
        t[w][i] = t[w][i-1] + timestep
        pos,mov2 = move(pos,ratearr,topslice,botslice)
        movs[int(mov2)] += 1
        pathw = np.row_stack((pathw, pos))
    tStart[w],tEnd[w],waittime[w],stepStart[w],stepEnd[w],permsteps[w],crosses[w] = zones(system,pathw[1:,2],t[w],topslice,botslice,SFstart,CCstart,IGstart,IGend)
    path.append(pathw)
mtStart = [(sum(tStart[:,m1])/walkers) for m1 in range (4)]
mtEnd = [(sum(tEnd[:,m2])/walkers) for m2 in range (4)]
msStart = [(sum(stepStart[:,m3])/walkers) for m3 in range (4)]
msEnd = [(sum(stepEnd[:,m4])/walkers) for m4 in range (4)]
mWaitTime = [(sum(waittime[:,m5])/walkers) for m5 in range (4)]
mPermSteps = [(sum(permsteps[:,m6])/walkers) for m6 in range (4)]
mCrosses = [(sum(crosses[:,m7])/walkers) for m7 in range (6)]
print("start times " + str(mtStart))
print("end times " + str(mtEnd))
print("start steps " + str(msStart))
print("end steps " + str(msEnd))
print("time means " + str(mWaitTime))
print("step means " + str(mPermSteps))
print("mean crosses " + str(mCrosses))
if (vacfversion != 0):
    path = np.array(path, dtype=int)
    sd = np.zeros([walkers,steps],dtype=int)
    dispx = np.zeros([walkers,steps],dtype=int)
    dispy = np.zeros([walkers,steps],dtype=int)
    dispz = np.zeros([walkers,steps],dtype=int)
    dispx2 = np.zeros([walkers,steps],dtype=int)
    dispy2 = np.zeros([walkers,steps],dtype=int)
    dispz2 = np.zeros([walkers,steps],dtype=int)
    velx = np.zeros([walkers,steps])
    vely = np.zeros([walkers,steps])
    velz = np.zeros([walkers,steps])
    for i2 in range (walkers):
        for j2 in range (steps):
            dispx[i2][j2] = (path[i2][j2][0] - path[i2][0][0]) * 100
            dispy[i2][j2] = (path[i2][j2][1] - path[i2][0][1]) * 100
            dispz[i2][j2] = (path[i2][j2][2] - path[i2][0][2]) * 100
            dispx2[i2][j2] = dispx[i2][j2] ** 2
            dispy2[i2][j2] = dispy[i2][j2] ** 2
            dispz2[i2][j2] = dispz[i2][j2] ** 2
            sd[i2][j2] = dispx2[i2][j2] + dispy2[i2][j2] + dispz2[i2][j2]
            velx[i2][j2] = (dispx[i2][j2]/(t[i2][j2])) if (j2 != 0) else 0 # do not do the sum, use t itself
            vely[i2][j2] = (dispy[i2][j2]/(t[i2][j2])) if (j2 != 0) else 0 # do not do the sum, use t itself
            velz[i2][j2] = (dispz[i2][j2]/(t[i2][j2])) if (j2 != 0) else 0 # do not do the sum, use t itself
    # Get vacfs for the full length of time (can skip this step)
    tvacfx, tvacfy, tvacfz, tvacfr = fullvacfs(vacfversion, vacfwindow, walkers, steps, velx, vely, velz)
    # Calculate the three zones vels
    # Find the max size of array needed for each portion
    stepmax = np.zeros([walkers,4],dtype=int)
    stepmax = [(np.amax(permsteps[:,zz]) if (np.amax(permsteps[:,zz]) > 0) else 0) for zz in range (4)]
    print("stepmax " + str(stepmax))
    # create velocity for each zone
    sfvels, ccvels, igvels, ilvels = zonevels(walkers, steps, stepStart, stepEnd, velx, vely, velz, stepmax)
    # You have arrays with vels and zeros when the vel is not needed, now make the vacf array for each zone
    sftvacf, cctvacf, igtvacf, iltvacf = zonevacfs(vacfversion, vacfwindow, walkers, stepmax, sfvels, ccvels, igvels, ilvels)
    # Now calculate the mean vacfs and D's for each zone (Also makes a file for each zone)
    msfvacf,mccvacf,migvacf,milvacf,sfD,ccD,igD,ilD = zonemvacfs(walkers, steps, stepStart, stepEnd, stepmax, sftvacf, cctvacf, igtvacf, iltvacf)
    # print(sfD,ccD,igD,ilD)
    # Now make file A for bookkeeping
    makefileA(version, steparr, walkerarr, t, sd, tvacfx, tvacfy, tvacfz, tvacfr)
    # Print the Dmatrix file
    l2_3 = open("test{0}/Dmatrix_{1}_{2}.dat".format(version,walkers,steps),'w')
    # l2_3.write("{0:4f} \n{1:4f} \n{2:4f} \n{3:4f} \n{4:4f} \n{5} \n".format(sfD[0],sfD[1],sfD[2],sfD[3],mWaitTime[0],mPermSteps[0]))
    # l2_3.write("{0:4f} \n{1:4f} \n{2:4f} \n{3:4f} \n{4:4f} \n{5} \n".format(ccD[0],ccD[1],ccD[2],ccD[3],mWaitTime[1],mPermSteps[1]))
    # l2_3.write("{0:4f} \n{1:4f} \n{2:4f} \n{3:4f} \n{4:4f} \n{5} \n".format(igD[0],igD[1],igD[2],igD[3],mWaitTime[2],mPermSteps[2]))
    # l2_3.write("{0:4f} \n{1:4f} \n{2:4f} \n{3:4f} \n{4:4f} \n{5} \n".format(ilD[0],ilD[1],ilD[2],ilD[3],mWaitTime[3],mPermSteps[3]))
    l2_3.write("{0:4f} \n{1:4f} \n{2} \n".format(sfD[3],mWaitTime[0],mPermSteps[0]))
    l2_3.write("{0:4f} \n{1:4f} \n{2} \n".format(ccD[3],mWaitTime[1],mPermSteps[1]))
    l2_3.write("{0:4f} \n{1:4f} \n{2} \n".format(igD[3],mWaitTime[2],mPermSteps[2]))
    l2_3.write("{0:4f} \n{1:4f} \n{2} \n".format(ilD[3],mWaitTime[3],mPermSteps[3]))
    l2_3.close()
if (vacfversion == 0):
    l2_3 = open("test{0}/Dmatrix_{1}_{2}.dat".format(version,walkers,steps),'w')
    l2_3.write("{0:4f} \n{1:4f} \n".format(mWaitTime[0],mPermSteps[0]))
    l2_3.write("{0:4f} \n{1:4f} \n".format(mWaitTime[1],mPermSteps[1]))
    l2_3.write("{0:4f} \n{1:4f} \n".format(mWaitTime[2],mPermSteps[2]))
    l2_3.write("{0:4f} \n{1:4f} \n".format(mWaitTime[3],mPermSteps[3]))
    l2_3.close()
# Print some trajs
# gettrajs(version, walkers, steps, sites, t, path, frameskip)

