"""
This code was my university project 
Code by: Reza Lotfi Navaei

"""

#===================================
#       Developing Flow in a
#       Rectangular  Channel
#         SMIPLE Algorithm
#     Energy & Species Solution
#===================================
import time
import numpy as np
import matplotlib.pyplot as plt


#============================ function definition =============================
def Initialize(T,Y1,Y2):

    for i in range(m+2):
       for j in range(n+2):
           T[i,j]=0
           Y1[i,j]=0
           Y2[i,j]=0
     
   
def FaceFlux(FW,FE,FS,FN,DW,DE,DS,DN,NS):
    
    # centeral scheme
    if NS==1:
        aW=DW+0.5*FW
        aE=DE-0.5*FE   
        aS=DS+0.5*FS
        aN=DN-0.5*FN
    
    # upwind scheme
    elif NS==2:
        aW=DW+max( FW,0.)
        aE=DE+max(-FE,0.)     
        aS=DS+max( FS,0.)
        aN=DN+max(-FN,0.)
    
    # hybrid scheme
    elif NS==3:
        aW=max( FW,DW+0.5*FW,0.)
        aE=max(-FE,DE-0.5*FE,0.)
        aS=max( FS,DS+0.5*FS,0.)
        aN=max(-FN,DN-0.5*FN,0.)
        
    return aW,aE,aS,aN
           
        
def ComputeNorm(z,z_old,L2):
    
    L2=0.
    num=0
    for i in range(1,m):
        for j in range(1,n):
            L2+=(z[i,j]-z_old[i,j])**2
            num+=1           
    L2=(L2/num)**0.5 
        
    return L2


def UpdateGhostCell_Temperature(T):
    
    # horizontal sides
    for i in range(m+2):
        T[i,0]=2*TB-T[i,1]
        T[i,n+1]=2*TT-T[i,n]
        # vertical sides
    for j in range(1,n+1):
        T[0,j]=TI
        T[m+1,j]=T[m,j]


def UpdateGhostCell_Species1(Y1):
    
    # horizontal sides
    for i in range(m+2):
        Y1[i,0]=Y1[i,1]
        Y1[i,n+1]=Y1[i,n]
        # vertical sides
    for j in range(1,n+1):
        if j<int(n/2)+1:          
            Y1[0,j]=Y1I
        else:
            Y1[0,j]=0.
          
        Y1[m+1,j]=Y1[m,j]
     
'''    # horizontal sides
    for i in range(m+2):
        Y1[i,0]=2*0-Y1[i,1]
        Y1[i,n+1]=2*0-Y1[i,n]
        # vertical sides
    for j in range(1,n+1):
        if j*dy<0.4*h:         
            Y1[0,j]=0
        elif j*dy>0.6*h:
            Y1[0,j]=0.
        else:
            Y1[0,j]=1.
              
        Y1[m+1,j]=Y1[m,j]
'''        
        
def UpdateGhostCell_Species2(Y2):
    
    # horizontal sides
    for i in range(m+2):
        Y2[i,0]=Y2[i,1]
        Y2[i,n+1]=Y2[i,n]
        # vertical sides
    for j in range(1,n+1):
        if j<int(n/2)+1:          
            Y2[0,j]=0
        else:
            Y2[0,j]=Y2I
            
        Y2[m+1,j]=Y2[m,j]
        

def ComputeTransportCoefficients(u,v,aW,aE,aS,aN,aP,Ct):
    
    for i in range(1,m+1):
        for j in range(1,n+1):
              
            FW=dy*u[i-1,j]
            FE=dy*u[i,j]
            FS=dx*v[i,j-1]
            FN=dx*v[i,j]
            
            DW=dy/(Re*Ct*dx)
            DE=dy/(Re*Ct*dx)
            DS=dx/(Re*Ct*dy)
            DN=dx/(Re*Ct*dy)
            
            # compute the flux at the cell faces
            aW[i,j],aE[i,j],aS[i,j],aN[i,j]=FaceFlux(FW,FE,FS,FN,DW,DE,DS,DN,NS)
    
            aP[i,j]=aW[i,j]+aE[i,j]+aS[i,j]+aN[i,j]+FE-FW+FN-FS 
 
        
def SolveTransportEquation(aW,aE,aS,aN,aP,U):
    
    for i in range(1,m+1):
        for j in range(1,n+1):       
            U[i,j]=aW[i,j]*U[i-1,j]+aE[i,j]*U[i+1,j]\
                  +aS[i,j]*U[i,j-1]+aN[i,j]*U[i,j+1]
            U[i,j]/=aP[i,j]
            
            
def Results(T,Y1,Y2,h):
        
    plt.figure(1)
    x=np.linspace(dx/2.,l-dx/2.,m)
    y=np.linspace(dy/2.,h-dy/2.,n)
    Y,X=np.meshgrid(y,x)
    plt.contourf(X,Y,T[1:m+1,1:n+1],50,cmap='jet')
    plt.axes().set_aspect('equal')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('T contours')
    plt.colorbar(orientation='horizontal')

    plt.figure(2)
    plt.contourf(X,Y,Y1[1:m+1,1:n+1],50,cmap='jet')
    plt.axes().set_aspect('equal')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Y1 contours')
    plt.colorbar(orientation='horizontal')

    plt.figure(3)
    plt.contourf(X,Y,Y2[1:m+1,1:n+1],50,cmap='jet')
    plt.axes().set_aspect('equal')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Y2 contours')
    plt.colorbar(orientation='horizontal')

    y=np.linspace(-dy/2.,h+dy/2.,n+2)
    y[0]=0.
    y[n+1]=h
   
    plt.figure(4)
    T[0,:]=100.
    T[m+1,:]=T[m,:]
    T[:,0]=TB
    T[:,n+1]=TT
    for i in range(9):
        plt.plot(T[i*int(m/8),:],y)
        plt.xlabel('T')
        plt.ylabel('y')
        plt.title('T profile')   
    plt.legend(['x=0','x=0.125L','x=0.25L','x=0.375L','x=0.5L'\
               ,'x=0.625L','x=0.75L','x=0.875L','x=L'],frameon=False)
    plt.show()
   
    
    
#================================ main program ================================
m=100                         # number of control volumes in the horizontal direction    
n=40                          # number of control volumes in the vertical direction

l=5.                          # lenght of the channel
h=1.                          # height of the channel
Re=50.                        # Reynolds number
Pr=1.                         # Prandtl number
Sc=0.47                       # Schmidt number

TI=100.                       # inlet temperature
TB=0.                         # bottom side temperature
TT=0.                         # top side temperature
Y1I=1                         # inlet value of species1
Y2I=1                         # inlet value of species2

alpha=0.8                     # pressure correction under relaxation factor
tol=1e-5                      # convergence criterion
iter_max=1e5
dx=l/m
dy=h/n

# matrix definition
u =np.zeros((m+1,n+2));
v =np.zeros((m+2,n+1));
T =np.zeros((m+2,n+2));
T_old=np.zeros((m+2,n+2))
Y1=np.zeros((m+2,n+2));
Y1_old=np.zeros((m+2,n+2))
Y2=np.zeros((m+2,n+2));
Y2_old=np.zeros((m+2,n+2))

aW=np.zeros((m+2,n+2))
aE=np.zeros((m+2,n+2))
aS=np.zeros((m+2,n+2))
aN=np.zeros((m+2,n+2))
aP=np.zeros((m+2,n+2))

# determine type of the scheme
print("\nChoose a scheme:")
print("  1-Centeral")
print("  2-Upwind")
print("  3-Hybrid")
NS=input("  ")
NS=int(NS)

# read the velocity data
# u
f=open("u.txt","r")
Line=f.readlines()
for i in range(m+1):
    for j in range(n+2):
            u[i,j]=float(Line[j+(n+2)*i])
f.close
 
# v
f=open("v.txt","r")
Line=f.readlines()
for i in range(m+2):
    for j in range(n+1):
            v[i,j]=float(Line[j+(n+1)*i])
f.close

# initialize the solution
Initialize(T,Y1,Y2)


#=============================== energy ================================
print('\n solve energy equation:')
tic=time.time()

# compute energy equation coefficients          
ComputeTransportCoefficients(u,v,aW,aE,aS,aN,aP,Pr)
  
it=0
L2T=1.
while L2T>tol and it<iter_max:
    it+=1
        
    # store the solution from the previous iteration
    T_old[:,:]=T[:,:]

    # update ghost cell values
    UpdateGhostCell_Temperature(T)

    # solve energy equation
    SolveTransportEquation(aW,aE,aS,aN,aP,T)
            
    # compute L2 norm of the error
    L2T=ComputeNorm(T,T_old,L2T)  
    
    # print iteration and L2 norm of errors
    print("{:10d}{:10.2e}".format(it,L2T))
           
    
#============================== Species1 ===============================
print('\n solve transport equation for species1:')

# compute species1 equation coefficients          
ComputeTransportCoefficients(u,v,aW,aE,aS,aN,aP,Sc)
  
it=0
L2Y1=1.
while L2Y1>tol and it<iter_max:
    it+=1
        
    # store the solution from the previous iteration
    Y1_old[:,:]=Y1[:,:]

    # update ghost cell values
    UpdateGhostCell_Species1(Y1)

    # solve energy equation
    SolveTransportEquation(aW,aE,aS,aN,aP,Y1)
            
    # compute L2 norm of the error
    L2Y1=ComputeNorm(Y1,Y1_old,L2Y1)  
    
    # print iteration and L2 norm of errors
    print("{:10d}{:10.2e}".format(it,L2Y1))
            
    
#============================== Species2 ===============================
print('\n solve transport equation for species2:')

# compute species2 coefficients          
ComputeTransportCoefficients(u,v,aW,aE,aS,aN,aP,Sc)
  
it=0
L2Y2=1.
while L2Y2>tol and it<iter_max:
    it+=1
        
    # store the solution from the previous iteration
    Y2_old[:,:]=Y2[:,:]

    # update ghost cell values
    UpdateGhostCell_Species2(Y2)

    # solve energy equation
    SolveTransportEquation(aW,aE,aS,aN,aP,Y2)
            
    # compute L2 norm of the error
    L2Y2=ComputeNorm(Y2,Y2_old,L2Y2)  
    
    # print iteration and L2 norm of errors
    print("{:10d}{:10.2e}".format(it,L2Y2))
    
    
toc=time.time()
#=============================== results ===============================
print('\n=================================================================\n')
print(' calculation completed\n')
print(' execution time:{:9.3f}s'.format(toc-tic))

Results(T,Y1,Y2,h)
