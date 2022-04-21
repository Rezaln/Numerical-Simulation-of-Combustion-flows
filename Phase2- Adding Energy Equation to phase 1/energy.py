"""
This code was my university project 
Code by: Reza Lotfi Navaei

"""
#===================================
#       Developing Flow in a
#       Rectangular  Channel
#         SMIPLE Algorithm
#===================================
import time
import numpy as np
import matplotlib.pyplot as plt



#============================ function definition =============================
def Initialize(u,v,p,T):
    
    for i in range(m+1):
        for j in range(n+2):
            u[i,j]=uI
 
    for i in range(m+2):
        for j in range(n+1):
            v[i,j]=0
           
    for i in range(m+2):
        for j in range(n+2):
            p[i,j]=0
            T[i,j]=0
     
        
def UpdateGhostCell_Momentum(u,v):
    
    # u
    # horizontal sides
    for i in range(m+1):
        u[i,0]=-u[i,1]
        u[i,n+1]=-u[i,n]
    # vertical sides
    for j in range(1,n+1):
        u[0,j]=uI
        u[m,j]=u[m-1,j]
    
    # v
    # horizontal sides
    for i in range(m+2):
        v[i,0]=0.
        v[i,n]=0.
    # vertical sides
    for j in range(1,n):
        v[0,j]=-v[1,j]
        v[m+1,j]=v[m,j]
 
    
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
           

def ComputeX_MomentumCoefficients(u,v,du,aXW,aXE,aXS,aXN,aXP,NS):
    
    for i in range(1,m):
        for j in range(1,n+1):
            
            FW=dy*0.5*(u[i,j]  +u[i-1,j]  )
            FE=dy*0.5*(u[i,j]  +u[i+1,j]  )
            FS=dx*0.5*(v[i,j-1]+v[i+1,j-1])
            FN=dx*0.5*(v[i,j]  +v[i+1,j]  )
            
            DW=dy/(Re*dx)
            DE=dy/(Re*dx)
            DS=dx/(Re*dy)
            DN=dx/(Re*dy)
            
            # compute the flux at the cell faces
            aXW[i,j],aXE[i,j],aXS[i,j],aXN[i,j]=FaceFlux(FW,FE,FS,FN,DW,DE,DS,DN,NS)
            
            aXP[i,j]=aXW[i,j]+aXE[i,j]+aXS[i,j]+aXN[i,j]+FE-FW+FN-FS
            du[i,j]=dy/aXP[i,j]
        
        
def ComputeY_MomentumCoefficients(u,v,dv,aYW,aYE,aYS,aYN,aYP,NS):
    
    for i in range(1,m+1):
        for j in range(1,n):
    
            FW=dy*0.5*(u[i-1,j]+u[i-1,j+1])
            FE=dy*0.5*(u[i,j]  +u[i,j+1]  )
            FS=dx*0.5*(v[i,j]  +v[i,j-1]  )
            FN=dx*0.5*(v[i,j]  +v[i,j+1]  )
    
            DW=dy/(Re*dx)
            DE=dy/(Re*dx)
            DS=dx/(Re*dy)
            DN=dx/(Re*dy)
            
            # compute the flux at the cell faces
            aYW[i,j],aYE[i,j],aYS[i,j],aYN[i,j]=FaceFlux(FW,FE,FS,FN,DW,DE,DS,DN,NS)

            aYP[i,j]=aYW[i,j]+aYE[i,j]+aYS[i,j]+aYN[i,j]+FE-FW+FN-FS  
            dv[i,j]=dx/aYP[i,j]  


def SolveX_Momentun(aXW,aXE,aXS,aXN,aXP,p,u):
    
    for i in range(1,m):
        for j in range(1,n+1):       
            u[i,j]=aXW[i,j]*u[i-1,j]+aXE[i,j]*u[i+1,j]\
                  +aXS[i,j]*u[i,j-1]+aXN[i,j]*u[i,j+1]\
                  +(p[i,j]-p[i+1,j])*dy
            u[i,j]/=aXP[i,j]
            

def SolveY_Momentun(aYW,aYE,aYS,aYN,aYP,p,v):
    
    for i in range(1,m+1):
        for j in range(1,n):       
            v[i,j]=aYW[i,j]*v[i-1,j]+aYE[i,j]*v[i+1,j]\
                  +aYS[i,j]*v[i,j-1]+aYN[i,j]*v[i,j+1]\
                  +(p[i,j]-p[i,j+1])*dx
            v[i,j]/=aYP[i,j]  
        
        
def ComputePressureCoefficients(du,dv,aW,aE,aS,aN,aP,Source):
    
    for i in range(1,m+1):    
        for j in range(1,n+1): 
            
            aW[i,j]=dy*du[i-1,j]
            aE[i,j]=dy*du[i,j]
            aS[i,j]=dx*dv[i,j-1]
            aN[i,j]=dx*dv[i,j]
            aP[i,j]=aW[i,j]+aE[i,j]+aS[i,j]+aN[i,j]
            
            Source[i,j]=dy*(u[i-1,j]-u[i,j])\
                       +dx*(v[i,j-1]-v[i,j])
                       
 
def SolvePressureCorrection(aW,aE,aS,aN,aP,Source,pc):
  
    pc[:,:]=0.
    for i in range(1,m+1):
        for j in range(1,n+1):       
            pc[i,j]=aW[i,j]*pc[i-1,j]+aE[i,j]*pc[i+1,j]\
                   +aS[i,j]*pc[i,j-1]+aN[i,j]*pc[i,j+1]+Source[i,j]
            pc[i,j]/=aP[i,j]


def UpdateVariables(du,dv,pc,u,v,p):
                    
    for i in range(1,m):
        for j in range(1,n+1):
            u[i,j]=u[i,j]+du[i,j]*(pc[i,j]-pc[i+1,j])
        
    for i in range(1,m+1):
        for j in range(1,n):
            v[i,j]=v[i,j]+dv[i,j]*(pc[i,j]-pc[i,j+1])
            
    for i in range(1,m+1):
        for j in range(1,n+1):
            p[i,j]=p[i,j]+alpha*pc[i,j]
       
        
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
        

def ComputeEnergyCoefficients(u,v,aW,aE,aS,aN,aP):
    
    for i in range(1,m+1):
        for j in range(1,n+1):
              
            FW=dy*u[i-1,j]
            FE=dy*u[i,j]
            FS=dx*v[i,j-1]
            FN=dx*v[i,j]
            
            DW=dy/(Re*Pr*dx)
            DE=dy/(Re*Pr*dx)
            DS=dx/(Re*Pr*dy)
            DN=dx/(Re*Pr*dy)
            
            # compute the flux at the cell faces
            aW[i,j],aE[i,j],aS[i,j],aN[i,j]=FaceFlux(FW,FE,FS,FN,DW,DE,DS,DN,NS)
    
            aP[i,j]=aW[i,j]+aE[i,j]+aS[i,j]+aN[i,j]+FE-FW+FN-FS 
 
        
def SolveEnergy(aW,aE,aS,aN,aP,T):
    
    for i in range(1,m+1):
        for j in range(1,n+1):       
            T[i,j]=aW[i,j]*T[i-1,j]+aE[i,j]*T[i+1,j]\
                  +aS[i,j]*T[i,j-1]+aN[i,j]*T[i,j+1]
            T[i,j]/=aP[i,j]
            
            
def Results(u,v,p,T,h):
        
    plt.figure(1)
    x=np.linspace(0.,l,m+1)
    y=np.linspace(-dy/2.,h+dy/2.,n+2)
    Y,X=np.meshgrid(y,x)
    Y[:,0]=0.
    Y[:,n+1]=h
    u[:,0]=0.
    u[:,n+1]=0.
    plt.contourf(X,Y,u,50,cmap='jet')
    plt.axes().set_aspect('equal')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('u contours')
    plt.colorbar(orientation='horizontal')
    
    plt.figure(2)
    x=np.linspace(-dx/2.,l+dx/2.,m+2)
    y=np.linspace(0.,h,n+1)
    Y,X=np.meshgrid(y,x)
    X[0,:]=0.
    X[m+1,:]=l
    v[0,:]=0.
    v[m+1,:]=v[m,:]
    plt.contourf(X,Y,v,50,cmap='jet')
    plt.axes().set_aspect('equal')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('v contours')
    plt.colorbar(orientation='horizontal')
    
    plt.figure(3)
    x=np.linspace(dx/2.,l-dx/2.,m)
    y=np.linspace(dy/2.,h-dy/2.,n)
    Y,X=np.meshgrid(y,x)
    plt.contourf(X,Y,p[1:m+1,1:n+1],50,cmap='jet')
    plt.axes().set_aspect('equal')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('p contours')
    plt.colorbar(orientation='horizontal')
  
    plt.figure(4)
    x=np.linspace(dx/2.,l-dx/2.,m)
    y=np.linspace(dy/2.,h-dy/2.,n)
    Y,X=np.meshgrid(y,x)
    plt.contourf(X,Y,T[1:m+1,1:n+1],50,cmap='jet')
    plt.axes().set_aspect('equal')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('T contours')
    plt.colorbar(orientation='horizontal')

    y=np.linspace(-dy/2.,h+dy/2.,n+2)
    y[0]=0.
    y[n+1]=h
        
    # calculate exact velocity profile at the outlet
    h=h/2.
    for j in range(n+2):
        u_exact[j]=1.5*uI*(2*h*y[j]-y[j]**2)/h**2
        
    plt.figure(5)
    plt.plot(u_exact,y,u[m,:],y)
    plt.xlabel('u')
    plt.ylabel('y')
    plt.title('u(y)')
    plt.legend(['exact','numerical'])
    
    # plot the error
    plt.figure(6)
    plt.plot(u_exact-u[m,:],y)
    plt.xlabel('error')
    plt.ylabel('y')
    plt.title('error')

    plt.figure(7)
    T[0,:]=100.
    T[m+1,:]=T[m,:]
    T[:,0]=TB
    T[:,n+1]=TT
    for i in range(1,11):
        plt.plot(i*int(m/10)*dx+T[i*int(m/10),:]/TI,y)
        plt.axes().set_aspect('equal')
        plt.xlabel('T')
        plt.ylabel('y')
        plt.ylim(0,h)
        plt.title('T profile') 
    plt.show()
   
    
    
#================================ main program ================================
m=100                         # number of control volumes in the horizontal direction    
n=40                          # number of control volumes in the vertical direction

l=20.                         # lenght of the channel
h=1.                          # height of the channel
Re=50.                        # Reynolds number
Pr=1.                         # Prandtl number
uI=1.                         # inlet velocity
TI=80.                        # inlet temperature
TB=50.                        # bottom side temperature
TT=50.                        # top side temperature
alpha=0.8                     # pressure correction under relaxation factor
tol=1e-4                      # convergence criterion
iter_max=1e5
dx=l/m
dy=h/n

# matrix definition
u =np.zeros((m+1,n+2));  u_old=np.zeros((m+1,n+2))
v =np.zeros((m+2,n+1));  v_old=np.zeros((m+2,n+1))
p =np.zeros((m+2,n+2));  p_old=np.zeros((m+2,n+2))
T =np.zeros((m+2,n+2));  T_old=np.zeros((m+2,n+2))
pc=np.zeros((m+2,n+2))

aXW=np.zeros((m+1,n+2)); aYW=np.zeros((m+2,n+1)); aW=np.zeros((m+2,n+2))
aXE=np.zeros((m+1,n+2)); aYE=np.zeros((m+2,n+1)); aE=np.zeros((m+2,n+2))
aXS=np.zeros((m+1,n+2)); aYS=np.zeros((m+2,n+1)); aS=np.zeros((m+2,n+2))
aXN=np.zeros((m+1,n+2)); aYN=np.zeros((m+2,n+1)); aN=np.zeros((m+2,n+2))
aXP=np.zeros((m+1,n+2)); aYP=np.zeros((m+2,n+1)); aP=np.zeros((m+2,n+2))
du =np.zeros((m+1,n+2)); dv =np.zeros((m+2,n+1)); 
Source=np.zeros((m+2,n+2))
u_exact=np.zeros(n+2)

# determine type of the scheme
print("\nChoose a scheme:")
print("  1-Central")
print("  2-Upwind")
print("  3-Hybrid")
NS=input("  ")
NS=int(NS)

# initialize the solution
Initialize(u,v,p,T)


#============================== main loop ==============================
tic=time.time()

print('\n solve Navier-Stokes equations:')
it=0
L2u=1.
L2v=1.
L2p=1.
while (L2u>tol or L2v>tol or L2p >tol) and it<iter_max:
    it+=1
    
    # store the solution from the previous iteration
    u_old[:,:]=u[:,:]
    v_old[:,:]=v[:,:]
    p_old[:,:]=p[:,:]

    # update ghost cell values
    UpdateGhostCell_Momentum(u,v)
        
    # compute x_momentun coefficients       
    ComputeX_MomentumCoefficients(u,v,du,aXW,aXE,aXS,aXN,aXP,NS)
    
    # compute y_momentun coefficients          
    ComputeY_MomentumCoefficients(u,v,dv,aYW,aYE,aYS,aYN,aYP,NS)
                 
    # solve x_momentum
    SolveX_Momentun(aXW,aXE,aXS,aXN,aXP,p,u)
                          
    # solve y_momentum 
    SolveY_Momentun(aYW,aYE,aYS,aYN,aYP,p,v)
            
    # pressure correction coefficients and source term
    ComputePressureCoefficients(du,dv,aW,aE,aS,aN,aP,Source) 
                       
    # solve pressure correction
    SolvePressureCorrection(aW,aE,aS,aN,aP,Source,pc)  
           
    # update u, v & p
    UpdateVariables(du,dv,pc,u,v,p)
                   
    # compute L2 norm of the error
    L2u=ComputeNorm(u,u_old,L2u)  
    L2v=ComputeNorm(v,v_old,L2v)  
    L2p=ComputeNorm(p,p_old,L2p)  
    
    # print iteration and L2 norm of errors
    print("{:10d}{:10.2e}{:10.2e}{:10.2e}".format(it,L2u,L2v,L2p))
    
    
#=============================== energy ================================
print('\n solve energy equation:')

# compute energy equation coefficients          
ComputeEnergyCoefficients(u,v,aW,aE,aS,aN,aP)
  
it=0
L2T=1.
while L2T>tol and it<iter_max:
    it+=1
        
    # store the solution from the previous iteration
    T_old[:,:]=T[:,:]

    # update ghost cell values
    UpdateGhostCell_Temperature(T)

    # solve energy equation
    SolveEnergy(aW,aE,aS,aN,aP,T)
            
    # compute L2 norm of the error
    L2T=ComputeNorm(T,T_old,L2T)  
    
    # print iteration and L2 norm of errors
    print("{:10d}{:10.2e}".format(it,L2T))
           
    
    
toc=time.time()
#=============================== results ===============================
print('\n=================================================================\n')
print(' calculation completed\n')
print(' execution time:{:9.3f}s'.format(toc-tic))
print(' u at the middle of the oulet:{:5.2f}\n'.format(u[m,int(h/(2*dy))]))

Results(u,v,p,T,h)