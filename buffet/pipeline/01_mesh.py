import numpy as np, math
N_AIR_HALF=200; N_WAKE=120; N_NORM=120
B=1.10; CHORD=1.0; RFF=25.0; WLEN=20.0; FCH=5e-6
def y0012(x): return 0.6*(0.2969*np.sqrt(x)-0.1260*x-0.3516*x**2+0.2843*x**3-0.1015*x**4)
bs=np.linspace(np.pi,0,N_AIR_HALF); xs=0.5*(1-np.cos(bs)); ys=y0012(xs)
ax=np.concatenate([xs[::-1],xs[1:]]); ay=np.concatenate([ys[::-1],-ys[1:]])
wx=np.linspace(CHORD,CHORD+WLEN,N_WAKE+1)[1:]
ix=np.concatenate([wx[::-1],ax,wx]); iy=np.concatenate([np.zeros(N_WAKE),ay,np.zeros(N_WAKE)])
def ofp(x,y,r=RFF):
    if x<=CHORD+1e-9:
        cx,cy=0.5,0; dx=x-cx; dy=y-cy; a=math.atan2(dy,dx)
        return cx+r*math.cos(a),cy+r*math.sin(a)
    s=1.0 if y>=0 else -1.0
    return x,s*r
ox=np.array([ofp(x,y)[0] for x,y in zip(ix,iy)]); oy=np.array([ofp(x,y)[1] for x,y in zip(ix,iy)])
eta=(B**np.arange(N_NORM)-1)/(B**N_NORM-1)
ef=FCH/RFF; eta=np.maximum(eta,ef*np.arange(N_NORM)/N_NORM)
Nx=ix.size; Ny=N_NORM
X=np.zeros((Nx,Ny)); Y=np.zeros((Nx,Ny))
for i in range(Nx):
    X[i,:]=ix[i]+eta*(ox[i]-ix[i]); Y[i,:]=iy[i]+eta*(oy[i]-iy[i])
def w(f,*a): f.write(' '.join(map(str,a))+'\n')
with open('mesh.su2','w') as f:
    w(f,'NDIME=',2); w(f,'NPOIN=',Nx*Ny)
    for j in range(Ny):
        for i in range(Nx): w(f,X[i,j],Y[i,j],i+j*Nx)
    ne=(Nx-1)*(Ny-1); w(f,'NELEM=',ne)
    for j in range(Ny-1):
        for i in range(Nx-1):
            n1=i+j*Nx; n2=(i+1)+j*Nx; n3=(i+1)+(j+1)*Nx; n4=i+(j+1)*Nx
            w(f,9,n1,n2,n3,n4)
    w(f,'NMARK=',2)
    w(f,'MARKER_TAG=','airfoil'); w(f,'MARKER_ELEMS=',Nx-1)
    for i in range(Nx-1): w(f,3,i,i+1)
    w(f,'MARKER_TAG=','farfield'); w(f,'MARKER_ELEMS=',Nx-1)
    j=Ny-1
    for i in range(Nx-1,0,-1): w(f,3,i+j*Nx,(i-1)+j*Nx)
print('wrote mesh.su2 npts',Nx*Ny,'nelem',ne)
