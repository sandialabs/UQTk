import numpy as npy
import matplotlib.pyplot as plt
from scipy import stats, mgrid, c_, reshape, random, rot90

npts=500
kdeSmall=1.e-30

def klFunc(x,y,z):
  dx = x[0,1]-x[0,0]
  dy = y[1,0]-y[0,0]
  psum=0.0
  z[npy.where(z<kdeSmall)]=kdeSmall
  for i in range(x.shape[0]-1):
    for j in range(x.shape[1]-1):
      psum = psum+0.25*dx*dy*(z[i,j]+z[i+1,j]+z[i,j+1]+z[i+1,j+1])
  #print k,psum
  z=z/psum
  ze = x.copy()
  for i in range(x.shape[0]):
    for j in range(x.shape[1]):
      ze[i,j] = -(1.0-x[i,j])**2-100.*(y[i,j]-x[i,j]*x[i,j])**2
  psum=0.0
  for i in range(x.shape[0]-1):
    for j in range(x.shape[1]-1):
      psum = psum+0.25*dx*dy*(npy.exp(ze[i,j])+npy.exp(ze[i+1,j])+npy.exp(ze[i,j+1])+npy.exp(ze[i+1,j+1]))
  #print k,psum
  ze=ze-npy.log(psum)
  pSol=x.copy();
  for i in range(x.shape[0]):
    for j in range(x.shape[1]):
      pSol[i,j] = z[i,j]*(npy.log(z[i,j])-ze[i,j])
  KLdiv=0.0
  for i in range(x.shape[0]-1):
    for j in range(x.shape[1]-1):
      KLdiv = KLdiv+0.25*dx*dy*(pSol[i,j]+pSol[i,j+1]+pSol[i+1,j]+pSol[i+1,j+1])
  print KLdiv
  return KLdiv;


d10k=[]
d100k=[]
d1000k=[]
for i in range(1,10+1):
    d10k.append(npy.genfromtxt('dens_10k_'+str(i)+'.dat'))
    d100k.append(npy.genfromtxt('dens_100k_'+str(i)+'.dat'))
    d1000k.append(npy.genfromtxt('dens_1000k_'+str(i)+'.dat'))

KL10k=[]
for k in range(10):
  x=reshape(d10k[k][:,0], (npts,npts))
  y=reshape(d10k[k][:,1], (npts,npts))
  z=reshape(d10k[k][:,2], (npts,npts))
  KL10k.append(klFunc(x,y,z))

KL100k=[]
for k in range(10):
  x=reshape(d100k[k][:,0], (npts,npts))
  y=reshape(d100k[k][:,1], (npts,npts))
  z=reshape(d100k[k][:,2], (npts,npts))
  KL100k.append(klFunc(x,y,z))

KL1000k=[]
for k in range(10):
  x=reshape(d1000k[k][:,0], (npts,npts))
  y=reshape(d1000k[k][:,1], (npts,npts))
  z=reshape(d1000k[k][:,2], (npts,npts))
  KL1000k.append(klFunc(x,y,z))

KLd=npy.array([KL10k,KL100k,KL1000k])

nTests=KLd.shape[0]
KLmn=[KLd[i].mean()  for i in range(nTests)]
print "Convergence rate based on the mean KL:"
for i in range(nTests-1):
  print npy.log(KLmn[i]/KLmn[i+1])/npy.log(10.0)


fig = plt.figure(figsize=(6,6))
ax=fig.gca(position=[0.16,0.10,0.80,0.85])

for i in range(3):
  plt.plot([i+1 for j in range(KLd.shape[1])],KLd[i],"ko",markersize=5)

plt.plot([1,2,3],KLmn,color="red", linewidth=2)
ax.set_xlim([0.8,3.2])
ax.xaxis.set_ticks([1,2,3])
ax.xaxis.set_ticklabels(["10k","100k","1000k"])

ax.set_ylim([0.005,0.15])
ax.set_yscale('log')

ax.set_xlabel(r"$N_{spl}$",fontsize=18,labelpad=10)
ax.set_ylabel("$D_k$",fontsize=18)

plt.savefig("KLdiv_rb.pdf")
plt.show()
