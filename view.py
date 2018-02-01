import numpy            as np
import matplotlib.pylab as plt


'Settings'
nzones = 2048
nssdw = 6 
alpha = 1.048		# 1.256 for n=6?

Esn = 1.00e+51
Mej = 2.00e+33 * 2.0
d0  = 1.00e-24
age = 3.15e+07 * 50

vt   = np.sqrt( (10.0*nssdw - 50) / (3.0*nssdw - 9.0) * Esn / Mej )
capA = (5.0*nssdw - 25) / (2.0*np.pi*nssdw) * Esn * vt**(-5.0)
Rsh  = alpha * (capA * vt**nssdw / d0)**(1.0/nssdw) * age**((nssdw-3.0)/nssdw)

d_plateau = capA / age**3
r_plateau = vt * age


data1 = np.loadtxt('x060a')
data2 = np.loadtxt('x060b')


dsc1 = d0
dsc2 = d0

rsc1 = Rsh
rsc2 = Rsh

zxc1 = data1[:,0]
zxc2 = data2[:,0]

zro1 = data1[:,1]
zro2 = data2[:,1]

rss1 = np.zeros(zxc1.size)
dss1 = np.zeros(zxc1.size)
pss1 = np.zeros(zxc1.size)
uss1 = np.zeros(zxc1.size)

rss2 = np.zeros(zxc2.size)
dss2 = np.zeros(zxc2.size)
pss2 = np.zeros(zxc2.size)
uss2 = np.zeros(zxc2.size)


i = zxc1.size-1
while zro1[i] < 1.1:
  i -= 1
zxc1 = zxc1 / zxc1[i]

i = zxc2.size-1
while zro2[i] < 1.1:
  i -= 1
zxc2 = zxc2 / zxc2[i]


zux1 = data1[:,3]
zux2 = data2[:,3]

usc1 = (zxc1[0]*Rsh) / age / zux1[0]
usc2 = (zxc2[0]*Rsh) / age / zux2[0]

psc1 = dsc1*usc1**2
psc2 = dsc2*usc2**2

for i in range(zxc1.size):
  rss1[i] = zxc1[i]*rsc1
  dss1[i] = data1[i,1]*dsc1
  pss1[i] = data1[i,2]*psc1
  uss1[i] = data1[i,3]*usc1



for i in range(zxc2.size):
  rss2[i] = zxc2[i]*rsc2
  dss2[i] = data2[i,1]*dsc2
  pss2[i] = data2[i,2]*psc2
  uss2[i] = data2[i,3]*usc2


print rss1[0], dss1[0], pss1[0], uss1[0]
print rss2[0], dss2[0], pss2[0], uss2[0]
print ' '


as1 = dss1[0]*uss1[0]**nssdw
as2 = dss2[0]*uss2[0]**nssdw

print as1, as2


xmin = 0.5*rss2[0]
xmax = 1.1*rss2[-1]

zxc = np.linspace(xmin, xmax, nzones)
zro1 = np.zeros(nzones)
zpr1 = np.zeros(nzones)
zux1 = np.zeros(nzones)

zro2 = np.zeros(nzones)
zpr2 = np.zeros(nzones)
zux2 = np.zeros(nzones)


for i in range(nzones):
  if zxc[i] < rss1[0]:
    zux1[i] = zxc[i] / age
    zro1[i] = min( d_plateau, as1 / (zux1[i]**nssdw) )
    zpr1[i] = pss1[0]
  elif zxc[i] > rss1[-1]: 
    zux1[i] = 0.0
    zro1[i] = d0
    zpr1[i] = pss1[-1]
  else:
    n = 1
    while rss1[n] < zxc[i]:
      n += 1
    nm = n-1

    frac = (zxc[i]-rss1[nm])/(rss1[n]-rss1[nm])
    zro1[i] = frac*dss1[n] + (1.0-frac)*dss1[nm]
    zpr1[i] = frac*pss1[n] + (1.0-frac)*pss1[nm]
    zux1[i] = frac*uss1[n] + (1.0-frac)*uss1[nm]


  if zxc[i] < rss2[0]:
    zux2[i] = zxc[i] / age
    zro2[i] = min( d_plateau, as2 / (zux2[i]**nssdw) )
    zpr2[i] = pss2[0]
  elif zxc[i] > rss2[-1]: 
    zux2[i] = 0.0
    zro2[i] = d0
    zpr2[i] = pss2[-1]
  else:
    n = 1
    while rss2[n] < zxc[i]:
      n += 1
    nm = n-1

    frac = (zxc[i]-rss2[nm])/(rss2[n]-rss2[nm])
    zro2[i] = frac*dss2[n] + (1.0-frac)*dss2[nm]
    zpr2[i] = frac*pss2[n] + (1.0-frac)*pss2[nm]
    zux2[i] = frac*uss2[n] + (1.0-frac)*uss2[nm]


'plot'
'''
dnorm = (np.log(zro)-np.max(np.log(zro)))/(np.max(np.log(zro))-np.min(np.log(zro)))+1
pnorm = zpr/np.max(zpr)
unorm = zux/np.max(zux)

plt.plot(zxc, dnorm, label=r'$\rho$')
plt.plot(zxc, pnorm, label=r'$P$'   )
plt.plot(zxc, unorm, label=r'$u$'   )
'''

plt.plot(np.log(zpr1))
plt.plot(np.log(zpr2), 'k')

#plt.legend(loc=6)
plt.show()
