# -*- coding: utf-8 -*-
"""
Created on Sun Jan 12 17:46:20 2020

@author: Yoon
"""

import sympy as sp

Rgt = sp.Symbol('R_gt', positive=True)
N = sp.symbols('N_Ca', function=True, real=True)
Ld, Lu, H0, cos1, cos2 =sp.symbols('L_d L_u H_0 cos(theta_s) cos(theta_d)' , Constant=True )
# eqn = Eq(6*(Ld/H0*(1-2/Rgt) +Lu/H0)*N +1.34*Rgt*N-cos1-cos2,0)
eqn_sol = sp.solve(6*(Ld/H0*(1-2/Rgt) +Lu/H0)*N +1.34*Rgt*N-cos1-cos2, Rgt)
# print(eqn_sol[0])
sp.init_printing()
eqn_sol

#%% sol1 sol2 png file
import sympy as sp

Rgt = sp.Symbol('R_gt', positive=True)
N = sp.symbols('N_Ca', function=True, positive=True)
Ld, Lu, H0, thetas, thetad =sp.symbols('L_d L_u H_0 theta_s theta_d' , Constant=True )
eqn_sol = sp.solve(6*(Ld/H0*(1-2/Rgt) +Lu/H0)*N +1.34*Rgt*N-sp.cos(thetas)-sp.cos(thetad), Rgt)
sp.init_printing()

#%%
import numpy as np
import matplotlib.pyplot as plt

k=0.02
x1=np.linspace(0,1,100)
x2=np.linspace(1.1, 10,100)
x=np.concatenate((x1, x2), axis=None)
x=np.array(x)
y=k+k/(x-1)
y2=k*np.ones(len(x))
plt.plot(x,y,x,y2,'--')
plt.xlim(0,10)
plt.xlabel('$R_{gt} $ ')
plt.ylabel('$N_{Ca} $ ')
plt.ylim(-2,2)
string='Weeping with $k=$ {}'.format(k)
plt.title(string)

#%%검산
import sympy as sp
sp.init_printing()
a,b =sp.symbols('a b', Constant=True)
expr=(a+b+(b**2)/a)**3
expr2=expr.expand()

#%%
import sympy as sp
sp.init_printing()
a,b =sp.symbols('alpha beta', Constant=True)
expr=(a+b+(b**2)/a)**2
expr.expand()

#%%
import sympy as sp
sp.init_printing()
x=sp.symbols('x')
a,b,c= sp.symbols('a b c', Constant=True)
expr=sp.Eq(a*x**3+b*x**2-c,0)
sol=sp.solve(expr,x)
sol[2]
expr1=a*sol[2]**3
expr2=b*sol[2]**2
expr3=expr1+expr2-c
expr3.simplify()

#%%
import numpy as np
import matplotlib.pyplot as plt

k=(1.6/6)*200/950
x1=np.linspace(0,2,100)
x2=np.linspace(2.001, 10,100)
x=np.concatenate((x1, x2), axis=None)
x=np.array(x)
y=k+2*k/(x-2)

t=(1.6/6)*100/950
x31=np.linspace(0,1,100)
x32=np.linspace(1.001,10,100)
x3=np.concatenate((x31, x32), axis=None)
x3=np.array(x3)
y3=t+t/(x3-1)

y2=k*np.ones(len(x))
plt.plot(x,y,'b',x3,y3,'r')
plt.legend(('Bead Break-up', 'Weeping'))
plt.xlim(0,10)
plt.xlabel('$R_{gt} $ ')
plt.ylabel('$N_{Ca} $ ')
plt.ylim(-2,2)
string='Operability Window when cosine term = 1.6'
plt.title(string)
plt.grid()

#%% c=1.6
import numpy as np
import matplotlib.pyplot as plt
r=np.array([2.79, 2.6, 2.5, 2.4, 2.3, 2.2, 2.1])
c=1.6
a=6.0*(950.0/279.0)*(1.0-2.0/r)
b=1.34*r
alpha=(b**2.0)*(c*a**5.0)**(-1.0/3.0)/9.0
beta=-b/3.0/a
thirdterm=(beta**2.0)/alpha
exact=(alpha+beta+thirdterm)**3.0
model=c/a
dev=np.stack((exact,model))


#%% c=1.0
import numpy as np
r=np.array([2.79, 2.6, 2.5, 2.4, 2.3, 2.2, 2.1])
c=1.0
a=6.0*(950.0/279.0)*(1.0-2.0/r)
b=1.34*r
alpha=(b**2.0)*(c*a**5.0)**(-1.0/3.0)/9.0
beta=-b/3.0/a
thirdterm=(beta**2.0)/alpha
exact=(alpha+beta+thirdterm)**3.0
model=c/a
dev=np.stack((exact,model))
dev

#%%
import numpy as np
import matplotlib.pyplot as plt
r=np.array([4,5,6,7,8,9,10])
c=1.0
a=6.0*(950.0/279.0)*(1.0-2.0/r)
b=1.34*r
alpha=(b**2.0)*(c*a**5.0)**(-1.0/3.0)/9.0
beta=-b/3.0/a
thirdterm=(beta**2.0)/alpha
exact=(alpha+beta+thirdterm)**3.0
model=c/a
plt.plot(r,exact,'r',r,model,'b')
#%% c=1.6 Bead Break-up1
import numpy as np
import matplotlib.pyplot as plt
r=np.array([3.0, 2.9, 2.79, 2.6, 2.5, 2.4, 2.3, 2.2, 2.1, 2.08])
c=1.6
a=6.0*(950.0/279.0)*(1.0-2.0/r)
b=1.34*r
alpha=(b**2.0)*(c*a**5.0)**(-1.0/3.0)/9.0
beta=-b/3.0/a
thirdterm=(beta**2.0)/alpha
exact=(alpha+beta+thirdterm)**3.0
model=c/a
dev=exact-model
plt.title(r'$ \cos{\theta_{d}}+\cos{\theta_{s}} = 1.6 $')
plt.plot(r,model, label='My model')
plt.plot(r,exact, 'r' ,linestyle='--', label='Exact Result')
plt.plot(r,dev, 'bo' ,label='Deviation')
plt.xlabel('$R_{gt} $ ')
plt.ylabel('$N_{Ca} $ ')
plt.legend()
plt.grid()
plt.show()

#%% c=1.0 Bead Break-up1
import numpy as np
import matplotlib.pyplot as plt
r=np.array([3.0, 2.9, 2.79, 2.6, 2.5, 2.4, 2.3, 2.2, 2.1, 2.08])
c=1.0
a=6.0*(950.0/279.0)*(1.0-2.0/r)
b=1.34*r
alpha=(b**2.0)*(c*a**5.0)**(-1.0/3.0)/9.0
beta=-b/3.0/a
thirdterm=(beta**2.0)/alpha
exact=(alpha+beta+thirdterm)**3.0
model=c/a
dev=exact-model
plt.title(r'$ \cos{\theta_{d}}+\cos{\theta_{s}} = 1.0 $')
plt.plot(r,model, label='My model')
plt.plot(r,exact, 'r' ,linestyle='--', label='Exact Result')
plt.plot(r,dev, 'bo' ,label='Deviation')
plt.xlabel('$R_{gt} $ ')
plt.ylabel('$N_{Ca} $ ')
plt.legend()
plt.grid()
plt.show()

#%% c=0.5 Bead Break-up1
import numpy as np
import matplotlib.pyplot as plt
r=np.array([2.79, 2.6, 2.5, 2.4, 2.3, 2.2, 2.1, 2.08])
c=0.5
a=6.0*(950.0/279.0)*(1.0-2.0/r)
b=1.34*r
alpha=(b**2.0)*(c*a**5.0)**(-1.0/3.0)/9.0
beta=-b/3.0/a
thirdterm=(beta**2.0)/alpha
exact=(alpha+beta+thirdterm)**3.0
model=c/a
dev=exact-model
plt.title(r'$ \cos{\theta_{d}}+\cos{\theta_{s}} = 0.5 $')
plt.plot(r,model, label='My model')
plt.plot(r,exact, 'r' ,linestyle='--', label='Exact Result')
plt.plot(r,dev, 'bo' ,label='Deviation')
plt.xlabel('$R_{gt} $ ')
plt.ylabel('$N_{Ca} $ ')
plt.legend()
plt.grid()
plt.show()




#%% Bead Break-up2
import numpy as np
import matplotlib.pyplot as plt
r=np.linspace(0.001, 1.5 , 20)

c=1.6
a=6.0*(950.0/279.0)*(1.0-2.0/r)
b=1.34*r
alpha = -(-c/a)**(1/3)
beta=-b/3.0/a
thirdterm=(beta**2.0)/alpha
exact=(alpha+beta+thirdterm)**3.0
model=c/a
dev=exact-model

fig=plt.figure(figsize=(11,5))
## fig.suptitle('Bead Break-up when $ R_{gt} < 2 $' , fontsize=18)

ax1=fig.add_subplot(1,2,1)
ax1.set_title(r'$ \cos{\theta_{d}}+\cos{\theta_{s}} = 1.6 $')
ax1.plot(r,model, label='My model')
ax1.plot(r,exact, 'r' ,linestyle='--', label='Exact Result')
ax1.plot(r,dev, 'bo' ,label='Deviation')
ax1.set_xlabel('$R_{gt} $ ' , fontsize=18)
ax1.set_ylabel('$N_{Ca} $ ' , fontsize=18)
plt.legend()
plt.grid()

c=1.0
a=6.0*(950.0/279.0)*(1.0-2.0/r)
b=1.34*r
alpha = -(-c/a)**(1/3)
beta=-b/3.0/a
thirdterm=(beta**2.0)/alpha
exact=(alpha+beta+thirdterm)**3.0
model=c/a
dev=exact-model

ax1=fig.add_subplot(1,2,2)
ax1.set_title(r'$ \cos{\theta_{d}}+\cos{\theta_{s}} = 1.0 $')
ax1.plot(r,model, label='My model')
ax1.plot(r,exact, 'r' ,linestyle='--', label='Exact Result')
ax1.plot(r,dev, 'bo' ,label='Deviation')
ax1.set_xlabel('$R_{gt} $ ' , fontsize=18)
ax1.set_ylabel('$N_{Ca} $ ' , fontsize=18)
plt.legend()
plt.grid()

fig.tight_layout()
plt.show()

#%% Weeping1
import numpy as np
import matplotlib.pyplot as plt
r=np.linspace(1.1, 5.0 , 20)

c=1.6
a=6.0*(950.0/279.0)*(2.0-2.0/r)         #Weeping
b=1.34*r
alpha = (b**2.0)*(c*a**5.0)**(-1.0/3.0)/9.0
beta=-b/3.0/a
thirdterm=(beta**2.0)/alpha
exact=(alpha+beta+thirdterm)**3.0
model=c/a
dev=exact-model

fig=plt.figure(figsize=(11,5))
## fig.suptitle('Bead Break-up when $ R_{gt} < 2 $' , fontsize=18)

ax1=fig.add_subplot(1,2,1)
ax1.set_title(r'$ \cos{\theta_{d}}+\cos{\theta_{s}} = 1.6 $')
ax1.plot(r,model, label='My model')
ax1.plot(r,exact, 'r' ,linestyle='--', label='Exact Result')
ax1.plot(r,dev, 'bo' ,label='Deviation')
ax1.set_xlabel('$R_{gt} $ ' , fontsize=18)
ax1.set_ylabel('$N_{Ca} $ ' , fontsize=18)
plt.legend()
plt.grid()

c=1.0
a=6.0*(950.0/279.0)*(2.0-2.0/r)          #Weeping
b=1.34*r
alpha = (b**2.0)*(c*a**5.0)**(-1.0/3.0)/9.0
beta=-b/3.0/a
thirdterm=(beta**2.0)/alpha
exact=(alpha+beta+thirdterm)**3.0
model=c/a
dev=exact-model

ax1=fig.add_subplot(1,2,2)
ax1.set_title(r'$ \cos{\theta_{d}}+\cos{\theta_{s}} = 1.0 $')
ax1.plot(r,model, label='My model')
ax1.plot(r,exact, 'r' ,linestyle='--', label='Exact Result')
ax1.plot(r,dev, 'bo' ,label='Deviation')
ax1.set_xlabel('$R_{gt} $ ' , fontsize=18)
ax1.set_ylabel('$N_{Ca} $ ' , fontsize=18)
plt.legend()
plt.grid()

fig.tight_layout()
plt.show()

#%% Weeping2
import numpy as np
import matplotlib.pyplot as plt
r=np.linspace(0.001, 1.0 , 20)

c=1.6
a=6.0*(950.0/279.0)*(2.0-2.0/r)         #Weeping
b=1.34*r
alpha = -(-c/a)**(1/3)
beta=-b/3.0/a
thirdterm=(beta**2.0)/alpha
exact=(alpha+beta+thirdterm)**3.0
model=c/a
dev=exact-model

fig=plt.figure(figsize=(11,5))
## fig.suptitle('Bead Break-up when $ R_{gt} < 2 $' , fontsize=18)

ax1=fig.add_subplot(1,2,1)
ax1.set_title(r'$ \cos{\theta_{d}}+\cos{\theta_{s}} = 1.6 $')
ax1.plot(r,model, label='My model')
ax1.plot(r,exact, 'r' ,linestyle='--', label='Exact Result')
ax1.plot(r,dev, 'bo' ,label='Deviation')
ax1.set_xlabel('$R_{gt} $ ' , fontsize=18)
ax1.set_ylabel('$N_{Ca} $ ' , fontsize=18)
plt.legend()
plt.grid()

c=1.0
a=6.0*(950.0/279.0)*(2.0-2.0/r)          #Weeping
b=1.34*r
alpha = -(-c/a)**(1/3)
beta=-b/3.0/a
thirdterm=(beta**2.0)/alpha
exact=(alpha+beta+thirdterm)**3.0
model=c/a
dev=exact-model

ax1=fig.add_subplot(1,2,2)
ax1.set_title(r'$ \cos{\theta_{d}}+\cos{\theta_{s}} = 1.0 $')
ax1.plot(r,model, label='My model')
ax1.plot(r,exact, 'r' ,linestyle='--', label='Exact Result')
ax1.plot(r,dev, 'bo' ,label='Deviation')
ax1.set_xlabel('$R_{gt} $ ' , fontsize=18)
ax1.set_ylabel('$N_{Ca} $ ' , fontsize=18)
plt.legend()
plt.grid()

fig.tight_layout()
plt.show()



#%% Weeping
import numpy as np
import matplotlib.pyplot as plt

r1=np.linspace(1.1, 5.0 , 20)

c=1.6


a=6.0*(950.0/279.0)*(2.0-2.0/r1)         #Weeping
b=1.34*r1
alpha1 = (b**2.0)*(c*a**5.0)**(-1.0/3.0)/9.0  #Weeping1
beta1=-b/3.0/a
thirdterm1=(beta1**2.0)/alpha1
exact1=(alpha1+beta1+thirdterm1)**3.0
model1=c/a
dev1=exact1-model1

r2=np.linspace(0.001, 1.0 , 20)
a=6.0*(950.0/279.0)*(2.0-2.0/r2)         #Weeping
b=1.34*r2
alpha2 = -(-c/a)**(1/3)                  #Weeping2
beta2=-b/3.0/a
thirdterm2=(beta2**2.0)/alpha2
exact2=(alpha2+beta2+thirdterm2)**3.0
model2=c/a
dev2=exact2-model2

r=np.concatenate((r2, r1 ), axis=None)
r=np.array(r)

model=np.concatenate((model2, model1 ), axis=None)
model=np.array(model)

exact=np.concatenate((exact2, exact1 ), axis=None)
exact=np.array(exact)

dev=np.concatenate((dev2, dev1), axis=None)
dev=np.array(dev)

fig=plt.figure(figsize=(11,5))
## fig.suptitle('Bead Break-up when $ R_{gt} < 2 $' , fontsize=18)

ax1=fig.add_subplot(1,2,1)
ax1.set_title(r'$ \cos{\theta_{d}}+\cos{\theta_{s}} = 1.6 $')
ax1.plot(r,model, label='My model')
ax1.plot(r,exact, 'r' ,linestyle='--', label='Exact Result')
ax1.plot(r,dev, 'bo' ,label='Deviation')
ax1.set_xlabel('$R_{gt} $ ' , fontsize=18)
ax1.set_ylabel('$N_{Ca} $ ' , fontsize=18)
plt.legend()
plt.grid()


c=1.0

a=6.0*(950.0/279.0)*(2.0-2.0/r1)         #Weeping
b=1.34*r1
alpha1 = (b**2.0)*(c*a**5.0)**(-1.0/3.0)/9.0  #Weeping1
beta1=-b/3.0/a
thirdterm1=(beta1**2.0)/alpha1
exact1=(alpha1+beta1+thirdterm1)**3.0
model1=c/a
dev1=exact1-model1

r2=np.linspace(0.001, 1.0 , 20)
a=6.0*(950.0/279.0)*(2.0-2.0/r2)         #Weeping
b=1.34*r2
alpha2 = -(-c/a)**(1/3)                  #Weeping2
beta2=-b/3.0/a
thirdterm2=(beta2**2.0)/alpha2
exact2=(alpha2+beta2+thirdterm2)**3.0
model2=c/a
dev2=exact2-model2

r=np.concatenate((r2, r1 ), axis=None)
r=np.array(r)

model=np.concatenate((model2, model1 ), axis=None)
model=np.array(model)

exact=np.concatenate((exact2, exact1 ), axis=None)
exact=np.array(exact)

dev=np.concatenate((dev2, dev1), axis=None)
dev=np.array(dev)

ax1=fig.add_subplot(1,2,2)
ax1.set_title(r'$ \cos{\theta_{d}}+\cos{\theta_{s}} = 1.0 $')
ax1.plot(r,model, label='My model')
ax1.plot(r,exact, 'r' ,linestyle='--', label='Exact Result')
ax1.plot(r,dev, 'bo' ,label='Deviation')
ax1.set_xlabel('$R_{gt} $ ' , fontsize=18)
ax1.set_ylabel('$N_{Ca} $ ' , fontsize=18)
plt.legend()
plt.grid()

fig.tight_layout()
plt.show()


#%% Bead Break-up
import numpy as np
import matplotlib.pyplot as plt

def jiconcen(r1,r2):
    r=np.concatenate((r2, r1 ), axis=None)
    result=np.array(r)
    return result

def jibreakup1(r,c):
    a=6.0*(950.0/279.0)*(1.0-2.0/r)
    b=1.34*r
    alpha=(b**2.0)*(c*a**5.0)**(-1.0/3.0)/9.0
    beta=-b/3.0/a
    thirdterm=(beta**2.0)/alpha
    exact=(alpha+beta+thirdterm)**3.0
    model=c/a
    dev=exact-model
    return exact, model, dev

def jibreakup2(r,c):
    a=6.0*(950.0/279.0)*(1.0-2.0/r)
    b=1.34*r
    alpha = -(-c/a)**(1/3)
    beta=-b/3.0/a
    thirdterm=(beta**2.0)/alpha
    exact=(alpha+beta+thirdterm)**3.0
    model=c/a
    dev=exact-model
    return exact, model, dev


r1=np.linspace(2.1, 5.0 , 20)
result1=jibreakup1(r1, 1.6)
r2=np.linspace(0.001, 2.0 , 20)
result2=jibreakup2(r2, 1.6)

r=jiconcen(r1,r2)
exact=jiconcen(result1[0],result2[0])
model=jiconcen(result1[1],result2[1])
dev=jiconcen(result1[2],result2[2])


fig=plt.figure(figsize=(11,5))

ax1=fig.add_subplot(1,2,1)
ax1.set_title(r'$ \cos{\theta_{d}}+\cos{\theta_{s}} = 1.6 $')
ax1.plot(r,model, label='My model')
ax1.plot(r,exact, 'r' ,linestyle='--', label='Exact Result')
ax1.plot(r,dev, 'bo' ,label='Deviation')
ax1.set_xlabel('$R_{gt} $ ' , fontsize=18)
ax1.set_ylabel('$N_{Ca} $ ' , fontsize=18)
plt.legend()
plt.grid()


result1=jibreakup1(r1, 1.0)
result2=jibreakup2(r2, 1.0)

r=jiconcen(r1,r2)
exact=jiconcen(result1[0],result2[0])
model=jiconcen(result1[1],result2[1])
dev=jiconcen(result1[2],result2[2])

ax1=fig.add_subplot(1,2,2)
ax1.set_title(r'$ \cos{\theta_{d}}+\cos{\theta_{s}} = 1.0 $')
ax1.plot(r,model, label='My model')
ax1.plot(r,exact, 'r' ,linestyle='--', label='Exact Result')
ax1.plot(r,dev, 'bo' ,label='Deviation')
ax1.set_xlabel('$R_{gt} $ ' , fontsize=18)
ax1.set_ylabel('$N_{Ca} $ ' , fontsize=18)
plt.legend()
plt.grid()

fig.tight_layout()
plt.show()