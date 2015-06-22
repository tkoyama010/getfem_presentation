import getfem as gf
import numpy as np
import scipy.io as io
import scipy.sparse.linalg as linalg
from scipy import io
from scipy.sparse import linalg
from scipy import sparse

print "parameters"

rho = 1.000e+00
Lambda = 1.000e+00
Mu = 1.000e+00

print "create a Mesh object"
d =  1.000e+00
x =  1.000e+00
y =  1.000e+00
z =  1.000e+00
m = gf.Mesh('cartesian Q1', np.arange(0., x+d, d), np.arange(0., y+d, d), np.arange(0., z+d, d))
m.set('optimize_structure')

print "create a MeshFem object"
mfu = gf.MeshFem(m,3) # displacement
print "assign the FEM"
mfu.set_fem(gf.Fem('FEM_QK(3,1)'))

print "build a MeshIm object"

mim = gf.MeshIm(m, gf.Integ('IM_HEXAHEDRON(5)'))

print "detect some boundary of the mesh"
P = m.pts()
ctop = (abs(P[0,:] - 0.) < 1e-6)
cbot = (abs(P[1,:] - 0.) < 1e-6)
pidtop = np.compress(ctop,range(0,m.nbpts()))
pidbot = np.compress(cbot,range(0,m.nbpts()))
ftop = m.faces_from_pid(pidtop)
fbot = m.faces_from_pid(pidbot)
print "create boundary region"
NEUMANN_BOUNDARY = 1
DIRICHLET_BOUNDARY = 2
m.set_region(NEUMANN_BOUNDARY,ftop)
m.set_region(DIRICHLET_BOUNDARY,fbot)

print "create a MeshFem object"

mfd = gf.MeshFem(m, 1) # data

print "assign the FEM"
mfd.set_fem(gf.Fem('FEM_QK(3,3)'))

print "assembly"
nbd = mfd.nbdof()
K = gf.asm_linear_elasticity(mim, mfu, mfd, np.repeat([Lambda], nbd), np.repeat([Mu], nbd))
M = gf.asm_mass_matrix(mim, mfu)*rho

K.save('hb', 'K.hb')
M.save('hb', 'M.hb')

print "solve"

A = io.hb_read('M.hb')
B = io.hb_read('K.hb')

K = B.todense()
np.savetxt("K.txt", K, fmt=' (%+.18e) ')
