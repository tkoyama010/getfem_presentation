import getfem as gf
import numpy as np
import scipy.io as io
import scipy.sparse.linalg as linalg
from scipy import io
from scipy.sparse import linalg

print "parameters"

E =  1.320e+11
Nu =  3.430e-01
rho =  8.960e+03

Lambda = E*Nu/((1+Nu)*(1-2*Nu))
Mu = E/(2*(1+Nu))

print "create a Mesh object"
d =  1.000e-03
x =  1.000e-02
y =  1.000e-01
z =  2.000e-02
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

(H,R) = gf.asm_dirichlet(DIRICHLET_BOUNDARY, mim, mfu, mfd, mfd.eval('[[1,0,0],[0,1,0],[0,0,1]]',globals(),locals()), mfd.eval('[0,0,0]'))
(N,U0) = H.dirichlet_nullspace(R)

Nt = gf.Spmat('copy',N)
Nt.transpose()
KK = Nt*K*N
MM = Nt*M*N

KK.save('hb', 'KK.hb')
MM.save('hb', 'MM.hb')

print "solve"

A = io.hb_read('MM.hb')
B = io.hb_read('KK.hb')

w, v = linalg.eigs(A,M=B)
omega = np.sqrt(1.000e+00/w)/2.000e+00/np.pi

print "post-processing"
np.savetxt("omega.txt", np.real(omega))
np.savetxt("v.txt", np.real(v))

sl = gf.Slice(('boundary',), mfu, 3)
sl.export_to_vtk('m.vtk', mfu, N*np.real(v[:,0]), 'Mode1', mfu, N*np.real(v[:,1]), 'Mode2', mfu, N*np.real(v[:,2]), 'Mode3', mfu, N*np.real(v[:,3]), 'Mode4', mfu, N*np.real(v[:,4]), 'Mode5', mfu, N*np.real(v[:,5]), 'Mode6')

