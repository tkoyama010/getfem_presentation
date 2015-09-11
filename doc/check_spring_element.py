# coding: utf-8
import getfem as gf
import numpy as np
m = gf.Mesh('cartesian', np.arange(0.0, 10.0, 1.0))
print m
mfu = gf.MeshFem(m,1)
mfu.set_fem(gf.Fem('FEM_PK(1,1)'))
mim = gf.MeshIm(m, gf.Integ('IM_GAUSS1D(1)'))
mfd = gf.MeshFem(m, 1)
mfd.set_fem(gf.Fem('FEM_PK(1,1)'))
nbd = mfd.nbdof()
Lambda = 1.0; Mu = 1.0
K = gf.asm_linear_elasticity(mim, mfu, mfd, np.repeat([Lambda], nbd), np.repeat([Mu], nbd))
K.full()
print 'Lambda = ', Lambda
print 'Mu     = ', Mu
print K.full()
print '2*Mu+Lambda = ', 2*Mu+Lambda
