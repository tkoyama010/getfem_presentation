
# coding: utf-8

# #はじめに
# 静的な場合の有限要素法問題における変位は以下のような連立方程式をとくことにより得られます。
# 
# $KU=F$
# 
# ただし、$K$は剛性行列、$U$は変位ベクトル、$F$は外力ベクトルです。この資料は三脚のメッシュをpythonでインポートし、境界条件を設定し剛性行列$K$と外力ベクトル$F$を作成して変位$U$を計算する過程を説明したものです。

# #モジュールのインポート
# まずはGetfem++のpythonインターフェースのモジュールをインポートします。さらに、境界条件設定に使うためnumpyをインポートします。

# In[1]:

import getfem as gf
import numpy as np


# #パラメータの設定
# まずはパラメータを設定します。この例では、このファイルからの相対パスが '../mesh/tripod.mesh'であるメッシュファイルを使用します。メッシュは四面体要素です。ヤング係数$E = 1.0\times10^3$、ポアソン比$\nu = 0.3$とします。

# In[2]:

file_msh = '../mesh/tripod.mesh'
E = 1e3
Nu = 0.3
Lambda = E*Nu/((1+Nu)*(1-2*Nu))
Mu = E/(2*(1+Nu))


# #メッシュの読み込み
# メッシュを外部ファイルから読み込みます。メッシュの入力データはGetfem++の独自のフォーマットになっていますが、解読は容易だと思います。

# In[3]:

m = gf.Mesh('load',file_msh)
m.set('optimize_structure')


# 次のコマンドで、mに設定したメッシュをParaviewのスクリプトでpng画像に打ち出し確認します。

# In[4]:

m.export_to_vtk('m.vtk')


# 以下のようにpythonファイルに書き込み実行してpng画像に打ち出します。複雑に見えますが、Paraviewのマクロの記録で作成したものですので、作成は難しくありません。

# In[5]:

get_ipython().run_cell_magic(u'writefile', u'plot.py', u'try: paraview.simple\nexcept: from paraview.simple import *\nparaview.simple._DisableFirstRenderCameraReset()\n\nm_vtk = LegacyVTKReader( FileNames=[\'m.vtk\'] )\n\nRenderView2 = CreateRenderView()\nRenderView2.CompressorConfig = \'vtkSquirtCompressor 0 3\'\nRenderView2.UseLight = 1\nRenderView2.LightSwitch = 0\nRenderView2.Background = [0.31999694819562063, 0.3400015259021897, 0.4299992370489052]\nRenderView2.CenterOfRotation = [9.392949104309082, 1.5, 0.0]\n\nAnimationScene1 = GetAnimationScene()\nAnimationScene1.ViewModules = RenderView2\n\nDataRepresentation2 = Show()\nDataRepresentation2.ScaleFactor = 8.50093994140625\nDataRepresentation2.ScalarOpacityUnitDistance = 8.145737909998818\nDataRepresentation2.EdgeColor = [0.0, 0.0, 0.5000076295109483]\n\nRenderView2.CameraPosition = [9.392949104309082, 1.5, 183.2819944361333]\nRenderView2.OrientationAxesVisibility = 0\nRenderView2.CameraClippingRange = [96.86482207477977, 292.5054971187886]\nRenderView2.RemoteRenderThreshold = 3.0\nRenderView2.Background = [1.0, 1.0, 1.0]\nRenderView2.CameraFocalPoint = [9.392949104309082, 1.5, 0.0]\nRenderView2.CenterAxesVisibility = 0\nRenderView2.CameraParallelScale = 57.398613649179126\nRenderView2.OrientationAxesLabelColor = [0.0, 0.0, 0.0]\n\nDataRepresentation2.EdgeColor = [0.0, 0.0, 0.0]\nDataRepresentation2.DiffuseColor = [0.0, 0.0, 0.0]\nDataRepresentation2.ColorArrayName = (\'POINT_DATA\', \'\')\nDataRepresentation2.AmbientColor = [0.0, 0.0, 0.0]\nDataRepresentation2.SelectionColor = [0.0, 0.0, 0.0]\nDataRepresentation2.BackfaceDiffuseColor = [0.0, 0.0, 0.0]\nDataRepresentation2.CubeAxesColor = [0.0, 0.0, 0.0]\nDataRepresentation2.Representation = \'Wireframe\'\n\na1_vtkEdgeFlags_PVLookupTable = GetLookupTableForArray( "vtkEdgeFlags", 1, RGBPoints=[0.0, 0.23, 0.299, 0.754, 0.5, 0.865, 0.865, 0.865, 1.0, 0.706, 0.016, 0.15], VectorMode=\'Magnitude\', NanColor=[0.25, 0.0, 0.0], ColorSpace=\'Diverging\', ScalarRangeInitialized=1.0 )\n\na1_vtkEdgeFlags_PiecewiseFunction = CreatePiecewiseFunction( Points=[0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0] )\n\nWriteImage(\'m1.png\')\nDataRepresentation2.ScalarOpacityFunction = a1_vtkEdgeFlags_PiecewiseFunction\nDataRepresentation2.LookupTable = a1_vtkEdgeFlags_PVLookupTable\n\na1_vtkEdgeFlags_PVLookupTable.ScalarOpacityFunction = a1_vtkEdgeFlags_PiecewiseFunction\n\n\n\nRenderView2.CameraViewUp = [0.0, 0.0, 1.0]\nRenderView2.CameraPosition = [9.392949104309082, -220.27121326772138, 0.0]\nRenderView2.CameraFocalPoint = [9.392949104309082, 1.5, 0.0]\nRenderView2.CameraClippingRange = [196.66850113504415, 253.90528146673722]\n\nWriteImage(\'m2.png\')\n\n\nRender()')


# In[6]:

get_ipython().system(u'python plot.py')


# 横と上から見たメッシュ図です。三脚のメッシュファイルであることがわかります。

# In[7]:

from IPython.core.display import Image
Image('m1.png')


# In[8]:

from IPython.core.display import Image
Image('m2.png')


# # FEMと立体求積法の設定
# 変数mにセットされたメッシュからそれぞれ変位応答格納用の変数と各節点データ操作用の変数を作成しておきます。GetFEM++の特徴として積分法の選択肢の多さがあります。

# In[9]:

mfu = gf.MeshFem(m,3) # displacement
mfd = gf.MeshFem(m,1) # data


# mfuとmfdにはそれぞれLagrange要素$Q_3$と$Q_1$が入ります。

# In[10]:

print mfu


# In[11]:

print mfd


# FEM手法として古典的なLagrange要素$P_k$を割り当てます。
# 
# | degree               | dimension            | d.o.f. number                        | class        | vectorial                | $\tau$-equivalent        | Polynomial   |
# |:--------------------:|:--------------------:|:------------------------------------:|:------------:|:------------------------:|:------------------------:|:------------:|
# | $K,0\leq K\leq255$   | $P,1\leq K\leq255$   | $\dfrac{\left(K+P\right)!}{K!P!}$    | $C^0$        | No$\left(Q=1\right)$     | Yes$\left(M=Id\right)$   | Yes          |
# 

# In[12]:

degree = 1
mfu.set_fem(gf.Fem('FEM_PK(3,%d)' % (degree,)));
mfd.set_fem(gf.Fem('FEM_PK(3,0)'))


# 立体求積法として15積分点・5次のtetrahedronを使用します。
# 
# <img src="getfemlistintmethodtetrahedron5.png">
# 

# In[13]:

mim = gf.MeshIm(m,gf.Integ('IM_TETRAHEDRON(5)'))


# # 境界条件の設定
# 最後に境界条件を設定します。今回は三脚の上端にNEUMANN条件を、下端にDIRICHLET条件を設定する。

# In[14]:

P = m.pts()
ctop = (abs(P[1,:] - 13) < 1e-6)
cbot = (abs(P[1,:] + 10) < 1e-6)
pidtop = np.compress(ctop,range(0,m.nbpts()))
pidbot = np.compress(cbot,range(0,m.nbpts()))
ftop = m.faces_from_pid(pidtop)
fbot = m.faces_from_pid(pidbot)


# 境界領域を作成します。

# In[15]:

NEUMANN_BOUNDARY = 1
DIRICHLET_BOUNDARY = 2
m.set_region(NEUMANN_BOUNDARY,ftop)
m.set_region(DIRICHLET_BOUNDARY,fbot)


# # 外力ベクトルと剛性マトリックスの組み立て
# 以上で設定した、条件を使用して外力ベクトルと剛性行列を設定します。

# In[16]:

help(gf.asm_boundary_source)


# In[17]:

help(gf.asm_linear_elasticity)


# In[18]:

nbd = mfd.nbdof()
F = gf.asm_boundary_source(NEUMANN_BOUNDARY, mim, mfu, mfd, np.repeat([[0],[-100],[0]],nbd,1))
K = gf.asm_linear_elasticity(mim, mfu, mfd, np.repeat([Lambda], nbd), np.repeat([Mu], nbd))


# # 固定条件の設定
# DIRICHLET条件(固定端条件)の設定をします。

# In[19]:

(H,R) = gf.asm_dirichlet(DIRICHLET_BOUNDARY, mim, mfu, mfd, mfd.eval('[[1,0,0],[0,1,0],[0,0,1]]'), mfd.eval('[0,0,0]'))
(N,U0) = H.dirichlet_nullspace(R)


# In[20]:

Nt = gf.Spmat('copy',N)
Nt.transpose()
KK = Nt*K*N
FF = Nt*(F-K*U0)


# In[21]:

# solve ...
P = gf.Precond('ildlt',KK)
UU = gf.linsolve_cg(KK,FF,P)
U = N*UU+U0


# # ポスト処理
# 以上で計算した変位をポスト処理結果として出力します。以下のコマンドで境界部分の変位の解析結果VTKファイルが出力されます。

# In[22]:

sl = gf.Slice(('boundary',), mfu, degree)
sl.export_to_vtk('tripod_ev.vtk', mfu, U, 'Displacement')


# In[23]:

get_ipython().run_cell_magic(u'writefile', u'plot.py', u'try: paraview.simple\nexcept: from paraview.simple import *\nparaview.simple._DisableFirstRenderCameraReset()\n\ntripod_ev_vtk = LegacyVTKReader( FileNames=[\'tripod_ev.vtk\'] )\n\nRenderView2 = CreateRenderView()\nRenderView2.CompressorConfig = \'vtkSquirtCompressor 0 3\'\nRenderView2.UseLight = 1\nRenderView2.LightSwitch = 0\nRenderView2.RemoteRenderThreshold = 3.0\nRenderView2.CenterOfRotation = [9.39224910736084, 1.5, 0.0]\n\nAnimationScene1 = GetAnimationScene()\nAnimationScene1.ViewModules = RenderView2\n\nDataRepresentation2 = Show()\nDataRepresentation2.ScaleFactor = 8.50093994140625\nDataRepresentation2.ScalarOpacityUnitDistance = 9.599809856069918\nDataRepresentation2.SelectionPointFieldDataArrayName = \'Displacement\'\nDataRepresentation2.EdgeColor = [0.0, 0.0, 0.0]\n\nRenderView2.CameraPosition = [9.39224910736084, 1.5, 221.76947835313635]\nRenderView2.OrientationAxesVisibility = 0\nRenderView2.CameraFocalPoint = [9.39224910736084, 1.5, 0.0]\nRenderView2.CameraClippingRange = [134.9674311526128, 331.57029329454673]\nRenderView2.Background = [1.0, 1.0, 1.0]\nRenderView2.CameraParallelScale = 57.398164620242895\nRenderView2.CenterAxesVisibility = 0\n\na3_Displacement_PVLookupTable = GetLookupTableForArray( "Displacement", 3, RGBPoints=[0.0, 0.23, 0.299, 0.754, 6.571885963503576, 0.865, 0.865, 0.865, 12.804577141990409, 0.706, 0.016, 0.15], VectorMode=\'Magnitude\', NanColor=[0.25, 0.0, 0.0], ColorSpace=\'Diverging\', ScalarRangeInitialized=1.0 )\n\na3_Displacement_PiecewiseFunction = CreatePiecewiseFunction( Points=[0.0, 0.0, 0.5, 0.0, 12.804577141990409, 1.0, 0.5, 0.0] )\n\nDataRepresentation2.Representation = \'Surface With Edges\'\nDataRepresentation2.ScalarOpacityFunction = a3_Displacement_PiecewiseFunction\nDataRepresentation2.ColorArrayName = (\'POINT_DATA\', \'Displacement\')\nDataRepresentation2.LookupTable = a3_Displacement_PVLookupTable\n\nWriteImage(\'tripod1.png\')\na3_Displacement_PVLookupTable.ScalarOpacityFunction = a3_Displacement_PiecewiseFunction\n\n\n\nRenderView2.CameraViewUp = [0.0, 0.0, 1.0]\nRenderView2.CameraPosition = [9.39224910736084, 223.26947835313635, 0.0]\nRenderView2.OrientationAxesVisibility = 0\nRenderView2.CameraClippingRange = [196.66678356960497, 253.9035205284334]\nRenderView2.CameraFocalPoint = [9.39224910736084, 1.5, 0.0]\nRenderView2.CenterAxesVisibility = 0\n\nWriteImage(\'tripod2.png\')\n\n\nRender()')


# In[24]:

get_ipython().system(u'python plot.py')


# In[25]:

from IPython.core.display import Image
Image('tripod1.png')


# In[26]:

from IPython.core.display import Image
Image('tripod2.png')


# In[26]:



