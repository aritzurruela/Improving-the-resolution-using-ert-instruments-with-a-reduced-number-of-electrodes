import pybert as pb
import pygimli as pg
import pygimli.meshtools
import matplotlib.pyplot as plt
import pygimli.meshtools as mt
import numpy as np

#inversion parameters

quality=34.5
maxCellArea=1
robustData=True
lam=20
paraDX=0.25
#creación grid para comparación

grid = pg.createGrid(x=pg.utils.grange(start=0.20, end=0.5, n=500),
                           y=-pg.utils.grange(0, 0.12, n=250, log=False))

grid2 = pg.createGrid(x=pg.utils.grange(start=0.33, end=0.35, n=500),
                           y=-pg.utils.grange(0.04, 0.06, n=250, log=False))

square=mt.createRectangle(start=[0.33,-0.04],end=[0.35,-0.06])

#MODELO

rhomap= [[1,150],[0,3000]]
background=mt.createWorld(start=[-1,0],end=[1,-1], area=1,marker=1)
circle=mt.createCircle(pos=[0.34,-0.05],radius=0.01,marker=0)
world = mt.mergePLC([background,circle])
mesh= mt.createMesh(world, quality=33, area=0.1,smooth=[1,2])

res_model=pg.solver.parseArgToArray(rhomap, mesh.cellCount(), mesh)

#INTERPOLACIÓN MODELO A GRID

nan=99.9
model_vector_nn = []
for pos in grid.cellCenters():
    cell = mesh.findCell(pos)
    if cell:
        model_vector_nn.append(res_model[cell.id()])
    else:
        model_vector_nn.append(nan) 
        
model_vector2 = []
for pos in grid2.cellCenters():
    cell = mesh.findCell(pos)
    if cell:
        model_vector2.append(res_model[cell.id()])
    else:
        model_vector2.append(nan)
        
model_data2=pb.pg.RVector(model_vector2)
        

#-------------------18x4 DIP-----------------------------------------

data_18x4_dip=pb.importData('D2_clean.bin')
data_18x4_dip.set('k', pb.geometricFactors(data_18x4_dip))
data_18x4_dip.set('rhoa', data_18x4_dip('u') / data_18x4_dip('i') * data_18x4_dip('k'))

#INVERSION

ert_18x4_dip =pb.ERTManager(data_18x4_dip)
#inv_18x4_dip =ert_18x4_dip.invert(quality=quality, maxCellArea=maxCellArea, robustData=robustData, lam=lam,paraDX=paraDX)
inv_18x4_dip =ert_18x4_dip.invert()
inversion_mesh_18x4_dip=ert_18x4_dip.paraDomain

#interpolación de datos inversión campo a grid

nan=99.9
inversion_vector_18x4_dip = []
for pos in grid.cellCenters():
    cell = inversion_mesh_18x4_dip.findCell(pos)
    if cell:
        inversion_vector_18x4_dip.append(inv_18x4_dip[cell.id()])
    else:
        inversion_vector_18x4_dip.append(nan) 
        
inversion_vector_18x4_dip2 = []
for pos in grid2.cellCenters():
    cell = inversion_mesh_18x4_dip.findCell(pos)
    if cell:
        inversion_vector_18x4_dip2.append(inv_18x4_dip[cell.id()])
    else:
        inversion_vector_18x4_dip2.append(nan) 
        
#comparación martorana

model_data_nn=pb.pg.RVector(model_vector_nn)
inversion_data_18x4_dip=pb.pg.RVector(inversion_vector_18x4_dip)
inversion_data_18x4_dip2=pb.pg.RVector(inversion_vector_18x4_dip2)
martorana_18x4_dip=np.log10(inversion_data_18x4_dip/model_data_nn)
martorana_18x4_dip2=np.log10(inversion_data_18x4_dip2/model_data2)

#valor abs
martorana_18x4_dip_abs=np.absolute(martorana_18x4_dip)
martorana_18x4_dip_abs2=np.absolute(martorana_18x4_dip2)
#stats
mean_18x4_dip=np.mean(martorana_18x4_dip_abs)
mean_18x4_dip2=np.mean(martorana_18x4_dip_abs2)

#------------------18x4 SCH---------------------

data_18x4_sch=pb.importData('S2_clean.bin')
data_18x4_sch.set('k', pb.geometricFactors(data_18x4_sch))
data_18x4_sch.set('rhoa', data_18x4_sch('u') / data_18x4_sch('i') * data_18x4_sch('k'))

#INVERSION

ert_18x4_sch =pb.ERTManager(data_18x4_sch)
inv_18x4_sch =ert_18x4_sch.invert()
inversion_mesh_18x4_sch=ert_18x4_sch.paraDomain

#interpolación de datos inversión campo a grid

nan=99.9
inversion_vector_18x4_sch = []
for pos in grid.cellCenters():
    cell = inversion_mesh_18x4_sch.findCell(pos)
    if cell:
        inversion_vector_18x4_sch.append(inv_18x4_sch[cell.id()])
    else:
        inversion_vector_18x4_sch.append(nan) 
        
inversion_vector_18x4_sch2 = []
for pos in grid2.cellCenters():
    cell = inversion_mesh_18x4_sch.findCell(pos)
    if cell:
        inversion_vector_18x4_sch2.append(inv_18x4_sch[cell.id()])
    else:
        inversion_vector_18x4_sch2.append(nan) 
#comparación martorana

inversion_data_18x4_sch=pb.pg.RVector(inversion_vector_18x4_sch)
inversion_data_18x4_sch2=pb.pg.RVector(inversion_vector_18x4_sch2)
martorana_18x4_sch=np.log10(inversion_data_18x4_sch/model_data_nn)
martorana_18x4_sch2=np.log10(inversion_data_18x4_sch2/model_data2)

#valor abs
martorana_18x4_sch_abs=np.absolute(martorana_18x4_sch)
martorana_18x4_sch_abs2=np.absolute(martorana_18x4_sch2)
#stats
mean_18x4_sch2=np.mean(martorana_18x4_sch_abs2)
mean_18x4_sch=np.mean(martorana_18x4_sch_abs)

#-------------18x1 CONC DIP------------


#importing data

data_CONC_dip=pb.importData('d_total.bin')
data_CONC_dip.set('k', pb.geometricFactors(data_CONC_dip))
data_CONC_dip.set('rhoa', data_CONC_dip('u') / data_CONC_dip('i') * data_CONC_dip('k'))

#inverting

ert_conc_dip =pb.ERTManager(data_CONC_dip)
#inv_conc_dip =ert_conc_dip.invert(quality=quality, maxCellArea=maxCellArea, robustData=robustData, lam=lam,paraDX=paraDX)
inv_conc_dip =ert_conc_dip.invert()
inversion_mesh_conc_dip=ert_conc_dip.paraDomain

#interpolación de datos inversión campo a grid

nan=99.9
inversion_vector_conc_dip = []
for pos in grid.cellCenters():
    cell = inversion_mesh_conc_dip.findCell(pos)
    if cell:
        inversion_vector_conc_dip.append(inv_conc_dip[cell.id()])
    else:
        inversion_vector_conc_dip.append(nan)

inversion_vector_conc_dip2 = []
for pos in grid2.cellCenters():
    cell = inversion_mesh_conc_dip.findCell(pos)
    if cell:
        inversion_vector_conc_dip2.append(inv_conc_dip[cell.id()])
    else:
        inversion_vector_conc_dip2.append(nan) 
        
#comparación martorana

inversion_data_conc_dip=pb.pg.RVector(inversion_vector_conc_dip)
inversion_data_conc_dip2=pb.pg.RVector(inversion_vector_conc_dip2)
martorana_conc_dip=np.log10(inversion_data_conc_dip/model_data_nn)
martorana_conc_dip2=np.log10(inversion_data_conc_dip2/model_data2)

#valor abs
martorana_conc_dip_abs=np.absolute(martorana_conc_dip)
martorana_conc_dip_abs2=np.absolute(martorana_conc_dip2)
#stats
mean_conc_dip=np.mean(martorana_conc_dip_abs)
mean_conc_dip2=np.mean(martorana_conc_dip_abs2)

#-------------18x1 CONC SCH-----------

data_conc_sch=pb.importData('s_total.bin')

data_conc_sch.set('k', pb.geometricFactors(data_conc_sch))
data_conc_sch.set('rhoa', data_conc_sch('u') / data_conc_sch('i') * data_conc_sch('k'))
#inverting

ert_conc_sch =pb.ERTManager(data_conc_sch)
inv_conc_sch =ert_conc_sch.invert()
inversion_mesh_conc_sch=ert_conc_sch.paraDomain

#interpolación de datos inversión campo a grid

nan=99.9
inversion_vector_conc_sch = []
for pos in grid.cellCenters():
    cell = inversion_mesh_conc_sch.findCell(pos)
    if cell:
        inversion_vector_conc_sch.append(inv_conc_sch[cell.id()])
    else:
        inversion_vector_conc_sch.append(nan) 

inversion_vector_conc_sch2 = []
for pos in grid2.cellCenters():
    cell = inversion_mesh_conc_sch.findCell(pos)
    if cell:
        inversion_vector_conc_sch2.append(inv_conc_sch[cell.id()])
    else:
        inversion_vector_conc_sch2.append(nan)        
#comparación martorana

inversion_data_conc_sch2=pb.pg.RVector(inversion_vector_conc_sch2)
inversion_data_conc_sch=pb.pg.RVector(inversion_vector_conc_sch)
martorana_conc_sch=np.log10(inversion_data_conc_sch/model_data_nn)
martorana_conc_sch2=np.log10(inversion_data_conc_sch2/model_data2)

#valor abs
martorana_conc_sch_abs=np.absolute(martorana_conc_sch)
martorana_conc_sch_abs2=np.absolute(martorana_conc_sch2)
#stats
mean_conc_sch=np.mean(martorana_conc_sch_abs)
mean_conc_sch2=np.mean(martorana_conc_sch_abs2)

#-----------------72x1 DIP -----------------

data_72x1_dip=pb.importData('72X1_DIP.bin')


data_72x1_dip.set('k', pb.geometricFactors(data_72x1_dip))
data_72x1_dip.set('rhoa', data_72x1_dip('u') / data_72x1_dip('i') * data_72x1_dip('k'))

#inversion

ert_72x1_dip =pb.ERTManager(data_72x1_dip)

#inv_72x1_dip =ert_72x1_dip.invert(quality=quality, maxCellArea=maxCellArea, robustData=robustData, lam=lam,paraDX=paraDX)
inv_72x1_dip =ert_72x1_dip.invert()
inversion_mesh_72x1_dip=ert_72x1_dip.paraDomain

#interpolación de datos inversión campo a grid

nan=99.9
inversion_vector_72x1_dip = []
for pos in grid.cellCenters():
    cell = inversion_mesh_72x1_dip.findCell(pos)
    if cell:
        inversion_vector_72x1_dip.append(inv_72x1_dip[cell.id()])
    else:
        inversion_vector_72x1_dip.append(nan) 

inversion_vector_72x1_dip2 = []
for pos in grid2.cellCenters():
    cell = inversion_mesh_72x1_dip.findCell(pos)
    if cell:
        inversion_vector_72x1_dip2.append(inv_72x1_dip[cell.id()])
    else:
        inversion_vector_72x1_dip2.append(nan)       
#comparación martorana

inversion_data_72x1_dip=pb.pg.RVector(inversion_vector_72x1_dip)
martorana_72x1_dip=np.log10(inversion_data_72x1_dip/model_data_nn)

inversion_data_72x1_dip2=pb.pg.RVector(inversion_vector_72x1_dip2)
martorana_72x1_dip2=np.log10(inversion_data_72x1_dip2/model_data2)

#valor abs
martorana_72x1_dip_abs=np.absolute(martorana_72x1_dip)
martorana_72x1_dip_abs2=np.absolute(martorana_72x1_dip2)
#stats
mean_72x1_dip=np.mean(martorana_72x1_dip_abs)
mean_72x1_dip2=np.mean(martorana_72x1_dip_abs2)

#-----------------72x1 SCH -----------------

data_72X1_sch=pb.importData('72X1_SCH.bin')

data_72X1_sch.set('k', pb.geometricFactors(data_72X1_sch))
data_72X1_sch.set('rhoa', data_72X1_sch('u') / data_72X1_sch('i') * data_72X1_sch('k'))

#inversion

ert_72x1_sch =pb.ERTManager(data_72X1_sch)

#inv_72x1_sch =ert_72x1_sch.invert(quality=quality, maxCellArea=maxCellArea, robustData=robustData, lam=lam,paraDX=paraDX)
inv_72x1_sch =ert_72x1_sch.invert()
inversion_mesh_72x1_sch=ert_72x1_sch.paraDomain

#interpolación de datos inversión campo a grid

nan=99.9
inversion_vector_72x1_sch = []
for pos in grid.cellCenters():
    cell = inversion_mesh_72x1_sch.findCell(pos)
    if cell:
        inversion_vector_72x1_sch.append(inv_72x1_sch[cell.id()])
    else:
        inversion_vector_72x1_sch.append(nan) 

inversion_vector_72x1_sch2 = []
for pos in grid2.cellCenters():
    cell = inversion_mesh_72x1_sch.findCell(pos)
    if cell:
        inversion_vector_72x1_sch2.append(inv_72x1_sch[cell.id()])
    else:
        inversion_vector_72x1_sch2.append(nan)         
#comparación martorana

inversion_data_72x1_sch=pb.pg.RVector(inversion_vector_72x1_sch)
martorana_72x1_sch=np.log10(inversion_data_72x1_sch/model_data_nn)

inversion_data_72x1_sch2=pb.pg.RVector(inversion_vector_72x1_sch2)
martorana_72x1_sch2=np.log10(inversion_data_72x1_sch2/model_data2)
#valor abs
martorana_72x1_sch_abs=np.absolute(martorana_72x1_sch)
martorana_72x1_sch_abs2=np.absolute(martorana_72x1_sch2)
#stats
mean_72x1_sch=np.mean(martorana_72x1_sch_abs)
mean_72x1_sch2=np.mean(martorana_72x1_sch_abs2)

#martorana vs martorana

martorana_18x4vs_conc_dip=np.greater(martorana_conc_dip_abs,martorana_18x4_dip_abs,dtype=object)
martorana_18x4vs_conc_sch=np.greater(martorana_conc_sch_abs,martorana_18x4_sch_abs,dtype=object)
martorana_conc_vs72x1_dip=np.greater(martorana_72x1_dip_abs,martorana_conc_dip_abs,dtype=object)
martorana_conc_vs72x1_sch=np.greater(martorana_72x1_sch_abs,martorana_conc_sch_abs,dtype=object)

mean_72x1_dd=np.mean(martorana_72x1_dip_abs)
mean_72x1_sch=np.mean(martorana_72x1_sch_abs)
mean_18x1_dd=np.mean(martorana_conc_dip_abs)
mean_18x1_sch=np.mean(martorana_conc_sch_abs)
mean_18x4_dd=np.mean(martorana_18x4_dip_abs)
mean_18x4_sch=np.mean(martorana_18x4_sch_abs)
stats=open('mitjes_martorana_sim2_ref.txt','w')
stats.write('la mitjana del perfil simple de 72 (dipol) és:'+str(mean_72x1_dd)+'\n')
stats.write('la mitjana del perfil simple de 72 (schlumb) és:'+str(mean_72x1_sch)+'\n')
stats.write('la mitjana del perfil simple de 18 (dipol) és:'+str(mean_18x4_dd)+'\n')
stats.write('la mitjana del perfil simple de 18 (schlumb) és:'+str(mean_18x4_sch)+'\n')
stats.write('la mitjana del perfil concatenat (dipol) és:'+str(mean_18x1_dd)+'\n')
stats.write('la mitjana del perfil concatenat (schlum) és:'+str(mean_18x1_sch)+'\n')
stats.close()


mean_72x1_dd_s=np.mean(martorana_72x1_dip_abs2)
mean_72x1_sch_s=np.mean(martorana_72x1_sch_abs2)
mean_18x1_dd_s=np.mean(martorana_conc_dip_abs2)
mean_18x1_sch_s=np.mean(martorana_conc_sch_abs2)
mean_18x4_dd_s=np.mean(martorana_18x4_dip_abs2)
mean_18x4_sch_s=np.mean(martorana_18x4_sch_abs2)
stats=open('mitjes_martorana_sim_s_ref.txt','w')
stats.write('la mitjana del perfil simple de 72 (dipol) és:'+str(mean_72x1_dd_s)+'\n')
stats.write('la mitjana del perfil simple de 72 (schlumb) és:'+str(mean_72x1_sch_s)+'\n')
stats.write('la mitjana del perfil simple de 18 (dipol) és:'+str(mean_18x4_dd_s)+'\n')
stats.write('la mitjana del perfil simple de 18 (schlumb) és:'+str(mean_18x4_sch_s)+'\n')
stats.write('la mitjana del perfil concatenat (dipol) és:'+str(mean_18x1_dd_s)+'\n')
stats.write('la mitjana del perfil concatenat (schlum) és:'+str(mean_18x1_sch_s)+'\n')
stats.close()

#FIGURES

#figura inversió

cmin=60
cmax=1000

fig, (ax) = plt.subplots(ncols=2,nrows=4,figsize=(20,20),dpi=300)
plt.subplots_adjust(hspace=1,wspace=1)
pg.show(world,fillRegion=False,ax=ax[0,0])
ert_18x4_dip.showResult(cMin=cmin,cMax=cmax,ax=ax[0,0])
pg.show(world,fillRegion=False,ax=ax[1,0])
ert_18x4_sch.showResult(cMin=cmin,cMax=cmax,ax=ax[1,0])
pg.show(world,fillRegion=False,ax=ax[2,0])
ert_conc_dip.showResult(cMin=cmin,cMax=cmax,ax=ax[2,0])
pg.show(world,fillRegion=False,ax=ax[3,0])
ert_conc_sch.showResult(cMin=cmin,cMax=cmax,ax=ax[3,0])
pg.show(world,fillRegion=False,ax=ax[0,1])
ert_72x1_dip.showResult(cMin=cmin,cMax=cmax,ax=ax[0,1])
pg.show(world,fillRegion=False,ax=ax[1,1])
ert_72x1_sch.showResult(cMin=cmin,cMax=cmax,ax=ax[1,1])
pg.show(world,ax=ax[2,1])
pg.show(mesh,rhomap,label='Resistivity $(\Omega$m)',ax=ax[3,1])
#ax[0,0].set(xlim=(20,50),ylim=(-10,0))
ax[0,0].set_title('18x4_dd')
ax[0,0].set_ylabel('Depth')
ax[0,0].set_xlabel('distance (m)')
#ax[1,0].set(xlim=(20,50),ylim=(-10,0))
ax[1,0].set_title('18x4_slm')
ax[1,0].set_ylabel('Depth')
ax[1,0].set_xlabel('distance (m)')
#ax[2,0].set(xlim=(20,50),ylim=(-10,0))
ax[2,0].set_title('18x1_dd_conc')
ax[2,0].set_ylabel('Depth')
ax[2,0].set_xlabel('distance (m)')
#ax[3,0].set(xlim=(20,50),ylim=(-10,0))
ax[3,0].set_title('18x1_slm_conc')
ax[3,0].set_ylabel('Depth')
ax[3,0].set_xlabel('distance (m)')
#ax[0,1].set(xlim=(20,50),ylim=(-10,0))
ax[0,1].set_title('72x1_dd')
ax[0,1].set_ylabel('Depth')
ax[0,1].set_xlabel('distance (m)')
#ax[1,1].set(xlim=(20,50),ylim=(-10,0))
ax[1,1].set_title('72x1_slm')
ax[1,1].set_ylabel('Depth')
ax[1,1].set_xlabel('distance (m)')
ax[3,1].set(xlim=(0,0.71),ylim=(-0.10,0))
ax[2,1].set_ylabel('Depth')
ax[2,1].set_xlabel('distance (m)')
ax[3,1].set(xlim=(0,0.71),ylim=(-0.10,0))
ax[3,1].set_title('model')
ax[3,1].set_ylabel('Depth')
ax[3,1].set_xlabel('distance (m)')

plt.savefig('inversion.png')


#figura mortoranes

fig, (ax) = plt.subplots(ncols=2,nrows=4,figsize=(20,20),dpi=300)
#plt.subplots_adjust(hspace=0.5)
pg.show(grid,martorana_18x4_dip,cMap='bwr',cMin=-1,cMax=1,ax=ax[0,0])
pg.show(grid,martorana_18x4_sch,cMap='bwr',cMin=-1,cMax=1,ax=ax[1,0])
pg.show(grid,martorana_conc_dip,cMap='bwr',cMin=-1,cMax=1,ax=ax[2,0])
pg.show(grid,martorana_conc_sch,cMap='bwr',cMin=-1,cMax=1,ax=ax[3,0])
pg.show(grid,martorana_72x1_dip,cMap='bwr',cMin=-1,cMax=1,ax=ax[0,1])
pg.show(grid,martorana_72x1_sch,cMap='bwr',cMin=-1,cMax=1,ax=ax[1,1])
pg.show(grid,martorana_conc_dip,cMap='bwr',cMin=-1,cMax=1,ax=ax[2,1])
pg.show(grid,martorana_conc_sch,cMap='bwr',cMin=-1,cMax=1,ax=ax[3,1])
ax[0,0].set_title('18x4_dd')
ax[0,0].set_ylabel('Depth')
ax[0,0].set_xlabel('distance (m)')
ax[0,0].text(0.34, -0.05,str(round(mean_18x4_dip,3)),verticalalignment='bottom', horizontalalignment='right',color='red',fontsize=30)
ax[1,0].set_title('18x4_slm')
ax[1,0].set_ylabel('Depth')
ax[1,0].set_xlabel('distance (m)')
ax[1,0].text(0.34, -0.05,str(round(mean_18x4_sch,3)),verticalalignment='bottom', horizontalalignment='right',color='red',fontsize=30)
ax[2,0].set_title('18x1_dd_conc')
ax[2,0].set_ylabel('Depth')
ax[2,0].set_xlabel('distance (m)')
ax[2,0].text(0.34, -0.05,str(round(mean_conc_dip,3)),verticalalignment='bottom', horizontalalignment='right',color='red',fontsize=30)
ax[3,0].set_title('18x1_slm_conc')
ax[3,0].set_ylabel('Depth')
ax[3,0].set_xlabel('distance (m)')
ax[3,0].text(0.34, -0.05,str(round(mean_conc_sch,3)),verticalalignment='bottom', horizontalalignment='right',color='red',fontsize=30)
ax[0,1].set_title('72x1_dd')
ax[0,1].set_ylabel('Depth')
ax[0,1].set_xlabel('distance (m)')
ax[0,1].text(0.34, -0.05,str(round(mean_72x1_dip,3)),verticalalignment='bottom', horizontalalignment='right',color='red',fontsize=30)
ax[1,1].set_title('72x1_slm')
ax[1,1].set_ylabel('Depth')
ax[1,1].set_xlabel('distance (m)')
ax[1,1].text(0.34, -0.05,str(round(mean_72x1_sch,3)),verticalalignment='bottom', horizontalalignment='right',color='red',fontsize=30)
ax[2,1].set_title('18x1_dd_conc')
ax[2,1].set_ylabel('Depth')
ax[2,1].set_xlabel('distance (m)')
ax[2,1].text(0.34, -0.05,str(round(mean_conc_dip,3)),verticalalignment='bottom', horizontalalignment='right',color='red',fontsize=30)
ax[3,1].set_title('18x1_slm_conc')
ax[3,1].set_ylabel('Depth')
ax[3,1].set_xlabel('distance (m)')
ax[3,1].text(0.34, -0.05,str(round(mean_conc_sch,3)),verticalalignment='bottom', horizontalalignment='right',color='red',fontsize=30)
plt.savefig('martorana_4x4_fieldb.png')

#grid pequeño
fig, (ax) = plt.subplots(ncols=2,nrows=4,figsize=(20,20),dpi=300)
#plt.subplots_adjust(hspace=0.5)
pg.show(grid2,martorana_18x4_dip2,cMap='bwr',cMin=-1,cMax=1,ax=ax[0,0])
pg.show(grid2,martorana_18x4_sch2,cMap='bwr',cMin=-1,cMax=1,ax=ax[1,0])
pg.show(grid2,martorana_conc_dip2,cMap='bwr',cMin=-1,cMax=1,ax=ax[2,0])
pg.show(grid2,martorana_conc_sch2,cMap='bwr',cMin=-1,cMax=1,ax=ax[3,0])
pg.show(grid2,martorana_72x1_dip2,cMap='bwr',cMin=-1,cMax=1,ax=ax[0,1])
pg.show(grid2,martorana_72x1_sch2,cMap='bwr',cMin=-1,cMax=1,ax=ax[1,1])
pg.show(grid2,martorana_conc_dip2,cMap='bwr',cMin=-1,cMax=1,ax=ax[2,1])
pg.show(grid2,martorana_conc_sch2,cMap='bwr',cMin=-1,cMax=1,ax=ax[3,1])
ax[0,0].set_title('18x4_dd')
ax[0,0].set_ylabel('Depth')
ax[0,0].set_xlabel('distance (m)')
ax[0,0].text(0.34, -0.05,str(round(mean_18x4_dip2,3)),verticalalignment='bottom', horizontalalignment='right',color='white',fontsize=30)
ax[1,0].set_title('18x4_slm')
ax[1,0].set_ylabel('Depth')
ax[1,0].set_xlabel('distance (m)')
ax[1,0].text(0.34, -0.05,str(round(mean_18x4_sch2,3)),verticalalignment='bottom', horizontalalignment='right',color='white',fontsize=30)
ax[2,0].set_title('18x1_dd_conc')
ax[2,0].set_ylabel('Depth')
ax[2,0].set_xlabel('distance (m)')
ax[2,0].text(0.34, -0.05,str(round(mean_conc_dip2,3)),verticalalignment='bottom', horizontalalignment='right',color='white',fontsize=30)
ax[3,0].set_title('18x1_slm_conc')
ax[3,0].set_ylabel('Depth')
ax[3,0].set_xlabel('distance (m)')
ax[3,0].text(0.34, -0.05,str(round(mean_conc_sch2,3)),verticalalignment='bottom', horizontalalignment='right',color='white',fontsize=30)
ax[0,1].set_title('72x1_dd')
ax[0,1].set_ylabel('Depth')
ax[0,1].set_xlabel('distance (m)')
ax[0,1].text(0.34, -0.05,str(round(mean_72x1_dip2,3)),verticalalignment='bottom', horizontalalignment='right',color='white',fontsize=30)
ax[1,1].set_title('72x1_slm')
ax[1,1].set_ylabel('Depth')
ax[1,1].set_xlabel('distance (m)')
ax[1,1].text(0.34, -0.05,str(round(mean_72x1_sch2,3)),verticalalignment='bottom', horizontalalignment='right',color='white',fontsize=30)
ax[2,1].set_title('18x1_dd_conc')
ax[2,1].set_ylabel('Depth')
ax[2,1].set_xlabel('distance (m)')
ax[2,1].text(0.34, -0.05,str(round(mean_conc_dip2,3)),verticalalignment='bottom', horizontalalignment='right',color='white',fontsize=30)
ax[3,1].set_title('18x1_slm_conc')
ax[3,1].set_ylabel('Depth')
ax[3,1].set_xlabel('distance (m)')
ax[3,1].text(0.34, -0.05,str(round(mean_conc_sch2,3)),verticalalignment='bottom', horizontalalignment='right',color='white',fontsize=30)
plt.savefig('martorana_4x4_field_smallb.png')

#figura martorana inversió

cmin=60
cmax=1500

fig, (ax) = plt.subplots(ncols=2,nrows=3,figsize=(20,15),dpi=300)
#plt.subplots_adjust(hspace=-2)
pg.show(world,fillRegion=False,ax=ax[0,0])
ert_18x4_dip.showResult(cMin=cmin,cMax=cmax,colorBar=False,ax=ax[0,0])
pg.show(world,fillRegion=False,ax=ax[1,0])
ert_conc_dip.showResult(cMin=cmin,cMax=cmax,colorBar=False,ax=ax[1,0])
pg.show(world,fillRegion=False,ax=ax[2,0])
ert_72x1_dip.showResult(cMin=cmin,cMax=cmax,colorBar=True,ax=ax[2,0])
pg.show(square,fillRegion=False,ax=ax[0,1])
pg.show(grid,martorana_18x4_dip,cMap='bwr',cMin=-2,cMax=2,ax=ax[0,1])
pg.show(square,fillRegion=False,ax=ax[1,1])
pg.show(grid,martorana_conc_dip,cMap='bwr',cMin=-2,cMax=2,ax=ax[1,1])
pg.show(square,fillRegion=False,ax=ax[2,1])
pg.show(grid,martorana_72x1_dip,cMap='bwr',label='Adj. index',cMin=-2,cMax=2,ax=ax[2,1])
ax[0,0].set(ylim=(-0.12,0))
ax[0,0].set(xlim=(0.2,0.5))
#ax[0,0].set_title('18x4_dd')
ax[0,0].set_ylabel('Depth')
ax[0,0].set_xlabel('distance (m)')
ax[1,0].set(ylim=(-0.12,0))
ax[1,0].set(xlim=(0.2,0.5))
#ax[1,0].set_title('conc_dd')
ax[1,0].set_ylabel('Depth')
ax[1,0].set_xlabel('distance (m)')
ax[2,0].set(ylim=(-0.12,0))
ax[2,0].set(xlim=(0.2,0.5))
#ax[2,0].set_title('72x1_dd')
ax[2,0].set_ylabel('Depth')
ax[2,0].set_xlabel('distance (m)')
ax[0,1].text(0.38, -0.07,str(round(mean_18x4_dip2,3)),size=20,verticalalignment='bottom', horizontalalignment='right',bbox={'facecolor': 'red', 'alpha': 0.2, 'pad': 5})
#ax[0,1].set(ylim=(-15,0))
#ax[0,1].set(xlim=(0,71))
#ax[0,1].set_title('martorana_18x1')
#ax[0,1].set_ylabel('Depth')
ax[0,1].set_xlabel('distance(m)')
ax[1,1].text(0.38, -0.07,str(round(mean_conc_dip2,3)),size=20,verticalalignment='bottom', horizontalalignment='right',bbox={'facecolor': 'red', 'alpha': 0.2, 'pad': 5})
#ax[1,1].set(ylim=(-15,0))
#ax[1,1].set(xlim=(0,71))
#ax[1,1].set_title('martorana_72')
#ax[1,1].set_ylabel('Depth')
ax[1,1].set_xlabel('distance (m)')
ax[2,1].text(0.38, -0.07,str(round(mean_72x1_dip2,3)),size=20,verticalalignment='bottom', horizontalalignment='right',bbox={'facecolor': 'red', 'alpha': 0.2, 'pad': 5})
#ax[2,1].set(ylim=(-15,0))
#ax[2,1].set(xlim=(0,71))
#ax[2,1].set_ylabel('Depth')
ax[2,1].set_xlabel('distance (m)')

plt.savefig('inversion_martorana_dd.png')

cmin=50
cmax=750

fig, (ax) = plt.subplots(ncols=2,nrows=3,figsize=(45,20),dpi=300)
#plt.subplots_adjust(hspace=-2)
pg.show(world,fillRegion=False,ax=ax[0,0])
ert_18x4_sch.showResult(cMin=cmin,cMax=cmax,colorBar=False,ax=ax[0,0])
pg.show(world,fillRegion=False,ax=ax[1,0])
ert_conc_sch.showResult(cMin=cmin,cMax=cmax,colorBar=False,ax=ax[1,0])
pg.show(world,fillRegion=False,ax=ax[2,0])
ert_72x1_sch.showResult(cMin=cmin,cMax=cmax,colorBar=True,ax=ax[2,0])
pg.show(world,fillRegion=False,ax=ax[0,1])
pg.show(grid2,martorana_18x4_sch2,cMap='bwr',cMin=-1,cMax=1,ax=ax[0,1])
pg.show(world,fillRegion=False,ax=ax[1,1])
pg.show(grid2,martorana_conc_sch2,cMap='bwr',cMin=-1,cMax=1,ax=ax[1,1])
pg.show(world,fillRegion=False,ax=ax[2,1])
pg.show(grid2,martorana_72x1_sch2,cMap='bwr',label='Adj. index',cMin=-2,cMax=2,ax=ax[2,1])
#ax[0,0].set(ylim=(-15,0))
#ax[0,0].set(xlim=(0,71))
ax[0,0].set_title('18x4_dd')
ax[0,0].set_ylabel('Depth')
ax[0,0].set_xlabel('distance (m)')
#ax[1,0].set(ylim=(-15,0))
#ax[1,0].set(xlim=(0,71))
ax[1,0].set_title('conc_dd')
ax[1,0].set_ylabel('Depth')
ax[1,0].set_xlabel('distance (m)')
#ax[2,0].set(ylim=(-15,0))
#ax[2,0].set(xlim=(0,71))
ax[2,0].set_title('72x1_dd')
ax[2,0].set_ylabel('Depth')
ax[2,0].set_xlabel('distance (m)')
ax[0,1].text(0.34, -0.05,str(round(mean_18x4_sch2,3)),size=30,verticalalignment='bottom', horizontalalignment='right',bbox={'facecolor': 'red', 'alpha': 0.2, 'pad': 5})
#ax[0,1].set(ylim=(-15,0))
#ax[0,1].set(xlim=(0,71))
ax[0,1].set_title('martorana_18x1')
ax[0,1].set_ylabel('Depth')
ax[0,1].set_xlabel('martorana_conc')
ax[1,1].text(0.34, -0.05,str(round(mean_conc_sch2,3)),size=30,verticalalignment='bottom', horizontalalignment='right',bbox={'facecolor': 'red', 'alpha': 0.2, 'pad': 5})
#ax[1,1].set(ylim=(-15,0))
#ax[1,1].set(xlim=(0,71))
ax[1,1].set_title('martorana_72')
ax[1,1].set_ylabel('Depth')
ax[1,1].set_xlabel('distance (m)')
ax[2,1].text(0.34, -0.05,str(round(mean_72x1_sch2,3)),size=30,verticalalignment='bottom', horizontalalignment='right',bbox={'facecolor': 'red', 'alpha': 0.2, 'pad': 5})
#ax[2,1].set(ylim=(-15,0))
#ax[2,1].set(xlim=(0,71))
ax[2,1].set_ylabel('Depth')
ax[2,1].set_xlabel('distance (m)')

plt.savefig('inversion_martorana_sch.png')

#figura martorana vs martorana

fig, (ax) = plt.subplots(ncols=2,nrows=2, figsize=(20,20),dpi=300)
pg.show(world,fillRegion=False,ax=ax[0,0])
pg.show(grid,martorana_18x4vs_conc_dip,cMap='bwr',cMin=-1,cMax=1,ax=ax[0,0])
pg.show(world,fillRegion=False,ax=ax[1,0])
pg.show(grid,martorana_18x4vs_conc_sch,cMap='bwr',cMin=-1,cMax=1,ax=ax[1,0])
pg.show(world,fillRegion=False,ax=ax[0,1])
pg.show(grid,martorana_conc_vs72x1_dip,cMap='bwr',cMin=-1,cMax=1,ax=ax[0,1])
pg.show(world,fillRegion=False,ax=ax[1,1])
pg.show(grid,martorana_conc_vs72x1_sch,cMap='bwr',cMin=-1,cMax=1,ax=ax[1,1])
ax[0,0].set_title('perfil simple vs. concatenació dipol')
ax[0,0].set_ylabel('Depth')
ax[0,0].set_xlabel('distance (m)')
ax[1,0].set_title('perfil simple vs. concatenació schlumberger')
ax[1,0].set_ylabel('Depth')
ax[1,0].set_xlabel('distance (m)')
ax[0,1].set_title('perfil 72 vs. concatenació dipol')
ax[0,1].set_ylabel('Depth')
ax[0,1].set_xlabel('distance (m)')
ax[1,1].set_title('perfil 72 vs. concatenació schlumberger')
ax[1,1].set_ylabel('Depth')
ax[1,1].set_xlabel('distance (m)')
plt.savefig('martorana_vs_martoranab.png')

fig, (ax) = plt.subplots(figsize=(10,5),dpi=300)
pg.show(mesh,rhomap,showMesh=True,label='Resistivity $(\Omega$m)',cMap='Reds',ax=ax)
ax.set(xlim=(0.20,0.48),ylim=(-0.10,0))
#ax.set_title('model')
ax.set_ylabel('Depth (m)')
ax.set_xlabel('distance (m)')
plt.savefig('model_red.png')
