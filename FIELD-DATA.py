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
lam=50
paraDX=0.25
#creación grid para comparación

grid = pg.createGrid(x=pg.utils.grange(start=30, end=40, n=500),
                           y=-pg.utils.grange(0, 4, n=250, log=False))

grid2 = pg.createGrid(x=pg.utils.grange(start=34, end=36, n=500),
                           y=-pg.utils.grange(0.5, 2.5, n=250, log=False))

#MODELO

rhomap= [[1,50],[0,200],[3,700],[2,75],[4,100]]
background=mt.createWorld(start=[-10,0],end=[81,-21], area=1,marker=1)
sup=mt.createPolygon([[-10,0],[30,0],[31.2,-1.5],[-10,-1.5]], isClosed=True,area=1,marker=2)
sup2=mt.createPolygon([[40,0],[81,0],[81,-1.5],[38.8,-1.5]], isClosed=True,area=1,marker=4)
sup3=mt.createPolygon([[30,0],[40,0],[39.6,-0.5],[30.4,-0.5]], isClosed=True,area=1,marker=4)
pol=mt.createPolygon([[30,0],[40,0],[38,-2.5],[32,-2.5]], isClosed=True,area=1,marker=0)
square= mt.createRectangle(start=[34,-0.5],end=[36,-2.5],area=1, marker=3)
world = mt.mergePLC([background,pol,square,sup,sup2,sup3])
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

data_18x4_dip=pb.importData('18X4_DIP_1A.bin')
data_18x4_dip.set('k', pb.geometricFactors(data_18x4_dip))
data_18x4_dip.set('rhoa', data_18x4_dip('u') / data_18x4_dip('i') * data_18x4_dip('k'))

#INVERSION

ert_18x4_dip =pb.ERTManager(data_18x4_dip)
inv_18x4_dip =ert_18x4_dip.invert(quality=quality, maxCellArea=maxCellArea, robustData=robustData, lam=lam,paraDX=paraDX)
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

data_18x4_sch=pb.importData('18X4 1SA_5187.bin')
data_18x4_sch.set('k', pb.geometricFactors(data_18x4_sch))
data_18x4_sch.set('rhoa', data_18x4_sch('u') / data_18x4_sch('i') * data_18x4_sch('k'))

#INVERSION

ert_18x4_sch =pb.ERTManager(data_18x4_sch)
inv_18x4_sch =ert_18x4_sch.invert(quality=quality, maxCellArea=maxCellArea, robustData=robustData, lam=lam,paraDX=paraDX)
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

data1=pb.importData('18X4_DIP_1A.bin')
data2=pb.importData('18X4_DIP_1B.bin')
data3=pb.importData('18X4_DIP_1C.bin')
data4=pb.importData('18X4_DIP_1D.bin')

#MODIFY ELECTRODE POSITION

for i in range(data2.sensorCount()):
    data2.setSensorPosition(i, [i*4+1, 0, 0])
    
for i in range(data3.sensorCount()):
    data3.setSensorPosition(i, [i*4+2, 0, 0])
    
for i in range(data4.sensorCount()):
    data4.setSensorPosition(i, [i*4+3, 0, 0])

#calculating rhoa

data1.set('k', pb.geometricFactors(data1))
data1.set('rhoa', data1('u') / data1('i') * data1('k'))

data2.set('k', pb.geometricFactors(data2))
data2.set('rhoa', data2('u') / data2('i') * data2('k'))

data3.set('k', pb.geometricFactors(data3))
data3.set('rhoa', data3('u') / data3('i') * data3('k'))

data4.set('k', pb.geometricFactors(data4))
data4.set('rhoa', data4('u') / data4('i') * data4('k'))

#concatenating data

data_CONC_dip = data1
data_CONC_dip.add(data2)
data_CONC_dip.add(data3)
data_CONC_dip.add(data4)

#inverting

ert_conc_dip =pb.ERTManager(data_CONC_dip)

inv_conc_dip =ert_conc_dip.invert(quality=quality, maxCellArea=maxCellArea, robustData=robustData, lam=lam,paraDX=paraDX)
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

#importing data

data1=pb.importData('18X4 1SA_5187.bin')
data2=pb.importData('18X4 1SB_5295.bin')
data3=pb.importData('18X4 1SC_5639.bin')
data4=pb.importData('18X4 1SD_5747.bin')

#MODIFY ELECTRODE POSITION

for i in range(data2.sensorCount()):
    data2.setSensorPosition(i, [i*4+1, 0, 0])
    
for i in range(data3.sensorCount()):
    data3.setSensorPosition(i, [i*4+2, 0, 0])
    
for i in range(data4.sensorCount()):
    data4.setSensorPosition(i, [i*4+3, 0, 0])


#calculating rhoa

data1.set('k', pb.geometricFactors(data1))
data1.set('rhoa', data1('u') / data1('i') * data1('k'))

data2.set('k', pb.geometricFactors(data2))
data2.set('rhoa', data2('u') / data2('i') * data2('k'))

data3.set('k', pb.geometricFactors(data3))
data3.set('rhoa', data3('u') / data3('i') * data3('k'))

data4.set('k', pb.geometricFactors(data4))
data4.set('rhoa', data4('u') / data4('i') * data4('k'))

#concatenating data

data_conc_sch = data1
data_conc_sch.add(data2)
data_conc_sch.add(data3)
data_conc_sch.add(data4)

#inverting

ert_conc_sch =pb.ERTManager(data_conc_sch)
inv_conc_sch =ert_conc_sch.invert(quality=quality, maxCellArea=maxCellArea, robustData=robustData, lam=lam,paraDX=paraDX)
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

data_72x1_dip=pb.importData('72x1_DIP.bin')


data_72x1_dip.set('k', pb.geometricFactors(data_72x1_dip))
data_72x1_dip.set('rhoa', data_72x1_dip('u') / data_72x1_dip('i') * data_72x1_dip('k'))

#inversion

ert_72x1_dip =pb.ERTManager(data_72x1_dip)

inv_72x1_dip =ert_72x1_dip.invert(quality=quality, maxCellArea=maxCellArea, robustData=robustData, lam=50,paraDX=paraDX)

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

data_72X1_sch=pb.importData('72X1SCH.bin')

data_72X1_sch.set('k', pb.geometricFactors(data_72X1_sch))
data_72X1_sch.set('rhoa', data_72X1_sch('u') / data_72X1_sch('i') * data_72X1_sch('k'))

#inversion

ert_72x1_sch =pb.ERTManager(data_72X1_sch)

inv_72x1_sch =ert_72x1_sch.invert(quality=quality, maxCellArea=maxCellArea, robustData=robustData, lam=50,paraDX=paraDX)

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
stats=open('mitjes_martorana_ALL.txt','w')
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
stats=open('mitjes_martorana_SPECIFIC.txt','w')
stats.write('la mitjana del perfil simple de 72 (dipol) és:'+str(mean_72x1_dd_s)+'\n')
stats.write('la mitjana del perfil simple de 72 (schlumb) és:'+str(mean_72x1_sch_s)+'\n')
stats.write('la mitjana del perfil simple de 18 (dipol) és:'+str(mean_18x4_dd_s)+'\n')
stats.write('la mitjana del perfil simple de 18 (schlumb) és:'+str(mean_18x4_sch_s)+'\n')
stats.write('la mitjana del perfil concatenat (dipol) és:'+str(mean_18x1_dd_s)+'\n')
stats.write('la mitjana del perfil concatenat (schlum) és:'+str(mean_18x1_sch_s)+'\n')
stats.close()
#FIGURES

#figura inversió

fig, (ax) = plt.subplots(ncols=2,nrows=4,figsize=(20,20),dpi=300)
plt.subplots_adjust(hspace=1,wspace=1)
pg.show(world,fillRegion=False,ax=ax[0,0])
ert_18x4_dip.showResult(cMin=15,cMax=125,ax=ax[0,0])
pg.show(world,fillRegion=False,ax=ax[1,0])
ert_18x4_sch.showResult(cMin=15,cMax=125,ax=ax[1,0])
pg.show(world,fillRegion=False,ax=ax[2,0])
ert_conc_dip.showResult(cMin=15,cMax=125,ax=ax[2,0])
pg.show(world,fillRegion=False,ax=ax[3,0])
ert_conc_sch.showResult(cMin=15,cMax=125,ax=ax[3,0])
pg.show(world,fillRegion=False,ax=ax[0,1])
ert_72x1_dip.showResult(cMin=15,cMax=125,ax=ax[0,1])
pg.show(world,fillRegion=False,ax=ax[1,1])
ert_72x1_sch.showResult(cMin=15,cMax=125,ax=ax[1,1])
pg.show(world,ax=ax[2,1])
pg.show(mesh,rhomap,label='Resistivity $(\Omega$m)',ax=ax[3,1])
ax[0,0].set(xlim=(20,50),ylim=(-10,0))
ax[0,0].set_title('18x4_dd')
ax[0,0].set_ylabel('Depth')
ax[0,0].set_xlabel('distance (m)')
ax[1,0].set(xlim=(20,50),ylim=(-10,0))
ax[1,0].set_title('18x4_slm')
ax[1,0].set_ylabel('Depth')
ax[1,0].set_xlabel('distance (m)')
ax[2,0].set(xlim=(20,50),ylim=(-10,0))
ax[2,0].set_title('18x1_dd_conc')
ax[2,0].set_ylabel('Depth')
ax[2,0].set_xlabel('distance (m)')
ax[3,0].set(xlim=(20,50),ylim=(-10,0))
ax[3,0].set_title('18x1_slm_conc')
ax[3,0].set_ylabel('Depth')
ax[3,0].set_xlabel('distance (m)')
ax[0,1].set(xlim=(20,50),ylim=(-10,0))
ax[0,1].set_title('72x1_dd')
ax[0,1].set_ylabel('Depth')
ax[0,1].set_xlabel('distance (m)')
ax[1,1].set(xlim=(20,50),ylim=(-10,0))
ax[1,1].set_title('72x1_slm')
ax[1,1].set_ylabel('Depth')
ax[1,1].set_xlabel('distance (m)')
ax[2,1].set(xlim=(20,50),ylim=(-10,0))
ax[2,1].set_ylabel('Depth')
ax[2,1].set_xlabel('distance (m)')
ax[3,1].set(xlim=(20,50),ylim=(-10,0))
ax[3,1].set_title('model')
ax[3,1].set_ylabel('Depth')
ax[3,1].set_xlabel('distance (m)')

plt.savefig('castellbisbal_invb.png')


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
ax[0,0].text(40, -4,str(round(mean_18x4_dip,3)),verticalalignment='bottom', horizontalalignment='right',color='red',fontsize=30)
ax[1,0].set_title('18x4_slm')
ax[1,0].set_ylabel('Depth')
ax[1,0].set_xlabel('distance (m)')
ax[1,0].text(40, -4,str(round(mean_18x4_sch,3)),verticalalignment='bottom', horizontalalignment='right',color='red',fontsize=30)
ax[2,0].set_title('18x1_dd_conc')
ax[2,0].set_ylabel('Depth')
ax[2,0].set_xlabel('distance (m)')
ax[2,0].text(40, -4,str(round(mean_conc_dip,3)),verticalalignment='bottom', horizontalalignment='right',color='red',fontsize=30)
ax[3,0].set_title('18x1_slm_conc')
ax[3,0].set_ylabel('Depth')
ax[3,0].set_xlabel('distance (m)')
ax[3,0].text(40, -4,str(round(mean_conc_sch,3)),verticalalignment='bottom', horizontalalignment='right',color='red',fontsize=30)
ax[0,1].set_title('72x1_dd')
ax[0,1].set_ylabel('Depth')
ax[0,1].set_xlabel('distance (m)')
ax[0,1].text(40, -4,str(round(mean_72x1_dip,3)),verticalalignment='bottom', horizontalalignment='right',color='red',fontsize=30)
ax[1,1].set_title('72x1_slm')
ax[1,1].set_ylabel('Depth')
ax[1,1].set_xlabel('distance (m)')
ax[1,1].text(40, -4,str(round(mean_72x1_sch,3)),verticalalignment='bottom', horizontalalignment='right',color='red',fontsize=30)
ax[2,1].set_title('18x1_dd_conc')
ax[2,1].set_ylabel('Depth')
ax[2,1].set_xlabel('distance (m)')
ax[2,1].text(40, -4,str(round(mean_conc_dip,3)),verticalalignment='bottom', horizontalalignment='right',color='red',fontsize=30)
ax[3,1].set_title('18x1_slm_conc')
ax[3,1].set_ylabel('Depth')
ax[3,1].set_xlabel('distance (m)')
ax[3,1].text(40, -4,str(round(mean_conc_sch,3)),verticalalignment='bottom', horizontalalignment='right',color='red',fontsize=30)
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
ax[0,0].text(35.7, -2.9,str(round(mean_18x4_dip2,3)),verticalalignment='bottom', horizontalalignment='right',color='white',fontsize=30)
ax[1,0].set_title('18x4_slm')
ax[1,0].set_ylabel('Depth')
ax[1,0].set_xlabel('distance (m)')
ax[1,0].text(35.7, -2.9,str(round(mean_18x4_sch2,3)),verticalalignment='bottom', horizontalalignment='right',color='white',fontsize=30)
ax[2,0].set_title('18x1_dd_conc')
ax[2,0].set_ylabel('Depth')
ax[2,0].set_xlabel('distance (m)')
ax[2,0].text(35.7, -2.9,str(round(mean_conc_dip2,3)),verticalalignment='bottom', horizontalalignment='right',color='white',fontsize=30)
ax[3,0].set_title('18x1_slm_conc')
ax[3,0].set_ylabel('Depth')
ax[3,0].set_xlabel('distance (m)')
ax[3,0].text(35.7, -2.9,str(round(mean_conc_sch2,3)),verticalalignment='bottom', horizontalalignment='right',color='white',fontsize=30)
ax[0,1].set_title('72x1_dd')
ax[0,1].set_ylabel('Depth')
ax[0,1].set_xlabel('distance (m)')
ax[0,1].text(35.7, -2.9,str(round(mean_72x1_dip2,3)),verticalalignment='bottom', horizontalalignment='right',color='white',fontsize=30)
ax[1,1].set_title('72x1_slm')
ax[1,1].set_ylabel('Depth')
ax[1,1].set_xlabel('distance (m)')
ax[1,1].text(35.7, -2.9,str(round(mean_72x1_sch2,3)),verticalalignment='bottom', horizontalalignment='right',color='white',fontsize=30)
ax[2,1].set_title('18x1_dd_conc')
ax[2,1].set_ylabel('Depth')
ax[2,1].set_xlabel('distance (m)')
ax[2,1].text(35.7, -2.9,str(round(mean_conc_dip2,3)),verticalalignment='bottom', horizontalalignment='right',color='white',fontsize=30)
ax[3,1].set_title('18x1_slm_conc')
ax[3,1].set_ylabel('Depth')
ax[3,1].set_xlabel('distance (m)')
ax[3,1].text(35.7, -2.9,str(round(mean_conc_sch2,3)),verticalalignment='bottom', horizontalalignment='right',color='white',fontsize=30)
plt.savefig('martorana_4x4_field_smallb.png')

#figura martorana inversió

cmin=15
cmax=125

fig, (ax) = plt.subplots(ncols=2,nrows=3,figsize=(20,15),dpi=300)
#plt.subplots_adjust(hspace=-2)
pg.show(world,fillRegion=False,ax=ax[0,0])
ert_18x4_dip.showResult(cMin=cmin,cMax=cmax,colorBar=False,ax=ax[0,0])
pg.show(world,fillRegion=False,ax=ax[1,0])
ert_conc_dip.showResult(cMin=cmin,cMax=cmax,colorBar=False,ax=ax[1,0])
pg.show(world,fillRegion=False,ax=ax[2,0])
ert_72x1_dip.showResult(cMin=cmin,cMax=cmax,colorBar=True,ax=ax[2,0])
pg.show(world,fillRegion=False,ax=ax[0,1])
pg.show(grid,martorana_18x4_dip,cMap='bwr',cMin=-1,cMax=1,ax=ax[0,1])
pg.show(world,fillRegion=False,ax=ax[1,1])
pg.show(grid,martorana_conc_dip,cMap='bwr',cMin=-1,cMax=1,ax=ax[1,1])
pg.show(world,fillRegion=False,ax=ax[2,1])
pg.show(grid,martorana_72x1_dip,cMap='bwr',label='Adj. index',cMin=-1,cMax=1,ax=ax[2,1])
ax[0,0].set(ylim=(-8,0))
ax[0,0].set(xlim=(25,45))
#ax[0,0].set_title('18x4_dd')
ax[0,0].set_ylabel('Depth')
ax[0,0].set_xlabel('distance (m)')
ax[1,0].set(ylim=(-8,0))
ax[1,0].set(xlim=(25,45))
#ax[1,0].set_title('conc_dd')
ax[1,0].set_ylabel('Depth')
ax[1,0].set_xlabel('distance (m)')
ax[2,0].set(ylim=(-8,0))
ax[2,0].set(xlim=(25,45))
#ax[2,0].set_title('72x1_dd')
ax[2,0].set_ylabel('Depth')
ax[2,0].set_xlabel('distance (m)')
ax[0,1].text(34.2, -2.5,str(round(mean_18x4_dip2,3)),verticalalignment='bottom',size=30,color='white')
#ax[0,1].text(35, -3,str(round(mean_18x4_dip,3)),verticalalignment='bottom',size=30,color='black')
#ax[0,1].set_title('martorana_18x1')
#ax[0,1].set_ylabel('Depth')
ax[0,1].set_xlabel('distance(m)')
ax[1,1].text(34.2, -2.5,str(round(mean_conc_dip2,3)),verticalalignment='bottom',size=30,color='white')
#ax[1,1].text(35, -3,str(round(mean_conc_dip,3)),verticalalignment='bottom',size=30,color='black')
#ax[1,1].set_title('martorana_72')
#ax[1,1].set_ylabel('Depth')
ax[1,1].set_xlabel('distance (m)')
ax[2,1].text(34.2, -2.5,str(round(mean_72x1_dip2,3)),verticalalignment='bottom',size=30,color='white')
#ax[2,1].text(35, -3,str(round(mean_72x1_dip,3)),verticalalignment='bottom',size=30,color='black')
#ax[2,1].set_ylabel('Depth')
ax[2,1].set_xlabel('distance (m)')

plt.savefig('inversion_martorana_dd.png')

fig, (ax) = plt.subplots(ncols=2,nrows=3,figsize=(45,20),dpi=300)
#plt.subplots_adjust(hspace=-2)
pg.show(world,fillRegion=False,ax=ax[0,0])
ert_18x4_sch.showResult(cMin=cmin,cMax=cmax,colorBar=False,ax=ax[0,0])
pg.show(world,fillRegion=False,ax=ax[1,0])
ert_conc_sch.showResult(cMin=cmin,cMax=cmax,colorBar=False,ax=ax[1,0])
pg.show(world,fillRegion=False,ax=ax[2,0])
ert_72x1_sch.showResult(cMin=cmin,cMax=cmax,colorBar=True,ax=ax[2,0])
pg.show(world,fillRegion=False,ax=ax[0,1])
pg.show(grid,martorana_18x4_dip,cMap='bwr',cMin=-1,cMax=1,ax=ax[0,1])
pg.show(world,fillRegion=False,ax=ax[1,1])
pg.show(grid,martorana_conc_dip,cMap='bwr',cMin=-1,cMax=1,ax=ax[1,1])
pg.show(world,fillRegion=False,ax=ax[2,1])
pg.show(grid,martorana_72x1_dip,cMap='bwr',label='Adj. index',cMin=-1,cMax=1,ax=ax[2,1])
ax[0,0].set(ylim=(-8,0))
ax[0,0].set(xlim=(20,50))
ax[0,0].set_title('18x4_dd')
ax[0,0].set_ylabel('Depth')
ax[0,0].set_xlabel('distance (m)')
ax[1,0].set(ylim=(-8,0))
ax[1,0].set(xlim=(20,50))
ax[1,0].set_title('conc_dd')
ax[1,0].set_ylabel('Depth')
ax[1,0].set_xlabel('distance (m)')
ax[2,0].set(ylim=(-8,0))
ax[2,0].set(xlim=(20,50))
ax[2,0].set_title('72x1_dd')
ax[2,0].set_ylabel('Depth')
ax[2,0].set_xlabel('distance (m)')
ax[0,1].text(35, -2,str(round(mean_18x4_sch2,3)),verticalalignment='bottom',size=30,color='white')
ax[0,1].text(35, -3,str(round(mean_18x4_sch,3)),verticalalignment='bottom',size=30,color='black')
ax[0,1].set_title('martorana_18x1')
ax[0,1].set_ylabel('Depth')
ax[0,1].set_xlabel('martorana_conc')
ax[1,1].text(35, -2,str(round(mean_conc_sch2,3)),verticalalignment='bottom',size=30,color='white')
ax[1,1].text(35, -3,str(round(mean_conc_sch,3)),verticalalignment='bottom',size=30,color='black')
ax[1,1].set_title('martorana_72')
ax[1,1].set_ylabel('Depth')
ax[1,1].set_xlabel('distance (m)')
ax[2,1].text(35, -2,str(round(mean_72x1_sch2,3)),verticalalignment='bottom',size=30,color='white')
ax[2,1].text(35, -3,str(round(mean_72x1_sch,3)),verticalalignment='bottom',size=30,color='black')
ax[2,1].set_ylabel('Depth')
ax[2,1].set_xlabel('distance (m)')

plt.savefig('inversion_martorana_sch.png')

#figura martorana vs martorana

fig, (ax) = plt.subplots(ncols=2,nrows=2, figsize=(20,20),dpi=300)
pg.show(world,fillRegion=False,ax=ax[0,0])
pg.show(grid,martorana_18x4vs_conc_dip,cMap='bwr',cMin=-2,cMax=2,ax=ax[0,0])
pg.show(world,fillRegion=False,ax=ax[1,0])
pg.show(grid,martorana_18x4vs_conc_sch,cMap='bwr',cMin=-2,cMax=2,ax=ax[1,0])
pg.show(world,fillRegion=False,ax=ax[0,1])
pg.show(grid,martorana_conc_vs72x1_dip,cMap='bwr',cMin=-2,cMax=2,ax=ax[0,1])
pg.show(world,fillRegion=False,ax=ax[1,1])
pg.show(grid,martorana_conc_vs72x1_sch,cMap='bwr',cMin=-2,cMax=2,ax=ax[1,1])
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

fig, (ax) = plt.subplots(figsize=(15,5),dpi=300)
pg.show(mesh,rhomap,showMesh=True,label='Resistivity $(\Omega$m)',cMap='Blues',ax=ax)
ax.set(xlim=(20,50),ylim=(-6,0))
#ax.set_title('model')
ax.set_ylabel('Depth (m)')
ax.set_xlabel('distance (m)')
plt.savefig('model_blue.png')

ert_conc_dip.saveResult(folder='results18x1dd')
ert_conc_sch.saveResult(folder='results18x1slm')
ert_18x4_dip.saveResult(folder='results18x4dd')
ert_18x4_sch.saveResult(folder='results18x4slm')
ert_72x1_dip.saveResult(folder='results72x1dd')
ert_72x1_sch.saveResult(folder='results72x1slm')
