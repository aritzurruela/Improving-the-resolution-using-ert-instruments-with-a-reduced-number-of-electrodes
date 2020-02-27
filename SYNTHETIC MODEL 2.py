# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 15:27:10 2019

@author: Aritz
"""
import numpy as np
import pygimli as pg
import pygimli.meshtools as mt
import matplotlib.pyplot as plt
import pybert as pb

#creación grid para comparación

grid = pg.createGrid(x=pg.utils.grange(start=0, end=71, n=500),
                           y=-pg.utils.grange(0, 15, n=250, log=False))

grid2= pg.createGrid(x=pg.utils.grange(start=50, end=60, n=500),
                           y=-pg.utils.grange(5.5, 9.5
                                              , n=250, log=False))
square=mt.createRectangle(start=[50,-5.5],end=[60,-9.5])

#parametros inversion

quality=34.5
maxCellArea=1
robustData=True
paraDX=0.25
lam=20

#PARAMETROS MODELADO

noise=0
rhomap= [[0,500],[1,50],[2,100],[3,2000]]
#CREACION MODELO GENERAL

background=mt.createWorld(start=[-10,0],end=[81,-20], area=1,marker=0)

pol0=mt.createPolygon([[-10,-2.5],[27.5,-2.5],[30,0],[-10,0]], isClosed=True,marker=1)
pol1=mt.createPolygon([[-10,-20],[10,-20],[27.5,-2.5],[-10,-2.5]], isClosed=True,marker=2)
pol2=mt.createPolygon([[60,0],[65,0],[45,-20],[40,-20]],isClosed=True, marker=3)

world = mt.mergePLC([background,pol0,pol1,pol2])
mesh= mt.createMesh(world, quality=33, area=0.1,smooth=[1,2])

#EXPORTACION DATOS DE RESISTIVIDAD DEL MODELO (VECTOR)

res_model=pg.solver.parseArgToArray(rhomap, mesh.cellCount(), mesh)

pg.show(mesh,res_model)

#INTERPOLACION DATOS MODELO A MESH INVERSION 

nan=99.9
model_vector = []
for pos in grid.cellCenters():
    cell = mesh.findCell(pos)
    if cell:
        model_vector.append(res_model[cell.id()])
    else:
        model_vector.append(nan)

nan=99.9
model_vector_s = []
for pos in grid2.cellCenters():
    cell = mesh.findCell(pos)
    if cell:
        model_vector_s.append(res_model[cell.id()])
    else:
        model_vector_s.append(nan)

model_data=pb.pg.RVector(model_vector)

model_data_s=pb.pg.RVector(model_vector_s)

#72x1 DD -------------------------------------------------------------

world_72_dd=world

scheme72_dd = pb.createData(elecs=pg.utils.grange(start=0, end=71, n=72),
                       schemeName='dd')
for pos in scheme72_dd.sensorPositions():
    world_72_dd.createNode(pos)
    world_72_dd.createNode(pos + pg.RVector3(0, -0.1))


mesh_72x1_dd= mt.createMesh(world_72_dd, quality=34)

#forward modelling

data_72x1_dd=pb.simulate(mesh_72x1_dd,res=rhomap,scheme=scheme72_dd,noise=noise, verbose=True)
#inversion

print('starting 72x1dd inversion')

ert_72x1_dd=pb.ERTManager(data_72x1_dd)

inversion_72x1_dd=ert_72x1_dd.invert(quality=quality, maxCellArea=0.3, robustData=robustData, lam=lam,paraDX=paraDX,verbose=True)

inversion_mesh_72x1_dd=ert_72x1_dd.paraDomain

print('interpolation of inversion data to grid')

#interpolation of inversion data to grid

nan=99.9
inversion_vector_72x1_dd = []
for pos in grid.cellCenters():
    cell = inversion_mesh_72x1_dd.findCell(pos)
    if cell:
        inversion_vector_72x1_dd.append(inversion_72x1_dd[cell.id()])
    else:
        inversion_vector_72x1_dd.append(nan)

nan=99.9
inversion_vector_72x1_dd_s = []
for pos in grid2.cellCenters():
    cell = inversion_mesh_72x1_dd.findCell(pos)
    if cell:
        inversion_vector_72x1_dd_s.append(inversion_72x1_dd[cell.id()])
    else:
        inversion_vector_72x1_dd_s.append(nan)

inversion_data_72x1_dd=pb.pg.RVector(inversion_vector_72x1_dd)

martorana_sim_72x1_dd=np.log10(inversion_data_72x1_dd/model_data)

inversion_data_72x1_dd_s=pb.pg.RVector(inversion_vector_72x1_dd_s)

martorana_sim_72x1_dd_s=np.log10(inversion_data_72x1_dd_s/model_data_s)

#72x1 SLM--------------------------------------------------------------------------------

world_72_sch=world

scheme72_slm = pb.createData(elecs=pg.utils.grange(start=0, end=71, n=72),
                       schemeName='slm')
for pos in scheme72_slm.sensorPositions():
    world_72_sch.createNode(pos)
    world_72_sch.createNode(pos + pg.RVector3(0, -0.1))
    
mesh_72x1_slm= mt.createMesh(world_72_sch)
    
#forward modelling

data_72x1_slm=pb.simulate(mesh_72x1_slm,res=rhomap,scheme=scheme72_slm,noiseLevel= noise, verbose=True)

#inversion

print('starting 72x1 sch inversion')

ert_72x1_slm=pb.ERTManager(data_72x1_slm)

inversion_72x1_slm=ert_72x1_slm.invert(quality=quality, maxCellArea=0.3, robustData=robustData, lam=lam,paraDX=paraDX)

inversion_mesh_72x1_slm=ert_72x1_slm.paraDomain

#interpolation of inversion data to grid

print('interpolation of inversion data to grid')

nan=99.9
inversion_vector_72x1_slm = []
for pos in grid.cellCenters():
    cell = inversion_mesh_72x1_slm.findCell(pos)
    if cell:
        inversion_vector_72x1_slm.append(inversion_72x1_slm[cell.id()])
    else:
        inversion_vector_72x1_slm.append(nan)

nan=99.9
inversion_vector_72x1_slm_s = []
for pos in grid2.cellCenters():
    cell = inversion_mesh_72x1_slm.findCell(pos)
    if cell:
        inversion_vector_72x1_slm_s.append(inversion_72x1_slm[cell.id()])
    else:
        inversion_vector_72x1_slm_s.append(nan)
        
inversion_data_72x1_slm=pb.pg.RVector(inversion_vector_72x1_slm)

martorana_sim_72x1_slm=np.log10(inversion_data_72x1_slm/model_data)

inversion_data_72x1_slm_s=pb.pg.RVector(inversion_vector_72x1_slm_s)

martorana_sim_72x1_slm_s=np.log10(inversion_data_72x1_slm_s/model_data_s)

#------------------------------------------------------------------------------------------------------

#18X4 DD

world_18x4_dip=world

scheme_18x4_dd = pb.createData(elecs=pg.utils.grange(start=0, end=71, n=18),
                       schemeName='dd')
for pos in scheme_18x4_dd.sensorPositions():
    world_18x4_dip.createNode(pos)
    world_18x4_dip.createNode(pos + pg.RVector3(0, -0.1))

   
mesh_18x4_dip=mt.createMesh(world_18x4_dip)
    
data_18x4_dd=pb.simulate(mesh_18x4_dip,res=rhomap,scheme=scheme_18x4_dd,noiseLevel= noise, verbose=True)

#18X4 SLM

world_18x4_sch=world

scheme_18x4_sch = pb.createData(elecs=pg.utils.grange(start=0, end=71, n=18),
                       schemeName='slm')
for pos in scheme_18x4_sch.sensorPositions():
    world_18x4_sch.createNode(pos)
    world_18x4_sch.createNode(pos + pg.RVector3(0, -0.1))
    
mesh_18x4_sch=mt.createMesh(world_18x4_sch)
    
data_18x4_slm=pb.simulate(mesh_18x4_sch,res=rhomap,scheme=scheme_18x4_sch,noiseLevel= noise, verbose=True)

#INVERSION
print('starting inversion of 18x4dd')
#18x4_dd

ert_18x4_dd=pb.ERTManager(data_18x4_dd)

inversion_18x4_dd=ert_18x4_dd.invert(quality=quality, maxCellArea=maxCellArea, robustData=robustData, lam=lam,paraDX=paraDX)

inversion_mesh_18x4_dd=ert_18x4_dd.paraDomain

#18x4_slm

print('starting inversion of 18x4sch')

ert_18x4_slm=pb.ERTManager(data_18x4_slm)

inversion_18x4_slm=ert_18x4_slm.invert(quality=quality, maxCellArea=maxCellArea, robustData=robustData, lam=lam,paraDX=paraDX)

inversion_mesh_18x4_slm=ert_18x4_slm.paraDomain
       
#pg.show(grid,model_vector_nn)  

#INTERPOLACION DATOS inversion A grid 18x4 dd

nan=99.9
inversion_vector_18x4_dd = []
for pos in grid.cellCenters():
    cell = inversion_mesh_18x4_dd.findCell(pos)
    if cell:
        inversion_vector_18x4_dd.append(inversion_18x4_dd[cell.id()])
    else:
        inversion_vector_18x4_dd.append(nan)

nan=99.9
inversion_vector_18x4_dd_s = []
for pos in grid2.cellCenters():
    cell = inversion_mesh_18x4_dd.findCell(pos)
    if cell:
        inversion_vector_18x4_dd_s.append(inversion_18x4_dd[cell.id()])
    else:
        inversion_vector_18x4_dd_s.append(nan)

inversion_data_18x4_dd=pb.pg.RVector(inversion_vector_18x4_dd)

inversion_data_18x4_dd_s=pb.pg.RVector(inversion_vector_18x4_dd_s)
    
#INTERPOLACION DATOS inversion A grid 18x4 slm

nan=99.9
inversion_vector_18x4_slm = []
for pos in grid.cellCenters():
    cell = inversion_mesh_18x4_slm.findCell(pos)
    if cell:
        inversion_vector_18x4_slm.append(inversion_18x4_slm[cell.id()])
    else:
        inversion_vector_18x4_slm.append(nan)

nan=99.9
inversion_vector_18x4_slm_s = []
for pos in grid2.cellCenters():
    cell = inversion_mesh_18x4_slm.findCell(pos)
    if cell:
        inversion_vector_18x4_slm_s.append(inversion_18x4_slm[cell.id()])
    else:
        inversion_vector_18x4_slm_s.append(nan)
        
inversion_data_18x4_slm=pb.pg.RVector(inversion_vector_18x4_slm)
inversion_data_18x4_slm_s=pb.pg.RVector(inversion_vector_18x4_slm_s)

#MARTORANA

martorana_sim_18x4_dd=np.log10(inversion_data_18x4_dd/model_data)

martorana_sim_18x4_slm=np.log10(inversion_data_18x4_slm/model_data)

martorana_sim_18x4_dd_s=np.log10(inversion_data_18x4_dd_s/model_data_s)

martorana_sim_18x4_slm_s=np.log10(inversion_data_18x4_slm_s/model_data_s)

#------------------------------------------------------------------------------------------

#18X1_dd

world1 = world

scheme1 = pb.createData(elecs=pg.utils.grange(start=0, end=71, n=18),
                       schemeName='dd')

for pos in scheme1.sensorPositions():
    world1.createNode(pos)
    world1.createNode(pos + pg.RVector3(0, -0.1))



    
mesh1= mt.createMesh(world1)

data1_dd=pb.simulate(mesh1,res=rhomap,scheme=scheme1,noiseLevel= noise, verbose=True)

#18x2_dd

world2 =world

scheme2 = pb.createData(elecs=pg.utils.grange(start=1, end=71, n=18),
                       schemeName='dd')

for pos in scheme2.sensorPositions():
    world2.createNode(pos)
    world2.createNode(pos + pg.RVector3(0, -0.1))

    
mesh2= mt.createMesh(world2)

data2_dd=pb.simulate(mesh2,res=rhomap,scheme=scheme2,noiseLevel= noise, verbose=True)

#18x3_dd

world3 = world

scheme3 = pb.createData(elecs=pg.utils.grange(start=1, end=71, n=18),
                       schemeName='dd')

for pos in scheme3.sensorPositions():
    world3.createNode(pos)
    world3.createNode(pos + pg.RVector3(0, -0.1))
    
mesh3= mt.createMesh(world3)

data3_dd=pb.simulate(mesh3,res=rhomap,scheme=scheme3,noiseLevel= noise, verbose=True)

#18x4_dd

world4 = world

scheme4 = pb.createData(elecs=pg.utils.grange(start=1, end=71, n=18),
                       schemeName='dd')

for pos in scheme4.sensorPositions():
    world4.createNode(pos)
    world4.createNode(pos + pg.RVector3(0, -0.1))
    
mesh4= mt.createMesh(world4)

data4_dd=pb.simulate(mesh4,res=rhomap,scheme=scheme4,noiseLevel= noise, verbose=True)

#CONCATENACION 18x1 dd

data_total_dd = data1_dd

data_total_dd.add(data2_dd)
data_total_dd.add(data3_dd)
data_total_dd.add(data4_dd)

#INVERSION 18x1 dd

ert_18x1_dd=pb.ERTManager(data_total_dd)

inversion_18x1_dd=ert_18x1_dd.invert(quality=quality, maxCellArea=maxCellArea, robustData=robustData, lam=lam,paraDX=paraDX)

inversion_mesh_18x1_dd=ert_18x1_dd.paraDomain

#INTERPOLACION DATOS INVERSION A MESH INVERSION 18x1 dd

nan=99.9
inversion_vector_18x1_dd = []
for pos in grid.cellCenters():
    cell = inversion_mesh_18x1_dd.findCell(pos)
    if cell:
        inversion_vector_18x1_dd.append(inversion_18x1_dd[cell.id()])
    else:
        inversion_vector_18x1_dd.append(nan)
        
nan=99.9
inversion_vector_18x1_dd_s = []
for pos in grid2.cellCenters():
    cell = inversion_mesh_18x1_dd.findCell(pos)
    if cell:
        inversion_vector_18x1_dd_s.append(inversion_18x1_dd[cell.id()])
    else:
        inversion_vector_18x1_dd_s.append(nan)

inversion_data_18x1_dd=pb.pg.RVector(inversion_vector_18x1_dd)

inversion_data_18x1_dd_s=pb.pg.RVector(inversion_vector_18x1_dd_s)

#MARTORANA

martorana_sim_18x1_dd=np.log10(inversion_data_18x1_dd/model_data)

martorana_sim_18x1_dd_s=np.log10(inversion_data_18x1_dd_s/model_data_s)

#-------------------------------------------------------------------------------------       
#18X1_slm

world1 = world

scheme1 = pb.createData(elecs=pg.utils.grange(start=0, end=71, n=18),
                       schemeName='slm')

for pos in scheme1.sensorPositions():
    world1.createNode(pos)
    world1.createNode(pos + pg.RVector3(0, -0.1))

    
mesh1= mt.createMesh(world1)

data1_slm=pb.simulate(mesh1,res=rhomap,scheme=scheme1,noiseLevel= noise, verbose=True)

#18x2_slm

world2 =world

scheme2 = pb.createData(elecs=pg.utils.grange(start=1, end=71, n=18),
                       schemeName='slm')

for pos in scheme2.sensorPositions():
    world2.createNode(pos)
    world2.createNode(pos + pg.RVector3(0, -0.1))

    
mesh2= mt.createMesh(world2)

data2_slm=pb.simulate(mesh2,res=rhomap,scheme=scheme2,noiseLevel= noise, verbose=True)

#18x3_slm

world3 = world

scheme3 = pb.createData(elecs=pg.utils.grange(start=2, end=71, n=18),
                       schemeName='slm')

for pos in scheme3.sensorPositions():
    world3.createNode(pos)
    world3.createNode(pos + pg.RVector3(0, -0.1))
    
mesh3= mt.createMesh(world3)

data3_slm=pb.simulate(mesh3,res=rhomap,scheme=scheme3,noiseLevel= noise, verbose=True)

#18x4_slm

world4 = world

scheme4 = pb.createData(elecs=pg.utils.grange(start=3, end=71, n=18),
                       schemeName='slm')

for pos in scheme4.sensorPositions():
    world4.createNode(pos)
    world4.createNode(pos + pg.RVector3(0, -0.1))
    
mesh4= mt.createMesh(world4)

data4_slm=pb.simulate(mesh4,res=rhomap,scheme=scheme4,noiseLevel= noise, verbose=True)

#CONCATENACION 18x1 slm

data_total_slm = data1_slm

data_total_slm.add(data2_slm)
data_total_slm.add(data3_slm)
data_total_slm.add(data4_slm)

#INVERSION 18x1 slm

ert_18x1_slm=pb.ERTManager(data_total_slm)

inversion_18x1_slm=ert_18x1_slm.invert(quality=quality, maxCellArea=maxCellArea, robustData=robustData, lam=lam,paraDX=paraDX)

inversion_mesh_18x1_slm=ert_18x1_slm.paraDomain

#INTERPOLACION DATOS INVERSION A MESH INVERSION 18x1 slm

nan=99.9
inversion_vector_18x1_slm = []
for pos in grid.cellCenters():
    cell = inversion_mesh_18x1_slm.findCell(pos)
    if cell:
        inversion_vector_18x1_slm.append(inversion_18x1_slm[cell.id()])
    else:
        inversion_vector_18x1_slm.append(nan)

nan=99.9
inversion_vector_18x1_slm_s = []
for pos in grid2.cellCenters():
    cell = inversion_mesh_18x1_slm.findCell(pos)
    if cell:
        inversion_vector_18x1_slm_s.append(inversion_18x1_slm[cell.id()])
    else:
        inversion_vector_18x1_slm_s.append(nan)

inversion_data_18x1_slm=pb.pg.RVector(inversion_vector_18x1_slm)
inversion_data_18x1_slm_s=pb.pg.RVector(inversion_vector_18x1_slm_s)

#MARTORANA

martorana_sim_18x1_slm=np.log10(inversion_data_18x1_slm/model_data)

martorana_sim_18x1_slm_s=np.log10(inversion_data_18x1_slm_s/model_data_s)
#------------------------------------------

# valor absoluto
martorana_sim_18x4_dd_abs=np.absolute(martorana_sim_18x4_dd)
martorana_sim_18x4_slm_abs=np.absolute(martorana_sim_18x4_slm)
martorana_sim_18x1_dd_abs=np.absolute(martorana_sim_18x1_dd)
martorana_sim_18x1_slm_abs=np.absolute(martorana_sim_18x1_slm)
martorana_sim_72x1_dd_abs=np.absolute(martorana_sim_72x1_dd)
martorana_sim_72x1_slm_abs=np.absolute(martorana_sim_72x1_slm)

martorana_sim_18x4_dd_s_abs=np.absolute(martorana_sim_18x4_dd_s)
martorana_sim_18x4_slm_s_abs=np.absolute(martorana_sim_18x4_slm_s)
martorana_sim_18x1_dd_s_abs=np.absolute(martorana_sim_18x1_dd_s)
martorana_sim_18x1_slm_s_abs=np.absolute(martorana_sim_18x1_slm_s)
martorana_sim_72x1_dd_s_abs=np.absolute(martorana_sim_72x1_dd_s)
martorana_sim_72x1_slm_s_abs=np.absolute(martorana_sim_72x1_slm_s)

#stats

# min_18x4= min(martorana_sim_18x4_dd_abs)
# max_18x4=max(martorana_sim_18x4_dd_abs)
# std_18x4=np.std(martorana_sim_18x4_dd_abs)
# mean_18x4=np.mean(martorana_sim_18x4_dd_abs)
# n = len(martorana_sim_18x4_dd_abs)

# stats=open('stats_martorana_sim_18x4_dd_abs.txt','w')
# stats.write('número de valors és : ' + str(n)+ '\n')
# stats.write('la mitjana és:'+str(mean_18x4)+'\n')
# stats.write('valor mínim és : ' + str(min_18x4)+ '\n')
# stats.write('valor màxim és : ' + str(max_18x4)+ '\n')
# stats.write('error standard és : ' + str(std_18x4)+ '\n')
# stats.close()

#martorana vs martorana

martorana_sim_18x4vs18x1_dd=np.greater(martorana_sim_18x1_dd_abs,martorana_sim_18x4_dd_abs,dtype=object)
martorana_sim_18x4vs18x1_slm=np.greater(martorana_sim_18x1_slm_abs,martorana_sim_18x4_slm_abs,dtype=object)
martorana_sim_18x1vs72x1_dd=np.greater(martorana_sim_72x1_dd_abs,martorana_sim_18x1_dd_abs,dtype=object)
martorana_sim_18x1vs72x1_slm=np.greater(martorana_sim_72x1_slm_abs,martorana_sim_18x1_slm_abs,dtype=object)

#STATS


mean_72x1_dd=np.mean(martorana_sim_72x1_dd_abs)
mean_72x1_sch=np.mean(martorana_sim_72x1_slm_abs)
mean_18x1_dd=np.mean(martorana_sim_18x1_dd_abs)
mean_18x1_sch=np.mean(martorana_sim_18x1_slm_abs)
mean_18x4_dd=np.mean(martorana_sim_18x4_dd_abs)
mean_18x4_sch=np.mean(martorana_sim_18x4_slm_abs)
stats=open('mitjes_martorana_sim2_ref.txt','w')
stats.write('la mitjana del perfil simple de 72 (dipol) és:'+str(mean_72x1_dd)+'\n')
stats.write('la mitjana del perfil simple de 72 (schlumb) és:'+str(mean_72x1_sch)+'\n')
stats.write('la mitjana del perfil simple de 18 (dipol) és:'+str(mean_18x4_dd)+'\n')
stats.write('la mitjana del perfil simple de 18 (schlumb) és:'+str(mean_18x4_sch)+'\n')
stats.write('la mitjana del perfil concatenat (dipol) és:'+str(mean_18x1_dd)+'\n')
stats.write('la mitjana del perfil concatenat (schlum) és:'+str(mean_18x1_sch)+'\n')
stats.close()


mean_72x1_dd_s=np.mean(martorana_sim_72x1_dd_s_abs)
mean_72x1_sch_s=np.mean(martorana_sim_72x1_slm_s_abs)
mean_18x1_dd_s=np.mean(martorana_sim_18x1_dd_s_abs)
mean_18x1_sch_s=np.mean(martorana_sim_18x1_slm_s_abs)
mean_18x4_dd_s=np.mean(martorana_sim_18x4_dd_s_abs)
mean_18x4_sch_s=np.mean(martorana_sim_18x4_slm_s_abs)
stats=open('mitjes_martorana_sim_s_ref.txt','w')
stats.write('la mitjana del perfil simple de 72 (dipol) és:'+str(mean_72x1_dd_s)+'\n')
stats.write('la mitjana del perfil simple de 72 (schlumb) és:'+str(mean_72x1_sch_s)+'\n')
stats.write('la mitjana del perfil simple de 18 (dipol) és:'+str(mean_18x4_dd_s)+'\n')
stats.write('la mitjana del perfil simple de 18 (schlumb) és:'+str(mean_18x4_sch_s)+'\n')
stats.write('la mitjana del perfil concatenat (dipol) és:'+str(mean_18x1_dd_s)+'\n')
stats.write('la mitjana del perfil concatenat (schlum) és:'+str(mean_18x1_sch_s)+'\n')
stats.close()



#figures 

#boxplots

dd_all=[martorana_sim_18x4_dd,martorana_sim_18x1_dd,martorana_sim_72x1_dd]

dd=[martorana_sim_18x4_dd_s_abs,martorana_sim_18x1_dd_s_abs,martorana_sim_72x1_dd_s_abs]
sch=[martorana_sim_18x4_slm_s_abs,martorana_sim_18x1_slm_s_abs,martorana_sim_72x1_slm_s_abs]

dds=[martorana_sim_18x4_dd_s,martorana_sim_18x1_dd_s,martorana_sim_72x1_dd_s]
schs=[martorana_sim_18x4_slm_s,martorana_sim_18x1_slm_s,martorana_sim_72x1_slm_s]
mean_dd=(mean_18x4_dd_s,mean_18x1_dd_s,mean_72x1_dd_s)

mean_sch=(mean_18x4_sch_s,mean_18x1_sch_s,mean_72x1_sch_s)

names=('profile 18','concatenation','profile 72')

fig,ax=plt.subplots(figsize=(20,20),dpi=300)
ax.hist(martorana_sim_18x4_dd,alpha=0.5,color='g', label='18',bins=100)
ax.hist(martorana_sim_18x1_dd,alpha=0.5,color='b', label='conc',bins=100)
ax.hist(martorana_sim_72x1_dd,alpha=0.5,color='r', label='72',bins=100)
ax.legend(loc='upper right')
plt.savefig('histogram_dd.png')

fig,ax=plt.subplots(figsize=(20,20),dpi=300)
ax.hist(martorana_sim_18x4_slm,alpha=0.5,color='g', label='18',bins=100)
ax.hist(martorana_sim_18x1_slm,alpha=0.5,color='b', label='conc',bins=100)
ax.hist(martorana_sim_72x1_slm,alpha=0.5,color='r', label='72',bins=100)
ax.legend(loc='upper right')
plt.savefig('histogram_sch.png')

fig,ax=plt.subplots(figsize=(20,20),dpi=300)
ax.hist(martorana_sim_18x4_dd_s,alpha=0.5,color='g', label='18',bins=100)
ax.hist(martorana_sim_18x1_dd_s,alpha=0.5,color='b', label='conc',bins=100)
ax.hist(martorana_sim_72x1_dd_s,alpha=0.5,color='r', label='72',bins=100)
ax.legend(loc='upper right')
ax.set(xlim=(-1,-0.8))
plt.savefig('histogram_dds.png')

fig,ax=plt.subplots(figsize=(20,20),dpi=300)
ax.hist(martorana_sim_18x4_slm_s,alpha=0.5,color='g', label='18',bins=100)
ax.hist(martorana_sim_18x1_slm_s,alpha=0.5,color='b', label='conc',bins=100)
ax.hist(martorana_sim_72x1_slm_s,alpha=0.5,color='r', label='72',bins=100)
ax.legend(loc='upper right')
ax.set(xlim=(-1,-0.8))
plt.savefig('histogram_schs.png')

fig,ax=plt.subplots(figsize=(20,20),dpi=300)
ax.set_title( 'absolute dipole distribution')
ax.boxplot(dd)
plt.savefig('boxplot dipole.png')

fig,ax=plt.subplots(figsize=(20,20),dpi=300)
ax.set_title( 'absolute schlumberger distribution')
ax.boxplot(sch)
plt.savefig('boxplot schlumberger.png')

fig,ax=plt.subplots(figsize=(20,20),dpi=300)
ax.set_title( 'specific dipole distribution')
ax.boxplot(dds,showfliers=False)
plt.savefig('boxplot dipole_s.png')

fig,ax=plt.subplots(figsize=(20,20),dpi=300)
ax.set_title( 'specific schlumberger distribution ')
ax.boxplot(schs,showfliers=False)
plt.savefig('boxplot schlumberger_s.png')

fig,ax=plt.subplots(figsize=(20,20),dpi=300)
ax.bar(names,mean_dd)
ax.set_title('dipole absolute average')
plt.savefig('mean_dd.png')

fig,ax=plt.subplots(figsize=(20,20),dpi=300)
ax.bar(names,mean_sch)
ax.set_title('schlumberger absolute average')
plt.savefig('mean_sch.png')

cmin=200
cmax=1500


fig, (ax) = plt.subplots(nrows=3,figsize=(20,20),dpi=300)
#plt.subplots_adjust(hspace=1,wspace=1)
pg.show(world,fillRegion=False,ax=ax[0])
ert_18x4_dd.showResult(cMin=cmin,cMax=cmax,colorBar=False,ax=ax[0])
pg.show(world,fillRegion=False,ax=ax[1])
ert_18x1_dd.showResult(cMin=cmin,cMax=cmax,colorBar=False,ax=ax[1])
pg.show(world,fillRegion=False,ax=ax[2])
ert_72x1_dd.showResult(cMin=cmin,cMax=cmax,ax=ax[2])
ax[0].set(ylim=(-15,0))
ax[0].set(xlim=(0,71))
ax[0].set_title('18x4_dd')
ax[0].set_ylabel('Depth')
ax[0].set_xlabel('distance (m)')
ax[1].set(xlim=(0,71))
ax[1].set(ylim=(-15,0))
ax[1].set_title('18x1_dd_conc')
ax[1].set_ylabel('Depth')
ax[1].set_xlabel('distance (m)')
ax[2].set(xlim=(0,71))
ax[2].set(ylim=(-15,0))
ax[2].set_title('72x1_dd')
ax[2].set_ylabel('Depth')
ax[2].set_xlabel('distance (m)')

plt.savefig('inversion_dipolo_ref.png')

fig, (ax) = plt.subplots(nrows=3,figsize=(20,20),dpi=300)
#plt.subplots_adjust(hspace=1,wspace=1)
pg.show(world,fillRegion=False,ax=ax[0])
ert_18x4_slm.showResult(cMin=cmin,cMax=cmax,colorBar=False,ax=ax[0])
pg.show(world,fillRegion=False,ax=ax[1])
ert_18x1_slm.showResult(cMin=cmin,cMax=cmax,colorBar=False,ax=ax[1])
pg.show(world,fillRegion=False,ax=ax[2])
ert_72x1_slm.showResult(cMin=cmin,cMax=cmax,ax=ax[2])
ax[0].set(ylim=(-15,0))
ax[0].set(xlim=(0,71))
ax[0].set_title('18x4_dd')
ax[0].set_ylabel('Depth')
ax[0].set_xlabel('distance (m)')
ax[1].set(ylim=(-15,0))
ax[1].set(xlim=(0,71))
ax[1].set_title('18x1_dd_conc')
ax[1].set_ylabel('Depth')
ax[1].set_xlabel('distance (m)')
ax[2].set(ylim=(-15,0))
ax[2].set(xlim=(0,71))
ax[2].set_title('72x1_dd')
ax[2].set_ylabel('Depth')
ax[2].set_xlabel('distance (m)')

plt.savefig('inversion_sch_ref.png')


fig, (ax) = plt.subplots(ncols=2,nrows=3,figsize=(32,15),dpi=300)
#plt.subplots_adjust(hspace=-2)
pg.show(world,fillRegion=False,ax=ax[0,0])
ert_18x4_dd.showResult(cMin=cmin,cMax=cmax,colorBar=False,ax=ax[0,0])
pg.show(world,fillRegion=False,ax=ax[1,0])
ert_18x1_dd.showResult(cMin=cmin,cMax=cmax,colorBar=False,ax=ax[1,0])
pg.show(world,fillRegion=False,ax=ax[2,0])
ert_72x1_dd.showResult(cMin=cmin,cMax=cmax,colorBar=True,ax=ax[2,0])
pg.show(square,fillRegion=False,ax=ax[0,1])
pg.show(grid,martorana_sim_18x4_dd,cMap='bwr',cMin=-1,cMax=1,ax=ax[0,1])
pg.show(square,fillRegion=False,ax=ax[1,1])
pg.show(grid,martorana_sim_18x1_dd,cMap='bwr',cMin=-1,cMax=1,ax=ax[1,1])
pg.show(square,fillRegion=False,ax=ax[2,1])
pg.show(grid,martorana_sim_72x1_dd,cMap='bwr',label='Adj. index',cMin=-1,cMax=1,ax=ax[2,1])
ax[0,0].set(ylim=(-15,0))
ax[0,0].set(xlim=(0,71))
#ax[0,0].set_title('18x4_dd')
ax[0,0].set_ylabel('Depth')
ax[0,0].set_xlabel('distance (m)')
ax[1,0].set(ylim=(-15,0))
ax[1,0].set(xlim=(0,71))
#ax[1,0].set_title('conc_dd')
ax[1,0].set_ylabel('Depth')
ax[1,0].set_xlabel('distance (m)')
ax[2,0].set(ylim=(-15,0))
ax[2,0].set(xlim=(0,71))
#ax[2,0].set_title('72x1_dd')
ax[2,0].set_ylabel('Depth')
ax[2,0].set_xlabel('distance (m)')
ax[0,1].text(65, -9.5,str(round(mean_18x4_dd_s,3)),size=20,verticalalignment='bottom', horizontalalignment='right',bbox={'facecolor': 'red', 'alpha': 0.2, 'pad': 5})
ax[0,1].set(ylim=(-15,0))
ax[0,1].set(xlim=(0,71))
#ax[0,1].set_title('martorana_18x1')
#ax[0,1].set_ylabel('Depth')
ax[0,1].set_xlabel('distance(m)')
ax[1,1].text(65, -9.5,str(round(mean_18x1_dd_s,3)),size=20,verticalalignment='bottom', horizontalalignment='right',bbox={'facecolor': 'red', 'alpha': 0.2, 'pad': 5})
ax[1,1].set(ylim=(-15,0))
ax[1,1].set(xlim=(0,71))
#ax[1,1].set_title('martorana_72')
#ax[1,1].set_ylabel('Depth')
ax[1,1].set_xlabel('distance (m)')
ax[2,1].text(65, -9.5,str(round(mean_72x1_dd_s,3)),size=20,verticalalignment='bottom', horizontalalignment='right',bbox={'facecolor': 'red', 'alpha': 0.2, 'pad': 5})
ax[2,1].set(ylim=(-15,0))
ax[2,1].set(xlim=(0,71))
#ax[2,1].set_ylabel('Depth')
ax[2,1].set_xlabel('distance (m)')

plt.savefig('inversion_martorana_dd_ref.png')

fig, (ax) = plt.subplots(ncols=2,nrows=3,figsize=(45,20),dpi=300)
#plt.subplots_adjust(hspace=-2)
pg.show(world,fillRegion=False,ax=ax[0,0])
ert_18x4_slm.showResult(cMin=cmin,cMax=cmax,colorBar=False,ax=ax[0,0])
pg.show(world,fillRegion=False,ax=ax[1,0])
ert_18x1_slm.showResult(cMin=cmin,cMax=cmax,colorBar=False,ax=ax[1,0])
pg.show(world,fillRegion=False,ax=ax[2,0])
ert_72x1_slm.showResult(cMin=cmin,cMax=cmax,colorBar=True,ax=ax[2,0])
pg.show(square,fillRegion=False,ax=ax[0,1])
pg.show(grid,martorana_sim_18x4_slm,cMap='bwr',cMin=-1,cMax=1,ax=ax[0,1])
pg.show(square,fillRegion=False,ax=ax[1,1])
pg.show(grid,martorana_sim_18x1_slm,cMap='bwr',cMin=-1,cMax=1,ax=ax[1,1])
pg.show(square,fillRegion=False,ax=ax[2,1])
pg.show(grid,martorana_sim_72x1_slm,cMap='bwr',label='Adj. index',cMin=-1,cMax=1,ax=ax[2,1])
ax[0,0].set(ylim=(-15,0))
ax[0,0].set(xlim=(0,71))
ax[0,0].set_title('18x4_slm')
ax[0,0].set_ylabel('Depth')
ax[0,0].set_xlabel('distance (m)')
ax[1,0].set(ylim=(-15,0))
ax[1,0].set(xlim=(0,71))
ax[1,0].set_title('conc_slm')
ax[1,0].set_ylabel('Depth')
ax[1,0].set_xlabel('distance (m)')
ax[2,0].set(ylim=(-15,0))
ax[2,0].set(xlim=(0,71))
ax[2,0].set_title('72x1_slm')
ax[2,0].set_ylabel('Depth')
ax[2,0].set_xlabel('distance (m)')
ax[0,1].text(65, -9.5,str(round(mean_18x4_sch_s,3)),size=20,verticalalignment='bottom', horizontalalignment='right',bbox={'facecolor': 'red', 'alpha': 0.2, 'pad': 5})
ax[0,1].set(ylim=(-15,0))
ax[0,1].set(xlim=(0,71))
ax[0,1].set_title('martorana_18x1')
ax[0,1].set_ylabel('Depth')
ax[0,1].set_xlabel('martorana_conc')
ax[1,1].text(65, -9.5,str(round(mean_18x1_sch_s,3)),size=20,verticalalignment='bottom', horizontalalignment='right',bbox={'facecolor': 'red', 'alpha': 0.2, 'pad': 5})
ax[1,1].set(ylim=(-15,0))
ax[1,1].set(xlim=(0,71))
ax[1,1].set_title('martorana_72')
ax[1,1].set_ylabel('Depth')
ax[1,1].set_xlabel('distance (m)')
ax[2,1].text(65, -9.5,str(round(mean_72x1_sch_s,3)),size=20,verticalalignment='bottom', horizontalalignment='right',bbox={'facecolor': 'red', 'alpha': 0.2, 'pad': 5})
ax[2,1].set(ylim=(-15,0))
ax[2,1].set(xlim=(0,71))
ax[2,1].set_ylabel('Depth')
ax[2,1].set_xlabel('distance (m)')

plt.savefig('inversion_martorana_sch_ref.png')

#figura amb tots els martoranes

fig, (ax) = plt.subplots(nrows=3,figsize=(30,20),dpi=300)
#plt.subplots_adjust(hspace=1,wspace=1)
pg.show(grid,martorana_sim_18x4_dd,cMap='bwr',cMin=-1,cMax=1,ax=ax[0])
pg.show(grid,martorana_sim_18x1_dd,cMap='bwr',cMin=-1,cMax=1,ax=ax[1])
pg.show(grid,martorana_sim_72x1_dd,cMap='bwr',cMin=-1,label='adj. index',cMax=1,ax=ax[2])
ax[0].set(ylim=(-15,0))
ax[0].set(xlim=(0,71))
ax[0].set_title('18x4_dd')
ax[0].set_ylabel('Depth')
ax[0].set_xlabel('distance (m)')
ax[1].set(xlim=(0,71))
ax[1].set(ylim=(-15,0))
ax[1].set_title('CONC')
ax[1].set_ylabel('Depth')
ax[1].set_xlabel('distance (m)')
ax[2].set(xlim=(0,71))
ax[2].set(ylim=(-15,0))
ax[2].set_title('72x1_dd')
ax[2].set_ylabel('Depth')
ax[2].set_xlabel('distance (m)')

plt.savefig('martorana_dipole_ref.png')

fig, (ax) = plt.subplots(nrows=3,figsize=(30,20),dpi=300)
#plt.subplots_adjust(hspace=1,wspace=1)
pg.show(grid,martorana_sim_18x4_slm,cMap='bwr',cMin=-1,cMax=1,ax=ax[0])
pg.show(grid,martorana_sim_18x1_slm,cMap='bwr',cMin=-1,cMax=1,ax=ax[1])
pg.show(grid,martorana_sim_72x1_slm,cMap='bwr',cMin=-1,cMax=1,label='adj. index',ax=ax[2])
ax[0].set(ylim=(-15,0))
ax[0].set(xlim=(0,71))
ax[0].set_title('18x4_dd')
ax[0].set_ylabel('Depth')
ax[0].set_xlabel('distance (m)')
ax[1].set(xlim=(0,71))
ax[1].set(ylim=(-15,0))
ax[1].set_title('CONC')
ax[1].set_ylabel('Depth')
ax[1].set_xlabel('distance (m)')
ax[2].set(xlim=(0,71))
ax[2].set(ylim=(-15,0))
ax[2].set_title('72x1_dd')
ax[2].set_ylabel('Depth')
ax[2].set_xlabel('distance (m)')

plt.savefig('martorana_sch_ref.png')

fig, (ax) = plt.subplots(ncols=2,nrows=4,figsize=(50,20),dpi=300)
plt.subplots_adjust(hspace=0.5)
pg.show(square,fillRegion=False,ax=ax[0,0])
pg.show(grid,martorana_sim_18x4_dd,cMap='bwr',label='Adjustment index',cMin=-1,cMax=1,ax=ax[0,0])
pg.show(square,fillRegion=False,ax=ax[1,0])
pg.show(grid,martorana_sim_18x4_slm,cMap='bwr',cMin=-1,cMax=1,ax=ax[1,0])
pg.show(square,fillRegion=False,ax=ax[2,0])
pg.show(grid,martorana_sim_18x1_dd,cMap='bwr',cMin=-1,cMax=1,ax=ax[2,0])
pg.show(square,fillRegion=False,ax=ax[3,0])
pg.show(grid,martorana_sim_18x1_slm,cMap='bwr',cMin=-1,cMax=1,ax=ax[3,0])
pg.show(square,fillRegion=False,ax=ax[0,1])
pg.show(grid,martorana_sim_72x1_dd,cMap='bwr',cMin=-1,cMax=1,ax=ax[0,1])
pg.show(square,fillRegion=False,ax=ax[1,1])
pg.show(grid,martorana_sim_72x1_slm,cMap='bwr',cMin=-1,cMax=1,ax=ax[1,1])
pg.show(square,fillRegion=False,ax=ax[2,1])
pg.show(grid,martorana_sim_18x1_dd,cMap='bwr',cMin=-1,cMax=1,ax=ax[2,1])
pg.show(square,fillRegion=False,ax=ax[3,1])
pg.show(grid,martorana_sim_18x1_slm,cMap='bwr',cMin=-1,cMax=1,ax=ax[3,1])
ax[0,0].set_title('18x4_dd')
ax[0,0].set_ylabel('Depth')
ax[0,0].set_xlabel('distance (m)')
ax[0,0].text(30, -4.5,str(round(mean_18x4_dd_s,3)),verticalalignment='bottom', horizontalalignment='right',bbox={'facecolor': 'red', 'alpha': 0.2, 'pad': 5})
ax[1,0].set_title('18x4_slm')
ax[1,0].set_ylabel('Depth')
ax[1,0].set_xlabel('distance (m)')
ax[1,0].text(30, -4.5,str(round(mean_18x4_sch_s,3)),verticalalignment='bottom', horizontalalignment='right',bbox={'facecolor': 'red', 'alpha': 0.2, 'pad': 5})
ax[2,0].set_title('18x1_dd_conc')
ax[2,0].set_ylabel('Depth')
ax[2,0].set_xlabel('distance (m)')
ax[2,0].text(30, -4.5,str(round(mean_18x1_dd_s,3)),verticalalignment='bottom', horizontalalignment='right',bbox={'facecolor': 'red', 'alpha': 0.2, 'pad': 5})
ax[3,0].set_title('18x1_slm_conc')
ax[3,0].set_ylabel('Depth')
ax[3,0].set_xlabel('distance (m)')
ax[3,0].text(30, -4.5,str(round(mean_18x1_sch_s,3)),verticalalignment='bottom', horizontalalignment='right',bbox={'facecolor': 'red', 'alpha': 0.2, 'pad': 5})
ax[0,1].set_title('72x1_dd')
ax[0,1].set_ylabel('Depth')
ax[0,1].set_xlabel('distance (m)')
ax[0,1].text(30, -4.5,str(round(mean_72x1_dd_s,3)),verticalalignment='bottom', horizontalalignment='right',bbox={'facecolor': 'red', 'alpha': 0.2, 'pad': 5})
ax[1,1].set_title('72x1_slm')
ax[1,1].set_ylabel('Depth')
ax[1,1].set_xlabel('distance (m)')
ax[1,1].text(30, -4.5,str(round(mean_72x1_sch_s,3)),verticalalignment='bottom', horizontalalignment='right',bbox={'facecolor': 'red', 'alpha': 0.2, 'pad': 5})
ax[2,1].set_title('18x1_dd_conc')
ax[2,1].set_ylabel('Depth')
ax[2,1].set_xlabel('distance (m)')
ax[2,1].text(30, -4.5,str(round(mean_18x1_dd_s,3)),verticalalignment='bottom', horizontalalignment='right',bbox={'facecolor': 'red', 'alpha': 0.2, 'pad': 5})
ax[3,1].set_title('18x1_slm_conc')
ax[3,1].set_ylabel('Depth')
ax[3,1].set_xlabel('distance (m)')
ax[3,1].text(30, -4.5,str(round(mean_18x1_sch_s,3)),verticalalignment='bottom', horizontalalignment='right',bbox={'facecolor': 'red', 'alpha': 0.2, 'pad': 5})
plt.savefig('martorana_4x4_ref.png')


fig, (ax) = plt.subplots(ncols=2,nrows=4,figsize=(50,20),dpi=300)
plt.subplots_adjust(hspace=0.5)
pg.show(grid2,martorana_sim_18x4_dd_s,cMap='bwr',cMin=-1,cMax=1,ax=ax[0,0])
pg.show(grid2,martorana_sim_18x4_slm_s,cMap='bwr',cMin=-1,cMax=1,ax=ax[1,0])
pg.show(grid2,martorana_sim_18x1_dd_s,cMap='bwr',cMin=-1,cMax=1,ax=ax[2,0])
pg.show(grid2,martorana_sim_18x1_slm_s,cMap='bwr',cMin=-1,cMax=1,ax=ax[3,0])
pg.show(grid2,martorana_sim_72x1_dd_s,cMap='bwr',cMin=-1,cMax=1,ax=ax[0,1])
pg.show(grid2,martorana_sim_72x1_slm_s,cMap='bwr',cMin=-1,cMax=1,ax=ax[1,1])
pg.show(grid2,martorana_sim_18x1_dd_s,cMap='bwr',cMin=-1,cMax=1,ax=ax[2,1])
pg.show(grid2,martorana_sim_18x1_slm_s,cMap='bwr',cMin=-1,cMax=1,ax=ax[3,1])
ax[0,0].set_title('18x4_dd')
ax[0,0].set_ylabel('Depth')
ax[0,0].set_xlabel('distance (m)')
ax[0,0].text(44, -4.5,str(round(mean_18x4_dd_s,3)),verticalalignment='bottom', horizontalalignment='right',bbox={'facecolor': 'red', 'alpha': 0.2, 'pad': 5})
ax[1,0].set_title('18x4_slm')
ax[1,0].set_ylabel('Depth')
ax[1,0].set_xlabel('distance (m)')
ax[1,0].text(44, -4.5,str(round(mean_18x4_sch_s,3)),verticalalignment='bottom', horizontalalignment='right',bbox={'facecolor': 'red', 'alpha': 0.2, 'pad': 5})
ax[2,0].set_title('18x1_dd_conc')
ax[2,0].set_ylabel('Depth')
ax[2,0].set_xlabel('distance (m)')
ax[2,0].text(44, -4.5,str(round(mean_18x1_dd_s,3)),verticalalignment='bottom', horizontalalignment='right',bbox={'facecolor': 'red', 'alpha': 0.2, 'pad': 5})
ax[3,0].set_title('18x1_slm_conc')
ax[3,0].set_ylabel('Depth')
ax[3,0].set_xlabel('distance (m)')
ax[3,0].text(44, -4.5,str(round(mean_18x1_sch_s,3)),verticalalignment='bottom', horizontalalignment='right',bbox={'facecolor': 'red', 'alpha': 0.2, 'pad': 5})
ax[0,1].set_title('72x1_dd')
ax[0,1].set_ylabel('Depth')
ax[0,1].set_xlabel('distance (m)')
ax[0,1].text(44, -4.5,str(round(mean_72x1_dd_s,3)),verticalalignment='bottom', horizontalalignment='right',bbox={'facecolor': 'red', 'alpha': 0.2, 'pad': 5})
ax[1,1].set_title('72x1_slm')
ax[1,1].set_ylabel('Depth')
ax[1,1].set_xlabel('distance (m)')
ax[1,1].text(44, -4.5,str(round(mean_72x1_sch_s,3)),verticalalignment='bottom', horizontalalignment='right',bbox={'facecolor': 'red', 'alpha': 0.2, 'pad': 5})
ax[2,1].set_title('18x1_dd_conc')
ax[2,1].set_ylabel('Depth')
ax[2,1].set_xlabel('distance (m)')
ax[2,1].text(44, -4.5,str(round(mean_18x1_dd_s,3)),verticalalignment='bottom', horizontalalignment='right',bbox={'facecolor': 'red', 'alpha': 0.2, 'pad': 5})
ax[3,1].set_title('18x1_slm_conc')
ax[3,1].set_ylabel('Depth')
ax[3,1].set_xlabel('distance (m)')
ax[3,1].text(44, -4.5,str(round(mean_18x1_sch_s,3)),verticalalignment='bottom', horizontalalignment='right',bbox={'facecolor': 'red', 'alpha': 0.2, 'pad': 5})
plt.savefig('martorana_4x4_s_ref.png')


#figura martorana vs martorana

fig, (ax) = plt.subplots(ncols=2,nrows=2, figsize=(50,20),dpi=300)
pg.show(world,fillRegion=False,ax=ax[0,0])
pg.show(grid,martorana_sim_18x4vs18x1_dd,cMap='bwr',cMin=-1,cMax=1,ax=ax[0,0])
pg.show(world,fillRegion=False,ax=ax[1,0])
pg.show(grid,martorana_sim_18x4vs18x1_slm,cMap='bwr',cMin=-1,cMax=1,ax=ax[1,0])
pg.show(world,fillRegion=False,ax=ax[0,1])
pg.show(grid,martorana_sim_18x1vs72x1_dd,cMap='bwr',cMin=-1,cMax=1,ax=ax[0,1])
pg.show(world,fillRegion=False,ax=ax[1,1])
pg.show(grid,martorana_sim_18x1vs72x1_slm,cMap='bwr',cMin=-1,cMax=1,ax=ax[1,1])
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
plt.savefig('martorana_vs_martorana_ref.png')

fig, (ax) = plt.subplots(figsize=(15,5),dpi=300)
pg.show(mesh,rhomap,showMesh=True,label='Resistivity $(\Omega$m)',cMap='Reds',ax=ax)
#ax.set(xlim=(10,60),ylim=(-14,0))
#ax.set_title('model')
ax.set_ylabel('Depth (m)')
ax.set_xlabel('distance (m)')
plt.savefig('model_red.png')

ert_18x1_dd.saveResult(folder='results18x1dd')
ert_18x1_slm.saveResult(folder='results18x1slm')
ert_18x4_dd.saveResult(folder='results18x4dd')
ert_18x4_slm.saveResult(folder='results18x4slm')
ert_72x1_dd.saveResult(folder='results72x1dd')
ert_72x1_slm.saveResult(folder='results72x1slm')