import networkx as nx
import math

'''

# Draw links uavs-uavs from linear solution
-----------------------------------------------------------------------
uavrange=125
# Collect data
theuavs=[]
#with open("finalres.csv") as f:
with open("out/G_0_coords") as f:
	lines=f.readlines()
	for x in lines :
		row=x.rstrip('\n').split(',')
		theuavs.append([float(row[0]),float(row[1])])

thegrounds=[]
with open("coords.csv") as f:
	lines=f.readlines()
	for x in lines :
		row=x.rstrip('\n').split(',')
		thegrounds.append([float(row[0]),float(row[1])])

with open('finallinksgroundsuavs.csv','w') as f:
	for u in theuavs :
		for g in thegrounds :
			euclid_dist = math.sqrt( (u[0]-g[0])**2 + (u[1]-g[1])**2 )
			if(euclid_dist < uavrange) :
				f.write(str(u[0])+","+str(u[1])+"\n")
				f.write(str(g[0])+","+str(g[1])+"\n\n")

--------------------------------------------------------------------------------------

# Draw different path for connectivity
---------------------------------------------------------------
g0=set()
with open("out/G_0_coords") as f:
	lines=f.readlines()
	for x in lines :
		row=x.rstrip('\n').split(',')
		g0.add((float(row[0]),float(row[1])))

g4=set()
with open("out/G_4_coords") as f:
	lines=f.readlines()
	for x in lines :
		row=x.rstrip('\n').split(',')
		g4.add((float(row[0]),float(row[1])))

''' Do same for g3, g5, etc.'''

onlyg3=g3-g0
with open('onlyg3','w') as f1:
	with open('onlyg3_coords','w') as f2:
		for ele in onlyg3 :
			x,y=ele
			f2.write(str(x)+","+str(y)+"\n")
			for u in theuavs :
				euclid_dist = math.sqrt( (x-u[0])**2 + (y-u[1])**2 )
				if(euclid_dist < uavrange) :
					f1.write(str(x)+","+str(y)+"\n")
					f1.write(str(u[0])+","+str(u[1])+"\n\n")

onlyg4=g4-g0
with open('onlyg4','w') as f1:
	with open('onlyg4_coords','w') as f2:
		for ele in onlyg4 :
			x,y=ele
			f2.write(str(x)+","+str(y)+"\n")
			for u in theuavs :
				euclid_dist = math.sqrt( (x-u[0])**2 + (y-u[1])**2 )
				if(euclid_dist < uavrange) :
					f1.write(str(x)+","+str(y)+"\n")
					f1.write(str(u[0])+","+str(u[1])+"\n\n")

# add dashed lines
------------------
uavsasset=set()
#with open("finalres.csv") as f:
with open("out/G_0_coords") as f:
	lines=f.readlines()
	for x in lines :
		row=x.rstrip('\n').split(',')
		uavsasset.add((float(row[0]),float(row[1])))


current=g0-(g0-g3) # then current=g0-onlyg4and so on and so forth
with open('withg3','w') as f:
	for ele in onlyg3 :
		x,y=ele
		for u in current :
			x2,y2=u
			euclid_dist = math.sqrt( (x-x2)**2 + (y-y2)**2 )
			if(euclid_dist <= uavrange) :
				f.write(str(x)+","+str(y)+"\n")
				f.write(str(x2)+","+str(y2)+"\n\n")

notonlyg5=g5-g0
current=g0-(g0-g5)
with open('onlyg5','w') as f:
	for ele in onlyg5 :
		x,y=ele
		for u in current :
			x2,y2=u
			euclid_dist = math.sqrt( (x-x2)**2 + (y-y2)**2 )
			if(euclid_dist <= uavrange) :
				f.write(str(x)+","+str(y)+"\n")
				f.write(str(x2)+","+str(y2)+"\n\n")

'''

uavrange=125
# Collect data
drones_coords=[]
theFiles=["G_"+str(i)+"_coords" for i in range(23)]
for aFile in theFiles :
	with open(aFile) as f:
		lines=f.readlines()
		buff=[]
		for x in lines :
			row=x.rstrip('\n').split(',')
			buff.append([float(row[0]),float(row[1])])
		drones_coords.append(buff)

Gs=[nx.Graph() for i in range(len(drones_coords))]
for l in range(len(drones_coords)) :
	Gs[l].add_nodes_from([0,len(drones_coords[l])-1])
	for i in range(len(drones_coords[l])) :
		for j in range(i+1,len(drones_coords[l])):
			euclid_dist = math.sqrt( (drones_coords[l][i][0]-drones_coords[l][j][0])**2 + (drones_coords[l][i][1]-drones_coords[l][j][1])**2 )
			if(euclid_dist < uavrange) :
				Gs[l].add_edge(i,j)

nodes_cuts=[]
for i in Gs :
	nodes_cuts.append(nx.minimum_node_cut(i))

for j in range(i+1,len(drones_coords[l])):
	for j in range(i+1,len(drones_coords[l])):
	euclid_dist = math.sqrt( (drones_coords[l][i][0]-drones_coords[l][j][0])**2 + (drones_coords[l][i][1]-drones_coords[l][j][1])**2 )
