import networkx as nx
import math

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
