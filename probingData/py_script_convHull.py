# Collect data
drones_coords=[]
aFile="ref_pbrlm.csv"
with open(aFile) as f:
	lines=f.readlines()
	i=0
	for x in lines :
		if (i<=3) :
			row=x.split('\n')
			i=i+1
			continue
		else:
			row=x.rstrip('\n').split(',')
			drones_coords.append([float(row[0]),float(row[1])])

# compute convex hull and plot
from numpy import array
points=array(drones_coords)
from scipy.spatial import ConvexHull, convex_hull_plot_2d
hull = ConvexHull(points)

import csv
with open('convHullcoords.csv', "w") as f:
	for ele in hull.vertices :
		f.write(str(drones_coords[ele][0])+","+str(drones_coords[ele][1])+"\n")

'''
# finding the frame for the grid
drones_coords=[]
with open("convHullcoords.csv") as f:
	lines=f.readlines()
	for x in lines :
		row=x.rstrip('\n').split(',')
		drones_coords.append([float(row[0]),float(row[1])])

minx=10000;miny=10000
maxx=0;maxy=0
for ele in drones_coords :
	if (ele[0] < minx) :
		minx=ele[0]
	if (ele[0] > maxx) :
		maxx=ele[0]
	if (ele[1] < miny) :
		miny=ele[1]
	if (ele[1] > maxy) :
		maxy=ele[1]

# found 20.6005730037, 23.4401850852 to 974.566876883, 967.101020767
# 974-20 = 954; 967-23= 944
# At the end chose to split with respect to whole map : 125*8 = 1000

# finding the coordinates for the grid without point (875,125)
with open('truegriddrawing.csv','w') as f:
	for move in [0,1] :
		for val in a :
			if (move==0) :
				if (val==875) :
					continue
				else :
					f.write(str(val)+",125\n"+str(val)+",875\n\n")
			if (move==1) :
				if (val==125) :
					continue
				else :
					f.write("125,"+str(val)+"\n"+"875,"+str(val)+"\n\n")
	f.write("125,125\n750,125\n\n875,250\n875,875")

'''