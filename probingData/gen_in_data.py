alldata=[]
import random
from random import randint

nuavs=500
uavsrange=125

alldata=[]
x0,xinf=0,10000
y0,yinf=0,10000

for i in range(nuavs) :
	alldata.append([round(random.uniform(x0,xinf),2),round(random.uniform(y0,yinf),2)])# round to 2nd member decimal

f=open("dummy"+str(nuavs),"w")
f.write(str(18)+"\n")# to correct, at stage useless, depicts number of available uavs
f.write(str(uavsrange)+"\n")
f.write(str(nuavs)+"\n")
f.write(str(xinf)+","+str(yinf)+"\n")

for point in alldata :
	f.write(str(point[0])+","+str(point[1])+"\n")
f.close()
