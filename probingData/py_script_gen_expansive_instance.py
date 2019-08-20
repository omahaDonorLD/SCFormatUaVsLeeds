import math
import random
import sympy as sym

def eucl_dist(a,b):
	return math.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2)

T=30
delta=10
uavsrange=125
sln=[]
## init
true_sln=[]
stack_sln=[]
true_sln.append([7.43344,12.22022])
stack_sln.append([7.43344,12.22022])
x=sym.symbols('x')
y=sym.symbols('y')
quadrants=[ [x-uavsrange-delta,x-uavsrange,y-uavsrange-delta,y+uavsrange] , [x-uavsrange,x+uavsrange+delta,y-uavsrange-delta,y-uavsrange] , [x+uavsrange,x+uavsrange+delta,y-uavsrange,y+uavsrange+delta] , [x-uavsrange-delta,x+uavsrange,y+uavsrange,y+uavsrange+delta] ]

k=0

while len(stack_sln) > 0 :
	# test all quadrants
	xk,yk=stack_sln.pop()
	for i in range(4):
		# for each quadrant, do 3 trials
		for j in range(3):# for each quadrant, do 3 trials
			gotit=False
			xstart=quadrants[i][0].subs( dict(x=xk) )
			xend=quadrants[i][1].subs( dict(x=xk) )
			ystart=quadrants[i][2].subs( dict(y=yk) )
			yend=quadrants[i][3].subs( dict(y=yk) )
			xk_1=random.uniform( xstart, xend)
			yk_1=random.uniform( ystart, yend )
			k=k+1
			#print(xstart,",",xend,",",ystart,",",yend,",",xk_1,",",yk_1," ite ", k)
			if xk_1 <= 0 or xk_1 >= 1000 or yk_1 <= 0 or yk_1 >= 1000 :# skip if out of limits
				continue
			
			gotit=True
			for uav in true_sln :
				if eucl_dist([xk_1,yk_1],uav) <= uavsrange :
					gotit=False
					break

			if gotit	 : # new candidate is good
				print("in if",xk_1,",",yk_1,"-",xk,",",yk,"===>",len(stack_sln))
				stack_sln.append([xk,yk])# put back current solution in solutions to visit as might be needed later
				stack_sln.append([xk_1,yk_1])# add new admissible solution
				true_sln.append([xk_1,yk_1])


print (true_sln)
print ("stack size",len(stack_sln))
