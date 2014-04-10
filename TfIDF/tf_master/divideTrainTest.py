"""
It takes a set of feature vector files and partition them according to their category

"""

import os
import mlpy
path="tf2/"

categories=[]
numFiles={}
for root,folder,files in os.walk(path):
	if files:
		for f in files:
			print f
			LINES=open(path+f,'r').readlines()[1:]
			for line in LINES:
				line=line.replace('	',',')
				category,rem=line.split(',',1)
				newCat=category.split('-')[1:-1]
				newpath='/'.join(newCat)
				if newpath not in numFiles:
					numFiles[newpath]=1
					print "creating >> ",newpath
					if not os.path.exists("top/"+newpath):
						os.makedirs('top/'+newpath)
					print "created >> ",newpath
				else:
					numFiles[newpath]+=1
				num = numFiles[newpath]
				fin=open('top/'+newpath+'/'+str(num),'w')
				print>>fin,line
