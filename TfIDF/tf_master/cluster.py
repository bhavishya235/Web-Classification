"""
It tries the clustering with k=4 clusters and goes on till 94.
The categories maintains the original known categories.

For each cluster:
	The output is "clustur" "highest category occuring in this clustur" "no of its occurences" "total clustur strength"


"""

import os
import mlpy
path="tf/files-for-clustur/"

cltuple=[]
catname=[]
catnum=[]
categories=[]
for root,folder,files in os.walk(path):
	if files:
		for f in files:
			try:
				print f
				trusted_training_set=open(path+f,'r').readlines()
				for line in trusted_training_set:
					category,rem=line.split(',',1)
					cltuple.append(eval("["+rem))
					newCat=category.split('-')[1]
					#Change above line to implement heirarchical clustering
					if newCat not in categories:
						categories.append(newCat)
					catname.append(newCat)
					catnum.append(categories.index(newCat))					 
			except:
				print "please correct later..."

def transpose(grid):
    return zip(*grid)

def removeBlankRows(grid):
    return [list(row) for row in grid if sum(row)!=0]


cltuple= removeBlankRows(transpose(removeBlankRows(transpose(cltuple))))
for it in range(0,len(categories)):
	print "no of doc in cat",it,"=",catnum.count(it)

print "vector len=",len(cltuple[0])
print "no of doc=",len(cltuple)
	
#computing K-Means with K = 2 (2 clusters)
for i in range(4,100,10):
	clusts={}
	for j in range(0,i):
		clusts[j]=[]
	kmeans = mlpy.Kmeans(k=i, init="plus", seed=0)
	cls=kmeans.compute(cltuple)
	print "steps=",kmeans.steps
	k=0
	for item in cls:
		clusts[item].append(catnum[k])
		#print item,catnum[k]
		k=k+1
	for j in range(0,i):
		maincat=max(set(clusts[j]), key=clusts[j].count)
		print j,maincat,clusts[j].count(maincat),len(clusts[j])
	#print clusts
	r=raw_input("paused...")
