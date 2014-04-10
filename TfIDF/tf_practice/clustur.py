import os
import mlpy
path="tf/files-for-clustur/"

cltuple=[]
for root,folder,files in os.walk(path):
	if files:
		for f in files:
			try:
				print f
				trusted_training_set=open(path+f,'r').readlines()
				for line in trusted_training_set:
					category,rem=line.split(',',1)
					cltuple.append(eval("["+rem))
			except:
				print "please correct later..."

def transpose(grid):
    return zip(*grid)

def removeBlankRows(grid):
    return [list(row) for row in grid if sum(row)!=0]

cltuple= removeBlankRows(transpose(removeBlankRows(transpose(cltuple))))
print len(cltuple[0])
#computing K-Means with K = 2 (2 clusters)
kmeans = mlpy.Kmeans(k=4, init="plus", seed=0)
kmeans.compute(cltuple)
print kmeans.means
print kmeans.steps
