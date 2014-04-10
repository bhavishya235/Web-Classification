# Iteratively run for each class
# This file is the master file which creates the training and test set based on the h-class defined.
# Create classifier for each level
# Calculate accuracy both level wise and averall

import os
from sklearn.ensemble import RandomForestClassifier 
path="tf/"

category = ["Arts"]
categoryDepth = len(category)

training_tuple = [[],[],[]]	# a multi-dimension array that stores the tuples for training for all heirarchies
totalVectors = 0
favVectors = 0
otherVectors = 0

#creating training feature vectors for all the heirarchical classes
for root,folder,files in os.walk(path+"training2/"):
	if files:
		for f in files:
			print f
			try:
				lines=open(path+"training/"+f,'r').readlines()
				for line in lines:
					totalVectors+=1
					lineFeature = []
					cat,rem=line.split(',',1)
					cat= cat.split('-')[1:categoryDepth+1]
					print cat
					catDepth=len(cat)
					res=''
					for i in range(0,min(categoryDepth,catDepth)):
						if cat[i]==category[i]:
							lineFeature="[1,"+rem
							training_tuple[i].append(eval(lineFeature))
							#training_tuple[i].append([totalVectors,1])
							res=res+',1'
							favVectors+=1
						else:
							lineFeature="[0,"+rem
							training_tuple[i].append(eval(lineFeature))
							#training_tuple[i].append([totalVectors,0])
							res=res+',0'
							otherVectors+=1
							break		# as deeper level will only be checked when above level satisfies
					print res
			except:
				print "please correct later..."

'''
for t in training_tuple[0]:
	print t

print "*"*50

for t in training_tuple[1]:
	print t

print "*"*50

for t in training_tuple[2]:
	print t
'''

print "favVectors=",favVectors
print "otherVectors=",otherVectors

X=[]
Y=[]
#clf = [svm.SVC(),svm.SVC(),svm.SVC()]
Forest = RandomForestClassifier(n_estimators = 50)

for i in range(categoryDepth):
	X.append([])
	Y.append([])
	for tup in training_tuple[i]:
		X[i].append(tup[1:])
		Y[i].append(tup[0])
	#clf[i].fit(X[i], Y[i])
	print "Training forest"
	Forest = Forest.fit(X[i],Y[i])

print "="*20, "Trained clf0"
print "="*20, "Prepairing Test Set"

#creating test feature vectors for all the heirarchical classes
totalTestVectors=0
test_tuple = []
testVectors = 0
otherTestVectors = 0

for root,folder,files in os.walk(path+"test/"):
	if files:
		for f in files:
			print f
			#try:
			lines=open(path+"test/"+f,'r').readlines()
			for line in lines:
				totalTestVectors+=1
				lineFeature = []
				cat,rem=line.split(',',1)
				cat= cat.split('-')[1:categoryDepth+1]
				print cat
				res=''
				if cat == category:
					lineFeature="[1,"+rem
					test_tuple.append(eval(lineFeature))
					#test_tuple.append([totalTestVectors,1])
					res=res+',1'
					testVectors+=1
					print res
				else:
					lineFeature="[0,"+rem
					test_tuple.append(eval(lineFeature))
					#test_tuple.append([totalVectors,0])
					res=res+',0'
					otherTestVectors+=1
					print res
			#except:
			#	print "please correct later..."

print "="*20, "Prepaired Test Set"

#classifying each vector heirarchicaly
tX=[]
tY=[]
sentence_index = [];
for tup in test_tuple:
	tX.append(tup[1:])
	tY.append(tup[0])
	sentence_index.append(tup[0])
	
total=correct=0.0
for tx,ty,sen in zip(tX,tY,sentence_index):
	for i in range(0,categoryDepth):
		total+=1
		ans = Forest.predict([tx])
		print i,ans,ty
		print type(ans)
		if ty==0:
			if ans==ty:
				correct+=1
			break
		else:
			if ans==ty:
				correct+=1
			else:
				break

print "% correct = ",correct*100.0/total
print "favVectors=",favVectors
print "otherVectors=",otherVectors
print "testVectors=",testVectors
print "otherTestVectors=",otherTestVectors