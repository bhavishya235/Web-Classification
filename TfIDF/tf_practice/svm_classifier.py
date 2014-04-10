from sklearn import svm
import os
path="tf/"
training_tuple=[]

for root,folder,files in os.walk(path+"training/"):
	if files:
		for f in files:
			try:
				print f
				trusted_training_set=open(path+"training/"+f,'r').readlines()
				for line in trusted_training_set:
					category,rem=line.split(',',1)
					if "Computers" in category:
						line="[1,"+rem
					else:
						line="[0,"+rem 
					training_tuple.append(eval(line))
			except:
				print "please correct later..."



X=[]
y=[]

for tup in training_tuple:
	X.append(tup[1:])
	y.append(tup[0])


clf = svm.SVC()
#print X
#print y

clf.fit(X, y)
print "="*20, "Trained"
test_tuple=[]


for root,folder,files in os.walk(path+"test/"):
	if files:
		for f in files:
			try:
				print f
				trusted_test_set=open(path+"test/"+f,'r').readlines()
				for line in trusted_test_set:
					category,rem=line.split(',',1)
					if "Computers" in category:
						line="[1,"+rem
					else:
						line="[0,"+rem
					test_tuple.append(eval(line))
			except: 
				print "correct later..."


tX=[]
ty=[]
sentence_index = [];
for tup in test_tuple:
	tX.append(tup[1:])
	ty.append(tup[0])
	sentence_index.append(tup[0])
	
total=correct=0.0
for tx,ty2,sen in zip(tX,ty,sentence_index):
	total+=1
	ans = clf.predict([tx])
	ty2=ty2+0.0
	print ans,ty2,sen
	if(ans==ty2): correct+=1

print "% correct = ",correct/total
