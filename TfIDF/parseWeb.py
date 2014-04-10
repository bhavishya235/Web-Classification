#The files in def path had =x20 and *x20 for seperating two webpages.
#This extracts the textual prime data from them and stores in writefolder.
import BeautifulSoup
import urllib
import os
import re
import Queue

defPath = "/home/tom/projects/UGP2/datafiles/dmoz/Top3"
WriteFolder="/home/tom/projects/TfIDF/AllVector"
os.chdir(defPath)
Dict=[]

def getTextData(html,path,i):
	try:
		soup = BeautifulSoup.BeautifulSoup(html)
		texts = soup.findAll(text=True)
		remWords=[['&[^ ]*;',''],['[ ]* ',' '],['\t',' ']]
		stopList=['===================','<!.*>','<?.*?>','\+{20}']
		def visible(element):
		    if element.parent.name in ['style', 'script', '[document]', 'head', 'title']:
		        return False
		        return False
		    for stop in stopList:
			if re.match(stop,str(element)):
				return False
		    return True

		visible_texts = filter(visible, texts)
		s= visible_texts
		st=""
		for p in s:
			if '\n' not in p:	#print st
			        st=st+str(p)
			st=st+" "
		for rem in remWords:
			st=re.sub(rem[0],rem[1],str(st))
		#print "here"
		#print str(path),i
		path=re.sub('/','-',path)
		path=path.split('-home-tom-projects-UGP2-datafiles-dmoz-')[1]
		print str(path)+str(i)
		fout=open(WriteFolder+"/"+path+'-'+str(i),'w')
		print>>fout,st

	except:
		print "*"*20+"error"+"*"*20

for root,folder,files in os.walk(defPath):
	if files:
		i=1
		for f in files:
			if f == "web1":
				ender="+"*20
				print root+'/'+f
				html=""
				filelines=open(root+"/web1",'r').readlines()
				#fout = open(root+"/web1",'w')
				#fout.close()
				for line in filelines:
					if ender in line:
						getTextData(html,root,i);
						i=i+1
						html=""
					else:
						html+=line
			else:
				os.system("rm "+root+'/'+f)
