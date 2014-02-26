import BeautifulSoup
import urllib
import os
import re
import Queue

defPath = "/home/tom/projects/UGP2/datafiles/dmoz/Top2"
os.chdir(defPath)

def getTextData(html):
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
		fout=open(root+"/web1",'a')
		print>>fout,st
	except:
		print "*"*20+"error"+"*"*20

for root,folder,files in os.walk(defPath):
	if files:
		for f in files:
			if f == "web1":
				ender="+"*20
				print root+'/'+f
				html=""
				filelines=open(root+"/web1",'r').readlines()
				fout = open(root+"/web1",'w')
				fout.close()
				for line in filelines:
					if ender in line:
						getTextData(html);
						html=""
					else:
						html+=line
			else:
				os.system("rm "+root+'/'+f)
