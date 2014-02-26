import sys
import urllib2,urllib
import urlparse
from urlparse import urljoin
import re
import os
from BeautifulSoup import BeautifulSoup
import Queue

#setting proxy connection
proxy = urllib2.ProxyHandler({'https': 'http://javeshg:Temporary@ironport1.iitk.ac.in:3128',
			      'http': 'http://javeshg:Temporary@ironport1.iitk.ac.in:3128'})
auth = urllib2.HTTPBasicAuthHandler()
opener = urllib2.build_opener(proxy, auth, urllib2.HTTPHandler,urllib2.HTTPSHandler,urllib2.HTTPRedirectHandler)
urllib2.install_opener(opener)
#proxy connection made

hdr = {'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.11 (KHTML, like Gecko) Chrome/23.0.1271.64 Safari/537.11',
       'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
       'Accept-Charset': 'ISO-8859-1,utf-8;q=0.7,*;q=0.3',
       'Accept-Encoding': 'none',
       'Accept-Language': 'en-US,en;q=0.8',
       'Connection': 'keep-alive'}


defPath = "/home/tom/projects/UGP2/datafiles/dmoz/Top1"
os.chdir(defPath)

def getData(root):
	os.chdir(root)
	print root
	fin = open("content",'r').read()
	soup = BeautifulSoup(fin)
	links = soup.findAll("link")
	count=0

	for l in links:
		if(count==20): break
		count=count+1
		url = l["r:resource"]
		print  url
		a=1
		req = urllib2.Request(url, headers=hdr)
		try:
			f = urllib2.urlopen(req)
			of = open("web1",'a')
			tt="="*20+"<" + url + ">"+"="*20+'\n' 
			of.write(tt)
			of.write(f.read())
			tt="+"*20+"<Over>"+"+"*20+'\n'
			of.write(tt)
			of.close()
		except urllib2.HTTPError, e:
			a=9
		except urllib2.URLError, e2:
			a=9

for root,folder,files in os.walk(defPath):
	#print root
	#print folder
	#print files

	if files:
		for f in files:
			if f == "content":
				print "="*20 + "Getting Content" + "="*21
				getData(root);
				print "="*26 + "Done" + "="*26
			else:
				print "*"*20 + "Unknown File >> "+ f + "*"*20
