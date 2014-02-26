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


defPath = "/home/tom/projects/UGP2/datafiles/YahooDirectories"
os.chdir(defPath)
levelEnd = "----level end----"

url_list = list()
cat_dir = list()
url_list.append("http://dir.yahoo.com/arts/")
cat_dir.append("arts")
category = {}
sitesList = {}

while url_list:
	curr_page_url = url_list.pop()
	curr_cat_dir = cat_dir.pop()
	
	print curr_page_url
	if curr_cat_dir == levelEnd:
		os.chdir('..')
	else:
		directory=curr_cat_dir
		if not os.path.exists(directory):
			os.makedirs(directory)
		os.chdir(directory)
		
		#fetching link
		curr_page_url = urlparse.urlsplit(curr_page_url)
		curr_page_url = curr_page_url.geturl()

		req = urllib2.Request(curr_page_url, headers=hdr)
		try:
				f = urllib2.urlopen(req)
		except urllib2.HTTPError, e:
			a=9
		
		except urllib2.URLError, e2:
			a=9

		#storing html page in temporary web1 file
		of = open('web1','w')
		of.write(f.read())
		of.close()

		#opening saved html page for parsing
		fin = open('web1','r')
		soup = BeautifulSoup(fin.read())
		#print(soup.prettify())

		#getting the <Top Categories> tag
		
		topCateg = soup.find("tr",text=re.compile("CATEGORIES"))

		#if categories
		if topCateg is not None:
			reqTable = topCateg.findParents('table')[0]

			#get all the required links from this table
			reqDivs = reqTable.findAll("div",{"class":"cat"})

			for div in reqDivs:
				reqLi = div.findAll("li")
				for rli in reqLi:
					reqA = rli.findAll("a")[0]
					category[reqA.find(text=True)] = urlparse.urljoin(curr_page_url,reqA['href'])+'/'

			if category:
				url_list.append(levelEnd)
				cat_dir.append(levelEnd)
				while category:
					k,v = category.popitem()
					print v,curr_page_url
					url_list.append(v)
					cat_dir.append(k)

			print '*'*80
		else:
			#top sites
			topSites = soup.find("tr",text=re.compile("SITE LISTINGS"))
			reqTable = topSites.findParents('table')[0]

			#get all the required links from this table
			reqDivs = reqTable.findAll("div",{"class":"st"})

			for div in reqDivs:
				reqLi = div.findAll("li")
				for rli in reqLi:
					reqA = rli.findAll("a")[0]
					sitesList[reqA.find(text=True)] = urlparse.urljoin(curr_page_url,reqA['href'])+'/'

			if sitesList:
				of = open('sites','w')
				while sitesList:
					k,v = sitesList.popitem()
					print>>of, k,"\t",v
