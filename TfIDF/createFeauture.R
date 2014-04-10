#Get the files in directory AllVector1/$i where $i varies from 1 to 72 and for the documents(around 100) present in this directory print their tf in file FeaturesTf$i
library("tm")
library("SnowballC")
#reading dictionary
train.dict = as.matrix(read.table("Dictionary.txt"))
train.dict = train.dict[,1]
for(i in 71:72){
	print(i)
	t <- paste("AllVector1/",i,sep='')
	(train.sample = Corpus(DirSource(t, encoding = "UTF-8")))

	#inspect(sample1)

	#remove extra white space
	train.sample <- tm_map(train.sample, stripWhitespace)

	#remove numbers
	train.sample <- tm_map(train.sample, removeNumbers)

	#remove punctuations
	train.sample <- tm_map(train.sample, removePunctuation)

	#convert to lower case
	train.sample <- tm_map(train.sample, tolower)

	#remove stopwords
	train.sample <- tm_map(train.sample, removeWords, stopwords("english"))

	#stemming
	train.sample <- tm_map(train.sample,stemDocument)

	#document-term matrix
	train.dtmTf <- DocumentTermMatrix(train.sample,list(dictionary = train.dict))

	#document-term matrix, normalized tf-idf
	train.dtmTfIdf <- DocumentTermMatrix(train.sample,control = list(weighting =function(x)weightTfIdf(x, normalize =TRUE),dictionary = train.dict))

	#convert to a matrix
	train.matTf = as.matrix(train.dtmTf)
	train.matTfIdf = as.matrix(train.dtmTfIdf)

	#save in a file
	write.table(x=train.matTf,file=paste("FeaturesTf",i,sep=''),quote=FALSE,sep="\t")
	write.table(x=train.matTfIdf,file=paste("FeaturesTfIdf",i,sep=''),quote=FALSE,sep="\t")
}

