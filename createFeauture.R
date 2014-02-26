library("tm")
library("SnowballC")
(train.sample = Corpus(DirSource("datafiles/trainingCorpus", encoding = "UTF-8")))

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

#reading dictionary
train.dict = as.matrix(read.table("datafiles/Dictionary.txt"))
train.dict = train.dict[,1]

#document-term matrix
train.dtmTf <- DocumentTermMatrix(train.sample,list(dictionary = train.dict))

#document-term matrix, normalized tf-idf
train.dtmTfIdf <- DocumentTermMatrix(train.sample,control = list(weighting =function(x)weightTfIdf(x, normalize =TRUE),dictionary = train.dict))

#convert to a matrix
train.matTf = as.matrix(train.dtmTf)
train.matTfIdf = as.matrix(train.dtmTfIdf)

#save in a file
write.table(x=train.matTf,file="datafiles/trainFeaturesTf.txt",quote=FALSE,sep="\t")
write.table(x=train.matTfIdf,file="datafiles/trainFeaturesTfIdf.txt",quote=FALSE,sep="\t")