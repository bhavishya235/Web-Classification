library("tm")
library("SnowballC")
(dict.sample = Corpus(DirSource("datafiles/dictionaryCorpus", encoding = "UTF-8")))
#sample1 = source("datafiles/sample1")

#inspect(sample1)

#remove extra white space
dict.sample <- tm_map(dict.sample, stripWhitespace)

#remove numbers
dict.sample <- tm_map(dict.sample, removeNumbers)

#remove punctuations
dict.sample <- tm_map(dict.sample, removePunctuation)

#convert to lower case
dict.sample <- tm_map(dict.sample, tolower)

#remove stopwords
dict.sample <- tm_map(dict.sample, removeWords, stopwords("english"))

#stemming
dict.sample <- tm_map(dict.sample,stemDocument)


#inspect(sample1_stem)

#document-term matrix
dict.dtmTf <- DocumentTermMatrix(dict.sample)

#document-term matrix, normalized tf-idf
#dict.dtmTfIdf <- DocumentTermMatrix(dict.sample,control = list(weighting =function(x)weightTfIdf(x, normalize =TRUE)))

#convert to a matrix
dict.matTf = as.matrix(dict.dtmTf)
#dict.matTfIdf = as.matrix(dict.dtmTfIdf)

#get col names i.e. all the words
dict.allWords = colnames(dict.matTf)

#get each word frequency
dict.wordFreq = as.matrix(colSums(dict.matTf))

#save in a file
write(x=t(dict.allWords), file="datafiles/Dictionary.txt", sep="\t")
write.table(x=dict.wordFreq,file="datafiles/DictionaryFrequency.txt",quote=FALSE,sep="\t")


