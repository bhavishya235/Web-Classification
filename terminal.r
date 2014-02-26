library("tm")
library("SnowballC")
(sample1 = Corpus(DirSource("datafiles", encoding = "UTF-8")))
inspect(sample1)
inspect(sample1)
sample1_spaces <- tm_map(sample1, stripWhitespace)
sample1_stop <- tm_map(sample1_spaces, removeWords, stopwords("english"))
sample1_stem <- tm_map(sample1_stop,stemDocument)
inspect(sample1_stem)
quit
exit
bye
quit("yes")
inspect()
x = 9
inspect(x)
q
q(save="yes")
source("createDictionary.R")
dict.mat
colSums(x=dict.mat)
dict.allWords
q(save="yes")
source("createFeauture.R")
train.matTf
d = Dictionary(c('also','euclidean'))
d = Dictionary(c('also','euclidean'))
inspect(DocumentTermMatrix(train.sample, list(dictionary = d)))
inspect(DocumentTermMatrix(train.sample, list(dictionary = c('d'))))
inspect(DocumentTermMatrix(train.sample, list(dictionary = c('also'))))
inspect(DocumentTermMatrix(train.sample, list(dictionary = c('also','euclidian'))))
dic = read("datafiles/Dictionary.txt")
dic = read.table("datafiles/Dictionary.txt")
dic
str(dic)
m= as.matrix(dict)
m= as.matrix(dic)
m
l = m[,1]
l
inspect(DocumentTermMatrix(train.sample, list(dictionary = l)))
quit(save = "feat.R")
quit(save = "terminal.r")
quit(save="terminal.r")
savehistory(file="terminal.r")
