

## load matlab library for the 'reshape' function code, why not?
#
library(matlab)
library(zoo)
library(stringr)


`%+%` <- function(a, b) paste(a, b, sep="")

## Compute Z normalization for the timeseries.
##  parameters:
##  ts - the timeseries to normalize
##
##  ts should be a matrix, we do care only about first row
#

findMtf = function(mtf){
  res = matrix(ncol = nchar(x)[1], nrow = length(x))
  res[] = 0
  tmp = str_locate_all(x,mtf)
  
  for(i in 1:length(tmp)){
    if(length(tmp[[i]]) == 0) next
    for(j in 1:nrow(tmp[[i]])){
      res[i,tmp[[i]][j,1]:tmp[[i]][j,2]] = 99
    }
  }
  return(res)
}

absExtrema = function(ww){
  ww[] = 0
  
  rangeStr = function(str){
    me = max(as.numeric(unlist(strsplit(str,""))),na.rm = T)
    ne = min(as.numeric(unlist(strsplit(str,""))),na.rm = T)
    return(c(as.character(me),as.character(ne)))
  }
  
  for(i in 1:nrow(ww)){
    me = rangeStr(x[i])[1]
    tmp = str_locate_all(x[i],perl(me %+% "+"))
    tmp = tmp[[1]]
    idx = which(tmp[,1] == 1:len.motif)
    if(length(idx) != 0){
      print("Max at start, Line: " %+% i)
      tmp = tmp[-idx ,]
      if(is.matrix(tmp) == F) tmp = t(as.matrix(tmp))
    }
    t.diff = tmp[,2] - tmp[,1]
    if(length(t.diff) == 0){
      print("No max, Line: " %+% i)
    } else if(length(t.diff) == 1){
      t.idx = 1
      ww[i,tmp[t.idx,1]:tmp[t.idx,2]] = 99
    } else{
      t.max = max(t.diff)
      t.idx = which(t.diff == t.max)
      for (j in 1:length(t.idx)){
        ww[i,tmp[t.idx[j],1]:tmp[t.idx[j],2]] = 99
      }
    }
    
    me = rangeStr(x[i])[2]
    tmp = str_locate_all(x[i],perl(me %+% "+"))
    tmp = tmp[[1]]
    idx = which(tmp[,1] == 1:len.motif)
    if(length(idx) != 0){
      print("Min at start, Line: " %+% i)
      tmp = tmp[-idx ,]
      if(is.matrix(tmp) == F) tmp = t(as.matrix(tmp))
    }
    t.diff = tmp[,2] - tmp[,1]
    if(length(t.diff) == 0){
      print("No min, Line: " %+% i)
    } else if(length(t.diff) == 1){
      t.idx = 1
      ww[i,tmp[t.idx,1]:tmp[t.idx,2]] = -99
    } else{
      t.max = max(t.diff)
      t.idx = which(t.diff == t.max)
      for (j in 1:length(t.idx)){
        ww[i,tmp[t.idx[j],1]:tmp[t.idx[j],2]] = -99
      }
    }
  }
  return(ww)
}

last2last = function(lag, symb, exchange){
  
  lag = abs(lag) + 1
  path = "https://www.google.com/finance/getprices?q=" %+% symb %+%
  "&x=" %+% exchange %+% "&" %+% "i=60&" %+% "p=" %+% lag %+% 
  "d&" %+% "f=d,o"
  
  xx = read.csv(path,skip=6,row.names=NULL)
  ## o.s = as.numeric(unlist(strsplit(readLines(path,n=7)[7],"="))[2])
  dd = as.numeric(unlist(strsplit(xx[1,1],"a"))[2])
  ## dd = dd + 60 * o.s
  
  xx.end = grep("a",xx[,1])[2]
  
  if(!is.na(xx.end)) xx = xx[1:(xx.end - 1),]
  
  xx[1,1] = dd
  xx[-1,1] = (as.numeric(xx[-1,1]))*60 + dd 
  
  xx[,1] = as.POSIXct(as.numeric(xx[,1]), origin = '1970-01-01')
  xx[,1] = as.POSIXct(format(xx[,1],tz="Asia/Jerusalem"))
  ts = zoo(x = xx[,2], order.by = as.POSIXct(xx[,1]))
  return(ts)
}

last4days = function(lag, mints, symb, exchange){
  
  lag = abs(lag) + 1
  mints = as.character(60 * mints)
  path = "https://www.google.com/finance/getprices?q=" %+% symb %+% "&x=" %+% exchange %+% "&i=" %+% mints %+% "&p=" %+% lag %+% "d&" %+% "f=d,o"
  
  xx = read.csv(path,skip=6,row.names=NULL)
  
  idx = which(diff(as.numeric(xx[,1])) > 1)+1
  idx = c(idx,which(is.na(diff(as.numeric(xx[,1]))))+1)
  idx = na.omit(idx)
  idx = unique(idx)
  
  dd = as.numeric(unlist(strsplit(xx[1,1],"a"))[2])
  dd = as.POSIXct(dd, origin = '1970-01-01')
  xx[,1] = seq(dd,Sys.time(),length.out = nrow(xx))
  
  ts = zoo(x = xx[,2], order.by = as.POSIXct(xx[,1]))
  ts = ts[complete.cases(ts)]
  return(ts)
}


last2sum = function(n, symb, exchange){
  if(abs(n) > 10) stop("N is Larger then abs(10) !!!")
  ts = last2last(n, symb, exchange)
  
  for(i in (n+1):0){
    ts = rbind(ts,last2last(i, symb, exchange))
  }
  return(ts)
}

cumsum2G25 = function(n){
  if(abs(n) > 10) stop("N is Larger then abs(10) !!!")
  ts = last2G25(n)
  
  for(i in (n+1):0){
    ts = rbind(ts,last2G25(i))
  }
  return(ts)
}

trimExtrema = function(nip){
  len = vector()
  for(i in 1:nrow(nip)){
    len[i] = max(rle(nip[i,])$
                   lengths[rle(nip[i,])$values == 99])/s.char[i]
  }
  idx = which(len >= 0.7)
  for(i in idx){
    ddx = which(nip[i,] == 99)
    nip[idx,ddx] = 0
  }
  
  for(i in 1:nrow(nip)){
    len[i] = max(rle(nip[i,])$
                   lengths[rle(nip[i,])$values == -99])/s.char[i]
  }
  idx = which(len >= 0.7)
  for(i in idx){
    ddx = which(nip[i,] == -99)
    nip[idx,ddx] = 0
  }
  return(nip)
}

padZeros = function(tmp){
  for(i in 1:length(tmp)){
    zeros = paste(rep("-",max(nchar(tmp))-nchar(tmp[i])),collapse="")
    tmp[i] = tmp[i] %+% zeros
  }
  return(tmp)
}

sax2Matrix = function(tmp){
  uu = as.numeric(unlist(strsplit(tmp[1],"")))
  for(i in 2:length(tmp)) uu = rbind(uu,as.numeric(unlist(strsplit(tmp[i],""))))
  uu = as.matrix(uu)
  return(uu)
}


znorm <- function(ts){
 ts.mean <- mean(ts[1,])
 ts.dev <- sd(ts[1,])
 (ts - ts.mean)/ts.dev
}

fact = function(x){
  x = as.integer(x)
  div = seq_len(abs(x))
  factors = div[x %% div == 0L]
  return(factors)
}

last2G25 = function(lag){
  
  lag = abs(lag) + 1
  
  path = "https://www.google.com/finance/getprices?q=T25&x=TLV&" %+%
    "i=60&" %+%
    "p=" %+% lag %+% "d&" %+%
    "f=d,o"
  
  xx = read.csv(path,skip=6,row.names=NULL)
  
  len = min(nrow(xx),451)
    
  dd = xx[1,1]
  dd = as.numeric(sub("a","",dd))
  dd = as.POSIXct(dd + 5*60*60, origin = '1970-01-01')
  
  dayIDX = seq.POSIXt(dd,by="60 sec",length.out=len)
  
  ts = as.zoo(xx[1:len,2])
    
  index(ts) = dayIDX
  return(ts)
}

lastL25 = function(pp)
{
  xx = read.table(pp,skip=1)
  xx = zoo(cbind(xx[,2]),order.by=xx[,1])
  base = as.numeric(readLines(pp,n=1))
  coredata(xx) = as.numeric(coredata(xx))
  index(xx) = as.POSIXct(as.character(index(xx)),format="%T",tz="Asia/Jerusalem")
  return(xx)
 }
 

last2Nasdaq = function(lag){
  
  lag = abs(lag) + 1
  
  path = "https://www.google.com/finance/getprices?q=IXIC&" %+%
    "i=60&" %+%
    "p=" %+% lag %+% "d&" %+%
    "f=d,o"
  
  
  xx = read.csv(path,skip=6,row.names=NULL)
  
  len = min(nrow(xx),451)
  
  dd = xx[1,1]
  dd = as.numeric(sub("a","",dd))
  dd = as.POSIXct(dd, origin = '1970-01-01')
  
  dayIDX = seq.POSIXt(dd,by="60 sec",length.out=len)
  
  ts = as.zoo(xx[1:len,2])
  
  index(ts) = dayIDX
  return(ts)
}

ts2word = function(pp, shrink = 4){
  xx = read.table(pp,skip=1)
  xx = zoo(cbind(xx[,2]),order.by=xx[,1])
  coredata(xx) = as.numeric(coredata(xx))
  index(xx) = as.POSIXct(as.character(index(xx)),format="%T",tz="Asia/Jerusalem")
  
  xx.data = as.matrix(t(coredata(xx)))
  xx.norm = znorm(xx.data)
  
  aSize = 9
  
  drop = ncol(xx.norm) %% shrink
  drop = length(xx.norm)-drop
  xx.norm = as.matrix(t(xx.norm[1:drop]))
  
  ## adist("aaaab","aaabaabbbb")
  ## adist("aaaab","aaaba")
  
  paaSize = ncol(xx.norm)/shrink
  xx.paa = paa(xx.norm, paaSize)
  ## xx.str = ts2string(xx.paa, aSize)
  xx.str = as.character(ts2num(xx.paa, aSize))
  ## print(xx.str)
  xx.word = paste(xx.str, sep="", collapse="")
  return(xx.word)
}

tsG2word = function(ts, shrink = 5){
  xx = ts
  xx.data = as.matrix(t(coredata(xx)))
  xx.norm = znorm(xx.data)
  
  aSize = 9
  
  drop = ncol(xx.norm) %% shrink
  drop = length(xx.norm)-drop
  xx.norm = as.matrix(t(xx.norm[1:drop]))
  
  ## adist("aaaab","aaabaabbbb")
  ## adist("aaaab","aaaba")
  
  paaSize = ncol(xx.norm)/shrink
  xx.paa = paa(xx.norm, paaSize)
  ## xx.str = ts2string(xx.paa, aSize)
  xx.str = as.character(ts2num(xx.paa, aSize))
  ## print(xx.str)
  xx.word = paste(xx.str, sep="", collapse="")
  return(xx.word)
}

tsDB2word = function(ts, shrink = 10){
  xx = ts
  xx.data = as.matrix(t(coredata(xx)))
  xx.norm = znorm(xx.data)
  
  aSize = 9
  
  drop = ncol(xx.norm) %% shrink
  drop = length(xx.norm)-drop
  xx.norm = as.matrix(t(xx.norm[1:drop]))
  
  ## adist("aaaab","aaabaabbbb")
  ## adist("aaaab","aaaba")
  
  paaSize = ncol(xx.norm)/shrink
  xx.paa = paa(xx.norm, paaSize)
  ## xx.str = ts2string(xx.paa, aSize)
  xx.str = as.character(ts2num(xx.paa, aSize))
  ## print(xx.str)
  xx.word = paste(xx.str, sep="", collapse="")
  return(xx.word)
}

##
## Compute PAA approximation for the timeseries with reduction
##  parameters:
##  ts - timeseries
##  ap - number of points in approximated timeseries
#
paa <- function(ts, ap){
 len <- ncol(ts)
 res <- ts
 if(len != ap){
  if( (len %% ap) == 0 ){
   res <- reshape(ts, len %/% ap, ap)
  }else{
   tmp <- matrix( rep(0, ap*len), ap, len)
   for(i in 1:ap){
    tmp[i, ] <- ts[1, ]
   }
   extended <- reshape(tmp, 1, ap*len)
   res <- reshape(extended, len, ap)
  }
 }
 matrix(colMeans(res), nrow=1, ncol=ap)
}


startG25idx = function(){
  path = "https://www.google.com/finance/getprices?q=T25&x=TLV&" %+%
  "i=60&" %+%
  "p=1d&" %+% 
  "f=d,o"  
  tmp = cbind(unlist(strsplit(readLines(path,n=8)[8],",")))
  tmp = tmp[1]
  tmp = as.numeric(sub("a","",tmp))
  tmp = as.POSIXct(tmp, origin = '1970-01-01')
  return(tmp)
}

lastG25 = function(){
  path = "https://www.google.com/finance/getprices?q=T25&x=TLV&" %+%
    "i=60&" %+%
    "p=1d&" %+% 
    "f=d,o"
  xx = read.csv(path,skip=7)
  xx = xx[,2]
  
  dayStart = startG25idx()
  dayIDX = seq.POSIXt(dayStart,by="60 sec",length.out=length(xx))
  
  ts = as.zoo(xx)
  index(ts) = dayIDX
  return(ts)
}
##
## Converts the specified resolution into the cut points
##
alphabet2cut <- function(alphabet_size){
 switch(alphabet_size,
  0.00,
  c(-Inf,  0.00),
  c(-Inf, -0.43,  0.43),
  c(-Inf, -0.67,  0.00,  0.67),
  c(-Inf, -0.84, -0.25,  0.25,  0.84),
  c(-Inf, -0.97, -0.43,  0.00,  0.43,  0.97),
  c(-Inf, -1.07, -0.57, -0.18,  0.18,  0.57,  1.07),
  c(-Inf, -1.15, -0.67, -0.32,  0.00,  0.32,  0.67,  1.15),
  c(-Inf, -1.22, -0.76, -0.43, -0.14,  0.14,  0.43,  0.76,  1.22),
  c(-Inf, -1.28, -0.84, -0.52, -0.25,  0.00,  0.25,  0.52,  0.84,  1.28),
  c(-Inf, -1.34, -0.91, -0.60, -0.35, -0.11,  0.11,  0.35,  0.60,  0.91, 1.34),
  c(-Inf, -1.38, -0.97, -0.67, -0.43, -0.21,  0.00,  0.21,  0.43,  0.67, 0.97, 1.38),
  c(-Inf, -1.43, -1.02, -0.74, -0.50, -0.29, -0.10,  0.10,  0.29,  0.50, 0.74, 1.02, 1.43),
  c(-Inf, -1.47, -1.07, -0.79, -0.57, -0.37, -0.18,  0.00,  0.18,  0.37, 0.57, 0.79, 1.07, 1.47),
  c(-Inf, -1.50, -1.11, -0.84, -0.62, -0.43, -0.25, -0.08,  0.08,  0.25, 0.43, 0.62, 0.84, 1.11, 1.5),
  c(-Inf, -1.53, -1.15, -0.89, -0.67, -0.49, -0.32, -0.16,  0.00,  0.16, 0.32, 0.49, 0.67, 0.89, 1.15, 1.53),
  c(-Inf, -1.56, -1.19, -0.93, -0.72, -0.54, -0.38, -0.22, -0.07,  0.07, 0.22, 0.38, 0.54, 0.72, 0.93, 1.19, 1.56),
  c(-Inf, -1.59, -1.22, -0.97, -0.76, -0.59, -0.43, -0.28, -0.14,  0.00, 0.14, 0.28, 0.43, 0.59, 0.76, 0.97, 1.22, 1.59),
  c(-Inf, -1.62, -1.25, -1.00, -0.80, -0.63, -0.48, -0.34, -0.20, -0.07, 0.07, 0.20, 0.34, 0.48, 0.63, 0.80, 1.00, 1.25, 1.62),
  c(-Inf, -1.64, -1.28, -1.04, -0.84, -0.67, -0.52, -0.39, -0.25, -0.13, 0.00, 0.13, 0.25, 0.39, 0.52, 0.67, 0.84, 1.04, 1.28, 1.64),
 )
}

##
## compute distance matrix for the alphabet size specified
##
distance_matrix <- function (alphabet_size){
 if(alphabet_size>1 && alphabet_size<20){
  cutlines <- alphabet2cut(alphabet_size)[2:alphabet_size]
  distance_matrix <- matrix(rep(0, alphabet_size*alphabet_size), byrow=T, nrow=alphabet_size, ncol=alphabet_size)
  i=1
  while(i <= alphabet_size){
    # the min_dist for adjacent symbols are 0, so we start with i+2
    j=i+2;
    while(j <= alphabet_size){
      # square the distance now for future use
      distance_matrix[i,j]=(cutlines[i]-cutlines[j-1])*(cutlines[i]-cutlines[j-1])
      # the distance matrix is symmetric
      distance_matrix[j,i] = distance_matrix[i,j]
      j=j+1;
    }
    i=i+1;
  }
  distance_matrix
 }
}

##
## Converts the specified resolution into the cut points
##
num2letter <- function(num){
  letters <- c("a",  "b",  "c",  "d",  "e",
               "f",  "g",  "h",  "i",  "j",
               "k",  "l",  "m",  "n",  "o",
               "p",  "q",  "r",  "s",  "t",
               "u",  "v",  "w",  "x",  "y",  "z")
  letters[num]
}


##
## Converts the timeseries into string
##
ts2string <- function(ts, aSize){
 cut_points <- alphabet2cut(aSize)
 res <- rep(0, ncol(ts))
 for(i in 1:ncol(ts)){
  res[i] = length(cut_points[cut_points<=ts[i]])
 }
 num2letter(res)
}

ts2dumb <- function(ts, aSize){
 cut_points <- alphabet2cut(aSize)
 res <- rep(0, ncol(ts))
 for(i in 1:ncol(ts)){
  res[i] = length(cut_points[cut_points<=ts[i]])
 }
 return(res)
}

ts2num <- function(ts, aSize){
 cut_points <- alphabet2cut(aSize)
 res <- rep(0, ncol(ts))
 for(i in 1:ncol(ts)){
  res[i] = length(cut_points[cut_points<=ts[i]])
 }
 as.vector(res)
}

##
## compute distance between strings
##
min_dist <- function(str1, str2, alphabet_size, compression_ratio){
 if(length(str1) != length(str2)){
  stop("error: the strings must have equal length");
 }else{    
  if(any(str1 > alphabet_size) | any(str2 > alphabet_size)){
   stop('error: some symbol(s) in the string(s) exceed(s) the alphabet size!');
  }else{
    dist_table <- distance_matrix(alphabet_size);
    dist <- 0;
    dist = sqrt(compression_ratio * sum(diag(dist_table[str1,str2])));
  }
 }
}

# str1 = mtf[i:(i+len.motif-1)]
# str2 = mtf[j:(j+len.motif-1)]

D.ist <- function(str1, str2){
  x = str1
  y = str2
  str1 = paste0(str1,collapse = "")
  str2 = paste0(str2,collapse = "")
  if(length(str1) != length(str2)) stop("error: the strings must have equal length")
  dist_table <- distance_matrix(9);
  dist <- 0;
  dist = sqrt(sum(diag(dist_table[x,y])))
}


##
## compute the Euclidean distance between the set of points
##
euclidean <- function(x, y){
 as.numeric(dist( rbind(x,y) ))
}

dis.cord = function(mtf,len.motif){
  bsf = 0
  bsf.loc = NA
  
  for(i in 1:(length(mtf)-len.motif+1)){
    nnd = Inf
    for(j in 1:(length(mtf)-len.motif+1)){
      if(abs(i-j) >= len.motif){
        tmp = D.ist(mtf[i:(i+len.motif-1)],mtf[j:(j+len.motif-1)])
        if(tmp < nnd){
          nnd = tmp
        }
      }
    }
    if(nnd > bsf){
      bsf = nnd
      bsf.loc = i
    }
  }
  return(bsf.loc)
}


find_plc = function(a,in.b){
  a = paste(as.character(a),sep="",collapse="")
  in.b = paste(as.character(in.b),sep="",collapse="")
  plc = str_locate_all(in.b,a)
  return(as.matrix(data.frame(plc)))
}

Motif2File = function(tmp, len.motif){
  
  split4me = function(t.p){
    t.v = c(rep(0,len.motif))
    for(i in 1:(length(t.p)/len.motif)){
      t.v = cbind(t.v,t.p[c((1:len.motif)+((i-1)*len.motif))])
    }
    t.v = t.v[,-1]
    return(t.v)
  }
  
  Un.Iq = unique(tmp$col.vec)
  Un.Iq = Un.Iq[-which(Un.Iq == "black")]
  t.v = c(rep(0,len.motif))
  for(i in 1:length(Un.Iq)){
    t.p = which(tmp$col.vec == Un.Iq[i])
    if(length(t.p) > len.motif)t.p = split4me(t.p)
    t.v = cbind(t.v,t.p)
  }
  t.v = t.v[,-1]
  t.v = t(t.v)
  tv.int = vector()
  
  if(length(Un.Iq) == 1){
    tv.int = as.numeric(paste(tmp[matrix(t(t.v))],collapse=""))
  }else{
    t.v = t.v[order(t.v[,1]),]
    for(i in 1:length(Un.Iq)){
      tv.int[i] = as.numeric(paste(tmp[matrix(t(t.v[i,]))],collapse=""))
    }
  }

  flat = as.numeric(paste(rep("1",len.motif),collapse = ""))
  idx = which(tv.int %% flat == 0)
  tv.int = tv.int[-idx]

  print(tv.int)
  
  write.table(tv.int, file = file.path(getwd(),"Motif-Int-List.tsv"), 
        sep = '\t', append = T,  col.names = F)
}

MagieNoire = function(mtf,len.motif){
  
  dc.mat = c(rep(0,len.motif))
  
  S.Needs = function(dc.mat,len.motif,mtf){
    dc.mat = na.omit(dc.mat)
    dc.mat = dc.mat[-1,]
    segs = find_plc(dc.mat,mtf)
    col.vec = c(rep("black",length(mtf)))
    col.vec[segs[1,1]:segs[1,2]] = brewer.pal(3, "Accent")[1]
    total = data.frame(cbind(mtf,col.vec))
    total$mtf = as.numeric(total$mtf)
    return(total)
  }
  
  motif.scan = function(mtf,len.motif){
    dc = dis.cord(mtf,len.motif)
    if(!is.na(dc)){
      dc.mat <<- rbind(dc.mat,mtf[dc:(dc+len.motif-1)])
      set1 = mtf[1:(dc-1)]
      set2 = mtf[(dc+len.motif):length(mtf)]
      if(length(set1) > (len.motif-1)) motif.scan(set1,len.motif)
      if(length(set2) > (len.motif-1)) motif.scan(set2,len.motif)
    }
  }
  
  motif.scan(mtf,len.motif)
  
  if(length(dc.mat) == len.motif)return(c(0))
  if(length(dc.mat) == 2 * len.motif)return(S.Needs(dc.mat,len.motif,mtf))
  
  dc.mat = na.omit(dc.mat)
  dc.mat = dc.mat[-1,]
  segs = matrix(ncol = 2)
  
  for(i in 1:nrow(dc.mat)){segs = rbind(segs,find_plc(dc.mat[i,],mtf))}
  
  segs = segs[-1,]
  
  col.vec = c(rep("black",length(mtf)))
  
  for(i in 1:nrow(segs)){
    col.vec[segs[i,1]:segs[i,2]] = brewer.pal(5, "Greys")[(i%%5)+1]
  }
  
  total = data.frame(cbind(mtf,col.vec))
  total$mtf = as.numeric(total$mtf)
  
  return(total)
}

find_plc = function(a,in.b){
  a = paste(as.character(a),sep="",collapse="")
  in.b = paste(as.character(in.b),sep="",collapse="")
  plc = str_locate_all(in.b,a)
  return(as.matrix(data.frame(plc)))
}

Motif2File = function(tmp, len.motif){
  
  split4me = function(t.p){
    t.v = c(rep(0,len.motif))
    for(i in 1:(length(t.p)/len.motif)){
      t.v = cbind(t.v,t.p[c((1:len.motif)+((i-1)*len.motif))])
    }
    t.v = t.v[,-1]
    return(t.v)
  }
  
  Un.Iq = unique(tmp$col.vec)
  Un.Iq = Un.Iq[-which(Un.Iq == "black")]
  t.v = c(rep(0,len.motif))
  for(i in 1:length(Un.Iq)){
    t.p = which(tmp$col.vec == Un.Iq[i])
    if(length(t.p) > len.motif)t.p = split4me(t.p)
    t.v = cbind(t.v,t.p)
  }
  t.v = t.v[,-1]
  t.v = t(t.v)
  tv.int = vector()
  
  if(length(Un.Iq) == 1){
    tv.int = as.numeric(paste(tmp[matrix(t(t.v))],collapse=""))
  }else{
    t.v = t.v[order(t.v[,1]),]
    for(i in 1:length(Un.Iq)){
      tv.int[i] = as.numeric(paste(tmp[matrix(t(t.v[i,]))],collapse=""))
    }
  }

  flat = as.numeric(paste(rep("1",len.motif),collapse = ""))
  idx = which(as.numeric(tv.int) %% flat == 0)
  tv.int = tv.int[-idx]
  
  write.table(tv.int, file = file.path(getwd(),"motif-int-list.tsv"), 
        sep = '\t', append = T,  col.names = F)
}

MagieNoire = function(yy,len.motif){
  
  dc.mat = c(rep(0,len.motif))
  
  S.Needs = function(dc.mat,len.motif,yy){
    dc.mat = na.omit(dc.mat)
    dc.mat = dc.mat[-1,]
    segs = find_plc(dc.mat,yy)
    col.vec = c(rep("black",length(yy)))
    col.vec[segs[1,1]:segs[1,2]] = "gray"
    total = data.frame(cbind(yy,col.vec))
    total$yy = as.numeric(total$yy)
    return(total)
  }
  
  motif.scan = function(yy,len.motif){
    dc = dis.cord(yy,len.motif)
    if(!is.na(dc)){
      dc.mat <<- rbind(dc.mat,yy[dc:(dc+len.motif-1)])
      set1 = yy[1:(dc-1)]
      set2 = yy[(dc+len.motif):length(yy)]
      if(length(set1) > (len.motif-1)) motif.scan(set1,len.motif)
      if(length(set2) > (len.motif-1)) motif.scan(set2,len.motif)
    }
  }
  
  motif.scan(yy,len.motif)
  
  if(length(dc.mat) == len.motif)return(c(0))
  if(length(dc.mat) == 2 * len.motif)return(S.Needs(dc.mat,len.motif,yy))
  
  dc.mat = na.omit(dc.mat)
  dc.mat = dc.mat[-1,]
  segs = matrix(ncol = 2)
  
  for(i in 1:nrow(dc.mat)){segs = rbind(segs,find_plc(dc.mat[i,],yy))}
  
  segs = segs[-1,]
  
  col.vec = c(rep("black",length(yy)))
  
  for(i in 1:nrow(segs)){
    col.vec[segs[i,1]:segs[i,2]] = rainbow(nrow(segs))[i]
  }
  
  total = data.frame(cbind(yy,col.vec))
  total$yy = as.numeric(total$yy)
  
  return(total)
}
