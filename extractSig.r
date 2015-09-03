
# --- header ---

require(signalextraction)

require(quantmod)

require(googleVis)

require(RColorBrewer)

source("alt_dll.r")

# --- foo ---

Sys.setenv(TZ = 'Asia/Jerusalem')

applyFit = function(ts,fit){
  last = index(fitted)[length(fitted)]
  out = ts[index(ts) > last]
  coredata(ts[index(ts) > last]) = outsamp(fit, out, sequel=T)
  return(ts)
}

extractSig = function(ts = NULL ,fit = NULL){  
  if(is.null(fit)){
    if(!is.null(ts)){
      fit = dfa(coredata(ts), quart = FALSE, d = 0, pb = 1/6.5, sb = 1/6,
                tpfilter = TRUE, lambda = 5, expweight = 1.5,
                pbd = 1.02, limamp = 3, i2 = TRUE, n.loops = 10,
                verbose = 1)
      fitted = zoo(fit$xf, order.by = index(ts))
    }else{
      fitted = ts
    }
    
  }else{
    fitted = applyFit(ts,fit)
  }
  
  ones = which(rle(sign(diff(coredata(fitted))))$lengths < 5)
  idx = rle(sign(diff(coredata(fitted))))$lengths[-ones]
  idx = cumsum(idx)
  idx = idx + 1
  sig = -sign(coredata(fitted[idx]) - coredata(fitted[idx-3]))
  pref = sig[-length(sig)] * diff(ts[idx])
  idx.gain = which(pref >= 0)
  idx.loss = which(pref < 0)
  stop.loss = quantile(pref[idx.loss],0.25)
  pref[idx.loss] = max(pref[idx.loss],stop.loss)
  plr = abs(mean(pref[idx.gain])/mean(pref[idx.loss]))
  
  return(list(ts = ts,fit = fit, fitted = fitted, idx = idx, sig = sig, 
              pref = pref, idx.gain = idx.gain, idx.loss = idx.loss,
              stop.loss = stop.loss, plr = plr))
}


# --- main ---

ticker = "3"

switch(ticker,
       '1' = {symb="T25"; exchange="TLV"},
       '2' = {symb=".IXIC"; exchange="INDEXNASDAQ"},
       '3' = {symb=".INX"; exchange="INDEXSP"},
       '4' = {symb="AAPL"; exchange="NASDAQ"},
       '5' = {symb="SX5E"; exchange="INDEXSTOXX"})

wrap = function(options){
  options = cbind(options)
  for(i in 1:length(options)){
    options[i] = "'" %+% options[i] %+% "'"
  }
  print(options)
  opt = "[" %+% options[1]
  for(i in 2:length(options)){
    opt = opt %+% "," %+% options[i]
  }
  opt = opt %+% "]"
  return(opt)
}

ts = last2last(0,symb,exchange)

out = extractSig(ts,NULL)

new_ts = extractSig(ts,out$fit)

stocks = data.frame(cbind(new_ts$ts,new_ts$fitted))

stocks = cbind(format(index(ts),"%T"),stocks)

colnames(stocks) = c("time", "original", "fitted")

colors = c("#f7f7f7", "#f1a340", "#998ec3")

colors = wrap(colors[2:3])

options = list(height=300, backgroundColor = colors[1], colors = colors)

chart = gvisLineChart(stocks, options = options)

chart$html$footer = NULL

# new_ts$idx[which(new_ts$idx > last(old_idx)[1]:length(new_ts$idx)]

plot(chart) # tag = "chart"
