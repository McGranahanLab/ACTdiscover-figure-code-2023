#wGII (Definition from Dissertation_Endersfelder): The chromosome based GII score for all autosomal chromosomes
#c=1,...,22 with segments 1,...,k corresponding copy numbers x1,...,xk and corresponding
#lengths l1,...,lk was defined as gc = (sum(li where xi<tl)+sum(li where xi>tg))/sum(li),
#where tl and tg are thresholds for copy number loss or gain (here tl = tg = ploidy).
#The wGII score for a sample was then defined as the average over all 22 (autosomal) chromosome based
#GII scores --> wGII = sum(gc)/22.

calculate.wGIIscore = function(cn.data, pro.Region = FALSE) {
  
  out.table = c()
  for(Sample in unique(cn.data$Patient)){
    
    data.sample = cn.data[cn.data$Patient == Sample,]
    
    #wGII score for every region if there is multiple region data for one sample
    if (pro.Region) {
      for(Region in unique(data.sample$region)){
        data.region = data.sample[data.sample$region == Region,]
        
        #GII score for every autosomal chromosome
        wGII = 0
        for(i in 1:22) {
          data.chr = data.region[data.region$chr == i,] 
          seg.gain = abs(data.chr$endpos[data.chr$cnTotal > round(data.chr$Ploidy)] - data.chr$startpos[data.chr$cnTotal > round(data.chr$Ploidy)])
          seg.loss = abs(data.chr$endpos[data.chr$cnTotal < round(data.chr$Ploidy)] - data.chr$startpos[data.chr$cnTotal < round(data.chr$Ploidy)])
          seg.total = abs(data.chr$endpos - data.chr$startpos)
          
          GII = (sum(seg.gain) + sum(seg.loss)) / sum(seg.total)
          wGII = wGII + GII
          
        }
        
        #weighted GII score
        wGII = wGII/22
        
        out       = c(data.region$Patient[1], data.region$region[1], wGII)
        out.table = rbind(out.table, out)
        colnames(out.table) = c('Patient', 'region', 'wGII')
        
      }
      
    }
    
    #wGII score for no multiple region data
    if (!pro.Region) {
      
      #GII score for every autosomal chromosome
      wGII = 0
      for (i in 1:22) {
        data.chr = data.sample[data.sample$chr == i,] 
        seg.gain = abs(data.chr$endpos[data.chr$cnTotal > round(data.chr$Ploidy)] - data.chr$startpos[data.chr$cnTotal > round(data.chr$Ploidy)])
        seg.loss = abs(data.chr$endpos[data.chr$cnTotal < round(data.chr$Ploidy)] - data.chr$startpos[data.chr$cnTotal < round(data.chr$Ploidy)])
        seg.total = abs(data.chr$endpos - data.chr$startpos)
        
        GII = (sum(seg.gain) + sum(seg.loss)) / sum(seg.total)
        wGII = wGII + GII
        
      }
      
      #weighted GII score
      wGII = wGII / 22
      
      out       = c(data.sample$Patient[1], wGII)
      out.table = rbind(out.table, out)
      colnames(out.table) = c('Patient', 'wGII')
      
    }
  }
  
  return(out.table)
  
}
