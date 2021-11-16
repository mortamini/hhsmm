#' @importFrom graphics axTicks legend lines par rect
#'
.add.states <- function(states,ht=0,greyscale=FALSE,leg=NA,J=length(unique(states)),time.scale=24,shift=0) {  
  J = length(unique(states))  
  
  if(greyscale) cols=c("#FFFFFF" ,"#F0F0F0" ,"#D9D9D9", "#BDBDBD" ,"#969696", "#737373", "#525252", "#252525")
  else cols = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3") #kind regards to RBrewerPal for these values  
  
  hats = rle(states)
  hats = list(intervals=cumsum(c(0,hats$lengths))/time.scale+shift,state=hats$values)
  for(ii in 1:length(hats$state)) 
    if(greyscale)  rect(hats$intervals[ii],ht,hats$intervals[ii+1],ht+(axTicks(2)[2]-axTicks(2)[1])/5,col=cols[hats$state[ii]],border=1)
  else rect(hats$intervals[ii],ht,hats$intervals[ii+1],ht+(axTicks(2)[2]-axTicks(2)[1])/5,col=cols[hats$state[ii]],border=cols[hats$state[ii]])
  if(any(!is.na(leg))) legend("topleft",legend=leg,fill=cols,bg="white")
}
