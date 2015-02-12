#' Function to Analyze Survey Scales for Possible Reduction
#'
#' This function provides a series of statistics to determine minimum scale length with properties similar to a full scale.
#' @param full.scale A dataset containing all items in the full scale, rescaled to range between 0 and 1; must not include NAs
#' @param retest Indicator for whether the scale was administered twice; in this case, responses from both administrations should be included in full.scale; default is FALSE
#' @param factors the number of factors you want to analyze in the principal components section; default is 2
#' @return \strong{abbreviated.scales} creates mean values of all subscales (including one-item scales and the full original scale)
#' @return \strong{abbreviated.scales.retest} if retest=TRUE, creates mean values of all subscales from second administration of the scale
#' @return \strong{statistics} provides information on each subscale, including which items the scale contains, its length, and all possible statistics calculated using this function (e.g., Cronbach's alpha cannot be calculated for the one-item scales, so these statistics are not provided here)
#' @return \strong{item.analysis} summarizes information on subscales of various length (see tables X and Y in Loose, et. al. (2015))
#' @return \strong{question.analysis} summarizes information on subscales containing each item from the original full scale (see tables X and Y in Loose, et. al. (2015))
#' @importFrom psy cronbach
#' @importFrom psych corr.test
#' @export
#' @author Krista Loose \email{loosek@@mit.edu}
#' @references Loose, Krista, Yue Hou, Adam Berinsky, and Cindy Kam.  2015.  "Strategies for Scale Reduction."  Working Paper.


ScaleReduce = function(full.scale, retest=FALSE, factors=2){

  requireNamespace(psych, quietly = TRUE)
  requireNamespace(psy, quietly = TRUE)
  ## determine length of full scale
  if(retest == TRUE){
    n = ncol(full.scale)/2
  }
  else{
    n = ncol(full.scale)
  }

  ## identify all possible subscales of original scale
  combinations = as.data.frame(matrix(NA, nrow=0, ncol=n))
  for(i in 1:n){
    cols = t(combn(n, i))
    for(j in 1:nrow(cols)){
      this.combination = cols[j,]
      no.NAs = n - length(this.combination)
      combinations[nrow(combinations)+1,] = c(this.combination, rep(NA, no.NAs))
    }
  }
  m = nrow(combinations)

  ## create values of scale for each respondent in data
  abbreviated.scales = matrix(NA, nrow=nrow(full.scale), ncol=m)
  for(i in 1:m){
    scale.columns = c(as.matrix(combinations[i, !is.na(combinations[i,])]))
    for(j in 1:nrow(full.scale)){
      abbreviated.scales[j,i] = sum(full.scale[j, scale.columns], na.rm=TRUE)/length(scale.columns)
    }
  }

  ## create values of scale for each respondent in re-test setting
  if(retest == TRUE){
    abbreviated.scales.retest = matrix(NA, nrow=nrow(full.scale), ncol=m)
    for(i in 1:m){
      scale.columns = c(as.matrix(combinations[i, !is.na(combinations[i,])] + n))
      for(j in 1:nrow(full.scale)){
        if(sum(is.na(full.scale[j, scale.columns]))==7) {abbreviated.scales.retest[j,] = NA}
        else{
          abbreviated.scales.retest[j,i] = sum(full.scale[j, scale.columns], na.rm=TRUE)/length(scale.columns)
        }
      }
    }
  }

  ## data.matrix will contain statistics on each possible subscale
  data.matrix = as.data.frame(matrix(NA, nrow=m, ncol=n))
  colnames(data.matrix) = colnames(full.scale)[1:n]
  for(i in 1:n){
    for(j in 1:m){
      data.matrix[j,i] = ifelse(i %in% combinations[j,], 1, 0)
    }
  }

  ## calculates the number of items in each subscale
  data.matrix$items = NA
  for(i in 1:nrow(data.matrix)){
    data.matrix$items[i] = sum(data.matrix[i,], na.rm=TRUE)
  }

  ## correlation with full scale
  data.matrix$correlations = NA
  for(i in 1:(m-1)){
    data.matrix$correlations[i] = corr.test(as.matrix(abbreviated.scales[,i]), as.matrix(abbreviated.scales[,m]))$r
  }

  ## Cronbach's alpha coefficient
  data.matrix$alphas = NA
  for(i in (n+1):m){
    alpha.columns = c(as.matrix(combinations[i, !is.na(combinations[i,])]))
    data.matrix$alphas[i] = cronbach(full.scale[, alpha.columns])$alpha
  }

  ## principal factor analysis
  ## code based on Rencher (2002, pp. 448-450)
  new.cols = (ncol(data.matrix)+1):(ncol(data.matrix)+(factors))
  data.matrix[,new.cols] = NA
  colnames(data.matrix)[new.cols] = paste("eigen", 1:factors, sep="")

  for(i in (n+1):m){
    ev.columns = c(as.matrix(combinations[i, !is.na(combinations[i,])]))
    R = cor(full.scale[,ev.columns], use="complete.obs")
    R1 = solve(R)
    for(j in 1:length(ev.columns)){
      R[j,j] = R[j,j] = 1 - (1/R1[j,j])
    }
    evs = eigen(R)$values
    for(j in 1:length(new.cols)){
      data.matrix[i, new.cols[j]] = evs[j]
    }
  }

  ## classification rate
  data.matrix$classification = NA
  high_full.scale = ifelse(abbreviated.scales[,m] > median(abbreviated.scales[,m]), 1, 0)
  for(i in 1:(m-1)){
    high_subscale = ifelse(abbreviated.scales[,i] > median(abbreviated.scales[,i]), 1, 0)
    data.matrix$classification[i] = 1 - (sum(high_full.scale != high_subscale)/length(high_full.scale))
  }

  ## test-retest statistics
  if(retest == TRUE){
    data.matrix$retest = NA
    for(i in 1:m){
      data.matrix$retest[i] = cor(abbreviated.scales[,i], abbreviated.scales.retest[,i], use="pairwise.complete.obs")
    }
  }

  ## create summary tables
  if(retest == TRUE){
  	number_of_metrics = 4 + factors
  }
  else{
  	number_of_metrics = 3 + factors
  }

  item.analysis = as.data.frame(matrix(NA, nrow=n, ncol=number_of_metrics*3))
  question.analysis = as.data.frame(matrix(NA, nrow=n, ncol=number_of_metrics*3))
  stats.function = function(x){
  	c(mean(x, na.rm=TRUE), quantile(x, c(0.025, 0.975), na.rm=TRUE))
  }

  for(j in 1:number_of_metrics){
  	metric.data = data.matrix[,(j+n+1)]
  	metric.cols = (j*3-2):(j*3)
  	metric.name = colnames(data.matrix)[(j+n+1)]
  	colnames(item.analysis)[metric.cols] = paste(metric.name, c("mean","lb","ub"), sep="_")
  	colnames(question.analysis)[metric.cols] = paste(metric.name, c("mean","lb","ub"), sep="_")
  	for(i in 1:n){
  		item.subset = metric.data[which(data.matrix$items==i)]
  		item.analysis[i,metric.cols] = stats.function(item.subset)
  		question.subset = metric.data[which(data.matrix[,i]==1)]
  		question.analysis[i,metric.cols] = stats.function(question.subset)
  	}
  }

  rownames(item.analysis) = paste("items", 1:n, sep="_")

  ## item scale correlation
  old.cols = ncol(data.matrix)
  data.matrix[,(old.cols+1):(old.cols+n)] = NA
  colnames(data.matrix)[(old.cols+1):(old.cols+n)] = paste("item_total", colnames(data.matrix)[1:n], sep="_")

  for(i in (n+1):m){
  	total.columns = c(as.matrix(combinations[i, !is.na(combinations[i,])]))
  	total.stats = alpha(full.scale[,total.columns])$item.stats
  	for(j in 1:length(total.columns)){
  		data.matrix[i,(total.columns[j]+old.cols)] = total.stats$r.drop[j]
  	}
  }

  reduced.matrix = data.matrix[,(n+2+number_of_metrics):(ncol(data.matrix))]
  item_totals = apply(reduced.matrix[,-127], 2, stats.function)
  colnames(item_totals) = colnames(data.matrix)[1:n]
  rownames(item_totals) = paste("item_total", c("mean","lb","ub"), sep="_")
  question.analysis = cbind(question.analysis, t(item_totals))

  ## return all data
  list.of.data = list(statistics = data.matrix,
                      abbreviated.scales = abbreviated.scales,
                      item.analysis = item.analysis,
                      question.analysis = question.analysis)
  if(retest == TRUE){
    list.of.data = c(list.of.data, abbreviated.scales.retest)
  }

  return(list.of.data)
}
