#' Orders data for each year
#' @description Orders data for each year in "yearlist" according to the max-min distance
#' ordering
#'
#' @param unordered_data Un-ordered data to order. Must contain a "years"
#' column. Data for each year will be ordered separately.
#' @param yearlist List of years to order.
#' @return Returns a dataframe the same size as "unordered_data"
#' @importFrom GPvecchia order_maxmin_exact
#' @export
order_yeardata=function(unordered_data,yearlist=2007:2016){
  ordered_augdata=unordered_data
  for(year in yearlist){
    yeardata=unordered_data[unordered_data$years==year,]
    myord = GPvecchia::order_maxmin_exact(cbind(yeardata$lon_rad,yeardata$lat_rad))
    cut = min(dim(yeardata)[1], 9)
    myord=c(myord[1], myord[-seq(1,cut)], myord[2:cut])
    yeardata_ordered=yeardata[myord,]
    ordered_augdata[unordered_data$years==year,]=yeardata_ordered
  }
  return(ordered_augdata)
}
