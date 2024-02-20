
String <- function(x="") {
  x <- as.character(paste(x, collapse=""))
  class(x) <- c("String","character")
  return(x)
}

"[.String" <- function(x,i,j,...,drop=TRUE) {
  unlist(strsplit(x,""))[i]
}
"[<-.String" <- function(x,i,j,...,value) {
  tmp <- x[]
  tmp[i] <- String(value)
  x <- String(tmp)
  x
}