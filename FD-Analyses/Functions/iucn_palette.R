iucn_palettes <- list(
  CO= "#000000",
  CR= "#FF0000",
  EN= "#FFA500",
  VU= "#FFFF00",
  NT= "#ADFF2F",
  LC= "#008000",
  DD= "#808080",
  NE= "#FFFFFF",
  All =  c("#FFFFFF", "#808080", "#008000", "#ADFF2F", "#FFFF00", "#FFA500", "#FF0000", "#000000")
)

iucn_palette <- function(category="All", exclude=NA){
  if(length(category)>0){
    if(length(category)==1){
      pal <- iucn_palettes[[category]]
    } else if(length(category)>1){
      pal <- iucn_palettes[grepl(paste(category, collapse="|"), names(iucn_palettes))]
      pal <- as.character(unlist(pal))
    }
  }
  if(is.na(exclude[1])==F){
    exclude_codes <- iucn_palette(category=exclude)
    pal <- pal[!pal %in% exclude_codes]
  }
  if (is.null(pal)){
    stop("Palette not found.")
  } else {
    #return(pal)
    type <- "discrete"
    n <- length(pal)
    out <- switch(type,
                  discrete = pal[1:n])
    
    structure(out, class = "palette", name = category)
  }
}