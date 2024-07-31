## Import R utilities from Jacques van Helden
export.plot <- function (file.prefix="PlotExport",
                         export.formats="pdf", # supported: postscript, jpg, png, bmp, pdf
                         width=11, # in inches
                         height=8, # in inches
                         horizontal=T,
                         ... ## Additional parameters are passed to the export method
) {
  
  ppi <- 72
  file.ext <- c(
    postscript = "ps",
    pdf = "pdf",
    ps = "ps",
    eps = "eps",
    jpeg="jpg",
    jpg="jpg",
    bmp="bmp",
    png="png",
    svg="svg",
    tiff="tiff")
  for (f in export.formats) {
    from.dev <- dev.cur();
    
    file.name <- paste(file.prefix,file.ext[f], sep=".")
    
    if ((f == "postscript") || (f == "ps")) {
      postscript(file.name,paper="special",width=width,height=height,horizontal=horizontal, ...)
    } else if (f == "eps") {
      postscript(file.name,paper="special",width=width,height=height,horizontal=horizontal,onefile=F, ...)
    } else if (f == "pdf") {
      pdf(file.name, paper="special",width=width,height=height, ...)
    } else if ((f == "jpg") || (f == "jpeg")) {
      jpeg(file.name,width=(width*ppi),height=(height*ppi),quality=100, ...)
    } else if (f == "png") {
      png(file.name,width=width*ppi,height=height*ppi, ...)
    } else if (f == "bmp") {
      bitmap(file.name,width=width*ppi,height=height*ppi, ...)
    } else if (f == "svg") {
      svg(file.name,width=width*ppi,height=height*ppi, ...)
    } else if (f == "tiff") {
      #tiff(filename = "Rplot%03d.tiff", width = 480, height = 480, units = "px", pointsize = 12, compression = c("none", "rle", "lzw", "jpeg", "zip"), bg = "white", res = NA,  ..., type = c("cairo", "Xlib", "quartz"), antialias)
      tiff(file.name,width=width*ppi,height=height*ppi, compression = 'none', ...)
    }
    else {
      print(paste("Error: format ", f, " is not supported", sep=""))
      return()
    }
    to.dev <- dev.cur()
    dev.set(which=from.dev)
    dev.copy(which=to.dev)
    dev.set(which=to.dev)
    dev.off()
    dev.set(which=from.dev) ## This is required because dev.off() returns to the first, not the last, device
  }
}
