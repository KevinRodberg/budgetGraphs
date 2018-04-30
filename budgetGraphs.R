list.of.packages <-c("reshape2","ggplot2","grDevices","rasterVis","maptools","future","listenv","classInt",
                     "rgdal","rgeos","tcltk2")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

#---------------------------------------------------------------------------------------------------
# The order of the following library functions are very important
# due to package:ggplt2 and package:latticeExtra fighting over the 'layer' object  
# and package:raster and package:future fighting over the 'value' object
#---------------------------------------------------------------------------------------------------

library(reshape2)
library(ggplot2)
library(grDevices)
library(rasterVis)
library(maptools)
library(future)
library(listenv)
library(classInt)

ggplotColours <- function(n = 6, h = c(0, 360) + 15) {
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

options(warn = -1)
options(scipen = 100000)

source_github <- function(u) {
  # load package
  require(RCurl)
  
  # read script lines from website
  script <- getURL(u, ssl.verifypeer = FALSE)
  script<-strsplit(script,"\r\n")
  # parase lines and evealuate in the global environement
  eval(parse(text = unlist(script)))
}

source_https<-function(u, unlink.tmp.certs = FALSE) {
  # load package
  require(RCurl)
  
  # read script lines from website using a security certificate
  if(!file.exists("cacert.pem")) download.file(url="http://curl.haxx.se/ca/cacert.pem", destfile = "cacert.pem")
  script <- getURL(u, followlocation = TRUE, cainfo = "cacert.pem")
  if(unlink.tmp.certs) unlink("cacert.pem")
  
  # parase lines and evealuate in the global environement
  eval(parse(text = script), envir= .GlobalEnv)
}


source_github("https://raw.githubusercontent.com/KevinRodberg/ReusableFunctions/master/rasterFuncs.R")
source_github("https://raw.githubusercontent.com/KevinRodberg/ReusableFunctions/master/gisFuncs.R")
source_github("https://raw.githubusercontent.com/KevinRodberg/ReusableFunctions/master/modflowDataFuncs.R")

#=================================================================
# Beginning of Script

MFmodel.Params <- defineMFmodel()
if (!exists("MFmodel")) {
  MFmodel <- 'x'
}

makePlots <- function(i, BasinNames, data, output) {
  wraps <-    strwrap(BasinNames$Name_HUC12[i],width = 80,initial = "",prefix = "\n",simplify = TRUE)
  # barplot(data[,i], main=paste0(wraps),sub=paste0('Basin Id: ',i),beside=T, ylim=c(0,80))
  oneRec <- data[, i]
  plotdat <- melt(oneRec)
  plotfile = paste0(outDir, "/Basin", i, "_Yr", year, ".jpg")
  labelNames = rownames(plotdat)
  
  if (output == "save") {

    ggplot(plotdat, aes(x = rownames(plotdat),y = value,fill = labelNames)) +
      geom_bar(stat = "identity") + theme(legend.position = "none") +
      theme(axis.text.x=element_text(angle=45,hjust=1)) +
      scale_x_discrete(limits=rownames(plotdat)) +
      ylim(-1, 70) + labs(x ="Budgets", fill = "Budget") +
      ggtitle(paste(unlist(wraps), collapse = ""), subtitle = paste0('Basin Id: ', i))
    ggsave(plotfile,width = 10,height = 7.5,units = "in",dpi = 300,device = "jpeg")
  }
  else{
    print(ggplot(plotdat, aes(x = rownames(plotdat),y = value,fill = labelNames)) +
            geom_bar(stat = "identity") + theme(legend.position = "none") +
            theme(axis.text.x=element_text(angle=45,hjust=1)) +
            scale_x_discrete(limits=rownames(plotdat)) +
            ylim(-1,70) + labs(x ="Budgets", fill = "Budget") +
            ggtitle(paste(unlist(wraps), collapse = ""), subtitle = paste0('Basin Id: ', i))
    )
  }
  
  return(0)
}

WMDbnd <- buildWMDbnd(HARNSP17ft)

plotRaster <- function(ras, colorRamp, atValues, rasName) {
  
  dir.create(file.path(outDir, "plots"), showWarnings = TRUE)
  pngFile = paste(outDir, "/plots/", rasName, ".png", sep = "")
  
  png(file = pngFile,width = 3000,height = 2400,units = "px",res = 300  )
  print(levelplot(ras,contour = FALSE,col.regions = colorRamp,at =  unique(atValues),
      main = rasName,scales = list(x = list(rot = 90)))+ 
        layer(sp.polygons(WMDbnd))
  )
  dev.off()
}

tkmessageBox(message = "Click OK and Provide Required I/O Info",
             icon = "warning",type = "ok")


if (!(MFmodel %in% (c('ECFTX', 'NPALM', 'LWCSIM')))) {
  MFmodel <- chooseModel()
  if ((!MFmodel %in% (c('ECFTX', 'NPALM', 'LWCSIM')))) {
    exit('Incorrect model choice')
  }
}
M = MFmodel.Params[MFmodel,]

if (!exists("SP_rng")) {
  SP_rng <- seq(1, M$nsp)
}
SPknt <- length(SP_rng)

# Calculate raster extents
xmax = M$xmin + (M$ncols * M$res)
ymax = M$ymin + (M$nrows * M$res)
modelras <-  raster(    resolution = M$res,
                        xmn = M$xmin,
                        xmx = xmax,
                        ymn = M$ymin,
                        ymx = ymax,
                        crs = HARNSP17ft  )


#Establish working directory for output:
setwd(paste0("//ad.sfwmd.gov/dfsroot/data/wsd/MOD/", MFmodel))
outDir <-  tk_choose.dir(default = getwd(), caption = "Select output directory")

setwd("//whqhpc01p/hpcc_shared/dbandara/CFWI/ECFTX/ET_RCH/ET_RECH_ECFTX")
inDir <-  tk_choose.dir(default = getwd(), caption = "Select input directory for budget files")

DaysPerMonth <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

#createBarCharts <-FALSE
createBarCharts <- TRUE
#onScreenOnly <- FALSE
onScreenOnly <- TRUE

plan(multiprocess)

BckGrndreadBud<-function(filename){
  budgetByYear <-  read.table(filename, header = TRUE, sep = ",")
  return(budgetByYear)
}

#years <- c(1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014)
#years <- c(2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014)
years <- c(2003)
#years <- c(2002,2003,2004)
budList <- listenv(NULL)
for (year in years) {
  print(paste("Reading Budgets in Background for year:", year))
  #filename <-paste0("//whqhpc01p/hpcc_shared/dbandara/CFWI/ECFTX/ET_RCH/ET_RECH_ECFTX/Output/040418_ReduceIRR_DESTO_IR_MTN/transient/budget_",year,".txt")
  filename <- paste0(inDir, "/budget_", year, ".txt")
  budList[[year]] <-future({BckGrndreadBud(filename)})
}
for (year in years){
  print(paste("Processing Budgets for year:", year))
  
  budgetByYear<-as.data.frame(future::value(budList[[year]]) )

  NameFile <-"//whqhpc01p/hpcc_shared/krodberg/CFWI/WatershedNameIDs.csv"
  BasinNames <- read.table(NameFile, header = TRUE, sep = ",")
  
  names(budgetByYear)
  budgetByYear$BASIN[budgetByYear$BASIN == -999] <- NA
  budgetByYear$BASIN[budgetByYear$BASIN == 0] <- NA
  budgetByYear$ERROR <-(budgetByYear$RAIN+budgetByYear$AGI+budgetByYear$LSI)-
    (budgetByYear$RUNOFF+budgetByYear$RECHARGE+budgetByYear$UZET)
  
 # Allterms <- names(budgetByYear)
  Allterms = c("RAIN", "AGI" ,"LSI","RUNOFF","INIABS","RECHARGE","SWC","IRR","UZET","GWETmax","RET","EXTD","DTW","PERV_AREA","IMPERV_AREA","DCIA","ERROR"  )
  meltBudgets <- melt(budgetByYear,
                      id.vars = c("XCOORD", "YCOORD", "ROW", "COL", "BASIN", "BASIN_NAME"),
                      measure.vars = Allterms )
  sumTerms = c("PERV_AREA", "IMPERV_AREA", "DCIA")
  meanTerms = sort(c("RAIN","AGI","LSI","RUNOFF","INIABS","RECHARGE","SWC","IRR","UZET","GWETmax","RET","EXTD","DTW","ERROR"))
  
  if (createBarCharts) {
    pltGrphs <- listenv(NULL)
    basinMeans <-dcast(subset(meltBudgets,variable %in% meanTerms),BASIN~variable,mean,na.rm = TRUE)
    #    basinSums <- dcast(subset(meltBudgets,variable%in%sumTerms) , BASIN~variable , sum,na.rm=TRUE)
    data <- t(basinMeans[, -1])
    colnames(data) <- basinMeans$BASIN
    labelNames = paste(rownames(data), '[in/year]')
    labelNames[labelNames == 'DTW [in/year]'] <- 'DTW [feet]'
    labelNames[labelNames == 'EXTD [in/year]'] <- 'EXTD [feet]'
    rownames(data) <- labelNames
    if (!onScreenOnly) {
      cat("\nCreating BarChart figure files\n")
      ix = 0
      ilen = length(unlist(basinMeans$BASIN))
      for (i in unlist(basinMeans$BASIN)) {
        ix = ix + 1
        cat(paste('\r', format(ix / ilen * 100, digits = 2, nsmall = 2), "%"))
        if (!is.na(i)) {
          pltGrphs[[i]] <- future({
            makePlots(i, BasinNames, data, "save")
          })
        }
      }
    }
    if (onScreenOnly){
      cat("\nCreating BarChart Plots\n")
      ix = 0
      ilen = length(unlist(basinMeans$BASIN))
      for (i in unlist(basinMeans$BASIN)) {
        ix = ix + 1
        cat(paste('\r', format(ix / ilen * 100, digits = 2, nsmall = 2), "%"))
        if (!is.na(i)) {
          x <- makePlots(i, BasinNames, data, "onscreen")
        }
      }
    }
  }
  
  printingOn = TRUE
  #  printingOn = FALSE
  
  if (printingOn) {
    cat(paste0("\nRaster plots are being created in background processes: \n",
        "//ad.sfwmd.gov/dfsroot/data/wsd/MOD/",MFmodel,"/WaterBudgets/plots \n" ) )
  } else{
    cat(paste0("\nRasters are being processed\n"))
  }
  
  # Initialize special list to receive 'futures' return values
  pltMaps <- listenv(NULL)
  
  for (iter in seq(1:length(meanTerms)))  {
    rasTerm = meanTerms[iter]
    rasNames <- sprintf("%s", paste0(rasTerm, '_', year))

    white = "#FFFFFF"
    
    # Change base color for raster plots by budget term
    BaseColorval = ggplotColours(length(meanTerms))[iter]
    colorRamp <- colorRampPalette(c(white, BaseColorval))
    
    # swap rows for columns with transpose function: t()
    Array4Ras <-t(array(budgetByYear[, rasTerm], dim = c(M$ncols, M$nrows)))
    ras <- raster(Array4Ras)
    crs(ras) <- HARNSP17ft
    e <- extent(c(M$xmin, xmax, M$ymin, ymax))
    extent(ras) <- e
    
    breaks.qt <- tryCatch({
      classIntervals(Array4Ras[Array4Ras > 0],n = 10,style = "quantile",intervalClosure = "right")
    }, error = function(cond) {
      message(paste("\nerror with ClassIntervals of budget=", rasTerm))
      message("Script will not be plotting raster")
      return(NA)
    })
    if (printingOn) {
      if (!is.na(breaks.qt)) {
        pltMaps[[iter]] <-future({
          plotRaster(ras, colorRamp, breaks.qt$brks, rasNames)
        })
      }
    }
    }
  rm(datSrc)
  gc(verbose = TRUE)
}
