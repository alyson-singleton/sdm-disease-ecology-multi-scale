
#load data
env_data_national <- list.files(path="~/Desktop/brazil_schisto_raster_folders/national_reduced_late_feb", pattern="tif", all.files=FALSE, full.names=TRUE,recursive=TRUE)
e_national <- raster::stack(env_data_national)
env_data_sp <- list.files(path="~/Desktop/brazil_schisto_raster_folders/sp_late_feb", pattern="tif", all.files=FALSE, full.names=TRUE,recursive=TRUE)
e_sp <- raster::stack(env_data_sp)
env_data_mg <- list.files(path="~/Desktop/brazil_schisto_raster_folders/mg_late_feb", pattern="tif", all.files=FALSE, full.names=TRUE,recursive=TRUE)
e_mg <- raster::stack(env_data_mg)

#---------------------------------------------------------------------------------------------------------
# bio4_temp_ (*100), temperature seasonality (standard deviation of the monthly mean temperatures)
#---------------------------------------------------------------------------------------------------------
plot(e_national[[5]])
#plot prep
min <- 0; max <- 50 #bio11_temp
zlim <- range(c(min,max))
breakpoints <- round(c(seq(min, max, (max-min)/8)),0)
colors2 <- c(RColorBrewer::brewer.pal(9, "Oranges"))[2:9] #BuGn, 
par(mar = c(0, 0, 0, 0))
#par(mfrow=c(2,3), mar = c(1,3,2,0))#, mai=c(0.05,0.05,0.05,0.05))
#def.par <- par(no.readonly = TRUE)
layout(matrix(c(1,2,1,3), ncol=2, nrow=2, byrow=T), widths=c(2,1), 
       heights=c(2,2))

layout.show(n=3)

image(e_national[[5]],# main="A",  
      main="", #ylab = "B. straminea\n",
      axes = FALSE, box = FALSE, legend = F,
      breaks = breakpoints, col = colors2,
      #xlab="Longitude", ylab="Latitude",
      cex.main = 2, cex.lab=2, zlim=zlim, ann=F, asp=1)
#title(ylab = "B. straminea\n", adj=0.5, line=-0.59, font.lab=4, cex.lab=2)
title(main="A",adj=0.01, cex.main=2, line = -3)

# legend("right",fill=rev(colors2), title= "Temp (C)",
#        legend = rev(breakpoints[2:9]/100))

image(e_mg[[5]],# main="A",  
      main = "", #ylab = "B. straminea\n",
      axes = FALSE, box = FALSE, legend = F,
      breaks = breakpoints, col = colors2,
      #xlab = "Longitude", ylab = "Latitude",
      cex.main = 2, cex.lab = 2, zlim = zlim, ann=F, asp=1)
#title(ylab = "B. straminea\n", adj = 0.5, line = -0.59, font.lab = 4, cex.lab = 2)
#title("B", adj = 0.01, cex.main = 2, line = -3)

image(e_sp[[5]],# main="A",  
      main = "", #ylab = "B. straminea\n",
      axes = FALSE, box = FALSE, legend = F,
      breaks = breakpoints, col = colors2,
      #xlab = "Longitude", ylab = "Latitude",
      cex.main = 2, cex.lab = 2, zlim = zlim, ann=F, asp=1)
#title(ylab = "B. straminea\n", adj = 0.5, line = -0.59, font.lab = 4, cex.lab = 2)
#title("C", adj = 0.01, cex.main = 2, line = -3)

#---------------------------------------------------------------------------------------------------------
# bio11_temp (0.1 // 273.15), average daily mean of coldest quarter
#---------------------------------------------------------------------------------------------------------
#plot prep
plot(e_national[[2]])

min <- 2800; max <- 3050 #bio11_temp
zlim <- range(c(min,max))
breakpoints <- round(c(seq(min, max, (max-min)/8)),0)
colors2 <- c(RColorBrewer::brewer.pal(9, "Oranges"))[2:9] #BuGn, 
par(mar = c(0, 0, 0, 0))
#par(mfrow=c(2,3), mar = c(1,3,2,0))#, mai=c(0.05,0.05,0.05,0.05))
#def.par <- par(no.readonly = TRUE)
layout(matrix(c(1,2,1,3), ncol=2, nrow=2, byrow=T), widths=c(2,1), 
       heights=c(2,2))

layout.show(n=3)

image(e_national[[2]],# main="A",  
     main="", #ylab = "B. straminea\n",
     axes = FALSE, box = FALSE, legend = F,
     breaks = breakpoints, col = colors2,
     #xlab="Longitude", ylab="Latitude",
     cex.main = 2, cex.lab=2, zlim=zlim, ann=F, asp=1)
#title(ylab = "B. straminea\n", adj=0.5, line=-0.59, font.lab=4, cex.lab=2)
title(main="B",adj=0.01, cex.main=2, line = -3)
# legend("right",fill=rev(colors2), title= "Temp (C)",
#        legend = rev(c("10", "13", "16", "19", "22", "26", "29", "32")))

image(e_mg[[2]],# main="A",  
     main = "", #ylab = "B. straminea\n",
     axes = FALSE, box = FALSE, legend = F,
     breaks = breakpoints, col = colors2,
     #xlab = "Longitude", ylab = "Latitude",
     cex.main = 2, cex.lab = 2, zlim = zlim, ann=F, asp=1)
#title(ylab = "B. straminea\n", adj = 0.5, line = -0.59, font.lab = 4, cex.lab = 2)
#title("B", adj = 0.01, cex.main = 2, line = -3)

image(e_sp[[2]],# main="A",  
     main = "", #ylab = "B. straminea\n",
     axes = FALSE, box = FALSE, legend = F,
     breaks = breakpoints, col = colors2,
     #xlab = "Longitude", ylab = "Latitude",
     cex.main = 2, cex.lab = 2, zlim = zlim, ann=F, asp=1)
#title(ylab = "B. straminea\n", adj = 0.5, line = -0.59, font.lab = 4, cex.lab = 2)
#title("C", adj = 0.01, cex.main = 2, line = -3)

#---------------------------------------------------------------------------------------------------------
# bio16_prec (0.1), mean monthly precipitation of the wettest quarter
#---------------------------------------------------------------------------------------------------------
#plot prep
plot(e_national[[3]])

min <- 0; max <- 18000 #bio11_temp
zlim <- range(c(min,max))
breakpoints <- round(c(seq(min, max, (max-min)/8)),0)
colors2 <- c(RColorBrewer::brewer.pal(9, "Oranges"))[2:9] #BuGn, 
par(mar = c(0, 0, 0, 0))
#par(mfrow=c(2,3), mar = c(1,3,2,0))#, mai=c(0.05,0.05,0.05,0.05))
#def.par <- par(no.readonly = TRUE)
layout(matrix(c(1,2,1,3), ncol=2, nrow=2, byrow=T), widths=c(2,1), 
       heights=c(2,2))

layout.show(n=3)

image(e_national[[3]],# main="A",  
      main="", #ylab = "B. straminea\n",
      axes = FALSE, box = FALSE, legend = F,
      breaks = breakpoints, col = colors2,
      #xlab="Longitude", ylab="Latitude",
      cex.main = 2, cex.lab=2, zlim=zlim, ann=F, asp=1)
#title(ylab = "B. straminea\n", adj=0.5, line=-0.59, font.lab=4, cex.lab=2)
title(main="C",adj=0.01, cex.main=2, line = -3)
# legend("right",fill=rev(colors2), title= "Precip (kg/m^2)",
#        legend = rev(breakpoints[2:9])/100)

image(e_mg[[3]],# main="A",  
      main = "", #ylab = "B. straminea\n",
      axes = FALSE, box = FALSE, legend = F,
      breaks = breakpoints, col = colors2,
      #xlab = "Longitude", ylab = "Latitude",
      cex.main = 2, cex.lab = 2, zlim = zlim, ann=F, asp=1)
#title(ylab = "B. straminea\n", adj = 0.5, line = -0.59, font.lab = 4, cex.lab = 2)
#title("B", adj = 0.01, cex.main = 2, line = -3)

image(e_sp[[3]],# main="A",  
      main = "", #ylab = "B. straminea\n",
      axes = FALSE, box = FALSE, legend = F,
      breaks = breakpoints, col = colors2,
      #xlab = "Longitude", ylab = "Latitude",
      cex.main = 2, cex.lab = 2, zlim = zlim, ann=F, asp=1)
#title(ylab = "B. straminea\n", adj = 0.5, line = -0.59, font.lab = 4, cex.lab = 2)
#title("C", adj = 0.01, cex.main = 2, line = -3)

#---------------------------------------------------------------------------------------------------------
# bio17_prec (0.1), mean monthly precipitation amount of the driest quarter
#---------------------------------------------------------------------------------------------------------
#plot prep
plot(e_national[[4]])

min <- 0; max <- 10000 #bio11_temp
zlim <- range(c(min,max))
breakpoints <- round(c(seq(min, max, (max-min)/8)),0)
colors2 <- c(RColorBrewer::brewer.pal(9, "Oranges"))[2:9] #BuGn, 
par(mar = c(0, 0, 0, 0))
#par(mfrow=c(2,3), mar = c(1,3,2,0))#, mai=c(0.05,0.05,0.05,0.05))
#def.par <- par(no.readonly = TRUE)
layout(matrix(c(1,2,1,3), ncol=2, nrow=2, byrow=T), widths=c(2,1), 
       heights=c(2,2))

layout.show(n=3)

image(e_national[[4]],# main="A",  
      main="", #ylab = "B. straminea\n",
      axes = FALSE, box = FALSE, legend = F,
      breaks = breakpoints, col = colors2,
      #xlab="Longitude", ylab="Latitude",
      cex.main = 2, cex.lab=2, zlim=zlim, ann=F, asp=1)
#title(ylab = "B. straminea\n", adj=0.5, line=-0.59, font.lab=4, cex.lab=2)
title(main="D",adj=0.01, cex.main=2, line = -3)
# legend("right",fill=rev(colors2), title= "Precip (kg/m^2)",
#        legend = rev(breakpoints[2:9])/100)

image(e_mg[[4]],# main="A",  
      main = "", #ylab = "B. straminea\n",
      axes = FALSE, box = FALSE, legend = F,
      breaks = breakpoints, col = colors2,
      #xlab = "Longitude", ylab = "Latitude",
      cex.main = 2, cex.lab = 2, zlim = zlim, ann=F, asp=1)
#title(ylab = "B. straminea\n", adj = 0.5, line = -0.59, font.lab = 4, cex.lab = 2)
#title("B", adj = 0.01, cex.main = 2, line = -3)

image(e_sp[[4]],# main="A",  
      main = "", #ylab = "B. straminea\n",
      axes = FALSE, box = FALSE, legend = F,
      breaks = breakpoints, col = colors2,
      #xlab = "Longitude", ylab = "Latitude",
      cex.main = 2, cex.lab = 2, zlim = zlim, ann=F, asp=1)
#title(ylab = "B. straminea\n", adj = 0.5, line = -0.59, font.lab = 4, cex.lab = 2)
#title("C", adj = 0.01, cex.main = 2, line = -3)

#---------------------------------------------------------------------------------------------------------
# clay (*100), clay percentage
#---------------------------------------------------------------------------------------------------------
plot(e_national[[6]])
#plot prep
min <- 0; max <- 70 #bio11_temp
zlim <- range(c(min,max))
breakpoints <- round(c(seq(min, max, (max-min)/8)),0)
colors2 <- c(RColorBrewer::brewer.pal(9, "Blues"))[2:9] #BuGn, 
par(mar = c(0, 0, 0, 0))
#par(mfrow=c(2,3), mar = c(1,3,2,0))#, mai=c(0.05,0.05,0.05,0.05))
#def.par <- par(no.readonly = TRUE)
layout(matrix(c(1,2,1,3), ncol=2, nrow=2, byrow=T), widths=c(2,1), 
       heights=c(2,2))

layout.show(n=3)

image(e_national[[6]],# main="A",  
      main="", #ylab = "B. straminea\n",
      axes = FALSE, box = FALSE, legend = F,
      breaks = breakpoints, col = colors2,
      #xlab="Longitude", ylab="Latitude",
      cex.main = 2, cex.lab=2, zlim=zlim, ann=F, asp=1)
#title(ylab = "B. straminea\n", adj=0.5, line=-0.59, font.lab=4, cex.lab=2)
title(main="E",adj=0.01, cex.main=2, line = -3)

# legend("right",fill=rev(colors2), title= "Clay (%)",
#        legend = rev(breakpoints[2:9]/100))

image(e_mg[[6]],# main="A",  
      main = "", #ylab = "B. straminea\n",
      axes = FALSE, box = FALSE, legend = F,
      breaks = breakpoints, col = colors2,
      #xlab = "Longitude", ylab = "Latitude",
      cex.main = 2, cex.lab = 2, zlim = zlim, ann=F, asp=1)
#title(ylab = "B. straminea\n", adj = 0.5, line = -0.59, font.lab = 4, cex.lab = 2)
#title("B", adj = 0.01, cex.main = 2, line = -3)

image(e_sp[[6]],# main="A",  
      main = "", #ylab = "B. straminea\n",
      axes = FALSE, box = FALSE, legend = F,
      breaks = breakpoints, col = colors2,
      #xlab = "Longitude", ylab = "Latitude",
      cex.main = 2, cex.lab = 2, zlim = zlim, ann=F, asp=1)
#title(ylab = "B. straminea\n", adj = 0.5, line = -0.59, font.lab = 4, cex.lab = 2)
#title("C", adj = 0.01, cex.main = 2, line = -3)

#---------------------------------------------------------------------------------------------------------
# hnd, height above nearest drainage
#---------------------------------------------------------------------------------------------------------
plot(e_national[[8]])
#plot prep
min <- 0; max <- 720 #bio11_temp
zlim <- range(c(min,800))
breakpoints <- c(0, 25, 50, 100, 200, 400, 800)#round(c(seq(min, max, (max-min)/4)),0)
colors2 <- c(RColorBrewer::brewer.pal(9, "Blues"))[c(2,3,5,7,8,9)] #BuGn, 
par(mar = c(0, 0, 0, 0))
#par(mfrow=c(2,3), mar = c(1,3,2,0))#, mai=c(0.05,0.05,0.05,0.05))
#def.par <- par(no.readonly = TRUE)
layout(matrix(c(1,2,1,3), ncol=2, nrow=2, byrow=T), widths=c(2,1), 
       heights=c(2,2))

layout.show(n=3)

image(e_national[[8]],# main="A",  
      main="", #ylab = "B. straminea\n",
      axes = FALSE, box = FALSE, legend = F,
      breaks = breakpoints, col = colors2,
      #xlab="Longitude", ylab="Latitude",
      cex.main = 2, cex.lab=2, zlim=zlim, ann=F, asp=1)
#title(ylab = "B. straminea\n", adj=0.5, line=-0.59, font.lab=4, cex.lab=2)
title(main="F",adj=0.01, cex.main=2, line = -3)

# legend("right",fill=rev(colors2), title= "HND (m)",
#        legend = rev(breakpoints[2:7]))

image(e_mg[[8]],# main="A",  
      main = "", #ylab = "B. straminea\n",
      axes = FALSE, box = FALSE, legend = F,
      breaks = breakpoints, col = colors2,
      #xlab = "Longitude", ylab = "Latitude",
      cex.main = 2, cex.lab = 2, zlim = zlim, ann=F, asp=1)
#title(ylab = "B. straminea\n", adj = 0.5, line = -0.59, font.lab = 4, cex.lab = 2)
#title("B", adj = 0.01, cex.main = 2, line = -3)

image(e_sp[[8]],# main="A",  
      main = "", #ylab = "B. straminea\n",
      axes = FALSE, box = FALSE, legend = F,
      breaks = breakpoints, col = colors2,
      #xlab = "Longitude", ylab = "Latitude",
      cex.main = 2, cex.lab = 2, zlim = zlim, ann=F, asp=1)
#title(ylab = "B. straminea\n", adj = 0.5, line = -0.59, font.lab = 4, cex.lab = 2)
#title("C", adj = 0.01, cex.main = 2, line = -3)

#---------------------------------------------------------------------------------------------------------
# pH
#---------------------------------------------------------------------------------------------------------
plot(e_national[[9]])
#plot prep
min <- 40; max <- 75 #bio11_temp
zlim <- range(c(min,max))
breakpoints <- round(c(seq(min, max, (max-min)/7)),0)
colors2 <- c(RColorBrewer::brewer.pal(9, "Blues"))[2:8] #BuGn, 
par(mar = c(0, 0, 0, 0))
#par(mfrow=c(2,3), mar = c(1,3,2,0))#, mai=c(0.05,0.05,0.05,0.05))
#def.par <- par(no.readonly = TRUE)
layout(matrix(c(1,2,1,3), ncol=2, nrow=2, byrow=T), widths=c(2,1), 
       heights=c(2,2))

layout.show(n=3)

image(e_national[[9]],# main="A",  
      main="", #ylab = "B. straminea\n",
      axes = FALSE, box = FALSE, legend = F,
      breaks = breakpoints, col = colors2,
      #xlab="Longitude", ylab="Latitude",
      cex.main = 2, cex.lab=2, zlim=zlim, ann=F, asp=1)
#title(ylab = "B. straminea\n", adj=0.5, line=-0.59, font.lab=4, cex.lab=2)
title(main="G",adj=0.01, cex.main=2, line = -3)

# legend("right",fill=rev(colors2), title= "pH",
#        legend = rev(breakpoints[2:8]/10))

image(e_mg[[9]],# main="A",  
      main = "", #ylab = "B. straminea\n",
      axes = FALSE, box = FALSE, legend = F,
      breaks = breakpoints, col = colors2,
      #xlab = "Longitude", ylab = "Latitude",
      cex.main = 2, cex.lab = 2, zlim = zlim, ann=F, asp=1)
#title(ylab = "B. straminea\n", adj = 0.5, line = -0.59, font.lab = 4, cex.lab = 2)
#title("B", adj = 0.01, cex.main = 2, line = -3)

image(e_sp[[9]],# main="A",  
      main = "", #ylab = "B. straminea\n",
      axes = FALSE, box = FALSE, legend = F,
      breaks = breakpoints, col = colors2,
      #xlab = "Longitude", ylab = "Latitude",
      cex.main = 2, cex.lab = 2, zlim = zlim, ann=F, asp=1)
#title(ylab = "B. straminea\n", adj = 0.5, line = -0.59, font.lab = 4, cex.lab = 2)
#title("C", adj = 0.01, cex.main = 2, line = -3)

#---------------------------------------------------------------------------------------------------------
# soil water percentage
#---------------------------------------------------------------------------------------------------------
plot(e_national[[10]])
#plot prep
min <- 10; max <- 80 #bio11_temp
zlim <- range(c(min,max))
breakpoints <- round(c(seq(min, max, (max-min)/7)),0)
colors2 <- c(RColorBrewer::brewer.pal(9, "Blues"))[2:8] #BuGn, 
par(mar = c(0, 0, 0, 0))
#par(mfrow=c(2,3), mar = c(1,3,2,0))#, mai=c(0.05,0.05,0.05,0.05))
#def.par <- par(no.readonly = TRUE)
layout(matrix(c(1,2,1,3), ncol=2, nrow=2, byrow=T), widths=c(2,1), 
       heights=c(2,2))

layout.show(n=3)

image(e_national[[10]],# main="A",  
      main="", #ylab = "B. straminea\n",
      axes = FALSE, box = FALSE, legend = F,
      breaks = breakpoints, col = colors2,
      #xlab="Longitude", ylab="Latitude",
      cex.main = 2, cex.lab=2, zlim=zlim, ann=F, asp=1)
#title(ylab = "B. straminea\n", adj=0.5, line=-0.59, font.lab=4, cex.lab=2)
title(main="H",adj=0.01, cex.main=2, line = -3)

# legend("right",fill=rev(colors2), title= "Soil Water (%)",
#        legend = rev(breakpoints[2:8]/100))

image(e_mg[[10]],# main="A",  
      main = "", #ylab = "B. straminea\n",
      axes = FALSE, box = FALSE, legend = F,
      breaks = breakpoints, col = colors2,
      #xlab = "Longitude", ylab = "Latitude",
      cex.main = 2, cex.lab = 2, zlim = zlim, ann=F, asp=1)
#title(ylab = "B. straminea\n", adj = 0.5, line = -0.59, font.lab = 4, cex.lab = 2)
#title("B", adj = 0.01, cex.main = 2, line = -3)

image(e_sp[[10]],# main="A",  
      main = "", #ylab = "B. straminea\n",
      axes = FALSE, box = FALSE, legend = F,
      breaks = breakpoints, col = colors2,
      #xlab = "Longitude", ylab = "Latitude",
      cex.main = 2, cex.lab = 2, zlim = zlim, ann=F, asp=1)
#title(ylab = "B. straminea\n", adj = 0.5, line = -0.59, font.lab = 4, cex.lab = 2)
#title("C", adj = 0.01, cex.main = 2, line = -3)

#---------------------------------------------------------------------------------------------------------
# ag mosaic (%)
#---------------------------------------------------------------------------------------------------------
plot(e_national[[1]])
#plot prep
min <- 0; max <- 1 #bio11_temp
zlim <- range(c(min,max))
breakpoints <- c(seq(min, max, (max-min)/8))
colors2 <- c(RColorBrewer::brewer.pal(9, "BuGn"))[2:9] #BuGn, 
par(mar = c(0, 0, 0, 0))
#par(mfrow=c(2,3), mar = c(1,3,2,0))#, mai=c(0.05,0.05,0.05,0.05))
#def.par <- par(no.readonly = TRUE)
layout(matrix(c(1,2,1,3), ncol=2, nrow=2, byrow=T), widths=c(2,1), 
       heights=c(2,2))

layout.show(n=3)

image(e_national[[1]],# main="A",  
      main="", #ylab = "B. straminea\n",
      axes = FALSE, box = FALSE, legend = F,
      breaks = breakpoints, col = colors2,
      #xlab="Longitude", ylab="Latitude",
      cex.main = 2, cex.lab=2, zlim=zlim, ann=F, asp=1)
#title(ylab = "B. straminea\n", adj=0.5, line=-0.59, font.lab=4, cex.lab=2)
title(main="I",adj=0.01, cex.main=2, line = -3)

# legend("right",fill=rev(colors2), title= "Crop Cover (%)",
#        legend = rev(breakpoints[2:8]))

image(e_mg[[1]],# main="A",  
      main = "", #ylab = "B. straminea\n",
      axes = FALSE, box = FALSE, legend = F,
      breaks = breakpoints, col = colors2,
      #xlab = "Longitude", ylab = "Latitude",
      cex.main = 2, cex.lab = 2, zlim = zlim, ann=F, asp=1)
#title(ylab = "B. straminea\n", adj = 0.5, line = -0.59, font.lab = 4, cex.lab = 2)
#title("B", adj = 0.01, cex.main = 2, line = -3)

image(e_sp[[1]],# main="A",  
      main = "", #ylab = "B. straminea\n",
      axes = FALSE, box = FALSE, legend = F,
      breaks = breakpoints, col = colors2,
      #xlab = "Longitude", ylab = "Latitude",
      cex.main = 2, cex.lab = 2, zlim = zlim, ann=F, asp=1)
#title(ylab = "B. straminea\n", adj = 0.5, line = -0.59, font.lab = 4, cex.lab = 2)
#title("C", adj = 0.01, cex.main = 2, line = -3)

#---------------------------------------------------------------------------------------------------------
# ag mosaic (%)
#---------------------------------------------------------------------------------------------------------
plot(e_national[[7]])
#plot prep
min <- 0; max <- 20000 #bio11_temp
zlim <- range(c(min,max))
breakpoints <- c(seq(min, max, (max-min)/5))
colors2 <- rev(c(RColorBrewer::brewer.pal(9, "BuGn"))[c(2,3,4,5,9)]) #BuGn, 
par(mar = c(0, 0, 0, 0))
#par(mfrow=c(2,3), mar = c(1,3,2,0))#, mai=c(0.05,0.05,0.05,0.05))
#def.par <- par(no.readonly = TRUE)
layout(matrix(c(1,2,1,3), ncol=2, nrow=2, byrow=T), widths=c(2,1), 
       heights=c(2,2))

layout.show(n=3)

image(e_national[[7]],# main="A",  
      main="", #ylab = "B. straminea\n",
      axes = FALSE, box = FALSE, legend = F,
      breaks = breakpoints, col = colors2,
      #xlab="Longitude", ylab="Latitude",
      cex.main = 2, cex.lab=2, zlim=zlim, ann=F, asp=1)
#title(ylab = "B. straminea\n", adj=0.5, line=-0.59, font.lab=4, cex.lab=2)
title(main="J",adj=0.01, cex.main=2, line = -3)

# legend("right",fill=rev(colors2)[c(1,2,3,4,5)], title= "Distance (km)",
#        legend = rev(breakpoints[c(2,3,4,5,6)]/1000))

image(e_mg[[7]],# main="A",  
      main = "", #ylab = "B. straminea\n",
      axes = FALSE, box = FALSE, legend = F,
      breaks = breakpoints, col = colors2,
      #xlab = "Longitude", ylab = "Latitude",
      cex.main = 2, cex.lab = 2, zlim = zlim, ann=F, asp=1)
#title(ylab = "B. straminea\n", adj = 0.5, line = -0.59, font.lab = 4, cex.lab = 2)
#title("B", adj = 0.01, cex.main = 2, line = -3)

image(e_sp[[7]],# main="A",  
      main = "", #ylab = "B. straminea\n",
      axes = FALSE, box = FALSE, legend = F,
      breaks = breakpoints, col = colors2,
      #xlab = "Longitude", ylab = "Latitude",
      cex.main = 2, cex.lab = 2, zlim = zlim, ann=F, asp=1)
#title(ylab = "B. straminea\n", adj = 0.5, line = -0.59, font.lab = 4, cex.lab = 2)
#title("C", adj = 0.01, cex.main = 2, line = -3)
