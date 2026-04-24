## Generic functions
## tom.hengl@opengeohub.org

## most probablye translation
sumcor <- function(x, a, b, lev, max=4, stats=FALSE){
  o1 = summary(as.factor(x[x[,a]==lev,b]))
  o2 = data.frame(count=o1, levels = attr(o1, "names"))
  o2 = o2[order(o2$count, decreasing=TRUE),][1:max,]
  o2$input = lev
  if(stats==FALSE){
    o2 = data.frame(input=lev, Class_1=o2$levels[1], Class_2=o2$levels[2], Class_3=o2$levels[3])
  }
  o2
}

## remove LOD sign
rem.LOD = function(x){ as.numeric(ifelse(x=="< LOD", 0, as.numeric(x)))}

## Overlay and extract values using a tiling system:
extract.tiled <- function(obj, tile.pol, path="/data/tt/LandGIS/grid250m", ID="ID", cpus=parallel::detectCores(), snowfall=TRUE){
  obj$row.index <- 1:nrow(obj)
  ov.c <- over(spTransform(obj, CRS(proj4string(tile.pol))), tile.pol)
  message("Done overlaying points and tiles.")
  ov.t <- which(!is.na(ov.c[,ID]))
  ## for each point get the tile name:
  ov.c <- data.frame(ID=paste(ov.c[ov.t,ID]), row.index=ov.t, stringsAsFactors = FALSE)
  names(ov.c)[1] = ID
  tiles.lst <- basename(dirname(list.files(path=path, pattern=glob2rx("*.rds$"), recursive=TRUE)))
  ov.c <- ov.c[paste(ov.c[,ID]) %in% sapply(tiles.lst, function(i){strsplit(i, "T")[[1]][2]}),]
  tiles <- levels(as.factor(paste(ov.c[,ID])))
  cov.c <- as.list(tiles)
  names(cov.c) <- tiles
  ## extract using snowfall
  if(snowfall==TRUE){
    require(snowfall)
    sfInit(parallel=TRUE, cpus=cpus)
    sfExport(list=c("obj", "path", "ov.c", "ID", "cov.c", ".extract.tile", "tile.pol"))
    sfLibrary(raster)
    sfLibrary(rgdal)
    ov.lst <- sfClusterApplyLB(1:length(cov.c), function(i){try(.extract.tile(i, x=obj, ID=ID, path=path, ov.c=ov.c, cov.c=cov.c, tile.pol=tile.pol), silent = TRUE)})
    snowfall::sfStop()
  } else {
    ov.lst <- pbmcapply::pbmclapply(1:length(cov.c), function(i){ try(.extract.tile(i, x=obj, ID=ID, path=path, ov.c=ov.c, cov.c=cov.c, tile.pol=tile.pol), silent = TRUE) }, mc.cores = cpus)
  }
  ## bind together:
  message("Done running overlay in parallel.")
  out <- dplyr::bind_rows(ov.lst)
  out <- plyr::join(obj@data, as.data.frame(out), type="left", by="row.index")
  return(out)
}

.extract.tile <- function(i, x, ID, path, ov.c, cov.c, tile.pol, filename=NULL){
  row.index <- ov.c$row.index[ov.c[,ID]==names(cov.c)[i]]
  pnts <- x[row.index,]
  pnts <- spTransform(pnts, CRS(proj4string(tile.pol)))
  m <- readRDS(paste0(path, "/T", names(cov.c)[i], "/T", names(cov.c)[i], ".rds"))
  out <- sp::over(y=m, x=pnts)
  out$band1 <- NULL
  out$row.index <- row.index
  xy <- data.frame(pnts@coords)
  names(xy) <- c("X","Y")
  out <- cbind(out, xy)
  if(is.null(filename)){
    return(out)
  } else {
    saveRDS(out, filename)
    return(filename)
  }
  gc()
}


## complete cases
complete.vars <- function(x, sel=c("w15l2", "w3cld","ksat_lab", "ksat_field"), check.coords=TRUE, coords=c("longitude_decimal_degrees", "latitude_decimal_degrees"), remove.duplicates = FALSE){
  sel.rows <- lapply(sel, function(i){!is.na(x[,i])})
  sel.rows <- do.call(cbind, sel.rows)
  sel.rows <- rowSums(sel.rows) > 0
  if(check.coords==TRUE){
    geo.rows <- !is.na(x[,coords[1]]) & !is.na(x[,coords[2]]) & x[,coords[1]] > -180 & x[,coords[1]] < 180 & !(x[,coords[1]] == 0 & x[,coords[2]] == 0) & x[,coords[2]] > -90 & x[,coords[2]] < 90
    sel.rows <- sel.rows & geo.rows
  }
  x <- x[which(sel.rows),]
  if(remove.duplicates==TRUE){
    x <- x[!duplicated(x[,c(sel, coords)]),]
  }
  return(x)
}

plot_gh <- function(pnts, out.pdf, world, lats, longs, crs_goode = "+proj=igh"){
  # https://wilkelab.org/practicalgg/articles/goode.html
  require(cowplot)
  require(sf)
  require(rworldmap)
  require(ggplot2)
  if(missing(world)){ world <- sf::st_as_sf(rworldmap::getMap(resolution = "low")) }
  if(missing(lats)){
    lats <- c(
      90:-90, # right side down
      -90:0, 0:-90, # third cut bottom
      -90:0, 0:-90, # second cut bottom
      -90:0, 0:-90, # first cut bottom
      -90:90, # left side up
      90:0, 0:90, # cut top
      90 # close
    )
  }
  if(missing(longs)){
    longs <- c(
      rep(180, 181), # right side down
      rep(c(80.01, 79.99), each = 91), # third cut bottom
      rep(c(-19.99, -20.01), each = 91), # second cut bottom
      rep(c(-99.99, -100.01), each = 91), # first cut bottom
      rep(-180, 181), # left side up
      rep(c(-40.01, -39.99), each = 91), # cut top
      180 # close
     )
  }
  goode_outline <-
    list(cbind(longs, lats)) %>%
    st_polygon() %>%
    st_sfc(
      crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    )
  # now we need to work in transformed coordinates, not in long-lat coordinates
  goode_outline <- st_transform(goode_outline, crs = crs_goode)
  # get the bounding box in transformed coordinates and expand by 10%
  xlim <- st_bbox(goode_outline)[c("xmin", "xmax")]*1.1
  ylim <- st_bbox(goode_outline)[c("ymin", "ymax")]*1.1
  # turn into enclosing rectangle
  goode_encl_rect <-
    list(
      cbind(
        c(xlim[1], xlim[2], xlim[2], xlim[1], xlim[1]),
        c(ylim[1], ylim[1], ylim[2], ylim[2], ylim[1])
      )
    ) %>%
    st_polygon() %>%
    st_sfc(crs = crs_goode)
  # calculate the area outside the earth outline as the difference
  # between the enclosing rectangle and the earth outline
  goode_without <- st_difference(goode_encl_rect, goode_outline)
  m <- ggplot(world) + geom_sf(fill = "gray80", color = "black", size = 0.5/.pt) +
    geom_sf(data = goode_without, fill = "white", color = "NA") +
    geom_sf(data = goode_outline, fill = NA, color = "gray30", size = 0.5/.pt) +
    cowplot::theme_minimal_grid() + theme(panel.background = element_rect(fill = "#56B4E950", color = "white", size = 1),  panel.grid.major = element_line(color = "gray30", size = 0.25)) +
    geom_sf(data = pnts, size = 0.8, shape = 21, fill = "yellow", color="black") +
    #geom_sf(data = pnts, size = 1, pch="+", color="black") +
    coord_sf(crs = crs_goode, xlim = 0.95*xlim, ylim = 0.95*ylim, expand = FALSE)
  ggsave(out.pdf, m, dpi=150, height = 5.35, width = 9)
}

cols2dms <- function(x,y,z,e){ifelse(is.na(e)|is.na(x), NA, as(sp::char2dms(paste(x, "d", y, "'", z, "\"", e, sep="")), "numeric"))}

saveRDS.gz <- function(object,file,threads=parallel::detectCores()) {
  con <- pipe(paste0("pigz -p",threads," > ",file),"wb")
  saveRDS(object, file = con)
  close(con)
}

readRDS.gz <- function(file,threads=parallel::detectCores()) {
  con <- pipe(paste0("pigz -d -c -p",threads," ",file))
  object <- readRDS(file = con)
  close(con)
  return(object)
}

## translate names
transvalues = function(df, out.name, in.name, fun.lst){
  if(!length(out.name)==length(in.name)){
    stop("Arguments 'out.name' and 'in.name' not equal length")
  }
  if(missing(fun.lst)){
    fun.lst = as.list(rep("x*1", length(out.name)))
  }
  ## https://stackoverflow.com/questions/61094854/storing-functions-in-an-r-list
  utility.fns = lapply(1:length(fun.lst), function(i){function(x){eval(parse(text = fun.lst[[i]]) )}})
  out <- as.data.frame(lapply(1:length(in.name), function(i){utility.fns[[i]](df[,in.name[i]])}))
  names(out) = out.name
  return(out)
}

add_s = function(x){
  y = c("*and$", "*bel$", "*tel$", "*hel$", "*uod$", "*yod$", "*hod$", "*mod$", "*ent$", "*ert$", "*ept$", "*ult$", "*cid$", "*bid$", "*lid$", "*gid$", "*hid$", "*rid$", "*sid$", "*alf$", "*oll$", "*ist$", "*sol$", "*ak$", "*zem$")
  s = sapply(y, function(i){grep(i, x)})
  if(sum(unlist(s))>0){
    x = paste0(x, "s")
  }
  x
}
