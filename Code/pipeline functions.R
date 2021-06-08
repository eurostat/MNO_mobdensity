# pop_gen_param <- function (tile.num = 10000, # number of tiles
#           base.tile.size = 100, # size of a single tile
#           city.num = 4, # number of single cities
#           city.size = 2000, # size of a single city polygon
#           hole.num = 20, # number of single holes
#           hole.size = 400) 
# {
#   nms <- names(formals(pop_gen_param))
#   lst <- sapply(nms, get, envir = environment(), simplify = FALSE)
#   class(lst) <- "pop_gen_param"
#   lst
# }


# Population generation function
pop_gen <- function(tile.num, base.tile.size, 
                    city.num, city.size, city.shape = "SQUARE", 
                    hole.num, hole.size, hole.shape = "SQUARE",
                    pop.dist.ls) {

  
  # save area parameters
  area.params <- list(tile.num = tile.num, 
                      base.tile.size = base.tile.size, 
                      city.num = city.num, 
                      city.size = city.size, 
                      city.shape = city.shape, 
                      hole.num = hole.num, 
                      hole.size = hole.size, 
                      hole.shape = hole.shape,
                      pop.dist.df = pop.dist.df)
  
  # calculate complete area size
  poly.size <- sqrt(tile.num) * base.tile.size
  
  # build basis polygon
  area.polygon <- st_polygon(list(rbind(c(0, 0), c(poly.size, 0), c(poly.size, poly.size), c(0, poly.size), c(0, 0)))) %>%
    st_sfc() %>%
    st_sf()
  
  # create bounding box area and corresponding raster object
  area.bbox <- st_bbox(area.polygon)
  area.raster <- create_raster(area.bbox, tile.size = base.tile.size)
  
  # specify elevation raster --> default at 0 currently
  area.elevation <- area.raster
  values(area.elevation) <- 0
  
  # retransform raster to sf for cell creation
  base.tiles <- st_as_sf(st_as_stars(area.raster))
  
  # build city geometry (1 obs.)
  cities <- st_sample(area.polygon, city.num) %>%
    st_buffer(dist = city.size, endCapStyle = city.shape) %>% # parameter if square or circle
    st_geometry() %>%
    st_union() %>%
    st_sf() %>%
    mutate(city = 1)
  
  # build hole geometry (1 obs.)
  holes <- st_sample(cities, hole.num) %>%
    st_buffer(dist = hole.size, endCapStyle = as.character(hole.shape)) %>% # parameter if square or circle
    st_geometry() %>%
    st_union() %>%
    st_sf() %>%
    mutate(hole = 1)
  
  # join tiles with cities and holes
  area.sf.helper <- base.tiles %>%
    st_join(cities) %>%
    st_join(holes) %>%
    mutate(type = case_when(city == 1 & is.na(hole) ~ "Urban",
                            city == 1 & hole == 1 ~ "Hole",
                            TRUE ~ "Rural")) %>%
    mutate(category = case_when(type %in% c("Urban", "Hole") ~ "Urban",
                                type == "Rural" ~ "Rural")) %>% 
    rownames_to_column("tile.id") %>%
    mutate(tile.id = as.integer(tile.id)) %>% 
    arrange(type)
  
  # find out how many tiles are of a certain type to match with respective pop-distribution function
  summary.area <- area.sf.helper %>% 
    st_drop_geometry() %>% 
    dplyr::select(tile.id, type) %>% 
    group_by(type) %>% 
    summarise(n = n()) %>% 
    left_join(pop.dist.df, by = "type") %>% 
    mutate(final = paste0(exprssion, ", n = ", n, ")")) %>% 
    arrange(type)
    
  # sample pop vector from input distribution function
  pop.helper <- unlist(map(summary.area$final, ~eval(parse(text = .x))))
  
  # append pop to the sf data frame, apply necessary rounding and create final sf dataframe
  area.sf <- area.sf.helper %>% 
    mutate(pop = pop.helper) %>% 
    mutate(pop = round(pop, 0)) %>%
    mutate(pop = if_else(pop < 0, 0, pop)) %>%
    mutate(centroid.geometry = st_centroid(.$geometry)) %>% 
    mutate(X.centroid = unlist(map(.$centroid.geometry, 1)),
           Y.centroid = unlist(map(.$centroid.geometry, 2))) %>% 
    dplyr::select(tile.id, type, category, pop, X.centroid, Y.centroid) %>%
    arrange(tile.id)
  
  # create non-sf data frame version
  area.df <- area.sf %>% 
    st_drop_geometry() %>% 
    arrange(tile.id) 
  
  # put everything into a output list
  final <- list(area.params = area.params,
                area.sf = area.sf,
                area.df = area.df,
                area.union = area.polygon,
                area.bbox = area.bbox,
                area.raster = area.raster,
                area.elevation = area.elevation)
  
  return(final)
}



# create tower positions with attached cells
create_cells <- function(area.sf, tower.dist, rotation.deg, jitter, small = FALSE, subscript, seed) {
  
  
  set.seed = seed   
  
  rotation = function(a){
    r = a * pi / 180 #degrees to radians
    matrix(c(cos(r), sin(r), -sin(r), cos(r)), nrow = 2, ncol = 2)
  } 
  
  layer_network_generate = function(x, tower.dist, rotation.deg){
    layer.geo <- x %>% 
      st_make_grid(cellsize = tower.dist, 
                   square = F, # hexagon 
                   flat_topped = T) %>%  # different cell size (qm)
      st_geometry()
    
    layer.centroid <- st_centroid(layer.geo)
    layer <- (layer.geo - layer.centroid) * rotation(rotation.deg) + layer.centroid # rotate by 35 degrees
    return(layer)
    
  }
  
  # create layer object, placing towers
  layer <- layer_network_generate(x = area.sf, tower.dist = tower.dist, rotation.deg = rotation.deg)
  
  # specify exact location of towers and labelling
  towers <- layer %>%
    st_centroid() %>% 
    st_jitter(jitter) %>%
    st_coordinates() %>%
    as_tibble() %>%
    dplyr::select(X.tow = X, Y.tow = Y) %>% 
    mutate(tower.id = paste0(subscript, ".", 1:n()))
  
  # create 3 cells per tower and labelling
  cells.unparam <- towers %>% 
    slice(rep(1:n(), each = 3)) %>%
    group_by(tower.id) %>%
    mutate(cell = paste(tower.id, "C", 1:3, sep = ".")) %>%
    ungroup() %>%
    mutate(cell.kind = subscript) %>% 
    mutate(intra.cell.number = str_sub(cell, -1)) %>% 
    mutate(small = small) %>% 
    mutate(rotation.deg = rotation.deg)
  
  return(cells.unparam)
  
}


# area.elevation <- sim.area$area.elevation
# area.sf <- sim.area$area.sf
# area.bbox <- sim.area$area.bbox
# cells.unparam <- MA.cells.unparam

# specify parameters of each cell
create_cellplan <- function(area.sf, area.bbox, area.elevation, cells.unparam, cell.param.mobloc) {
  
  crs.set <- st_crs(area.sf)
  
  # specify parameters if missing (NAs will automatically be replaced with mobloc default value)
  if (missing(cell.param.mobloc)) {
    specified.param <- mobloc_param()
  } else {
    specified.param <- cell.param.mobloc
  }
  
  # prepare cells for validation
  cellplan.unval <- cells.unparam %>% 
    mutate(direction.int = case_when(small == FALSE & intra.cell.number == "1" ~ rotation.deg + 0,
                                     small == FALSE & intra.cell.number == "2" ~ rotation.deg + 120,
                                     small == FALSE & intra.cell.number == "3" ~ rotation.deg + 240,
                                     TRUE ~ NA_real_)) %>% 
    mutate(direction = case_when(is.na(direction.int) ~ NA_real_,
                                 direction.int > 360 ~ 360 - direction.int, # direction needs to be [0; 360]
                                 direction.int <= 360 ~ direction.int),
           # tilt, etc will be changed if parameters are set in cell.param
           tilt = NA,
           beam_h = NA,
           beam_v = NA) %>%
    st_as_sf(coords = c("X.tow", "Y.tow"), crs = crs.set) %>%
    st_crop(area.bbox) %>%
    st_set_agr("aggregate") %>% 
    st_intersection(area.sf) %>% 
    dplyr::select(cell, small, direction, tilt, beam_h, beam_v)
  
  
  # validate cellplan and add antenna position indication through offset
  cellplan.val <- validate_cellplan(cp = cellplan.unval, 
                                    param = specified.param,
                                    elevation = area.elevation
                                    # region = region,
                                    # envir = envir
                                    ) 
    # mutate(move_cells_into_prop_direction(., offset = 100)) %>% 
    # mutate(x.offset = st_coordinates(.)[, 1],
    #        y.offset = st_coordinates(.)[, 2]) %>% 
    # st_drop_geometry()
  
  # put everything into a list
  final <- list(cellplan.val = cellplan.val,
                cell.param.mobloc = specified.param)
  
  return(final)
  
}


smart_round <- function(x, digits = 0) {
  up <- 10 ^ digits
  x <- x * up
  y <- floor(x)
  indices <- tail(order(x - y), round(sum(x)) - sum(y))
  y[indices] <- y[indices] + 1
  y / up
}


# implement selected iterations option with input vector

# New MLE iteration function by Matyas
EM_est <- function(c.vec.dt, P.dt, a.vec.dt, n.iter, selected.range, ldt = 10^-04, message = T) {
  
  cdt <- c.vec.dt
  pij <- cdt[P.dt, on = "i"]
  pij <- pij[c > 0] # remove those lines where c==0 because it will create 0 division
  # adt <- data.table(a = a.vec)
  # tiles <- adt[, .(j = 1:.N, u = a)]
  tiles <- a.vec.dt
  keep <- a.vec.dt # base dataframe for the selected iterations
  
  for(m in 1:(n.iter)){
    
    if(message == T) {
      cat(format(Sys.time()), paste0("---- calculating u", m), "----\n")
    }
    
    cols <- c("j", paste0("u"))
    ju <- tiles[, cols, with = F]
    setnames(ju, c("j", "u"))
    pij <- ju[pij, on = "j"]
    denom <- pij[, .(sum_pik_uk = sum(u * pij)), by = i]
    pij <- denom[pij, on = "i"]
    faktor <- pij[, .(f = sum(c * pij / sum_pik_uk)), by = j]
    faktor.adj <- faktor[, f := fifelse(test = {is.na(f) | is.nan(f) | is.infinite(f)}, 1, f)] # if else to assure that the posterior is 1 to secure the same estimand value after ldt
    pij[, c("u", "sum_pik_uk") := NULL]
    tiles <- faktor.adj[tiles, on = "j"]
    # tiles <- eval(parse(text = paste0("tiles[, u := u * f]")))
    tiles <- eval(parse(text = paste0("tiles[,  u := fifelse(u * f < ldt, 0, u * f)]")))
    tiles[, "f" := NULL]
    
    if(m %in% selected.range) {
      keep <- tiles[keep, on = "j"]
      keep <- eval(parse(text = paste0("keep[, u", m, ":= u]")))
      keep[, "u" := NULL]
    }
    
  }
  
  # final <- list(tiles = tiles,
  #               keep = keep)
  
  return(keep)
  
}



# # New MLE iteration function by Matyas
# EM_est <- function(c.vec, P.dt, a.vec, n.iter, ldt = 10^-04) {
# 
# 
#   cdt <- data.table(c = c.vec)
#   cdt <- cdt[, .(i = 1:.N, c = c)]
#   pag <- cdt[P.dt, on = "i"]
#   pag <- pag[c > 0] # remove those lines where c==0 because it will create 0 division
#   adt <- data.table(a = a.vec)
#   tiles <- adt[, .(j = 1:.N, u0 = a)]
# 
#   for(m in 0:(n.iter - 1)){
#     cat(format(Sys.time()), paste0("---- calculating u", m + 1), "----\n")
#     cols <- c("j", paste0("u", m))
#     ju <- tiles[, cols, with = F]
#     setnames(ju, c("j", "u"))
#     pag <- ju[pag, on = "j"]
#     denom <- pag[, .(sum_pik_uk = sum(u * pag)), by = i]
#     pag <- denom[pag, on = "i"]
#     faktor <- pag[, .(f = sum(c * pag / sum_pik_uk)), by = j]
#     faktor.adj <- faktor[, f := fifelse(test = {is.na(f) | is.nan(f) | is.infinite(f)}, 1, f)] # if else to assure that the posterior is 1 to secure the same estimand value after ldt
#     pag[, c("u", "sum_pik_uk") := NULL]
#     tiles <- faktor.adj[tiles, on = "j"]
#     tiles <- eval(parse(text = paste0("tiles[, u", m + 1, " := u", m, "* f]")))
#     tiles <- eval(parse(text = paste0("tiles[,  u", m + 1, " := fifelse(u", m + 1, " < ldt, 0, u", m + 1, ")]")))
#     tiles[, "f" := NULL]
#   }
# 
#   return(tiles)
# 
# }



DF_est <- function(c.vec.dt, P.star.spm, a.supertile.vec){
  
  c.vec <- c(c.vec.dt)$c
  A.spm <- .sparseDiagonal(n = length(a.supertile.vec), x = a.supertile.vec)
  
  Y <- P.star.spm %*% A.spm %*% t(P.star.spm) 
  Y1 <- VCA::MPinv(Y) %*% (c.vec - P.star.spm %*% a.supertile.vec)
  u <- as.vector(A.spm %*% t(P.star.spm) %*% Y1 + a.supertile.vec)
  
  return(u)
}


DF_est_iterated <- function(c.vec.dt, P.star.spm, a.supertile.vec, P.dt, DF.threshold = 1, selected.range.DF, n.iter.DF){
  
  c.vec <- c(c.vec.dt)$c
  
  DF.tiles <- data.table(j = as.numeric(names(a.supertile.vec)),
                         u = a.supertile.vec)
  keep <- data.table(j = as.numeric(names(a.supertile.vec)),
                     prior = a.supertile.vec)
  
  for(m.DF in 1:(n.iter.DF)){
    
    m.DF = 1

    
    # time/iteration indication
    cat(format(Sys.time()), paste0("---- calculating u", m.DF), "----\n")
    
    # specify prior (A)
    U.spm <- .sparseDiagonal(n = length(a.supertile.vec), x = DF.tiles$u)
    
    # calculate DF
    Y <- P.star.spm %*% U.spm %*% t(P.star.spm)
    Y1 <- VCA::MPinv(Y) %*% (c.vec - P.star.spm %*% a.supertile.vec)
    u.unadj <- as.vector(U.spm %*% t(P.star.spm) %*% Y1 + a.supertile.vec)
    
    # save as datatable
    DF.tiles <- data.table(j = as.numeric(names(a.supertile.vec)),
                           u = u.unadj)
    
    # Transforming based on dynamic DF.threshold
    # DF.tiles <- eval(parse(text = paste0("DF.tiles[, u := fifelse(u < (DF.threshold / m.DF), DF.threshold / m.DF, u)]")))
    DF.tiles[, u := fifelse(u < (DF.threshold / m.DF), DF.threshold / m.DF, u)]
    
    # Renormalizing with 1 EM iteration
    u.dt <- EM_est(c.vec.dt = c.vec.dt,
                   P.dt = P.star.oracle.supertile.dt,
                   a.vec.dt = DF.tiles,
                   selected.range = c(1,2, 3),
                   n.iter = 3)
    setnames(u.dt, c("j", "u.prior", "u"))
    u.dt.final <- u.dt[, c("j", "u")]
    DF.tiles <- DF.tiles[, "j"]
    DF.tiles <- u.dt.final[DF.tiles, on = "j"]
    
    # define datatable object "keep" which contains only the iterations that were indicated in DF.selected.range
    if(m.DF %in% selected.range.DF) {
      keep <- DF.tiles[keep, on = "j"]
      keep <- eval(parse(text = paste0("keep[, u", m.DF, ":= u]")))
      keep[, "u" := NULL]
    }
    
  }
  
  return(keep)
}



# paramter: sim.area, cellplan.combined, signal.strength.comb.dt, c.vec.df, method, offset


# Voronoi
# aggregating the antennas to towers (and corresponding values) and identifying problematic tower locations (outside of focus area)

### 1 ####


VOR_est <- function(area, cellplan.combined, signal.strength.comb.dt, C.vec.df, seed = c("tower", "cell.offset", "cell.hotpoint"), offset = 20) {
  
  ### helper 
  
  base.tile.size <- area$area.params$base.tile.size
  crs.set <- st_crs(area$area.sf)
  cellplan.reduced <- cellplan.combined %>% 
    dplyr::select(cell, direction)
  
  if (seed == "tower") {
    
    ### tower ###
    seed.object.unadj <- C.vec.df %>% 
      left_join(cellplan.combined, by = "cell") %>% 
      mutate(tower = str_extract(cell, "[A-Z]+.[:digit:]+")) %>% 
      group_by(tower) %>% 
      mutate(phones.sum = sum(phones.sum)) %>% 
      distinct(tower, .keep_all = T) %>% 
      ungroup() %>% 
      dplyr::select(seed = tower, X.tow = x, Y.tow = y, phones.sum) %>%
      st_as_sf(coords = c("X.tow", "Y.tow")) %>% 
      st_sf(crs = crs.set) %>%  # optional
      mutate(within.fa = lengths(st_within(., area$area.union))) # find seeds outside the focus area
    
  } else if (seed == "cell.offset") {
    
    #### offset ####
    cellplan.offset <- move_cells_into_prop_direction(st_as_sf(cellplan.combined), offset = offset)
    
    seed.object.unadj <- C.vec.df %>% 
      left_join(cellplan.offset, by = "cell") %>% 
      dplyr::select(seed = cell, phones.sum, geometry) %>%
      st_sf() %>% 
      st_sf(crs = crs.set) %>%
      mutate(within.fa = lengths(st_within(., area$area.union))) # find seeds outside the focus area
    
  } else if (seed == "cell.hotpoint") {
    
    ### Hotpoint ###
    cell.max.location <- signal.strength.comb.dt %>% 
      as_tibble() %>% 
      rename(tile.id = rid) %>% 
      group_by(cell) %>% 
      mutate(max.cell.s = max(dBm)) %>% # find the max SIGNAL STRENGTH per cell
      filter(dBm == max.cell.s) %>%  # filter the rows where corresponds to the cell specific max.s to retain the tile.id
      ungroup() %>% 
      group_by(tile.id) %>% 
      mutate(count.same.centroids = n()) %>% # tiles that act as hotpoint for multiple cells
      ungroup() %>% 
      mutate(same.centroids = case_when(count.same.centroids > 1 ~ 1, # mark identical hotpoints for further adjustment into respective direction
                                        TRUE ~ 0)) %>% 
      left_join(area$area.sf, by = "tile.id") %>% 
      left_join(cellplan.reduced, by = "cell") %>% # join to receive direction info
      dplyr::select(cell, same.centroids, direction, geometry) %>% 
      st_sf(crs = crs.set) %>% 
      st_centroid() %>% 
      mutate(X = unlist(map(.$geometry, 1)),
             Y = unlist(map(.$geometry, 2))) %>% 
      mutate(X.adj = case_when(same.centroids == "1" ~ X + SIN(direction) * 10, # move point in cell defined direction with offset of 10m if identical hotpoints
                               TRUE ~ X),
             Y.adj = case_when(same.centroids == "1" ~ Y + COS(direction) * 10,
                               TRUE ~ Y)) %>% 
      st_drop_geometry() %>% 
      st_as_sf(coords = c("X.adj", "Y.adj"), crs = crs.set) %>% 
      dplyr::select(-X, -Y) # add again for consistency check plot if the adjutsment was right
    
    # cell.max.location %>% 
    #   ggplot() +
    #   geom_sf() +
    #   geom_point(aes(x = X, y = Y, color = factor(same.centroids)))
    
    
    seed.object.unadj <- C.vec.df %>% 
      left_join(cell.max.location, by = "cell") %>% 
      dplyr::select(seed = cell, phones.sum, geometry) %>%
      mutate(seed = as.character(seed)) %>% 
      st_sf(crs = crs.set) %>% 
      mutate(within.fa = lengths(st_within(., area$area.union))) # find seeds outside the focus area
    
  } else if (seed == "cell.barycenter") {
    
    ### barycenter ###
    
    # calculate the centroid of all tiles and save as point coordinates
    tile.centroid <- area$area.sf %>% 
      st_centroid() %>% 
      mutate(centroid.X = unlist(map(.$geometry, 1)),
             centroid.Y = unlist(map(.$geometry, 2))) %>% 
      dplyr::select(tile.id, centroid.X, centroid.Y) %>% 
      st_drop_geometry()
    
    # join signal strength with tile centroid. calculate weighted 2d mean of x and y coordinate with mean(sig.dom.) as weights
    cell.max.location <- signal.strength.comb.dt %>% 
      as_tibble() %>% 
      rename(tile.id = rid) %>% 
      left_join(tile.centroid, by = "tile.id") %>% 
      group_by(cell) %>% 
      summarise(X.barycenter = weighted.mean(centroid.X, s),
                Y.barycenter = weighted.mean(centroid.Y, s)) %>% # calculate for each cell the barycenter location
      ungroup() %>% 
      st_as_sf(coords = c("X.barycenter", "Y.barycenter"), crs = crs.set)
    
    
    seed.object.unadj <- C.vec.df %>% 
      left_join(cell.max.location, by = "cell") %>% 
      dplyr::select(seed = cell, phones.sum, geometry) %>%
      mutate(seed = as.character(seed)) %>% 
      st_sf(crs = crs.set) %>% 
      mutate(within.fa = lengths(st_within(., area$area.union))) # find seeds outside the focus area
    
  }
  
  
  # saving seed IDs of problematic seeds (2)
  
  ### 2 ####
  
  seed.outside <- seed.object.unadj %>% 
    filter(within.fa == 0) %>%
    st_drop_geometry() %>% 
    dplyr::select(seed) %>% 
    deframe()
  
  # filtering unproblematic seeds
  seed.inside <- seed.object.unadj %>% 
    filter(within.fa == 1)
  
  if (length(seed.outside) >= 1) {
    
    seed.object.adj <- seed.object.unadj %>% 
      filter(within.fa == 0) %>%
      st_nearest_points(., area$area.union) %>%
      st_cast("POINT") %>%
      .[seq(2, length(.), 2)] %>% #
      st_as_sf() %>%
      mutate(seed = names(seed.outside)) %>%
      mutate(geometry = x) %>%
      st_sf(sf_column_name = "geometry") %>%
      dplyr::select(-x) %>%
      bind_rows(seed.inside)
    
  } else {
    
    seed.object.adj <- seed.object.unadj
  }
  
  
  # Finding respective nearest point on focus area border for every problematic seed location
  
  
  sum(unlist(st_intersects(seed.object.adj, area$area.union))) ==
    length(seed.object.unadj$seed) # check if all are within now
  
  # Using the seed object to calculate the Voronoi regions and their spatial densities per Voronoi region
  
  ### 3 ####
  seed.voronoi.est <- seed.object.adj %>%  
    st_geometry() %>% 
    st_union() %>% 
    st_voronoi() %>% 
    st_collection_extract(type = "POLYGON") %>%
    st_sf() %>% # check if crs is listed
    st_join(seed.object.unadj) %>%  # rejoin with seed object to retain seed id
    st_intersection(area$area.union) %>% 
    mutate(vor.area = row_number()) %>% 
    mutate(vor.area.size = as.numeric(st_area(.$geometry))) %>% 
    mutate(vor.est = phones.sum / vor.area.size)
  
  Voronoi.regions.plot <- seed.voronoi.est %>%
    ggplot() +
    # geom_sf(aes(fill = phones.sum), color = "blue") +
    geom_sf(color = "blue") +
    theme(text = element_text(size = 13))
    # scale_fill_viridis_c("Phones") +
  
  # Joining the regions with the tile specific data
  seed.voronoi.tile <- seed.voronoi.est %>% 
    st_join(area$area.sf) %>% # & re-connect the data items
    st_set_agr("aggregate") %>% # clean up
    group_by(tile.id) %>% 
    mutate(count = n()) %>% 
    ungroup()
  
  # identifiying tiles intersecting with multiple Voronoi regions
  seed.multiple <- seed.voronoi.tile %>%
    st_drop_geometry() %>% 
    filter(count > 1) %>% 
    distinct(tile.id) %>% 
    deframe()
  
  # calculate area within competing voronoi regions of "multiple" tiles
  seed.intersect.tiles <- area$area.sf %>% 
    filter(tile.id %in% seed.multiple) %>%
    st_intersection(seed.voronoi.est) %>% 
    # st_collection_extract(type = "POLYGON") %>% # select the polygons
    mutate(amount.tiles = as.numeric(st_area(.$geometry)) / base.tile.size^2) # checked if it adds up to 1
  
  # final datatset to calculate spatial density
  seed.voronoi.final <- seed.intersect.tiles %>% 
    st_drop_geometry() %>% 
    dplyr::select(tile.id, seed, amount.tiles) %>% 
    right_join(seed.voronoi.tile, by = c("tile.id", "seed")) %>% 
    mutate(amount.tiles = case_when(is.na(amount.tiles) ~ 1,
                                    TRUE ~ amount.tiles)) %>% 
    group_by(tile.id) %>% 
    summarise(u.VOR = weighted.mean(x = vor.est, w = amount.tiles) * base.tile.size^2)
  # should result in same length as the raw tiles object and the sum of the voronoi est corrected should resemble the sum of the c.vec
  
  
  return(list(seed.voronoi.final = seed.voronoi.final,
              Voronoi.regions.plot = Voronoi.regions.plot))
  
}






###################
### Plots #########
###################

map_density <- function(data, var, label, pointsize = 1.9, pixels = c(900, 900)) {
  
  colors <- c("white", "light grey", "light blue", "blue", "light green", "yellow", "orange", "red", "#654321")
  var.label <- paste(label)
  
  plot <- data %>% 
    ggplot(aes(x = X.centroid, y = Y.centroid)) +
    # geom_sf(aes_string(fill = var), color = "transparent") +
    geom_scattermore(aes_string(color = var), pointsize = pointsize, pixels = pixels) +
    # scico::scale_fill_scico(palette = "bilbao", limits = c(0, 3.48), direction = 1) +
    scale_color_manual(values = colors, drop = F, name = label) +
    coord_sf() +
    theme_minimal() +
    theme(text = element_text(size = 13)) +
    labs(x = "", y = "") +
    guides(colour = guide_legend(override.aes = list(shape = 15, size = 5)))
    # theme(axis.title.x = element_blank(),
    #       axis.text.x = element_blank(),
    #       axis.ticks.x = element_blank(),
    #       axis.title.y = element_blank(),
    #       axis.text.y = element_blank(),
    #       axis.ticks.y = element_blank()) +
  
  
  return(plot)
}




r2d <- function(x) x * 180 / pi
d2r <- function(x) x / 180 * pi

COS <- function(x) cos(d2r(x))
SIN <- function(x) sin(d2r(x))
TAN <- function(x) tan(d2r(x))

ACOS <- function(x) r2d(acos(x))
ASIN <- function(x) r2d(asin(x))
ATAN <- function(x) r2d(atan(x))
ATAN2 <- function(y, x) r2d(atan2(y, x))



db2s <- function(dBm, midpoint, steepness) {
  scale <- (dBm - midpoint) * steepness
  1 / (1 + exp(1)^(-scale))
}


dBW2dBm <- function(dBW) {
  dBW + 30
}


dBW2W <- function(dBW) {
  10^(dBW / 10)
}


W2dBW <- function(W) {
  10 * log10(W)
}


W2dBm <- function(W) {
  dBW2dBm(W2dBW(W))
}


dBm2dBW <- function(dBm) {
  dBm - 30
}


dBm2W <- function(dBm) {
  dBW2W(dBm - 30)
}


### Signal Parameter plots (theoretical)
########################################

sig_param_plots <- function(param.df, range.max = 20000, base_size = 11) {
  distance <- fill <- xmin <- xmax <- ymin <- ymax <- NULL
  
  # define the helper parameters for the plot data frame creation
  cell.kind.unique <- n_distinct(param.df$cell.kind)
  range.total <- rep(seq(10, range.max, by = 10), cell.kind.unique)
  length.range.total <- rep(length(range.total) / cell.kind.unique, cell.kind.unique)
  
  # use the helpers to construct plot dataframe and join with the input params
  df <- tibble(cell.kind = factor(rep(param.df$cell.kind, length.range.total)),
               distance = range.total) %>% 
    left_join(param.df, by = "cell.kind") %>% 
    mutate(dBm = W2dBm(W)) %>% 
    mutate(distance.log10 = log10(distance)) %>% 
    mutate(dBm = distance2dB(distance, ple, W)) %>% 
    mutate(sig.dom = db2s(dBm, 
                          midpoint = midpoint, 
                          steepness = steepness)) %>% 
    mutate(below.dominance.th = case_when(sig.dom >= dominance.th ~ "Above", 
                                          sig.dom < dominance.th ~ "Below"))
  
  minor.breaks <- rep(1:9, 21) * (10^rep(-10:10, each = 9))
  
  
  strength.distance.plot <- ggplot(df, aes(x = distance, y = dBm)) + 
    geom_line(aes(color = label), size = 1.4) +
    scale_x_log10(labels = scales::trans_format("log10", 
                                                scales::math_format(10^.x)),
                  minor_breaks = minor.breaks) +
    annotation_logticks(sides = "b") +
    labs(title = "Distance vs. Signal strength", 
         x = "log10(Distance (m))",
         y = "Signal strength (dBm)",
         color = "Cell Kind") +
    theme_bw(base_size = base_size) + 
    theme(panel.grid.major = element_line("grey85"))
  
  strength.dominance.plot <- ggplot(df) + 
    geom_line(aes(x = dBm, y = sig.dom, color = label, 
                  alpha = below.dominance.th, linetype = below.dominance.th), size = 1.4) +
    scale_alpha_discrete(range = c(1, 0.4), guide = F) +
    scale_linetype_discrete() +
    labs(title = "Signal dominance vs. strength", 
         x = "Signal strength (dBm)",
         y = "Signal dominance",
         color = "Cell Kind",
         linetype = "Dominance Threshold") +
    guides(color = guide_legend(order = 1), 
           linetype = guide_legend(order = 2)) +
    theme_bw(base_size = base_size) + 
    theme(panel.grid.major = element_line("grey85"))
  
  dominance.distance.plot <- ggplot(df) + 
    geom_line(aes(x = distance, y = sig.dom, color = label, 
                  alpha = below.dominance.th, linetype = below.dominance.th), size = 1.4) +
    scale_alpha_discrete(range = c(1, 0.4), guide = F) +
    scale_linetype_discrete() +
    labs(title = "Distance vs. Signal dominance", 
         x = "Distance (m)",
         y = "Signal dominance",
         color = "Cell Kind",
         linetype = "Dominance Threshold") +
    guides(color = guide_legend(order = 1), 
           linetype = guide_legend(order = 2)) +
    theme_bw(base_size = base_size) + 
    theme(panel.grid.major = element_line("grey85"))
  
  table.gen <- param.df %>% 
    dplyr::select(label, W, ple) %>% 
    mutate(`Transmit Power (dBm)` = W2dBm(W)) %>% 
    dplyr::select(-W, `Path Loss Exponent (pos)` = ple) %>% 
    pivot_longer(-label, names_to = "Gen.Parameter", values_to = "Value") %>% 
    pivot_wider(id_cols = Gen.Parameter, names_from = label, values_from = Value) %>% 
    tableGrob(rows = NULL)
  
  table.mod <- param.df %>% 
    dplyr::select(label, Midpoint = midpoint, Steepness = steepness, `Dominance Threshold` = dominance.th) %>% 
    pivot_longer(-label, names_to = "Mod.Parameter", values_to = "Value") %>% 
    pivot_wider(id_cols = Mod.Parameter, names_from = label, values_from = Value) %>% 
    tableGrob(rows = NULL)
  
  final <- arrangeGrob(strength.distance.plot, table.gen, table.mod, strength.dominance.plot, dominance.distance.plot,
                       layout_matrix = rbind(c(1, 1, 2, 2),
                                             c(1, 1, 3, 3),
                                             c(4, 4, 5, 5),
                                             c(4, 4, 5, 5)))
  
  return(list(final = final,
              strength.distance.plot = strength.distance.plot,
              strength.dominance.plot = strength.dominance.plot,
              dominance.distance.plot = dominance.distance.plot,
              table.gen = table.gen,
              table.mod = table.mod))
  
}

### Type plot
#############

type_plot <- function(data) {
  
  data %>% 
    ggplot() +
    geom_sf(aes(fill = type, color = type)) +
    scale_fill_manual(values = c("white", "grey", "black"),
                      breaks = c("Hole", "Rural", "Urban")) +
    scale_color_manual(values = c("white", "grey", "black"),
                       breaks = c("Hole", "Rural", "Urban"),
                       guide = F) +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    labs(title = "Simulated area by type",
         fill = "Type")
}


### Pop dist plot
#############

pop_plot <- function(data) {
  
  data %>% 
    ggplot() +
    geom_sf(aes(fill = pop, color = pop)) +
    scale_fill_gradient(low = "white", high = "black") +
    scale_color_gradient(low = "white", high = "black", 
                         guide = F) +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    labs(title = "Simulated area by population density",
         fill = "Population")
} 


### pop summary results
###########

pop_summary_results <- function(data) {
  
  data %>% 
    st_drop_geometry() %>% 
    group_by(type) %>% 
    summarise(n.type = n(),
              prop.type = n() / length(data$area.sf$tile.id),
              mean.pop = mean(pop), 
              sd.pop = sd(pop), 
              min.pop = min(pop), 
              max.pop = max(pop),
              sum.pop = sum(pop)) %>% 
    mutate_if(is.numeric, round, 2)
    # left_join(pop.dist.df, by = "type")
  
  }
  
  

### Pop Density plots 
#############

# custom_ecdf_prep <- function(data) {
#   dat <- data %>% 
#     mutate(pop.plot = pop + 1) %>%  
#     arrange(pop.plot) %>%  
#     mutate(prob = 1 / n()) %>%  
#     mutate(cum.prob = cumsum(prob)) %>%  
#     mutate(cum.prob.comp = 1 - cum.prob) %>%  
#     mutate(log10.cum.prob.comp = log10(cum.prob.comp)) %>% 
#     mutate(log10.pop = log10(pop.plot)) %>%  
#     mutate(cum.prob.comp = 1 - cum.prob)
#   
#   return(dat)
# }

density_plots <- function(data) {
  
  custom_ecdf_prep <- function(data) {
    dat <- data %>% 
      mutate(pop.plot = pop + 1) %>%  
      arrange(pop.plot) %>%  
      mutate(prob = 1 / n()) %>%  
      mutate(cum.prob = cumsum(prob)) %>%  
      mutate(cum.prob.comp = 1 - cum.prob) %>%  
      mutate(log10.cum.prob.comp = log10(cum.prob.comp)) %>% 
      mutate(log10.pop = log10(pop.plot)) %>%  
      mutate(cum.prob.comp = 1 - cum.prob)
    
    return(dat)
  }
  
  minor.breaks <- rep(1:9, 21) * (10^rep(-10:10, each = 9))
  
  pop.emp.dist <- data %>% 
    ungroup() %>% 
    custom_ecdf_prep()
  
  ECCDF.df <- pop.emp.dist %>% 
    dplyr::select(log10.cum.prob.comp, log10.pop, type) %>%
    mutate(log10.cum.prob.comp = round(log10.cum.prob.comp, 3)) %>% # effective plot sample --> faster plotting excluding overplot
    distinct()
  
  ECDF.df <- pop.emp.dist %>% 
    dplyr::select(cum.prob.comp, pop.plot, type) %>%
    mutate(cum.prob.comp = round(cum.prob.comp, 3)) %>% # effective plot sample --> faster plotting excluding overplot
    distinct()
  
  ECCDF.pop.plot <- ECDF.df %>%   
    ggplot() + 
    geom_point(aes(x = pop.plot, y = cum.prob.comp
                   # color = type
                   )) + 
    # geom_hline(yintercept = -0.3010300, linetype = "dotted") + 
    # geom_hline(yintercept = -1, linetype = "dotted") + 
    # geom_text(x = 1.8, y = -0.15, label = "50% of the data") + 
    # geom_text(x = 1.8, y = -0.8, label = "90% of the data") + 
    scale_color_ptol() + 
    scale_y_log10(labels = scales::trans_format("log10", 
                                                scales::math_format(10^.x)),
                  minor_breaks = minor.breaks) +
    scale_x_log10(labels = scales::trans_format("log10", 
                                                scales::math_format(10^.x)),
                  minor_breaks = minor.breaks) +
    annotation_logticks(sides = "lb") +
    labs(y = "ECCDF", x = "Mobile phones",  
         colour = "") + 
    theme(legend.position = "bottom",
          text = element_text(size = 13))
  
  ECDF.pop.plot <- ECDF.df %>%   
    ggplot() + 
    geom_point(aes(x = pop.plot, y = cum.prob.comp
                   # color = type
                   )) + 
    scale_color_ptol() +
    xlim(0, 30) +
    labs(title = "", y = "", x = "") +
    theme(legend.position = "none",
          plot.margin = unit(c(-0.5, 0, 0, -0.5), "cm")) 
  
  
  combined <- ECCDF.pop.plot +
    annotation_custom(ggplotGrob(ECDF.pop.plot), 
                      xmin = 0, xmax = 1.5, 
                      ymin = -3, ymax = -1.5)
  
  # return(list(ECDF = ECDF.pop.plot,
  #             ECCDF = ECCDF.pop.plot))
  return(combined)
  
}


scatter_density <- function(point, estimator.name){
  
  
  
  
  base.data.ls <- point %>% 
    filter(estimator %in% estimator.name) %>% 
    group_split(scale)
  
  complete.plot <- base.data.ls %>% 
    map(~ggplot(., aes(x = pop, y = estimate)) +
          geom_pointdensity(size = 0.4) +
          scale_color_viridis(guide = F) +
          geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
          # lims(x = c(0, 5), y = c(0, 5)) +
          labs(title = "",
               x = "",
               y = "") +
          facet_wrap(vars(scale)) +
          theme(plot.margin = unit(c(-1, -1, -1, -1), "lines"),
                axis.text = element_text(size = 4),
                strip.text.x = element_text(size = 8)))
  
  complete.plot.final <- arrangeGrob(complete.plot[[1]], complete.plot[[2]], complete.plot[[3]], complete.plot[[4]],
                                     # padding = 2,
                                     layout_matrix = rbind(c(1, 2, 3, 4)))
  # ggsave("2d_density/d.png", complete.plot.final, device = "png")
  # grid.arrange(complete.plot[[1]], complete.plot[[2]], complete.plot[[3]], complete.plot[[4]],
  #             padding = 0,
  #             layout_matrix = rbind(c(1, 2),
  #                                   c(3, 4)))
  
  
  zoom.xy <- base.data.ls %>% 
    map(~ggplot(., aes(x = pop, y = estimate)) +
          geom_pointdensity(size = 0.4) +
          scale_color_viridis(guide = F) +
          geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
          lims(x = c(0, 5), y = c(0, 5)) +
          labs(title = "",
               x = "",
               y = "") +
          facet_wrap(vars(scale)) +
          theme(plot.margin = unit(c(-1, -1, -1, -1), "lines"),
                axis.text = element_text(size = 4),
                strip.text.x = element_text(size = 8)))
  
  zoom.xy.final <- arrangeGrob(zoom.xy[[1]], zoom.xy[[2]], zoom.xy[[3]], zoom.xy[[4]],
                               # padding = 0,
                               top = textGrob("Zoom on XY", gp = gpar(fontsize = 10)),
                               layout_matrix = rbind(c(1, 2, 3, 4)))
  
  minor.breaks <- rep(1:9, 21) * (10^rep(-10:10, each = 9))
  
  log.both <- base.data.ls %>% 
    map(~ggplot(., aes(x = log10(pop + 1), y = log10(estimate + 1))) +
          geom_pointdensity(size = 0.4) +
          scale_color_viridis(guide = F) +
          geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
          scale_y_log10(labels = scales::trans_format("log10", 
                                                      scales::math_format(10^.x)),
                        minor_breaks = minor.breaks) +
          scale_x_log10(labels = scales::trans_format("log10", 
                                                      scales::math_format(10^.x)),
                        minor_breaks = minor.breaks) +
          annotation_logticks(sides = "lb") +
          labs(title = "",
               x = "",
               y = "") +
          facet_wrap(vars(scale)) +
          theme(plot.margin = unit(c(-1, -1, -1, -1), "lines"),
                axis.text = element_text(size = 4),
                strip.text.x = element_text(size = 8)))
  
  log.both.final <- arrangeGrob(log.both[[1]], log.both[[2]], log.both[[3]], log.both[[4]],
                                # padding = 0,
                                top = textGrob("Joint Density log10", gp = gpar(fontsize = 10)),
                                layout_matrix = rbind(c(1, 2, 3, 4)))
  
  scatter.density <- arrangeGrob(complete.plot.final, zoom.xy.final, log.both.final,
                                 # padding = 0,
                                 layout_matrix = rbind(c(1),
                                                       c(2),
                                                       c(3)))
  
  return(scatter.density)
  
}




P_equalizer <- function(P.long.df, signal.strength.comb.dt) {
  
  P.model.equal.start <- P.long.df %>% 
    left_join(signal.strength.comb.dt, by = c("tile.id.num" = "rid", "cell")) 
  
  P.model.summary <- P.model.equal.start %>% 
    mutate(dummy = 1) %>% 
    group_by(tile.id.num) %>% 
    summarise(TH.05 = sum(dummy[s > 0.5]),
              TH.025 = sum(dummy[s > 0.25 & s <= 0.5])) %>% 
    mutate(TH.025 = case_when(TH.05 != 0 ~ 0,
                              TRUE ~ TH.025))
  
  under.any.TH <- P.model.summary %>% 
    filter(TH.05 == 0 & TH.025 ==  0) %>% 
    distinct(tile.id.num) %>% 
    deframe()
  
  TH.05 <- P.model.summary %>% 
    filter(TH.05 >= 1 & TH.025 == 0) %>% 
    distinct(tile.id.num) %>% 
    deframe()
  
  TH.025 <- P.model.summary %>% 
    filter(TH.05 == 0 & TH.025 >= 1) %>% 
    distinct(tile.id.num) %>% 
    deframe()
  
  
  # calculate new p for each of the three groups
  P.model.equal.u.TH <- P.model.equal.start %>% 
    filter(tile.id.num %in% under.any.TH) %>% 
    mutate(pij.new = 1) %>% 
    group_by(tile.id.num) %>% 
    mutate(pij.equal = pij.new / sum(pij.new)) %>% 
    ungroup()
  
  P.model.equal.TH.05 <- P.model.equal.start %>% 
    filter(tile.id.num %in% TH.05) %>%
    filter(s > 0.5) %>%  
    mutate(pij.new = 1) %>% 
    group_by(tile.id.num) %>% 
    mutate(pij.equal = pij.new / sum(pij.new)) %>% 
    ungroup()
  
  P.model.equal.TH.025 <- P.model.equal.start %>% 
    filter(tile.id.num %in% TH.025) %>% 
    filter(s > 0.25) %>%  
    mutate(pij.new = 1) %>% 
    group_by(tile.id.num) %>% 
    mutate(pij.equal = pij.new / sum(pij.new)) %>% 
    ungroup()
  
  P.model.equal.final <- bind_rows(P.model.equal.u.TH, P.model.equal.TH.05, P.model.equal.TH.025)
  
  
  
  return(P.model.equal.final)
  
}


custom_ecdf_prep <- function(data) {
  dat <- data %>% 
    mutate(pop.plot = values + 1) %>%  
    arrange(pop.plot) %>%  
    mutate(prob = 1 / n()) %>%  
    mutate(cum.prob = cumsum(prob)) %>%  
    mutate(cum.prob.comp = 1 - cum.prob) %>%  
    mutate(log10.cum.prob.comp = log10(cum.prob.comp)) %>% 
    mutate(log10.pop = log10(pop.plot)) %>%  
    mutate(cum.prob.comp = 1 - cum.prob)
  
  return(dat)
}


## dev to cell
c.vec.sampler <- function(x) {
  data.table(sample(x = as.character(x$cell), size = mean(x$pop),
                    replace = T, prob = x$pij))
}
