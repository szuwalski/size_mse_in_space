

costs_to_fish <-
  function(land_mask,
           port_lat,
           port_lon,
           lat,
           lon){
    
    #=====================================================
    #===CALCULATE COSTS TO A PATCH FROM THE PORT FOR COSTS
    #=====================================================
    
    
    cost <- raster(
      nrow = length(lat),
      ncol = length(lon),
      xmn = min(lon),
      xmx = max(lon),
      ymn = min(lat),
      ymx = max(lat),
      crs = "+proj=utm"
    )
    
    
    # have to think this value of cost 
    # need to talk with cody
    cost[] <- 1
    for (x in 1:nrow(cost))
      for (y in 1:ncol(cost))
        if (land_mask[x, y] == 0)
          cost[x, y] <- 100000
    
    trCost <- transition(1 / cost, min, directions = 8)
    trCostC <- geoCorrection(trCost, type = "c")
    trCostR <- geoCorrection(trCost, type = "r")
    
    ## Create three points (representing three points in time series)
    pnts <- cbind(x = c(port_lon, -155, -160),
                  y = c(port_lat, 52, 58))
    costDistance(trCostC, pnts)
    
    ## Display results for one set of points
    plot(cost)
    plot(SpatialPoints(pnts),
         add = TRUE,
         pch = 20,
         col = "red")
    plot(shortestPath(trCostC, pnts[1, ], pnts[2, ], output = "SpatialLines"),
         add = TRUE)
    plot(shortestPath(trCostC, pnts[1, ], pnts[3, ], output = "SpatialLines"),
         add = TRUE)
    plot(shortestPath(trCostC, pnts[2, ], pnts[3, ], output = "SpatialLines"),
         add = TRUE)
    
    #==find the distance from harbor to all points
    distance_map <- matrix(ncol = length(lon), nrow = length(lat))
    colnames(distance_map) <- lon
    rownames(distance_map) <- (lat)
    for (x in 1:length(lon))
      for (y in 1:length(lat))
      {
        if (land_mask[y, x] != 0)
        {
          pts <-
            cbind(x = c(port_lon, as.numeric(colnames(distance_map)[x])),
                  y = c(port_lat,  as.numeric(rownames(distance_map)[y])))
          distance_map[y, x] <- costDistance(trCostC, pts[1, ], pts[2, ])
          if (distance_map[y, x] > 10000)
            distance_map[y, x] <- NA
        }
      }
    filled.contour(
      x = lon,
      y = rev(lat),
      g(distance_map * land_mask),
      plot.axes = c(
        map(add = TRUE, fill = T, col = 'grey'),
        points(
          y = port_lat,
          x = port_lon,
          pch = 16,
          col = 'red'
        )
      )
    )
    #write.csv(distance_map,'dist.csv')
    
    
    return(list("distance_map" = distance_map))
  }
