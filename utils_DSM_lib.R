# ////////////////////////////////////////////////////////////////////////////
#   // ***** LOAD LANDSAT COMPOSITE ***** //
#   ////////////////////////////////////////////////////////////////////////////
#
# Assumes the image is a Landsat image

maskCloudsAndSuch = function(img, cloudThresh  = 20){
  #Bust clouds
  cs = ee$Algorithms$Landsat$simpleCloudScore(img)$select('cloud')$gt(cloudThresh)
  #Make sure all or no bands have data
  numberBandsHaveData = img$mask()$reduce(ee$Reducer$sum())
  allOrNoBandsHaveData = numberBandsHaveData$eq(0)$Or(numberBandsHaveData$gte(7))

  #If it's Landsat 5- defringe by nibbling away at the fringes
  allBandsHaveData = allOrNoBandsHaveData

  #.focal_min(1,'square','pixels',8)

  #Make sure no band is just under zero
  allBandsGT = img$reduce(ee$Reducer$min())$gt(-0.001)
  return(img$mask(img$mask()$And(cs$Not())$And(allBandsHaveData)$And(allBandsGT)))
}


getComp = function(compositeArea, year, compositingPeriod, startJulian, endJulian, exportToDrive, roiName, crs) {

  if(class(compositeArea)[[1]] == "sf"){
  compositeArea <- compositeArea %>% sf::st_transform(crs = 4326, proj4string = "+init=epsg:4326")
  bb <- sf::st_bbox(compositeArea)
  compositeArea <- ee$Geometry$Rectangle(bb)
  }

    # // **** USER EDITABLE VARIABLES **** //
    #   //Composite parameters
    buffer_distance = 1
    cloudThresh = 20 # Specify cloud threshold (0-100)- lower number masks out more clouds
    possibleSensors = c('L5','L7','L8') # Specify which sensors to pull from- supports L5,L7, and L8

    reducerPercentile = 50 # Reducer for compositing
    #
    # //////////////////////////////////////////////////////
    #   //Globals
    # //////////////////////////////////////////////////////
    #   // Establish composite area based on input geometry.

    reducer = ee$Reducer$percentile(list(reducerPercentile))

    #reducer for compositing
    prioritizeL5 = FALSE # Binary true or false for prioritizing L5, L8, then L7 in compositing
    studyArea = compositeArea$buffer(buffer_distance)

    fB=studyArea
    # Map$centerObject(studyArea)
    # //////////////////////////////////////////////////////
    #   //Use data mask from Hansen's Global Forest Change as a water mask
    #Use data mask from Hansen's Global Forest Change as a water mask

    forestChangeImage = ee$Image('UMD/hansen/global_forest_change_2019_v1_7')
    mskW = forestChangeImage$select('datamask')
    mskW = mskW$eq(1)
    # Run a focal_mode convolution on the image.
    maskFocalMode = mskW$focal_mode()
    # Further smooth the image via focal_max
    watermask = maskFocalMode$focal_max(5, "square", "pixels", 5 )

    #Band combinations for each sensor corresponding to final selected corresponding bands
    sensor_band_dict = ee$Dictionary(list(L8 = ee$List(c(1,2,3,4,5,9,6)),
                                          L7 = ee$List(c(0,1,2,3,4,5,7)),
                                          L5 = ee$List(c(0,1,2,3,4,5,6)),
                                          L4 = ee$List(c(0,1,2,3,4,5,6))))

    #band names
    bandNames = ee$List(c('blue','green','red','nir','swir1','temp','swir2'))
    STD_NAMES = c('blue','green','red','nir','swir1','temp','swir2')
    bandNumbers = c(0,1,2,3,4,5,6)




  #Define dates
  y1Image = year
  y2Image = year + compositingPeriod

  startDate = ee$Date$fromYMD(ee$Number(year),1,1)$advance(startJulian,'day')
  endDate = ee$Date$fromYMD(ee$Number(year)$add(ee$Number(compositingPeriod)),1,1)$advance(endJulian,'day')


  #Get the images for composite and shadow model
  if(!is.null(possibleSensors)){
    l5s = getCollection('L5',startDate,endDate,startJulian,endJulian)
  }
  else{l5s = getCollection('L5',ee$Date('1000-01-01'),ee$Date('1001-01-01'),0,365)}

  if(!is.null(possibleSensors)){
    l7s = getCollection('L7',startDate,endDate,startJulian,endJulian)
  }
  else{l7s = getCollection('L7',ee$Date('1000-01-01'),ee$Date('1001-01-01'),0,365)}

  if(!is.null(possibleSensors)){
    l8s = getCollection('L8',startDate,endDate,startJulian,endJulian)
  }
  else{l8s = getCollection('L8',ee$Date('1000-01-01'),ee.Date('1001-01-01'),0,365)}

  ls = ee$ImageCollection(l5s$merge(l7s)$merge(l8s))
  composite = ls$reduce(reducer)$select(bandNumbers,bandNames)$mask(watermask)
  composite = composite$mask(composite$mask()$And(watermask)$clip(studyArea))

  #fullName = paste0(roiName,'_',y1Image$toString(),'_' ,y2Image$toString(),'_',startJulian$toString(),'_',endJulian$toString(),'_Composite')

  #Set up our final composite with bands we'd like to include
  composite = composite$select(c('blue','green','red','nir','swir1','swir2'))$multiply(10000)$int16()$clip(fB)

  return(composite)

}

# Basic shadow masking using sum of specified bands
# Tends to include hill shadows and water

maskShadows = function(img){
  shadowThresh = 0.1
shadowSumBands = c('nir','swir1','swir2')
ss = img$select(shadowSumBands)$reduce(ee$Reducer$sum())
  return(img$mask(img$mask()$And(ss$gt(shadowThresh))))
}

# Function to handle empty collections that will cause subsequent processes to fail
# If the collection is empty, will fill it with an empty image
fillEmptyCollections = function(inCollection,dummyImage){
  dummyCollection = ee$ImageCollection(dummyImage$mask(ee$Image(0)))
  imageCount = inCollection$toList(1)$length()
  return(ee$ImageCollection(ee$Algorithms$If(imageCount$gt(0),inCollection,dummyCollection)))

}
#Helper function to get images from a specified sensor

getCollection = function(sensor,startDate,endDate,startJulian,endJulian){

  #Names of collections to look in
  #Add _L1T for L1T imagery
  #TOA is computed on both the L1G or L1T
  collection_dict = list(L8 = "LANDSAT/LC08/C01/T1_TOA",
                         L7 = "LANDSAT/LE07/C01/T1_TOA",
                         L5 = "LANDSAT/LT05/C01/T1_TOA",
                         L4 = "LANDSAT/LT04/C01/T1_TOA")
  sensor_band_dict = ee$Dictionary(list(L8 = ee$List(c(1,2,3,4,5,9,6)),
                                        L7 = ee$List(c(0,1,2,3,4,5,7)),
                                        L5 = ee$List(c(0,1,2,3,4,5,6)),
                                        L4 = ee$List(c(0,1,2,3,4,5,6))))
  bandNames = ee$List(c('blue','green','red','nir','swir1','temp','swir2'))
  collectionName = collection_dict[sensor][[1]]

  #Start with an un-date-confined collection of iamges
  WOD = ee$ImageCollection(collectionName)

  #Pop off an image to serve as a template if there are no images in the date range
  dummy = ee$Image(WOD$first())

  #Filter by the dates
  ls = WOD$filterDate(startDate,endDate)$filter(ee$Filter$calendarRange(startJulian,endJulian))

  #Fill the collection if it's empty
  ls = fillEmptyCollections(ls,dummy)

  #Clean the collection up- clouds, fringes....
  ls = ls$map(rgee::ee_utils_pyfunc(maskCloudsAndSuch))$select(sensor_band_dict$get(sensor),bandNames)$map(rgee::ee_utils_pyfunc(maskShadows))
  return(ls)
}
#
# ////////////////////////////////////////////////////////////////////////////
#   // ***** LOAD SPECTRAL BANDS ***** //
#   ////////////////////////////////////////////////////////////////////////////

getNDVI = function (inImage) {
   #Calculate NDVI = (nir - red) / (nir + red)
  # bring in composite from ex_1 and calc NDVI

  ndvi <- inImage$normalizedDifference(c('nir','red'))$rename('ndvi')

  # now add this to composite as band

  composite_ndvi <- inImage$addBands(ndvi)
  }


# /////////////////////////////////////////////////////////////////////////////////////
#   // 2.2 Calculate mineral/geologic indices
getMineralIndices = function(inImage, crs){

  # // Clay Minerals = swir1 / swir2
  clayIndex = inImage$select('swir1')$divide(inImage$select('swir2'))$rename('clayIndex')

  # // Ferrous Minerals = swir / nir
   ferrousIndex = inImage$select('swir1')$divide(inImage$select('nir'))$rename('ferrousIndex')

  # // Carbonate Index = (red - green) / (red + green)
   carbonateIndex = inImage$normalizedDifference(c('red','green'))$rename('carbonateIndex')

  # // Rock Outcrop Index = (swir1 - green) / (swir1 + green)
   rockOutcropIndex = inImage$normalizedDifference(c('swir1','green'))$rename('rockOutcropIndex')

  # // Add bands
  outStack = inImage$addBands(c(clayIndex, ferrousIndex, carbonateIndex, rockOutcropIndex))

  # // if(mineralIndices.contains('clayIndex')){
  #   //   composite = composite.addBands(clayIndex);
  #   // } else if(mineralIndices.contains('ferrousIndex')) {
  #     //   composite = composite.addBands(ferrousIndex);
  #     // } else if(mineralIndices.contains('carbonateIndex')) {
  #       //   composite = composite.addBands(carbonateIndex);
  #       // } else if(mineralIndices.contains('rockOutcropIndex')) {
  #         //   composite = composite.addBands(rockOutcropIndex);
  #         // }
  #
  # // Display mineral indices
  # //////////////////////////////////////////////////////////////////////////////////
  #   // Calculate visual params using min and max of image
  #
  # // get dictionary of min values for each band
 #  Mineralviz_min = (outStack.reduceRegion({
 #    reducer: ee.Reducer.min(),
 #    scale: 10,
 #    crs: crs,
 #    bestEffort: true
 #  }));
 #  # get dictionary of max values for each band
 # Mineralviz_max = (outStack.reduceRegion({
 #    reducer: ee.Reducer.max(),
 #    scale: 10,
 #    crs: crs,
 #    bestEffort: true
 #  }));
 #  // Add to map
  # // Map.addLayer(outStack.select("clayIndex"), {min: Mineralviz_min.getNumber('clayIndex').getInfo(), max: Mineralviz_max.getNumber('clayIndex').getInfo()}, "clayIndex");
  # // Map.addLayer(outStack.select("ferrousIndex"), {min: Mineralviz_min.getNumber('ferrousIndex').getInfo(), max: Mineralviz_max.getNumber('ferrousIndex').getInfo()}, "ferrousIndex");
  # // Map.addLayer(outStack.select("carbonateIndex"), {min: Mineralviz_min.getNumber('carbonateIndex').getInfo(), max: Mineralviz_max.getNumber('carbonateIndex').getInfo()}, "carbonateIndex");
  # // Map.addLayer(outStack.select("rockOutcropIndex"), {min: Mineralviz_min.getNumber('rockOutcropIndex').getInfo(), max: Mineralviz_max.getNumber('rockOutcropIndex').getInfo()}, "rockOutcropIndex");

}

# /////////////////////////////////////////////////////////////////////////////////////
#   // 2.3 Tasseled cap image transformation
#
# // add tassel cap function
getTasseledCap = function(image) {

  bands = ee$List(c('blue','green','red','nir','swir1','swir2'))

  # //Crist 1985 coeffs - TOA refl (http://www.gis.usu.edu/~doug/RS5750/assign/OLD/RSE(17)-301.pdf)
  coefficients = ee$Array(list(c(0.2043, 0.4158, 0.5524, 0.5741, 0.3124, 0.2303),
                               c(-0.1603, -0.2819, -0.4934, 0.7940, -0.0002, -0.1446),
                               c(0.0315, 0.2021, 0.3102, 0.1594, -0.6806, -0.6109),
                               c(-0.2117, -0.0284, 0.1302, -0.1007, 0.6529, -0.7078),
                               c(-0.8669, -0.1835, 0.3856, 0.0408, -0.1132, 0.2272),
                               c(0.3677, -0.8200, 0.4354, 0.0518, -0.0066, -0.0104)))
  # // Make an Array Image, with a 1-D Array per pixel.
   arrayImage1D = image$select(bands)$toArray()

  # // Make an Array Image with a 2-D Array per pixel, 6x1.
   arrayImage2D = arrayImage1D$toArray(1)

   componentsImage = ee$Image(coefficients)$matrixMultiply(arrayImage2D)$arrayProject(list(0))$arrayFlatten(list(c('brightness', 'greenness', 'wetness', 'fourth', 'fifth', 'sixth')))$float()
meta = componentsImage$getInfo()
  return(image$addBands(componentsImage))
}


getTC = function(inImage, crs){

  # // Apply Tasseled cap transformation.
   TCall = getTasseledCap(inImage)
  # // Choose desired bands -- New composite will add 3 new layers: brightness, greenness, wetness
  TCBands = TCall$select(c("brightness", "greenness", "wetness"))
  # // Stack onto Landsat image composite
   outStack = inImage$addBands(TCBands)

  return(outStack)
}

# // Function to calculate all spectral indices/image enhancement (NDVI, mineral, and tasseled cap transformation)
getSpectralIndices = function (inImage, crs){
  stackNDVI = getNDVI(inImage)
  stackNDVIMinerals = getMineralIndices(stackNDVI, crs)
  stackNDVIMineralsTC = getTC(stackNDVIMinerals, crs)

  return(stackNDVIMineralsTC)
}

#
# ////////////////////////////////////////////////////////////////////////////
#   // ***** LOAD TOPO DERIVS ***** //
#   ////////////////////////////////////////////////////////////////////////////

getTopo = function(compositeArea, image){

  if(class(compositeArea)[[1]] == "sf"){
    compositeArea <- compositeArea %>% sf::st_transform(crs = 4326, proj4string = "+init=epsg:4326")
    bb <- sf::st_bbox(compositeArea)
    compositeArea <- ee$Geometry$Rectangle(bb)
  }

    # //set some params
    resolution = 10
    # //var projection = 'EPSG:5070';//26912';

topoBands = getNEDTopography(compositeArea)$clip(compositeArea)
outStack = image$addBands(topoBands)


return(outStack)

}
 # // Function to add USGS 1/3 arc second topography and derive slope, aspect, & hillshade

    getNEDTopography = function(compositeArea){

      if(class(compositeArea)[[1]] == "sf"){
        compositeArea <- compositeArea %>% sf::st_transform(crs = 4326, proj4string = "+init=epsg:4326")
        bb <- sf::st_bbox(compositeArea)
        compositeArea <- ee$Geometry$Rectangle(bb)
      }
      #Import NED elevation data
      elevation = ee$Image('USGS/NED')

      # Clip to area of interest
      elevation = elevation$clip(compositeArea)

      # Calculate slope and aspect
      topo = ee$Algorithms$Terrain(elevation)

      # Calculate % slope
      slopeDeg = topo$select('slope')
      slopeRads = slopeDeg$multiply(base::pi)$divide(ee$Number(180))
      slopeTan = slopeRads$tan()
      slopePCT = slopeTan$multiply(ee$Number(100))$rename('slopePCT')

      # Add 8-direction aspect
      aspect = topo$select('aspect')
      aspectRad = aspect$multiply(base::pi)$divide(180)
      aspectSin = aspectRad$sin()$rename('sin')
      aspectCos = aspectRad$cos()$rename('cos')
      aspect_8 = (aspect$multiply(8)$divide(360))$add(1)$floor()$uint8()$rename('aspect_8')

      # Option: add increased bit depth to values.
      #aspectSin = aspectSin.multiply(10000).int32()

      #Add 3 equally-spaced sun azimuth hillshades
      hill_1 = ee$Terrain$hillshade(elevation,30)$rename('hill_1')
      hill_2 = ee$Terrain$hillshade(elevation,150)$rename('hill_2')
      hill_3 = ee$Terrain$hillshade(elevation,270)$rename('hill_3')

      #Compile selected topography bands
      topo = topo$select('elevation')$addBands(c(slopePCT,
                                                 aspectSin,
                                                 aspectCos))

      return(topo)
    }
# ////////////////////////////////////////////////////////////////////////////
#                     // ***** LOAD CLIMATE DATA ***** //
# ////////////////////////////////////////////////////////////////////////////

  getClimate = function(compositeArea, year, compositingPeriod, crs){

  clim = ee$ImageCollection("OREGONSTATE/PRISM/AN81d")


  if(class(compositeArea)[[1]] == "sf"){
    compositeArea <- compositeArea %>% sf::st_transform(crs = 4326, proj4string = "+init=epsg:4326")
    bb <- sf::st_bbox(compositeArea)
    compositeArea <- ee$Geometry$Rectangle(bb)
  }

  summer_startMonth = 6 #//Calendar month 1-12
  summer_startDay = 20 #//Day of month
  summer_length = 94 #//Number of days of the season
  winter_startMonth = 12
  winter_startDay = 21
  winter_length = 89


  # // Calculate multi-year mean of total seasonal precipitation: summer, winter, seasonal difference
  winter_precip = season_climate(clim,'ppt', winter_startMonth, winter_startDay, winter_length,  year, compositingPeriod)$mean()$rename("Winter Precip")
  summer_precip = season_climate(clim,'ppt', summer_startMonth, summer_startDay, summer_length,  year, compositingPeriod)$mean()$rename("Summer Precip")

  # // Calculate multi-year mean of mean seasonal temperature: summer, winter, seasonal difference
  winter_tmean = season_climate(clim,'tmean', winter_startMonth, winter_startDay, winter_length,  year, compositingPeriod)$mean()$rename("Winter Temp")
  summer_tmean = season_climate(clim,'tmean', summer_startMonth, summer_startDay, summer_length,  year, compositingPeriod)$mean()$rename("Summer Temp")

  # // Calculate seasonal climate differences
  difference_precip = summer_precip$subtract(winter_precip)$rename("Seasonal Precip Difference")
  difference_tmean = summer_tmean$subtract(winter_tmean)$rename("Seasonal Temp Difference")

  # // Stack climate bands together
  climateBands = winter_precip$addBands(c(summer_precip,
                                          difference_precip,
                                          winter_tmean,
                                          summer_tmean,
                                          difference_tmean))


  # // Clip to ROI boundaries
  climateBands = climateBands$clip(compositeArea)

  return(climateBands)

  }



  # // Function inputs are band, seasonal dates, season length
  season_climate = function(image, band, month_start, day_start, season_length, year, compositingPeriod){

    clim = image
    # // Specify year range for calculating multi-year average
    # // Same as year range used for image composite
    start_year = year
    end_year = year + compositingPeriod
    years = ee$List$sequence(start_year, end_year)

    # // Calculate seasonal climate values for each year
    return(ee$ImageCollection(years$map(rgee::ee_utils_pyfunc(function(yr){

      # //////////////////////////////////////////////////////////////////////
      #   // For seasonal temp ('tmean' band) -- caluculates MEAN of mean daily temperature per season, per year
      # //////////////////////////////////////////////////////////////////////
      if(band == 'tmean'){
        bandName = clim$select('tmean')
        # // filter to dates of interest
        return(bandName$filter(ee$Filter$date(
          ee$Date$fromYMD(yr,month_start,day_start),
          # // Adds season length to starting date to get ending date
          ee$Date$fromYMD(yr,month_start,day_start)$advance(season_length,'day')))$mean()$set('system:time_start', ee$Date$fromYMD(yr, month_start,1)$millis()))
      }
      #
      # //////////////////////////////////////////////////////////////////////
      #   // For seasonal precip ('ppt' band) -- calculates SUM of total daily precip precipitation per season, per year
      # //////////////////////////////////////////////////////////////////////
      if(band == 'ppt'){
        bandName = clim$select('ppt')
        # // filter to dates of interest
        return(bandName$filter(ee$Filter$date(
          ee$Date$fromYMD(yr,month_start,day_start),
          # // Adds season length to starting date to get ending date
          ee$Date$fromYMD(yr,month_start,day_start)$advance(season_length,'day')))$sum()$set('system:time_start', ee$Date$fromYMD(yr, month_start,1)$millis()))
      }

    })
    ))) #// end iterating function over years
  }
