
# load library
library(exploreRGEE)
library(rgee)
library(mapedit)
ee_Initialize()

# create your AOI
shapes <- mapedit::drawFeatures()

# this creates a GEE geom with an sf object
shapes_boundary <- exploreRGEE:::sf_setup(shapes)$geom

#Composite area of interest
compositeArea = shapes_boundary  #use a polygon drawn on the editmap
roiName = 'MT_random' #Give the study area a descriptive name.
buffer_distance = 1 #Distance to buffer composite area in meters. Must be > 0.

#Composite time period
year = 2019 # Start year for composite
startJulian = 152 # Starting Julian Date
endJulian = 252 # Ending Julian date
compositingPeriod = 0 # Number of years into the future to include

#Composite parameters
cloudThresh = 20 # Specify cloud threshold (0-100)- lower number masks out more clouds
possibleSensors = c('L5','L7','L8') # Specify which sensors to pull from- supports L5,L7, and L8

reducerPercentile = 50 # Reducer for compositing

reducer = ee$Reducer$percentile(list(reducerPercentile))

#reducer for compositing
prioritizeL5 = FALSE # Binary true or false for prioritizing L5, L8, then L7 in compositing
studyArea = compositeArea$buffer(buffer_distance)

fB=studyArea
Map$centerObject(studyArea)

 ############# Now all the functions to get a final function getComp() #########################

#Use data mask from Hansen's Global Forest Change as a water mask
forestChangeImage = ee$Image('UMD/hansen/global_forest_change_2019_v1_7')
mskW = forestChangeImage$select('datamask')
mskW = mskW$eq(1)
# Run a focal_mode convolution on the image.
maskFocalMode = mskW$focal_mode()
# Further smooth the image via focal_max
watermask = maskFocalMode$focal_max(5, "square", "pixels", 5 )

#Names of collections to look in
#Add _L1T for L1T imagery
#TOA is computed on both the L1G or L1T
collection_dict = list(L8 = "LANDSAT/LC08/C01/T1_TOA",
                       L7 = "LANDSAT/LE07/C01/T1_TOA",
                       L5 = "LANDSAT/LT05/C01/T1_TOA",
                       L4 = "LANDSAT/LT04/C01/T1_TOA")

#Band combinations for each sensor corresponding to final selected corresponding bands
sensor_band_dict = ee$Dictionary(list(L8 = ee$List(c(1,2,3,4,5,9,6)),
                                      L7 = ee$List(c(0,1,2,3,4,5,7)),
                                      L5 = ee$List(c(0,1,2,3,4,5,6)),
                                      L4 = ee$List(c(0,1,2,3,4,5,6))))

#band names
bandNames = ee$List(c('blue','green','red','nir','swir1','temp','swir2'))
STD_NAMES = c('blue','green','red','nir','swir1','temp','swir2')
bandNumbers = c(0,1,2,3,4,5,6)

# Assumes the image is a Landsat image
maskCloudsAndSuch = function(img){
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
# Basic shadow masking using sum of specified bands
# Tends to include hill shadows and water
shadowThresh = 0.1
shadowSumBands = c('nir','swir1','swir2')
maskShadows = function(img){
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

getComp = function(year,compositingPeriod, startJulian,endJulian){

  #Define dates
  y1Image = year
  y2Image = year + compositingPeriod

  startDate = ee$Date$fromYMD(ee$Number(year),1,1)$advance(startJulian,'day')
  endDate = ee$Date$fromYMD(ee$Number(year)$add(ee$Number(compositingPeriod)),1,1)$advance(endJulian,'day')


  #Helper function to get images from a specified sensor
  getCollection = function(sensor,startDate,endDate,startJulian,endJulian){

    collectionName = collection_dict[sensor][[1]]

    #Start with an un-date-confined collection of iamges
    WOD = ee$ImageCollection(collectionName)$filterBounds(fB)

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

  #Set up our final composite with bands we'd like to include
  composite = composite$select(c('blue','green','red','nir','swir1','swir2'))$multiply(10000)$int16()$clip(fB)

  return(composite)

}

# now that we created the function we can use it!

composite = getComp(2016,4,startJulian,endJulian)

Map$addLayer(composite$clip(compositeArea),
             visParams = list(min = 0, max = 3000, bands = c('red', 'green', 'blue'), gamma = 1.5))

#or using exploreRGEE

composite %>% exploreRGEE::viz(min = 0, max = 3000,
                               band = c('nir', 'green', 'blue'),
                               gamma = 1.5,
                               user_shape = shapes)
