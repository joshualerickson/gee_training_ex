# /*
#   Day 2, Exercise 4.1: Run a Random Forest Soil Classification in GEE
# NRCS Soil Mapping/Classification in Earth Engine
# Author:
#   Juliette Bateman
# juliette.bateman@usda.gov
# RedCastle Resources
# USFS Geospatial Technology and Applications Center (GTAC)
# Last Modified:
#
#   Summary:
#   This script uses:
#   a cloud free landsat composite for an user specified AOI,
# spectral indices/image enhancements (NDVI, mineral indices, tasseled cap transformation),
# topographic derivates (% slope, aspect cos & sin, elevation),
# climatatic data (seasonal precipitation and temperature).
# to run a random forest soil classification and validation in Essex, VT.
# */


#   /////////////////////////////////////////////////////////////////////////////////////
#   // 1. Prepare Landsat Composite & Predictor Layers
# /////////////////////////////////////////////////////////////////////////////////////
#
#
  # // **** User Editable Variables **** //

# we'll need to read-in Juliette's pedons and boundary
library(exploreRGEE)
library(rgee)
ee_Initialize(gcs = TRUE)

TX_pedons <- ee$FeatureCollection('users/julbateman/NRCS_GEE_DSM/Texas/TX_RyanFlats_pedons')
TX_boundary <- ee$FeatureCollection('users/julbateman/NRCS_GEE_DSM/Texas/TX_RyanFlats_boundary')
tx <- ee_as_sf(TX_boundary)
#check
Map$addLayer(TX_boundary) + Map$addLayer(TX_pedons)

year = 2019 #// Start year for composite
startJulian = 100 #// Starting Julian Date
endJulian = 272 #// Ending Julian date
compositingPeriod = 1 #// Number of years into the future to include
compositeArea = TX_boundary$geometry()$bounds() #// ROI shapefile/asset or polygon
roiName = 'RyanFlats_TX' #// Give the study area a descriptive name.
exportToDrive = 'no' #// Option to export landsat composite to drive
crs = 'EPSG:32613' #// EPSG number for output projection. 32613 = WGS84/UTM Zone 13N. For more info- http://spatialreference.org/ref/epsg/



# Use function to load Landsat composite
new_image = getComp(compositeArea, year, compositingPeriod, startJulian, endJulian, exportToDrive, roiName, crs)
Map$addLayer(new_image, visParams = list(bands = c('red', 'green', 'blue'),
                                         min = 0,
                                         max = 3000))

# // Load in other predictive layers, stack them all
inSpectral = getSpectralIndices(new_image, crs)

Map$addLayer(inSpectral, visParams = list(bands = c('carbonateIndex', 'ferrousIndex',
                                                    'clayIndex'),
                                         min = 0,
                                         max = 3))

inSpectralTopo = getTopo(compositeArea, inSpectral)

Map$addLayer(inSpectralTopo, visParams = list(bands = c('slopePCT'),
                                          min = 0,
                                          max = 2))

inSpectralTopoClimate = getClimate(TX_boundary, year, compositingPeriod, crs)
inSpectralTopoClimate = inSpectralTopoClimate$addBands(inSpectralTopo)
Map$addLayer(inSpectralTopoClimate, visParams = list(bands = c('Summer Temp'),
                                              min = 20,
                                              max = 30))

# // Reproject predictor layers to the same projection and scale
inStack = inSpectralTopoClimate$reproject(crs, NULL, 30)
inStack = inSpectralTopo$reproject(crs, NULL, 30)
# /////////////////////////////////////////////////////////////////////////////////////
#   // 2. Prepare Training/Validation Data
# /////////////////////////////////////////////////////////////////////////////////////
#
#   // Load in training data, separate a 70%/30% for training/validation
data = ee$FeatureCollection(TX_pedons, 'geometry')
datawithColumn = data$randomColumn('random')
split = 0.7 # separate 70% for training, 30% for validation
trainingData = datawithColumn$filter(ee$Filter$lt('random', split))
validationData = datawithColumn$filter(ee$Filter$gte('random', split))

# /////////////////////////////////////////////////////////////////////////////////////
# // 3. Random Forest Classification
# /////////////////////////////////////////////////////////////////////////////////////
#
# // Select predictor layers to include in classification

bands = inStack$bandNames() # //All bands are included here

# // Intersect training points with predictor layers to get training data
training = inStack$select(bands)$sampleRegions(
  collection = trainingData,
  properties = list('Class_num'),
  scale = 30
  )

# // Make RF classifier, train it
classifier = ee$Classifier$smileRandomForest(100, NULL, 1, 0.5, NULL, 0)$setOutputMode('CLASSIFICATION')$train(
    features = training,
    classProperty = 'Class_num',
    inputProperties = bands
    )

# // Classifiy image
classified = inStack$select(bands)$classify(classifier)

#
# /////////////////////////////////////////////////////////////////////////////////////
# // 4. Add classifcation to Map and Create a Legend
# /////////////////////////////////////////////////////////////////////////////////////
#
# // **** Display final classification **** //
# // Make palette for 7 soil types
palette = c(
  '8B0000', #// Barlite 1
  'FFA500', #// Berrend 2
  '4169E1', #// Chilimol 3
  '006400', #// Marfa 4
  '808000', #// Murray 5
  'BC8F8F', #// Musquiz 6
  'FFFF00' #// Phantom 7
)

stats <- classified$reduceRegion(reducer = ee$Reducer$max(),
       geometry = TX_boundary,
       scale = 100
)

print(stats$getInfo())
# // Display the classification

Map$addLayer(classified, visParams = list(min = 1, max = 7, band = 'classification'),
             'classification')

classified %>% exploreRGEE::viz(min = 1, max = 7, user_shape = tx, band = 'classification')

class_ras <- ee_as_raster(
  classified,
  region = TX_boundary$geometry(),
  dsn = 'class_ras',
  via = 'gcs',
  container = 'joshualerickson',
  scale = 30,
  lazy = TRUE
)

# // **** Make a Legend **** //
# // Create the panel for the legend items.



# /////////////////////////////////////////////////////////////////////////////////////
#   // 5. Create Accuracy Assessment Statistics & Figures
# /////////////////////////////////////////////////////////////////////////////////////

  # // Get accuracy assessment statistics, print them
trainAccuracy = classifier$confusionMatrix()

conf_mat <- trainAccuracy$getInfo()
conf_mat <- as.data.frame(matrix(unlist(conf_mat), nrow=length(unlist(conf_mat[1]))))

conf_mat

# // **** Variable Importance Histogram **** //
#   // Get variable importance
dict = classifier$explain()

# //print("Explain:", dict);
variableImportance = ee$Feature(NULL, ee$Dictionary(dict)$get('importance'))

library(tidyverse)
variableImportance <- variableImportance$getInfo()$properties %>% data.frame() %>%
  pivot_longer(everything())

# now plot
library(wesanderson)
variableImportance %>% mutate(name = fct_reorder(name, value)) %>%
  ggplot(aes(name, value, fill = value)) +
  geom_col() + coord_flip() + scale_fill_gradientn(colours = wes_palette('Zissou1', n = 11, type = 'continuous')) +
  labs(y = 'Importance Scores', x = 'Features', title = 'Variable Importance Plot',
       fill = 'Importance Scores', subtitle = 'predicting pedons in southwest TX') + theme_bw()



# // **** Predicted vs Observed Scatterplot **** //

#   // Get predicted classification points in same location as training data
predicted = classified$sampleRegions(collection = trainingData, geometries = TRUE)

# // Separate the predicted (classification) and observed (Class_num) properties
sample = predicted$select(c('classification', 'Class_num'))
sample = sample$getInfo()

predicted = map_df(sample$features, ~.$properties['classification'])
observed = map_df(sample$features, ~.$properties['Class_num'])
results = data.frame(predicted = predicted$classification, observed = observed$Class_num)
# // Create chart, print it

results %>% ggplot(aes(observed, predicted)) +
  geom_point() + geom_line() + theme_bw()

# /////////////////////////////////////////////////////////////////////////////////////
#   // 6. Validation
# /////////////////////////////////////////////////////////////////////////////////////

  # // Sample predicted classification points for validation
validation = classified$sampleRegions(
  collection = validationData,
  properties = list('Class_num'),
  scale = 30
)

# // Compare validation data to classification
testAccuracy = validation$errorMatrix('Class_num', 'classification')

ta <- testAccuracy$getInfo()


conf_mat_acc <- as.data.frame(matrix(unlist(ta), nrow=length(unlist(ta[1]))))
conf_mat_acc


# // **** Predicted vs Observed Scatterplot for validation data!!!!! **** //

#   // Get predicted classification points in same location as training data
predicted = classified$sampleRegions(collection = validationData, geometries = TRUE)

# // Separate the predicted (classification) and observed (Class_num) properties
sample = predicted$select(c('classification', 'Class_num'))
sample = sample$getInfo()

predicted = map_df(sample$features, ~.$properties['classification'])
observed = map_df(sample$features, ~.$properties['Class_num'])
results = data.frame(predicted = predicted$classification, observed = observed$Class_num)
# // Create chart, print it

results %>% ggplot(aes(observed, predicted)) +
  geom_point() + geom_abline() + theme_bw()
