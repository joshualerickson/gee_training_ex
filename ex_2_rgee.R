# excercise 2
library(exploreRGEE)
library(rgee)
ee_Initialize()

# bring in composite from ex_1 and calc NDVI

ndvi <- composite$normalizedDifference(c('nir','red'))$rename('ndvi')

# now add this to composite as band

composite_ndvi <- composite$addBands(ndvi)

#check meta data for band name etc.

meta <- composite_ndvi$getInfo()


#now visualize

composite_ndvi %>% exploreRGEE::viz(min = 0, max = 1,
                                    palette = 'RdYlGn',
                                    user_shape = shapes,
                                    band = 'ndvi')
