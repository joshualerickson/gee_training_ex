

library(exploreRGEE)
library(rgee)
library(mapedit)
ee_Initialize()

# create your AOI
shapes <- mapedit::drawFeatures()

# this creates a GEE geom with an sf object
shapes_boundary <- exploreRGEE:::sf_setup(shapes)$geom


# 2. CALCULATE TOPOGRAPHIC DERIVATIVES, STACK LAYERS

 # IMPORTANT*****aspect layer will have holes where there is zero slope

# Function to add USGS 1/3 arc second topography and derive slope, aspect, & hillshade

getNEDTopography = function(compositeArea){

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

topo = getNEDTopography(shapes_boundary)

  # Display topo derivatives
Map$addLayer(topo$select("elevation"),visParams =  list(min= 200, max= 1200), "Elevation")
Map$addLayer(topo$select("slopePCT"),visParams =  list(min= 15, max= 55), "% Slope")
Map$addLayer(topo$select("sin"),visParams =  list(min= -1, max= 1), "Aspect: sin")
Map$addLayer(topo$select("cos"),visParams =  list(min= -1, max= 1), "Aspect: cos")


# climate section ---------------------------------------------------------

#now work with climate data

# /*
#   Day 1, Exercise 3.2: Calculate Climate Data in GEE
# NRCS Soil Mapping/Classification in Earth Engine
# Summary:
#   This script calculates climate variables for a user-specified AOI using PRISM climate data.
# Climate outputs include: total winter precip, total summer precip,
# mean winter max temp, mean summer max temp,
# and seasonal differences in precip and temp.
# */

#   /////////////////////////////////////////////////////////////////////////////////////
#   // 1. INPUT PREP
# /////////////////////////////////////////////////////////////////////////////////////
#
  # // Specify composite area
# // Specify composite time period - same range will be used for seasonal data
year = 2019 #// Start year for composite
compositingPeriod = 0 #// Number of years into the future to include

# /////////////////////////////////////////////////////////////////////////////////////
#   // 2. CALCULATE CLIMATE DATA
# /////////////////////////////////////////////////////////////////////////////////////
#
#   // Import PRISM daily climate data, ~5km resolution
# // May want to use different climate inputs at small spatial scales

clim = ee$ImageCollection("OREGONSTATE/PRISM/AN81d")
#
# // Specify seasonal date ranges.
# // Seasonal date ranges are based on following dates:
#   //winter: 12/21 - 3/19
# // - Winter is about 89 days long.
# //spring: 3/20 - 6/19
# // -- Spring is 93 days.
# //summer:6/20 - 9/21
# // -- Summer is 94 days.
# //fall: 9/22 - 12/2
# // -- Autumn is 90 days.

summer_startMonth = 6 #//Calendar month 1-12
summer_startDay = 20 #//Day of month
summer_length = 94 #//Number of days of the season
winter_startMonth = 12
winter_startDay = 21
winter_length = 89

# // Function inputs are band, seasonal dates, season length
season_climate = function(band, month_start, day_start, season_length, year, compositingPeriod){

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

# // Calculate multi-year mean of total seasonal precipitation: summer, winter, seasonal difference
winter_precip = season_climate('ppt', winter_startMonth, winter_startDay, winter_length, 2019, 1)$mean()$rename("Winter Precip")
summer_precip = season_climate('ppt', summer_startMonth, summer_startDay, summer_length, 2019, 1)$mean()$rename("Summer Precip")

# // Calculate multi-year mean of mean seasonal temperature: summer, winter, seasonal difference
 winter_tmean = season_climate('tmean', winter_startMonth, winter_startDay, winter_length, 2019, 1)$mean()$rename("Winter Temp")
summer_tmean = season_climate('tmean', summer_startMonth, summer_startDay, summer_length, 2019, 1)$mean()$rename("Summer Temp")

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
climateBands = climateBands$clip(shapes_boundary)
#
# /////////////////////////////////////////////////////////////////////////////////////
#   // 3. DISPLAY
# /////////////////////////////////////////////////////////////////////////////////////
#
  # // Visualize precip layers
Map$addLayer(climateBands$select('Winter Precip'), visParams = list(min=130, max=450), 'winter precip')
Map$addLayer(climateBands$select('Summer Precip'), visParams = list(min=130, max=450), 'summer precip')
Map$addLayer(climateBands$select('Seasonal Precip Difference'), visParams = list(min=50, max=150), 'precip difference')
# // Visualize temp layers
Map$addLayer(climateBands$select('Winter Temp'), visParams = list(min=-10, max=10), 'winter temp')
Map$addLayer(climateBands$select('Summer Temp'), visParams = list(min=15, max=30), 'summer temp')
Map$addLayer(climateBands$select('Seasonal Temp Difference'),visParams = list(min=-10, max=10), 'temp difference')

