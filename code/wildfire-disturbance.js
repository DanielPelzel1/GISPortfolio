/*  This script applies a disturbance index developed by Mildrexler, et al., (2007 - Ecological Applications)
to the Landsat 5,7,8, 9 time series to identify deviations of year X against normal values derived from 
preceding years imagery.

The disturbance index relies on a ratio between yearly maximum NDVI and yearly maximum surface temperature
                      (LST / NDVI)

I have converted all simple ratios to normalized difference ratios

(LST/NDVI) for year N is compared to the average (LST/NDVI) of preceding years imagery.
The original publication used a simple ration of year N (LST/NDVI) and the average (LST/NDVI)

(LST/NDVI) for year N / average LST/NDVI

I have modified this by generating a normalized difference between year N and average

(YearN - Ave) / (YearN + Ave)


*/

//Specify the year to map  
//  The start year is the Year - 10 so the selected year is compared to the previous 10 years.
var Year = 2020
var startYear = Year - 20
var Years = ee.List.sequence(startYear, Year);

// Identify a region to map
// var fc =ee.Geometry.Rectangle([-117.129, 32.552, -106.626, 43.138]); // Utah, Central Strip
// Load a Fusion Table of state boundaries and filter.
var fc = Tiger.filter(ee.Filter.or(
         ee.Filter.eq('NAME', 'Utah'),ee.Filter.eq('NAME', 'Nevada')))

//  This sets up the selection and renaming of bands so all images will have the same spectral bands
var l89bands = ['SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7','ST_B10'];
var l57bands = ['SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','SR_B7','ST_B6'];
var renameTo = ['B','G','R','NIR','SWIR1','SWIR2','Thermal'];

//**********************************************************************************************
//************* FUNCTIONS **********************************************************************
//**********************************************************************************************
// This example demonstrates the use of the Landsat 8 Collection 2, Level 2
// QA_PIXEL band (CFMask) to mask unwanted pixels.

function maskClouds(image) {
  // Bit 0 - Fill
  // Bit 1 - Dilated Cloud
  // Bit 2 - Cirrus
  // Bit 3 - Cloud
  // Bit 4 - Cloud Shadow
  var qaMask = image.select('QA_PIXEL').bitwiseAnd(parseInt('11111', 2)).eq(0);
  var saturationMask = image.select('QA_RADSAT').eq(0);

  // Apply the scaling factors to the appropriate bands.
  var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);
  var thermalBands = image.select('ST_B.*').multiply(0.00341802).add(149.0);

  // Replace the original bands with the scaled ones and apply the masks.
  return image.addBands(opticalBands, null, true)
      .addBands(thermalBands, null, true)
      .updateMask(qaMask)
      .updateMask(saturationMask);
}

// NDVI Functions for L8 and L7 (the L7 function is used for L4&5)
function addDist(image) {
  var ndvi = image.normalizedDifference(['NIR', 'R']).add(1).multiply(25).rename('NDVI'); //since I added 1 to the NDVI, then multiply by 25 to put it in the same numerical range as LST
  // var ndviMsk = ndvi.updateMask(ndvi.gte(0)).unmask(0)  //weird tings happen with negative NDVI's so turn them all to 0.  the line above corrects for this
  var lst = image.select('Thermal').subtract(273.15).rename('LST'); // Temperature is currently in Kelvin, convert to C.
  var dist = ((lst.subtract(ndvi)).divide(lst.add(ndvi))).multiply(50).rename('DIST')
  return image.addBands([ndvi,lst,dist]);
}
function renameL89(image){
  return image.select(l89bands).rename(renameTo) 
}
function renameL57(image){
  return image.select(l57bands).rename(renameTo)
}
//**********************************************************************************************
//**********************************************************************************************
//*****************Get the Image Collections and merge into one********************************
// Load the Landsat 9 ImageCollection.  Mask clouds and add NDVI.
var l9 = L9SR.filter(ee.Filter.calendarRange(startYear, Year, 'year'))
            .filter(ee.Filter.calendarRange(6, 8, 'month'))
            .filterBounds(fc).map(maskClouds).map(renameL89).map(addDist)

// Load the Landsat 8 ImageCollection.  Mask clouds and add NDVI.
var l8 = L8SR.filter(ee.Filter.calendarRange(startYear, Year, 'year'))
            .filter(ee.Filter.calendarRange(6, 8, 'month'))
            .filterBounds(fc).map(maskClouds).map(renameL89).map(addDist)

// Load the Landsat 7 ImageCollection.  Mask clouds and add NDVI.  MINIMIZE use of L7
var l7 = L7SR.filterDate('2012-03', '2013-05')
            .filter(ee.Filter.calendarRange(6, 8, 'month'))
            .filterBounds(fc).map(maskClouds).map(renameL57).map(addDist)

// Load the Landsat 5 ImageCollection.  Mask clouds and add NDVI.
var l5 = L5SR.filter(ee.Filter.calendarRange(startYear, Year, 'year'))
            .filter(ee.Filter.calendarRange(6, 8, 'month'))
            .filterBounds(fc).map(maskClouds).map(renameL57).map(addDist)

var Merged = l5.merge(l7).merge(l8).merge(l9).sort("system:time_start").select('NDVI','LST','DIST')
// print('Merged',Merged.size())

//**********************************************************************************************
//**********************************************************************************************
//***************** CALCULATE YEARLY AVERAGES **************************************************
// Group images by year and reduce to a single image
// the result is an ImageCollection with one image for each year.

//  MODIFY TO SEPARATELY SELECT MAX LST AND MAX NDVI FOR EACH YEAR.  I BELIEVE THE .MAX() FUNCTION
//   DOES THIS, BUT NEED TO VALIDATE.
var byYear = ee.ImageCollection.fromImages(
      Years.map(function (Y) {
        return Merged.filter(ee.Filter.calendarRange(Y, Y, 'year'))
                    .max()  // take the max values of NDVI, LST, and DIST
                    .set('year', Y);
}));
print('byYear',byYear)

// Select all images prior to selected year and reduce those to the mean.
var average= byYear.filterMetadata('year','less_than',Year).mean().clip(fc);
print('Previous years',byYear.filterMetadata('year','less_than',Year))

// Extract from the NDVI and LST yearly image composites - the year that I want to compare against the norm of earlier years
//  Convert the image collection to an image (there should only be one image in the collection.
var theYear = byYear.filterMetadata('year','equals',Year).first().clip(fc);
print('theYear',theYear)

// Check on the min/max of the surface temperature ('LST') band to see if it makes sense
var minMax = theYear.reduceRegion({
    reducer: ee.Reducer.minMax(),
    geometry: fc,
    scale: 90,
    maxPixels: 1e15
})
print('Min/Max',minMax)

//  Calculate the normalized difference of the NDVI and LST between the selected year and the previous years average.
var disturbance = (((theYear.select('LST').divide(theYear.select('NDVI'))).subtract((average.select('LST').divide(average.select('NDVI')))))
                  .divide(((theYear.select('LST').divide(theYear.select('NDVI'))).add((average.select('LST').divide(average.select('NDVI'))))))).rename('Disturbance');
print('disturbance',disturbance)

//  The visualization parameters ramps most negative values as green to most positive values as red
var vizParams = {min:-0.3,
                max:0.3,
                palette: ['0f6e28','4dc96e','e0e6d1','d16969','ff0000']}

//**********************************************************************************************
//***************** MAKE A PRETTY MAP SECTION **************************************************
//**********************************************************************************************
var HillShade = ee.Terrain.hillshade(ee.Image('USGS/3DEP/10m').clip(fc).multiply(4),135, 45).divide(255);

var disturbanceSHD = disturbance.visualize(vizParams).multiply(HillShade);

//  ADD WATER
var NLCD = ee.Image('USGS/NLCD/NLCD2016').select('landcover').clip(fc)
var Water = NLCD.updateMask(NLCD.eq(11)); // mask it
//  ADD County Boundaries
var Counties = ee.FeatureCollection('TIGER/2018/Counties')
    .filter(ee.Filter.or(
        ee.Filter.eq('STATEFP', '49'),ee.Filter.eq('STATEFP', '32')));
var boundary = ee.Image().toByte()
    .paint(Counties, 3, 1); // Outline using color 3, width 5.

//Blend the Shaded NDVI with the water and then state boundaries.
var shaded_water = disturbanceSHD.visualize({gain:1}).blend(Water.visualize({palette: ['789BE2']})).blend(boundary)

// Map.centerObject(fc,7)
Map.addLayer(disturbance,vizParams,'disturbance',false);
Map.addLayer(shaded_water,{},'Shaded');
Map.addLayer(theYear,{bands: 'LST'}, 'TheYear LST',false)
Map.addLayer(theYear,{bands: 'NDVI'}, 'TheYear NDVI',false)

//**********************************************************************************************
//**********************************************************************************************
//***************************** CREATE QUERY PANEL *********************************************
// Create a panel to hold our widgets.
var panel = ui.Panel();
panel.style().set('width', '500px');

// Create an intro panel with labels.
var intro = ui.Panel([
  ui.Label({
    value: 'Disturbance Index Inspector',
    style: {fontSize: '20px', fontWeight: 'bold'}
  }),
  ui.Label('Click a point on the map to inspect.')
]);
panel.add(intro);

// Create panels to hold lon/lat values.
var lon = ui.Label();
var lat = ui.Label();
panel.add(ui.Panel([lon, lat], ui.Panel.Layout.flow('horizontal')));

// Register a callback on the default map to be invoked when the map is clicked.
Map.onClick(function(coords) {
  // Update the lon/lat panel with values from the click event.
  lon.setValue('lon: ' + coords.lon.toFixed(2)),
  lat.setValue('lat: ' + coords.lat.toFixed(2));

  // Add a red dot for the point clicked on.
  var point = ee.Geometry.Point(coords.lon, coords.lat);
  var dot = ui.Map.Layer(point, {color: 'FFFF00'});
  Map.layers().set(5, dot);  // redo this to query the number of images on the screen and then add 1

  // Create an NDVI temporal chart.
  var distChart = ui.Chart.image.series(byYear.select(['DIST','LST','NDVI']), point, ee.Reducer.mean(), 30,'year'); 
  distChart.setOptions({
    title: 'Disturbance Index Over Time\n'+ 'lon:' + coords.lon.toFixed(2) + ' lat:' + coords.lat.toFixed(2),
    pointSize: 4,
    vAxis: {title: 'LST/NDVI'},
    hAxis: {title: 'Date', gridlines: {count: 12}},
  });
  panel.widgets().set(20, distChart);

});

Map.style().set('cursor', 'crosshair');

// Add the panel to the ui.root.
ui.root.insert(0, panel);
