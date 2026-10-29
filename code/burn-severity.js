////////////////////////////////////////////////////////
//  WILDFIRE BURN SEVERITY  //
//////////////////////////////////////////////////////// 



/////////////////////////
// USER DEFINED FIRES //
/////////////////////////
var fire1_ID = 'UT3771811283720170617';
var fire2_ID = 'UT4060611070920200821';

// Import dNBR assets from Part 1 and rename band to 'dNBR'
var fire1 = ee.Image('projects/danielpelzel1/assets/UT3771811283720170617_dNBR_Asset').rename('dNBR');
var fire2 = ee.Image('projects/danielpelzel1/assets/UT4060611070920200821_dNBR_Asset').rename('dNBR');

/////////////////////////////
// FUNCTION: CATEGORICAL BURN SEVERITY //
/////////////////////////////
function addCategorical(image) {
  var dNBR_Cat = image.expression(
      "(b('dNBR') > 659) ? 4" +       // High severity
      ": (b('dNBR') > 269) ? 3" +     // Moderate severity
      ": (b('dNBR') > 100) ? 2" +     // Low severity
      ": (b('dNBR') > -100) ? 1" +    // Unburned
      ": 0"                           // Enhanced Regrowth
  ).rename('dNBR_Cat')
   .updateMask(image.abs().add(1));
  
  return image.addBands(dNBR_Cat);
}

// Apply to both fires
fire1 = addCategorical(fire1);
fire2 = addCategorical(fire2);

/////////////////////////////
// VISUALIZATION SETTINGS //
/////////////////////////////
var dNBR_vis = {min: -100, max: 660, palette: ['2AEE20','2AB123','C9E029','FDC30A','D84E17']};
var cat_vis = {min: 0, max: 4, palette: ['2AEE20','2AB123','C9E029','FDC30A','D84E17']};

// Map both fires
Map.centerObject(fire1, 9);
Map.addLayer(fire1.select('dNBR'), dNBR_vis, 'Fire 1 Continuous');
Map.addLayer(fire2.select('dNBR'), dNBR_vis, 'Fire 2 Continuous');

Map.addLayer(fire1.select('dNBR_Cat'), cat_vis, 'Fire 1 Categorical');
Map.addLayer(fire2.select('dNBR_Cat'), cat_vis, 'Fire 2 Categorical');

/////////////////////////////
// FIRE STATISTICS //
/////////////////////////////
var firePerimeters = ee.FeatureCollection("projects/ee-dougramsey/assets/MTBS/MTBS_Perimeters_v2024");

function fireStats(image, fireID) {
  var perimeter = firePerimeters.filter(ee.Filter.eq('Event_ID', fireID)).first();
  
  var mean_val = image.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: perimeter.geometry(),
    scale: 30
  });
  
  var hist = image.reduceRegion({
    reducer: ee.Reducer.histogram(),
    geometry: perimeter.geometry(),
    scale: 30
  });
  
  var hist_chart = ui.Chart.image.histogram({
    image: image.select('dNBR'),
    region: perimeter.geometry(),
    scale: 30,
    minBucketWidth: 10
  });
  
  print('Fire', fireID, 'Mean dNBR:', mean_val);
  print('Fire', fireID, 'Histogram:', hist);
  print('Fire', fireID, 'Histogram Chart:', hist_chart);
}

// Compute stats for both fires
fireStats(fire1, fire1_ID);
fireStats(fire2, fire2_ID);

/////////////////////////////
// DISTANCE TO NON-HIGH-SEVERITY //
/////////////////////////////
function addDistance(image) {
  var dNBR_HighSev = image.expression(
      "(b('dNBR_Cat') == 4) ? 0 : 1"
  ).rename('dNBR_HighSev');
  
  image = image.addBands(dNBR_HighSev);
  
  var dNBR_Distance = image.select('dNBR_HighSev')
                            .distance({
                              kernel: ee.Kernel.euclidean({radius:2000, units:'meters'}), 
                              skipMasked:true
                            })
                            .rename('dNBR_Distance')
                            .updateMask(image.select('dNBR').abs().add(1))
                            .round();
  
  return image.addBands(dNBR_Distance);
}

// Apply distance maps
fire1 = addDistance(fire1);
fire2 = addDistance(fire2);

/////////////////////////////
// READY: FIRE LAYERS AND STATS //
/////////////////////////////
print('Fire 1 final image with categorical + distance:', fire1);
print('Fire 2 final image with categorical + distance:', fire2);
