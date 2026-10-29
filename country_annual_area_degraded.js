// ---------------- 1. User settings ----------------

var COUNTRIES_ASSET =
  'projects/ee-naisoild/assets/global_countries_china_verified';

var STATUS_ASSET_ROOT = 'projects/ee-naisoild/assets/GlobalD/GRASS_Diagnostics_';

var START_YEAR = 2016;
var END_YEAR = 2024;
var SCALE = 1000;

var countries = ee.FeatureCollection(COUNTRIES_ASSET);

var SAFE_DATA_BOUNDS = ee.Geometry.Rectangle(
  [-180, -60, 180, 88.0], 'EPSG:4326', false
);
var safeCountries = countries
  .filterBounds(SAFE_DATA_BOUNDS)
  .map(function(feature) {
  var safeGeometry = feature.geometry().intersection(
    SAFE_DATA_BOUNDS, ee.ErrorMargin(100)
  );
  return feature.setGeometry(safeGeometry);
  });


// ---------------- 2. Calculate one year ----------------

///////////// 3 types of status
// function summarizeYear(year) {
//   var status = ee.Image(STATUS_ASSET_ROOT + String(year))
//     // .select('SDG1531_Status_'+String(year))
//     .select('SDG1531_StatusPath_'+String(year))
//     .rename('status')
//     .toInt8()
//     .clip(SAFE_DATA_BOUNDS);

//   var pixelAreaKm2 = ee.Image.pixelArea().divide(1e6);
//   var areaBands = ee.Image.cat([
//     pixelAreaKm2.updateMask(status.eq(2))
//       .rename('improved_area_km2'),
//     pixelAreaKm2.updateMask(status.eq(3))
//       .rename('degraded_area_km2'),
//     pixelAreaKm2.updateMask(status.eq(1))
//       .rename('stable_area_km2')
//   ]).toDouble();

//   var result = areaBands.reduceRegions({
//     collection: safeCountries,
//     reducer: ee.Reducer.sum(),
//     scale: SCALE,
//     crs: status.projection(),
//     tileScale: 4,
//     maxPixelsPerRegion: 1e9
//   });

//   return result.map(function(feature) {
//     var improved = ee.Algorithms.If(
//       feature.get('improved_area_km2'),
//       feature.get('improved_area_km2'), 0
//     );
//     var degraded = ee.Algorithms.If(
//       feature.get('degraded_area_km2'),
//       feature.get('degraded_area_km2'), 0
//     );
//     var stable = ee.Algorithms.If(
//       feature.get('stable_area_km2'),
//       feature.get('stable_area_km2'), 0
//     );

//     return feature
//       .set('year', year)
//       .set('improved_area_km2', improved)
//       .set('degraded_area_km2', degraded)
//       .set('stable_area_km2', stable)
//       .setGeometry(null);
//   });
// }

///////////// 7 types of status (pathway)
var TYPE_DEFS = [
  {
    code: 1,
    field: 'persistent_degradation_area_km2'
  },
  {
    code: 2,
    field: 'recent_degradation_area_km2'
  },
  {
    code: 3,
    field: 'baseline_degradation_area_km2'
  },
  {
    code: 4,
    field: 'stability_area_km2'
  },
  {
    code: 5,
    field: 'baseline_improvement_area_km2'
  },
  {
    code: 6,
    field: 'recent_improvement_area_km2'
  },
  {
    code: 7,
    field: 'persistent_improvement_area_km2'
  }
];
function summarizeYear(year) {
  var statusPath = ee.Image(
    STATUS_ASSET_ROOT + String(year)
  )
    .select('SDG1531_StatusPath_'+String(year))
    .rename('status_path')
    .toInt8()
    .clip(SAFE_DATA_BOUNDS);

  var pixelAreaKm2 = ee.Image.pixelArea().divide(1e6);

  var areaImages = TYPE_DEFS.map(function(typeDef) {
    return pixelAreaKm2
      .updateMask(statusPath.eq(typeDef.code))
      .rename(typeDef.field);
  });

  var areaBands = ee.Image.cat(areaImages).toDouble();

  var result = areaBands.reduceRegions({
    collection: safeCountries,
    reducer: ee.Reducer.sum(),
    scale: SCALE,
    crs: statusPath.projection(),
    tileScale: 4,
    maxPixelsPerRegion: 1e9
  });

  return result.map(function(feature) {
    var outputProperties = {
      year: year
    };

    TYPE_DEFS.forEach(function(typeDef) {
      var value = feature.get(typeDef.field);

      outputProperties[typeDef.field] = ee.Algorithms.If(
        ee.Algorithms.IsEqual(value, null),
        0,
        value
      );
    });
    return feature
      .set(outputProperties)
      .setGeometry(null);
  });
}


// ---------------- 3. Combine all years ----------------
var annualTables = [];
for (var year = START_YEAR; year <= END_YEAR; year++) {
  annualTables.push(summarizeYear(year));
}

var allYears = ee.FeatureCollection(annualTables).flatten();

print('Country-year result sample', allYears.limit(20));
print('Number of country-year records', allYears.size());


// ---------------- 4. Export one CSV ----------------

// Export.table.toDrive({
//   collection: allYears,
//   description: 'country_SDG1531_area_2016_2024',
//   fileNamePrefix: 'country_SDG1531_area_2016_2024',
//   fileFormat: 'CSV'
// });
Export.table.toDrive({
  collection: allYears,
  description: 'country_SDG1531_status_paths_2016_2024',
  fileNamePrefix: 'country_SDG1531_status_paths_2016_2024',
  fileFormat: 'CSV'
});