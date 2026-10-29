var study_area = 
    /* color: #d63000 */
    /* shown: false */
    ee.Geometry.Polygon(
        [[[-180, 88],
          [0, 88],
          [180, 88],
          [180, -60],
          [0, -60],
          [-180, -60]]], null, false);

// -----------------------------------------------------

var LC_ASSET_ROOT = 'users/wxwang/GlobalD/';
var OUT_ASSET_ROOT = 'projects/ee-naisoild/assets/GlobalD/';
var LAI_COLLECTION = 'users/wxwang/Modis_LAI_1km';

var TARGET_YEAR = 2024;  //  2016-2024
var EXPORT_TAG = 'GRASS';

var BASELINE_START_YEAR = 2000;
var BASELINE_END_YEAR = 2015;
var LP_WINDOW_YEARS = 16;
var LP_RECENT_YEARS = 3;

var SCALE = 1000;
var SOC_THRESHOLD = 0.10;
var PERFORMANCE_THRESHOLD = 0.50;

if (TARGET_YEAR < 2016 || TARGET_YEAR > 2024) {
  throw new Error('TARGET_YEAR must be between 2016 and 2024 for the current LC/LAI asset range.');
}

var fromValues = [1,2,3,4,5,6,7,8,9,10];
var toValues   = [10,20,30,40,50,60,70,80,90,100];

var tempmoist = ee.Image('projects/ee-tufangbobo/assets/Ecosystem/Temp_Moist');
var ecoRegions = ee.FeatureCollection('RESOLVE/ECOREGIONS/2017').filterBounds(study_area);


// ----------------- 1. Common functions -----------------

function getLC(year) {
  return ee.Image(LC_ASSET_ROOT + 'lc_' + String(year) + '_1km')
    .remap(fromValues, toValues)
    .rename('landcover')
    .toInt16();
}

var lc2000_ = getLC(BASELINE_START_YEAR);
var lc2015_ = getLC(BASELINE_END_YEAR);
var landMask = lc2015_.neq(60).selfMask();

function asClass(img) {
  return img.rename('class').toInt8();
}

function applyStatusMatrix(baselineAssessment, periodAssessment, bandName, g) {
  var b = baselineAssessment.rename('baseline');
  var p = periodAssessment.rename('period');

  return p.expression(
    'p == 3 ? 3 : (p == 2 ? 2 : b)',
    {'b': b, 'p': p}
  ).rename(bandName).clip(g).toInt8();
}

function characterizeStatus(baselineAssessment, periodAssessment, year, g) {
  var b = baselineAssessment.rename('baseline');
  var p = periodAssessment.rename('period');

  return p.expression(
    'b == 3 && p == 3 ? 1 : ' +
    '((b == 1 || b == 2) && p == 3 ? 2 : ' +
    '(b == 3 && p == 1 ? 3 : ' +
    '(b == 1 && p == 1 ? 4 : ' +
    '(b == 2 && p == 1 ? 5 : ' +
    '((b == 3 || b == 1) && p == 2 ? 6 : 7)))))',
    {'b': b, 'p': p}
  ).rename('SDG1531_StatusPath_' + String(year)).clip(g).toInt8();
}


// ----------------------- 2. LC -------------------------

function classifyLCTransition(lcStartRaw, lcEndRaw, g) {
  var lcStart = lcStartRaw.where(lcStartRaw.eq(40), ee.Image(30)).clip(g);
  var lcEnd = lcEndRaw.where(lcEndRaw.eq(40), ee.Image(30)).clip(g);

  var change = lcStart.neq(lcEnd).selfMask().toInt();
  var stable = asClass(lcStart.eq(lcEnd).selfMask().toInt());

  var forest = asClass(change.updateMask(lcStart.eq(20))
    .expression('end != 20 ? 3 : 1', {'end': lcEnd})
    .updateMask(change.gt(0).and(lcStart.eq(20))));

  var grassland = asClass(change.updateMask(lcStart.eq(30))
    .expression('end == 20 || end == 50 || end == 10 ? 2 : 3', {'end': lcEnd})
    .updateMask(change.gt(0).and(lcStart.eq(30))));

  var cropland = asClass(change.updateMask(lcStart.eq(10))
    .expression('end == 20 || end == 50 ? 2 : 3', {'end': lcEnd})
    .updateMask(change.gt(0).and(lcStart.eq(10))));

  var wetlands = asClass(change.updateMask(lcStart.eq(50))
    .expression('end != 50 ? 3 : 1', {'end': lcEnd})
    .updateMask(change.gt(0).and(lcStart.eq(50))));

  var settlements = asClass(change.updateMask(lcStart.eq(80))
    .expression('end != 80 ? 2 : 1', {'end': lcEnd})
    .updateMask(change.gt(0).and(lcStart.eq(80))));

  var others = asClass(change.updateMask(
      lcStart.eq(60).or(lcStart.eq(70)).or(lcStart.eq(90)).or(lcStart.eq(100))
    )
    .expression('end == 80 ? 3 : 2', {'end': lcEnd})
    .updateMask(change.gt(0).and(
      lcStart.eq(60).or(lcStart.eq(70)).or(lcStart.eq(90)).or(lcStart.eq(100))
    )));

  return ee.ImageCollection([
    stable, forest, grassland, cropland, wetlands, settlements, others
  ]).mosaic().rename('class').clip(g).toInt8();
}

function lcBaselineAssessment(g) {
  return classifyLCTransition(lc2000_, lc2015_, g)
    .rename('LC_baseline').updateMask(landMask).toInt8();
}

function lcPeriodAssessment(year, g) {
  return classifyLCTransition(lc2015_, getLC(year), g)
    .rename('LC_period_' + String(year)).updateMask(landMask).toInt8();
}


// --------------------- 3. SOC -----------------------

function toUNCCD7(lc) {
  return ee.Image(0)
    .where(lc.eq(20), 1)
    .where(lc.eq(30).or(lc.eq(40)).or(lc.eq(70)), 2)
    .where(lc.eq(10), 3)
    .where(lc.eq(50), 4)
    .where(lc.eq(80), 5)
    .where(lc.eq(90).or(lc.eq(100)), 6)
    .where(lc.eq(60), 7)
    .rename('unccd_lc').toInt8();
}

function croplandSOCFactor() {
  return ee.Image(0.58)
    .where(tempmoist.eq(1).or(tempmoist.eq(2))
      .or(tempmoist.eq(10)).or(tempmoist.eq(11))
      .or(tempmoist.eq(13)).or(tempmoist.eq(14))
      .or(tempmoist.eq(16)).or(tempmoist.eq(17)), 0.80)
    .where(tempmoist.eq(3).or(tempmoist.eq(9))
      .or(tempmoist.eq(12)).or(tempmoist.eq(15)), 0.69)
    .where(tempmoist.eq(7), 0.48)
    .where(tempmoist.eq(8), 0.64)
    .rename('f');
}

function socConversionFactor(lcStartRaw, lcEndRaw, g) {
  var s = toUNCCD7(lcStartRaw).clip(g);
  var e = toUNCCD7(lcEndRaw).clip(g);
  var f = croplandSOCFactor().clip(g);

  return s.expression(
    '(s == 1 || s == 2) ? ' +
      '(e == 3 ? f : ((e == 5 || e == 6) ? 0.1 : 1)) : ' +
    '(s == 3 ? ' +
      '((e == 1 || e == 2) ? 1 / f : ' +
       '(e == 4 ? 1 / 0.71 : ((e == 5 || e == 6) ? 0.1 : 1))) : ' +
    '(s == 4 ? ' +
      '(e == 3 ? 0.71 : ((e == 5 || e == 6) ? 0.1 : 1)) : ' +
    '((s == 5 || s == 6) ? ' +
      '((e == 1 || e == 2 || e == 3 || e == 4) ? 2 : 1) : 1)))',
    {'s': s, 'e': e, 'f': f}
  ).rename('conversion_factor').clip(g);
}

function yearsSinceFirstLCChange(lcStartRaw, startYear, endYear, g) {
  var changeYears = [];

  for (var y = startYear + 1; y <= endYear; y++) {
    changeYears.push(
      ee.Image.constant(y)
        .toInt16()
        .updateMask(getLC(y).neq(lcStartRaw))
        .rename('change_year')
        .clip(g)
    );
  }

  var firstYear = ee.ImageCollection.fromImages(changeYears)
    .min()
    .rename('change_year')
    .toInt16()
    .clip(g);

  return ee.Image.constant(endYear)
    .toInt16()
    .subtract(firstYear)
    .add(1)
    .unmask(0)
    .max(0)
    .rename('years_since_change')
    .toInt16()
    .clip(g);
}

function classifySOCTransition(lcStartRaw, lcEndRaw, startYear, endYear, g) {
  var valid = toUNCCD7(lcStartRaw).gt(0)
    .and(toUNCCD7(lcEndRaw).gt(0))
    .and(lcStartRaw.neq(60))
    .selfMask().clip(g);
  var conversionFactor = socConversionFactor(lcStartRaw, lcEndRaw, g);

  var yearsSinceChange = yearsSinceFirstLCChange(
    lcStartRaw, startYear, endYear, g
  );

  var rate = conversionFactor.subtract(1)
    .multiply(yearsSinceChange.divide(20))
    .updateMask(valid)
    .clip(g);

  return rate.expression(
    'r < -thr ? 3 : (r > thr ? 2 : 1)',
    {
      'r': rate,
      'thr': ee.Number(SOC_THRESHOLD)
    }
  ).rename('class').clip(g).updateMask(valid).toInt8();
}

function socBaselineAssessment(g) {
  return classifySOCTransition(
    lc2000_, lc2015_, BASELINE_START_YEAR, BASELINE_END_YEAR, g
  )
    .rename('SOC_baseline').updateMask(landMask).toInt8();
}

function socPeriodAssessment(year, g) {
  return classifySOCTransition(
    lc2015_, getLC(year), BASELINE_END_YEAR, year, g
  )
    .rename('SOC_period_' + String(year)).updateMask(landMask).toInt8();
}


// ---------------------- 4. LP-LAI -------------------------

var prodCol = ee.ImageCollection(LAI_COLLECTION)
  .map(function(img) {
    var year = ee.Number(ee.Algorithms.If(
      img.propertyNames().contains('year'),
      img.get('year'),
      ee.Date(img.get('system:time_start')).get('year')
    ));

    return img.select([0])
      .rename('prod')
      .toFloat()
      .set('year', year)
      .set('system:time_start', ee.Date.fromYMD(year, 1, 1).millis());
  })
  .filter(ee.Filter.gte('year', BASELINE_START_YEAR))
  .filter(ee.Filter.lte('year', 2024));

function filterProd(startYear, endYear) {
  return prodCol
    .filter(ee.Filter.gte('year', startYear))
    .filter(ee.Filter.lte('year', endYear));
}

function mannKendallTheilSen(imageCollection, bandName) {
  var ordered = imageCollection.sort('year');
  var n = ee.Number(ordered.size());
  var list = ordered.toList(n);

  var s = ee.ImageCollection.fromImages(
    ee.List.sequence(0, n.subtract(2)).map(function(i) {
      i = ee.Number(i);
      var earlier = ee.Image(list.get(i)).select(bandName);
      return ee.ImageCollection.fromImages(
        ee.List.sequence(i.add(1), n.subtract(1)).map(function(j) {
          var later = ee.Image(list.get(ee.Number(j))).select(bandName);
          var difference = later.subtract(earlier);
          return difference.gt(0)
            .subtract(difference.lt(0))
            .rename('S');
        })
      ).sum();
    })
  ).sum().rename('S');

  var varS = n.multiply(n.subtract(1)).multiply(n.multiply(2).add(5)).divide(18);
  var z = s.expression(
    'S > 0 ? (S - 1) / sd : (S < 0 ? (S + 1) / sd : 0)',
    {'S': s, 'sd': varS.sqrt()}
  ).rename('z');

  var xy = ordered.map(function(img) {
    var x = ee.Image.constant(ee.Number(img.get('year'))).rename('year').toFloat();
    return x.addBands(img.select(bandName).rename('prod').toFloat());
  });
  var sen = xy.select(['year', 'prod']).reduce(ee.Reducer.sensSlope());

  return sen.select('slope').rename('slope').addBands(z);
}

function performanceReference(windowMean, referenceLC, g) {
  var classes = [10, 20, 30, 40, 50, 70, 90, 100];
  var refs = [];

  classes.forEach(function(code) {
    var classMean = windowMean.updateMask(referenceLC.eq(code));
    var stats = classMean.reduceRegions({
      collection: ecoRegions.filterBounds(g),
      reducer: ee.Reducer.percentile([90], ['p90']),
      scale: SCALE,
      tileScale: 4
    });
    refs.push(stats.reduceToImage(['p90'], ee.Reducer.first())
      .updateMask(referenceLC.eq(code)).rename('p90'));
  });

  return ee.ImageCollection.fromImages(refs).mosaic().rename('p90').clip(g);
}

function classifyLPWindow(endYear, g) {
  var startYear = endYear - LP_WINDOW_YEARS + 1;
  var recentStartYear = endYear - LP_RECENT_YEARS + 1;
  var historyEndYear = recentStartYear - 1;

  var trendCol = filterProd(startYear, endYear);
  var trendOut = mannKendallTheilSen(trendCol, 'prod').clip(g);

  var output_trend = trendOut.expression(
    's > 0 && z >= 1.96 ? 2 : (s < 0 && z <= -1.96 ? 3 : 1)',
    {
      's': trendOut.select('slope'),
      'z': trendOut.select('z')
    }
  ).rename('trend').clip(g).updateMask(landMask).toInt8();

  var stateBase = filterProd(startYear, historyEndYear);
  var stateTarget = filterProd(recentStartYear, endYear);

  var base_mean = stateBase.mean().rename('prod').clip(g);
  var base_std = stateBase.reduce(ee.Reducer.stdDev()).rename('prod').clip(g);
  var target_mean = stateTarget.mean().rename('prod').clip(g);

  var z_score = target_mean.subtract(base_mean)
    .divide(base_std.divide(ee.Number(LP_RECENT_YEARS).sqrt()))
    .updateMask(base_std.gt(0))
    .clip(g);

  var output_state = z_score.expression(
    'z < -1.96 ? 3 : (z > 1.96 ? 2 : 1)',
    {'z': z_score}
  ).rename('state').clip(g).updateMask(landMask).toInt8();

  var window_mean = trendCol.mean().rename('prod').clip(g);
  var target_perf = performanceReference(window_mean, lc2015_, g);
  var target_ratio = window_mean.divide(target_perf)
    .updateMask(target_perf.gt(0))
    .rename('ratio');

  var output_performance = target_ratio.expression(
    'r < threshold ? 3 : 1',
    {'r': target_ratio, 'threshold': ee.Number(PERFORMANCE_THRESHOLD)}
  ).rename('performance').clip(g).updateMask(landMask).toInt8();

  return output_trend.expression(
    't == 3 ? 3 : ' +
    '(t == 2 ? (s == 3 && p == 3 ? 3 : 2) : ' +
    '(s == 3 ? 3 : 1))',
    {
      't': output_trend,
      's': output_state,
      'p': output_performance
    }
  ).rename('class').clip(g).updateMask(landMask).toInt8();
}

function lpBaselineAssessment(g) {
  return classifyLPWindow(BASELINE_END_YEAR, g)
    .rename('LP_baseline').toInt8();
}

function lpPeriodAssessment(year, g) {
  var lcTarget = getLC(year);
  var lpMask = lcTarget.lte(50).or(lcTarget.eq(70)).or(lcTarget.eq(90)).clip(g);
  return classifyLPWindow(year, g)
    .rename('LP_period_' + String(year)).updateMask(lpMask).toInt8();
}


// ----------------- 5. 1OAO, Status matrix and characterization -----------------

function combine1OAO(lcImg, socImg, lpImg, bandName, g) {
  var lc = lcImg.rename('class');
  var soc = socImg.rename('class');
  var lp = lpImg.rename('class');

  var degraded = lc.eq(3).unmask(0)
    .or(soc.eq(3).unmask(0))
    .or(lp.eq(3).unmask(0));

  var improved = lc.eq(2).unmask(0)
    .or(soc.eq(2).unmask(0))
    .or(lp.eq(2).unmask(0));

  var valid = lc.mask().gt(0).unmask(0)
    .or(soc.mask().gt(0).unmask(0))
    .or(lp.mask().gt(0).unmask(0));

  return ee.Image(1)
    .where(improved, 2)
    .where(degraded, 3)
    .rename(bandName)
    .clip(g)
    .updateMask(valid)
    .toInt8();
}

var lcBaseline = lcBaselineAssessment(study_area);
var socBaseline = socBaselineAssessment(study_area);
var lpBaseline = lpBaselineAssessment(study_area);

var lcPeriod = lcPeriodAssessment(TARGET_YEAR, study_area);
var socPeriod = socPeriodAssessment(TARGET_YEAR, study_area);
var lpPeriod = lpPeriodAssessment(TARGET_YEAR, study_area);

var baselineAssessment = combine1OAO(
  lcBaseline, socBaseline, lpBaseline,
  'SDG1531_BaselineAssessment_2015', study_area
);

var periodAssessment = combine1OAO(
  lcPeriod, socPeriod, lpPeriod,
  'SDG1531_PeriodAssessment_' + String(TARGET_YEAR), study_area
);

var finalImg = applyStatusMatrix(
  baselineAssessment, periodAssessment,
  'SDG1531_Status_' + String(TARGET_YEAR), study_area
);

var statusPath = characterizeStatus(
  baselineAssessment, periodAssessment, TARGET_YEAR, study_area
);

var counterbalance = ee.Image(0)
  .where(statusPath.eq(1).or(statusPath.eq(2)), 1)
  .where(statusPath.eq(6).or(statusPath.eq(7)), 2)
  .rename('LDN_Counterbalance_' + String(TARGET_YEAR))
  .updateMask(finalImg.mask())
  .toInt8();

var diagnosticStack = ee.Image.cat([
  lcBaseline, socBaseline, lpBaseline,
  lcPeriod, socPeriod, lpPeriod,
  baselineAssessment, periodAssessment,
  finalImg, statusPath, counterbalance,
  lc2015_.rename('land_type_2015')
]);


// ----------------- 6. export -----------------

// Export.image.toAsset({
//   image: lcPeriod,
//   description: 'SDG1531_LC_period_' + String(TARGET_YEAR),
//   assetId: OUT_ASSET_ROOT + EXPORT_TAG + '_LC_period_' + String(TARGET_YEAR),
//   region: study_area,
//   scale: SCALE,
//   maxPixels: 1e13,
//   pyramidingPolicy: {'.default': 'mode'}
// });

// Export.image.toAsset({
//   image: socPeriod,
//   description: 'SDG1531_SOC_period_' + String(TARGET_YEAR),
//   assetId: OUT_ASSET_ROOT + EXPORT_TAG + '_SOC_period_' + String(TARGET_YEAR),
//   region: study_area,
//   scale: SCALE,
//   maxPixels: 1e13,
//   pyramidingPolicy: {'.default': 'mode'}
// });

// Export.image.toAsset({
//   image: lpPeriod,
//   description: 'SDG1531_LP_LAI_period_' + String(TARGET_YEAR),
//   assetId: OUT_ASSET_ROOT + EXPORT_TAG + '_LP_LAI_period_' + String(TARGET_YEAR),
//   region: study_area,
//   scale: SCALE,
//   maxPixels: 1e13,
//   pyramidingPolicy: {'.default': 'mode'}
// });



// Export.image.toAsset({
//   image: finalImg,
//   description: 'SDG1531_Status_' + String(TARGET_YEAR),
//   assetId: OUT_ASSET_ROOT + EXPORT_TAG + '_Status_' + String(TARGET_YEAR),
//   region: study_area,
//   scale: SCALE,
//   maxPixels: 1e13,
//   pyramidingPolicy: {'.default': 'mode'}
// });

// Export.image.toAsset({
//   image: statusPath,
//   description: 'SDG1531_StatusPath_' + String(TARGET_YEAR),
//   assetId: OUT_ASSET_ROOT + EXPORT_TAG + '_StatusPath_' + String(TARGET_YEAR),
//   region: study_area,
//   scale: SCALE,
//   maxPixels: 1e13,
//   pyramidingPolicy: {'.default': 'mode'}
// });

Export.image.toAsset({
  image: diagnosticStack,
  description: 'SDG1531_Diagnostics_' + String(TARGET_YEAR),
  assetId: OUT_ASSET_ROOT + EXPORT_TAG + '_Diagnostics_' + String(TARGET_YEAR),
  region: study_area,
  scale: SCALE,
  maxPixels: 1e13,
  pyramidingPolicy: {'.default': 'mode'}
});

print('Baseline assessment', baselineAssessment);
print('Period assessment', periodAssessment);
print('Final status', finalImg);
print('Seven-class status pathway', statusPath);
