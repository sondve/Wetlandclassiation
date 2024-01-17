function mapping(year){
var roi = AOI.geometry();
Map.addLayer(roi, {'color':'grey'}, 'StudyArea'); 
Map.centerObject(roi,5);

/*** 获取ALOS DEM数据集 ***/ 
var alos_demCol = ee.ImageCollection('JAXA/ALOS/AW3D30/V3_2').select('DSM');
var alos_elevation = alos_demCol.mosaic().clip(roi).unmask(0);
var proj = alos_demCol.first().select(0).projection();
var alos_slope = ee.Terrain.slope(alos_elevation.setDefaultProjection(proj));
var terrains = alos_elevation.addBands(alos_slope);

/*** Remove cloud  ***/
function cloudMask(image){
    var scored = ee.Algorithms.Landsat.simpleCloudScore(image);
    var mask = scored.select(['cloud']).lte(20);
    var qa = image.select('QA_PIXEL');
    var cloud_shadow = qa.bitwiseAnd(1 << 4);
    mask=mask.and(cloud_shadow.not());
    image = image.addBands(mask.rename('CloudMask'));
    image = image.addBands(scored.select('cloud').multiply(-1).add(100).rename('NotCloudScore'));
    return image.updateMask(mask) ;
}
/*** Blue:0.45 - 0.51 μm, Green:0.53 - 0.59 μm, Red:0.64 - 0.67 μm, Near infrared:0.85 - 0.88 μm
     Shortwave infrared 1:1.57 - 1.65 μm;Shortwave infrared 2:2.11 - 2.29 μm ***/
var l8rawName = ee.List(['B2', 'B3', 'B4', 'B5', 'B6', 'B7']);
/*** Blue:0.45 - 0.52 μm, Green:0.52 - 0.60 μm, Red:0.63 - 0.69 μm, Near infrared:0.77 - 0.90 μm
     Shortwave infrared 1:1.55 - 1.75 μm;Shortwave infrared 2:2.08 - 2.35 μm ***/
var l57rawName = ee.List(['B1', 'B2', 'B3', 'B4', 'B5', 'B7']);
var bands_index = ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'NDVI', 'LSWI', 'mNDWI', 'LTideI'];//, 'glcm'];
var Wetland_MaxExtentMask = Wetland_MaxExtent.filterBounds(roi).mosaic().clip(roi);

function Landsat_vis(image) {
  // 计算不同的指数
  var ndwi = image.expression('(bg-bnir)/(bg+bnir)', {'bg': image.select('B3'), 'bnir': image.select('B5')});
  var ndvi = image.expression('(bnir-br)/(br+bnir)', {'br': image.select('B4'), 'bnir': image.select('B5')}).rename('NDVI');
  var lswi = image.expression('(bnir-bswir)/(bswir+bnir)', {'bnir': image.select('B5'), 'bswir': image.select('B6')});
  var evi = image.expression('(b5-b4)*2.5/(b5+6*b4-7.5*b2+10000.0)', {'b2': image.select('B2'), 'b4': image.select('B4'), 'b5': image.select('B5')});
  var mndwi = image.expression('(bg-bswir)/(bg+bswir)', {'bg': image.select('B3'), 'bswir': image.select('B6')});
  var LtideI = image.expression('(bnir-bmax)/(bmax+bnir+1000)', {'bmax': image.select([1,2,5]).reduce(ee.Reducer.max()), 'bnir': image.select('B5')}).rename('LTideI');
  return image.addBands(ndwi.rename('NDWI')).addBands(ndvi.rename('NDVI')).addBands(lswi.rename('LSWI'))
              .addBands(evi.rename('EVI')).addBands(mndwi.rename('mNDWI')).addBands(LtideI.rename('LTideI'));
}
//Create and rescale a grayscale image for GLCM 
//----- Define the GLCM indices used in input for the PCA
function GLCM(image){
  var glcm_bands= ["gray_asm","gray_contrast","gray_corr","gray_ent","gray_var","gray_idm","gray_savg"];
  var gray = image.expression(
        '(0.3 * NIR) + (0.59 * R) + (0.11 * G)', {
        'NIR': image.select('B5'),
        'R': image.select('B4'),
        'G': image.select('B3'),
  }).rename('gray');
  
  // the glcmTexture size (in pixel) can be adjusted considering the spatial resolution and the object textural characteristics
  var glcm = gray.unitScale(0, 0.30).multiply(100).toInt().glcmTexture({size: 2});
  //--- Before the PCA the glcm bands are scaled
  var imageglcm_bands = glcm.select(glcm_bands);
  // calculate the min and max value of an imageglcm_bands
  var minMax = imageglcm_bands.reduceRegion({
    reducer: ee.Reducer.minMax(),
    geometry: roi,
    scale: 30,
    maxPixels: 1e13,
  }); 
  var glcm = ee.ImageCollection.fromImages(
    imageglcm_bands.bandNames().map(function(name){
      name = ee.String(name);
      var band = imageglcm_bands.select(name);
      return band.unitScale(ee.Number(minMax.get(name.cat('_min'))), ee.Number(minMax.get(name.cat('_max'))));
  })).toBands().rename(imageglcm_bands.bandNames());
  
  //---- Apply the PCA
  // The code relating to the PCA was adapted from the GEE documentation  https://developers.google.com/earth-engine/guides/arrays_eigen_analysis
  // Get some information about the input to be used later.
  //var scale = glcm.projection().nominalScale();
  var scale = 30;
  var bandNames = glcm.bandNames();
  // Mean center the data to enable a faster covariance reducer and an SD stretch of the principal components.
  var meanDict = glcm.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: roi, 
      scale: scale,
      maxPixels: 1e13
  });
  var means = ee.Image.constant(meanDict.values(bandNames));
  var centered = glcm.subtract(means);
  
  // This helper function returns a list of new band names.
  var getNewBandNames = function(prefix) {
    var seq = ee.List.sequence(1, bandNames.length());
    return seq.map(function(b) {
      return ee.String(prefix).cat(ee.Number(b).int());
    });
  };
  
  // This function accepts mean centered imagery, a scale and a region in which to perform the analysis. 
  // It returns the Principal Components (PC) in the region as a new image.
  var getPrincipalComponents = function(centered, scale, region) {
    // Collapse the bands of the image into a 1D array per pixel.
    var arrays = centered.toArray();
    // Compute the covariance of the bands within the region.
    var covar = arrays.reduceRegion({
      reducer: ee.Reducer.centeredCovariance(),
      geometry: region,
      scale: scale, 
      maxPixels: 1e13
    });
    // Get the 'array' covariance result and cast to an array.
    // This represents the band-to-band covariance within the region.
    var covarArray = ee.Array(covar.get('array'));
    // Perform an eigen analysis and slice apart the values and vectors.
    var eigens = covarArray.eigen();
    // This is a P-length vector of Eigenvalues.
    var eigenValues = eigens.slice(1, 0, 1);
    // This is a PxP matrix with eigenvectors in rows.
    var eigenVectors = eigens.slice(1, 1);
    // Convert the array image to 2D arrays for matrix computations.
    var arrayImage = arrays.toArray(1);
    // Left multiply the image array by the matrix of eigenvectors.
    var principalComponents = ee.Image(eigenVectors).matrixMultiply(arrayImage);
    // Turn the square roots of the Eigenvalues into a P-band image.
    var sdImage = ee.Image(eigenValues.sqrt())
      .arrayProject([0]).arrayFlatten([getNewBandNames('sd')]);
    // Turn the PCs into a P-band image, normalized by SD.
    return principalComponents
      // Throw out an an unneeded dimension, [[]] -> [].
      .arrayProject([0])
      // Make the one band array image a multi-band image, [] -> image.
      .arrayFlatten([getNewBandNames('pc')])
      // Normalize the PCs by their SDs.
      .divide(sdImage);
  };
  // Get the PCs at the specified scale and in the specified region
  var glcm_pc1 = getPrincipalComponents(centered, scale, roi).select('pc1').rename('glcm');
  return image.addBands(glcm_pc1);
}
// 选择Landsat数据集
var dataset = ee.ImageCollection("LANDSAT/LT05/C02/T1_TOA").filterBounds(roi).map(cloudMask).select(l57rawName, l8rawName)// Landsat5
  .merge(ee.ImageCollection("LANDSAT/LE07/C02/T1_TOA").filterBounds(roi).map(cloudMask).select(l57rawName, l8rawName))// Landsat7
  .merge(ee.ImageCollection("LANDSAT/LC08/C02/T1_TOA").filterBounds(roi).map(cloudMask).select(l8rawName))// Landsat8
  .filterDate('2000-01-01', '2020-12-31').map(Landsat_vis);//.map(GLCM);
//var year = 2020;
var epoch_dataset = dataset.filterDate(year + '-01-01', year + '-12-31');
/*** when year=2000, the real timescale is 2000-2001. ***/
// var epoch_dataset = dataset.filterDate('2000-01-01', '2001-12-31');
var epoch_highTide = epoch_dataset.select(['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'mNDWI', 'LTideI']).qualityMosaic('mNDWI');
var epoch_lowrTide = epoch_dataset.select(['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'mNDWI', 'LTideI']).qualityMosaic('LTideI');
var epoch_percentiles = epoch_dataset.select(bands_index).reduce(ee.Reducer.percentile([10, 30, 50, 70, 90]));
var epoch_textures = epoch_percentiles.select(['B5_p10', 'B5_p30', 'B5_p50', 'B5_p70', 'B5_p90']).multiply(10000).toInt().glcmTexture(3).select([8, 1, 3, 4, 2]);

var train_features = epoch_highTide.addBands(epoch_lowrTide).addBands(epoch_percentiles).addBands(epoch_textures).addBands(terrains).clip(roi);
var local_trainsamples = stable_sample_library1.filterBounds(roi); 
var train_sampleFeatures = train_features.sampleRegions({collection: local_trainsamples, scale: 30, tileScale: 16, });
//****************生成一个随机数构成的列，均匀分布在0-1之间********************//
var withRandom = train_sampleFeatures.randomColumn('random');//样本点随机的排列
var split = 0.8; 
var trainingPartition = withRandom.filter(ee.Filter.lt('random', split));//筛选80%的样本作为训练样本
var testingPartition = withRandom.filter(ee.Filter.gte('random', split));//筛选20%的样本作为测试样本

var Wetland_rfmodel = ee.Classifier.smileRandomForest({numberOfTrees: 100}).train(trainingPartition, 'landcover', train_features.bandNames());
var Wetland_map = train_features.updateMask(Wetland_MaxExtentMask).classify(Wetland_rfmodel).unmask(0).clip(roi);
  return Wetland_map;
}

exports.mapping = mapping;
