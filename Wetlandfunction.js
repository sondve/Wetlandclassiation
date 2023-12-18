// data转为image函数
var CreateBands = function(classValues, classNames, image, year) {
  return ee.Image(classValues.map(function(classValues, index) {
    var band = image.eq(classValues);
    return band.rename(classNames[index]);
  })).set('year', year);
}

var createLandCoverSample = function(img, roi, classValue, pointNum) {
  // 为图像设置 "system:footprint" 属性
  var imgNeibor = ee.Image(img).set("system:footprint", roi);
  // 对图像进行分层采样
  var pointSample = imgNeibor.selfMask().stratifiedSample({
    numPoints: pointNum,
    region: roi,
    scale: 30,
    seed: 0,
    geometries: true
  });
  // 设置 landcover 属性
  var labeledSample = pointSample.map(function(fea) {
    return fea.set('landcover', classValue);
  });
  return labeledSample;
}

// 生成月度合成影像函数
var generateMonthlyComposite = function(image, start_date, roi) {
  var monthcomposite = ee.List.sequence(0, 1 * 6).map(function(n) {
    var start = ee.Date(start_date).advance(n, 'month');
    var end = start.advance(1, 'month');
    return image.filterDate(start, end).select("NDVI").median().clip(roi);
  }).flatten();
  return monthcomposite;
};

// 从随机设置的样本点中采样数据
var reclassSamplefromSample = function(classSample, image) {
  var ReclassSample = image.sampleRegions({
  collection: classSample,
  properties: ['landcover'],
  scale: 30,
  tileScale: 16,
  geometries: true
});
 return ReclassSample;
}


// 定义函数计算样本点每个波段的均值和标准差
var calculateMeanAndStdDev = function(sample, image) {
  var imgBand = image.bandNames();
  var mean = sample.reduceColumns({
    reducer: ee.Reducer.mean().repeat(7),
    selectors: imgBand
  });
  mean = ee.List(ee.Dictionary(mean).get('mean'));

  var stdDev = sample.reduceColumns({
    reducer: ee.Reducer.stdDev().repeat(7),
    selectors: imgBand
  });
  stdDev = ee.List(ee.Dictionary(stdDev).get('stdDev'));

  var meanStdDevPlus = ee.List([0,1,2,3,4,5,6]).map(function(item){
    var MeanStdDevPlus = ee.Number(mean.get(item)).add(stdDev.get(item));
    return MeanStdDevPlus;
  });
  
  var meanStdDevMinus = ee.List([0,1,2,3,4,5,6]).map(function(item){
  var MeanStdDevMinus = ee.Number(mean.get(item)).subtract(stdDev.get(item));
  return MeanStdDevMinus;
  });

  return {
    mean: mean,
    stdDev: stdDev,
    meanStdDevPlus: meanStdDevPlus,
    meanStdDevMinus: meanStdDevMinus
  }
}

// 添加 flag 属性，用于筛选 Grassland_Sample_flat 样本
var  selectsample = function(sample1, sample_stats, image){
  var imgBand = image.bandNames();
  var Sample2 = sample1.map(function(fea){
  var feaValue = ee.List(ee.Feature(fea).toDictionary(imgBand).values());
  var flag = ee.List([0,1,2,3,4,5,6]).map(function(item){
    var a_b = ee.Number(sample_stats.meanStdDevMinus.get(item)).subtract(feaValue.get(item));
    var c_b = ee.Number(sample_stats.meanStdDevPlus.get(item)).subtract(feaValue.get(item));
    return ee.Number(a_b).multiply(c_b);
  });
  return fea.set('flag0', flag.get(0)).set('flag1', flag.get(1)).set('flag2', flag.get(2))
            .set('flag3', flag.get(3)).set('flag4', flag.get(4)).set('flag5', flag.get(5))
            .set('flag6', flag.get(6));
  });
  // 筛选出 flag 为负值的样本点，即 sample 样本本点
  var Sample3 = Sample2.filter(ee.Filter.lt("flag0", 0))
                            .filter(ee.Filter.lt("flag1", 0))
                            .filter(ee.Filter.lt("flag2", 0))
                            .filter(ee.Filter.lt("flag3", 0))
                            .filter(ee.Filter.lt("flag4", 0))
                            .filter(ee.Filter.lt("flag5", 0))
                            .filter(ee.Filter.lt("flag6", 0));
  return Sample3;
}


exports.CreateBands = CreateBands;
exports.createLandCoverSample = createLandCoverSample;
exports.generateMonthlyComposite = generateMonthlyComposite;
exports.reclassSamplefromSample = reclassSamplefromSample;
exports.calculateMeanAndStdDev = calculateMeanAndStdDev;
exports.selectsample = selectsample;
