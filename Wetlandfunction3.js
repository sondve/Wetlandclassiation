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

exports.CreateBands = CreateBands;
exports.createLandCoverSample = createLandCoverSample;
exports.reclassSamplefromSample = reclassSamplefromSample;
