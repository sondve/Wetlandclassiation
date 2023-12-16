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

function maskL8sr(image) {
  var cloudShadowBitMask = 1 << 3;
  var cloudsBitMask = 1 << 5;

  var qa = image.select('pixel_qa');

  var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
      .and(qa.bitwiseAnd(cloudsBitMask).eq(0));

  return image.updateMask(mask).divide(10000)
      .copyProperties(image, ["system:time_start"]); //      .select("B[0-9]*")
}

var mosaicByDate = function(imcol){
  // convert the ImageCollection into List
  var imlist = imcol.toList(imcol.size());
  // print(imlist)
  
  // Obtain the distinct image dates from the ImageCollection
  var unique_dates = imlist.map(function(im){
    return ee.Image(im).date().format("YYYY-MM-dd");
  }).distinct();
  // print(unique_dates);
  
  // mosaic the images acquired on the same date
  var mosaic_imlist = unique_dates.map(function(d){
    d = ee.Date(d);
    var im = imcol.filterDate(d, d.advance(1, "day")).mosaic();
    //print(im)
    // return the mosaiced same-date images and set the time properties
    return im.set(
      "system:time_start", d.millis(), 
      "system:id", d.format("YYYY-MM-dd")
      );
  });
  return ee.ImageCollection(mosaic_imlist);
} 

// 生成月度合成影像
var monthcomposite = function(image){
  var monthcomposite = ee.List.sequence(0, 1*6).map(function(n) {
  var start = ee.Date(year+'-04-01').advance(n, 'month');
  var end = start.advance(1, 'month');
  return image.filterDate(start, end).select("NDVI").median().clip(roi);
  }).flatten();
  // print(monthcomposite);
  // 将月度合成影像合并为图像集合
  return imageCol = ee.ImageCollection.fromImages(monthcomposite);
}


exports.CreateBands = CreateBands;
exports.createLandCoverSample = createLandCoverSample;
exports.maskL8sr = maskL8sr;
exports.mosaicByDate = mosaicByDate;
exports.monthcomposite = monthcomposite;
