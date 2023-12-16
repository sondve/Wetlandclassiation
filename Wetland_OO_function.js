// data转为image函数
var createGLCBands = function(classValues, classNames, image, year) {
  return ee.Image(classValues.map(function(classValues, index) {
    var band = image.eq(classValues);
    return band.rename(classNames[index]);
  })).set('year', year);
}

exports.createGLCBands = createGLCBands;
