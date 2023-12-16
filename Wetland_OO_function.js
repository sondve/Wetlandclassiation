// data转为image函数
function createGLCBands(classValues, classNames, image, year) {
  return ee.Image(classValues.map(function(classValues, index) {
    var band = image.eq(classValues);
    return band.rename(classNames[index]);
  })).set('year', year);
}


var generateGrid = function(xmin, ymin, xmax, ymax, dx, dy, marginx, marginy, opt_proj) {
    var proj = opt_proj || 'EPSG:4326';
    
    dx = ee.Number(dx);
    dy = ee.Number(dy);
  
    var xx = ee.List.sequence(xmin, ee.Number(xmax).subtract(ee.Number(dx).multiply(0.1)), dx);
    var yy = ee.List.sequence(ymin, ee.Number(ymax).subtract(ee.Number(dy).multiply(0.1)), dy);
    
    var cells = xx.map(function(x) {
      return yy.map(function(y) {
        var x1 = ee.Number(x).subtract(marginx);
        var x2 = ee.Number(x).add(ee.Number(dx)).add(marginx);
        var y1 = ee.Number(y).subtract(marginy);
        var y2 = ee.Number(y).add(ee.Number(dy)).add(marginy);
        
        var coords = ee.List([x1, y1, x2, y2]);
        var rect = ee.Algorithms.GeometryConstructors.Rectangle(coords, proj, false);
      
        var nx = x1.add(dx.multiply(0.5)).subtract(xmin).divide(dx).floor();
        var ny = y1.add(dy.multiply(0.5)).subtract(ymin).divide(dy).floor();
      
        return ee.Feature(rect)
          .set({ 
            nx: nx.format('%d'),
            ny: ny.format('%d'),
          });
          // .set({cell_id: x1.format('%.3f').cat('_').cat(y1.format('%.3f')) })
      });
    }).flatten();
  
    return ee.FeatureCollection(cells);
  }; 
