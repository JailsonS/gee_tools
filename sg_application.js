// config
var region = geometry
var startDate = '2024-01-01'
var endDate = '2024-12-30'
var polynomialOrder = 2
var windowSize = 3
var halfWindow = (windowSize - 1)/2


function getFractions(image) {
      // default endmembers
      var ENDMEMBERS = [
          [0.0119,0.0475,0.0169,0.625,0.2399,0.0675], // GV
          [0.1514,0.1597,0.1421,0.3053,0.7707,0.1975], // NPV
          [0.1799,0.2479,0.3158,0.5437,0.7707,0.6646], // Soil
          [0.4031,0.8714,0.79,0.8989,0.7002,0.6607] // Cloud
      ]

      var outBandNames = ['gv', 'npv', 'soil', 'cloud']
      
      var fractions = ee.Image(image)
          .select(['blue', 'green', 'red', 'nir', 'swir1', 'swir2'])
          .unmix(ENDMEMBERS) 
          .max(0);
          //.multiply(100) 
          //.byte() ;
      
      fractions = fractions.rename(outBandNames);
      
      var summed = fractions.expression('b("gv") + b("npv") + b("soil")');
      
      var shade = summed 
          .subtract(1.0) 
          .abs() 
          //.byte() 
          .rename("shade");

      fractions = fractions.addBands(shade);
      
      return image.addBands(fractions);
}


function getNdfi(image){
      var summed = image.expression('b("gv") + b("npv") + b("soil")')
  
      var gvs = image.select("gv").divide(summed).rename("gvs");
  
      var npvSoil = image.expression('b("npv") + b("soil")');
  
      var ndfi = ee.Image.cat(gvs, npvSoil) 
          .normalizedDifference() 
          .rename('ndfi');
  
      // rescale NDFI from 0 to 200 \
      //ndfi = ndfi.expression('byte(b("ndfi") * 100 + 100)');
  
      image = image.addBands(gvs);
      image = image.addBands(ndfi.clamp(-1, 1));
  
      return ee.Image(image);
}


function scaleFactorBands(collection) {

    collection = collection.map(function(image) { 
        var tmStart = image.get('system:time_start');
        var tmEnd = image.get('system:time_end');
        
        
        var optical_bands = image.select('SR_B.').multiply(0.0000275).add(-0.2)
        var thermal_bands = image.select('ST_B.*').multiply(0.00341802).add(149.0)
    
        return image.addBands(optical_bands, null, true).addBands(thermal_bands, null, true).selfMask()
            .copyProperties(image)
            .set('system:time_start', tmStart)
            .set('system:time_end', tmEnd);
    })

    return collection
}


function getCollection(dateStart, dateEnd, cloudCover, roi) {
  
    var collection = null;
    
    var l5 = null;
    var l7 = null;
    var l8 = null;
    var l9 = null;
    
    var bands = [
      'blue',
      'green',
      'red',
      'nir',
      'swir1',
      'swir2',
      'pixel_qa',
      'tir'
    ]
    

    l5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
            .filter('CLOUD_COVER <= ' + cloudCover)
            .filterBounds(roi)
            .filterDate(dateStart, dateEnd)
            .map(function(img){return img.set('time', img.get('system:time_start'))});
            
    
            
    l5 = scaleFactorBands(l5) 
            .select(
              ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7', 'QA_PIXEL', 'ST_B6'], 
              bands
            )
        
    

    
    l7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
            .filter('CLOUD_COVER <= ' + cloudCover)
            .filterBounds(roi)
            .filterDate(dateStart, dateEnd)
            .map(function(img){return img.set('time', img.get('system:time_start'))});
            
    l7 = scaleFactorBands(l7)  
            .select(
              ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7', 'QA_PIXEL', 'ST_B6'], 
              bands
            )

    
    
    
    
    
    l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
            .filter('CLOUD_COVER <= ' + cloudCover)
            .filterBounds(roi)
            .filterDate(dateStart, dateEnd)
            .map(function(img){return img.set('time', img.get('system:time_start'))});
    
    l8 = scaleFactorBands(l8)
            .select(
              ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7', 'QA_PIXEL', 'ST_B10'], 
              bands
            )

    
    
    
    
            
    l9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
            .filter('CLOUD_COVER <= ' + cloudCover)
            .filterBounds(roi)
            .filterDate(dateStart, dateEnd)
            .map(function(img){return img.set('time', img.get('system:time_start'))});
            
    l9 = scaleFactorBands(l9)
            .select(
              ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7', 'QA_PIXEL', 'ST_B10'], 
              bands
            )
            
    
    collection = l5.merge(l7).merge(l8).merge(l9);
    
    
    
    
    
    return collection;

}


function removeCloud(image) {

    var cloudShadowBitMask = 1 << 3;
    var cloudsBitMask = 1 << 4;

    // Get the pixel QA band.
    var qa = image.select('pixel_qa');

    var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
                .and(qa.bitwiseAnd(cloudsBitMask).eq(0));

    return image.updateMask(mask).copyProperties(image);
}

// =================================================================================================


var timeWindow = {

    addTimeRadians: function(image, unit) {
        
        var date = ee.Date(image.get('system:time_start'));
        var years = date.difference(ee.Date('1970-01-01'), unit);
        var timeRadians = ee.Image(years.multiply(2 * Math.PI));
        
        return image.addBands(timeRadians.rename('t').float());
    
    },

    addDateBand: function(image) {
    
        return image.addBands(image.metadata('system:time_start').divide(1e18).rename('time'));
    
    },

    getCollectionFitted: function(collection, winSize, band, unit) {
        
        function smoother(item) {

            function applyFit(img){
                return img.select('time').multiply(fit.select('scale')).add(fit.select('offset'))
                        .set('system:time_start',img.get('system:time_start')).rename('fitted');
            }
            
            var itemDate = ee.List(item).get(0);
            var itemBand = ee.List(item).get(1);
            var itemUnit = ee.List(item).get(2);
            
            var t = ee.Date(itemDate);
            var bandName = itemBand;
            var unit = itemUnit;
            
            var window = dataset.filterDate(t.advance(-windowSize, unit), t.advance(windowSize, unit));
            
            var fit = window.select(['time', bandName]).reduce(ee.Reducer.linearFit());
            
                
            return window.map(applyFit).toList(window.size());
        }
        
        function getMeanWindow(item) {
            var itemDate = ee.List(item).get(0);
            var itemBand = ee.List(item).get(1);
            var itemUnit = ee.List(item).get(2);
            
            var t = ee.Date(itemDate);
            
            return fittedCollection.filterDate(t.advance(-windowSize, itemUnit),t.advance(windowSize, itemUnit))
                .mean().set('system:time_start',t.millis()).rename('fitted');
        }
        
        var dataset = collection.map(this.addDateBand);
        
      
        var dates = ee.List(dataset.aggregate_array('system:time_start'));
        var bandsList = ee.List.repeat(band, dates.size());
        var unitList = ee.List.repeat(unit, dates.size());
        
        
        
        var list = dates.zip(bandsList).zip(unitList).map(function(el) {
            return ee.List(el).flatten();
        });
        
        var windowSize = winSize;
        
        var fittedCollection = ee.ImageCollection(list.map(smoother).flatten());
        
        var smoothedCol = ee.ImageCollection(list.map(getMeanWindow));
        
        return smoothedCol;
    
    }

};



/**
 * 
 * @param {image} image
 * 
 */
function setVariables(image){

    var timestamp = ee.Date(image.get('system:time_start'));

    var ddiff = timestamp.difference(ee.Date(startDate), 'hour');

    var constant = ee.Image(1).toFloat().rename('constant');
    
    image = image.set('date', timestamp)

    var features = image.addBands(constant)

    for(var x=1; x <= polynomialOrder;x++) {
        
        if (x == 1) {
            features = features.addBands(ee.Image(ddiff).toFloat().rename('t' + x.toString()))
        } else {
            features = features.addBands(ee.Image(ddiff).pow(ee.Image(x)).toFloat().rename('t' + x.toString()))
        }
        

    }

    return features
}

/**
 * 
 * @param {imageCollection} collection - collection of indices (ie. 'NDFI')
 * @param {geometry} region 
 * @param {string} startDate 
 * @param {string} endDate 
 * @param {int} polynomialOrder 
 * @param {int} windowSize 
 * 
 */
function applyFilter(collection, region, startDate, endDate, polynomialOrder, windowSize, targetVar) {

    function getCoefFit(i) {
        // Obtém um subconjunto da matriz
        var subarray = array.arraySlice(imageAxis, ee.Number(i).int(), ee.Number(i).add(windowSize).int());
        var predictors = subarray.arraySlice(bandAxis, 2, 2 + polynomialOrder + 1);
        var response = subarray.arraySlice(bandAxis, 0, 1); // Índice de vegetação
        

        // Verifica se a matriz pode ser invertida usando matrixPseudoInverse()
        var coeff = predictors.matrixPseudoInverse().matrixMultiply(response);

        coeff = coeff.arrayProject([0]).arrayFlatten(coeffFlattener);
        
        return coeff;
    }

    // Define the axes of variation in the collection array.
    var imageAxis = 0;
    var bandAxis = 1;

    var coeffFlattener = ['constant']
    var indepSelectors = ['constant']

    for(var x=1; x <= polynomialOrder; x++) {
        coeffFlattener.push('x' + x.toString())
        indepSelectors.push('t' + x.toString())
    }

    coeffFlattener = [coeffFlattener];


    // Add predictors for SG fitting, using date difference
    var temporalCollection = collection
        .filterBounds(region)
        .filterDate(startDate, endDate)
        .select(targetVar)
        .map(setVariables);
    
    

    // Step 3: convert to array type and list
    var array = temporalCollection.toArray();
    var listCollection = temporalCollection.toList(temporalCollection.size());
    

    
    // it process portions of images to smooth 
    var runLength = ee.List.sequence(0, temporalCollection.size().subtract(windowSize));
    
  
    
    
    // Run the SG solver over the series, and return the smoothed image version
    var sgSeries = runLength.map(function(i) {
        var ref = ee.Image(listCollection.get(ee.Number(i).add(halfWindow)));
        var fitted = getCoefFit(i).multiply(ref.select(indepSelectors)).reduce(ee.Reducer.sum())
        return fitted.rename('fitted').copyProperties(ref)
            .set('system:time_start', ref.get('system:time_start'))
            .set('system:time_end', ref.get('system:time_end'))
            .set('system:index', ref.get('system:index'))
    });
    
    return [listCollection, sgSeries]
}


// =================================================================================================



//var landsatColMosaic = ee.ImageCollection('LANDSAT/COMPOSITES/C02/T1_L2_32DAY')
//      .filterDate(startDate, endDate)
//      .filterBounds(region)
//      .map(getFractions)
//      .map(getNdfi)
//      .map(function(img){return img.unmask(-1)})
//      .map(function(img){return img.clip(region)});


var landsatCol = getCollection(startDate, endDate, 100, region)
      .map(removeCloud)
      //.map(function(img){return img.unmask(-1)})
      .map(getFractions)
      .map(getNdfi)
      .select(['ndfi', 'gv'])
      .map(function(img){return img.clip(region)});
      




// apply Time Window Filter
var timeWinCollection = timeWindow.getCollectionFitted(landsatCol, 1, 'ndfi', 'month')




// apply SG filter
var result = applyFilter(landsatCol, region, startDate, endDate, polynomialOrder, windowSize, ['ndfi', 'gv']);
var fitColSg = result[1];
var fitSgCollection = ee.ImageCollection(fitColSg);




// =================================================================================================

Map.setOptions('satellite')

var chart = ui.Chart.image.series({
  imageCollection: landsatCol.select('ndfi').merge(timeWinCollection.select(['fitted'], ['time_win'])),
  region: pt,
  reducer: ee.Reducer.mean(),
  scale: 30, // Ajuste a escala para reduzir o número de pixels processados
  xProperty: 'system:time_start'
}).setOptions({
  title: 'Dados Originais vs. Ajustados',
  hAxis: {title: 'Tempo'},
  vAxis: {title: 'Valor do Índice'},
  interpolateNulls: false,
  lineWidth: 2,
  series: {
    0: {color: 'black', lineDashStyle: [1, 1], pointSize: 2},
    1: {color: 'red', lineDashStyle: [1, 1], pointSize: 2},
    //2: {color: 'blue'}
  }
}).setChartType('LineChart');

print(chart);




