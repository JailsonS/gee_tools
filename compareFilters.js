var palettes = require('users/gena/packages:palettes');
var palettesMb = require('users/mapbiomas/modules:Palettes.js');
var ColorRamp = require('users/joaovsiqueira1/packages:ColorRamp.js');


// config
var region = geometry
var startDate = '2009-07-01'
var endDate = '2011-08-30'

var polynomialOrder = 3
var windowSize = 9
var halfWindow = ee.Number(windowSize).divide(2).floor();

// regular filter
var daysInterval = 30;



var years = [
  //1988,1989,1990,
  //1991,1992,1993,1994,1995,1996,1997,
  //1998,1999,2000,
  2001,2002,2003,2004,2005,
  2006,2007,2008,2009,2010,2011,2012,
  2013,2014,2015,2016,2017,2018,2019, 
  2020,2021,2022,2023,2024
];


// =================================================================================================


var assetP1 = 'projects/ee-mapbiomas-imazon/assets/degradation/dam-frequency-c2';
var assetP2 = 'projects/ee-simex/assets/degradation/dam-frequency-c2';


// =================================================================================================

function getNdvi(image) {
    
    var ndvi = image.expression('(nir - red) / (nir + red)', {
      'nir': image.select('nir'),
      'red': image.select('red')
    }).rename('ndvi')
    
    return image.addBands(ndvi)
    
}


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



var twLinearFit = {
  
    // main goal is smooth a time series

    addDateBand: function(image) {
    
        return image.addBands(image.metadata('system:time_start').divide(1e18).rename('time'));
    
    },

    getCollectionFitted: function(collection, winSize, band, unit) {
        
        function smoother(item) {

            function applyFit(img){
                return img.select('time').multiply(fit.select('scale')).add(fit.select('offset'))
                        .set('system:time_start', img.get('system:time_start'))
                        .rename('fitted');
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
                .mean().set('system:time_start',t.millis()).rename(band);
        }
        
        var dataset = collection.map(this.addDateBand);
        
      
        //var dates = ee.List(dataset.aggregate_array('system:time_start'));
        var dates = ee.List.sequence(
            ee.Date(startDate).millis(),
            ee.Date(endDate).millis(),
            ee.Number(daysInterval).multiply(24).multiply(60).multiply(60).multiply(1000)
        );


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

var twLinearInterpolateFit = {

      // main goal is fill gaps by using linear interpolation
      
      addDateBand: function(image) {

          // The time image doesn't have a mask
          var timeImage = image.metadata('system:time_start').toFloat().rename('time');
           
          // We set the mask of the time band to be the same as the first band of the image
          var timeImageMasked = timeImage.updateMask(image.mask().select(0));
          
          return image.addBands(timeImageMasked).toFloat();

      },

      timeWindowFilter: function(collection, daysTimeFilter) {
          
          var dataset = collection.map(this.addDateBand);
          
          var daysMillis = ee.Number(daysTimeFilter).multiply(1000*60*60*24);
          
          // filter definitions
          var maxDiff = ee.Filter.maxDifference({difference: daysMillis,
            leftField: 'system:time_start',
            rightField: 'system:time_start'
          });
          
          var ltEqFilter = ee.Filter.lessThanOrEquals({leftField:'system:time_start',rightField:'system:time_start'});
          var gtEqFilter = ee.Filter.greaterThanOrEquals({leftField: 'system:time_start',rightField: 'system:time_start'});
          
          
          // combine filters
          var maxDiffAndLt = ee.Filter.and(maxDiff, ltEqFilter);
          var maxDiffAndGt = ee.Filter.and(maxDiff, gtEqFilter);
          
          
          // join filters 
          var joinAfter = ee.Join.saveAll({matchesKey: 'after', ordering: 'system:time_start', ascending: false});
          var joinBefore = ee.Join.saveAll({matchesKey: 'before', ordering: 'system:time_start', ascending: false});
          
          
          //apply joins
          var resAfter = joinAfter.apply({'primary': dataset,'secondary': dataset,'condition': maxDiffAndLt});
          var resBefore = joinBefore.apply({'primary': resAfter,'secondary': resAfter,'condition': maxDiffAndGt});
          
          return resBefore;
          
      },
      
      getCollectionFitted: function(collection, daysInterval, daysTimeFilter) {
        

          function getInterpolation(image){
              
              image = ee.Image(image)
            
              var beforeImages = ee.List(image.get('before'));
              var y1 = ee.ImageCollection.fromImages(beforeImages).mosaic();
              
              var afterImages = ee.List(image.get('after'));
              var y2 = ee.ImageCollection.fromImages(afterImages).mosaic();
                
                
                
              // get variables
              var t1 = y1.select('time').rename('t1');
              var t2 = y2.select('time').rename('t2');
              var t = image.metadata('system:time_start').rename('t');   
              
    
              // aplicamos a interpolação linear
              var timeRatio = t.subtract(t1).divide(t2.subtract(t1));  // (t - t1) / (t2 - t1)
              var interpolated = y1.add(y2.subtract(y1).multiply(timeRatio)); // y = y1 + (y2-y1) * ratio
    
              
              return image.unmask(interpolated).copyProperties(image, ['system:time_start'])
          }



          var totalDays = ee.Date(endDate).difference(ee.Date(startDate), 'day');
          var daysToInterpolate = ee.List.sequence(0, totalDays, daysInterval);

          var interpolatedImages = daysToInterpolate.map(function(day) {
            var image = ee.Image().rename('ndfi').set({
              'system:index': ee.Number(day).format('%d'),
              'system:time_start': ee.Date(startDate).advance(day, 'day').millis(),
              'type': 'interpolated'
            })
            return image
          });
          
          var initerpCol = ee.ImageCollection.fromImages(interpolatedImages)

          var collectionToInterpolate = collection.merge(initerpCol);
        

          var result = this.timeWindowFilter(collectionToInterpolate, daysTimeFilter)
              .map(getInterpolation);

          return ee.ImageCollection(result);  
      }



};

var sgFilter = {
  
  
  /**
   * 
   * @param {image} image
   * 
   */
  setVariables: function (image){
    
      image = ee.Image(image)
  
      var timestamp = ee.Date(image.get('system:time_start'));
      var ddiff = timestamp.difference(ee.Date(startDate), 'hour');
      var constant = ee.Image(1).toFloat().rename('constant');
      
      image = image
          .set('system:time_start', image.get('system:time_start'))
          .set('system:time_end', image.get('system:time_end'))
  
      var features = image.addBands(constant)
  
      for(var x=1; x <= polynomialOrder;x++) {
          
          if (x == 1) {
              features = features.addBands(ee.Image(ddiff).toFloat().rename('t'))
          } else {
              features = features.addBands(ee.Image(ddiff).pow(ee.Image(x)).toFloat().rename('t' + x.toString()))
          }
      }
  
      return features
  },
  
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
  applyFilter: function (collection, region, startDate, endDate, polynomialOrder, windowSize, targetVar) {
  
      function getCoefFit(i) {
          // Obtém um subconjunto da matriz
          var subarray = array.arraySlice(imageAxis, ee.Number(i).int(), ee.Number(i).add(windowSize).int());
          //var predictors = subarray.arraySlice(bandAxis, 2, 2 + polynomialOrder + 1);
          var predictors = subarray.arraySlice(bandAxis, 1, 1 + polynomialOrder + 1);
          var response = subarray.arraySlice(bandAxis, 0, 1); // Índice de vegetação
          
  
          // Verifica se a matriz pode ser invertida usando matrixPseudoInverse()
          //var coeff = predictors.matrixPseudoInverse().matrixMultiply(response);
          var coeff = predictors.matrixSolve(response);
  
          coeff = coeff.arrayProject([0]).arrayFlatten(coeffFlattener);
          
          return coeff;
      }
  
      // Define the axes of variation in the collection array.
      var imageAxis = 0;
      var bandAxis = 1;
  
      var coeffFlattener = ['constant']
      var indepSelectors = ['constant']
  
      for(var x=1; x <= polynomialOrder; x++) {
        
          if(x == 1) {
              coeffFlattener.push('x')
              indepSelectors.push('t')
          } else {
              coeffFlattener.push('x' + x.toString())
              indepSelectors.push('t' + x.toString())
          }
        

      }
  
      coeffFlattener = [coeffFlattener];
  
  
      // Add predictors for SG fitting, using date difference
      var temporalCollection = collection.select(targetVar)
          .map(this.setVariables)
          .sort('system:time_start');
      
  
      // Step 3: convert to array type and list
      var array = temporalCollection.toArray();
      var listCollection = temporalCollection.toList(temporalCollection.size());
      
  
      
      // it process portions of images to smooth 
      var runLength = ee.List.sequence(0, temporalCollection.size().subtract(windowSize));
      
    
      
      
      // Run the SG solver over the series, and return the smoothed image version
      var sgSeries = runLength.map(function(i) {
          var ref = ee.Image(listCollection.get(ee.Number(i).add(halfWindow)));
          var fitted = getCoefFit(i).multiply(ref.select(indepSelectors)).reduce(ee.Reducer.sum())
          return fitted.rename(targetVar[0] + '_fitted').copyProperties(ref)
              .set('system:time_start', ref.get('system:time_start'))
              .set('system:time_end', ref.get('system:time_end'))
              .set('system:index', ref.get('system:index'))
      });
      
      return ee.ImageCollection(sgSeries)
  }



}



// =================================================================================================

var dam = ee.ImageCollection(assetP1).merge(ee.ImageCollection(assetP2));
var disturbance = dam.select('freq_dam');

var disturbanceAll = disturbance.filter(ee.Filter.eq('year', 2010));
    disturbanceAll = disturbanceAll.reduce(ee.Reducer.sum());


// =================================================================================================


var landsatCol = ee.ImageCollection('LANDSAT/COMPOSITES/C02/T1_L2_32DAY')
      .filterDate(startDate, endDate)
      .filterBounds(region)
      .map(getFractions)
      .map(getNdfi)
      .map(function(img){return img.clip(region)});


//var landsatCol = getCollection(startDate, endDate, 100, region)
//      .map(removeCloud)
//      .map(getFractions)
//      .map(getNdfi)
//      .map(getNdvi)
//      .map(function(img){return img.clip(region)})
//      //.map(function(image){return image.unmask(0).clip(region)});



// linear fit
var linearCollection = twLinearFit.getCollectionFitted(landsatCol, 3, 'ndfi', 'month')


var regCollection = twLinearInterpolateFit.getCollectionFitted(linearCollection.select('ndfi'), daysInterval, 45)
    .select('ndfi')
    .map(function(image){return image.unmask(1).clip(region)})




// apply sg
var SgCollection = sgFilter.applyFilter(regCollection, region, startDate, endDate, polynomialOrder, windowSize, ['ndfi'])
    .map(function(image){return image.updateMask(image.neq(0))})



// =================================================================================================

var vis = {
    'min': 1,
    'max': 5, // o max ´é em  torno de 40, mas para melhorar a vis, foi setado 12
    'palette': palettes.matplotlib.inferno[7].slice(1),
    'format': 'png'
};

var imgEx1 = linearCollection.filterDate('2010-01-01', '2010-01-30').select('ndfi');
var imgEx2 = SgCollection.filterDate('2010-01-01', '2010-01-30').select('ndfi_fitted');

Map.addLayer(imgEx1, {
  palette:['red','orange', 'yellow', 'green'],
  min:-1,
  max:1
});

Map.addLayer(imgEx2, {
  palette:['red','orange', 'yellow', 'green'],
  min:-1,
  max:1
});
/*
Map.addLayer(imgExReg, {
  palette:['red','orange', 'yellow', 'green'],
  min:-1,
  max:1
});
*/

Map.addLayer(disturbanceAll, vis, 'dam', false);

// =================================================================================================

Map.setOptions('satellite')

// Função para arredondar timestamps para um mesmo dia
function roundTime(image) {
  var timestamp = image.date().format('YYYY-MM-dd'); // Arredonda para o dia
  return image.set('system:time_start', ee.Date(timestamp).millis());
}

// Aplica a função em todas as coleções
var landsatCol = landsatCol.map(roundTime);
var SgCollection = SgCollection.map(roundTime);
var regCollection = regCollection.map(roundTime);
var linearCollection = linearCollection.map(roundTime);

// Recria o gráfico com as datas alinhadas
var chart = ui.Chart.image.series({
  imageCollection: landsatCol.select(['ndfi'], ['Orig'])
    .merge(SgCollection.select(['ndfi_fitted'], ['SG Filter']))
    .merge(regCollection.select(['ndfi'], ['Reg Col']))
    .merge(linearCollection.select(['ndfi'], ['Lin Col'])),
  region: pt,
  reducer: ee.Reducer.mean(),
  scale: 30,
  xProperty: 'system:time_start'
}).setOptions({
  title: 'Série Temporal com Datas Alinhadas',
  interpolateNulls: true, // Para conectar os pontos
  lineWidth: 2,
  series: {
    0: {color: 'black', pointSize: 2},
    1: {color: 'red', pointSize: 2},
    2: {color: 'blue', pointSize: 2},
    3: {color: 'orange', pointSize: 2}
  }
}).setChartType('LineChart');

print(chart)


// =================================================================================================

var visParam =  {
    palette:['red','orange', 'yellow', 'green'],
    min:-1,
    max:1 
} 

var visualizeImage = function(image) {
  return image.visualize(visParam).clip(region).selfMask()
}

var visCollectionRegular = landsatCol.select(['ndfi'])
  .map(visualizeImage)

var visualizeSgFiltered = SgCollection.select(['ndfi_fitted'])
  .map(visualizeImage)

Export.video.toDrive({
  collection: visCollectionRegular,
  description: 'Regular_Time_Series',
  folder: 'earthengine',
  fileNamePrefix: 'regular',
  framesPerSecond: 2,
  dimensions: 800,
  region: region})
  
Export.video.toDrive({
  collection: visualizeSgFiltered,
  description: 'Filtered_Time_Series',
  folder: 'earthengine',
  fileNamePrefix: 'sg_filtered',
  framesPerSecond: 2,
  dimensions: 800,
  region: region})


// estabelecer limite de área de geometria