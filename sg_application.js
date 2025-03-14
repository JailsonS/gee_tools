// config
var region = geometry
var startDate = '2021-01-01'
var endDate = '2024-12-30'
var polynomialOrder = 3
var windowSize = 7
var halfWindow = (windowSize - 1)/2

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


// =================================================================================================



var landsatCol = ee.ImageCollection('LANDSAT/COMPOSITES/C02/T1_L2_32DAY')
      .filterDate(startDate, endDate)
      .filterBounds(region)
      .map(getFractions)
      .map(getNdfi)
      .map(function(img){return img.unmask(-1)})
      .map(function(img){return img.clip(region)});


var result = applyFilter(landsatCol, region, startDate, endDate, polynomialOrder, windowSize, ['ndfi', 'gv']);


var srcCol = result[0];
var fitCol = result[1];

// convert data

var srcCollection = ee.ImageCollection(srcCol);
var fitCollection = ee.ImageCollection(fitCol);


var srcImage = srcCol.slice(1).iterate(
    function(current, previous){
        return ee.Image(current).addBands(ee.Image(previous))
    },
    srcCol.get(0)
);

var fitImage = fitCol.iterate(
    function(current, previous){
        return ee.Image(current).addBands(ee.Image(previous))
    },
    fitCol.get(0)
);




// =================================================================================================



var chart = ui.Chart.image.series({
  imageCollection: srcCollection.select('ndfi').merge(fitCollection.select('fitted')),
  region: pt,
  reducer: ee.Reducer.mean(),
  scale: 30, // Ajuste a escala para reduzir o número de pixels processados
  xProperty: 'system:time_start'
}).setOptions({
  title: 'Dados Originais vs. Ajustados',
  hAxis: {title: 'Tempo'},
  vAxis: {title: 'Valor do Índice'},
  lineWidth: 2,
  series: {
    0: {color: 'blue', lineDashStyle: [4, 4]},
    1: {color: 'red'}
  }
}).setChartType('LineChart');

print(chart);




