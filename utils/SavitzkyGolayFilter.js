
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


// Solve 
function getLocalFit(i) {
    // Get a slice corresponding to the window_size of the SG smoother
    var subarray = array.arraySlice(imageAxis, ee.Number(i).int(), ee.Number(i).add(window_size).int())
    var predictors = subarray.arraySlice(bandAxis, 2, 2 + order + 1)
    var response = subarray.arraySlice(bandAxis, 0, 1); // vegetation indice
    var coeff = predictors.matrixSolve(response)
  
    coeff = coeff.arrayProject([0]).arrayFlatten(coeffFlattener)
    return coeff  
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
exports.applyFilter = function(collection, region, startDate, endDate, polynomialOrder, windowSize) {

    function getLocalFit(i) {
        // Get a slice corresponding to the window_size of the SG smoother
        var subarray = array.arraySlice(imageAxis, ee.Number(i).int(), ee.Number(i).add(windowSize).int())
        var predictors = subarray.arraySlice(bandAxis, 2, 2 + polynomialOrder + 1)
        var response = subarray.arraySlice(bandAxis, 0, 1); // vegetation indice
        var coeff = predictors.matrixSolve(response)
            coeff = coeff.arrayProject([0]).arrayFlatten(coeffFlattener);
        
        return coeff  
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


    // Step 1 - Add predictors for SG fitting, using date difference
    var temporalCollection = collection
        .filterBounds(region)
        .filterDate(startDate, endDate)
        .map(setVariables);
    
    

    // Step 2: Set up Savitzky-Golay smoothing
    var halfWindow = (windowSize - 1)/2


    // Step 3: convert to array type and list
    var array = temporalCollection.toArray();
    var listCollection = temporalCollection.toList(temporalCollection.size());
    
    // it process portions of images to smooth 
    var runLength = ee.List.sequence(0, temporalCollection.size().subtract(windowSize));
    

    // Run the SG solver over the series, and return the smoothed image version
    var sgSeries = runLength.map(function(i) {
        var ref = ee.Image(listCollection.get(ee.Number(i).add(halfWindow)));
        var fitted = getLocalFit(i).multiply(ref.select(indepSelectors)).reduce(ee.Reducer.sum())
        return fitted.rename('fitted').copyProperties(ref)
            .set('system:time_start', ref.get('system:time_start'))
            .set('system:time_end', ref.get('system:time_end'))
            .set('system:index', ref.get('system:index'))
    });
    
    return [temporalCollection, sgSeries]
}