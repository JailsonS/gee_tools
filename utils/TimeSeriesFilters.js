
/**
 * 
 * add harmonic model
 * 
*/
exports.HantModel = {
    //NDVIt = β0 + β1t + β2cos(2πωt) + β3sin(2πωt) + …
    
    // independents = [t, constant]
    getfitValues: function (n_cycles, dependent, independents, timeDif, collection) {
        
        var harmonicTrend = this.getCof(n_cycles, dependent, independents, timeDif, collection);
        
        var harmonicTrendCoefficients = harmonicTrend[1].select('coefficients')
            .arrayProject([0])
            .arrayFlatten([ harmonicTrend[2] ]);
            
        function fitValues(image){
            return image.addBands(image.select(harmonicTrend[2])
                .multiply(harmonicTrendCoefficients)
                .reduce('sum')
                .rename('fitted'));
        }
        
        var fittedHarmonic = harmonicTrend[0].map(fitValues);
        
        return [fittedHarmonic, harmonicTrendCoefficients];
            
    },
    
    getCof: function (n_cycles, dependent, independents, timeDif, collection) {
        
        var harmonicFrequencies = ee.List.sequence(1, n_cycles);
        
        function bandNames(base, list) { 
            return ee.List(list).map( function(i) {
                return ee.String(base).cat(ee.Number(i).int());
            });
        }
        
        var cosName = bandNames('cos_', harmonicFrequencies);
        var sinName = bandNames('sin_', harmonicFrequencies);
        
        independents = cosName.cat(sinName).cat(independents);
        
        function addDependents (image) {
            
            var years = image.date().difference(timeDif, 'year');
            var timeRadians = ee.Image(years.multiply(2 * Math.PI)).rename('t');
            var constant = ee.Image(1);
            
            return image.addBands(constant).addBands(timeRadians.float());
            
        }
        
        function addHarmonics (freqs) {
            return function (image) {
                
                var frequence = ee.Image.constant(freqs);
                var time = ee.Image(image).select('t');
                
                var sin = time.multiply(frequence).sin().rename(sinName);
                var cos = time.multiply(frequence).cos().rename(cosName);
                
                return image.addBands(cos).addBands(sin);
            };
        }
        
        var harmonicCollection = collection.map(addDependents).map( addHarmonics(harmonicFrequencies) );
        
        var harmonicTrend = harmonicCollection
            .select(independents.add(dependent))
            .reduce(ee.Reducer.linearRegression(independents.length(), 1));
        
        return [harmonicCollection, harmonicTrend, independents];
    } 

};




/**
* 
* add time window filter 
* 
*/
exports.LinearFit = {
  
    // main goal is smooth a time series

    addDateBand: function(image) {
    
        return image.addBands(image.metadata('system:time_start').divide(1e18).rename('time'));
    
    },

    getCollectionFitted: function(collection, startDate, endDate, winSize, band, unit) {
        
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
                .mean().set('system:time_start',t.millis()).rename('fitted');
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
        
        var smoothedCol = ee.ImageCollection(list.map(getMeanWindow)).map(
            function(image) {
                return image.clamp(-1, 1).set('system:time_start', image.get('system:time_start'))
            }
        );
        
        return smoothedCol;
    
    }
};


/**
* 
* add linear interpolation
* 
*/

exports.LinearInterpolation = {

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
    
    getCollectionFitted: function(collection, startDate, endDate, daysInterval, daysTimeFilter) {
      

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



/**
* 
* SG Filter
* 
*/

exports.SGFilter =  {

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


    applyFilter: function (collection, polynomialOrder, windowSize, targetVar) {

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









