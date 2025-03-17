
/**
 * 
 * add harmonic model
 * 
*/
exports.hantModel = {
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
exports.timeWindow = {

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
* add gap fill filter with time window
* 
*/

exports.gapFill = {

    _collectionWithTimeBand: function(collection) {
    var colWithTimeBand = collection.map(function(image) {
        // The time image doesn't have a mask
        var timeImage = image.metadata('time_start').rename('timestamp');
        
        // We set the mask of the time band to be the same as the first band of the image
        var timeImageMasked = timeImage.updateMask(image.mask().select(0));
        
        return image.addBands(timeImageMasked).toFloat();
    }); 
    
    return colWithTimeBand;
    },

    timeWindowFilter: function(days, collection) {
        
        var collectionWithTimeBand = this._collectionWithTimeBand(collection);


        var millis = ee.Number(days).multiply(1000*60*60*24);
        
        // -------------------------------------------------------------------------------------------
        
        var maxDiff = ee.Filter.maxDifference({
            difference: millis,
            leftField: 'time_start',
            rightField: 'time_start'
        });
        
        // We need a lessThanOrEquals filter to find all images after a given image
        // This will compare the given image's timestamp against other images' timestamps
        var lessEqFilter = ee.Filter.lessThanOrEquals({
            leftField: 'time_start',
            rightField: 'time_start'
        });

        // We need a greaterThanOrEquals filter to find all images before a given image
        // This will compare the given image's timestamp against other images' timestamps
        var greaterEqFilter = ee.Filter.greaterThanOrEquals({
            leftField: 'time_start',
            rightField: 'time_start'
        });
        
        // apply joins
        // For the first join, we need to match all images that are after the given image.
        // To do this we need to match 2 conditions
        // 1. The resulting images must be within the specified time-window of target image
        // 2. The target image's timestamp must be lesser than the timestamp of resulting images
        // Combine two filters to match both these conditions
        
        var timeWindowFilter = ee.Filter.and(maxDiff, lessEqFilter);
        
        // This join will find all images after, sorted in descending order
        // This will gives us images so that closest is last
        // -------------------------------------------------------------------------------------------
        var join1 = ee.Join.saveAll({
            matchesKey: 'after',
            ordering: 'time_start',
            ascending: false
        });
        
        // apply join
        var joinRes = join1.apply({
            primary: collectionWithTimeBand,
            secondary: collectionWithTimeBand,
            condition: timeWindowFilter
        });
        
        var timeWindowFilter2 = ee.Filter.and(maxDiff, greaterEqFilter);
        
        var join2 = ee.Join.saveAll({
            matchesKey: 'before',
            ordering: 'time_start',
            ascending: false
        });
        
        // apply join
        var joinRes2 = join2.apply({
            primary: joinRes,
            secondary: joinRes,
            condition: timeWindowFilter2
        });
        
        return joinRes2;
    
    },

    _interpolateFilter: function(image) {
        image = ee.Image(image);
        
        // We get the list of before and after images from the image property
        // Mosaic the images so we a before and after image with the closest unmasked pixel
        
        var beforeImages = ee.List(image.get('before'));
        var beforeMosaic = ee.ImageCollection.fromImages(beforeImages).mosaic();
        
        var afterImages = ee.List(image.get('after'));
        var afterMosaic = ee.ImageCollection.fromImages(afterImages).mosaic();
        
        // Interpolation formula
        // y = y1 + (y2-y1)*((t – t1) / (t2 – t1))
        // y = interpolated image
        // y1 = before image
        // y2 = after image
        // t = interpolation timestamp
        // t1 = before image timestamp
        // t2 = after image timestamp
        
        // We first compute the ratio (t – t1) / (t2 – t1)
        
        // Get image with before and after times
        var t1 = beforeMosaic.select('timestamp').rename('t1');
        var t2 = afterMosaic.select('timestamp').rename('t2');
        
        var t = image.metadata('time_start').rename('t');
        
        var timeImage = ee.Image.cat([t1, t2, t]);
        
        var timeRatio = timeImage.expression('(t - t1) / (t2 - t1)', {
        't': timeImage.select('t'),
        't1': timeImage.select('t1'),
        't2': timeImage.select('t2'),
        });
        
        // Compute an image with the interpolated image y
        var interpolated = beforeMosaic
        .add((afterMosaic.subtract(beforeMosaic).multiply(timeRatio)));
        // Replace the masked pixels in the current image with the average value
        var result = image.unmask(interpolated);
        
        return result.copyProperties(image, ['time_start']);
    },

    applyInterpolate: function(collection) {
        return ee.ImageCollection(collection.map(this._interpolateFilter));
    },

  
};














