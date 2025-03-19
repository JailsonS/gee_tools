var palettes = require('users/gena/packages:palettes');
var damCore = require('users/engsoaresfilho/temporal_eval:dam_code');

// parâmetros para o DAM
var defaultParams = {
    'tresh_dam_min': -0.250, // limiar para degradação 
    'tresh_dam_max': -0.095, // limiar para degradação 
    'tresh_df_min': -0.250, // limiar para desmatamento
    'time_window': 3
}

var listParams = [
    [2024, defaultParams]
]

// ======================================================================================================


var assetP1 = 'projects/ee-mapbiomas-imazon/assets/degradation/dam-frequency-c2';
var assetP2 = 'projects/ee-simex/assets/degradation/dam-frequency-c2';
var assetLulc = 'projects/mapbiomas-raisg/public/collection3/mapbiomas_raisg_panamazonia_collection3_integration_v2';


// ======================================================================================================

var damOrg = ee.ImageCollection(assetP1).merge(ee.ImageCollection(assetP2));

// ======================================================================================================

var interaface = {
  
    layers:{
    
      before:  ui.Map({ 
        style: {
          border: '2px solid black',
          width: '50%' 
        }
      }),
      
      after:  ui.Map({ 
        style: {
          border: '2px solid black',
          width: '50%' 
        }
      }),
      
    },
    
    widgets: {
        
        runButton: ui.Button({
            "label": "Executar",
            "onClick": function (button) {
                var disabled = button.getDisabled();
                if (!disabled) {
                    interaface.run();
                }
            },
            "style": {
                'padding': '1px',
                'stretch': 'horizontal',
                'position': 'top-left'
            }
        }),
        
        exportButton: ui.Button({
            "label": "Exportar",
            "onClick": function (button) {
                var disabled = button.getDisabled();
                if (!disabled) {
                    interaface.exportImage();
                }
            },
            "style": {
                'padding': '1px',
                'stretch': 'horizontal',
                'position': 'top-left'
            }
        }),
        
        init: function() {
          
          var dt = interaface.layers.before.drawingTools();
              dt.clear();
          
          var defaultGeom = ee.Geometry.Polygon([
              [
                [
                  -55.027232999772195,
                  -12.347199463422708
                ],
                [
                  -53.97941195485032,
                  -12.347199463422708
                ],
                [
                  -53.97941195485032,
                  -11.676936280491203
                ],
                [
                  -55.027232999772195,
                  -11.676936280491203
                ],
                [
                  -55.027232999772195,
                  -12.347199463422708
                ]
              ]
          ])
          
          dt.addLayer([defaultGeom], 'region', 'black', false);
          
          interaface.layers.before.add(this.runButton);
          interaface.layers.before.add(this.exportButton);
          
        }
        
    },
    
    run: function() {
      
        this.clearLayers();
        
        var dt = interaface.layers.before.drawingTools();
        var roi = dt.layers().get(0).toGeometry();
        
        this.clearGeometries();
        
        dt.addLayer([roi], 'region', 'black', false);
        
        interaface.layers.before.add(interaface.widgets.runButton);
        interaface.layers.before.add(interaface.widgets.exportButton);
        
        var mapLeft = interaface.layers.before;
        var mapRight = interaface.layers.after;
        
        mapLeft.setOptions('satellite')
        mapRight.setOptions('satellite')
        
        listParams.forEach(function(params){
          
            var year = params[0];
          
            var band = 'classification_' + (year - 1).toString()
            
            if (year >= 2020) { band = 'classification_2019' }
            
            var lulc = ee.Image(assetLulc).select(band)
          
            var roi = dt.layers().get(0).toGeometry();
            
            
            var styled = ee.FeatureCollection([ee.Feature(roi)]).style({
              color: 'red',        // Cor da borda
              width: 2,            // Espessura da borda
              fillColor: '00000000' // Transparente no interior
            });
              
            var dam = damCore.getDAM(params, roi).updateMask(lulc.eq(3));
            
            var damReg = dam.select('freq_dam')
            var damDefReg = dam.select('freq_dam_df').gt(0)
            
            var damDefOrig = damOrg.select('freq_dam_df').filter(ee.Filter.eq('year', params[0])).sum().gt(0);
            var damOrig = damOrg.select('freq_dam').filter(ee.Filter.eq('year', params[0])).sum();
            
            
            
            mapLeft.centerObject(roi, 11)
          
            mapLeft.addLayer(damReg, {
             bands:['freq_dam'],
              min:1, max:10,
              palette:palettes.cmocean.Thermal[7]
            }, 'DEGRADAÇÃO - REG');
            
            mapLeft.addLayer(damDefReg.selfMask(), {
              min:0, max:1,
              bands:['freq_dam_df'],
              palette:['red']
            }, 'DESMATAMENTO - REG'); 
            
            
            mapRight.addLayer(damOrig, {
              min:1, max:10,
              palette:palettes.cmocean.Thermal[7]
            }, 'DEGRADAÇÃO');
            
            mapRight.addLayer(damDefOrig, {
              min:0, max:1,
              palette: ['red']
            }, 'DESMATAMENTO');       
      
      
            mapLeft.addLayer(styled, {}, 'ROI');
            mapRight.addLayer(styled, {}, 'ROI'); 


            
        });
        
    },
    
    exportImage: function() {
      
        var layers = interaface.layers.before.layers();
        var region = interaface.layers.before.drawingTools().layers().get(0).toGeometry();
        
        var deforestationDAM = layers.get(0).getEeObject();
        var degradationnDAM = layers.get(1).getEeObject();
        
        
        Export.image.toDrive({
          image: deforestationDAM, 
          description: 'deforestationDAM', 
          folder: 'dam', 
          fileNamePrefix: 'deforestationDAM', 
          region:region , 
          scale: 30, 
          maxPixels: 1e13
        });
        
        
        Export.image.toDrive({
          image: degradationnDAM, 
          description: 'degradationnDAM', 
          folder: 'dam', 
          fileNamePrefix: 'degradationnDAM', 
          region:region , 
          scale: 30, 
          maxPixels: 1e13
        });
      
    },
    
    clearLayers: function(){
        
        interaface.layers.before.clear();
        interaface.layers.after.clear();
      
    },
    
    clearGeometries: function() {
        
      var dt = interaface.layers.before.drawingTools();
          dt.clear();
    
    },
    
    initializeLayers: function() {
        var dt = interaface.layers.before.drawingTools();
        var roi = dt.layers().get(0).toGeometry();
        
        this.clearLayers();

        
        interaface.layers.before.add(interaface.widgets.runButton);
        interaface.layers.before.add(interaface.widgets.exportButton);
        
        var mapLeft = interaface.layers.before;
        var mapRight = interaface.layers.after;
        
        mapLeft.setOptions('satellite')
        mapRight.setOptions('satellite')
        
        listParams.forEach(function(params){
          
            var year = params[0];
          
            var band = 'classification_' + (year - 1).toString()
            
            if (year >= 2020) { band = 'classification_2019' }
            
            var lulc = ee.Image(assetLulc).select(band)
          
            var roi = dt.layers().get(0).toGeometry();
            
            var styled = ee.FeatureCollection([ee.Feature(roi)]).style({
              color: 'red',        // Cor da borda
              width: 2,            // Espessura da borda
              fillColor: '00000000' // Transparente no interior
            });
              
            
            
            var dam = damCore.getDAM(params, roi).updateMask(lulc.eq(3));
            
            var damReg = dam.select('freq_dam')
            var damDefReg = dam.select('freq_dam_df').gt(0)
            
            var damDefOrig = damOrg.select('freq_dam_df').filter(ee.Filter.eq('year', params[0])).sum().gt(0);
            var damOrig = damOrg.select('freq_dam').filter(ee.Filter.eq('year', params[0])).sum();
            
            
            
            mapLeft.centerObject(roi, 11)
          
            mapLeft.addLayer(damReg, {
             bands:['freq_dam'],
              min:1, max:10,
              palette:palettes.cmocean.Thermal[7]
            }, 'DEGRADAÇÃO - REG');
            
            mapLeft.addLayer(damDefReg.selfMask(), {
              min:0, max:1,
              bands:['freq_dam_df'],
              palette:['red']
            }, 'DESMATAMENTO - REG'); 
            
            
            mapRight.addLayer(damOrig, {
              min:1, max:10,
              palette:palettes.cmocean.Thermal[7]
            }, 'DEGRADAÇÃO');
            
            mapRight.addLayer(damDefOrig, {
              min:0, max:1,
              palette: ['red']
            }, 'DESMATAMENTO');       
      
      
            mapLeft.addLayer(styled, {}, 'ROI');
            mapRight.addLayer(styled, {}, 'ROI'); 
            
        });
        
    },
    
    init: function(){
        
        var linked = ui.Map.Linker([interaface.layers.before, interaface.layers.after])
        
        var splitPanels = ui.SplitPanel({
            firstPanel: interaface.layers.before,
            secondPanel: interaface.layers.after,
            orientation: 'horizontal',
            wipe: true,
            style: {stretch: 'both'}
        })
        
        interaface.widgets.init();
        
        interaface.initializeLayers(); 
        
        ui.root.widgets().reset([splitPanels]);
    },
}

// ======================================================================================================

interaface.init();

// ======================================================================================================

var vis = {
    'min': 1,
    'max': 5, // o max ´é em  torno de 40, mas para melhorar a vis, foi setado 12
    'palette': palettes.matplotlib.inferno[7].slice(1),
    'format': 'png'
};

// ======================================================================================================





