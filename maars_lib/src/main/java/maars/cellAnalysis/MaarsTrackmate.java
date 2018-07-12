package maars.cellAnalysis;

import fiji.plugin.trackmate.Model;
import fiji.plugin.trackmate.SelectionModel;
import fiji.plugin.trackmate.Settings;
import fiji.plugin.trackmate.TrackMate;
import fiji.plugin.trackmate.detection.LogDetectorFactory;
import fiji.plugin.trackmate.features.spot.*;
import fiji.plugin.trackmate.visualization.hyperstack.HyperStackDisplayer;
import ij.IJ;
import ij.ImagePlus;
import maars.utils.ImgUtils;

import java.util.HashMap;
import java.util.Map;

import static fiji.plugin.trackmate.detection.DetectorKeys.*;

public class MaarsTrackmate {

   private Settings settings;

   public MaarsTrackmate(ImagePlus img, double radius, double quality, boolean gaussian_blur) {
      if (gaussian_blur){
         IJ.run(img, "Gaussian Blur 3D...", "x=1.2 y=1.2 z=2.4");
         while (IJ.macroRunning()){
            IJ.wait(200);
         }
      }

      settings = new Settings();
      settings.setFrom(ImgUtils.zProject(img, img.getCalibration()));

      // Computer different features (in order)

      settings.addSpotAnalyzerFactory(new SpotIntensityAnalyzerFactory<>());
      settings.addSpotAnalyzerFactory(new SpotContrastAndSNRAnalyzerFactory<>());
      settings.addSpotAnalyzerFactory(new SpotMorphologyAnalyzerFactory<>());
      settings.addSpotAnalyzerFactory(new SpotRadiusEstimatorFactory<>());
      settings.addSpotAnalyzerFactory(new SpotContrastAnalyzerFactory<>());

      // Set up detection parameters.

      settings.detectorFactory = new LogDetectorFactory<>();
      Map<String, Object> detectorSettings = new HashMap<>();
      detectorSettings.put(KEY_DO_SUBPIXEL_LOCALIZATION, true);
      detectorSettings.put(KEY_RADIUS, radius);
      detectorSettings.put(KEY_TARGET_CHANNEL, DEFAULT_TARGET_CHANNEL);
      detectorSettings.put(KEY_THRESHOLD, quality);
      detectorSettings.put(KEY_DO_MEDIAN_FILTERING, true);
      settings.detectorSettings = detectorSettings;
   }

   /**
    * Take parameters in the constructor then initalize trakemate object to get
    * unfiltered spots.
    *
    * @return Model a Trackmate style data structure
    */
   public Model doDetection() {
      TrackMate trackmate = new TrackMate(settings);

      trackmate.execDetection();

      trackmate.execInitialSpotFiltering();

      trackmate.computeSpotFeatures(false);

      return trackmate.getModel();
   }

   public static void executeTrackmate(ImagePlus img, double spotRadius, double quality) {
      assert img != null;
      MaarsTrackmate tmTest = new MaarsTrackmate(img.duplicate(), spotRadius, quality, true);
      Model model = tmTest.doDetection();
      model.getSpots().setVisible(true);
      SelectionModel selectionModel = new SelectionModel(model);
      HyperStackDisplayer displayer = new HyperStackDisplayer(model, selectionModel, tmTest.settings.imp);
      IJ.run(tmTest.settings.imp, "Enhance Contrast", "saturated=0.35");
      displayer.render();
      displayer.refresh();
   }
}
