package maars.cellAnalysis;

import fiji.plugin.trackmate.Logger;
import fiji.plugin.trackmate.detection.DetectorKeys;
import fiji.plugin.trackmate.detection.LogDetectorFactory;
import fiji.plugin.trackmate.features.FeatureFilter;
import fiji.plugin.trackmate.features.edges.EdgeTargetAnalyzer;
import fiji.plugin.trackmate.features.edges.EdgeTimeLocationAnalyzer;
import fiji.plugin.trackmate.features.edges.EdgeVelocityAnalyzer;
import fiji.plugin.trackmate.features.spot.*;
import fiji.plugin.trackmate.features.track.*;
import fiji.plugin.trackmate.tracking.LAPUtils;
import fiji.plugin.trackmate.tracking.TrackerKeys;
import fiji.plugin.trackmate.tracking.oldlap.SimpleLAPTrackerFactory;
import ij.IJ;
import ij.ImagePlus;
import ij.plugin.Duplicator;
import maars.agents.DefaultSetOfCells;
import fiji.plugin.trackmate.Model;
import fiji.plugin.trackmate.Settings;
import fiji.plugin.trackmate.TrackMate;
import maars.utils.ImgUtils;

public class TrackmateAnalyzer implements Runnable {

    private ImagePlus[] targetImgs;
    private double radius;
    private double quality;
    private String channel;
    private DefaultSetOfCells soc;
    private double minMitosisDuration;

    public TrackmateAnalyzer(ImagePlus fluoImage, DefaultSetOfCells soc, String ch, double rad, double qual,
                             double minMitosisDuration){
        this.targetImgs = new ImagePlus[soc.size()];
        for (int i = 1 ; i <= soc.size(); i++){
            fluoImage.setRoi(soc.getCell(i).getCellShapeRoi());
            Duplicator duplicator = new Duplicator();
            ImagePlus temp_img = duplicator.run(fluoImage);
            temp_img.setRoi(ImgUtils.centerCroppedRoi(soc.getCell(i).getCellShapeRoi()));
            targetImgs[i-1] = temp_img;
        }
        this.radius = rad;
        this.quality = qual;
        this.channel = ch;
        this.soc = soc;
        this.minMitosisDuration = minMitosisDuration;
    }

    @Override
    public void run() {
        for (int i = 0; i < targetImgs.length; i++) {
            Model model = new Model();
            model.setLogger(Logger.DEFAULT_LOGGER);

            Settings settings = new Settings();
            settings.setFrom(targetImgs[i]);

            settings.detectorFactory = new LogDetectorFactory<>();
            settings.detectorSettings.put(DetectorKeys.KEY_DO_SUBPIXEL_LOCALIZATION, true);
            settings.detectorSettings.put(DetectorKeys.KEY_RADIUS, this.radius);
            settings.detectorSettings.put(DetectorKeys.KEY_THRESHOLD, quality);
            settings.detectorSettings.put(DetectorKeys.KEY_DO_MEDIAN_FILTERING, true);

            settings.addSpotAnalyzerFactory(new SpotIntensityAnalyzerFactory<>());
            settings.addSpotAnalyzerFactory(new SpotContrastAndSNRAnalyzerFactory<>());
            settings.addSpotAnalyzerFactory(new SpotMorphologyAnalyzerFactory<>());
            settings.addSpotAnalyzerFactory(new SpotRadiusEstimatorFactory<>());
            settings.addSpotAnalyzerFactory(new SpotContrastAnalyzerFactory<>());


//            FeatureFilter filter1 = new FeatureFilter("QUALITY", quality, true);
//            settings.addSpotFilter(filter1);

            settings.addEdgeAnalyzer(new EdgeVelocityAnalyzer());
            settings.addEdgeAnalyzer(new EdgeTargetAnalyzer());
            settings.addEdgeAnalyzer(new EdgeTimeLocationAnalyzer());

            settings.trackerFactory = new SimpleLAPTrackerFactory();
            settings.trackerSettings = LAPUtils.getDefaultLAPSettingsMap();
//            settings.trackerSettings.put(TrackerKeys.KEY_ALLOW_TRACK_SPLITTING, true);
//            settings.trackerSettings.put(TrackerKeys.KEY_ALLOW_TRACK_MERGING, true);
//            settings.trackerSettings.put(TrackerKeys.KEY_MERGING_MAX_DISTANCE, 1d);
            settings.trackerSettings.put(TrackerKeys.KEY_LINKING_MAX_DISTANCE, 1.5d);
            settings.trackerSettings.put(TrackerKeys.KEY_GAP_CLOSING_MAX_DISTANCE, 1.5d);
            settings.trackerSettings.put(TrackerKeys.KEY_GAP_CLOSING_MAX_FRAME_GAP, 10);

            settings.addTrackAnalyzer(new TrackDurationAnalyzer());
            settings.addTrackAnalyzer(new TrackBranchingAnalyzer());
            settings.addTrackAnalyzer(new TrackIndexAnalyzer());
            settings.addTrackAnalyzer(new TrackLocationAnalyzer());
            settings.addTrackAnalyzer(new TrackSpeedStatisticsAnalyzer());
            settings.addTrackAnalyzer(new TrackSpotQualityFeatureAnalyzer());

            FeatureFilter filter2 = new FeatureFilter(TrackDurationAnalyzer.TRACK_DURATION,
                    minMitosisDuration, true);
            settings.addTrackFilter(filter2);

            TrackMate trackmate = new TrackMate(model, settings);

            boolean ok = trackmate.checkInput();
            if (!ok) {
                IJ.error(trackmate.getErrorMessage());
            }

            ok = trackmate.process();
            if (!ok) {
                IJ.log(trackmate.getErrorMessage() + "-- cell # : " + (i + 1) + "_" + this.channel);
            }else{
                soc.addPotentialMitosisCell(i+1);
                soc.getCell(i+1).putModel(channel, trackmate.getModel());
            }
        }
    }
}
