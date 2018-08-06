package maars.cellAnalysis;

import fiji.plugin.trackmate.action.ExportTracksToXML;
import fiji.plugin.trackmate.detection.LogDetectorFactory;
import fiji.plugin.trackmate.features.FeatureFilter;
import fiji.plugin.trackmate.features.edges.EdgeTargetAnalyzer;
import fiji.plugin.trackmate.features.edges.EdgeTimeLocationAnalyzer;
import fiji.plugin.trackmate.features.edges.EdgeVelocityAnalyzer;
import fiji.plugin.trackmate.features.spot.*;
import fiji.plugin.trackmate.features.track.*;
import fiji.plugin.trackmate.io.TmXmlWriter;
import fiji.plugin.trackmate.tracking.LAPUtils;
import fiji.plugin.trackmate.tracking.sparselap.SparseLAPTrackerFactory;
import fiji.plugin.trackmate.visualization.hyperstack.HyperStackDisplayer;
import ij.IJ;
import ij.ImagePlus;
import ij.plugin.Duplicator;
import maars.agents.DefaultSetOfCells;
import fiji.plugin.trackmate.Model;
import fiji.plugin.trackmate.Settings;
import fiji.plugin.trackmate.TrackMate;
import fiji.plugin.trackmate.SelectionModel;
import fiji.plugin.trackmate.Logger;

import java.io.File;
import java.io.IOException;

public class TrackmateAnalyzer implements Runnable {

    private ImagePlus[] targetImgs;
    private Duplicator duplicator = new Duplicator();
    private double radius;
    private double quality;
    private double maxNbSpot;
    private String channel;

    public TrackmateAnalyzer(ImagePlus fluoImage, DefaultSetOfCells soc, String ch, int maxNb,
                             double rad, double qual){
        this.targetImgs = new ImagePlus[soc.size()];
        for (int i = 1 ; i <= soc.size(); i++){
            fluoImage.setRoi(soc.getCell(i).getCellShapeRoi());
            targetImgs[i-1] = duplicator.run(fluoImage);
        }
        this.radius = rad;
        this.quality = qual;
        this.maxNbSpot = maxNb;
        this.channel = ch;
    }

    @Override
    public void run() {
        for (int i = 0; i < targetImgs.length; i++) {
            Model model = new Model();
//            model.setLogger(Logger.IJ_LOGGER);

            Settings settings = new Settings();
            settings.setFrom(targetImgs[i]);

            settings.detectorFactory = new LogDetectorFactory<>();
            settings.detectorSettings.put("DO_SUBPIXEL_LOCALIZATION", true);
            settings.detectorSettings.put("RADIUS", this.radius);
            settings.detectorSettings.put("THRESHOLD", 0.);
            settings.detectorSettings.put("DO_MEDIAN_FILTERING", true);

            settings.addSpotAnalyzerFactory(new SpotIntensityAnalyzerFactory<>());
            settings.addSpotAnalyzerFactory(new SpotContrastAndSNRAnalyzerFactory<>());
            settings.addSpotAnalyzerFactory(new SpotMorphologyAnalyzerFactory<>());
            settings.addSpotAnalyzerFactory(new SpotRadiusEstimatorFactory<>());
            settings.addSpotAnalyzerFactory(new SpotContrastAnalyzerFactory<>());


            FeatureFilter filter1 = new FeatureFilter("QUALITY", quality, true);
            settings.addSpotFilter(filter1);

            settings.addEdgeAnalyzer(new EdgeVelocityAnalyzer());
            settings.addEdgeAnalyzer(new EdgeTargetAnalyzer());
            settings.addEdgeAnalyzer(new EdgeTimeLocationAnalyzer());

            settings.trackerFactory = new SparseLAPTrackerFactory();
            settings.trackerSettings = LAPUtils.getDefaultLAPSettingsMap();
            settings.trackerSettings.put("ALLOW_TRACK_SPLITTING", true);
            settings.trackerSettings.put("ALLOW_TRACK_MERGING", true);

            settings.addTrackAnalyzer(new TrackDurationAnalyzer());
            settings.addTrackAnalyzer(new TrackBranchingAnalyzer());
            settings.addTrackAnalyzer(new TrackIndexAnalyzer());
            settings.addTrackAnalyzer(new TrackLocationAnalyzer());
            settings.addTrackAnalyzer(new TrackSpeedStatisticsAnalyzer());
            settings.addTrackAnalyzer(new TrackSpotQualityFeatureAnalyzer());

//            FeatureFilter filter2 = new FeatureFilter("TRACK_DISPLACEMENT", 10, true);
//            settings.addTrackFilter(filter2);
            TrackMate trackmate = new TrackMate(model, settings);

            boolean ok = trackmate.checkInput();
            if (!ok) {
                IJ.error(trackmate.getErrorMessage());
            }

            ok = trackmate.process();
            if (!ok) {
                IJ.error(trackmate.getErrorMessage());
            }

            SelectionModel selectionModel = new SelectionModel(model);
            HyperStackDisplayer displayer = new HyperStackDisplayer(model, selectionModel, targetImgs[i]);
            displayer.render();
            displayer.refresh();

            File outputFile = new File("/media/tongli/0ABC6EF952B52BF5/tmp/" + i + "_" + channel);
            TmXmlWriter writer = new TmXmlWriter(outputFile);
            writer.appendModel(model);
//            writer.appendSettings(settings);
            try {
                writer.writeToFile();
            } catch (IOException e) {
                e.printStackTrace();
            }

        }
    }
}
