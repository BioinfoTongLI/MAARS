package org.micromanager.cellstateanalysis;

import static fiji.plugin.trackmate.detection.DetectorKeys.DEFAULT_TARGET_CHANNEL;
import static fiji.plugin.trackmate.detection.DetectorKeys.KEY_DO_MEDIAN_FILTERING;
import static fiji.plugin.trackmate.detection.DetectorKeys.KEY_DO_SUBPIXEL_LOCALIZATION;
import static fiji.plugin.trackmate.detection.DetectorKeys.KEY_RADIUS;
import static fiji.plugin.trackmate.detection.DetectorKeys.KEY_TARGET_CHANNEL;
import static fiji.plugin.trackmate.detection.DetectorKeys.KEY_THRESHOLD;

import java.util.HashMap;
import java.util.Map;

import fiji.plugin.trackmate.Model;
import fiji.plugin.trackmate.Settings;
import fiji.plugin.trackmate.TrackMate;
import fiji.plugin.trackmate.detection.LogDetectorFactory;
import fiji.plugin.trackmate.features.spot.SpotContrastAnalyzerFactory;
import fiji.plugin.trackmate.features.spot.SpotContrastAndSNRAnalyzerFactory;
import fiji.plugin.trackmate.features.spot.SpotIntensityAnalyzerFactory;
import fiji.plugin.trackmate.features.spot.SpotMorphologyAnalyzerFactory;
import fiji.plugin.trackmate.features.spot.SpotRadiusEstimatorFactory;
import ij.ImagePlus;
import net.imglib2.type.numeric.real.FloatType;

public class MaarsTrackmate {

	private Model model;
	private Settings settings;

	public MaarsTrackmate(ImagePlus img, double radius) {
		img.deleteRoi();
		model = new Model();

		settings = new Settings();
		settings.setFrom(img);

		// Computer different features (in order)

		settings.addSpotAnalyzerFactory(new SpotRadiusEstimatorFactory<FloatType>());
		settings.addSpotAnalyzerFactory(new SpotContrastAnalyzerFactory<FloatType>());
		settings.addSpotAnalyzerFactory(new SpotIntensityAnalyzerFactory<FloatType>());
		settings.addSpotAnalyzerFactory(new SpotMorphologyAnalyzerFactory<FloatType>());
		settings.addSpotAnalyzerFactory(new SpotContrastAndSNRAnalyzerFactory<FloatType>());

		// Set up detection parameters.

		settings.detectorFactory = new LogDetectorFactory<FloatType>();
		Map<String, Object> detectorSettings = new HashMap<String, Object>();
		detectorSettings.put(KEY_DO_SUBPIXEL_LOCALIZATION, true);
		detectorSettings.put(KEY_RADIUS, radius);
		detectorSettings.put(KEY_TARGET_CHANNEL, DEFAULT_TARGET_CHANNEL);
		// TODO to figure out what value to use
		detectorSettings.put(KEY_THRESHOLD, (double) 1);
		detectorSettings.put(KEY_DO_MEDIAN_FILTERING, true);
		settings.detectorSettings = detectorSettings;
	}

	/**
	 * Take parameters in the constructor then initalize trakemate object to get
	 * unfiltered spots.
	 */
	public Model doDetection() {
		TrackMate trackmate = new TrackMate(model, settings);

		trackmate.execDetection();

		trackmate.execInitialSpotFiltering();

		trackmate.computeSpotFeatures(false);

		return trackmate.getModel();
	}
}