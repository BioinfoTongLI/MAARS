package org.micromanager.maars;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.regex.Pattern;

import org.micromanager.cellstateanalysis.FluoAnalyzer;
import org.micromanager.cellstateanalysis.SetOfCells;
import org.micromanager.internal.utils.ReportingUtils;
import org.micromanager.utils.FileUtils;

import ij.IJ;
import ij.ImagePlus;
import ij.measure.Calibration;
import ij.plugin.frame.RoiManager;
import mmcorej.CMMCore;

/**
 * @author Tong LI, mail: tongli.bioinfo@gmail.com
 * @version Nov 22, 2015
 */
public class MAARSNoAcq {
	private PrintStream curr_err;
	private PrintStream curr_out;

	public MAARSNoAcq(CMMCore mmc, MaarsParameters parameters, SetOfCells soc) {
		ExecutorService es = null;
		// Start time
		long start = System.currentTimeMillis();

		// Acquisition path arrangement
		ExplorationXYPositions explo = new ExplorationXYPositions(mmc, parameters);

		for (int i = 0; i < explo.length(); i++) {
			System.out.println("x : " + explo.getX(i) + " y : " + explo.getY(i));
			String xPos = String.valueOf(Math.round(explo.getX(i)));
			String yPos = String.valueOf(Math.round(explo.getY(i)));
			String pathToSegDir = FileUtils
					.convertPath(parameters.getSavingPath() + "/movie_X" + xPos + "_Y" + yPos + "/");
			String pathToSegMovie = FileUtils.convertPath(pathToSegDir + "MMStack.ome.tif");
			ImagePlus segImg = null;
			try {
				segImg = IJ.openImage(pathToSegMovie);
			} catch (Exception e) {
				e.printStackTrace();
				System.out.println("Invalid path");
			}
			// --------------------------segmentation-----------------------------//
			MaarsSegmentation ms = new MaarsSegmentation(parameters, xPos, yPos);
			ms.segmentation(segImg);
			if (ms.roiDetected()) {
				// from Roi initialize a set of cell
				soc.loadCells(xPos, yPos);
				soc.setRoiMeasurementIntoCells(ms.getRoiMeasurements());
				// Get the focus slice of BF image
				Calibration bfImgCal = segImg.getCalibration();
				// ----------------start acquisition and analysis --------//
				try {
					PrintStream ps = new PrintStream(pathToSegDir + "CellStateAnalysis.LOG");
					curr_err = System.err;
					curr_out = System.err;
					System.setOut(ps);
					System.setErr(ps);
				} catch (FileNotFoundException e) {
					e.printStackTrace();
				}
				int frame = 0;
				String pathToFluoDir = parameters.getSavingPath() + "/movie_X" + xPos + "_Y" + yPos + "_FLUO/";
				String[] listAcqNames = new File(pathToFluoDir).list();
				String pattern = "(\\d+)(_)(\\w+)";
				int frameCounter = 0;
				for (String acqName : listAcqNames) {
					if (Pattern.matches(pattern, acqName)) {
						frameCounter++;
					}
				}
				String channels = parameters.getUsingChannels();
				String[] arrayChannels = channels.split(",", -1);
				frameCounter = frameCounter / arrayChannels.length;
				int nThread = Runtime.getRuntime().availableProcessors();
				es = Executors.newFixedThreadPool(nThread);
				Future<?> future = null;
				while (frame < frameCounter) {
					for (String channel : arrayChannels) {
						ReportingUtils.logMessage("Analysing channel " + channel + "_" + frame);
						String[] id = new String[] { xPos, yPos, String.valueOf(frame), channel };
						soc.addAcqID(id);
						String pathToFluoMovie = parameters.getSavingPath() + "/movie_X" + xPos + "_Y" + yPos + "_FLUO/"
								+ frame + "_" + channel + "/MMStack.ome.tif";
						ImagePlus fluoImage = IJ.openImage(pathToFluoMovie);
						future = es
								.submit(new FluoAnalyzer(fluoImage, bfImgCal, soc, channel,
										Integer.parseInt(parameters.getChMaxNbSpot(channel)), Double
												.parseDouble(parameters.getChSpotRaius(channel)),
										frame, Double.parseDouble(
												parameters.getFluoParameter(MaarsParameters.TIME_INTERVAL))));
					}
					frame++;
				}
				try {
					future.get();
				} catch (InterruptedException e1) {
					e1.printStackTrace();
				} catch (ExecutionException e1) {
					e1.printStackTrace();
				}
			}
			// RoiManager.getInstance().reset();
			// RoiManager.getInstance().close();
		}
		es.shutdown();
		try {
			es.awaitTermination(120, TimeUnit.MINUTES);
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		System.setErr(curr_err);
		System.setOut(curr_out);
		if (soc.size() != 0) {
			long startWriting = System.currentTimeMillis();
			soc.saveSpots();
			soc.saveFeatures();
			soc.saveCroppedImgs();
			ReportingUtils.logMessage("it took " + (double) (System.currentTimeMillis() - startWriting) / 1000
					+ " sec for writing results");
		}
		System.out.println("it took " + (double) (System.currentTimeMillis() - start) / 1000 + " sec for analysing");
	}
}
