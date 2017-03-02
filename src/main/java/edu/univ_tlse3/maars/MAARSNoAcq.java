package edu.univ_tlse3.maars;

import edu.univ_tlse3.cellstateanalysis.Cell;
import edu.univ_tlse3.cellstateanalysis.FluoAnalyzer;
import edu.univ_tlse3.cellstateanalysis.SetOfCells;
import edu.univ_tlse3.display.SOCVisualizer;
import edu.univ_tlse3.gui.MaarsMainDialog;
import edu.univ_tlse3.utils.FileUtils;
import edu.univ_tlse3.utils.IOUtils;
import edu.univ_tlse3.utils.ImgUtils;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.Concatenator;
import ij.plugin.frame.RoiManager;

import java.awt.*;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.*;
import java.util.regex.Pattern;

/**
 * @author Tong LI, mail: tongli.bioinfo@gmail.com
 * @version Nov 22, 2015
 */
public class MAARSNoAcq implements Runnable {
   private PrintStream curr_err;
   private PrintStream curr_out;
   private MaarsParameters parameters;
   private SetOfCells soc_;
   private String rootDir;
   private ArrayList<String> arrayChannels = new ArrayList<>();
   private SOCVisualizer socVisualizer_;
   private ExecutorService es_;
   public boolean skipAllRestFrames = false;
   private CopyOnWriteArrayList<Map<String, Future>> tasksSet_;

   public MAARSNoAcq(MaarsParameters parameters, SOCVisualizer socVisualizer,
                     ExecutorService es, CopyOnWriteArrayList<Map<String, Future>> tasksSet,
                     SetOfCells soc) {
      this.parameters = parameters;
      tasksSet_ = tasksSet;
      rootDir = parameters.getSavingPath();
      soc_ = soc;
      socVisualizer_ = socVisualizer;
      es_ = es;
   }

   private ArrayList<String[]> getAcqPositions() {
      ArrayList<String[]> acqPos = new ArrayList<>();
      String[] listAcqNames = new File(rootDir).list();
      String pattern = "(X)(\\d+)(_)(Y)(\\d+)(_FLUO)";
      assert listAcqNames != null;
      for (String acqName : listAcqNames) {
         if (Pattern.matches(pattern, acqName)) {
            acqPos.add(new String[]{acqName.split("_", -1)[0].substring(1),
                    acqName.split("_", -1)[1].substring(1)});
         }
      }
      return acqPos;
   }

   @Override
   public void run() {
      // Start time
      long start = System.currentTimeMillis();
      for (String[] pos : getAcqPositions()) {
         if (skipAllRestFrames) {
            break;
         }
         soc_.reset();
         String xPos = pos[0];
         String yPos = pos[1];
         IJ.log("x : " + xPos + " y : " + yPos);
         String pathToSegDir = FileUtils.convertPath(rootDir + "/X" + xPos + "_Y" + yPos);
         String pathToFluoDir = pathToSegDir + "_FLUO/";
         String pathToSegMovie = FileUtils.convertPath(pathToSegDir + "/_1/_1_MMStack_Pos0.ome.tif");
         //update saving path
         parameters.setSavingPath(pathToSegDir);
         Boolean skipSegmentation = Boolean.parseBoolean(parameters.getSkipSegmentation());
         ImagePlus segImg = null;
         try {
            segImg = IJ.openImage(pathToSegMovie);
            parameters.setCalibration(String.valueOf(segImg.getCalibration().pixelWidth));
         } catch (Exception e) {
            IOUtils.printErrorToIJLog(e);
         }
         // --------------------------segmentation-----------------------------//
         MaarsSegmentation ms = new MaarsSegmentation(parameters, segImg);
         Future future;
         if (!skipSegmentation) {
            future = es_.submit(ms);
            try {
               future.get();
            } catch (InterruptedException | ExecutionException e) {
               IOUtils.printErrorToIJLog(e);
            }
         }
         if (ms.roiDetected()) {
            soc_.reset();
            // from Roi.zip initialize a set of cell
            soc_.loadCells(pathToSegDir);
            ResultsTable rt;
            if (!skipSegmentation) {
               rt = ms.getRoiMeasurements();
            }else{
               IJ.open(pathToSegDir + File.separator +"BF_Results.csv");
               rt = ResultsTable.getResultsTable();
               ResultsTable.getResultsWindow().close(false);
            }
            soc_.setRoiMeasurementIntoCells(rt);

            Calibration bfImgCal = null;
            assert segImg != null;
            bfImgCal = segImg.getCalibration();
            segImg = null;
            // ----------------start acquisition and analysis --------//
            try {
               PrintStream ps = new PrintStream(pathToSegDir + "/CellStateAnalysis.LOG");
               curr_err = System.err;
               curr_out = System.err;
               System.setOut(ps);
               System.setErr(ps);
            } catch (FileNotFoundException e) {
               IOUtils.printErrorToIJLog(e);
            }
            String[] listAcqNames = new File(pathToFluoDir).list();
            String pattern = "(\\w+)(_)(\\d+)";
            ArrayList<Integer> arrayImgFrames = new ArrayList<>();
            assert listAcqNames != null;
            for (String acqName : listAcqNames) {
               if (Pattern.matches(pattern, acqName)) {
                  String current_channel = acqName.split("_", -1)[0];
                  String current_frame = acqName.split("_", -1)[1];
                  if (!arrayChannels.contains(current_channel)) {
                     arrayChannels.add(current_channel);
                  }
                  if (!arrayImgFrames.contains(Integer.parseInt(current_frame))) {
                     arrayImgFrames.add(Integer.parseInt(current_frame));
                  }
               }
            }
            Collections.sort(arrayImgFrames);
            Concatenator concatenator = new Concatenator();
            concatenator.setIm5D(true);
            ImagePlus concatenatedFluoImgs = null;
            //TODO
            Boolean saveRam = null;
            HashMap<String, HashMap<Integer,ImagePlus>> allCroppedImgs = null;
            if (Runtime.getRuntime().maxMemory()/1024/1024/1024>16){
               saveRam = false;
            }else{
               saveRam = true;
               allCroppedImgs = new HashMap<>();
            }
            for (Integer current_frame : arrayImgFrames) {
               Map<String, Future> channelsInFrame = new HashMap<>();
               for (String channel : arrayChannels) {
                  IJ.log("Analysing channel " + channel + "_" + current_frame);
                  String pathToFluoMovie = pathToFluoDir + channel + "_" + current_frame + "/" + channel + "_" + current_frame + "_MMStack_Pos0.ome.tif";
                  ImagePlus fluoImage = IJ.openImage(pathToFluoMovie);
                  ImagePlus zProjectedFluoImg = ImgUtils.zProject(fluoImage);
                  zProjectedFluoImg.setTitle(fluoImage.getTitle() + "_" + channel + "_projected");
                  zProjectedFluoImg.setCalibration(fluoImage.getCalibration());
                  future = es_.submit(new FluoAnalyzer(zProjectedFluoImg.duplicate(), bfImgCal, soc_, channel,
                          Integer.parseInt(parameters.getChMaxNbSpot(channel)),
                          Double.parseDouble(parameters.getChSpotRaius(channel)),
                          Double.parseDouble(parameters.getChQuality(channel)), current_frame, socVisualizer_,
                          parameters.useDynamic()));
                  channelsInFrame.put(channel, future);
                  if (!saveRam) {
                     for (int i = 1; i <= fluoImage.getStack().size(); i++) {
                        fluoImage.getStack().setSliceLabel(channel + "_" + current_frame, i);
                     }
                     if (concatenatedFluoImgs == null) {
                        concatenatedFluoImgs = fluoImage;
                     } else {
                        concatenatedFluoImgs = concatenator.concatenate(concatenatedFluoImgs, fluoImage, false);
                     }
                  }else{
                     if (!allCroppedImgs.keySet().contains(channel)){
                        allCroppedImgs.put(channel, new HashMap<>());
                        for (Cell c : soc_){
                           int cellNb = c.getCellNumber();
                           if (!allCroppedImgs.get(channel).keySet().contains(cellNb)){
                              allCroppedImgs.get(channel).put(cellNb,
                                      ImgUtils.cropImgWithRoi(fluoImage, c.getCellShapeRoi()));
                           }
                        }
                     }else{
                        for (Cell c : soc_){
                           int cellNb = c.getCellNumber();
                           allCroppedImgs.get(channel).put(cellNb, concatenator.concatenate(allCroppedImgs.get(channel).get(cellNb),
                                   ImgUtils.cropImgWithRoi(fluoImage, c.getCellShapeRoi()), false));
                        }
                     }
                  }
                  fluoImage = null;
                  zProjectedFluoImg = null;
               }
               tasksSet_.add(channelsInFrame);
               if (skipAllRestFrames){
                  break;
               }
            }
            MaarsMainDialog.waitAllTaskToFinish(tasksSet_);
            if (!skipAllRestFrames) {
               RoiManager.getInstance().reset();
               RoiManager.getInstance().close();
               double timeInterval = Double.parseDouble(parameters.getFluoParameter(MaarsParameters.TIME_INTERVAL));
               if (soc_.size() != 0) {
                  long startWriting = System.currentTimeMillis();
                  concatenatedFluoImgs.getCalibration().frameInterval = timeInterval / 1000;
                  if (!saveRam) {
                     MAARS.saveAll(soc_, concatenatedFluoImgs, pathToFluoDir, parameters.useDynamic(),
                             arrayChannels.size(), arrayImgFrames.size());
                  }else{
                     MAARS.saveAll(soc_, allCroppedImgs, pathToFluoDir, parameters.useDynamic(),
                             arrayChannels.size(), arrayImgFrames.size());
                  }
                  IJ.log("it took " + (double) (System.currentTimeMillis() - startWriting) / 1000
                          + " sec for writing results");
                  if (parameters.useDynamic()) {
                     if (IJ.isWindows()) {
                        pathToSegDir = FileUtils.convertPathToLinuxType(pathToSegDir);
                     }
                     MAARS.analyzeMitosisDynamic(soc_, parameters, pathToSegDir, true);
                  }
               }else if (soc_.size() == 0){
                  try {
                     org.apache.commons.io.FileUtils.deleteDirectory(new File(pathToSegDir));
                  } catch (IOException e) {
                     IOUtils.printErrorToIJLog(e);
                  }
               }
            }
            concatenatedFluoImgs = null;
         }
      }
      System.setErr(curr_err);
      System.setOut(curr_out);
      if (!skipAllRestFrames){
         IJ.log("it took " + (double) (System.currentTimeMillis() - start) / 1000 + " sec for analysing all fields");
      }
      System.gc();
   }
}