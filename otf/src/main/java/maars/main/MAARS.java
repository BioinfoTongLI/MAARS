package maars.main;

import maars.acquisition.MAARS_mda;
import maars.cellAnalysis.FluoAnalyzer;
import maars.agents.SetOfCells;
import maars.display.SOCVisualizer;
import maars.gui.MaarsMainDialog;
import maars.utils.FileUtils;
import maars.io.IOUtils;
import maars.utils.ImgUtils;
import ij.IJ;
import ij.ImagePlus;
import ij.plugin.frame.RoiManager;
import mmcorej.CMMCore;
import org.micromanager.internal.MMStudio;
import org.micromanager.internal.utils.ReportingUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;
import java.util.concurrent.*;

/**
 * Main MAARS program
 *
 * @author Tong LI, mail: tongli.bioinfo@gmail.com
 * @version Nov 21, 2015
 */
public class MAARS implements Runnable {
   private static boolean stop_ = false;
   private PrintStream curr_err;
   private PrintStream curr_out;
   private MMStudio mm;
   private CMMCore mmc;
   private MaarsParameters parameters;
   private ArrayList<SetOfCells> socList_;
   private ArrayList<SOCVisualizer> socVisualizerList_;
   private CopyOnWriteArrayList<Map<String, Future>> tasksSet_;

   /**
    * * Constructor
    *
    * @param mm                MMStudio object (gui)
    * @param mmc               CMMCore object (core)
    * @param parameters        MAARS parameters object
    * @param socVisualizerList list of set of cell visualizer
    * @param tasksSet          tasks to be terminated
    * @param socList           list of set of cell
    */
   public MAARS(MMStudio mm, CMMCore mmc, MaarsParameters parameters, ArrayList<SOCVisualizer> socVisualizerList,
                CopyOnWriteArrayList<Map<String, Future>> tasksSet,
                ArrayList<SetOfCells> socList) {
      this.mmc = mmc;
      this.parameters = parameters;
      socList_ = socList;
      this.mm = mm;
      tasksSet_ = tasksSet;
      socVisualizerList_ = socVisualizerList;
   }

   public void interrupt() {
      stop_ = true;
   }

   @Override
   public void run() {
      String BF = "BF";
      String FLUO = "FLUO";
      // Start time
      long start = System.currentTimeMillis();
      parameters.setCalibration(String.valueOf(mm.getCachedPixelSizeUm()));
      ArrayList<String> arrayChannels = new ArrayList<>();
      Collections.addAll(arrayChannels, parameters.getUsingChannels().split(",", -1));
//        for (SetOfCells soc : socList_){
//           soc.reset();
//        }
      String savingPath = FileUtils.convertPath(parameters.getSavingPath());
      //acquisition

      ExecutorService es = Executors.newSingleThreadExecutor();
      HashMap<Integer, ImagePlus[]> segImgs = null;
      try {
         segImgs = es.submit(new MAARS_mda(
               parameters.getSegmentationParameter(MaarsParameters.PATH_TO_BF_ACQ_SETTING),
               savingPath, BF)).get();
      } catch (InterruptedException | ExecutionException e) {
         IOUtils.printErrorToIJLog(e);
      }
      if (segImgs == null) {
         IJ.log("No images acquired");
         return;
      }
      //update saving path
      parameters.setSavingPath(savingPath + File.separator + BF + "_1");
      MaarsSegmentation ms;
      ArrayList<MaarsSegmentation> arrayMs = new ArrayList<>();
      for (Integer posNb : segImgs.keySet()) {
         ImagePlus segImg = segImgs.get(posNb)[0];
         // --------------------------segmentation-----------------------------//
         ms = new MaarsSegmentation(parameters, segImg, posNb);
         arrayMs.add(ms);
         try {
            es.submit(ms).get();
         } catch (InterruptedException | ExecutionException e) {
            IOUtils.printErrorToIJLog(e);
         }
      }
      es.shutdown();
      try {
         es.awaitTermination(3, TimeUnit.MINUTES);
      } catch (InterruptedException e) {
         IOUtils.printErrorToIJLog(e);
      }
      parameters.setSavingPath(savingPath);
      // from Roi initialize a set of cell
      for (Integer posNb : segImgs.keySet()) {
         SetOfCells soc = socList_.get(posNb);
         soc.reset();
         soc.loadCells(savingPath + File.separator + BF + "_1", posNb);
         soc.setRoiMeasurementIntoCells(arrayMs.get(posNb).getRoiMeasurements());
      }
      // ----------------start acquisition and analysis --------//
      redirectLog(savingPath);
      int frame = 0;
      Boolean do_analysis = Boolean.parseBoolean(parameters.getFluoParameter(MaarsParameters.DO_ANALYSIS));
      double fluoTimeInterval = Double.parseDouble(parameters.getFluoParameter(MaarsParameters.TIME_INTERVAL));
      ExecutorService es1 = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
      if (parameters.useDynamic()) {
         // being dynamic acquisition
         double startTime = System.currentTimeMillis();
         double timeLimit = Double.parseDouble(parameters.getFluoParameter(MaarsParameters.TIME_LIMIT)) * 60
               * 1000;
         while (System.currentTimeMillis() - startTime <= timeLimit) {
            double beginAcq = System.currentTimeMillis();
            if (stop_) {
               es1.shutdownNow();
               break;
            }
            // Section to acquire bf images between fluos, can be skipped
            if (frame != 0) {
               try {
                  segImgs = es1.submit(new MAARS_mda(parameters.getSegmentationParameter(MaarsParameters.PATH_TO_BF_ACQ_SETTING),
                        savingPath, BF)).get();
               } catch (InterruptedException | ExecutionException e) {
                  IOUtils.printErrorToIJLog(e);
               }
               parameters.setSavingPath(savingPath + File.separator + BF + "_" + String.valueOf(frame + 1));
               for (Integer posNb : segImgs.keySet()) {
                  ImagePlus segImg = segImgs.get(posNb)[0];
                  //update saving path
//                          parameters.setSavingPath(savingPath + File.separator + BF + "_"+posNb);
                  // --------------------------segmentation-----------------------------//
                  ms = new MaarsSegmentation(parameters, segImg, posNb);
                  try {
                     es1.submit(ms).get();
                  } catch (InterruptedException | ExecutionException e) {
                     IOUtils.printErrorToIJLog(e);
                  }
               }
            }

            Map<String, Future> channelsInFrame = new HashMap<>();
            HashMap<Integer, ImagePlus[]> fluos = null;
            try {
               fluos = es1.submit(new MAARS_mda(
                     parameters.getFluoParameter(MaarsParameters.PATH_TO_FLUO_ACQ_SETTING),
                     savingPath, FLUO)).get();
            } catch (InterruptedException | ExecutionException e) {
               IOUtils.printErrorToIJLog(e);
            }
            for (Integer posNb : fluos.keySet()) {
               for (ImagePlus chImp:fluos.get(posNb)){
                  //TODO
                  if (do_analysis) {
                     String channel = chImp.getTitle();
                     Future future2 = es1.submit(new FluoAnalyzer(chImp, chImp.getCalibration(), socList_.get(posNb), channel,
                           Integer.parseInt(parameters.getChMaxNbSpot(channel)),
                           Double.parseDouble(parameters.getChSpotRaius(channel)),
                           Double.parseDouble(parameters.getChQuality(channel)), frame, socVisualizerList_.get(posNb), parameters.useDynamic()));
                     channelsInFrame.put(channel, future2);
                  }
               }
            }

            tasksSet_.add(channelsInFrame);
            frame++;
            double acqTook = System.currentTimeMillis() - beginAcq;
            if (fluoTimeInterval > acqTook) {
               try {
                  Thread.sleep((long) (fluoTimeInterval - acqTook));
               } catch (InterruptedException e) {
                  IOUtils.printErrorToIJLog(e);
               }
            } else {
               IJ.log("Attention : acquisition before took longer than " + fluoTimeInterval
                     / 1000 + " s : " + acqTook);
            }
         }
         IJ.log("Acquisition Done, proceeding to post-analysis");
      } else {
         // being static acquisition
         Map<String, Future> channelsInFrame = new HashMap<>();
         HashMap<Integer, ImagePlus[]> fluos = null;
         try {
            fluos = es1.submit(new MAARS_mda(
                  parameters.getFluoParameter(MaarsParameters.PATH_TO_FLUO_ACQ_SETTING),
                  savingPath, FLUO)).get();
         } catch (InterruptedException | ExecutionException e) {
            IOUtils.printErrorToIJLog(e);
         }
         for (Integer posNb : fluos.keySet()) {
            for (ImagePlus chImp:fluos.get(posNb)){
               String channel = chImp.getTitle();
               if (do_analysis) {
                  Future future2 = es1.submit(new FluoAnalyzer(chImp, chImp.getCalibration(), socList_.get(posNb), channel,
                        Integer.parseInt(parameters.getChMaxNbSpot(channel)),
                        Double.parseDouble(parameters.getChSpotRaius(channel)),
                        Double.parseDouble(parameters.getChQuality(channel)), frame, socVisualizerList_.get(posNb), parameters.useDynamic()));
                  channelsInFrame.put(channel, future2);
               }
            }
         }
         tasksSet_.add(channelsInFrame);
      }
      parameters.setSavingPath(savingPath);
      es1.shutdown();
      try {
         es1.awaitTermination(60, TimeUnit.MINUTES);
      } catch (InterruptedException e) {
         e.printStackTrace();
      }
      Maars_Interface.waitAllTaskToFinish(tasksSet_);
      RoiManager.getInstance().reset();
      RoiManager.getInstance().close();
      for (Integer posNb : segImgs.keySet()) {
         SetOfCells soc = socList_.get(posNb);
         if (do_analysis) {
            long startWriting = System.currentTimeMillis();
            String pathToFluoDir = savingPath + File.separator + FLUO + "_1";
            ImagePlus mergedImg = ImgUtils.loadFullFluoImgs(pathToFluoDir);
            mergedImg.getCalibration().frameInterval = fluoTimeInterval / 1000;
            IOUtils.saveAll(soc, mergedImg, pathToFluoDir, parameters.useDynamic(), arrayChannels);
            if (parameters.useDynamic()) {
               if (IJ.isWindows()) {
                  savingPath = FileUtils.convertPathToLinuxType(savingPath + File.separator + BF + "_1");
               }
               Maars_Interface.analyzeMitosisDynamic(soc, parameters, savingPath);
            }
            ReportingUtils.logMessage("it took " + (double) (System.currentTimeMillis() - startWriting) / 1000
                  + " sec for writing results");
         }
//                RemoteNotification.mailNotify("tongli.bioinfo@gmail.com");
         mmc.setAutoShutter(true);
         System.setErr(curr_err);
         System.setOut(curr_out);
         soc.reset();
      }
      IJ.log("it took " + (double) (System.currentTimeMillis() - start) / 1000 + " sec for analysing all fields");
      System.gc();
      IJ.showMessage("MAARS: Done!");
      MaarsMainDialog.okMainDialogButton.setEnabled(true);
   }

   private void redirectLog(String savingPath) {
      try {
         PrintStream ps = new PrintStream(savingPath + File.separator + "FluoAnalysis.LOG");
         curr_err = System.err;
         curr_out = System.err;
         System.setOut(ps);
         System.setErr(ps);
      } catch (FileNotFoundException e) {
         IOUtils.printErrorToIJLog(e);
      }
   }
}