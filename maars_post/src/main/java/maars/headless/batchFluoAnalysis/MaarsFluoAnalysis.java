package maars.headless.batchFluoAnalysis;

import ij.IJ;
import ij.ImagePlus;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.Duplicator;

import maars.agents.Cell;
import maars.agents.DefaultSetOfCells;
import maars.cellAnalysis.FluoAnalyzer;
import maars.cellAnalysis.PythonPipeline;
import maars.cellAnalysis.TrackmateAnalyzer;
import maars.display.SOCVisualizer;
import maars.io.IOUtils;
import maars.main.MaarsParameters;
import maars.main.Maars_Interface;
import maars.utils.FileUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.*;
import java.util.concurrent.CopyOnWriteArrayList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 * Created by tong on 30/06/17.
 */
public class MaarsFluoAnalysis implements Runnable{
   public static final String MITODIRNAME = "Mitosis";
   private static Logger logger = LoggerFactory.getLogger(MaarsFluoAnalysis.class);
   private MaarsParameters parameters_;
   private String[] usingChannels;
   private Duplicator duplicator = new Duplicator();
   private SOCVisualizer visualizer_;
   private ImagePlus[] processedChs_;
   private String posName_;
   private PrintStream curr_err = null;
   private PrintStream curr_out = null;
   private String segAnaDir;
   public static int METHOD = -1;
   private DefaultSetOfCells soc_;
   public MaarsFluoAnalysis(DefaultSetOfCells soc, ImagePlus[] impChs, String posName, MaarsParameters parameters, SOCVisualizer visualizer){
      soc_ = soc;
      processedChs_ = impChs;
      parameters_ = parameters;
      usingChannels = parameters_.getUsingChannels().split(",", -1);
      visualizer_ = visualizer;
      posName_ = posName;
      segAnaDir = FileUtils.convertPath(parameters.getSavingPath()) + File.separator +
            parameters.getSegmentationParameter(MaarsParameters.SEG_PREFIX) + Maars_Interface.SEGANALYSIS_SUFFIX;
      try {
         PrintStream ps = new PrintStream(parameters.getSavingPath() + File.separator + "FluoAnalysis.LOG");
         curr_err = System.err;
         curr_out = System.out;
         System.setOut(ps);
         System.setErr(ps);
      } catch (FileNotFoundException e) {
         IOUtils.printErrorToIJLog(e);
      }
   }
   @Override
   public void run() {
      String currentPosPrefix = segAnaDir + posName_ + File.separator;
      String currentZipPath = currentPosPrefix + "ROI.zip";
      if (FileUtils.exists(currentZipPath)) {
         // from Roi.zip initialize a set of cell
         soc_.loadCells(currentZipPath);
         IJ.open(currentPosPrefix + "Results.csv");
         ResultsTable rt = ResultsTable.getResultsTable();
         ResultsTable.getResultsWindow().close(false);
         soc_.addRoiMeasurementIntoCells(rt);
         // ----------------start acquisition and analysis --------//

         CopyOnWriteArrayList<Map<String, Future>> tasksSet = new CopyOnWriteArrayList<>();
         // main analysis step
         switch (METHOD){
            case 0:
               processStackedImg(processedChs_, parameters_, soc_, tasksSet);
               break;
            default:
               processStackedImg(processedChs_, parameters_, soc_, visualizer_, tasksSet);
         }
         Maars_Interface.waitAllTaskToFinish(tasksSet);
         if (soc_.size() != 0) {
            IOUtils.saveAll(METHOD, soc_, processedChs_, parameters_.getSavingPath() + File.separator, parameters_.useDynamic(),
                  usingChannels, posName_, parameters_.getFluoParameter(MaarsParameters.FLUO_PREFIX));
            analyzeMitosisDynamic(soc_, parameters_, processedChs_[0].getCalibration());
         }
      System.gc();
      }
      System.setErr(curr_err);
      System.setOut(curr_out);
   }


   private void processStackedImg(ImagePlus[] processedChs, MaarsParameters parameters, DefaultSetOfCells soc,
                                  CopyOnWriteArrayList<Map<String, Future>> tasksSet) {
      ExecutorService es = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
      Map<String, Future> chAnalysisTasks = new HashMap<>();
      for (int j = 0; j < processedChs.length; j++) {
         String channel = usingChannels[j];
         IJ.log("Processing channel " + channel );
         Future future = es.submit(new TrackmateAnalyzer(processedChs[j],
                 soc, channel, Double.parseDouble(parameters.getChSpotRaius(channel)),
                 Double.parseDouble(parameters.getChQuality(channel)),
                 Double.parseDouble(parameters.getMinimumMitosisDuration())));
         chAnalysisTasks.put(channel, future);
      }
      tasksSet.add(chAnalysisTasks);
   }


   private void processStackedImg(ImagePlus[] processedChs, MaarsParameters parameters, DefaultSetOfCells soc,
                                       SOCVisualizer socVisualizer, CopyOnWriteArrayList<Map<String, Future>> tasksSet){
      ExecutorService es = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());

      for (int frame = 1; frame <= processedChs[0].getDimensions()[4]; frame++) {
         Map<String, Future> chAnalysisTasks = new HashMap<>();
         for (int j = 0; j < processedChs.length; j++) {
            String channel = usingChannels[j];
            IJ.log("Processing channel " + channel + "_" + frame);
            ImagePlus slicedImg = duplicator.run(processedChs[j], frame, frame);
            Future future = es.submit(new FluoAnalyzer(slicedImg, slicedImg.getCalibration(),
                  soc, channel, Integer.parseInt(parameters.getChMaxNbSpot(channel)),
                  Double.parseDouble(parameters.getChSpotRaius(channel)),
                  Double.parseDouble(parameters.getChQuality(channel)), frame, socVisualizer,
                  parameters.useDynamic()));
            chAnalysisTasks.put(channel, future);
         }
         tasksSet.add(chAnalysisTasks);
      }
      es.shutdown();
   }

   private static void findAbnormalCells(String mitoDir,
                                         DefaultSetOfCells soc,
                                         HashMap slopeChanges) {
      assert FileUtils.exists(mitoDir);
      PrintWriter out = null;
      try {
         out = new PrintWriter(mitoDir + File.separator + "abnormalCells.txt");
      } catch (FileNotFoundException e) {
         IOUtils.printErrorToIJLog(e);
      }
      assert out != null;
      for (Object cellNb : slopeChanges.keySet()) {
         int cellNbInt = Integer.parseInt(String.valueOf(cellNb));
         int anaBOnsetFrame = Integer.valueOf(((String[]) slopeChanges.get(cellNb))[1]);
         int lastAnaphaseFrame = Integer.valueOf(((String[]) slopeChanges.get(cellNb))[3]);
         Cell cell = soc.getCell(cellNbInt);
         ArrayList<Integer> spotInBtwnFrames = cell.getSpotInBtwnFrames();
         if (spotInBtwnFrames.size() > 0) {
            Collections.sort(spotInBtwnFrames);
            int laggingTimePoint = spotInBtwnFrames.get(spotInBtwnFrames.size() - 1);
            if (laggingTimePoint > anaBOnsetFrame && laggingTimePoint < lastAnaphaseFrame) {
               String laggingMessage = "Lagging :" + cellNb + "_lastLaggingTimePoint_" + laggingTimePoint + "_anaBonset_" + anaBOnsetFrame;
               out.println(laggingMessage);
               logger.info(laggingMessage);
               IJ.openImage(mitoDir + File.separator + "croppedImgs"
                     + File.separator + cellNb + "_GFP.tif").show();
            }
         }
         //TODO to show unaligned cell
         if (cell.unalignedSpotFrames().size() > 0) {
            String unalignKtMessage = "Unaligned : Cell " + cellNb + " detected with unaligned kinetochore(s)";
            logger.info(unalignKtMessage);
            out.println(unalignKtMessage);
         }
      }
      out.close();
      IJ.log("lagging detection finished");
   }

   public static void analyzeMitosisDynamic(DefaultSetOfCells soc, MaarsParameters parameters, Calibration calib) {
      IJ.log("Start python analysis");
      String pos = soc.getPosLabel();
      String pathToRoot = parameters.getSavingPath() + File.separator;
      String mitoDir = pathToRoot + parameters.getFluoParameter(MaarsParameters.FLUO_PREFIX)+"_"+MITODIRNAME
            + File.separator + pos + File.separator;
      FileUtils.createFolder(mitoDir);
      switch (METHOD){
         case 0:
            //TODO
            break;
         default:
            String[] mitosis_cmd = new String[]{PythonPipeline.getPythonDefaultPathInConda(), MaarsParameters.DEPS_DIR +
                    PythonPipeline.ANALYSING_SCRIPT_NAME, pathToRoot, Double.toString(calib.pixelWidth), pos,
                    Double.toString(calib.frameInterval)};
            ArrayList<String> cmds = new ArrayList<>();
            cmds.add(String.join(" ", mitosis_cmd));
            String bashPath = mitoDir + "pythonAnalysis.sh";
            FileUtils.writeScript(bashPath, cmds);
            IJ.log("Script saved. If it fails, you can still run it manually afterward.");
            PythonPipeline.runPythonScript(mitosis_cmd, mitoDir + "mitosisDetection_log.txt");
//      if (parameters.useDynamic()) {
//         findAbnormalCells(mitoDir, soc, getMitoticCellNbs(mitoDir));
//      }
      }
   }
}
