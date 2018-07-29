package maars.headless.batchFluoAnalysis;

import ij.IJ;
import ij.ImagePlus;
import ij.measure.ResultsTable;
import ij.plugin.Duplicator;

import loci.formats.FormatException;
import maars.agents.Cell;
import maars.agents.DefaultSetOfCells;
import maars.cellAnalysis.FluoAnalyzer;
import maars.cellAnalysis.PythonPipeline;
import maars.display.SOCVisualizer;
import maars.io.IOUtils;
import maars.main.MaarsParameters;
import maars.main.Maars_Interface;
import maars.utils.FileUtils;
import maars.utils.ImgUtils;
import org.apache.commons.io.filefilter.WildcardFileFilter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.*;
import java.util.concurrent.CopyOnWriteArrayList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicBoolean;

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
   public MaarsFluoAnalysis(ImagePlus[] impChs, String posName, MaarsParameters parameters, SOCVisualizer visualizer){
      processedChs_ = impChs;
      parameters_ = parameters;
      usingChannels = parameters_.getUsingChannels().split(",", -1);
      visualizer_ = visualizer;
      posName_ = posName;
      try {
         PrintStream ps = new PrintStream(parameters_.getSavingPath() + File.separator + "FluoAnalysis.LOG");
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
      AtomicBoolean stop = new AtomicBoolean(false);

      DefaultSetOfCells soc;
      String segAnaDir = FileUtils.convertPath(parameters_.getSavingPath()) + File.separator +
            parameters_.getSegmentationParameter(MaarsParameters.SEG_PREFIX) + Maars_Interface.SEGANALYSIS_SUFFIX;
      soc = new DefaultSetOfCells(posName_);
      String currentPosPrefix = segAnaDir + posName_ + File.separator;
      String currentZipPath = currentPosPrefix + "ROI.zip";
      if (FileUtils.exists(currentZipPath)) {
         // from Roi.zip initialize a set of cell
         soc.loadCells(currentZipPath);
         IJ.open(currentPosPrefix + "Results.csv");
         ResultsTable rt = ResultsTable.getResultsTable();
         ResultsTable.getResultsWindow().close(false);
         soc.addRoiMeasurementIntoCells(rt);
         // ----------------start acquisition and analysis --------//

         CopyOnWriteArrayList<Map<String, Future>> tasksSet = new CopyOnWriteArrayList<>();
         processStackedImg(processedChs_, parameters_, soc, visualizer_,
                 tasksSet, stop);
         Maars_Interface.waitAllTaskToFinish(tasksSet);
         if (!stop.get() && soc.size() != 0) {
            long startWriting = System.currentTimeMillis();
            ArrayList<String> arrayChannels = new ArrayList<>();
            Collections.addAll(arrayChannels, usingChannels);
            FileUtils.createFolder(parameters_.getSavingPath() + File.separator + parameters_.getFluoParameter(MaarsParameters.FLUO_PREFIX)
                  +Maars_Interface.FLUOANALYSIS_SUFFIX);
            IOUtils.saveAll(soc, processedChs_, parameters_.getSavingPath() + File.separator, parameters_.useDynamic(),
                  arrayChannels, posName_, parameters_.getFluoParameter(MaarsParameters.FLUO_PREFIX));
            IJ.log("It took " + (double) (System.currentTimeMillis() - startWriting) / 1000
                  + " sec for writing results");
            analyzeMitosisDynamic(soc, parameters_, processedChs_[0].getCalibration().pixelWidth);
         }
      soc.reset();
      System.gc();
      }
      System.setErr(curr_err);
      System.setOut(curr_out);
   }

   private void processStackedImg(ImagePlus[] processedChs, MaarsParameters parameters, DefaultSetOfCells soc,
                                       SOCVisualizer socVisualizer, CopyOnWriteArrayList<Map<String, Future>> tasksSet,
                                         AtomicBoolean stop) {
      ExecutorService es = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
      int totalChannel = Integer.parseInt(processedChs[0].getStringProperty("SizeC"));
      int totalFrame = Integer.parseInt(processedChs[0].getStringProperty("SizeT"));
      for (int i = 1; i <= totalFrame; i++) {
         Map<String, Future> chAnalysisTasks = new HashMap<>();
         for (int j = 0; j < totalChannel; j++) {
            String channel = usingChannels[j];
            IJ.log("Processing channel " + channel + "_" + i);
            ImagePlus slicedImg = duplicator.run(processedChs[j], i, i);
            Future future = es.submit(new FluoAnalyzer(slicedImg, slicedImg.getCalibration(),
                  soc, channel, Integer.parseInt(parameters.getChMaxNbSpot(channel)),
                  Double.parseDouble(parameters.getChSpotRaius(channel)),
                  Double.parseDouble(parameters.getChQuality(channel)), i, socVisualizer,
                  parameters.useDynamic()));
            chAnalysisTasks.put(channel, future);
         }
         tasksSet.add(chAnalysisTasks);
         if (stop.get()) {
            break;
         }
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
         cell.setAnaBOnsetFrame(anaBOnsetFrame);
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
      assert out != null;
      out.close();
      IJ.log("lagging detection finished");
   }

   static HashMap getMitoticCellNbs(String mitoDir) {
      return FileUtils.readTable(mitoDir + File.separator + "slopeChanges.csv");
   }

   public static void analyzeMitosisDynamic(DefaultSetOfCells soc, MaarsParameters parameters, double calib) {
      IJ.log("Start python analysis");
      String pos = soc.getPosLabel();
      String pathToRoot = parameters.getSavingPath() + File.separator;
      String mitoDir = pathToRoot + parameters.getFluoParameter(MaarsParameters.FLUO_PREFIX)+"_"+MITODIRNAME
            + File.separator + pos + File.separator;
      FileUtils.createFolder(mitoDir);

      String[] mitosis_cmd = new String[]{PythonPipeline.getPythonDefaultPathInConda(), MaarsParameters.DEPS_DIR +
            PythonPipeline.ANALYSING_SCRIPT_NAME, pathToRoot, Double.toString(calib), pos};
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
