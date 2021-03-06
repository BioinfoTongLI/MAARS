package maars.headless.batchFluoAnalysis;

import ij.ImagePlus;
import maars.agents.DefaultSetOfCells;
import maars.display.SOCVisualizer;
import maars.main.MaarsParameters;
import maars.main.Maars_Interface;
import maars.utils.FileUtils;
import maars.utils.ImgUtils;
import net.imagej.ops.AbstractOp;
import org.apache.commons.io.filefilter.WildcardFileFilter;
import org.scijava.plugin.Attr;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import java.io.File;
import java.io.FilenameFilter;
import java.util.Map;
import java.util.Objects;

@Plugin(type=BatchFluoAnalysis.class, name = BatchFluoAnalysis.NAME,
      attrs = { @Attr(name = "aliases", value = BatchFluoAnalysis.ALIASES) })
public class DefaultBatchFluoAnalysis extends AbstractOp implements BatchFluoAnalysis{

   @Parameter
   private String d;

   @Parameter
   private String configName;

   @Parameter
   private String suffix;
   
   @Parameter
   private String acq;
   
   @Parameter
   private int met;

   @Override
   public void run() {
      Maars_Interface.copyDeps();
      System.out.println(d);
      MaarsParameters parameter = MaarsParameters.fromFile(d + File.separator  + configName);
      parameter.setSavingPath(d);
      parameter.setSegmentationParameter(MaarsParameters.SEG_PREFIX, "BF_" + acq);
      parameter.setFluoParameter(MaarsParameters.FLUO_PREFIX, "FLUO_" + acq);
      parameter.save(d);
      String[] usingChannels = parameter.getUsingChannels().split(",", -1);

      String fluoDir = FileUtils.convertPath(parameter.getSavingPath()) + File.separator +
              parameter.getFluoParameter(MaarsParameters.FLUO_PREFIX);
      String imgPath = Objects.requireNonNull(new File(fluoDir).listFiles(
              (FilenameFilter) new WildcardFileFilter("*." + suffix)))[0].getAbsolutePath();
      Map<Integer, String> serieNbPos = ImgUtils.populateSeriesImgNames(imgPath);


      for (int serie : serieNbPos.keySet()) {
         SOCVisualizer visualizer = new SOCVisualizer(d + "|" + serieNbPos.get(serie), usingChannels);
         visualizer.setVisible(true);
         DefaultSetOfCells soc = new DefaultSetOfCells(serieNbPos.get(serie));
         String processedImgFolder = fluoDir + "_processed_" + serieNbPos.get(serie);
         FileUtils.createFolder(processedImgFolder);
   
         ImagePlus[] processedChs = ImgUtils.getChImages(usingChannels, imgPath, serie, processedImgFolder,
               true, true);
         MaarsFluoAnalysis.METHOD = met;
         Thread th = new Thread(new MaarsFluoAnalysis(soc, processedChs, serieNbPos.get(serie),
                 parameter, visualizer));
         th.start();
         try {
            th.join();
         } catch (InterruptedException e) {
            e.printStackTrace();
         }
//         TODO to find a way to save these intermedia data
//         soc.reset();
//         visualizer.clear();
      }
   }
}
