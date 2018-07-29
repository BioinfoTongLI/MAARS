package maars.headless.batchFluoAnalysis;

import ij.ImagePlus;
import loci.formats.FormatException;
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
import java.io.IOException;
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

   @Override
   public void run() {
      Maars_Interface.copyDeps();
      System.out.println(d);
      MaarsParameters parameter = MaarsParameters.fromFile(d + File.separator  + configName);
      parameter.setSavingPath(d);
      parameter.save(d);
      String[] usingChannels = parameter.getUsingChannels().split(",", -1);

      String fluoDir = FileUtils.convertPath(parameter.getSavingPath()) + File.separator +
              parameter.getFluoParameter(MaarsParameters.FLUO_PREFIX);
      String imgPath = Objects.requireNonNull(new File(fluoDir).listFiles(
              (FilenameFilter) new WildcardFileFilter("*." + suffix)))[0].getAbsolutePath();
      Map<Integer, String> serieNbPos = ImgUtils.populateSeriesImgNames(imgPath);
      SOCVisualizer visualizer;

      for (int serie : serieNbPos.keySet()) {
         visualizer = new SOCVisualizer(usingChannels);
         visualizer.setVisible(true);
         String processedImgFolder = fluoDir + "_processed_" + serieNbPos.get(serie);
         FileUtils.createFolder(processedImgFolder);

         ImagePlus concatenatedFluoImgs = null;
         try {
            concatenatedFluoImgs = ImgUtils.lociImport(imgPath, serie);
         } catch (IOException | FormatException e) {
            e.printStackTrace();
         }
         assert concatenatedFluoImgs != null;

         ImagePlus[] processedChs = ImgUtils.preprocessChs(concatenatedFluoImgs, usingChannels,
                 processedImgFolder, true, true,
                 Double.parseDouble(parameter.getFluoParameter(MaarsParameters.TIME_INTERVAL)) / 1000);

         Thread th = new Thread(
                 new MaarsFluoAnalysis(processedChs, serieNbPos.get(serie),
                 parameter, visualizer));
         th.start();
         try {
            th.join();
         } catch (InterruptedException e) {
            e.printStackTrace();
         }
         visualizer.clear();
         visualizer.setVisible(true);
      }
   }
}
