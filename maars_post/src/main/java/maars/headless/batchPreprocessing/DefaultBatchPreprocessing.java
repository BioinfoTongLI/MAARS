package maars.headless.batchPreprocessing;

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

@Plugin(type=BatchPreprocessing.class, name = BatchPreprocessing.NAME,
      attrs = { @Attr(name = "aliases", value = BatchPreprocessing.ALIASES) })
public class DefaultBatchPreprocessing extends AbstractOp implements BatchPreprocessing{

   @Parameter
   private String d;

   @Parameter
   private String configName;

   @Parameter
   private String suffix;
   
   @Parameter
   private String acq;
   
   @Parameter
   private boolean doGaussianBlur;
   
   @Parameter
   private boolean doAlignment;
   
   
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
         String processedImgFolder = fluoDir + "_processed_" + serieNbPos.get(serie);
         FileUtils.createFolder(processedImgFolder);
         ImgUtils.getChImages(usingChannels, imgPath, serie, processedImgFolder, doGaussianBlur, doAlignment);
         System.out.println("Images are preprocessing and saved into : \\n" + processedImgFolder);
      }
   }
}
