package maars.io;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.plugin.Duplicator;
import maars.utils.FileUtils;
import maars.utils.ImgUtils;

import java.io.File;

public class MAARSImgSaver {
   private static final String croppedImgs = "croppedImgs";
   private String croppedImgDir;
   private Duplicator duplicator = new Duplicator();
   private String[] chNames;
   public boolean exists;

   MAARSImgSaver(String pathToFluoDir, String[] chNames) {
      this.croppedImgDir = pathToFluoDir + File.separator + croppedImgs + File.separator;
      exists = FileUtils.createFolder(croppedImgDir);
      this.chNames = chNames;
   }

   public void saveImgs(ImagePlus[] entireField, Roi roi, int cellNb){
      ImagePlus croppedImg;
      for (int j = 1; j <= chNames.length; j++) {
         entireField[j-1].setRoi(roi);
         croppedImg = duplicator.run(entireField[j-1], 1, entireField[j-1].getNFrames());
         IJ.run(croppedImg, "Grays", "");
         croppedImg.setRoi(ImgUtils.centerCroppedRoi(roi));
         IJ.run(croppedImg, "Enhance Contrast", "saturated=0.35");
         IJ.run(croppedImg, "Clear Outside", "stack");
         String pathToCroppedImg = croppedImgDir + String.valueOf(cellNb) + "_" + chNames[j-1] + ".tif";
         IJ.saveAsTiff(croppedImg, pathToCroppedImg);
      }
   }


//   public void exportChannelBtf(Boolean splitChannel, Set<String> arrayChannels) {
//      if (splitChannel) {
//         for (String channel : arrayChannels) {
//            final String btfPath = pathToFluoDir + File.separator + channel + ".ome.btf";
//            if (!FileUtils.exists(btfPath)) {
//               ImageStack currentStack = new ImageStack(mergedFullFieldImg.getWidth(),
//                       mergedFullFieldImg.getHeight());
//               for (int j = 1; j <= mergedFullFieldImg.getImageStack().getSize(); j++) {
//                  if (mergedFullFieldImg.getStack().getSliceLabel(j).equals(channel)) {
//                     currentStack.addSlice(mergedFullFieldImg.getStack().getProcessor(j));
//                  }
//               }
//               final String macroOpts = "outfile=[" + btfPath
//                       + "] splitz=[0] splitc=[0] splitt=[0] compression=[Uncompressed]";
//               ImagePlus currentChImp = new ImagePlus(channel, currentStack);
//               currentChImp.setCalibration(mergedFullFieldImg.getCalibration());
//               lociExporter.setup(macroOpts, currentChImp);
//               lociExporter.run(null);
//            }
//         }
//      } else {
//         String btfPath = pathToFluoDir + "merged.ome.btf";
//         if (!FileUtils.exists(btfPath)) {
//            final String macroOpts = "outfile=[" + btfPath
//                    + "] splitz=[0] splitc=[0] splitt=[0] compression=[Uncompressed]";
//            lociExporter.setup(macroOpts, mergedFullFieldImg);
//            lociExporter.run(null);
//         }
//      }
//   }
}
