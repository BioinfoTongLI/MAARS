package maars.utils;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.ChannelSplitter;
import ij.plugin.RoiScaler;
import loci.formats.FormatException;
import loci.formats.ImageReader;
import loci.formats.MetadataTools;
import loci.formats.meta.IMetadata;
import loci.plugins.BF;
import loci.plugins.in.ImporterOptions;
import maars.segmentPombe.SegPombeParameters;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


/**
 * @author Tong LI, mail: tongli.bioinfo@gmail.com
 * @version Nov 25, 2015
 */
public class ImgUtils {

   /**
    * change unit of "cm" to "micron"
    *
    * @param img img with calibration to change
    */
   public static void unitCmToMicron(ImagePlus img) {
      img.getCalibration().setUnit("micron");
      img.getCalibration().pixelWidth = img.getCalibration().pixelWidth * 10000;
      img.getCalibration().pixelHeight = img.getCalibration().pixelHeight * 10000;
   }

   /**
    * re-calculate the position of ROI. make it adapt to the cropped image
    *
    * @param roi roi to process
    * @return processed ROI
    */
   public static Roi centerCroppedRoi(Roi roi) {
      int[] newXs = roi.getPolygon().xpoints;
      int[] newYs = roi.getPolygon().ypoints;
      int nbPoints = roi.getPolygon().npoints;
      for (int i = 0; i < nbPoints; i++) {
         newXs[i] = newXs[i] - (int) roi.getXBase();
         newYs[i] = newYs[i] - (int) roi.getYBase();
      }
      float[] newXsF = new float[nbPoints];
      float[] newYsF = new float[nbPoints];
      for (int i = 0; i < nbPoints; i++) {
         newXsF[i] = (float) newXs[i];
         newYsF[i] = (float) newYs[i];
      }
      return new PolygonRoi(newXsF, newYsF, Roi.POLYGON);
   }

   /**
    * Calculate rescale factor.
    *
    * @param cal1 calibration of imagePlus
    * @param cal2 calibration of imagePlus
    * @return factors of correction for x and y
    */
   public static double[] getRescaleFactor(Calibration cal1, Calibration cal2) {
      double[] factors = new double[2];
      if (cal1.equals(cal2)) {
         factors[0] = 1;
         factors[1] = 1;
      } else {
         factors[0] = cal1.pixelWidth / cal2.pixelWidth;
         factors[1] = cal1.pixelHeight / cal2.pixelHeight;
      }
      return factors;
   }

   /**
    * @param oldRoi  :roi to rescale
    * @param factors : double[] where first one is a factor to change width and
    *                second one is a factor to change height
    * @return rescaled ROI
    */
   @Deprecated
   public static Roi rescaleRoi(Roi oldRoi, double[] factors) {
      Roi roi = RoiScaler.scale(oldRoi, factors[0], factors[1], true);
      roi.setName("rescaledRoi");
      return roi;
   }

   /**
    * Method to shrink the image scale if it is bigger than supposed to the
    * maximum width an height that are in pixels
    *
    * @param maxHeight new height
    * @param maxWidth  new width
    */
   @Deprecated
   public static void changeScale(ImagePlus img, int maxWidth, int maxHeight, SegPombeParameters parameters) {
      int newWidth;
      int newHeight;
      int newMaxParticleSize;
      int newMinParticleSize;
      Calibration newCal = new Calibration();
      newCal.setUnit("micron");
      IJ.log("Before Width : " + String.valueOf(img.getWidth()) + ", Before Height : "
            + String.valueOf(img.getHeight()));
      if (img.getWidth() > maxWidth) {
         IJ.log("Image width is greater than maximum width allowed");

         newWidth = maxWidth;
         newHeight = img.getHeight() * maxWidth / img.getWidth();

         newMinParticleSize = (int) parameters.getMinParticleSize() * maxWidth / img.getWidth();
         newMaxParticleSize = (int) parameters.getMaxParticleSize() * maxWidth / img.getWidth();

         if (img.getCalibration().scaled()) {

            newCal.pixelWidth = parameters.getScale(SegPombeParameters.WIDTH) * img.getWidth() / maxWidth;
            newCal.pixelHeight = parameters.getScale(SegPombeParameters.HEIGHT) * img.getWidth() / maxWidth;
            newCal.pixelDepth = parameters.getScale(SegPombeParameters.DEPTH);
         }

         IJ.log("New values : w = " + newWidth + " h = " + newHeight);
         if (newHeight > maxHeight) {
            IJ.log("New height is still greater than maximum height allowed");
            newHeight = maxHeight;
            newWidth = img.getWidth() * maxHeight / img.getHeight();

            newMinParticleSize = (int) parameters.getMinParticleSize() * maxHeight / img.getHeight();
            newMaxParticleSize = (int) parameters.getMaxParticleSize() * maxHeight / img.getHeight();

            if (img.getCalibration().scaled()) {
               newCal.pixelWidth = parameters.getScale(SegPombeParameters.WIDTH) * img.getHeight() / maxHeight;
               newCal.pixelHeight = parameters.getScale(SegPombeParameters.HEIGHT) * img.getHeight() / maxHeight;
               newCal.pixelDepth = parameters.getScale(SegPombeParameters.DEPTH);
            }

            IJ.log("New values : w = " + newWidth + " h = " + newHeight);

         }

         rescale(newWidth, newHeight, newMinParticleSize, newMaxParticleSize, newCal, parameters, img);
      } else {
         if (img.getHeight() > maxHeight) {
            IJ.log("Image height is greater than maximum width allowed");

            newHeight = maxHeight;
            newWidth = img.getWidth() * maxHeight / img.getHeight();

            newMinParticleSize = (int) parameters.getMinParticleSize() * maxHeight / img.getHeight();
            newMaxParticleSize = (int) parameters.getMaxParticleSize() * maxHeight / img.getHeight();

            if (img.getCalibration().scaled()) {
               newCal.pixelWidth = parameters.getScale(SegPombeParameters.WIDTH) * img.getHeight() / maxHeight;
               newCal.pixelHeight = parameters.getScale(SegPombeParameters.HEIGHT) * img.getHeight() / maxHeight;
               newCal.pixelDepth = parameters.getScale(SegPombeParameters.DEPTH);
            }

            IJ.log("New values : w = " + newWidth + " h = " + newHeight);

            if (newWidth > maxWidth) {
               IJ.log("New Width is still greater than maximum height allowed");

               newWidth = maxWidth;
               newHeight = img.getHeight() * maxWidth / img.getWidth();

               if (img.getCalibration().scaled()) {
                  newCal.pixelWidth = parameters.getScale(SegPombeParameters.WIDTH) * img.getWidth() / maxWidth;
                  newCal.pixelHeight = parameters.getScale(SegPombeParameters.HEIGHT) * img.getWidth() / maxWidth;
                  newCal.pixelDepth = parameters.getScale(SegPombeParameters.DEPTH);
               }

               IJ.log("New values : w = " + newWidth + " h = " + newHeight);

               newMinParticleSize = (int) parameters.getMinParticleSize() * maxWidth / img.getWidth();
               newMaxParticleSize = (int) parameters.getMaxParticleSize() * maxWidth / img.getWidth();

            }

            rescale(newWidth, newHeight, newMinParticleSize, newMaxParticleSize, newCal, parameters, img);
         }
      }
   }

   /**
    * Method to change size and scale of the image to analyze : need to compute
    * parameters before so it can be coherent
    */
   @Deprecated
   public static void rescale(int newWidth, int newHeight, int newMinParticleSize, int newMaxParticleSize,
                        Calibration newCal, SegPombeParameters parameters, ImagePlus img) {

      parameters.setMinParticleSize(newMinParticleSize);
      parameters.setMaxParticleSize(newMaxParticleSize);

      IJ.log("min area = " + newMinParticleSize + " max area = " + newMaxParticleSize);

      ImageStack newImgStack = new ImageStack(newWidth, newHeight);

      for (int slice = 0; slice < img.getNSlices(); slice++) {
         img.setZ(slice);
         newImgStack.addSlice(img.getProcessor().resize(newWidth, newHeight));
      }

      ImagePlus newImagePlus = new ImagePlus("rescaled_" + img.getTitle(), newImgStack);
      newImagePlus.setCalibration(newCal);

      parameters.setImageToAnalyze(newImagePlus);

      checkImgUnitsAndScale(img, parameters);
   }

   /**
    * Method to check if the image is scaled and if the unit matches 'micron'
    */
   public static Boolean checkImgUnitsAndScale(ImagePlus img, SegPombeParameters parameters) {
      IJ.log("Check if image is scaled");
      if (img.getCalibration().scaled()) {

         if (img.getCalibration().getUnit().equals("cm")) {
            img.getCalibration().setUnit("micron");
            img.getCalibration().pixelWidth = img.getCalibration().pixelWidth * 10000;
            img.getCalibration().pixelHeight = img.getCalibration().pixelHeight * 10000;
         }

         IJ.log("Check if unit of calibration is micron");
         if (img.getCalibration().getUnit().equals("micron")
               || img.getCalibration().getUnit().equals("µm")) {
            double[] scale = new double[3];
            scale[SegPombeParameters.WIDTH] = img.getCalibration().pixelWidth;
            scale[SegPombeParameters.HEIGHT] = img.getCalibration().pixelHeight;
            scale[SegPombeParameters.DEPTH] = img.getCalibration().pixelDepth;
            parameters.setScale(scale);
            IJ.log("Get and set calibration as scale");
            return true;
         } else {
            IJ.error("Wrong scale unit",
                  "The scale of your image must be in microns.\nTo change it you can go to Properties ...");
         }
      } else {
         IJ.error("No scale",
               "You must set a scale to your image\nif you want to enter measures in micron.\nYou can go to Properties ...");
      }
      return false;
   }

   /**
    * Allows to convert width or height or depth of a size in microns to a size
    * in pixels. To choose if you want to convert width or height you can use
    * CellsBoundaries constants : WIDTH and HEIGHT and DEPTH. Return an int
    * width or height in pixels.
    */
   public static int convertMicronToPixel(double micronSize, int widthOrHeightOrDepth, SegPombeParameters parameters) {
      return (int) Math.round(micronSize / parameters.getScale(widthOrHeightOrDepth));
   }

   /**
    * Allows to convert width or height or depth of a size in pixels to a size
    * in microns. To choose if you want to convert width or height you can use
    * CellsBoundaries constants : WIDTH and HEIGHT and DEPTH.
    *
    * @param pixelSize            a size in pixels
    * @param widthOrHeightOrDepth Dimension
    * @return Return an double width or height in microns.
    */
   public static double convertPixelToMicron(int pixelSize, int widthOrHeightOrDepth, SegPombeParameters parameters) {
      return parameters.getScale(widthOrHeightOrDepth) * pixelSize;
   }

   public static ImagePlus lociImport(String tiffPath, int serie_number) throws IOException, FormatException {
      return lociImport(tiffPath, serie_number, true);
   }
   
   public static ImagePlus lociImport(String tiffPath){
      IJ.run("Bio-Formats Importer", "open="+tiffPath+" color_mode=Default open_files view=Hyperstack stack_order=XYCZT");
      ImagePlus imp = IJ.getImage();
      imp.hide();
      return imp;
   }

   public static ImagePlus lociImport(String tiffPath, int serie_number, Boolean virtual) throws IOException, FormatException {
      ImporterOptions options = new ImporterOptions();
      options.setVirtual(virtual);
      options.setGroupFiles(false);
      options.setMustGroup(false);
      options.setId(tiffPath);
      options.setStackFormat("Hyperstack");
      options.clearSeries();
      options.setSeriesOn(serie_number, true);
      return BF.openImagePlus(options)[0];
   }

   public static HashMap<Integer, String> populateSeriesImgNames(String pathToTiffFile){
      Pattern pattern = Pattern.compile(".*_(.*)");

      HashMap<Integer, String> seriesImgNames = new HashMap<>();
      int seriesCount = 0;
      try (ImageReader reader = new ImageReader()) {
         IMetadata omexmlMetadata = MetadataTools.createOMEXMLMetadata();
         reader.setMetadataStore(omexmlMetadata);
         try {
            reader.setId(pathToTiffFile);
         } catch (FormatException | IOException e) {
            e.printStackTrace();
         }
         seriesCount = reader.getSeriesCount();
         for (int i = 0; i < seriesCount; i++) {
            reader.setSeries(i);
            String name = omexmlMetadata.getImageName(i); // this is the image name stored in the file
            String pos = null;
            if (name.equals("")){
               pos = getPosNameFromFileName(pathToTiffFile);
            }else{
               Matcher matcher = pattern.matcher(name);
               if (matcher.find()) {
                  pos = matcher.group(1);
               }else{
                  new FormatException("Can not load series names from metadata");
               }
            }
            assert pos !=null;
            seriesImgNames.put(i, pos);
         }
      } catch (IOException e) {
         e.printStackTrace();
      }
      assert seriesCount !=0 ;
      IJ.log(seriesCount + " series registered");
      return seriesImgNames;
   }

   public static String getPosNameFromFileName(String filePath){
      String sep;
      if (IJ.isWindows()){
         sep = "\\\\";
      }else{
         sep = File.separator;
      }
      String[] splits = filePath.split(sep);
      String fileName = splits[splits.length-1];
      Pattern pattern = Pattern.compile(".*_(.*).ome.tif*");
      Matcher matcher = pattern.matcher(fileName);
      String pos = null;
      if (matcher.find()) {
         pos = matcher.group(1);
      }
      assert pos != null;
      return pos;
   }
   
   public static HashMap<Integer, String> populateSeriesImgNames(File[] fileNames){
      HashMap<Integer, String> seriesImgNames = new HashMap<>();
      int counter = 0;
      for (File f : fileNames){
         seriesImgNames.put(counter, f.getAbsolutePath());
         counter++;
      }
      return seriesImgNames;
   }

   public static ImagePlus[] alignChannels(ImagePlus[] imps, String savingPath, String[] arrayChannels) {
      int nbOfImgs = imps.length;
      ImagePlus[] projectedImps = new ImagePlus[nbOfImgs];
      for (int i = 0 ;  i < nbOfImgs; i++) {
         IJ.log("Z-projecting...");
         IJ.run(imps[i], "Z Project...", "projection=[Max Intensity] all");
         projectedImps[i] = IJ.getImage();
         projectedImps[i].hide();
      }

      int imgWidth = projectedImps[1].getWidth();
      int imgHeight = projectedImps[1].getWidth();

      int cropMagnitude = 2;
      int newW = Math.floorDiv(imgWidth, cropMagnitude);
      int newH = Math.floorDiv(imgHeight, cropMagnitude);
      int selectionOriginX = Math.floorDiv(imgWidth, 2) - newW / 2;
      int selectionOriginY = Math.floorDiv(imgHeight, 2) - newH / 2;

      String tranfoFileName = "transformation_from_" + arrayChannels[1] + ".txt";
      IJ.run(projectedImps[1], "Align slices in stack...", "method=5 windowsizex=" +
              newW + " windowsizey=" + newH + " x0=" + selectionOriginX + " y0=" + selectionOriginY +
              " swindow=0 subpixel=false itpmethod=0 ref.slice=1 show=true");
      projectedImps[1].killRoi();
      IJ.saveAs("Results", savingPath + File.separator + tranfoFileName);
      ResultsTable rt = ResultsTable.getResultsTable();
      for (int i = 0; i < projectedImps[0].getNFrames()-1; i++){
         projectedImps[0].setSlice((int) rt.getValue("Slice", i));
         projectedImps[0].getProcessor().translate(
                 (int) rt.getValue("dX", i),
                 (int) rt.getValue("dY", i));
      }
      rt.reset();
      IJ.run("Clear Results", "");
      for (ImagePlus im : projectedImps) {
         im.hide();
      }
      return projectedImps;
   }

   public static ImagePlus[] blurChannels(ImagePlus[] chs){
      IJ.log("Denoising the image with 3D-blurring...");
      for (ImagePlus im : chs){
         IJ.run(im, "Gaussian Blur 3D...", "x=0.85 y=0.85 z=1.7");
      }
      return chs;
   }

   public static ImagePlus[] preprocessChs(ImagePlus concatenatedFluoImgs, String[] usingChannels,
                                           String processedImgFolder, boolean gaussian_blur, boolean align,
                                           boolean saveIntermedia){
      int totalChannel = Integer.parseInt(concatenatedFluoImgs.getStringProperty("SizeC"));
      double interval = concatenatedFluoImgs.getCalibration().frameInterval;

      ImagePlus[] processedChs = ChannelSplitter.split(concatenatedFluoImgs);
      concatenatedFluoImgs.close();
      if (gaussian_blur) {
         processedChs = ImgUtils.blurChannels(processedChs);
         for (int i = 0; i < totalChannel; i++) {
            processedChs[i].getCalibration().frameInterval = interval;
            if (saveIntermedia) {
               IJ.saveAsTiff(processedChs[i], processedImgFolder + File.separator + usingChannels[i]
                       + "_denoised");
            }
         }
      }
      if (align) {
         processedChs = ImgUtils.alignChannels(processedChs, processedImgFolder, usingChannels);
         for (int i = 0; i < totalChannel; i++) {
            processedChs[i].getCalibration().frameInterval = interval;
            if (saveIntermedia) {
               IJ.saveAsTiff(processedChs[i], processedImgFolder + File.separator + usingChannels[i] + "_aligned");
            }
         }
      }
      for (int i = 0; i < totalChannel; i++) {
         processedChs[i].getCalibration().frameInterval = interval;
         IJ.saveAsTiff(processedChs[i], processedImgFolder + File.separator + usingChannels[i] + "_final");
      }
      System.gc();
      return processedChs;
   }

   public static ImagePlus[] getChImages(String[] usingChannels, String imgPath, int serie,
                                         String processedImgFolder, boolean doGaussianBlue, boolean doAlignment) {
      ImagePlus[] processedChs;
      String firstChImgPath = processedImgFolder + File.separator + usingChannels[0] + "_final.tif";
      String secChImgPath = processedImgFolder + File.separator + usingChannels[1] + "_final.tif";
      if (FileUtils.exists(firstChImgPath) &&
              FileUtils.exists(secChImgPath)) {
         IJ.log("Images are previously preprocessed... if you want to redo analysis, please delete them.");
         processedChs = new ImagePlus[]{
                 IJ.openImage(firstChImgPath),
                 IJ.openImage(secChImgPath)
         };
      } else {
         ImagePlus concatenatedFluoImgs = null;
         try {
            concatenatedFluoImgs = ImgUtils.lociImport(imgPath, serie);
         } catch (IOException | FormatException e) {
            e.printStackTrace();
            IJ.error("Can not load fluo images...");
         }
         assert concatenatedFluoImgs != null;
         processedChs = ImgUtils.preprocessChs(concatenatedFluoImgs, usingChannels,
                 processedImgFolder, doGaussianBlue, doAlignment, true);
      }
      return processedChs;
   }
}
