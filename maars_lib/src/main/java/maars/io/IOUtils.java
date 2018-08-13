package maars.io;

import fiji.plugin.trackmate.Model;
import fiji.plugin.trackmate.io.TmXmlWriter;
import ij.IJ;
import ij.ImagePlus;
import maars.agents.Cell;
import maars.agents.DefaultSetOfCells;
import maars.main.Maars_Interface;
import maars.utils.FileUtils;


import java.io.*;
import java.util.concurrent.CopyOnWriteArrayList;

/**
 * Created by tongli on 27/12/2016.
 */
public class IOUtils {
   public static void printErrorToIJLog(Exception exception) {
      StringWriter sw = new StringWriter();
      PrintWriter ps = new PrintWriter(sw);
      exception.printStackTrace(ps);
      IJ.error(sw.toString());
      try {
         sw.close();
      } catch (IOException e1) {
         e1.printStackTrace();
      }
      ps.close();
   }

   public static void serializeSoc(String pathToFluoDir, DefaultSetOfCells soc) {
      File f = new File(pathToFluoDir + "SetOfCell.serialize");
      ObjectOutputStream objOut = null;
      try {
         objOut = new ObjectOutputStream(new BufferedOutputStream(
               new FileOutputStream(f)));
         objOut.writeObject(soc);
         objOut.flush();

         IJ.log("Set of cel object is serialized.");
      } catch (IOException i) {
         IJ.log((i.getMessage()));
      } finally {
         if (objOut != null) {
            try {
               objOut.close();
            } catch (IOException e) {
               IOUtils.printErrorToIJLog(e);
            }
         }
      }
   }

   public static void saveAll(int method, DefaultSetOfCells soc, ImagePlus[] processStack, String pathToDir,
                              Boolean useDynamic, String[] arrayChannels, String posNb, String prefix) {
      long startWriting = System.currentTimeMillis();
      String analysisOutputFolder = pathToDir + prefix +Maars_Interface.FLUOANALYSIS_SUFFIX;
      FileUtils.createFolder(analysisOutputFolder);
      IJ.log("Saving information of each cell on disk");
      String dest = analysisOutputFolder + posNb + File.separator;
      FileUtils.createFolder(dest);
//        TODO
      CopyOnWriteArrayList<Integer> cellIndex = soc.getPotentialMitosisCell();
      MAARSImgSaver imgSaver = new MAARSImgSaver(dest, arrayChannels);
      MAARSSpotsSaver spotSaver = new MAARSSpotsSaver(dest);
      MAARSGeometrySaver geoSaver = new MAARSGeometrySaver(dest);
      for (int i : cellIndex) {
         Cell cell = soc.getCell(i);
//        for (Cell cell : soc){
         switch (method){
            case 0:
               String tracksDir = dest + "tracks" + File.separator;
               FileUtils.createFolder(tracksDir);
               for (int j = 1; j <= arrayChannels.length; j++) {
                  File outputFile = new File(tracksDir + i + "_" + arrayChannels[j - 1] + ".xml");
                  TmXmlWriter writer = new TmXmlWriter(outputFile);
                  Model currentModel = cell.getModel(arrayChannels[j - 1]);
                  if (currentModel != null) {
                     writer.appendModel(currentModel);
                     try {
                        writer.writeToFile();
                     } catch (IOException e) {
                        e.printStackTrace();
                     }
                  }
               }
               break;
            default:
               geoSaver.save(cell);
               spotSaver.save(cell);
         }
      }
      if (!imgSaver.exists) {
         for (int i : cellIndex){
            Cell cell = soc.getCell(i);
            imgSaver.saveImgs(processStack, cell.getCellShapeRoi(), i);
         }
      }
      if (useDynamic) {
         IOUtils.serializeSoc(dest, soc);
      }
      IJ.log("It took " + (double) (System.currentTimeMillis() - startWriting) / 1000
              + " sec for writing results");
   }
}
