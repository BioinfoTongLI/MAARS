package maars.main;

/*
  This class stores all the parameters that are needed to run MAARS

  @author Tong LI, mail: tongli.bioinfo@gmail.com
 * @version Nov 10, 2015
 */

import ij.IJ;
import maars.io.IOUtils;
import maars.utils.FileUtils;
import org.jdom2.Document;
import org.jdom2.Element;
import org.jdom2.JDOMException;
import org.jdom2.input.SAXBuilder;
import org.jdom2.output.Format;
import org.jdom2.output.XMLOutputter;

import java.awt.*;
import java.io.*;

/*
 * MaarsParameters reads a configuration file written as a XML,
 * then allows to access values thanks to all constants defined.
 *
 * SEGMENTATION_PARAMETERS
 *    |
 *    +-----> TOLERANCE
 *    +-----> PATH_TO_BF_ACQ_SETTING
 *    +-----> CHANNEL
 *    +-----> FRAME_NUMBER
 *    +-----> RANGE_SIZE_FOR_MOVIE
 *    +-----> STEP
 *    +-----> CELL_THICKNESS
 *    +-----> MINIMUM_CELL_AREA
 *    +-----> MEAN_GREY_VALUE
 *    +-----> SOLIDITY
 *    +-----> FILTER_MEAN_GREY_VALUE
 *    +-----> FILTER_SOLIDITY
 *    +-----> MAXIMUM_CELL_AREA
 *
 * PATH_TO_POSITION_LIST
 *   
 * FLUO_ANALYSIS_PARAMETERS
 *    |
 *    +-----> PATH_TO_FLUO_ACQ_SETTING
 *    +-----> USING
 *    +-----> DYNAMIC
 *    +-----> PROJECTED
 *    
 * GENERAL_ACQUISITION_PARAMETERS
 *    |
 *    +-----> BATCH_MODE
 *    +-----> SAVING_PATH
 *    +-----> CHANNEL_GROUP
 *    +-----> DEFAULT_CHANNEL_PARAMATERS
 * 					|
 * 					+-----> channel name
 * 								|
 * 								+-----> COLOR
 * 								+-----> EXPOSURE
 * 								+-----> SHUTTER
 * 							   +-----> SPOT_RADIUS
 *    				         +-----> MAXIMUM_NUMBER_OF_SPOT
 *    				         +-----> QUALITY
 * MITOSIS_DETECTION_PARAMETERS
 *    |
 *    +-----MINIMUM_DURATION
 *    +-----DETECTION_CHANNEL
 * @author Tong LI && Marie
 *
 */
public class MaarsParameters {

   public static final String STEP = "STEP";
   public static final String DYNAMIC = "DYNAMIC";
   public static final String PATH_TO_POSITION_LIST = "PATH_TO_POSITION_LIST";
   public static final String MINIMUM_CELL_AREA = "MINIMUM_CELL_AREA";
   public static final String MAXIMUM_CELL_AREA = "MAXIMUM_CELL_AREA";
   public static final String FILTER_MEAN_GREY_VALUE = "FILTER_MEAN_GREY_VALUE";
   public static final String MEAN_GREY_VALUE = "MEAN_GREY_VALUE";
   public static final String FILTER_SOLIDITY = "FILTER_SOLIDITY";
   public static final String SOLIDITY = "SOLIDITY";
   public static final String CHANNEL = "CHANNEL";
   public static final String PROJECTED = "PROJECTED";
   public static final String PATH_TO_BF_ACQ_SETTING = "PATH_TO_BF_ACQ_SETTING";
   public static final String PATH_TO_FLUO_ACQ_SETTING = "PATH_TO_FLUO_ACQ_SETTING";
   public static final String SEG_PREFIX = "SEG_PREFIX";
   public static final String FLUO_PREFIX = "FLUO_PREFIX";
   static final String CELL_THICKNESS = "CELL_THICKNESS";
   private static final String SEGMENTATION_PARAMETERS = "SEGMENTATION_PARAMETERS";
   private static final String FLUO_ANALYSIS_PARAMETERS = "FLUO_ANALYSIS_PARAMETERS";
   private static final String MITOSIS_DETECTION_PARAMETERS = "MITOSIS_DETECTION_PARAMETERS";
   private static final String MINIMUM_DURATION = "MINIMUM_DURATION";
   private static final String DETECTION_CHANNEL = "DETECTION_CHANNEL";
   private static final String SPOT_RADIUS = "SPOT_RADIUS";
   private static final String MAXIMUM_NUMBER_OF_SPOT = "MAXIMUM_NUMBER_OF_SPOT";
   private static final String QUALITY = "QUALITY";
   private static final String SAVING_PATH = "SAVING_PATH";
   private static final String USING = "USING";
   private static final String GENERAL_ACQUISITION_PARAMETERS = "GENERAL_ACQUISITION_PARAMETERS";
   private static final String DEFAULT_CHANNEL_PARAMATERS = "DEFAULT_CHANNEL_PARAMATERS";
   public static final String FOCUS= "FOCUS";
   public static final String DIRECTION = "DIRECTION";
   public static final String DEPS_DIR = IJ.getDirectory("macros");
   public static final String DEFAULT_CONFIG_NAME = "maars_config.xml";
   private Document doc;
   private Element root;

   /**
    * Constructor of Element need path to configuration file
    *
    * @param defaultParametersStream input stream containing xml file information
    */
   public MaarsParameters(InputStream defaultParametersStream) {

      final SAXBuilder sb = new SAXBuilder();
      try {
         try {
            doc = sb.build(defaultParametersStream);
         } catch (IOException e) {
            IOUtils.printErrorToIJLog(e);
         }
      } catch (JDOMException e) {
         IOUtils.printErrorToIJLog(e);
      }
      root = (Element) doc.getContent(0);
   }

   /**
    * empty
    */
   public MaarsParameters() {
   }

   /**
    * The few following colors are return as Color object : GREEN, CYAN, RED,
    * BLUE, WHITE NB : return GRAY if unknown color
    *
    * @param colorName name of the color
    * @return Color
    */
   public static Color getColor(String colorName) {
      if (colorName.equals("GREEN")) {
         return Color.GREEN;
      } else {
         if (colorName.equals("CYAN")) {
            return Color.CYAN;
         } else {
            if (colorName.equals("RED")) {
               return Color.RED;
            } else {
               if (colorName.equals("BLUE")) {
                  return Color.BLUE;
               } else {
                  if (colorName.equals("WHITE")) {
                     return Color.WHITE;
                  } else {
                     return Color.GRAY;
                  }
               }
            }
         }
      }
   }

   // Getter

   public static int getChNb(MaarsParameters parameters) {
      String channelsString = parameters.getUsingChannels();
      String[] arrayChannels = channelsString.split(",", -1);
      return arrayChannels.length;
   }

   public static String[] getChArray(MaarsParameters parameters) {
      String channelsString = parameters.getUsingChannels();
      return channelsString.split(",", -1);
   }

   /**
    * Write the parameters into the configuration file
    *
    * @throws IOException error than can not write xml file
    */
   public void save() throws IOException {
      doc.setContent(root);
      XMLOutputter xmlOutput = new XMLOutputter();
      xmlOutput.setFormat(Format.getPrettyFormat());
      xmlOutput.output(doc, new FileWriter(DEFAULT_CONFIG_NAME));
   }

   /**
    * * Write the parameters into the configuration file
    *
    * @param path path to save
    */
   public void save(String path) {
      doc.setContent(root);
      XMLOutputter xmlOutput = new XMLOutputter();
      xmlOutput.setFormat(Format.getPrettyFormat());
      try {
         xmlOutput.output(doc, new FileWriter(path + File.separator + DEFAULT_CONFIG_NAME));
      } catch (IOException e) {
         e.printStackTrace();
      }
   }

   /**
    * @return path to position list file
    */
   public String getPathToPositionList() {
      return root.getChildText(PATH_TO_POSITION_LIST);
   }

   /**
    * @param pathToPositionList path to positionlist file
    */
   public void setPathToPositionList(String pathToPositionList) {
      root.getChild(PATH_TO_POSITION_LIST).setText(pathToPositionList);
   }

   /**
    * @return analysis with dynamic or not
    */
   public boolean useDynamic() {
      return Boolean.parseBoolean(root.getChild(FLUO_ANALYSIS_PARAMETERS).getChildText(DYNAMIC));
   }

   /**
    * @return saving folder of MAARS output
    */
   public String getSavingPath() {
      return root.getChild(GENERAL_ACQUISITION_PARAMETERS).getChildText(SAVING_PATH);
   }

   /**
    * set saving path
    *
    * @param path : corresponding value of parameter
    */
   public void setSavingPath(String path) {
      root.getChild(GENERAL_ACQUISITION_PARAMETERS).getChild(SAVING_PATH).setText(path);
   }

   /**
    * @param parameter name of fluo parameter
    * @return time limit of fluorescence acquisition for one acquisition
    */
   public String getFluoParameter(final String parameter) {
      return root.getChild(FLUO_ANALYSIS_PARAMETERS).getChildText(parameter);
   }

   /**
    * @param parameter name of fluo parameter
    * @return time limit of fluorescence acquisition for one acquisition
    */
   public String getSegmentationParameter(final String parameter) {
      return root.getChild(SEGMENTATION_PARAMETERS).getChildText(parameter);
   }

   /**
    * @param ch: GFP, CFP, DAPI, TXRED
    * @return MAXIMUM_NUMBER_OF_SPOT of corresponding channel
    */
   public String getChMaxNbSpot(String ch) {
      return root.getChild(GENERAL_ACQUISITION_PARAMETERS).getChild(DEFAULT_CHANNEL_PARAMATERS).getChild(ch)
            .getChildText(MAXIMUM_NUMBER_OF_SPOT);
   }

   /**
    * @param ch: GFP, CFP, DAPI, TXRED
    * @return SPOT_RADIUS of corresponding channel
    */
   public String getChSpotRaius(String ch) {
      return root.getChild(GENERAL_ACQUISITION_PARAMETERS).getChild(DEFAULT_CHANNEL_PARAMATERS).getChild(ch).getChildText(SPOT_RADIUS);
   }

   /**
    * @param ch: GFP, CFP, DAPI, TXRED
    * @return SPOT_RADIUS of corresponding channel
    */
   public String getChQuality(String ch) {
      return root.getChild(GENERAL_ACQUISITION_PARAMETERS).getChild(DEFAULT_CHANNEL_PARAMATERS).getChild(ch).getChildText(QUALITY);
   }

   //////////// Setters

   /**
    * @return get channels used for fluo analysis
    */
   public String getUsingChannels() {
      return root.getChild(FLUO_ANALYSIS_PARAMETERS).getChildText(USING);
   }

   /**
    * set channels to USING channel
    *
    * @param channels channels that are being using for acquisitions
    */
   public void setUsingChannels(String channels) {
      root.getChild(FLUO_ANALYSIS_PARAMETERS).getChild(USING).setText(channels);
   }

   /**
    * @return the parameter minimum mitosis duration
    */
   public String getMinimumMitosisDuration() {
      return root.getChild(MITOSIS_DETECTION_PARAMETERS).getChildText(MINIMUM_DURATION);
   }

   /**
    * @param duration minimum mitosis duration (sec)
    */
   public void setMinimumMitosisDuration(String duration) {
      root.getChild(MITOSIS_DETECTION_PARAMETERS).getChild(MINIMUM_DURATION).setText(duration);
   }

   /**
    * @return the channel to use for mitosis detection
    */
   public String getDetectionChForMitosis() {
      return root.getChild(MITOSIS_DETECTION_PARAMETERS).getChildText(DETECTION_CHANNEL);
   }

   /**
    * @param chForMitosis the channel to detection mitosis
    */
   public void setDetectionChForMitosis(String chForMitosis) {
      root.getChild(MITOSIS_DETECTION_PARAMETERS).getChild(DETECTION_CHANNEL).setText(chForMitosis);
   }

   /**
    * @param projected whether or not project z stack
    */
   public void setProjected(String projected) {
      root.getChild(FLUO_ANALYSIS_PARAMETERS).getChild(PROJECTED).setText(projected);
   }

   /**
    * set segmentation parameter
    *
    * @param parameter : static final String of MaarsParameters
    * @param value     : corresponding value of parameter
    */
   public void setSegmentationParameter(String parameter, String value) {
      root.getChild(SEGMENTATION_PARAMETERS).getChild(parameter).setText(value);
   }

   /**
    * set fluo analysis parameter
    *
    * @param parameter : static final String of MaarsParameters
    * @param value     : corresponding value of parameter
    */
   public void setFluoParameter(String parameter, String value) {
      root.getChild(FLUO_ANALYSIS_PARAMETERS).getChild(parameter).setText(value);
   }

   /**
    * @param ch        GFP, CFP, DAPI, TXRED
    * @param maxNbSpot maximum number of spot for corresponding channel
    */
   public void setChMaxNbSpot(String ch, String maxNbSpot) {
      root.getChild(GENERAL_ACQUISITION_PARAMETERS).getChild(DEFAULT_CHANNEL_PARAMATERS).getChild(ch).getChild(MAXIMUM_NUMBER_OF_SPOT)
            .setText(maxNbSpot);
   }

   /**
    * @param ch         GFP, CFP, DAPI, TXRED
    * @param spotRaidus spotRaidus for corresponding channel
    */
   public void setChSpotRaius(String ch, String spotRaidus) {
      root.getChild(GENERAL_ACQUISITION_PARAMETERS).getChild(DEFAULT_CHANNEL_PARAMATERS).getChild(ch).getChild(SPOT_RADIUS)
            .setText(spotRaidus);
   }

   /**
    * @param ch      GFP, CFP, DAPI, TXRED
    * @param quality quality of spots for corresponding channel
    */
   public void setChQuality(String ch, String quality) {
      root.getChild(GENERAL_ACQUISITION_PARAMETERS).getChild(DEFAULT_CHANNEL_PARAMATERS).getChild(ch).getChild(QUALITY).setText(quality);
   }

   /**
    * @param segChannel bright field channel name
    */
   public void setSegChannel(String segChannel) {
      root.getChild(SEGMENTATION_PARAMETERS).getChild(CHANNEL).setText(segChannel);
   }

   /**
    * @param root the dataset of this class
    */
   private void setRoot(Element root) {
      this.root = root;
   }

   /**
    * duplicate this object
    *
    * @return the a duplicate version of this class
    */
   public MaarsParameters duplicate() {
      Element newRoot = root.clone();
      MaarsParameters newParams = new MaarsParameters();
      newParams.setRoot(newRoot);
      return newParams;
   }

   /**
    * @param path to config file
    */
   public static MaarsParameters fromFile(String path){
      InputStream inStream = null;
      if (FileUtils.exists(path)) {
         try {
            inStream = new FileInputStream(path);
         } catch (FileNotFoundException e) {
            e.printStackTrace();
         }
      } else {
         inStream = FileUtils.getInputStreamOfScript("maars_default_config.xml");
      }
      return new MaarsParameters(inStream);
   }
}