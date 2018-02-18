import ij.ImagePlus;

import loci.formats.FormatException;
import loci.formats.IFormatReader;
import loci.formats.in.OMETiffReader;
import loci.plugins.BF;
import loci.plugins.in.ImporterOptions;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;

@RunWith(Parameterized.class)
public class  BFSegTest {
   private static Logger logger = LoggerFactory.getLogger(BFSegTest.class);
   @Parameterized.Parameters
   public static Collection<Object[]> prepareFiles() throws IOException, FormatException {
      String root = System.getProperty("user.dir") + "/src/main/resources/";

      ImporterOptions options = new ImporterOptions();
      options.setVirtual(true);
      options.setId(root + "BF.tif");
      options.setStackFormat("Hyperstack");
      ImagePlus[] bfimg = BF.openImagePlus(options);

      options.setId(root + "FLUO.tif");
      ImagePlus[] fluoimg = BF.openImagePlus(options);
      return Arrays.asList(new Object[][] {{bfimg[0], fluoimg[0]}});
   }
   @Parameterized.Parameter
   public ImagePlus bfimg;

   @Parameterized.Parameter(1)
   public ImagePlus fluoimg;

//   @Test
//   public void show() throws ImgIOException, IOException, FormatException, io.scif.FormatException {
//      String path = "/Volumes/Macintosh/curioData/screening/20_10_17_2/BF_1/BF_1_MMStack_wt.ome.tif";
//      File oneExample = new File(path);
//      final SCIFIO scifio = new SCIFIO();
//
//      SCIFIOConfig config = new SCIFIOConfig();
//      config = config.imgOpenerSetOpenAllImages(true);
//      config = config.imgOpenerSetImgModes(SCIFIOConfig.ImgMode.AUTO);
//      config = config.checkerSetOpen(true);
//
//      final Format format =
//            scifio.format().getFormat(oneExample.getAbsolutePath(), config);
//      final Parser parser = format.createParser();
//      final Metadata meta = parser.parse(oneExample);
//      BioFormatsFormat.Reader bfreader = new BioFormatsFormat.Reader();
//      bfreader.setMetadata(meta);

//      Reader reader = scifio.initializer().initializeReader(path, config);
//      List<Dataset> datasets =  scifio.datasetIO().openAll();


//      ImgPlus<FloatType>  imp = dataset.typedImg(new FloatType());



//      logger.info(me + "");
//      for (int i = 0; i < meta.getImageCount(); i++) {
//         System.out.println("image #" + i + " dimensions:");
//         for (int a = 0; a < meta.get(i).getAxes().size(); a++) {
//            final AxisType axisType = meta.get(i).getAxis(a).type();
//            final long axisLength = meta.get(i).getAxisLength(a);
//            System.out.println("\t" + axisLength + " : " + axisType);
//         }
//      }

//      IFormatReader reader = new OMETiffReader();

//      bfreader.openPlane(0,0,config).getLengths();
//      logger.info(bfreader.getCurrentFile());
//      bfFormat.addReader(IFormatReader.class);
//      Reader reader = bfFormat.createReader();
//      ImageReader reader = bfFormat.createImageReader();
//      logger.info(reader.setSource(););

//      reader.setSource(oneExample, config);
//      reader.setId(oneExample);
//      logger.info(reader.getFormatName() + "");

//      Reader reader = new TIFFFormat.Reader<>();
//      reader.setContext(new Context());
//      reader.setSource("/Volumes/Macintosh/curioData/screening/22_11_17/BF_1/BF_1_MMStack_wt_1.ome.tif", config);
//      logger.info(reader.getImageCount() + "");
//      ImgOpener imgOpener = new ImgOpener();
//      List<SCIFIOImgPlus<FloatType>> imgs = imgOpener.openImgs(reader, new FloatType(), config);
//      ImageJFunctions.show( imgs.get(0) );
//      for (int i = 0; i < imgs.size(); i++) {
//         final SCIFIOImgPlus<?> testImg = imgs.get(i);
//         assertEquals(2560, testImg.dimension(0));
//         assertEquals(2160, testImg.dimension(1));
//         assertEquals(34, testImg.dimension(2));
//      }

//      ImageJFunctions.show( bfimg );
//      ImageJFunctions.show( fluoimg );

//      OpService ops = new DefaultOpService();
//      FinalInterval finalInterval = new FinalInterval(0,0,5,5);
//
//      RandomAccessibleInterval<UnsignedByteType> view = Views.interval(bfimg.getImg(), finalInterval);
//      ImageJFunctions.show( view );
//      logger.info((int) ops.run(Ops.Math.Add.class, 1,2) + "");
//      int x = (int) ops.run("math.add", 1,1);
//      System.out.println(x);
//      RandomAccessibleInterval
//      logger.info(x + "");
//      ImageJFunctions.show( img. );
//      try {
//         Thread.sleep(100000);
//      } catch (InterruptedException e) {
//         e.printStackTrace();
//      }
//   }

   @Test
   public void bfloadTest() throws IOException, FormatException {
//      String oneExample = "/Volumes/Macintosh/curioData/screening/20_10_17_2/BF_1";
//      File[] listFiles = new File(oneExample).listFiles((FilenameFilter) new WildcardFileFilter("*.tif"));
//      logger.info(listFiles[0] + "");
      String oneExample = "/Volumes/Macintosh/curioData/screening/20_10_17_2/BF_1/BF_1_MMStack_wt.ome.tif";
      IFormatReader reader = new OMETiffReader();
      reader.setId(oneExample);
      for (String s :reader.getUsedFiles()) logger.info(s);
//      logger.info(reader.getSeriesCount() + "");
//      assertEquals(4, reader.getSeriesCount());
//
//      for (int i = 0; i < 4; i ++) {
//         reader.setSeries(i);
//         String[] files = reader.getSeriesUsedFiles();
//         for (String s : files) {
//            logger.info(s);
//         }
//      }
//
//      ImporterOptions options = new ImporterOptions();
//      options.setVirtual(true);
//      options.setId(oneExample);
//      options.setStackFormat("Hyperstack");
//      options.clearSeries();
//      options.setSeriesOn(3, true);
//      ImagePlus[] imps = BF.openImagePlus(options);
//      imps[0].show();
//      try {
//         Thread.sleep(100000);
//      } catch (InterruptedException e) {
//         e.printStackTrace();
//      }
   }

}