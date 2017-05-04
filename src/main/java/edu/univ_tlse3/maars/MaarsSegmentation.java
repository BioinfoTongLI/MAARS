package edu.univ_tlse3.maars;

import edu.univ_tlse3.segmentPombe.ParametersProcessor;
import edu.univ_tlse3.segmentPombe.SegPombe;
import edu.univ_tlse3.segmentPombe.SegPombeParameters;
import ij.IJ;
import ij.ImagePlus;
import ij.measure.ResultsTable;

/**
 * Class to segment a multiple z-stack bright field image then find and record
 * cell shape Rois
 *
 * @author Tong LI
 */
public class MaarsSegmentation implements Runnable {
    private MaarsParameters parameters;
    private ResultsTable rt;
    private ImagePlus img_;

    /**
     * * Constructor :
     *
     * @param parameters : MAARS parameters (see class MaarsParameters)
     * @param img           image to segment
     */
    public MaarsSegmentation(MaarsParameters parameters, ImagePlus img) {
        img_ = img;
        this.parameters = parameters;
    }

    ResultsTable getRoiMeasurements() {
        return this.rt;
    }

    @Override
    public void run() {
        IJ.log("Prepare parameters for segmentation...");
        SegPombeParameters segPombeParam = new SegPombeParameters();

        segPombeParam.setImageToAnalyze(img_);
        segPombeParam.setSavingPath(parameters.getSavingPath());

        segPombeParam.setFilterAbnormalShape(
                Boolean.parseBoolean(parameters.getSegmentationParameter(MaarsParameters.FILTER_SOLIDITY)));

        segPombeParam.setFiltrateWithMeanGrayValue(
                Boolean.parseBoolean(parameters.getSegmentationParameter(MaarsParameters.FILTER_MEAN_GREY_VALUE)));
        segPombeParam.getImageToAnalyze().getCalibration().pixelDepth = Double
                .parseDouble(parameters.getSegmentationParameter(MaarsParameters.STEP));
        // Calibrate parameters
        ParametersProcessor process = new ParametersProcessor(segPombeParam);

        process.checkImgUnitsAndScale();
        process.changeScale(
                Integer.parseInt(parameters.getSegmentationParameter(MaarsParameters.NEW_MAX_WIDTH_FOR_CHANGE_SCALE)),
                Integer.parseInt(parameters.getSegmentationParameter(MaarsParameters.NEW_MAX_HEIGTH_FOR_CHANGE_SCALE)));

        segPombeParam = process.getParameters();

        segPombeParam.setSigma(
                (int) Math.round(Double.parseDouble(parameters.getSegmentationParameter(MaarsParameters.CELL_SIZE))
                        / Double.parseDouble(parameters.getSegmentationParameter(MaarsParameters.STEP))));

        segPombeParam.setMinParticleSize((int) Math
                .round(Double.parseDouble(parameters.getSegmentationParameter(MaarsParameters.MINIMUM_CELL_AREA))
                        / segPombeParam.getImageToAnalyze().getCalibration().pixelWidth)
                / segPombeParam.getImageToAnalyze().getCalibration().pixelHeight);

        segPombeParam.setMaxParticleSize((int) Math
                .round(Double.parseDouble(parameters.getSegmentationParameter(MaarsParameters.MAXIMUM_CELL_AREA))
                        / segPombeParam.getImageToAnalyze().getCalibration().pixelWidth)
                / segPombeParam.getImageToAnalyze().getCalibration().pixelHeight);

        segPombeParam.setSolidityThreshold(
                Double.parseDouble(parameters.getSegmentationParameter(MaarsParameters.SOLIDITY)));

        segPombeParam.setMeanGreyValueThreshold(
                Double.parseDouble(parameters.getSegmentationParameter(MaarsParameters.MEAN_GREY_VALUE)));
        IJ.log("Done.");
        // Main segmentation process
        IJ.log("Begin segmentation...");
        SegPombe segPombe = new SegPombe(segPombeParam);
        segPombe.createCorrelationImage();
        segPombe.convertCorrelationToBinaryImage();
        segPombe.analyseAndFilterParticles();
        segPombe.showAndSaveResultsAndCleanUp();
        IJ.log("Segmentation done");
        this.rt = segPombe.getRoiMeasurements();
    }
}
