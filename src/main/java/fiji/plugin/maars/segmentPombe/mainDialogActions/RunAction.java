package fiji.plugin.maars.segmentPombe.mainDialogActions;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.measure.Calibration;

import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JTextField;

import fiji.plugin.maars.segmentPombe.SegPombe;
import fiji.plugin.maars.segmentPombe.SegPombeMainDialog;
import fiji.plugin.maars.segmentPombe.SegPombeParameters;
import fiji.plugin.maars.utils.FileUtils;

public class RunAction implements ActionListener {

	private SegPombeMainDialog mainDialog;
	private SegPombeParameters parameters;
	private ImagePlus imgToAnalysis;

	// Parameters of the algorithm
	private SegPombe segPombe;
	private int sigma;
	private int maxWidth;
	private int maxHeight;
	private double minParticleSize;
	private double maxParticleSize;
	private float zFocus;
	private double solidityThreshold;
	private double meanGrayValThreshold;

	// Variable for the actionPerformed : if this is true that means units has
	// been checked
	private boolean unitsChecked;

	public RunAction(SegPombeMainDialog mainDialog, SegPombeParameters parameters) {
		this.mainDialog = mainDialog;
		this.parameters = parameters;
		// String savingPath = cB.getFileNameField().getText();
		// this.parameters.setSavingPath(savingPath);
		// this.parameters.setWillFiltrateUnusualShape(cB.getFilterUnususalCkb().getState());
		// this.parameters.setWillFiltrateWithMeanGrayValue(cB.getFilterWithMeanGreyValueCkb().getState());
		// this.parameters.setWillSaveBinaryImg(cB.getWillSaveBinaryImgCkb().getState());
		// this.parameters.setWillSaveCorrelationImg(cB.getWillSaveCorrelationImgCkb().getState());
		// this.parameters.setWillSaveDataFrame(cB.getWillSaveDataFrameCkb().getState());
		// this.parameters.setWillSaveFocusImage(cB.getWillSaveFocusImageCkb().getState());
		// this.parameters.setWillSaveRoi(cB.getwillSaveRoiCkb().getState());
		// this.parameters.setWillShowBinaryImg(cB.getWillShowBinaryImgCkb().getState());
		// this.parameters.setWillShowCorrelationImg(cB.getWillShowCorrelationImgCkb().getState());
		// this.parameters.setWillShowDataFrame(cB.getWillShowDataFrameCkb().getState());
		// this.parameters.setWillChangeScale(cB.getWillChangeScaleCkb().getState());
		// this.parameters.setWillShowFocusImage(cB.getWillShowFocusImageCkb().getState());
		unitsChecked = false;
		this.imgToAnalysis = parameters.getImageToAnalyze();

	}

	/**
	 * Method to get the state of all the result Option checkbox and return
	 * false if none of them are selected
	 */
	public boolean checkResultOptions() {

		return (mainDialog.getShowCorrelationImgCkb().getState() || mainDialog.getSaveBinaryImgCkb().getState()
				|| mainDialog.getShowDataFrameCkb().getState() || mainDialog.getSaveCorrelationImgCkb().getState()
				|| mainDialog.getSaveBinaryImgCkb().getState() || mainDialog.getSaveDataFrameCkb().getState()
				|| mainDialog.getShowFocusImageCkb().getState() || mainDialog.getSaveFocusImageCkb().getState()
				|| mainDialog.getSaveRoiCkb().getState());
	}

	/**
	 * Allows to convert width or height or depth of a size in microns to a size
	 * in pixels. To choose if you want to convert width or height you can use
	 * CellsBoundaries constants : WIDTH and HEIGHT and DEPTH. Return an int
	 * width or height in pixels.
	 */
	public int convertMicronToPixelSize(double micronSize, int widthOrHeightOrDepth) {
		int pixelSize = (int) Math.round(micronSize / parameters.getScale(widthOrHeightOrDepth));
		return pixelSize;
	}

	/**
	 * Allows to convert width or height or depth of a size in pixels to a size
	 * in microns. To choose if you want to convert width or height you can use
	 * CellsBoundaries constants : WIDTH and HEIGHT and DEPTH. Return an double
	 * width or height in microns.
	 */
	public double convertPixelToMicronSize(int pixelSize, int widthOrHeightOrDepth) {
		double micronSize = (double) parameters.getScale(widthOrHeightOrDepth) * pixelSize;

		return micronSize;
	}

	/**
	 * Method to check if the image is scaled and if the unit matches 'micron'
	 */
	public void checkUnitsAndScale() {
		System.out.println("Check if image is scaled");
		if (imgToAnalysis.getCalibration().scaled()) {

			if (imgToAnalysis.getCalibration().getUnit().equals("cm")) {
				imgToAnalysis.getCalibration().setUnit("micron");
				imgToAnalysis.getCalibration().pixelWidth = imgToAnalysis.getCalibration().pixelWidth * 10000;
				imgToAnalysis.getCalibration().pixelHeight = imgToAnalysis.getCalibration().pixelHeight * 10000;
			}

			System.out.println("Check if unit of calibration is micron");
			if (imgToAnalysis.getCalibration().getUnit().equals("micron")
					|| imgToAnalysis.getCalibration().getUnit().equals("µm")) {
				double[] scale = new double[3];
				scale[SegPombeParameters.WIDTH] = imgToAnalysis.getCalibration().pixelWidth;
				scale[SegPombeParameters.HEIGHT] = imgToAnalysis.getCalibration().pixelHeight;
				scale[SegPombeParameters.DEPTH] = imgToAnalysis.getCalibration().pixelDepth;

				parameters.setScale(scale);
				unitsChecked = true;
				System.out.println("Get and set calibration as scale");
			} else {
				IJ.error("Wrong scale unit",
						"The scale of your image must be in microns.\nTo change it you can go to Properties ...");
			}
		} else {
			IJ.error("No scale",
					"You must set a scale to your image\nif you want to enter measures in micron.\nYou can go to Properties ...");
		}
	}

	/**
	 * Method to change the image scale if it is bigger than supposed to the
	 * maximum width an height are in pixels
	 */
	public void changeScale(int maxWidth, int maxHeight) {
		int newWidth;
		int newHeight;
		int newMaxParticleSize;
		int newMinParticleSize;
		Calibration newCal = new Calibration();
		newCal.setUnit("micron");
		ImagePlus img = imgToAnalysis;
		if (img.getWidth() > maxWidth) {
			System.out.println("Image width is greater than maximum width allowed");

			newWidth = maxWidth;
			newHeight = (int) img.getHeight() * maxWidth / img.getWidth();

			newMinParticleSize = (int) minParticleSize * maxWidth / img.getWidth();
			newMaxParticleSize = (int) maxParticleSize * maxWidth / img.getWidth();

			if (img.getCalibration().scaled()) {

				newCal.pixelWidth = parameters.getScale(SegPombeParameters.WIDTH) * img.getWidth() / maxWidth;
				newCal.pixelHeight = parameters.getScale(SegPombeParameters.HEIGHT) * img.getWidth() / maxWidth;
				newCal.pixelDepth = parameters.getScale(SegPombeParameters.DEPTH);
			}

			System.out.println("New values : w = " + newWidth + " h = " + newHeight);
			if (newHeight > maxHeight) {
				System.out.println("New height is still greater than maximum height allowed");
				newHeight = maxHeight;
				newWidth = (int) img.getWidth() * maxHeight / img.getHeight();

				newMinParticleSize = (int) minParticleSize * maxHeight / img.getHeight();
				newMaxParticleSize = (int) maxParticleSize * maxHeight / img.getHeight();

				if (img.getCalibration().scaled()) {
					newCal.pixelWidth = parameters.getScale(SegPombeParameters.WIDTH) * img.getHeight() / maxHeight;
					newCal.pixelHeight = parameters.getScale(SegPombeParameters.HEIGHT) * img.getHeight() / maxHeight;
					newCal.pixelDepth = parameters.getScale(SegPombeParameters.DEPTH);
				}

				System.out.println("New values : w = " + newWidth + " h = " + newHeight);

			}

			rescale(newWidth, newHeight, newMinParticleSize, newMaxParticleSize, newCal);
		} else {
			if (img.getHeight() > maxHeight) {
				System.out.println("Image height is greater than maximum width allowed");

				newHeight = maxHeight;
				newWidth = (int) img.getWidth() * maxHeight / img.getHeight();

				newMinParticleSize = (int) minParticleSize * maxHeight / img.getHeight();
				newMaxParticleSize = (int) maxParticleSize * maxHeight / img.getHeight();

				if (img.getCalibration().scaled()) {
					newCal.pixelWidth = parameters.getScale(SegPombeParameters.WIDTH) * img.getHeight() / maxHeight;
					newCal.pixelHeight = parameters.getScale(SegPombeParameters.HEIGHT) * img.getHeight() / maxHeight;
					newCal.pixelDepth = parameters.getScale(SegPombeParameters.DEPTH);
				}

				System.out.println("New values : w = " + newWidth + " h = " + newHeight);

				if (newWidth > maxWidth) {
					System.out.println("New Width is still greater than maximum height allowed");

					newWidth = maxWidth;
					newHeight = (int) img.getHeight() * maxWidth / img.getWidth();

					if (img.getCalibration().scaled()) {
						newCal.pixelWidth = parameters.getScale(SegPombeParameters.WIDTH) * img.getWidth() / maxWidth;
						newCal.pixelHeight = parameters.getScale(SegPombeParameters.HEIGHT) * img.getWidth() / maxWidth;
						newCal.pixelDepth = parameters.getScale(SegPombeParameters.DEPTH);
					}

					System.out.println("New values : w = " + newWidth + " h = " + newHeight);

					newMinParticleSize = (int) minParticleSize * maxWidth / img.getWidth();
					newMaxParticleSize = (int) maxParticleSize * maxWidth / img.getWidth();

				}

				rescale(newWidth, newHeight, newMinParticleSize, newMaxParticleSize, newCal);
			}
		}
	}

	/**
	 * Method to change size and scale of the image to analyze : need to compute
	 * parameters before so it can be coherent
	 */
	public void rescale(int newWidth, int newHeight, int newMinParticleSize, int newMaxParticleSize,
			Calibration newCal) {

		minParticleSize = newMinParticleSize;
		maxParticleSize = newMaxParticleSize;

		System.out.println("min area = " + newMinParticleSize + " max area = " + newMaxParticleSize);

		ImageStack newImgStack = new ImageStack(newWidth, newHeight);

		for (int slice = 0; slice < imgToAnalysis.getNSlices(); slice++) {
			imgToAnalysis.setZ(slice);
			newImgStack.addSlice(imgToAnalysis.getProcessor().resize(newWidth, newHeight));
		}

		ImagePlus newImagePlus = new ImagePlus("rescaled_" + imgToAnalysis.getTitle(), newImgStack);
		newImagePlus.setCalibration(newCal);

		parameters.setImageToAnalyze(newImagePlus);

		checkUnitsAndScale();
	}

	/**
	 * Action performed when run Button is triggered. It checks all parameters
	 * then run Algorithm.
	 */
	public void actionPerformed(ActionEvent e) {
		
		Boolean parametersChecked = false;
		
		while (!parametersChecked){
			if (!FileUtils.isFilenameValid(mainDialog.getFileNameField().getText())){
				IJ.error("No file", "You must select a file to process");
				mainDialog.getFileNameField().setBackground(Color.ORANGE);
			};
		}
		
			// Check if none of the result ckeckBox is selected : in this case,
			// the user would not get any result
			System.out.println("Check if none of the result ckeckBox is selected");
			boolean thereIsAResult = checkResultOptions();
			if (!thereIsAResult) {
				IJ.error("No result possible", "You have not selected a way to see your results");
			} else {
				// check if typical cell size field has been filled
				System.out.println("check if typical cell Z size field has been filled");
				if (mainDialog.getTypicalSizeTf().getText().isEmpty()) {
					IJ.error("Missing parameter", "You need to fill the typical cell Z size field");
				} else {

					// if the unit chosen is a micron it must be converted
					System.out.println("Check if one of the unite used is micron");
					if (SegPombeParameters.MICRONS == mainDialog.getTypicalSizeUnitCombo().getSelectedIndex()
							|| SegPombeParameters.MICRONS == mainDialog.getmaxWidthUnitCombo().getSelectedIndex()
							|| SegPombeParameters.MICRONS == mainDialog.getMaxHeightUnitCombo().getSelectedIndex()
							|| SegPombeParameters.MICRONS == mainDialog.getMinParticleSizeUnitCombo().getSelectedIndex()
							|| SegPombeParameters.MICRONS == mainDialog.getMaxParticleSizeUnitCombo()
									.getSelectedIndex()) {
						checkUnitsAndScale();
					}

					// Get cell typical size and check if the user did not make
					// any mistake while filling the field
					double typicalCellSize = 0;
					System.out.println("Try to get typical size");
					try {
						typicalCellSize = Double.parseDouble(mainDialog.getTypicalSizeTf().getText());
					} catch (NumberFormatException nfe) {
						IJ.error("Wrong parameter", "The typical cell size is supposed to be a number");
					}
					// Check if the cell size is a size
					if (typicalCellSize <= 0) {
						IJ.error("Wrong parameter", "The typical cell size must be a positive, not null value");
					} else {

						// Convert size into pixels
						if (SegPombeParameters.MICRONS == cB.getSizeComUnit().getSelectedIndex() && unitsChecked) {
							sigma = convertMicronToPixelSize(typicalCellSize, SegPombeParameters.DEPTH);
							System.out.println("typical size is in micron, convert it in pixels : " + sigma);
						} else {
							if (SegPombeParameters.PIXELS == cB.getSizeComUnit().getSelectedIndex()) {
								sigma = (int) typicalCellSize;
							}
						}

						// Read minimum cell area
						System.out.println("Try to read minimum particle size");
						double minParticleSizeTemp = 0;

						try {
							minParticleSizeTemp = Double.parseDouble(cB.getMinParticleSizeField().getText());
						} catch (NumberFormatException nfe) {
							IJ.error("Wrong parameter", "The minimum area is supposed to be a number");
						}

						if (minParticleSizeTemp <= 0) {
							IJ.error("Wrong parameter", "The minimum area must be a positive not null value");
						} else {
							// Covert minimum area if needed to
							if (SegPombeParameters.MICRONS == cB.getMinParticleSizeComboUnit().getSelectedIndex()
									&& unitsChecked) {
								minParticleSize = minParticleSizeTemp
										* convertMicronToPixelSize(1, SegPombeParameters.WIDTH)
										* convertMicronToPixelSize(1, SegPombeParameters.HEIGHT);
								System.out.println("Cell Area is in micron, convert it in pixels : " + minParticleSize);
							} else {
								if (SegPombeParameters.PIXELS == cB.getMinParticleSizeComboUnit().getSelectedIndex()) {
									minParticleSize = minParticleSizeTemp;
								}
							}

							// Read maximum cell area
							System.out.println("Try to read maximum particle size");
							double maxParticleSizeTemp = 0;

							try {
								maxParticleSizeTemp = Double.parseDouble(cB.getMaxParticleSizeField().getText());
							} catch (NumberFormatException nfe) {
								IJ.error("Wrong parameter", "The maximum area is supposed to be a number");
							}

							if (maxParticleSizeTemp <= 0) {
								IJ.error("Wrong parameter", "The maximum area must be a positive not null value");
							} else {
								// Covert maximum area if needed to
								if (SegPombeParameters.MICRONS == cB.getMaxParticleSizeComboUnit().getSelectedIndex()
										&& unitsChecked) {
									maxParticleSize = maxParticleSizeTemp
											* convertMicronToPixelSize(1, SegPombeParameters.WIDTH)
											* convertMicronToPixelSize(1, SegPombeParameters.HEIGHT);
									System.out.println(
											"Cell Area is in micron, convert it in pixels : " + maxParticleSize);
								} else {
									if (SegPombeParameters.PIXELS == cB.getMaxParticleSizeComboUnit()
											.getSelectedIndex()) {
										maxParticleSize = maxParticleSizeTemp;
									}
								}

								// If the user chose to change the scale
								System.out.println("Check if user wants to change scale");
								if (cB.getScaleCkb().getState()) {
									// Check if the user entered correct values
									// then get the values and convert them if
									// necessary
									double newMaxWidth = 0;
									System.out.println("Try to read width value");
									try {
										newMaxWidth = Double.parseDouble(cB.getMaxWTextField().getText());
									} catch (NumberFormatException nfe) {
										IJ.error("Wrong parameter", "The maximum width is supposed to be a number");
									}

									if (newMaxWidth <= 0) {
										IJ.error("Wrong parameter",
												"The maximum width must be a positive not null value");
									} else {

										double newMaxHeight = 0;
										System.out.println("Try to read height value");
										try {
											newMaxHeight = Double.parseDouble(cB.getMaxHTextField().getText());
										} catch (NumberFormatException nfe) {
											IJ.error("Wrong parameter",
													"The maximum height is supposed to be a number");
										}

										if (newMaxHeight <= 0) {
											IJ.error("Wrong parameter",
													"The maximum height must be a positive not null value");
										} else {

											if (SegPombeParameters.MICRONS == cB.getMaxWComboUnit().getSelectedIndex()
													&& unitsChecked) {

												maxWidth = convertMicronToPixelSize(newMaxWidth,
														SegPombeParameters.WIDTH);
												System.out.println(
														"Width value is in micron, convert it in pixel : " + maxWidth);
											} else {
												if (SegPombeParameters.PIXELS == cB.getMaxWComboUnit()
														.getSelectedIndex()) {
													maxWidth = (int) newMaxWidth;
												}
											}
											if (SegPombeParameters.MICRONS == cB.getMaxHComboUnit().getSelectedIndex()
													&& unitsChecked) {

												maxHeight = convertMicronToPixelSize(newMaxHeight,
														SegPombeParameters.HEIGHT);
												System.out.println("Height value is in micron, convert it in pixel : "
														+ maxHeight);
											} else {
												if (SegPombeParameters.PIXELS == cB.getMaxHComboUnit()
														.getSelectedIndex()) {
													maxHeight = (int) newMaxHeight;
												}
											}
											// Then we can change scale
											System.out.println("Change scale");
											changeScale(maxWidth, maxHeight);

										}
									}
								}

								System.out.println("Check if user wants to precise z focus");
								if (cB.getPreciseZFocusCheckbox().getState()) {

									float zf = (imgToAnalysis.getNSlices() / 2);

									System.out.println("Try to read z focus value");
									try {
										zf = Float.parseFloat(cB.getPreciseZFocusTextField().getText());
									} catch (NumberFormatException nfe) {
										IJ.error("Wrong parameter",
												"The slice corresponding to focus\nis supposed to be a number\n1 <= focus slice <= "
														+ imgToAnalysis.getNSlices()
														+ "\nProgram will continue with default focus slice");
										cB.getPreciseZFocusCheckbox().setState(false);
									}
									if (zf <= 0) {
										IJ.error("Wrong parameter",
												"If you want to precise focus slice,\nplease fill correctly the field\n1 <= focus slice <= "
														+ imgToAnalysis.getNSlices()
														+ "\nProgram will continue with default focus slice");
										cB.getPreciseZFocusCheckbox().setState(false);
									}

									zFocus = zf - 1;

								} else {
									zFocus = (imgToAnalysis.getNSlices() / 2) - 1;
								}

								System.out.println("Check if user want to filter shape using solidity");
								if (cB.getFilterUnususalCkb().getState()) {
									double sldtThrshld = 0;

									System.out.println("Try to read solidity threshold value");
									try {
										sldtThrshld = Double.parseDouble(cB.getSolidityField().getText());
									} catch (NumberFormatException nfe) {
										IJ.error("Wrong parameter",
												"The solidity threshold value\nis supposed to be a number\nProgram will continue without filter");
										cB.getFilterUnususalCkb().setState(false);
									}
									if (sldtThrshld <= 0 || sldtThrshld > 1) {
										IJ.error("Wrong parameter",
												"If you want to filter unusual shapes,\nplease fill correctly the field\n0 < soliditythreshold <= 1\nProgram will continue without filter");
										cB.getFilterUnususalCkb().setState(false);
									} else {
										solidityThreshold = sldtThrshld;
									}

								} else {
									solidityThreshold = 0;
								}

								System.out.println("Check if user want to filter background using mean gray value");
								if (cB.getFilterWithMeanGreyValueCkb().getState()) {
									double mnGryVlThrshld = 0;

									System.out.println("Try to read solidity threshold value");
									try {
										mnGryVlThrshld = Double
												.parseDouble(cB.getMeanGreyValueThresholdField().getText());
									} catch (NumberFormatException nfe) {
										IJ.error("Wrong parameter",
												"The mean gray value threshold value\nis supposed to be a number\nProgram will continue without filter");
										cB.getFilterWithMeanGreyValueCkb().setState(false);
									}
									if (mnGryVlThrshld == 0) {
										IJ.error("Wrong parameter",
												"If you want to filter unusual background,\nplease fill correctly the field\nProgram will continue without filter");
										cB.getFilterWithMeanGreyValueCkb().setState(false);
									} else {
										meanGrayValThreshold = mnGryVlThrshld;
									}

								} else {
									meanGrayValThreshold = 0;
								}

								// Then Algorithm can be run
								System.out.println("Aaaaand ACTION!");

								cellBoundId = new CellsBoundariesIdentification(parameters, sigma, minParticleSize,
										maxParticleSize, cB.getDirection(), zFocus, solidityThreshold,
										meanGrayValThreshold, false, true);
								cellBoundId.identifyCellsBoundaries();
							}
						}
					}
				}
			}
		}
	}
}
