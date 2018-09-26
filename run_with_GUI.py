#@OpService ops
#@String(choices={"fluoConfigurator", "batchSegmentation", "batchPreprocessing", "batchFluoAnalysis", "guiSeg"}, style="radioButtonVertical") method
#@File(label="Select a directory", style="directory") d
#@String(value="tif") suffix
#@String(value="45_1") acq
#@Boolean doGaussian
#@Boolean doAlign
#@Integer(label="Fluo-Analysis method", value=0) met
if method=="fluoConfigurator":
    ops.run(method, d.getPath(), "maars_config.xml")
elif method == "batchSegmentation":
	ops.run(method, d.getPath(), "maars_config.xml", suffix, acq, False)
elif method == "batchPreprocessing":
	ops.run(method, d.getPath(), "maars_config.xml", suffix, acq, doGaussian, doAlign)
elif method == "guiSeg":
	from maars.gui import SegPombeMainDialog
	seg_dia = SegPombeMainDialog()
	seg_dia.showDialog()
else: 
    ops.run(method, d.getPath(), "maars_config.xml", suffix, acq, met)

