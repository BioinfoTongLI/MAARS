#@OpService ops
#@String(choices={"fluoConfigurator", "batchSegmentation", "batchPreprocessing", "batchFluoAnalysis"}, style="radioButtonVertical") method
#@File(label="Select a directory", style="directory") d
#@String(value="tif") suffix
#@Integer(label="Fluo-Analysis method", value=0) met
if method=="fluoConfigurator":
    ops.run(method, d.getPath(), "maars_config.xml")
elif method == "batchSegmentation":
	ops.run(method, d.getPath(), "maars_config.xml", suffix, False)
elif method == "batchPreprocessing":
	ops.run(method, d.getPath(), "maars_config.xml", suffix)
else:
    ops.run(method, d.getPath(), "maars_config.xml", suffix, met)
