#@OpService ops
params = {
	"/media/tongli/transfer/screening-2018/20180921_mph1_mad3" : ["45_1", "85_1"]
	}
suffix = "tif"

for p in params:
    for acq in params[p]:
    	print(p, acq)
        ops.run("batchSegmentation", p, "maars_config.xml", suffix, acq, False)

for p in params:
    for acq in params[p]:
    	print(p, acq)
        ops.run("batchPreprocessing", p, "maars_config.xml", suffix, acq, False, True)
        ops.run("batchFluoAnalysis", p, "maars_config.xml", suffix, acq, "1")