#@OpService ops
params = {
	"/media/tongli/0ABC6EF952B52BF5/screening_nikon/calibrated/20180918_mph1_bub3" : ["45", "85"], 
	"/media/tongli/0ABC6EF952B52BF5/screening_nikon/calibrated/20180918_wt_mad3" : ["85" , "140"]}
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