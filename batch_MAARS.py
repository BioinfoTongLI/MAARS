#@OpService ops
default_set_1 = "45_1"
default_set_2 = "85_1"
params = {
	"/media/tongli/transfer/screening-2018/20180920_wt_bub1" : [default_set_1],
}
suffix = "tif"
do_seg = False
do_fluo = True

if do_seg:
	for p in params:
	    for acq in params[p]:
	    	print(p, acq)
	        ops.run("batchSegmentation", p, "maars_config.xml", suffix, acq, False)

if do_fluo:
	for p in params:
	    for acq in params[p]:
	    	print(p, acq)
	        ops.run("batchPreprocessing", p, "maars_config.xml", suffix, acq, True, True)
	        ops.run("batchFluoAnalysis", p, "maars_config.xml", suffix, acq, "1")