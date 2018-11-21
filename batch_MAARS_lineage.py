#@OpService ops
ids = []
n_round = 3
n_pos= 4
suffix = "tif"
do_seg = False
do_fluo = True
if do_fluo:
	for i in range(2, n_round):
	    for j in range(2, n_pos):
	    	ids.append("Round" + str(i) +"_Pos" + str(j) + "_1")
else:
	for i in range(n_round):
	    for j in range(n_pos):
	        for l in range(1, n_round+1):
	            ids.append("Round" + str(i) +"_Pos" + str(j) + "_" + str(l))
params = {
	"/media/tongli/Lineage_2018/Tests/20181106" : ids,
}

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
	        ops.run("batchFluoAnalysis", p, "maars_config.xml", suffix, acq, "0")