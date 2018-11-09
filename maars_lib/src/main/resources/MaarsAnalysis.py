# coding: utf-8
# !/usr/bin/env python3
import os
import sys

from maarsanalyzer import Plotting
from maarsanalyzer.MaarsAnalyzer import MaarsAnalyzer
sys.path.append(os.path.abspath(".."))
from maarsanalyzer.analysers import ComputeFeatures
from maarsanalyzer.io_helpers import Loader
from maarsanalyzer.io_helpers.Consts import Constants
import pandas as pd
from pathlib import Path
idx = pd.IndexSlice

if __name__ == '__main__':
    args = Loader.set_attributes_from_cmd_line()
    root = Path(args.baseDir)
    params = Loader.load_parameters(root)
    analyser = MaarsAnalyzer(params, pos=args.pos, calibration=args.calibration,
                             to_ch5=args.ch5, interval=args.interval)

    paths_to_features = analyser.get_paths_to_feature()
    existing_paths = Loader.validate_paths(paths_to_features)
    mitoFilter = analyser.get_mitosis_filter(existing_paths)
    path_to_mito_features = existing_paths[mitoFilter]
    path_to_mito_kt_features = [Path(str(path).replace("CFP", "GFP")) for path in path_to_mito_features]
    df_sp_lens = analyser.get_raw_elongations(path_to_mito_features)
    df_kt_lens = analyser.get_raw_elongations(path_to_mito_kt_features)
    df_sp_lens.to_hdf(root / Constants.MITO_DIR / "mitotic_elongations.h5", key="elongation")
    df_kt_lens.to_hdf(root / Constants.MITO_DIR / "mitotic_kt_elongations.h5", key="elongation")
    # plotting#################
    Plotting.plot_elong(df_sp_lens, root / Constants.MITO_DIR)
    Plotting.plot_elong(df_kt_lens, root / Constants.MITO_DIR)

    path_to_mito_pole_spots = [Path.joinpath(p.parent.parent, analyser.csts.SPOTS, p.stem + '.xml')
                               for p in path_to_mito_features]
    validated_pole_paths = Loader.validate_paths(path_to_mito_pole_spots)
    poles = Loader.load_spots(validated_pole_paths)

    path_to_mito_kt_spots = [Path.joinpath(p.parent, p.stem.split("_")[0] + "_" + params[Constants.KT_CH] + ".xml")
                             for p in path_to_mito_pole_spots]
    validated_kt_paths = Loader.validate_paths(path_to_mito_kt_spots)
    kts = Loader.load_spots(validated_kt_paths)

    original_dots_features_df = ComputeFeatures.merge_kt_spb(poles.unstack(level=[0, 1]), kts.unstack(level=[0, 1]))
    original_dots_features_df.to_hdf(root / Constants.MITO_DIR / "original_dots_features.h5", "original_dots_features")

    cell_rois_features = pd.read_csv(analyser.csts.ROI_FEATURES, index_col=0)

    mito_cell_numbers = [int(x) for x in original_dots_features_df.index.levels[-2]]
    mito_cell_roi_features = cell_rois_features.loc[mito_cell_numbers]
    mito_cell_roi_features.to_hdf(root / Constants.MITO_DIR / "mito_cell_roi_features.h5", key="mito_cell_roi_features")
    geos = ComputeFeatures.compute_geometries(original_dots_features_df, mito_cell_roi_features)
    # print(geos.groupby("Cell").plot())
    # import matplotlib.pyplot as plt
    # plt.show()
    # time_table = analyser.get_time_table(df_sp_lens)
    # print(time_table)
    # import matplotlib.pyplot as plt
    # for c in sp_geos.keys():
    #     sp_geos[c].plot()
    # plt.show()

    geos.to_hdf(root / Constants.MITO_DIR / "calculated_geos.h5", key="geos")
    analyser.shutdown()
    print("Done")
