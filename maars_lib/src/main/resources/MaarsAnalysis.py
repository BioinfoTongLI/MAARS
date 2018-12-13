# coding: utf-8
# !/usr/bin/env python3

from maarsanalyzer.classification.maarsClassic import is_in_mitosis
from maarsanalyzer.MaarsAnalyzer import MaarsAnalyzer
from maarsanalyzer.analysers import ComputeFeatures
from maarsanalyzer.io_helpers import Loader
import pandas as pd
from pathlib import Path
import numpy as np
idx = pd.IndexSlice

if __name__ == '__main__':

    args = Loader.set_attributes_from_cmd_line()
    root = Path(args.baseDir)
    params = Loader.load_parameters(root)
    analyser = MaarsAnalyzer(params, pos=args.pos, calibration=args.calibration, interval=args.interval)
    analyser.cell_rois_features.to_hdf(analyser.csts.MITO_DIR / "roi_features.h5", key="roi_features")

    spots = Loader.load_all_spots(analyser.csts.FLUO_SPOTS).swaplevel(0, 2).sort_index()
    spots.index.names = ["TimePoint", "Channel", "Cell", "ID_dot"]
    spots.to_hdf(analyser.csts.MITO_DIR / "original_spots.h5", key="spots")
    print("Spots saved.")
    # spots = pd.read_hdf(analyser.csts.MITO_DIR / "original_spots.h5")

    xy = spots[["POSITION_X", "POSITION_Y"]].astype(np.float)
    gb_names = ["TimePoint", "Cell"]
    gb_t_cell = xy.groupby(gb_names)
    geos = analyser.apply_parallel(gb_t_cell, ComputeFeatures.evaluate)
    geos.index.names = gb_names

    geos.to_hdf(analyser.csts.MITO_DIR / "geos.h5", key="geos")
    print("Geos saved.")
    # geos = pd.read_hdf(analyser.csts.MITO_DIR / "geos.h5")

    sp_lens = geos["sp_len"].unstack(level=[1]).astype(np.float)

    cells_has_sp_mask = sp_lens.apply(lambda one_series: one_series.any())
    cells_has_sp_mask.to_hdf(analyser.csts.MITO_DIR / "has_sp_mask.h5", key="has_sp_mask")

    sp_mito_filter = sp_lens.apply(is_in_mitosis,
                                   p=analyser.csts.DYNAMIC_P_THRESHOLD,
                                   min_seg_len=analyser.csts.MIN_SEG_LEN)
    sp_mito_filter.to_hdf(analyser.csts.MITO_DIR / "mitosis_mask.h5", key="mitosis")
    print("Mitosis mask saved")

    mito_elong = sp_lens.loc[:, sp_mito_filter]

    kt_area = geos["kt_area"].unstack(level=[1]).astype(np.float)
    exist_kt_area = kt_area.apply(lambda one_series: one_series.any())
    exist_kt_area.to_hdf(analyser.csts.MITO_DIR / "has_kt_cloud.h5", key="has_kt_cloud")
    print("Kt cloud mask saved.")
    kt_area_data = sp_lens.loc[:, exist_kt_area]

    has_kt_area_and_sp_elong = exist_kt_area & sp_mito_filter
    has_kt_area_and_sp_elong.to_hdf(analyser.csts.MITO_DIR / "has_kt_cloud_and_mitosis.h5",
                                    key="has_kt_cloud_and_mitosis")

    has_kt_area_or_sp_elong = exist_kt_area | sp_mito_filter
    has_kt_area_or_sp_elong.to_hdf(analyser.csts.MITO_DIR / "has_kt_area_or_sp_elong.h5",
                                   key="has_kt_area_or_sp_elong")

    # plt.plot(sp_lens.loc[:, has_kt_area_or_sp_elong])
    # plt.show()

    # Plotting.plot_elong(sp_lens, analyser.csts.MITO_DIR)

    # time_table = analyser.get_time_table(df_sp_lens)
    print("Done")