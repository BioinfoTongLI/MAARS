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
    analyser = MaarsAnalyzer(params, pos=args.pos, calibration=args.calibration, interval=args.interval)

    spots = Loader.load_all_spots(analyser.csts.FLUO_SPOTS).swaplevel(0, 2).sort_index()
    spots.index.names = ["TimePoint", "Channel", "Cell", "ID_dot"]
    spots.to_hdf(analyser.csts.MITO_DIR / "original_spots.h5", key="spots")

    # spots = pd.read_hdf(analyser.csts.MITO_DIR / "original_spots.h5")
    # spots.index.names = ["TimePoint", "Channel", "Cell", "dot_id"]

    xy = spots[["POSITION_X", "POSITION_Y"]].astype(np.float)
    gb_t_cell = xy.groupby(["TimePoint", "Cell"])
    geos = gb_t_cell.apply(ComputeFeatures.evaluate, roi_features=analyser.cell_rois_features)
    geos.to_hdf(analyser.csts.MITO_DIR / "geos.h5", key="geos")

    geos = pd.read_hdf(analyser.csts.MITO_DIR / "geos.h5")
    # print(geos.columns)

    sp_lens = geos["sp_len"].unstack(level=[1])

    cells_has_sp_mask = sp_lens.apply(lambda one_series: one_series.any())
    cells_has_sp_mask.to_hdf(analyser.csts.MITO_DIR / "has_sp_mask.h5", key="has_sp_mask")

    sp_mito_filter = sp_lens.apply(is_in_mitosis,
                                   p=analyser.csts.DYNAMIC_P_THRESHOLD,
                                   min_seg_len=analyser.csts.MIN_SEG_LEN)
    sp_mito_filter.to_hdf(analyser.csts.MITO_DIR / "mitosis_mask.h5", key="mitosis")

    mito_elong = sp_lens.loc[:, sp_mito_filter]
    # plt.plot(mito_elong)
    # plt.show()

    kt_area = geos["kt_area"].unstack(level=[1])
    exist_kt_area = kt_area.apply(lambda one_series: one_series.any())
    exist_kt_area.to_hdf(analyser.csts.MITO_DIR / "has_kt_cloud.h5", key="has_kt_cloud")

    kt_area_data = sp_lens.loc[:, exist_kt_area]
    plt.plot(kt_area_data)
    plt.show()

    has_kt_area_and_sp_elong = exist_kt_area & sp_mito_filter
    has_kt_area_and_sp_elong.to_hdf(analyser.csts.MITO_DIR / "has_kt_cloud_and_mitosis.h5",
                                    key="has_kt_cloud_and_mitosis")

    has_kt_area_or_sp_elong = exist_kt_area | sp_mito_filter
    has_kt_area_or_sp_elong.to_hdf(analyser.csts.MITO_DIR / "has_kt_area_or_sp_elong.h5",
                                   key="has_kt_area_or_sp_elong")

    # Plotting.plot_elong(sp_lens, analyser.csts.MITO_DIR)

    # time_table = analyser.get_time_table(df_sp_lens)
    print("Done")
