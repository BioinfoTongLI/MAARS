# coding: utf-8
# !/usr/bin/env python3
import os
import sys

from maarsanalyzer import Plotting
from maarsanalyzer.MaarsAnalyzer import MaarsAnalyzer
sys.path.append(os.path.abspath(".."))
from maarsanalyzer.analysers import ComputeFeatures
from maarsanalyzer.io_helpers import Writer, Loader
from maarsanalyzer.io_helpers.Consts import Constants
import pandas as pd
import numpy as np
from pathlib import Path

if __name__ == '__main__':
    args = Loader.set_attributes_from_cmd_line()
    root = Path(args.baseDir)
    params = Loader.load_parameters(root)
    analyser = MaarsAnalyzer(params, pos=args.pos, calibration=args.calibration, to_ch5=args.ch5)
    paths_to_features = analyser.get_paths_to_feature()
    existing_paths = Loader.validate_paths(paths_to_features)
    mitoFilter = analyser.get_mitosis_filter(existing_paths)
    path_to_mito_features = existing_paths[mitoFilter]

    path_to_mito_pole_spots = np.array([Path.joinpath(p.parent.parent, analyser.csts.SPOTS, p.stem + '.xml')
                                        for p in path_to_mito_features])
    validated_pole_paths = Loader.validate_paths(path_to_mito_pole_spots)
    dict_id_poles = Loader.load_spots(validated_pole_paths)

    path_to_mito_kt_spots = np.array(
        [Path.joinpath(p.parent, p.stem.split("_")[0] + "_" + params[Constants.KT_CH] + ".xml")
         for p in path_to_mito_pole_spots])
    validated_kt_paths = Loader.validate_paths(path_to_mito_kt_spots)
    dict_id_kts = Loader.load_spots(validated_kt_paths)

    df_sp_lens = analyser.get_elongations(path_to_mito_features)

    # plotting#################
    Plotting.plot_elong(df_sp_lens, root / Constants.MITO_DIR)

    valid_cells = ComputeFeatures.cells_contain_kts(dict_id_poles, dict_id_kts)
    print(valid_cells)
    dict_id_two_poles, dict_id_kts_of_two_poles = \
        ComputeFeatures.remove_single_poles(dict_id_poles, dict_id_kts, valid_cells)

    all_cell_rois = pd.read_csv(analyser.csts.ROI_FEATURES)
    valid_cell_rois = all_cell_rois.loc[valid_cells.astype(int)-1]
    valid_cell_rois.set_index(" ", inplace=True)
    valid_cell_rois.index = valid_cell_rois.index.astype(str)

    calculated_dots = ComputeFeatures.calculate_position_of_dots(dict_id_two_poles, dict_id_kts_of_two_poles,
                                                                 valid_cells)
    calculated_dots[""] = 0
    calculated_dots.set_index("", append=True, inplace=True)
    sp_geos = ComputeFeatures.calculate_sp_geos(dict_id_two_poles, valid_cell_rois,
                                                calculated_dots["sp_center_x"], calculated_dots["sp_center_y"])
    # time_table = analyser.get_time_table(sp_geos)
    # print(time_table)
    # import matplotlib.pyplot as plt
    # for c in sp_geos.keys():
    #     sp_geos[c].plot()
    # plt.show()

    fluo_dots = pd.merge(pd.concat(dict_id_poles).unstack(level=0), pd.concat(dict_id_kts).unstack(level=0),
                         left_index=True, right_index=True)
    for col in fluo_dots:
        if col[0] != "name":
            fluo_dots[col] = fluo_dots[col].astype(np.float16)
        else:
            fluo_dots[col] = fluo_dots[col].astype(str)

    unstacked_calculated_dots = calculated_dots.astype(np.str).unstack(level=0).astype(np.float16)
    all_dots = pd.concat([fluo_dots, unstacked_calculated_dots], axis=1, join="outer", sort=True)
    all_dots.index.set_names(['time', '#'], inplace=True)
    all_dots.to_hdf(root / Constants.MITO_DIR / "dots.h5", key="dots")

    kt_geos = ComputeFeatures.calculate_kt_geos(dict_id_kts_of_two_poles, valid_cells,
                                                calculated_dots["kt_center_x"], calculated_dots["kt_center_y"],
                                                calculated_dots.iloc[:, 4:10], calculated_dots.iloc[:, 10:16])
    tmp_idx = sp_geos.index
    sp_geos.index = sp_geos.index.set_levels(
        [tmp_idx.levels[0].astype(int), tmp_idx.levels[1]])
    for col in kt_geos:
        kt_geos[col] = kt_geos[col].astype(np.float64)
    pd.concat([sp_geos, kt_geos], axis=1, join="outer", sort=True).\
        to_hdf(root / Constants.MITO_DIR / "geos.h5", key="geos")
    df_sp_lens.to_hdf(root / Constants.MITO_DIR / "elongations.h5", key="elongation")
    print("Done")
