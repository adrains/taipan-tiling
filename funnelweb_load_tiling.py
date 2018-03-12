"""File for loading in the results of tiling
"""
import pickle
import cPickle

pkl_file_name = "/Users/adamrains/results/tiling/180222_2311_21_fw_tiling.pkl"
#pkl_file_name = "/Users/adamrains/Code/taipan-tiling/results/170712_1107_42_fw_tiling.pkl"

pkl_file = open(pkl_file_name, "rb")
(tiling, remaining_targets, run_settings) = cPickle.load(pkl_file)