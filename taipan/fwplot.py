"""Plotting functions to visualise/record the results of a FunnelWeb tiling run
"""
import taipan.core as tp
import taipan.fwtiling as fwtl
import matplotlib.patches as mpatches
import numpy as np
from collections import OrderedDict
from collections import Counter
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
from itertools import cycle
import glob

def plot_tile_pos(fig, gs, tiling):
    """
    """
    # Plot an Aitoff projection of the tile centre positions on-sky
    ax_tile_pos = fig.add_subplot(gs, projection='aitoff')
    ax_tile_pos.grid(True)
    
    # Count the number of tiles per field
    coords = Counter(["%f_%f" % (tile.ra, tile.dec) for tile in tiling])
    coords = np.array([[float(key.split("_")[0]), float(key.split("_")[1]), 
                      coords[key]] for key in coords.keys()])
    
    ax_tile_pos_plt = ax_tile_pos.scatter(np.radians(coords[:,0] - 180.), 
                          np.radians(coords[:,1]), c=coords[:,2], marker='o',
                          lw=0, s=9, cmap="rainbow")
    ax_tile_pos.set_title('Tile centre positions', y=1.1)
    ax_tile_pos.set_axisbelow(True)
    
    # Colour bar
    ax_tile_pos_plt.set_clim([np.min(coords[:,2]), np.max(coords[:,2])])
    cbar = plt.colorbar(ax_tile_pos_plt, orientation='horizontal')
    cbar.set_label("# Tiles")
    
    ax_tile_pos.tick_params(axis='both', which='major', labelsize=5)
    ax_tile_pos.tick_params(axis='both', which='minor', labelsize=5)
    
    return ax_tile_pos
    
def plot_box_and_whisker(fig, gs, gs_range, targets_per_tile, 
                         standards_per_tile, guides_per_tile):
    """
    """
    ax_bw = fig.add_subplot(gs[1,0])
    ax_bw.boxplot([targets_per_tile, standards_per_tile, guides_per_tile], 
                vert=False)
    ax_bw.set_yticklabels(['T', 'S', 'G'])
    ax_bw.set_title('Box-and-whisker plots of number of assignments')
    
    return ax_bw
    
def plot_tile_vs_completeness(fig, gs, gs_range, targets_per_tile, tiling, 
                              run_settings):
    """
    """
    ax_compl = fig.add_subplot(gs[2,0])
    targets_per_tile_sorted = sorted(targets_per_tile, key=lambda x: -1.*x)
    xpts = np.asarray(range(len(targets_per_tile_sorted))) + 1
    ypts = [np.sum(targets_per_tile_sorted[:i+1]) for i in xpts]
    ax_compl.plot(xpts, ypts, 'k-', lw=.9)
    ax_compl.plot(len(tiling), np.sum(targets_per_tile), 'ro',
             label='No. of tiles: %d' % len(tiling))
    ax_compl.hlines(run_settings["num_targets"], ax_compl.get_xlim()[0], 
               ax_compl.get_xlim()[1], lw=.75, colors='k', linestyles='dashed', 
               label='100% completion')
    ax_compl.hlines(0.975 * run_settings["num_targets"], ax_compl.get_xlim()[0], 
               ax_compl.get_xlim()[1], lw=.75, colors='k', linestyles='dashdot', 
               label='97.5% completion')
    ax_compl.hlines(0.95 * run_settings["num_targets"], ax_compl.get_xlim()[0], 
               ax_compl.get_xlim()[1], lw=.75, colors='k', linestyles='dotted', 
               label='95% completion')
    ax_compl.legend(loc='lower right', 
               title='Time to %3.1f comp.: %dm:%2.1fs' % (
               run_settings["tiling_completeness"] * 100., 
               int(np.floor(run_settings["mins_to_complete"])), 
               run_settings["mins_to_complete"] % 60.))
    ax_compl.set_title('Completeness progression')
    ax_compl.set_xlabel('No. of tiles')
    ax_compl.set_ylabel('No. of assigned targets')
    ax_compl.set_xlim([0, 1.05*len(tiling)])
    ax_compl.set_ylim([0, 1.05*run_settings["num_targets"]])
    
    return ax_compl
    
def plot_settings_table(fig, gs, run_settings):
    """
    """
    ax_tab = fig.add_subplot(gs)
    col_labels = ("Parameter", "Value")
    ax_tab.axis("off")
    
    # Prep table
    fmt_table = []
    
    for key, value in zip(run_settings.keys(), run_settings.values()):
        # Show only file name (not full path) for space constraints
        if (type(value) is str or type(value) is np.string_) and "/" in value:
            fmt_table.append([key, value.split("/")[-1]])
        # Format any short lists as strings, don't show long lists
        elif type(value) is list:
            if len(value) <= 10:
                fmt_table.append([key, str(value)])
        # Leave all other values as is (but don't plot the description, as it
        # is already at the top of the plot and likely to be long
        elif key != "description":
            fmt_table.append([key, value])
    
    settings_tab = ax_tab.table(cellText=fmt_table, colLabels=col_labels, 
                             loc="center")
    #if plot_other:
        #settings_tab.set_fontsize(7)
        #settings_tab.scale(1.0, 0.6)
    #else:
    settings_tab.auto_set_font_size(False)
    settings_tab.set_fontsize(9.5)
        #settings_tab.scale(2.2, 1.2)
    
    ax_tab.set_title("Run Settings & Overview", y=0.96)
    
    return ax_tab    
    
def plot_targets_per_tile(fig, gs, targets_by_mag_range, targets_per_tile,
                          run_settings, colour, i, label, unique_targets_range,
                          tiles_by_mag_range):
    """
    """
    # Plot a histogram of the number of targets per tile
    ax_tar_pt = fig.add_subplot(gs)
    ax_tar_pt.hist(targets_by_mag_range[i], 
                 bins=np.arange(0, max(targets_per_tile)+1, 1), 
                 color=colour, align='right', label=label)
    ax_tar_pt.vlines(run_settings["TARGET_PER_TILE"], ax_tar_pt.get_ylim()[0],  
                   ax_tar_pt.get_ylim()[1], linestyles='dashed', colors='k', 
                   label='Ideally-filled tile')
    ax_tar_pt.legend(loc='upper center')
    ax_tar_pt.set_xlabel('No. of targets per tile')
    ax_tar_pt.set_ylabel('Frequency')
    ax_tar_pt.set_yscale('log')
    ax_tar_pt.set_xlim(0, max(targets_per_tile) + 1)
    ax_tar_pt.xaxis.set_major_locator(ticker.MultipleLocator(10))
    
    ax_tar_pt.tick_params(axis='both', which='major', labelsize=5)
    ax_tar_pt.tick_params(axis='both', which='minor', labelsize=5)
    
    if len(targets_by_mag_range[i]) > 0:
        tile_mean = np.mean(targets_by_mag_range[i])
        tile_median = np.median(targets_by_mag_range[i])
    else:
        tile_mean = 0
        tile_median = 0
    
    # Now plot text describing the number of tiles, targets, stats, and
    # completion target for all histogram plots (with the exception of 
    # skipping completion target on the total histogram)
    if i > 0:
        completeness = run_settings["completeness_targets"][i-1] * 100
        ax_tar_pt.text(0.5, 0.4, "%5.2f %% Completion" % completeness,
                     ha="center", transform=ax_tar_pt.transAxes)
    
    ax_tar_pt.text(0.5, 0.5, "Mean: %i, Median: %i" % (tile_mean, 
                                                     tile_median),
                 ha="center", transform=ax_tar_pt.transAxes)
    
    ax_tar_pt.text(0.5, 0.6, 
                 "{:,} Unique Targets".format(unique_targets_range[i]),
                 ha="center", transform=ax_tar_pt.transAxes)
                 
    ax_tar_pt.text(0.5, 0.7, 
                 "{:,} Tiles".format(len(tiles_by_mag_range[i])),
                 ha="center", transform=ax_tar_pt.transAxes)
                 
    return ax_tar_pt
    
def plot_standards_per_tile(fig, gs, standards_by_mag_range, 
                            standards_per_tile, run_settings, i, colour, 
                            label):
    """
    """
    # Plot a histogram of the number of standards per tile
    ax_std_pt = fig.add_subplot(gs)
    ax_std_pt.hist(standards_by_mag_range[i], 
                 bins=max(standards_per_tile), color=colour, 
                 align='right', label=label)
    ax_std_pt.vlines(run_settings["STANDARDS_PER_TILE"], 
                   ax_std_pt.get_ylim()[0], ax_std_pt.get_ylim()[1], 
                   linestyles='dashed', colors='k', 
                   label='Ideally-filled tile')
    ax_std_pt.vlines(run_settings["STANDARDS_PER_TILE_MIN"],
                   ax_std_pt.get_ylim()[0], ax_std_pt.get_ylim()[1],
                   linestyles='dotted',  colors='k', 
                   label='Minimum standards per tile')
    ax_std_pt.set_xlabel('No. of standards per tile')
    ax_std_pt.set_ylabel('Frequency')
    ax_std_pt.legend(loc='upper center')
    ax_std_pt.set_xlim(0, max(standards_per_tile) + 1)
    #ax_std_pt.xaxis.set_major_locator(ticker.MultipleLocator(2))
    
    ax_std_pt.tick_params(axis='both', which='major', labelsize=5)
    ax_std_pt.tick_params(axis='both', which='minor', labelsize=5)
    
    if len(standards_by_mag_range[i]) > 0:
        standard_mean = np.mean(standards_by_mag_range[i])
        standard_median = np.median(standards_by_mag_range[i])
    else:
        standard_mean = 0
        standard_median = 0
    
    ax_std_pt.text(ax_std_pt.get_xlim()[1]/2, ax_std_pt.get_ylim()[1]/2,
                 "Mean: %i, Median: %i" % (standard_mean, standard_median), 
                 ha="center")    
                 
    return ax_std_pt
    
def plot_guides_per_tile(fig, gs, guides_by_mag_range, guides_per_tile,
                         run_settings, i, colour, label):
    """
    """
    # Plot a histogram of the number of guides per tile
    ax_gde_pt = fig.add_subplot(gs)
    ax_gde_pt.hist(guides_by_mag_range[i], bins=max(guides_per_tile), 
                 color=colour, align='right', label=label)
    ax_gde_pt.vlines(run_settings["GUIDES_PER_TILE"], ax_gde_pt.get_ylim()[0], 
                   ax_gde_pt.get_ylim()[1], linestyles='dashed', colors='k', 
                   label='Ideally-filled tile')
    ax_gde_pt.vlines(run_settings["GUIDES_PER_TILE_MIN"], 
                   ax_gde_pt.get_ylim()[0], ax_gde_pt.get_ylim()[1],
                   linestyles='dotted', colors='k', 
                   label='Minimum guides per tile')
    ax_gde_pt.set_xlabel('No. of guides per tile')
    ax_gde_pt.set_ylabel('Frequency')
    ax_gde_pt.legend(loc='upper center')
    ax_gde_pt.set_xlim(0, max(guides_per_tile) + 1)
    ax_gde_pt.xaxis.set_major_locator(ticker.MultipleLocator(1))
    
    ax_gde_pt.tick_params(axis='both', which='major', labelsize=5)
    ax_gde_pt.tick_params(axis='both', which='minor', labelsize=5)
    
    if len(guides_by_mag_range[i]) > 0:
        guides_mean = np.mean(guides_by_mag_range[i])
        guides_median = np.median(guides_by_mag_range[i])
    else:
        guides_mean = 0
        guides_median = 0
    
    ax_gde_pt.text(ax_gde_pt.get_xlim()[1]/2, ax_gde_pt.get_ylim()[1]/2,
                 "Mean: %i, Median: %i" % (guides_mean, guides_median), 
                                           ha="center") 
    
    return ax_gde_pt 
    
    
def plot_sky_per_tile(fig, gs, sky_by_mag_range, sky_per_tile, run_settings, 
                      i, colour, label):
    # Plot a histogram of the number of sky per tile
    ax_sky_pt = fig.add_subplot(gs)
    ax_sky_pt.hist(sky_by_mag_range[i], bins=max(sky_per_tile), 
                 color=colour, align='right', label=label)
    ax_sky_pt.vlines(run_settings["SKY_PER_TILE"], ax_sky_pt.get_ylim()[0], 
                   ax_sky_pt.get_ylim()[1], linestyles='dashed', colors='k', 
                   label='Ideally-filled tile')
    ax_sky_pt.vlines(run_settings["SKY_PER_TILE_MIN"], 
                   ax_sky_pt.get_ylim()[0], ax_sky_pt.get_ylim()[1],
                   linestyles='dotted', colors='k', 
                   label='Minimum sky per tile')
    ax_sky_pt.set_xlabel('No. of sky per tile')
    ax_sky_pt.set_ylabel('Frequency')
    ax_sky_pt.legend(loc='upper center')
    ax_sky_pt.set_xlim(0, max(sky_per_tile) + 1)
    ax_sky_pt.xaxis.set_major_locator(ticker.MultipleLocator(1))
    
    if len(sky_by_mag_range[i]) > 0:
        sky_mean = np.mean(sky_by_mag_range[i])
        sky_median = np.median(sky_by_mag_range[i])
    else:
        sky_mean = 0
        sky_median = 0
    
    ax_sky_pt.text(ax_sky_pt.get_xlim()[1]/2, ax_sky_pt.get_ylim()[1]/2,
                 "Mean: %i, Median: %i" % (sky_mean, sky_median), 
                                           ha="center") 

    return ax_sky_pt 
    
def plot_target_histogram_grid(fig, gs, gs_xy, targets_by_mag_range, 
                               targets_per_tile, unique_targets_range, 
                               tiles_by_mag_range, standards_by_mag_range, 
                               standards_per_tile, guides_by_mag_range, 
                               guides_per_tile, sky_by_mag_range, sky_per_tile,
                               run_settings, tile_count_labels):
    """
    """
    # Initialise axes lists, and colour format
    ax_tar_pt = []
    ax_std_pt = []
    ax_gde_pt = []
    ax_sky_pt = []
    plt_colours = ["MediumSeaGreen","SkyBlue","Gold","Orange", "Tomato"]
    colour_cycler = cycle(plt_colours)
    
    for i, label in enumerate(tile_count_labels): 
        # Match the colours for each magnitude range
        colour = next(colour_cycler)
        
        # Number of Targets per Tile
        ax_tar_pt.append(plot_targets_per_tile(fig, gs[gs_xy[0]+i, gs_xy[1]+0], 
                                               targets_by_mag_range, 
                                               targets_per_tile, run_settings, 
                                               colour, i, label, 
                                               unique_targets_range,
                                               tiles_by_mag_range))
         
        # Number of standards per tile
        ax_std_pt.append(plot_standards_per_tile(fig, 
                                                 gs[gs_xy[0]+i, gs_xy[1]+1], 
                                                 standards_by_mag_range,
                                                 standards_per_tile, 
                                                 run_settings, i, colour, 
                                                 label)) 
        
        # Number of guides per tile
        ax_gde_pt.append(plot_guides_per_tile(fig, gs[gs_xy[0]+i, gs_xy[1]+2], 
                                              guides_by_mag_range, 
                                              guides_per_tile, run_settings, i,
                                              colour, label))
    
        # Number of sky per tile
        ax_sky_pt.append(plot_sky_per_tile(fig, gs[gs_xy[0]+i, gs_xy[1]+3], 
                         sky_by_mag_range, sky_per_tile, run_settings, i, 
                         colour, label))

    # Set plot titles
    ax_tar_pt[0].set_title("Stars per Tile")
    ax_std_pt[0].set_title("Standards per Tile")
    ax_gde_pt[0].set_title("Guides per Tile")
    ax_sky_pt[0].set_title("Sky per Tile")
    
    return [ax_tar_pt, ax_std_pt, ax_gde_pt, ax_sky_pt]


def plot_priority_histogram(fig, gs, priorities_by_mag_range, run_settings):
    """
    """
    priorities = [0, 1, 2, 3, 4, 5]
    
    for mag_range in xrange(len(run_settings["mag_ranges"])):
        ax_priority = fig.add_subplot(gs)
        ax_priority.hist(sky_by_mag_range[i], bins=max(sky_per_tile), 
                     color=colour, align='right', label=label)
    
        ax_sky_pt.set_xlabel("Priority Level")
        ax_sky_pt.set_ylabel("Number of Targets Observed")
        ax_sky_pt.legend(loc='upper right')
        ax_sky_pt.set_xlim(0, 6)
        ax_sky_pt.xaxis.set_major_locator(ticker.MultipleLocator(1))
    

    
    ax_sky_pt.text(ax_sky_pt.get_xlim()[1]/2, ax_sky_pt.get_ylim()[1]/2,
                 "Mean: %i, Median: %i" % (sky_mean, sky_median), 
                                           ha="center") 

    return ax_sky_pt 


def plot_tiling(tiling, run_settings):
    """Function to plot an overview of a FunnelWeb tiling run.
    
    Plots:
    - Tile centres on sky using aitoff projection
    - Box & whisker plots of target/standard/guide assignments
    - Completeness progression (targets as a function of tiles)
    - Table of run_settings
    - Histograms of targets per tile for all tiles, plus each magnitude bin
    - Histograms of standards per tile for all tiles, plus each magnitude bin
    - Histograms of guides per tile for all tiles, plus each magnitude bin
    
    Parameters
    ----------
    tiling: list
        The list of TaipanTiles from a tiling run.
    
    run_settings: OrderedDict
        OrderedDict containing input settings and results from the tiling run.
        
    plot_other: boolean
        Whether to plot box-and-whisker and completeness target plots, or 
        instead plot the variable table in a larger format.
    """
    # Initialise plot, use GridSpec to have a NxM grid, write title
    plt.clf()
    gs = gridspec.GridSpec(7,7)
    gs.update(wspace=0.2)
    fig = plt.gcf()
    #fig.set_size_inches(8.3, 11.7)
    fig.set_size_inches(28.8, 24.) # A4 Paper
    plt.suptitle(run_settings["run_id"] + ": " + run_settings["description"], 
                 fontsize=10)
    
    # Separate targets, standards, and guides by magnitude range
    tiles_by_mag_range = [tiling]
    targets_by_mag_range = []
    standards_by_mag_range = []
    guides_by_mag_range = []
    sky_by_mag_range = []
    priorities_by_mag_range = []
    
    all_priorities = Counter({0:0, 1:0, 2:0, 3:0, 4:0, 5:0})
    
    # Initialise legend label list
    tile_count_labels = ["All"]
        
    for mrange in run_settings["mag_ranges"]:
        tiles_for_range = []
        targets_for_range = []
        standards_for_range = []
        guides_for_range = []
        sky_for_range = []
        priorities_for_range = {0:0, 1:0, 2:0, 3:0, 4:0, 5:0}
        
        for tile in tiling:
            if tile.mag_min == mrange[0] and tile.mag_max == mrange[1]:
                tiles_for_range.append(tile)
                targets_for_range.append(tile.count_assigned_targets_science())
                standards_for_range.append(
                    tile.count_assigned_targets_standard())
                guides_for_range.append(tile.count_assigned_targets_guide())
                sky_for_range.append(tile.count_assigned_targets_sky())
                
                tile_priorities = tile.count_target_priorities()
                for level in tile_priorities.keys():
                    priorities_for_range[level] += tile_priorities[level]
        
        tiles_by_mag_range.append(tiles_for_range)
        targets_by_mag_range.append(targets_for_range)
        standards_by_mag_range.append(standards_for_range)
        guides_by_mag_range.append(guides_for_range)
        sky_by_mag_range.append(sky_for_range)
        priorities_by_mag_range.append(priorities_for_range)
        all_priorities += Counter(priorities_for_range)
        
        tile_count_labels.append("Mag range: %s-%s" % (mrange[0], mrange[1]))
    
    targets_per_tile = run_settings["targets_per_tile"]
    standards_per_tile = run_settings["standards_per_tile"]
    guides_per_tile = run_settings["guides_per_tile"]
    sky_per_tile = run_settings["sky_per_tile"]
    
    # Insert the non-magnitude range specific numbers
    targets_by_mag_range.insert(0, [t.count_assigned_targets_science() 
                                  for t in tiling]) 
    standards_by_mag_range.insert(0, [t.count_assigned_targets_standard() 
                                      for t in tiling]) 
    guides_by_mag_range.insert(0, [t.count_assigned_targets_guide() 
                                   for t in tiling])
    sky_by_mag_range.insert(0, [t.count_assigned_targets_sky() 
                                for t in tiling])      
    
    priorities_by_mag_range.insert(0, all_priorities)  
                                
    # Now count the number of unique targets per magnitude range
    unique_targets_range =  []
    for mag_range in tiles_by_mag_range:
        unique, tot, dup = fwtl.count_unique_science_targets(mag_range, True) 
        unique_targets_range.append(unique)                   
    
    # Tile positions
    ax_tile_pos = plot_tile_pos(fig, gs[0:2,0:2], tiling)

    # Table of Run Settings
    # Plot a table for referencing the run settings/results
    ax_tab = plot_settings_table(fig, gs[2:5,0:2], run_settings)
    
    # Create plots for all tiles, as well as each magnitude range
    ax_tgt_hists =  plot_target_histogram_grid(fig, gs, [0,2], 
                                               targets_by_mag_range, 
                                               targets_per_tile, 
                                               unique_targets_range, 
                                               tiles_by_mag_range, 
                                               standards_by_mag_range, 
                                               standards_per_tile, 
                                               guides_by_mag_range,
                                               guides_per_tile, 
                                               sky_by_mag_range, sky_per_tile,
                                               run_settings, tile_count_labels)
    # -------------------------------------------------------------------------
    # Priorities
    # -------------------------------------------------------------------------
    priorities = [0, 1, 2, 3, 4, 5]
    
    #import pdb
    #pdb.set_trace()
    for range_i, label in enumerate(tile_count_labels):
        priority_count = [priorities_by_mag_range[range_i][i] for i in priorities]
        ax_priority = fig.add_subplot(gs[6,2+range_i])
        ax_priority.bar(priorities, priority_count, 
                     color="DarkOrchid", label=label)
    
        ax_priority.hlines(50000, ax_priority.get_xlim()[0], 
               ax_priority.get_xlim()[1], lw=.75, colors='k', linestyles='dashed', 
               label="50,000")
        ax_priority.set_xlabel("Priority Level")
        ax_priority.set_ylabel("Number of Targets Observed")
        ax_priority.legend(loc='upper center')
        ax_priority.set_xlim(0, 6)
        ax_priority.xaxis.set_major_locator(ticker.MultipleLocator(1))
    
        ax_priority.tick_params(axis='both', which='major', labelsize=7)
        ax_priority.tick_params(axis='both', which='minor', labelsize=7)
        
        ax_priority.grid(True)
        
        ax_priority.set_yscale('log')
    
    #ax_sky_pt.text(ax_sky_pt.get_xlim()[1]/2, ax_sky_pt.get_ylim()[1]/2,
                 #"Mean: %i, Median: %i" % (sky_mean, sky_median), 
                                           #ha="center") 
    
    # -------------------------------------------------------------------------
    
    # Save plot
    name = "results/" + run_settings["run_id"] + "_tiling_run_overview.pdf"
    fig.savefig(name, fmt='pdf')
    
    
def create_tiling_visualisation(tiling, run_settings, increment=1000):
    """
    For increasingly larger slices of the tiling set, run the plotting code.
    
    To create a video from these:
    ffmpeg -framerate 4 -pattern_type glob -i "*.png" -c:v libx264 
        -pix_fmt yuv420p xx.mp4
    
    Parameters
    ----------
    tiling: list
        The list of TaipanTiles from a tiling run.
    run_settings: OrderedDict
        OrderedDict containing input settings and results from the tiling run.
    increment: int
        The number of new tiles to include in each plot/frame (i.e. 
        increment=1000 means that each iteration of the loop will create a plot
        with 1000 more tiles than the last).
    """
    # Create the list of tile increments, ensuring the final number is the 
    # complete plot
    steps = list(np.arange(increment, len(tiling), increment))
    
    if steps[-1] != len(tiling):
        steps.append(len(tiling))
    
    # Generate a plot for each increment
    for frame, step in enumerate(steps):
        plot_tiling(tiling[:step], run_settings)
        
        name = "results/visualisation/%s_overview_f%04i.png" % (
                                                        run_settings["run_id"], 
                                                        frame)
        plt.savefig(name)