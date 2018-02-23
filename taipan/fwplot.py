"""Plotting functions to visualise/record the results of a FunnelWeb tiling run
"""
import taipan.core as tp
import taipan.fwtiling as fwtl
import taipan.fwpriorities as fwpri
import matplotlib.patches as mpatches
import numpy as np
from collections import OrderedDict
from collections import Counter
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
from itertools import cycle
from matplotlib.backends.backend_pdf import PdfPages
import glob

def plot_tile_pos(fig, gs, tiling, axis_siz=5):
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
                          lw=0, s=30, cmap="rainbow")
    ax_tile_pos.set_title('Tile centre positions', fontsize=8, y=1.0)
    ax_tile_pos.set_axisbelow(True)
    
    
    # Colour bar
    ax_tile_pos_plt.set_clim([np.min(coords[:,2]), np.max(coords[:,2])])
    cbar = plt.colorbar(ax_tile_pos_plt, orientation='horizontal', pad=0.05)
    cbar.set_label("# Tiles", fontsize=axis_siz*2)
    
    #ax_tile_pos.tick_params(axis='both', which='major', labelsize=7)
    #ax_tile_pos.tick_params(axis='both', which='minor', labelsize=7)
    
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
    ax_compl.hlines(0.975 * run_settings["num_targets"], 
                    ax_compl.get_xlim()[0], ax_compl.get_xlim()[1], lw=.75, 
                    colors='k', linestyles='dashdot', label='97.5% completion')
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
    settings_tab.set_fontsize(4)
    #settings_tab.scale(2.2, 1.2)
    
    ax_tab.set_title("Run Settings & Overview", fontsize=8, y=1.075)
    
    return ax_tab    
    
    
def plot_targets_per_tile(fig, gs, target_count_by_range, targets_per_tile,
                          run_settings, colour, i, label, unique_targets_range,
                          tiles_by_mag_range, leg_siz=4, axis_siz=5, 
                          plot_xy_label=(False, False)):
    """
    """
    # Plot a histogram of the number of targets per tile
    ax_tar_pt = fig.add_subplot(gs)
    ax_tar_pt.hist(target_count_by_range[i], 
                 bins=np.arange(0, max(targets_per_tile)+1, 1), 
                 color=colour, align='right', label=label)
    ax_tar_pt.vlines(run_settings["TARGET_PER_TILE"], ax_tar_pt.get_ylim()[0],  
                   ax_tar_pt.get_ylim()[1], linestyles='dashed', colors='k', 
                   label='Ideally-filled tile', linewidth=1)
    ax_tar_pt.legend(loc='upper center', prop={'size': leg_siz})
    if plot_xy_label[0]: 
        ax_tar_pt.set_xlabel('No. of targets per tile', fontsize=axis_siz)
    if plot_xy_label[1]: 
        ax_tar_pt.set_ylabel('Frequency', fontsize=axis_siz)
    ax_tar_pt.set_yscale('log')
    ax_tar_pt.set_xlim(0, max(targets_per_tile) + 1)
    ax_tar_pt.xaxis.set_major_locator(ticker.MultipleLocator(20))
    
    ax_tar_pt.tick_params(axis='both', which='major', labelsize=5)
    ax_tar_pt.tick_params(axis='both', which='minor', labelsize=5)
    
    if len(target_count_by_range[i]) > 0:
        tile_mean = np.mean(target_count_by_range[i])
        tile_median = np.median(target_count_by_range[i])
    else:
        tile_mean = 0
        tile_median = 0
    
    # Now plot text describing the number of tiles, targets, stats, and
    # completion target for all histogram plots (with the exception of 
    # skipping completion target on the total histogram)
    if i > 0:
        completeness = run_settings["completeness_targets"][i-1] * 100
        ax_tar_pt.text(0.5, 0.4, "%5.2f %% Completion" % completeness,
                     ha="center", transform=ax_tar_pt.transAxes, 
                     fontsize=leg_siz)
    
    ax_tar_pt.text(0.5, 0.5, "Mean: %i, Median: %i" % (tile_mean, 
                                                     tile_median),
                 ha="center", transform=ax_tar_pt.transAxes, fontsize=leg_siz)
    
    ax_tar_pt.text(0.5, 0.6, 
                 "{:,} Unique Targets".format(unique_targets_range[i]),
                 ha="center", transform=ax_tar_pt.transAxes, fontsize=leg_siz)
                 
    ax_tar_pt.text(0.5, 0.7, 
                 "{:,} Tiles".format(len(tiles_by_mag_range[i])),
                 ha="center", transform=ax_tar_pt.transAxes, fontsize=leg_siz)
                 
    return ax_tar_pt
    
def plot_standards_per_tile(fig, gs, standard_count_by_range, 
                            standards_per_tile, run_settings, i, colour, 
                            label, leg_siz=4, axis_siz=5, 
                            plot_xy_label=(False, False)):
    """
    """
    # Plot a histogram of the number of standards per tile
    ax_std_pt = fig.add_subplot(gs)
    ax_std_pt.hist(standard_count_by_range[i], 
                 bins=max(standards_per_tile), color=colour, 
                 align='right', label=label)
    ax_std_pt.vlines(run_settings["STANDARDS_PER_TILE"], 
                   ax_std_pt.get_ylim()[0], ax_std_pt.get_ylim()[1], 
                   linestyles='dashed', colors='k', 
                   label='Ideally-filled tile', linewidth=1)
    ax_std_pt.vlines(run_settings["STANDARDS_PER_TILE_MIN"],
                   ax_std_pt.get_ylim()[0], ax_std_pt.get_ylim()[1],
                   linestyles='dotted',  colors='k', 
                   label='Minimum standards per tile', linewidth=1)
    if plot_xy_label[0]: 
        ax_std_pt.set_xlabel('No. of standards per tile', fontsize=axis_siz)
    if plot_xy_label[1]: 
        ax_std_pt.set_ylabel('Frequency', fontsize=axis_siz)
    ax_std_pt.legend(loc='upper center', prop={'size': leg_siz})
    ax_std_pt.set_xlim(0, max(standards_per_tile) + 1)
    #ax_std_pt.xaxis.set_major_locator(ticker.MultipleLocator(2))
    
    ax_std_pt.tick_params(axis='both', which='major', labelsize=5)
    ax_std_pt.tick_params(axis='both', which='minor', labelsize=5)
    
    if len(standard_count_by_range[i]) > 0:
        standard_mean = np.mean(standard_count_by_range[i])
        standard_median = np.median(standard_count_by_range[i])
    else:
        standard_mean = 0
        standard_median = 0
    
    ax_std_pt.text(ax_std_pt.get_xlim()[1]/2, ax_std_pt.get_ylim()[1]/2,
                 "Mean: %i, Median: %i" % (standard_mean, standard_median), 
                 ha="center", fontsize=leg_siz)    
                 
    return ax_std_pt
    
def plot_guides_per_tile(fig, gs, guide_count_by_range, guides_per_tile,
                         run_settings, i, colour, label, leg_siz=4, axis_siz=5, 
                         plot_xy_label=(False, False)):
    """
    """
    # Plot a histogram of the number of guides per tile
    ax_gde_pt = fig.add_subplot(gs)
    ax_gde_pt.hist(guide_count_by_range[i], bins=max(guides_per_tile), 
                 color=colour, align='right', label=label)
    ax_gde_pt.vlines(run_settings["GUIDES_PER_TILE"], ax_gde_pt.get_ylim()[0], 
                   ax_gde_pt.get_ylim()[1], linestyles='dashed', colors='k', 
                   label='Ideally-filled tile', linewidth=1)
    ax_gde_pt.vlines(run_settings["GUIDES_PER_TILE_MIN"], 
                   ax_gde_pt.get_ylim()[0], ax_gde_pt.get_ylim()[1],
                   linestyles='dotted', colors='k', 
                   label='Minimum guides per tile', linewidth=1)
    if plot_xy_label[0]:                
        ax_gde_pt.set_xlabel('No. of guides per tile', fontsize=axis_siz)
    if plot_xy_label[1]: 
        ax_gde_pt.set_ylabel('Frequency', fontsize=axis_siz)
    ax_gde_pt.legend(loc='upper center', prop={'size': leg_siz})
    ax_gde_pt.set_xlim(0, max(guides_per_tile) + 1)
    ax_gde_pt.xaxis.set_major_locator(ticker.MultipleLocator(1))
    
    ax_gde_pt.tick_params(axis='both', which='major', labelsize=5)
    ax_gde_pt.tick_params(axis='both', which='minor', labelsize=5)
    
    if len(guide_count_by_range[i]) > 0:
        guides_mean = np.mean(guide_count_by_range[i])
        guides_median = np.median(guide_count_by_range[i])
    else:
        guides_mean = 0
        guides_median = 0
    
    ax_gde_pt.text(ax_gde_pt.get_xlim()[1]/2, ax_gde_pt.get_ylim()[1]/2,
                 "Mean: %i, Median: %i" % (guides_mean, guides_median), 
                                           ha="center", fontsize=leg_siz) 
    
    return ax_gde_pt 
    
    
def plot_sky_per_tile(fig, gs, sky_count_by_range, sky_per_tile, run_settings, 
                      i, colour, label, leg_siz=4, axis_siz=5, 
                      plot_xy_label=(False, False)):
    # Plot a histogram of the number of sky per tile
    ax_sky_pt = fig.add_subplot(gs)
    ax_sky_pt.hist(sky_count_by_range[i], bins=max(sky_per_tile), 
                 color=colour, align='right', label=label)
    ax_sky_pt.vlines(run_settings["SKY_PER_TILE"], ax_sky_pt.get_ylim()[0], 
                   ax_sky_pt.get_ylim()[1], linestyles='dashed', colors='k', 
                   label='Ideally-filled tile', linewidth=1)
    ax_sky_pt.vlines(run_settings["SKY_PER_TILE_MIN"], 
                   ax_sky_pt.get_ylim()[0], ax_sky_pt.get_ylim()[1],
                   linestyles='dotted', colors='k', 
                   label='Minimum sky per tile', linewidth=1)
    if plot_xy_label[0]:               
        ax_sky_pt.set_xlabel('No. of sky per tile', fontsize=axis_siz)
    if plot_xy_label[1]:
        ax_sky_pt.set_ylabel('Frequency', fontsize=axis_siz)
    ax_sky_pt.legend(loc='upper center', prop={'size': leg_siz})
    ax_sky_pt.set_xlim(0, max(sky_per_tile) + 1)
    ax_sky_pt.xaxis.set_major_locator(ticker.MultipleLocator(1))
    
    ax_sky_pt.tick_params(axis='both', which='major', labelsize=5)
    ax_sky_pt.tick_params(axis='both', which='minor', labelsize=5)
    
    if len(sky_count_by_range[i]) > 0:
        sky_mean = np.mean(sky_count_by_range[i])
        sky_median = np.median(sky_count_by_range[i])
    else:
        sky_mean = 0
        sky_median = 0
    
    ax_sky_pt.text(ax_sky_pt.get_xlim()[1]/2, ax_sky_pt.get_ylim()[1]/2,
                 "Mean: %i, Median: %i" % (sky_mean, sky_median), 
                                           ha="center", fontsize=leg_siz) 

    return ax_sky_pt 
    
def plot_target_histogram_grid(fig, gs, gs_xy, target_count_by_range, 
                               targets_per_tile, unique_targets_range, 
                               tiles_by_mag_range, standard_count_by_range, 
                               standards_per_tile, guide_count_by_range, 
                               guides_per_tile, sky_count_by_range, 
                               sky_per_tile, run_settings, tile_count_labels):
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
        
        # Determine whether to plot x axis labels
        plot_x_label = (i + 1 == len(tile_count_labels))
        
        # Number of Targets per Tile
        ax_tar_pt.append(plot_targets_per_tile(fig, gs[gs_xy[0]+i, gs_xy[1]+0], 
                                               target_count_by_range, 
                                               targets_per_tile, run_settings, 
                                               colour, i, label, 
                                               unique_targets_range,
                                               tiles_by_mag_range,
                                               plot_xy_label=
                                               (plot_x_label, True)))
         
        # Number of standards per tile
        ax_std_pt.append(plot_standards_per_tile(fig, 
                                                 gs[gs_xy[0]+i, gs_xy[1]+1], 
                                                 standard_count_by_range,
                                                 standards_per_tile, 
                                                 run_settings, i, colour, 
                                                 label, plot_xy_label=
                                                 (plot_x_label, False)))
        
        # Number of guides per tile
        ax_gde_pt.append(plot_guides_per_tile(fig, gs[gs_xy[0]+i, gs_xy[1]+2], 
                                              guide_count_by_range, 
                                              guides_per_tile, run_settings, i,
                                              colour, label, plot_xy_label=
                                              (plot_x_label, False)))
    
        # Number of sky per tile
        ax_sky_pt.append(plot_sky_per_tile(fig, gs[gs_xy[0]+i, gs_xy[1]+3], 
                         sky_count_by_range, sky_per_tile, run_settings, i, 
                         colour, label, plot_xy_label=(plot_x_label, False)))

    # Set plot titles
    ax_tar_pt[0].set_title("Stars per Tile", fontsize=8)
    ax_std_pt[0].set_title("Standards per Tile", fontsize=8)
    ax_gde_pt[0].set_title("Guides per Tile", fontsize=8)
    ax_sky_pt[0].set_title("Sky per Tile", fontsize=8)
    
    return [ax_tar_pt, ax_std_pt, ax_gde_pt, ax_sky_pt]


def plot_priority_histogram(fig, gs, gs_xy, priority_count_by_range, 
                            run_settings, tile_count_labels, leg_siz=5, 
                            axis_siz=5):
    """
    """
    priorities = [0, 1, 2, 3, 4, 5]
    
    ax_priority = []
    
    for range_i, label in enumerate(tile_count_labels):
        count = [priority_count_by_range[range_i][i] for i in priorities]
        ax_priority.append(fig.add_subplot(gs[gs_xy[0]+range_i, gs_xy[1]]))
        ax_priority[range_i].bar(priorities, count, color="DarkOrchid", 
                        label=label)
        if range_i == 0:
            ax_priority[range_i].hlines(50000, 
                ax_priority[range_i].get_xlim()[0], 
                ax_priority[range_i].get_xlim()[1], lw=.75, colors='k', 
                linestyles='dashed', label="50,000", linewidth=1)
        #ax_priority[range_i].set_xlabel("Priority Level", fontsize=axis_siz)
        ax_priority[range_i].set_ylabel("Number of Targets Observed", 
                                        fontsize=axis_siz)
        ax_priority[range_i].legend(loc='upper center', prop={'size': leg_siz})
        ax_priority[range_i].set_xlim(0, 6)
        ax_priority[range_i].xaxis.set_major_locator(ticker.MultipleLocator(1))
        ax_priority[range_i].tick_params(axis='both', which='major', 
                                         labelsize=7)
        ax_priority[range_i].tick_params(axis='both', which='minor', 
                                         labelsize=7)
        #ax_priority[range_i].grid(True)
        ax_priority[range_i].set_yscale('log')
        
        ax_priority[range_i].tick_params(axis='both', which='major', 
                                         labelsize=5)
        ax_priority[range_i].tick_params(axis='both', which='minor', 
                                         labelsize=5)
        
        ax_priority[range_i].text(0.6, 0.5, "# Priority 3: %i" % count[3],
                 ha="center", transform=ax_priority[range_i].transAxes, 
                 fontsize=leg_siz)
    
        ax_priority[range_i].text(0.6, 0.6, "# Priority 4: %i" % count[4],
                 ha="center", transform=ax_priority[range_i].transAxes, 
                 fontsize=leg_siz)
                 
        ax_priority[range_i].text(0.6, 0.7, "# Priority 5: %i" % count[5],
                 ha="center", transform=ax_priority[range_i].transAxes, 
                 fontsize=leg_siz)
        
    ax_priority[0].set_title("Priority Completion", fontsize=8)
    
    ax_priority[-1].set_xlabel("Priority Level", fontsize=axis_siz)
        
    return ax_priority 


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
    # -------------------------------------------------------------------------
    # Preparation and Analysis
    # -------------------------------------------------------------------------
    # Separate targets, standards, and guides by magnitude range
    tiles_by_mag_range = [tiling]
    targets_by_mag_range = []
    
    target_count_by_range = []
    standard_count_by_range = []
    guide_count_by_range = []
    sky_count_by_range = []
    priority_count_by_range = []
    
    all_priorities = Counter({0:0, 1:0, 2:0, 3:0, 4:0, 5:0})
    
    # Initialise legend label list
    tile_count_labels = ["All"]
        
    for mrange in run_settings["mag_ranges"]:
        tiles_for_range = []
        targets_for_range = []
        
        target_count = []
        standard_count = []
        guide_count = []
        sky_count = []
        priority_count = {0:0, 1:0, 2:0, 3:0, 4:0, 5:0}
        
        for tile in tiling:
            if tile.mag_min == mrange[0] and tile.mag_max == mrange[1]:
                # Collect the tile and target objects for this mag range
                tiles_for_range.append(tile)
                targets_for_range.extend(tile.get_assigned_targets_science())
                
                # Count assigned targets, standards, guides, and sky
                target_count.append(tile.count_assigned_targets_science())
                standard_count.append(
                    tile.count_assigned_targets_standard())
                guide_count.append(tile.count_assigned_targets_guide())
                sky_count.append(tile.count_assigned_targets_sky())
                
                # Count the priorities
                #tile_priorities = tile.count_target_priorities()
                #for level in tile_priorities.keys():
                    #priority_count[level] += tile_priorities[level]
        
        # Tile and Target object
        tiles_by_mag_range.append(tiles_for_range)
        targets_by_mag_range.append(targets_for_range)
        
        # Count Assigned
        target_count_by_range.append(target_count)
        standard_count_by_range.append(standard_count)
        guide_count_by_range.append(guide_count)
        sky_count_by_range.append(sky_count)
        priority_count_by_range.append(priority_count)
        
        # Priorities
        tile_priorities = fwpri.count_unique_priorities(targets_for_range)
        for level in tile_priorities.keys():
            priority_count[level] += tile_priorities[level]
        all_priorities += Counter(priority_count)
        
        tile_count_labels.append("Mag range: %s-%s" % (mrange[0], mrange[1]))
    
    targets_per_tile = run_settings["targets_per_tile"]
    standards_per_tile = run_settings["standards_per_tile"]
    guides_per_tile = run_settings["guides_per_tile"]
    sky_per_tile = run_settings["sky_per_tile"]
    
    # Insert the non-magnitude range specific numbers
    target_count_by_range.insert(0, [t.count_assigned_targets_science() 
                                  for t in tiling]) 
    standard_count_by_range.insert(0, [t.count_assigned_targets_standard() 
                                      for t in tiling]) 
    guide_count_by_range.insert(0, [t.count_assigned_targets_guide() 
                                   for t in tiling])
    sky_count_by_range.insert(0, [t.count_assigned_targets_sky() 
                                for t in tiling])      
    
    priority_count_by_range.insert(0, all_priorities)  
                                
    # Now count the number of unique targets per magnitude range
    unique_targets_range =  []
    for mag_range in tiles_by_mag_range:
        unique, tot, dup = fwtl.count_unique_science_targets(mag_range, True) 
        unique_targets_range.append(unique)                   
    
    # -------------------------------------------------------------------------
    # Plotting
    # -------------------------------------------------------------------------
    name = "results/" + run_settings["run_id"] + "_tiling_run_overview.pdf"
    with PdfPages(name) as pdf:
        # ---------------------------------------------------------------------
        # Page 1
        # ---------------------------------------------------------------------
        plt.clf()
        gs = gridspec.GridSpec(6,6)
        gs.update(wspace=0.2)
        fig = plt.gcf()
        fig.set_size_inches(11.7,8.3) # A4 Paper
        #fig.set_size_inches(28.8, 24.) 
        plt.suptitle("[Page 1] %s: %s" % (run_settings["run_id"], 
                     run_settings["description"]), fontsize=10)
    
        # Tile positions
        ax_tile_pos = plot_tile_pos(fig, gs[:,:], tiling)
        
        pdf.savefig()
        plt.close()
        
        # ---------------------------------------------------------------------
        # Page 2
        # ---------------------------------------------------------------------
        # Assigned target, standard, guide, and sky histograms
        gs = gridspec.GridSpec(5,4)
        gs.update(wspace=0.3)
        fig = plt.gcf()
        fig.set_size_inches(8.3, 11.7) # A4 Paper
        plt.suptitle("[Page 2]", fontsize=10)
        
        # Create plots for all tiles, as well as each magnitude range
        ax_tgt_hists =  plot_target_histogram_grid(fig, gs, [0,0], 
                                                   target_count_by_range, 
                                                   targets_per_tile, 
                                                   unique_targets_range, 
                                                   tiles_by_mag_range, 
                                                   standard_count_by_range, 
                                                   standards_per_tile, 
                                                   guide_count_by_range,
                                                   guides_per_tile, 
                                                   sky_count_by_range, 
                                                   sky_per_tile,
                                                   run_settings, 
                                                   tile_count_labels)
        
        pdf.savefig()
        plt.close()
        
        # ---------------------------------------------------------------------
        # Page 3
        # ---------------------------------------------------------------------
        # Priority histograms
        gs = gridspec.GridSpec(5,3)
        gs.update(wspace=0.3)
        fig = plt.gcf()
        fig.set_size_inches(8.3, 11.7) # A4 Paper
        plt.suptitle("[Page 3]", fontsize=10)        
                         
        ax_priority = plot_priority_histogram(fig, gs, [0,0],
                                              priority_count_by_range, 
                                              run_settings, tile_count_labels)
                                              
        # Table of Run Settings
        # Plot a table for referencing the run settings/results
        ax_tab = plot_settings_table(fig, gs[0:,1:], run_settings)
                                              
        pdf.savefig()
        plt.close()
        
        # ---------------------------------------------------------------------
        # Page 4
        # ---------------------------------------------------------------------
        # All target map
        plt_colours = ["SkyBlue","Gold","Orange", "Tomato"]
        
        fig = plt.gcf()
        fig.set_size_inches(11.7, 8.3)
        ax_target_pos = fig.add_subplot(111, projection='aitoff')
        ax_target_pos.grid(True)
        plt.suptitle("[Page 4]", fontsize=10)   
        
        for range_i, label in enumerate(tile_count_labels[1:]):
            print range_i, label
            ra_deg = [target.ra for target in targets_by_mag_range[range_i]]
            target_ras = np.radians(np.array(ra_deg) - 180)
            dec_deg = [target.dec for target in targets_by_mag_range[range_i]]
            target_decs = np.radians(np.array(dec_deg))

            ax_target_pos.plot(target_ras, target_decs, ".", label=label, 
                               markersize=0.1)
            
            ax_target_pos.set_title("All Science Targets", fontsize=8)
            ax_target_pos.set_xlabel("RA")
            ax_target_pos.set_ylabel("DEC")
            legend = ax_target_pos.legend(loc="best", markerscale=100)
             
        pdf.savefig(fig)
        plt.close()
    
    
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