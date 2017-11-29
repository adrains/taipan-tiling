"""Quick script to test assigning sky fibres to dark sky coordinates.
"""
import cPickle
import numpy as np
import taipan.core as tp
import taipan.tiling as tl
import taipan.fwtiling as fwtl
from astropy.table import Table
import funnelweb_tiling_settings as fwts
import matplotlib.pylab as plt

# Ensure taipan.core defaults are updated for FW
tp.SKY_PER_TILE = fwts.script_settings["SKY_PER_TILE"]
tp.SKY_PER_TILE_MIN = fwts.script_settings["SKY_PER_TILE_MIN"]
tp.QUAD_RADII = fwts.script_settings["QUAD_RADII"]
tp.QUAD_PER_RADII = fwts.script_settings["QUAD_PER_RADII"]

if len(tp.QUAD_PER_RADII) != len(tp.QUAD_RADII)-1:
    raise ValueError('QUAD_PER_RADII must have one less element than '
                 'QUAD_RADII')
if sum(tp.QUAD_PER_RADII) != tp.SKY_PER_TILE:
    raise UserWarning('The number of defined tile quandrants does not match '
                      'SKY_PER_TILE. These are meant to be the same!')

tp.FIBRES_PER_QUAD = []
for i in range(len(tp.QUAD_RADII))[:-1]:
    theta = 360. / tp.QUAD_PER_RADII[i]
    for j in range(tp.QUAD_PER_RADII[i]):
        tp.FIBRES_PER_QUAD.append(
            [k for k in tp.BUGPOS_OFFSET.keys() if
             k not in tp.FIBRES_GUIDE and
             tp.QUAD_RADII[i+1] <= tp.BUGPOS_OFFSET[k][0] < tp.QUAD_RADII[i] and
             j*theta <= tp.BUGPOS_OFFSET[k][1] < (j+1)*theta])

def visualise_tile(tile):
    """
    """
    sky_fibres = [fibre for fibre in tile.fibres if tile.fibres[fibre] == "sky"]
    
    sky_fibre_pos = [tp.BUGPOS_MM[fibre] for fibre in sky_fibres]
    
    for fibre in sky_fibre_pos:
        plt.plot(fibre[0]*tp.ARCSEC_PER_MM, fibre[1]*tp.ARCSEC_PER_MM, 
                 ".", label="Sky")
        patrol_radius = plt.Circle((fibre[0]*tp.ARCSEC_PER_MM, 
                                 fibre[1]*tp.ARCSEC_PER_MM), 
                                 tp.PATROL_RADIUS, color="r", alpha=0.2)
        plt.gcf().gca().add_artist(patrol_radius)
    
    # Plot the field
    field = plt.Circle((0,0), tp.TILE_RADIUS, 
                                 color="b", alpha=0.2)
    plt.gcf().gca().add_artist(field)


# Import the results of the tiling run
print "Importing tiling run"
pkl_file_name = "results/172711_1533_46_fw_tiling.pkl"

pkl_file = open(pkl_file_name, "rb")
(tiling, remaining_targets, run_settings) = cPickle.load(pkl_file)

try:
    if dark_coords:
        print "Already imported dark coordinates"
except:
    # Import the sky catalogue and generate sky "targets"
    print "Importing dark sky catalogue"
    dark_sky_fits = "/Users/adamrains/Catalogues/skyfibers_v17_gaia_ucac4_final.fits"
    
    dark_sky_table = Table.read(dark_sky_fits)
    
    dark_coords = [tp.FWTarget(-1*int(star["pkey_id"]), star["ra"], star["dec"], 
                   priority=run_settings["priority_normal"], mag=25, difficulty=1, 
                   science=False, standard=False, guide=False, sky=True) 
                   for star in dark_sky_table 
                   if (run_settings["ra_min"] < star["ra"] < run_settings["ra_max"])
                   and (run_settings["dec_min"] < star["dec"] < run_settings["dec_max"])]
                   
    burn = [t.compute_usposn() for t in dark_coords]
    tp.compute_target_difficulties(dark_coords)

all_removed_targets = []

# All imported, now test assigning sky fibres - test with first tile
for tile_i, tile in enumerate(tiling):
    # Gather sky fibres
    sky_fibres = [fibre for fibre in tile.fibres if tile.fibres[fibre] == "sky"]

    # Get the local dark sky coordinates
    nearby_dark_coords = tp.targets_in_range(tile.ra, tile.dec, dark_coords, 
                                             1*tp.TILE_RADIUS) 
    
    successful_assignments = 0
    
    print "Assigning fibres for tile %4i (%4i total coords)" % (tile_i, 
                                                                len(nearby_dark_coords)),
    """
    # Assign sky fibres
    for fibre in sky_fibres:
        # Note that if a sky fibre fails to be assigned, it was reset to empty and has
        # lost any memory of being a "sky" fibre
        
        remaining_coords, former_tgt = tile.assign_fibre(fibre, nearby_dark_coords, 
                              check_patrol_radius=True, 
                              check_tile_radius=run_settings["check_tile_radius"],
                              recompute_difficulty=run_settings["recompute_difficulty"],
                              order_closest_secondary=True,
                              method=run_settings["tile_unpick_method"],
                              combined_weight=run_settings["combined_weight"],
                              sequential_ordering=(2,1,0))
                              
        if len(remaining_coords) < len(nearby_dark_coords):
            successful_assignments += 1  
            
        print "%2i/%i sky fibres successfully assigned" % (successful_assignments, 
                                                       len(sky_fibres))                        
        """
    # Test sky assigning algorithm (modified from guide assignment)
    removed_targets = tile.assign_sky_fibres(nearby_dark_coords, 
                          target_method=run_settings["tile_unpick_method"], 
                          combined_weight=run_settings["combined_weight"],
                          sequential_ordering=run_settings["sequential_ordering"],
                          check_tile_radius=True, rank_sky=False)
                          
    all_removed_targets.extend(removed_targets)
                              
    print "%2i sky fibres successfully assigned" % tile.count_assigned_targets_sky(),
    print "with %2i science fibres removed" % len(removed_targets)  
    
    