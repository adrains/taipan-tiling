"""Module for generating FunnelWeb's list of priorities associated with each 
science target.
"""
from astropy.io import fits
import numpy as np
import pylab as pl
from astropy.table import Table
import random
from collections import Counter

def count_unique_priorities(target_list):
    """Counts the number of targets at each priority level for a set of tiles.
    
    Returns
    -------
    priority_count: Counter
        collections.Counter object, containing a dictionary with keys of
        each priority level, and values corresponding to each time it 
        appears.
    """
    target_list = set(target_list)
    target_priorities = [target.priority for target in target_list]
    priority_count = Counter(target_priorities)
    
    return priority_count
    

def gen_normal_priorities(normal, file="fw_input_catalogue.fits"):
    """Function to generate a full set of uniform priorities (i.e. no star is
    prioritised higher than any other).
    
    Parameters
    ----------
    normal: int
        The normal (i.e. default) priority for an unobserved star.
    file: string
        Pathway to the survey input catalogue. 
    """
    # Load fits file and two requisite tables
    fits_file = fits.open(file, mode="update")
    fw_input_cat = fits_file[1].data

    # Get the IDs from the input catalogue
    fw_ids = [row[0] for row in fw_input_cat]

    # Create a default list of priorities
    priorities = normal * np.ones(len(fw_ids)).astype('int')

    # A note on formats:
    # J --> 32 bit integer
    # K --> 64 bit integer
    id_col = fits.ColDefs([fits.Column(name="Gaia_ID", format="K", 
                                       array=fw_ids)])
    priority_col = fits.ColDefs([fits.Column(name="Priority", format="J", 
                                             array=priorities)])

    hdu = fits.BinTableHDU.from_columns(id_col + priority_col)

    # Save the results            
    hdu.writeto("fw_priorities.fits")

    # Close
    fits_file.close()


def create_reduced_input_catalogue(file="fw_input_catalogue.fits", frac=100):
    """
    """
    # Load fits file and two requisite tables
    fits_file = fits.open(file, mode="update")
    fw_input_cat = fits_file[1].data
    
    reduced_input_cat = fw_input_cat[::frac]
    
    new_name = file.replace(".fits", "_reduced_x%i.fits" % frac)

    hdu = fits.BinTableHDU(data=reduced_input_cat)
    hdu.writeto(new_name)

    # Close
    fits_file.close()

    


def gen_stat_priorities(normal=2, tabdata=None, stars_per_level=[50,50,50],
                        file="fw_input_catalogue.fits"):
    """Code to generate priorities for the FunnelWeb input catalogue, 
    assigning N stars per elevated priority, and all others to normal.
    
    Non-elevated priority levels for FunnelWeb are:
        0 - Desperate
        1 - Marginal
        2 - Normal
        3 - Preferred
        4 - Required
        5 - Urgent
        
    With 0 and 1 being for stars that have already been observed. This script
    assigns the priorities as follows:
        0 - 0 stars
        1 - 0 stars
        2 - Total - 3N stars
        3 - N stars
        4 - N stars
        5 - N stars
        
    The purpose of this is to facilitate simulations on FunnelWeb's 
    completeness as a function of priority level and magnitude range. By
    selecting 3N stars randomly, we each magnitude range should be 
    proportionally represented.
    
    Note: Code is *hella* slow - thankfully it doesn't need to be run often!
    
    Parameters
    ----------
    normal: int
        The normal (i.e. default) priority for an unobserved star.
    tabdata: fits.io.table
        Survey input catalogue in fits form.
    file: string
        Pathway to the survey input catalogue.
        
    Returns
    -------
    fw_ids: list
        A complete list of all IDs from the input catalogue.
    priorities: list
        A complete list, with corresponding order, of target priorities.
    priorities_main: list
        List of priorities associated with the main FunnelWeb survey.
    fw_main_coords: list
        Coordinates in RA and DEC of main survey targets, order corresponding
        to priorities_main.        
    """
    # Priorities - ignore 0 and 1 as these will be reserved for already
    # observed targets
    priority_levels = [3, 4, 5]
    
    # Limits
    ra_min = 0
    ra_max = 360
    dec_min = -90
    dec_max = 0
    gal_lat_limit = 10

    # Load fits file
    if not tabdata:
        print "Loading survey input catalogue..."
        tabdata = Table.read(file)
    
    # Initialise
    fw_ids_main = []
    fw_main_coords = []
    fw_ids_other = []
    
    # Collect all the stars from the main survey (i.e. satisfying the above 
    # limits) in one list, and everything else in another 
    for star in tabdata:
        # Only consider targets which satisfy RA/DEC/b restrictions
        if ((ra_min < star["RA_ep2015"] < ra_max)
            and (dec_min < star["Dec_ep2015"] < dec_max)
            and (np.abs(star['b']) > gal_lat_limit) 
            and star["Gaia_G_mag"] <= 30):
            # Target is acceptable
            fw_ids_main.append(star["Gaia_ID"])
            fw_main_coords.append([star["RA_ep2015"], star["Dec_ep2015"]])
            
        else:
            # Add to total list
            fw_ids_other.append(star["Gaia_ID"])

    # Generate a random selection of *main survey* indices
    num_upvoted = sum(stars_per_level)
    random_indices = random.sample(range(len(fw_ids_main)), num_upvoted)
    assert len(random_indices) == len(set(random_indices))
    
    # Initialise the priorities for the main survey to be the default
    priorities_main = normal * np.ones(len(fw_ids_main)).astype('int')
    
    # For each priority level, assign the associated initial priority to 
    # stars_per_level stars using the randomly selected IDs
    for level, priority in enumerate(priority_levels):
        range_start = sum(stars_per_level[:level])
        range_end = sum(stars_per_level[:level+1])

        for star_i in random_indices[range_start:range_end]:
            priorities_main[star_i] = priority
        
    # Create a default list of priorities for the remaining stars
    priorities_other = normal * np.ones(len(fw_ids_other)).astype('int')

    # Join the two separate lists and save
    fw_ids = fw_ids_main + fw_ids_other
    priorities = np.append(priorities_main, priorities_other)

    # Save the results 
    # A note on formats:
    # J --> 32 bit integer
    # K --> 64 bit integer
    id_col = fits.ColDefs([fits.Column(name="Gaia_ID", format="K", 
                                       array=fw_ids)])
    priority_col = fits.ColDefs([fits.Column(name="Priority", format="J", 
                                             array=priorities)])

    hdu = fits.BinTableHDU.from_columns(id_col + priority_col)

    string_rep = ""
    for level in stars_per_level: 
        string_rep += str(level) + "-"

    hdu.writeto("fw_stat_%s_priorities.fits" % string_rep[:-1], overwrite=True)
    
    # Plot the coordinates of the high priority stars for sanity checking
    """
    pl.figure()
    pl.subplot(projection="aitoff")
    pl.xlabel("RA")
    pl.ylabel("DEC")
    pl.grid(True)
    
    # Ensure each priority range has a different set of points
    points = ["rx", "b.", "go"]
    
    for point, priority in enumerate(priority_levels):
        ra = [fw_main_coords[i][0] for i, x in enumerate(priorities_main) 
              if priorities_main[i] == priority] 
        dec = [fw_main_coords[i][1] for i, x in enumerate(priorities_main) 
               if priorities_main[i] == priority]    

        pl.plot(ra, dec, points[point], label=priority)
        
    pl.legend(loc="best")
    """
    return fw_ids, priorities, priorities_main, fw_main_coords