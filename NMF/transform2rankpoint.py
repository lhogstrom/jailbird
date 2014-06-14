#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.5.4/bin/python
# <Script name>
# Author: Angela Brooks
# Program Completion Date:
# Description:
# Modification Date(s):
# Copyright (c) 2011, Angela Brooks. anbrooks@gmail.com
# All rights reserved.


import sys
import optparse 
import os
import pdb

import cmap.io.gct as gct
import cmap

import numpy
#############
# CONSTANTS #
#############

#################
# END CONSTANTS #
#################


###########
# CLASSES #
###########
class OptionParser(optparse.OptionParser):
    """
    Adding a method for required arguments.
    Taken from:
    http://www.python.org/doc/2.3/lib/optparse-extending-examples.html
    """
    def check_required(self, opt):
        option = self.get_option(opt)

        # Assumes the option's 'default' is set to None!
        if getattr(self.values, option.dest) is None:
            print "%s option not supplied" % option
            self.print_help()
            sys.exit(1)


###############
# END CLASSES #
###############
 
########
# MAIN #	
########
def main():
	
    opt_parser = OptionParser()
   
    # Add Options. Required options should have default=None
    opt_parser.add_option("-i",
                          dest="input",
                          type="string",
                          help="""Input GCT(X) file. This should be a pairwise
                                  comparison matrix""",
                          default=None)
    opt_parser.add_option("-o",
                          dest="output",
                          type="string",
                          help="""Output GCTX file where the values are
                                  transformed to a rankpoint matrix. Scale of
                                  -100 to 100""",
                          default=None)

    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("-i")
    opt_parser.check_required("-o")

    gct_object = gct.GCT(options.input) 

    gct_object.read()

    rid = gct_object.get_rids()
    cid = gct_object.get_cids()

    num_sigs = len(rid)

    if len(cid) != num_sigs:
        print "Expecting row and column length to be the same"
        opt_parser.print_help()
        sys.exit(1)       
    
    # empty matrix for transformed dat
    new_mat = []
    
    for raw_vals in gct_object.matrix:
        neg_vals = [0]
        pos_vals = [0]

        for this_val in raw_vals:
            if this_val < 0:
                neg_vals.append(abs(this_val))
            if this_val > 0:    
                pos_vals.append(this_val)

        neg_vals.sort()
        pos_vals.sort()

        len_neg_vals = len(neg_vals)
        len_pos_vals = len(pos_vals) 

        neg_val2rank_tuples = zip(neg_vals, range(len_neg_vals))
        neg_val2rank = dict(neg_val2rank_tuples)

        pos_val2rank_tuples = zip(pos_vals, range(len_pos_vals))
        pos_val2rank = dict(pos_val2rank_tuples)

        a_neg, b_neg = getSlopeIntercept(len_neg_vals)
        a_plus, b_plus = getSlopeIntercept(len_pos_vals)

        transformed_vals = []
        for this_val in raw_vals:
            if this_val == 0.0:
                tranformed_vals.append(0.0)
            elif this_val < 0:
                transformed_vals.append(-(a_neg + b_neg*neg_val2rank[abs(this_val)]))
            else: # this_val > 0:
                transformed_vals.append(a_plus + b_plus*pos_val2rank[this_val])

        new_mat.append(transformed_vals)
               
    # Make new gct object and write
    new_gct_object = gct.GCT()
    new_gct_object.build(numpy.array(new_mat), rid, cid)

    new_gct_object.write(options.output)

			
    sys.exit(0)

############
# END_MAIN #
############

#############
# FUNCTIONS #
#############
def formatDir(i_dir):
    i_dir = os.path.realpath(i_dir)
    if i_dir.endswith("/"):
        i_dir = i_dir.rstrip("/")
    return i_dir

def formatLine(line):
    line = line.replace("\r","")
    line = line.replace("\n","")
    return line

def getSlopeIntercept(num_vals):
    num_vals -= 1 # accounts for 0
    # Translated values will go from 0-100
    b = 100.00/num_vals
    a = 100 - b*num_vals

    return a, b
    
#################
# END FUNCTIONS #	
#################	
if __name__ == "__main__": main()
