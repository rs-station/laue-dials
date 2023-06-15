#!/usr/bin/env python
"""
This script combines integrated MTZ files by image into one single MTZ file.
"""

import sys
import logging
import time

import libtbx.phil
import reciprocalspaceship as rs
from dials.util import show_mail_handle_errors
from dials.util.options import ArgumentParser
from tqdm import trange

logger = logging.getLogger("laue-dials.command_line.combine_mtzs")

help_message = """

This program takes a set of MTZ files and combines them into a single output
MTZ file.

The output is an MTZ file consisting of all the data in the input MTZ files
suitable for merging in software such as careless.

Examples::

    laue.combine_mtzs [options] *.mtz
"""

# Set the phil scope
phil_scope = libtbx.phil.parse(
    """

output {
  filename = 'total_integrated.mtz'
    .type = str
    .help = "The output MTZ filename."

  log = 'laue.combine_mtzs.log'
    .type = str
    .help = "The log filename."
  }
""",
    process_includes=True,
)

working_phil = phil_scope.fetch(sources=[phil_scope])


@show_mail_handle_errors()
def run(args=None, *, phil=working_phil):
    # Parse arguments
    usage = "laue.combine_mtzs [options] *.mtz"

    parser = ArgumentParser(
        usage=usage,
        phil=phil,
        read_reflections=False,
        read_experiments=False,
        check_format=False,
        epilog=help_message,
    )

    params, options, files = parser.parse_args(
        args=args, show_diff_phil=True, ignore_unhandled=True, return_unhandled=True
    )

    # Configure logging                                                         
    console = logging.StreamHandler(sys.stdout)                                 
    fh = logging.FileHandler(params.output.log, mode="w", encoding="utf-8")
    loglevel = logging.INFO                                                     
                                                                                
    logger.addHandler(fh)                                                       
    logger.addHandler(console)                                                  
    logging.captureWarnings(True)                                               
    warning_logger = logging.getLogger("py.warnings")                           
    warning_logger.addHandler(fh)                                               
    warning_logger.addHandler(console)                                          
    dials_logger = logging.getLogger("dials")                                   
    dials_logger.addHandler(fh)                                                 
    dials_logger.addHandler(console)                                            
    dxtbx_logger = logging.getLogger("dxtbx")                                   
    dxtbx_logger.addHandler(fh)                                                 
    dxtbx_logger.addHandler(console)                                            
    xfel_logger = logging.getLogger("xfel")                                     
    xfel_logger.addHandler(fh)                                                  
    xfel_logger.addHandler(console)                                             
                                                                                
    logger.setLevel(loglevel)                                                   
    dials_logger.setLevel(loglevel)                                             
    dxtbx_logger.setLevel(loglevel)                                             
    xfel_logger.setLevel(loglevel)                                              
    fh.setLevel(loglevel)                                                       

    # Log diff phil
    diff_phil = parser.diff_phil.as_str()
    if diff_phil != "":
        logger.info("The following parameters have been modified:\n")
        logger.info(diff_phil)

    # Print help if no input                                        
    if len(files) == 0:
        parser.print_help()                                         
        return                                                      

    # Get initial time for process
    start_time = time.time()      

    # Loop over input files
    logger.info("Beginning combination.")
    total_integrated_mtz = rs.read_mtz(files[0])  # First still
    for i in trange(1, len(files)):
        try:
            # Get MTZ data
            img_mtz = rs.read_mtz(files[i])
            img_mtz["BATCH"] = i + 1

            # Add MTZ data to total
            total_integrated_mtz = rs.concat([total_integrated_mtz, img_mtz])
        except:
            # Print integration failure to user
            i + 1
            logger.info("Image {j:06d} could not be integrated.")
            continue

    # Fix data type of BATCH
    total_integrated_mtz["BATCH"] = total_integrated_mtz["BATCH"].infer_mtz_dtype()

    # Save reflections
    logger.info("Saving combined MTZ data to %s", params.output.filename)
    total_integrated_mtz.write_mtz(params.output.filename)

    # Final logs                                                                
    logger.info("")                                                             
    logger.info(                                                                
        "Time Taken for Total Processing = %f seconds", time.time() - start_time
    )                                                                           

if __name__ == "__main__":
    run()
