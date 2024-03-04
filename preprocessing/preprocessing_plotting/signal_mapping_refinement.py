#### IMPORTS
import logging
from pathlib import Path
import pod5 as p5
import polars as pl
import plotnine as p9
import remora
from remora import io, refine_signal_map, util
import matplotlib.pyplot as plt

#### Since some package versions required different versions, I could not use EXACTLY the code in the .ipynb notebooks from remora
#### This file is for visualization after signal sequence refinement of reads.
#### https://github.com/nanoporetech/remora/blob/master/notebooks/basic_read_plotting.ipynb

#Read in aligned BAM
bam_fh = io.ReadIndexedBam("aligned_basecalls.bam") 

#### Select a specific read ID
read_id = "0024a438-fcf8-41bd-bd8c-2b9f7601acc3"
bam_read = bam_fh.get_first_alignment(read_id)


#################
###LOAD IN DATA##
#################

with p5.Reader("DEAMERNANOPORE_20161117_FNFAB43577_MN16450_sequencing_run_MA_821_R9_4_NA12878_11_17_16_0.pod5") as reader:

    # Read the selected read from the pod5 file
    # next() is required here as Reader.reads() returns a Generator
    pod5_read = next(reader.reads(selection=[read_id]))
    io_read = io.Read.from_pod5_and_alignment(pod5_read, bam_read)
    print(f"Basecalls length: {io_read.seq_len}")
    print(f"Reference mapping length: {io_read.ref_seq_len}")
    print(f"Reference location: {io_read.ref_reg}")

    # Load in kmer Table (from remora repository, not sure if I am using the right one now)
    level_table = "9mer_levels_v1.txt"
    sig_map_refiner = refine_signal_map.SigMapRefiner(
        kmer_model_filename=level_table,
        do_rough_rescale=True,
        scale_iters=0,
        do_fix_guage=True,
    )

    # Plot the sequence from a read along with the expected signal levels derived from basecalls/reference
    start_base, end_base = 1_000, 1_050
    start_of_basecalls = io_read.extract_basecall_region(
        start_base=start_base, end_base=end_base
    )
    model_levels = sig_map_refiner.extract_levels(util.seq_to_int(io_read.seq))

    # for figure output
    fig, ax = start_of_basecalls.plot_on_base_coords(
        levels=model_levels[start_base:end_base]
    )

    # Save the figure to an image file (e.g., PNG, PDF, SVG)
    fig.savefig('basecall_anchored_expected_signal.png')  # Change the file extension as needed (e.g., 'output_figure.pdf')

    # plot the sequence from a read along with expected signal levels derived from reference sequence
    start_of_mapping = io_read.extract_ref_reg(
    io_read.ref_reg.adjust(end_adjust=50 - io_read.ref_reg.len)
    )
    ref_model_levels = sig_map_refiner.extract_levels(
        util.seq_to_int(io_read.ref_seq)
    )

    # for figure output
    fig, ax = start_of_mapping.plot_on_base_coords(
        levels=model_levels[start_base:end_base]
    )

    # Save the figure to an image file (e.g., PNG, PDF, SVG)
    fig.savefig('reference_anchored_expected_signal.png')  # Change the file extension as needed (e.g., 'output_figure.pdf')

    #############################
    # PERFORM SIGNAL REFINEMENT #
    #############################

    # silence DEBUG messages from signal refinement
    logging.getLogger("Remora").setLevel(logging.INFO)
    #io_read.set_refine_signal_mapping(sig_map_refiner, ref_mapping=False) #perform signal mapping refinement on basecall mapping
    io_read.set_refine_signal_mapping(sig_map_refiner, ref_mapping=True) #perform signal mapping refinement on reference mapping

    start_of_basecalls = io_read.extract_basecall_region(
        start_base=start_base, end_base=end_base)

    # for figure output
    fig, ax = start_of_basecalls.plot_on_base_coords(
        levels=model_levels[start_base:end_base]
    )

    # Save the figure to an image file (e.g., PNG, PDF, SVG)
    fig.savefig('reference_anchored_refined_signal.png')  # Change the file extension as needed (e.g., 'output_figure.pdf')
