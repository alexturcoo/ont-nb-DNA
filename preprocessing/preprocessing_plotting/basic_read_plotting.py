#### IMPORTS
from pathlib import Path
import pod5 as p5
import plotnine as p9
import remora
from remora import io
import numpy as np
import matplotlib.pyplot as plt

#### Since some package versions required different versions, I could not use EXACTLY the code in the .ipynb notebooks from remora
#### This file is for basic visualization of reads.
#### https://github.com/nanoporetech/pod5-file-format/blob/master/python/pod5/README.md
#### https://github.com/nanoporetech/remora/blob/master/notebooks/basic_read_plotting.ipynb

#### Theme set for plotnine - trouble with plotnine right now
#p9.theme_set(p9.theme_minimal() + p9.theme(figure_size=(10, 3)))

#Read in aligned BAM
bam_fh = io.ReadIndexedBam("aligned_basecalls.bam") 


################################
##POD 5 AND BAM INPUTS##########
################################

#### PRINT THE FIRST n Reads in the multi pod5 file
count = 0
with p5.Reader("DEAMERNANOPORE_20161117_FNFAB43577_MN16450_sequencing_run_MA_821_R9_4_NA12878_11_17_16_0.pod5") as reader:
    for read_record in reader.reads():
        print(read_record.read_id)
        count += 1
        if count >= 5:
            break

#### Select a specific read ID from output of previous command
read_id = "0024a438-fcf8-41bd-bd8c-2b9f7601acc3"
bam_read = bam_fh.get_first_alignment(read_id)


################################
### BASECALL-ANCHORED PLOTTING##
################################

with p5.Reader("DEAMERNANOPORE_20161117_FNFAB43577_MN16450_sequencing_run_MA_821_R9_4_NA12878_11_17_16_0.pod5") as reader:

    # Read the selected read from the pod5 file
    # next() is required here as Reader.reads() returns a Generator
    pod5_read = next(reader.reads(selection=[read_id]))
    io_read = io.Read.from_pod5_and_alignment(pod5_read, bam_read)
    print(f"Basecalls length: {io_read.seq_len}")
    print(f"Reference mapping length: {io_read.ref_seq_len}")
    print(f"Reference location: {io_read.ref_reg}")

    # Extract begining and ending of the read - for extracting sections of basecalls
    start_of_basecalls = io_read.extract_basecall_region(end_base=50)
    end_of_basecalls = io_read.extract_basecall_region(
        start_base=io_read.seq_len - 50
    )

    # Assuming start_of_basecalls.plot_on_base_coords() returns a Matplotlib figure
    fig, ax = start_of_basecalls.plot_on_base_coords()

    # Save the figure to an image file (e.g., PNG, PDF, SVG)
    fig.savefig('basecall_anchored.png')  # Change the file extension as needed (e.g., 'output_figure.pdf')

    # Extract beginning and ending - for extracting sections of reference alignment
    start_of_mapping = io_read.extract_ref_reg(
        io_read.ref_reg.adjust(end_adjust=50 - io_read.ref_reg.len)
    )
    end_of_mapping = io_read.extract_ref_reg(
        io_read.ref_reg.adjust(start_adjust=io_read.ref_reg.len - 50)
    )

    # Assuming start_of_basecalls.plot_on_base_coords() returns a Matplotlib figure
    fig, ax = end_of_mapping.plot_on_base_coords()

    # Save the figure to an image file (e.g., PNG, PDF, SVG)
    fig.savefig('reference_anchored.png')  # Change the file extension as needed (e.g., 'output_figure.pdf')
    
    # PLOT A READ'S SIGNAL DATA AGAINST TIME USING JUST POD5 PACKAGE FUNCTIONs (more basic plot than above)
    # Get the signal data and sample rate
    #sample_rate = pod5_read.run_info.sample_rate
    #signal = pod5_read.signal

    # Compute the time steps over the sampling period
    #time = np.arange(len(signal)) / sample_rate

    # Plot using matplotlib
    #plt.plot(time, signal)
    #plt.savefig('test_signal.png')
    

    



