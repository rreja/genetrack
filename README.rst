GeneTrack
=========

Introduction
-------------

Genetrack is a generic, strand independent peak-detection algorithm implemented in C programming language.
It may be used to predict peaks on any data that is represented as measurements on
each strand of a genome.

Installation
------------

We have binary versions available at https://github.com/bcclib/genetrack/bin .
Download the version for your platform and run it. Alternatively you may
recompile from source code::

    $ git clone git@github.com:bcclib/genetrack.git
    $ cd genetrack
    $ make test

Detailed usage
---------------

When run from command line the program prints help::

    $ genetrack

The most common invocation is that to predict peaks, see::

    $ genetrack peaks -h
    
In the simplest case the the program is invoked as::

    $ genetrack peaks inputfile.gff > peaks.gff
    
Several options may be use to tune the peak prediction parameters::

    Usage: ./genetrack peaks [options] input.gff > output.gff
    
    Options:
        -w <wigprefix>  also generate coverage files in wiggle format
        -s <number>     kernel width (smoothing) parameter, default = 5
        -x <integer>    minimum peak to peak distance in basepairs, default = 10
        -m <number>     minimum peak height, default = 2.5
        -b              input files in BED format
        -d              generate more messages
    
    Note: if no input file is specified the program will read the standard input

The following options are of interest

- -w <wigprefix> : Also produce the output in WIG_ format. This option also generates two wiggle files for each strand.
  The wiggle files contain the smoothed coverage that was used to detect peaks.

- -s <number> : Bandwith (smoothing) parameter. Typically s=5 is used for fine-grain peak calls and s=20
  for coarse-grain peak calls.

- -x <number> : Minimum peak-to-peak distances. It may be used to forbid calling peaks that are too close too close
  to one another. The predicted peaks will also be of this width around the midpoint. Note that
  the smoothing parameter also affects the distance over which a second peak may form.

- -m <number> : Minimum peak height for a peak to be predicted. Default is 2.5, the program will not output
  peaks smaller than this value.

- -b : if the input file is in BED_ format. The output will still be produced as GFF_

- -d : generates more messages during processing

Examples
--------

Run examples::

    $ genetrack peaks input.gff > peaks.gff
    
equivalent to::

    $ genetrack peaks < input.gff > peaks.gff
    
Or with options::

    $ genetrack -w out -x 10 -m 0 -f input.gff > peaks.gff

Notes
-----

The peak height can be thought of as a moving average of the reads(tags) that fall within a moving window.
The readcounts reported by GeneTrack is the number of reads that fall within this window. By default the
window size is ``3 * sigma`` the factor ``3`` may be changed via the ``-f`` parameter.
For example if ``sigma=5`` then the readcounts reported in the output will be computed over the range of
``(midpoint - 15, midpoint + 15)`` of the peak. Note that this range is unrelated to do with the size of
the exclusion zone or the width of the peak itself!

In ideal conditions where all reads indicate the same binding location the peak height will be equal to the readcount for
any (realistic) sigma value. The measurement errors or translational settings of the binding will reduce the peak
height while keeping the read count identical. We recommend  using the difference between the readcount and peak height
as an indicator of the ``fuzziness`` or uncertainty of the binding event. This value is more robust to parameter changes
than most other measures.

.. _GFF: http://genome.ucsc.edu/FAQ/FAQformat#format3
.. _BED: http://genome.ucsc.edu/FAQ/FAQformat.html#format1
.. _WIG: http://genome.ucsc.edu/FAQ/FAQformat.html#format6
