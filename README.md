GCcalc
======

Calculate circos format data of GC content and GC skew from genome data. 


Usage
=====

Usage: GCcalc.py -f input.fa [-w 1000] [-s 1000]

Options:
				-h, --help            show this help message and exit

				-f FASTA, --file=FASTA
                        Input Fasta format file
  			-w WINDOWSIZE, --window=WINDOWSIZE
                        default:1000 WindowSize to calculate
  			-s STEPSIZE, --step=STEPSIZE
                        default:1000 StepSize for slide widows
