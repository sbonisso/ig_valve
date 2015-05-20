ig_valve
=============
IgVALVE - Ig Vdj Antibody Labeling Validation and Evaluation

Validates antibody VDJ labeling of reads in either supervised (to simulated ground truth) or unsupervised (to real data) modes. Will measure:

* gene-segments
* allelic variants
* V/D/J partitions
* CDR3 distribution
* clone distribution

##### Use

valve executable contains two subcommands: unsupervised and supervised

```
$ ig_valve help
Commands:
  ig_valve help [COMMAND]          # Describe available commands or one specific...
  ig_valve supervised [options]    # compare predictions to simulated dataset wi...
  ig_valve unsupervised [options]  # compare predictions on dataset without labels
```

supervised requires a ground truth file, and a single prediction file, or a file with paths to multiple prediction files. 
unsupervised requires a file with paths to multiple prediction files.

Both subcommands will require a type of comparison to perform (-y), one of gene/allele/partition/cdr3/clone.
 
 