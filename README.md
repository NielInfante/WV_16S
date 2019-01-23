# 16S Analysis
Templates for 16S analysis


These are a series of R scripts to help with 16S analysis using Dada2. These are not ready to run, you'll have to make changes to directories and condition names.

Steps that have been done are marked with a  :white_check_mark: while those that have yet to be started are marking with a  :x:. Steps in progress are marked by a :large_orange_diamond:.

## Steps

First, set up your directories. There should be sub directories called Data, filtered, deseq, vegan, and pics. In Data there should be the trained taxonomic data.

### :x: <!--:large_orange_diamond: :white_check_mark:--> Start with doDada.R

* Point to wwhere the raw files are
* There should be a file called meta that contains the meta data
* You may need to do a bit of work to get the meta data to agree with the file names
* Look at the quality and change parameters in the filterAndTrim step
* Make sure error rates look good. If not, change the filtering parameters
* Set min and max length of merged sequence
* Make sure the taxonomy files are in place, and pick old or new assign taxonomy

This will create a few QC figures, and save filtered files, seqtab, taxa, and ultimately a phyloseq object.


### :x: <!--:large_orange_diamond: :white_check_mark:-->  Do Filtering with doPhyloseq.R

* Set working directory
* Step through to see how you want to filter data

Output a filtered phyloseq object

### :x: <!--:large_orange_diamond: :white_check_mark:--> Do differential abundance with doDeseq.R

* Set working directory
* Set output directory
* Filter the samples you want
* Set the output prefix, intGroup, contrast
* If there are too many zeros, you may need to use the alternate code to load the dds object

Outputs many pictures, and a list of significantly different taxa

### :x: <!--:large_orange_diamond: :white_check_mark:--> Do Beta diversity with doVeganBeta.R

### :x: <!--:large_orange_diamond: :white_check_mark:--> Do Alpha analysis with doAlpha.R

### :x: <!--:large_orange_diamond: :white_check_mark:--> Make pretty pictures with drawPics.R

