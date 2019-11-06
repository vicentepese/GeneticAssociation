#!bin/bash

# Go to data
cd ../Data  

# Filter based on maximum and minimum position of the referred study 
#awk '{ if($3 >= 25997458) print $0}' CHR6_new > CHR6_filt
awk '{ if($3 <= 33188223) print $0}' CHR6_filt > CHR6_filtdef
# Check 
#awk '{ if($3 -ge 33188223 || $3 -le 25997458) print $0;}' CHR6_filt > CHR6_check