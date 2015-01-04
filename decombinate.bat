REM this is the directory of the files to be processed. However, note that the files must ALSO be inside Decombinator folder.
REM set directory="C:\New folder\test"
REM set directory ="C:\Users\Benny Chain\Dropbox\R\10_12_2014\output"

REM cycle though this file and DECOMBINATE
REM note the Decombinator commands
REM -i is input file -o is output folder
REM -b True means extract barcode from bs1 to be1; and from bs2 to be2. Note that first base pair is 0, NOT 1

FOR /f %%A IN ('DIR /B *.fastq')DO decombinatorV2_2.py -i %%A  -o %%A -b True -bs1 0 -be1 6 -bs2 6 -be2 6 -p True -c True -of True -f True
