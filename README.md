# stitcher.py
[![DOI](https://zenodo.org/badge/258522635.svg)](https://zenodo.org/badge/latestdoi/258522635)

_stitcher.py_ reconstructs molecules from Smart-seq3 data processed with zUMIs => 2.6.0 https://github.com/sdparekh/zUMIs and outputs a .bam file which can be used for further analysis. See the full github page for more information regarding the computational analysis described in the Smart-seq3 article here https://github.com/sandberg-lab/Smart-seq3.

## System Requirements for stitcher.py

_stitcher.py_  is a python3 script with dependencies:
```
pandas
numpy
pysam
joblib
pygtrie
portion
```
No further installation is needed.

## Usage

usage: stitcher.py [-h] [-i input] [-o output] [-g gtf] [-iso iso] [-jun jun]
                   [-t threads] [--single-end] [--cells cells]
                   [--contig contig] [-v]

**arguments:**
```
  -h, --help       show this help message and exit
  --i input        Input .bam file
  --o output       Output .bam file
  --g gtf          gtf file with gene information
  -iso iso, --isoform iso json file with isoform information
  -jun jun, --junction jun json file with exon-exon structure
  --t threads      Number of threads
  --cells cells    List of cell barcodes to stitch molecules (text file, one cell barcode per line).
  --contig contig  Restrict stitching to contig
  -v, --version    show program's version number and exit
```
## Input

_stitcher.py_ takes .bam file processed with zUMIs => 2.6.0 together with a gtf file and a custom json file to reconstruct molecules for the genes in the gtf file.

As optional parameter, the user can specify the number of threads used for parallelization, the cells to process, and restrict the reconstruction to a given contig.

You can find pre-processed .json files here: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4548731.svg)](https://doi.org/10.5281/zenodo.4548731)

## Output 

_stitcher.py_ writes its results to a .bam file as the reads are being processed. Some of the fields have a slightly different interpretation than usual. The algorithm does not handle insertions at the moment, and remove those before stitching the molecule. The query name is in the format "cell:gene:umi". The D character in the CIGAR string indicates missing coverage. The MAPQ is always 255. 

In some cases the UMI reads contain conflicting information regarding the splicing of the molecule. This could either be due to differences in mapping or due to UMI collisions. In those cases, _stitcher.py_ prefer the unspliced option and writes two additional tags which contain the number of bases which have conflicting information and the intervals themselves. 

The .bam file contain many additional custom tags:

```
NR : Number of reads used to stitch.
ER : Number of reads covering an exon.
IR : Number of reads covering an intron.
BC : Cell barcode (same as zUMIs).
XT : Gene barcode (same as zUMIs < 2.6.0).
UB : UMI (same as zUMIs).
EL : Locations of read ends (strand dependent).
NC : Conflict in the reconstruction, the number of conflicting bases.
IL : Conflict in the reconstruction, the intervals where there is a conflict. Written as start1,end1,start2,end2,...
CT : List of transcripts compatible with the molecule
```
If the molecule is entirely in an intron (this may happen due to internal priming events), the CT tag is simply 'intronic'. Additionally, if there is no transcript consistent with the positions the molecule cover, the CT tag is an empty string.

## Example 

```
python3 stitcher.py --i smartseq3_file.bam --o smartseq3_molecules.bam --g Mus_musculus.GRCm38.91.chr.clean.gtf --isoform mm10_unique_intervals_for_isoforms.json --t 10 --contig chr1 --cells cells.txt
```
## gtf_to_json.py

_gtf_to_json.py_ is a helper script which takes the gtf file you are using as an input and writes the json file you need for _stitcher.py_ as output. There is a database file written as an intermediary file required by the gffutils package. Use this script if you have a custom gtf file or want to be extra careful.

You can find pre-processed .json files here: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4548731.svg)](https://doi.org/10.5281/zenodo.4548731)


### Usage 
gtf_to_json.py [-h] [-g gtf] [-d db] [-j json] [-t threads]
**arguments:**
```
  -h, --help            show this help message and exit
  -g gtf, --gtf gtf     Input gtf file
  -d db, --db db        Intermediary database (db) file
  -ji json, --json_intervals json Output json file for coverage
  -jr json, --json_refskip json Output json file for refskip
  -t threads, --threads threads Number of threads
```
### Example
```
python3 gtf_to_json.py -g Mus_musculus.GRCm38.91.chr.clean.gtf -d Mus_musculus.GRCm38.91.chr.clean.db -ji Mus_musculus.GRCm38.91.chr.clean.intervals.json -jr Mus_musculus.GRCm38.91.chr.clean.refskip.json -t 5

```
