##Description
The script extract subsequences from a given fasta file by a given
coordinate txt. All extracted subsequences resulting in a output fasta file.

##Requirements
- Python >3.5
--> install instruction for [python](https://docs.python.org/3/using/unix.html#getting-and-installing-the-latest-version-of-python)
- python package Biopython ([biopython.org](http://biopython.org))
pip install biopython

##Run
Options:
-s Source fasta file *(required)*
-c coordinate file *(required)*
-o output fasta *(required)* *(will be created by the script)*

call:
python sub_seq_extrctr.py -s *fasta_file*.fasta -c *coordinate_file*.txt -o *output.fasta*

##Example

*example.fasta*
<pre><code>>Aca_clpP-cp
ATGCCCATTGGTGTTCCAAAAGTACCTTTTCGGAGTCCTGGAGAGGAAGATGCAGCTTGGGTTGACATATAGTGCGACTT
>Cta_clpP-cp
ATGCCCATTGGTGTTCCAAAAGTTCCTTTTCGGATCCCTGGAGAGGAAGATGCAGTTTGGGTCGACGTATAGTGCGACTT
GTCGGACATATCGGGTTATACGGAATTTCCCCGTTATCTCTCCTGATGAGATATTTCCTGCTTCATTCAGGATCGATTCA
</code></pre>

*coordinate_example.txt* 
<pre><code># sequence_name start stop
Ath_clpP-cp 3   8
# comment
Cta_clpP-cp 1   5
</code></pre>

*output.fasta*
<pre><code>>Ath_clpP-cp <fragment start: 3 end: 8>
GCCTAT
>Cta_clpP-cp <fragment start: 1 end: 5>
ATGCC
</code></pre>
