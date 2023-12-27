# ATG_CRISPR_Score_Calculator

This program calculates CFD and MIT scores for SpyCas9 off-targets. Written specifically in the C++ language to achieve optimal performance. Calculations are adjusted to allow scoring off-targets with RNA/DNA bulges, in particular, to allow scoring CasOffinder-bulge output.


### Background:
 **CFD**: Cutting Frequence Determination (CFD) scores range between 0 to 1 and represent the percentage activity of *Streptococcus pyogenes* Cas9 (SpyCas9) and its CRISPR guide RNA (crRNA) at a specific off-target sequence. CFD scoring system take number, position, and type and mismatches between crRNA::DNA into account. The final score is adjustment for non-conical PAMs.

**MIT**: takes a value between 0 to 1 and represents the likelihood  of a CRISPR guide RNA (crRNA) cutting a specific off-target sequence. A higher MIT Score suggests higher off-target editing possibility. MIT scores takes number, position, and distance of mismatches between crRNA::DNA into account. Adjustment for non-conical PAMs is done using CFD Scores for non-conical PAMs.

 ### Command line arguments:

 You can run this program in two modes: `singlePair` and `pairList`.

 ### Usage: 
**singlePair**: `./ATG_CRISPR_Score_Calculator singlePair <crRNASequence> <DNASequence> `

Example:
```
./ATG_CRISPR_Score_Calculator singlePair ATCGATGCTGATGCTAGATAAGG ACCGATGCTGATGCTAGATAAGG 

output:
CFD     0.857143        MIT     1
```

**pairList**: `ATG_CRISPR_Score_Calculator pairList <pairListFilePath>`
Example: 

```
./ATG_CRISPR_Score_Calculator pairList data/crRNA_DNA_Pair_Examples.txt

output:
crRNA   DNA     CFD     MIT
ATCGATGCTGATGCTAGATAAGG ATCGATGCTGATGCTAGATAAGG 1       1
ATCGATGCTGATGCTAGATAAGG AcCGATGCTGATGCTAGATAAGG 0.857143        1
ATCGATGCTGATGCTAGATAAGG ATCGATGCTGATGCTAGATtAGG 0.6     0.417
ATCGATGCTGATGCTAGATAAGG ATCGATGCTGATGCTAGATAAGa 0.0694444       0.0694444
ATCGATGCTGATGCTAGATAAGG A-CGATGCTGATGCTAGATAAGG 0.953125        1
AT-CGATGCTGATGCTAGATAAGG        ATaCGATGCTGATGCTAGATAAGG        0.727273        1
```


Note: crRNA and DNA squences must be at least 23 nucleotides or longer and end with their PAMs. 

### Requirements: 
 - **CFD_Scores.txt**: A valid cfdScoreLine in the `CFD_Scores.txt` file is composed of a Label as string and a decimal number as double, in a tab-delimited format (e.g.,: `rU:dG_2   0.857143`). Labels should also include the CFD scores for PAMs. See `data/CFD_Scores.txt` file for details. 
- **crRNA** sequence: case-insensitive letters: A,T/U,C,G,-, 23 nucleotides or longer, ending with PAM. 
- **DNA** sequence: case-insensitive letters: A,T,C,G,-, 23 nucleotides or longer (but must have equal length to crRNA), ending with PAM. 

### Installation
```
git clone https://github.com/Ghahfarokhi/ATG_CRISPR_Score_Calculator.git
cd ATG_CRISPR_Score_Calculator
make ATG_CRISPR_Score_Calculator
./ATG_CRISPR_Score_Calculator -h
```

**References**: 
- **CFD**: Doench, J., et al. (2016) Nat Biotechnol 34, 184–191. 
- **MIT**: Hsu, P., et al. (2013) Nat Biotechnol 31, 827–832.  


*Note*: CFD Scores in the original python code published by Doench J. et al., contained calculation for mismatches only (RNA/DNA bulges were not included). Here, I extracted the CFD scores for DNA and RNA bulges from Supplementary Table 19 and listed them in the `data/CFD_Scores.txt` file. However, it has to be used with cautions as the combination of mismatches and bulges CFD scores is not peer-reviewed, and is vaguely described. The output is identical to the original python code if the provided crRNA and DNA don't contain RNA/DNA bulges. 

**Developer**: Amir.Taheri.Ghahfarokhi@Gmail.com 

**Initial release**: 2023-12-27