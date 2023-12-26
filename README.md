# CFD Scores Calculator in C++ 
 CFD (Cutting Frequency Determination) score for a given CRISPR guide RNA (crRNA) and DNA sequence is used to estimate the likelihood of a particular guide RNA inducing a double-strand break at a specific genomic targeted by *Streptococcus pyogenes* Cas9 (SpyCas9).

 ### Command line arguments:
      
    arg_1: cfdScores_file. 
    arg_2: crRNA sequence (case-insensitive letters: A,T/U,C,G,-), with PAM.
    arg_3: DNA sequence (case-insensitive letters: A,T,C,G,-), with PAM.

Note: crRNA and DNA squences must be at least 23 nucleotides or longer and end with their PAMs. 

### Requirements: 
 - **CFD_Scores.txt**: A valid cfdScoreLine in the CFD_Scores.txt file is composed of a Label as string and a decimal number as double, in a tab-delimited format (e.g.,: `rU:dG_2   0.857143`). Labels should also include the CFD scores for PAMs. See `data/CFD_Scores.txt` file for details. 

### Installation
```
git clone https://github.com/Ghahfarokhi/CFD_Scores_Calculator.git
cd CFD_Scores_Calculator
make cfdScoresCalculator
```

### Usage

**Example usage**:
```
./cfdScoresCalculator data/CFD_Scores.txt ATCGATGCTGATGCTAGATAAGG ACCGATGCTGATGCTAGATAAGG 
```

**Output**: 0.857143

**Reference**: 
- Doench, J., et al. (2016) Nat Biotechnol 34, 184–191.
- *Note*: CFD Scores in the original python code published by Doench J. et al., contained calculation for mismatches only (RNA/DNA bulges were not included). Here, I extracted the CFD scores for DNA and RNA bulges from Supplementary Table 19 and listed them in the `data/CFD_Scores.txt` file. However, it has to be used with cautions as the combination of mismatches and bulges CFD scores is not peer-reviewed, and is vaguely described. The output is identical to the original python code if the provided crRNA and DNA don't contain RNA/DNA bulges. 

**Developer**: Amir.Taheri.Ghahfarokhi@Gmail.com 

**Initial release**: 2023-12-25