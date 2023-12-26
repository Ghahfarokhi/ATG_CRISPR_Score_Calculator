#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>

// Function to get the reverse complement of a DNA sequence
std::string revComp(const std::string& DNASequence) {
    std::string complement;
    for (char base : DNASequence) {
        switch (base) {
            case 'A':
                complement += 'T';
                break;
            case 'T':
                complement += 'A';
                break;
            case 'C':
                complement += 'G';
                break;
            case 'G':
                complement += 'C';
                break;
            default:
                complement += base;
        }
    }
    return complement;
}

void printHelp(){
    const char *helpText = "\nThis program takes the following command line arguments:\n\
      arg_1: cfdScores_file. \n\
      arg_2: crRNA sequence (case-insensitive letters: A,T/U,C,G,-), with PAM.\n\
      arg_3: DNA sequence (case-insensitive letters: A,T,C,G,-), with PAM.\n\n\
      Note: crRNA and DNA squences must be at least 23 nucleotides or longer and end with their PAMs. \n\n\
    Usage example: \n\
    cfdScores CFD_Scores.txt ATCGATGCTGATGCTAGATAAGG ACCGATGCTGATGCTAGATAAGG \n\n\
    A valid cfdScoreLine in the CFD_Scores.txt file is composed of a Label as string and a decimal number as double, tab-delimited. \n\
    Example: rA:dA_10	0.882352941 \n\
    Labels should also include the CFD scores for PAMs. \n\n\
    CFD Score Reference: Doench, J., et al. (2016) Nat Biotechnol 34, 184–191. \n\n\
    Developer: Amir.Taheri.Ghahfarokhi@Gmail.com \n\
    Release: 2023-12-25\n\n";
    std::cout << helpText << std::endl;
}

int main(int argc, char *argv[]) {
    if (argc == 2 && std::string(argv[1]) == "-h") {
        printHelp();
        return 0;
    }
    
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <cfdScores_file> <crRNA_sequence> <DNA_sequence>" << std::endl;
        return 1;
    }

    std::string cfdScoresFilePath(argv[1]);
    std::ifstream cfdScoresFile(cfdScoresFilePath);

    if (!cfdScoresFile.is_open()) {
        std::cerr << "Error opening file: " << cfdScoresFilePath << std::endl;
        return 1;
    }

    std::map<std::string, double> cfdScores;
    std::string cfdScoreLine;

    while (std::getline(cfdScoresFile, cfdScoreLine)) {
        // Skip cfdScoreLines that start with "Label" or other headers
        if (cfdScoreLine.find("Label") != std::string::npos || cfdScoreLine.find("Error") != std::string::npos) {
            continue;
        }

        std::istringstream iss(cfdScoreLine);
        std::string label;
        double cfd;
        if (iss >> label >> cfd) {
            cfdScores[label] = cfd;
        } else {
            std::cerr << "Error parsing cfdScoreLine: " << cfdScoreLine << std::endl;
        }
    }

    // Get crRNA and DNA sequences from command cfdScoreLine arguments
    std::string crRNASequence(argv[2]);
    std::string DNASequence(argv[3]);

    // Extract the rightmost 2 characters of the DNA as PAM
    std::string PAM = DNASequence.substr(DNASequence.length() - 2);

    // Update crRNA and DNA with the exact 20 characters upstream of the last 3 characters
    if (crRNASequence.length() >= 23 && DNASequence.length() >= 23) {
        crRNASequence = crRNASequence.substr(crRNASequence.length() - 23, 20);
        DNASequence = DNASequence.substr(DNASequence.length() - 23, 20);
    } else {
        std::cerr << "Error: crRNA and DNA sequences must be at least 23 nucleotides long." << std::endl;
        return 1;
    }

    // Replace "T" with "U" in crRNA and transform both DNA and crRNA to uppercase
    for (char &ch : crRNASequence) {
        if (ch == 'T') {
            ch = 'U';
        }
        ch = std::toupper(ch);
    }

    for (char &ch : DNASequence) {
        ch = std::toupper(ch);
    }

    double product = 1.0;
    product *= cfdScores[PAM];

    if (product > 0) {
        // Perform pairwise comparison and product calculation
        for (int i = 0; i < crRNASequence.length(); ++i) {
            // Construct label for each pair with underscore and position
            std::string label = "r" + std::string(1, crRNASequence[i]) + ":d" + revComp(std::string(1, DNASequence[i])) + "_" + std::to_string(i + 1);

            // Retrieve CFD score from the map
            double cfdScore = cfdScores[label];

            // Multiply the product
            product *= cfdScore;
        }
    }
    
    // Print the final CFD score
    std::cout << product << std::endl;

    return 0;
}
