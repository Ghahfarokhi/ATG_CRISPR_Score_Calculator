#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>

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

std::map<std::string, double> cfdScoresFunc(std::string cfdScoresFilePath){
    // Read the cfdScoresFile that will be used in both singlePair and pairList modes.
    std::ifstream cfdScoresFile(cfdScoresFilePath);
    if (!cfdScoresFile.is_open()) {
        std::cerr << "Error opening file: " << cfdScoresFilePath << std::endl;
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

    return cfdScores;
}

double singlePairFunc (std::map<std::string, double> cfdScores, std::string crRNASequence, std::string DNASequence) {
    // Replace "T" with "U" in crRNA and transform both DNA and crRNA to uppercase
    for (char &ch : crRNASequence) {
        ch = std::toupper(ch);
        if (ch == 'T') {
            ch = 'U';
        }
    }

    for (char &ch : DNASequence) {
        ch = std::toupper(ch);
    }

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

    double cfdProduct = 1.0;
    cfdProduct *= cfdScores[PAM];

    if (cfdProduct > 0) {
        // Perform pairwise comparison and product calculation
        for (int i = 0; i < crRNASequence.length(); ++i) {
            // Construct label for each pair with underscore and position
            std::string label = "r" + std::string(1, crRNASequence[i]) + ":d" + revComp(std::string(1, DNASequence[i])) + "_" + std::to_string(i + 1);

            // Retrieve CFD score from the map
            double cfdScore = cfdScores[label];

            // Multiply the product
            cfdProduct *= cfdScore;
        }
    }

    return cfdProduct;
};

void printHelp(){
    const char *helpText = "\nCutting Frequency Determination (CFD) Calculator\n\
    Usage: \n\
        cfdScoresCalculator -h \n\
        OR\n\
        cfdScoresCalculator singlePair <cfdScoresFilePath> <crRNASequence> <DNASequence> \n\
        OR\n\
        cfdScoresCalculator pairList <cfdScoresFilePath> <pairListFilePath>\n\n\
    Examples:\n\
        ./cfdScoresCalculator singlrPair data/CFD_Scores.txt ATCGATGCTGATGCTAGATAAGG ACCGATGCTGATGCTAGATAAGG \n\
        ./cfdScoresCalculator pairList data/CFD_Scores.txt data/crRNA_DNA_Pair_Examples.txt \n\n\
    Note: crRNA and DNA squences provided in singlePair and pairList run modes must be at least 23 nucleotides or longer and end with their PAMs.\n\n\
    Reference: Doench, J., et al. (2016) Nat Biotechnol 34, 184–191. \n\n\
    Developer: Amir.Taheri.Ghahfarokhi@Gmail.com \n\
    Release: 2023-12-26\n\n";
    std::cout << helpText << std::endl;
}

int main(int argc, char *argv[]) {
    // Define variables that will be used to store command line arguments
    std::string runMode, cfdScoresFilePath, crRNASequence, DNASequence, pairListFilePath;
    double cfdProduct = 1.0;

    // Check the provided arguments
    if (argc == 2 && std::string(argv[1]) == "-h") {
        printHelp();
        return 0;
    } else if (argc == 5 && std::string(argv[1]) == "singlePair") {
        runMode = std::string(argv[1]);
        cfdScoresFilePath = std::string(argv[2]);
        crRNASequence = std::string(argv[3]);
        DNASequence = std::string(argv[4]);
    } else if (argc == 4 && std::string(argv[1]) == "pairList") {
        runMode = std::string(argv[1]);
        cfdScoresFilePath = std::string(argv[2]);
        pairListFilePath = std::string(argv[3]);
    } else {
        printHelp();
        return 1;
    }

    // Read cfdScores from cfdScoresFilePath
    std::map<std::string, double> cfdScores;
    cfdScores = cfdScoresFunc(cfdScoresFilePath);

    if (runMode == "singlePair") {

        cfdProduct = singlePairFunc(cfdScores, crRNASequence, DNASequence);
        // Print the final CFD score
        std::cout << cfdProduct << std::endl;

        return 0;

    } else if (runMode == "pairList") {
        
        // Read pairList from pairListFilePath
        std::ifstream pairListFile(pairListFilePath);
        if (!pairListFile.is_open()) {
            std::cerr << "Error opening file: " << pairListFilePath << std::endl;
            return 1;
        }

        // Process each line in the pairList file
        std::string pairLine;
        bool isFirstLine = true; 
        while (std::getline(pairListFile, pairLine)) {
            std::istringstream iss(pairLine);
            // Skip the header line
            if (isFirstLine) {
                isFirstLine = false;
                std::cout << pairLine << "\t" << "CFD" << std::endl;
                continue;
            }
                
            // Split the line into crRNA and DNA
            if (std::getline(iss, crRNASequence, '\t') && std::getline(iss, DNASequence, '\t')) {
                cfdProduct = singlePairFunc(cfdScores, crRNASequence, DNASequence);
                std::cout << crRNASequence << "\t" << DNASequence << "\t" << cfdProduct << std::endl;
            } else {
                std::string errorText;
                errorText = "Error parsing pairLine: " + pairLine + " in " + pairListFilePath;
                std::cerr << errorText << std::endl;
                return 1;
            }
        }

        return 0;
    }
}
