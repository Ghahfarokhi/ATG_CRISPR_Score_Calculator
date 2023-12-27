#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>

std::map<int, double> mitWeights;
std::map<std::string, double> cfdScores;
bool includeDistance = true;

std::string revComp(std::string DNA) {
    std::string complement;
    for (char base : DNA) {
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

std::map<std::string, std::string> prepareRnaDnaPam (std::string RNA, std::string DNA){
    std::map<std::string, std::string> RnaDnaPam;
    // Replace "U" with "T" in crRNA and transform both DNA and crRNA to uppercase
    for (char &ch : RNA) {
        ch = std::toupper(ch);
        if (ch == 'U') {
            ch = 'T';
        }
    }

    for (char &ch : DNA) {
        ch = std::toupper(ch);
    }

    // Extract the rightmost 2 characters of the DNA as PAM
    std::string PAM = DNA.substr(DNA.length() - 2);

    // Update crRNA and DNA with the exact 20 characters upstream of the last 3 characters
    if (RNA.length() >= 23 && DNA.length() >= 23) {
        RNA = RNA.substr(RNA.length() - 23, 20);
        DNA = DNA.substr(DNA.length() - 23, 20);
    } else {
        std::cerr << "Error: crRNA and DNA sequences must be at least 23 nucleotides long." << std::endl;
    }
    RnaDnaPam["RNA"]= RNA;
    RnaDnaPam["DNA"]= DNA;
    RnaDnaPam["PAM"]= PAM;

    return RnaDnaPam;
}

std::map<int, double> importMITWeights (){
    std::map<int, double> mitWeights;
    std::ifstream mitWeightsFile("data/MIT_Weights.txt");
    if (!mitWeightsFile.is_open()) {
        std::cerr << "Error opening data/MIT_Weights.txt" << std::endl;
    }
    int pos;
    std::string pam;
    std::string line;
    double w;
    
    while (std::getline(mitWeightsFile, line)) {
        std::istringstream iss(line);
        if (iss >> pos >> w) {
            mitWeights[pos] = w;
        } else {
            std::cerr << "Error reading data/MIT_Weights.txt: "  << line << std::endl;
        }
    }

    return mitWeights;
}

std::map<std::string, double> importCFDScores(){
    // Read the cfdScoresFile that will be used in both singlePair and pairList modes.
    std::ifstream cfdScoresFile("data/CFD_Scores.txt");
    if (!cfdScoresFile.is_open()) {
        std::cerr << "Error opening file: data/CFD_Scores.txt" << std::endl;
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

double calculateMITProduct(std::string RNA, std::string DNA, std::string PAM, bool includeDistance) {
    double mitProduct = 1.0;
    mitProduct *= cfdScores[PAM];
    
    if (mitProduct > 0) {
        std::vector<int> mismatchIndices;
        for (int i = 0; i < RNA.length(); ++i) {
            if (RNA[i] != DNA[i]) {
                mismatchIndices.push_back(i);
            }
        } 

        if (mismatchIndices.empty()) {
            return mitProduct;
        } else {
            int m = mismatchIndices.size();
            double d = (m == 1) ? 19.0 : (static_cast<double>(*std::max_element(mismatchIndices.begin(), mismatchIndices.end()) - *std::min_element(mismatchIndices.begin(), mismatchIndices.end()))) / (m - 1);
            double t1 = 1.0;
            for (int index : mismatchIndices) {
                t1 *= (1.0 - mitWeights[index+1]);
            }
            double t2 = 1.0 / (m * m);
            double t3 = 1.0 / ((19.0 - d) / 19.0 * 4.0 + 1.0);
            if (includeDistance){
                mitProduct = mitProduct * t1 * t2 * t3;
            } else {
                mitProduct = mitProduct * t1 * t2;
            }
        }
    }
    
    return mitProduct;
}

double calculateCFDProduct (std::string RNA, std::string DNA, std::string PAM) {
    double cfdProduct = 1.0;
    cfdProduct *= cfdScores[PAM];
    if (cfdProduct > 0) {
        // Perform pairwise comparison and product calculation
        for (int i = 0; i < RNA.length(); ++i) {
            // Construct label for each pair with underscore and position
            std::string label = "r" + std::string(1, RNA[i]) + ":d" + revComp(std::string(1, DNA[i])) + "_" + std::to_string(i + 1);
            // Retrieve CFD score from the map
            cfdProduct = cfdProduct * cfdScores[label];
        }
    }

    return cfdProduct;
};

void printHelp(){
    const char *helpText = "\nATG_CRISPR_Score_Calculator (CFD & MIT)\n\
Usage: \n\
    ATG_CRISPR_Score_Calculator -h \n\
    OR\n\
    ATG_CRISPR_Score_Calculator singlePair <crRNASequence> <DNASequence> \n\
    OR\n\
    ATG_CRISPR_Score_Calculator pairList <pairListFilePath>\n\n\
Examples:\n\
    ./ATG_CRISPR_Score_Calculator singlrPair ATCGATGCTGATGCTAGATAAGG ACCGATGCTGATGCTAGATAAGG \n\
    ./ATG_CRISPR_Score_Calculator pairList data/crRNA_DNA_Pair_Examples.txt \n\n\
Note: crRNA and DNA squences provided in singlePair and pairList run modes must be at least 23 nucleotides \
or longer and end with their PAMs.\n\n\
\nCFD: Cutting Frequence Determination (CFD) scores range between 0 to 1 and \
represent the percentage activity of a CRISPR guide RNA (crRNA) at a specific off-target sequence. \
CFD scoring system take number, position, and type and mismatches between crRNA::DNA into account. \
The final score is adjustment for non-conical PAMs. Note: CFD scores for RNA/DNA bulges provided in \
the Supplementary Table 19 of the reference publication, however, these scores were not used in the \
Supplementary Python code. In this program, these adjustment factors are also used, however, use this \
feature with cautiously as these have been described vaguely and combinatorial use of mismatches and bulges \
CFD scores is not peer-reviewed.\n\
\nMIT: takes a value between 0 to 1 and represents the likelihood  of a CRISPR guide RNA (crRNA) \
cutting a specific off-target sequence. A higher MIT Score suggests higher off-target editing possibility. \
MIT scores takes number, position, and distance of mismatches between crRNA::DNA into account. \
Adjustment for non-conical PAMs is done using CFD Scores for non-conical PAMs.\n\
\nReferences: \n\
    CFD: Doench, J., et al. (2016) Nat Biotechnol 34, 184–191. \n\
    MIT: Hsu, P., et al. (2013) Nat Biotechnol 31, 827–832. \n\
\nDeveloper: Amir.Taheri.Ghahfarokhi@Gmail.com \n\
Release: 2023-12-27\n\n";
    std::cout << helpText << std::endl;
}

int main(int argc, char *argv[]) {
    // Define variables that will be used to store command line arguments
    std::string runMode, crRNASequence, DNASequence, pairListFilePath;

    // Check the provided arguments
    if (argc == 2 && std::string(argv[1]) == "-h") {
        printHelp();
        return 0;
    } else if (argc == 4 && std::string(argv[1]) == "singlePair") {
        runMode = std::string(argv[1]);
        crRNASequence = std::string(argv[2]);
        DNASequence = std::string(argv[3]);
    } else if (argc == 3 && std::string(argv[1]) == "pairList") {
        runMode = std::string(argv[1]);
        pairListFilePath = std::string(argv[2]);
    } else {
        std::cerr << "\nError in provided arguments!\n" << std::endl;
        printHelp();
        return 1;
    }

    std::map<std::string, std::string> RnaDnaPam;
    cfdScores = importCFDScores();
    mitWeights = importMITWeights();

    double cfdProduct, mitProduct;
    if (runMode == "singlePair") {
        RnaDnaPam = prepareRnaDnaPam(crRNASequence, DNASequence);
        cfdProduct = calculateCFDProduct(RnaDnaPam["RNA"], RnaDnaPam["DNA"], RnaDnaPam["PAM"]);
        mitProduct = calculateMITProduct(RnaDnaPam["RNA"], RnaDnaPam["DNA"], RnaDnaPam["PAM"], includeDistance);
        std::cout << "CFD" << "\t" << cfdProduct << "\t" << "MIT" << "\t" << mitProduct << std::endl;

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
                std::cout << pairLine << "\t" << "CFD" << "\t" << "MIT" << std::endl;
                continue;
            }
                
            // Split the line into crRNA and DNA
            if (std::getline(iss, crRNASequence, '\t') && std::getline(iss, DNASequence, '\t')) {
                RnaDnaPam = prepareRnaDnaPam(crRNASequence, DNASequence);
                cfdProduct = calculateCFDProduct(RnaDnaPam["RNA"], RnaDnaPam["DNA"], RnaDnaPam["PAM"]);
                mitProduct = calculateMITProduct(RnaDnaPam["RNA"], RnaDnaPam["DNA"], RnaDnaPam["PAM"], includeDistance);
                std::cout << crRNASequence << "\t" << DNASequence << "\t" << cfdProduct << "\t" << mitProduct << std::endl;
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
