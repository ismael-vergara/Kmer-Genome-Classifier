/*
 * Metodología de la Programación: Kmer5
 * Curso 2023/2024
 */

#include <cstring>

#include "KmerCounter.h"

/** 
 * @file LEARN.cpp
 * @author Silvia Acid Carrillo <acid@decsai.ugr.es>
 * @author Andrés Cano Utrera <acu@decsai.ugr.es>
 * @author Luis Castillo Vidal <L.Castillo@decsai.ugr.es>
 * @author Javier Martínez Baena <jbaena@ugr.es>
 * 
 * Created on 22 December 2023, 10:00
 */

/**
 * Shows help about the use of this program in the given output stream
 * @param outputStream The output stream where the help will be shown (for example,
 * cout, cerr, etc) 
 */
void showEnglishHelp(std::ostream& outputStream) {
    outputStream << "ERROR in LEARN parameters" << std::endl;
    outputStream << "Run with the following parameters:" << std::endl;
    outputStream << "LEARN [-t|-b] [-k kValue] [-n nucleotidesSet] [-p profileId] [-o outputFilename] <file1.dna> [<file2.dna> <file3.dna> .... ]" << std::endl;
    outputStream << std::endl;
    outputStream << "Parameters:" << std::endl;
    outputStream << "-t|-b: text mode or binary mode for the output file (-t by default)" << std::endl;
    outputStream << "-k kValue: number of nucleotides in a kmer (5 by default)" << std::endl;
    outputStream << "-n nucleotidesSet: set of possible nucleotides in a kmer (ACGT by default). "
            << "Note that the characters should be provided in uppercase" << std::endl;
    outputStream << "-p profileId: profile identifier (unknown by default)" << std::endl;
    outputStream << "-o outputFilename: name of the output file (output.prf by default)" << std::endl;
    outputStream << "<file1.dna> <file2.dna> <file3.dna> ....: names of the input files (at least one is mandatory)" << std::endl;
    outputStream << std::endl;
    outputStream << "This program learns a profile model from a set of " <<
            "input DNA files <file1.dna> <file2.dna> <file3.dna> ...." << std::endl;
    outputStream << std::endl;
}

/**
 * This program learns a Profile model from a set of input DNA files (file1.dna,
 * file2.dna, ...). The learned Profile object is then zipped (kmers with any 
 * missing nucleotide or with frequency equals to zero will be removed) 
 * and ordered by frequency and saved in 
 * the file outputFilename (or output.prf if the output file is not provided).
 * 
 * Running sintax:
 * > LEARN [-t|-b] [-k kValue] [-n nucleotidesSet] [-p profileId] [-o outputFilename] <file1.dna> [<file2.dna> <file3.dna> ....]
 * 
 * Running example:
 * > LEARN -k 2 -p bug -o /tmp/unknownACGT.prf ../Genomes/unknownACGT.dna
 * 
 * > cat /tmp/unknownACGT.prf
MP-KMER-T-1.0
bug
7
GG 2
AC 1
AG 1
AT 1
CC 1
GA 1
TA 1
 * 
 * @param argc The number of command line parameters
 * @param argv The vector of command line parameters (cstrings)
 * @return 0 If there is no error; a value > 0 if error
 */
int main(int argc, char *argv[]) {
    // Process the main() arguments
    if (argc < 2) {
        showEnglishHelp(std::cerr);
        return 1;
    }

    int pos = 1;
    bool continues = true;
    char mode = 't';
    int kValue = 5;
    std::string nucleotidesSet = "ACGT";
    std::string profileId;
    std::string outputFilename = "output.prf";
    while (pos < argc && continues) {
        if (argv[pos][0] == '-') {
            if (strlen(argv[pos]) == 2) {
                switch(argv[pos][1]) {
                    case 't': 
                        pos++;
                        break;
                    
                    case 'b':
                        mode = 'b';
                        pos++;
                        break;
                        
                    case 'k':
                        if (pos + 1 < argc) {
                            kValue = atoi(argv[pos + 1]);
                            pos += 2;
                        }
                        
                        else {
                            showEnglishHelp(std::cerr);
                            return 1;
                        }
                        
                        break;

                    case 'n':
                        if (pos + 1 < argc) {
                            nucleotidesSet = argv[pos + 1];
                            pos += 2;
                        }
                        
                        else {
                            showEnglishHelp(std::cerr);
                            return 1;
                        }
                        
                        break;

                    case 'p':
                        if (pos + 1 < argc) {
                            profileId = argv[pos + 1];
                            pos += 2;
                        }
                        
                        else {
                            showEnglishHelp(std::cerr);
                            return 1;
                        }
                        
                        break;

                    case 'o': 
                        if (pos + 1 < argc) {
                            outputFilename = argv[pos + 1];
                            pos += 2;
                        }
                        
                        else {
                            showEnglishHelp(std::cerr);
                            return 1;                            
                        }
                        
                        break;

                    default:
                        showEnglishHelp(std::cerr);
                        return 1;
                        break;
                }
            }    
            
            else {
                showEnglishHelp(std::cerr);
                return 1;
            }
        }

        else
            continues = false;
    }

    if (pos == argc) {
        showEnglishHelp(std::cerr);
        return 1;
    }

    // Loop to calculate the kmer frecuencies of the input genome files using a KmerCounter object
    KmerCounter kmerCounter(kValue, nucleotidesSet);
    
    while (pos < argc) {
        KmerCounter tmp(kValue, nucleotidesSet);
        tmp.calculateFrequencies(argv[pos]);
        kmerCounter += tmp;
        pos++;
    }

    // Obtain a Profile object from the KmerCounter object
    Profile profile = kmerCounter.toProfile();
    if (!profileId.empty())
        profile.setProfileId(profileId);

    // Zip the Profile object
    profile.zip(true);

    // Sort the Profile object
    profile.sort();

    // Save the Profile object in the output file
    profile.save(outputFilename.c_str(), mode);

    return 0;
}