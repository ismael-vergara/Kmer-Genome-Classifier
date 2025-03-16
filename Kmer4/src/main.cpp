/*
 * Metodología de la Programación: Kmer4
 * Curso 2023/2024
 */

/**
 * @file main.cpp
 * @author Silvia Acid Carrillo <acid@decsai.ugr.es>
 * @author Andrés Cano Utrera <acu@decsai.ugr.es>
 * @author Luis Castillo Vidal <L.Castillo@decsai.ugr.es>
 * @author Javier Martínez Baena <jbaena@ugr.es>
 *
 * Created on 17 November 2023, 12:45
 */

#include <iostream>

#include "Profile.h"

/**
 * Shows help about the use of this program in the given output stream
 * @param outputStream The output stream where the help will be shown (for example,
 * cout, cerr, etc) 
 */
void showEnglishHelp(std::ostream& outputStream) {
    outputStream << "ERROR in Kmer4 parameters" << std::endl;
    outputStream << "Run with the following parameters:" << std::endl;
    outputStream << "kmer4 [-t min|max] <file1.prf> <file2.prf> [ ... <filen.prf>]" << std::endl;
    outputStream << std::endl;
    outputStream << "Parameters:" << std::endl;
    outputStream << "-t min | -t max: search for minimun distances or maximum distances (-t min by default)" << std::endl;
    outputStream << "<file1.prf>: source profile file for computing distances" << std::endl;
    outputStream << "<file2.prf> [ ... <filen.prf>]: target profile files for computing distances" << std::endl;  
    outputStream << std::endl;
    outputStream << "This program computes the distance from profile <file1.prf> to the rest" << std::endl;
    outputStream << std::endl;
}

int PosMinMax (const Profile& profile, const Profile* profiles, int n_input_profiles, bool (*Compare)(double n1, double n2)) {
    
    int pos_min_max = 0;
    double min_max_distance = profile.getDistance(profiles[0]);
    for (int i = 1; i < n_input_profiles; i++) {
        double distance = profile.getDistance(profiles[i]);
        if (Compare(min_max_distance, distance)) {
            min_max_distance = distance;
            pos_min_max = i; 
        }
    }
    
    return pos_min_max;
}

bool Minimum (double n1, double n2) {
    
    return n1 > n2;
}

bool Maximum (double n1, double n2) {
    
    return n1 < n2;
}

/**
 * This program reads an undefined number of Profile objects from the set of 
 * files passed as parameters to main(). All the Profiles object, except the 
 * first one, must be stored in a dynamic array of Profile objects. Then, 
 * for each Profile in the dynamic array, this program prints to the 
 * standard output the name of the file of that Profile and the distance from 
 * the first Profile to the current Profile. 
 * Finally, the program should print in the standard output, the name of 
 * the file with the Profile with the minimum|maximum  distance to the Profile 
 * of the first file and its profile identifier.
 * 
 * At least, two Profile files are required to run this program.
 * 
 * This program assumes that the profile files are already normalized and 
 * sorted by frequency. This is not checked in this program. Unexpected results
 * will be obtained if those conditions are not met.
 * 
 * Running sintax:
 * > kmer4 [-t min|max] <file1.prf> <file2.prf> [  ... <filen.prf>] 
 * 
 * Running example:
 * > kmer4 ../Genomes/human1.prf ../Genomes/worm1.prf ../Genomes/mouse1.prf 
Distance to ../Genomes/worm1.prf: 0.330618
Distance to ../Genomes/mouse1.prf: 0.224901
Nearest profile file: ../Genomes/mouse1.prf
Identifier of the nearest profile: mus musculus
 * 
 * Running example:
 * > kmer4 -t max ../Genomes/human1.prf ../Genomes/worm1.prf ../Genomes/mouse1.prf 
Distance to ../Genomes/worm1.prf: 0.330618
Distance to ../Genomes/mouse1.prf: 0.224901
Farthest profile file: ../Genomes/worm1.prf
Identifier of the farthest profile: worm
 */
int main(int argc, char* argv[]) {
    
    // Process the main() arguments
    bool min = true;
    int pos_first_file;
    if (argc < 3) {
        showEnglishHelp(std::cerr);
        return 1;
    }
    
    if (argv[1][0] == '-') {
        if (std::string(argv[1]) == "-t") {
            if (std::string(argv[2]) == "max") 
                min = false;
            
            else if (std::string(argv[2]) != "min") {
                showEnglishHelp(std::cerr);
                return 1;
            }

            pos_first_file = 3;
        } 
        
        else {
            showEnglishHelp(std::cerr);
            return 1;            
        }
    }   
        
    else pos_first_file = 1; 
 
    // Allocate a dynamic array of Profiles
    int n_input_profiles = argc - pos_first_file - 1;
    Profile* profiles;
    profiles = new Profile[n_input_profiles];
    
    // Load the input Profiles
    Profile profile;
    profile.load(argv[pos_first_file]);
    
    for (int i = 0; i < n_input_profiles; i++)
        profiles[i].load(argv[pos_first_file + 1 + i]);   
            
    // Calculate and print the distance from the first Profile to the rest
    for (int i = 0; i < n_input_profiles; i++) 
        std::cout << "Distance to " << argv[pos_first_file + 1 + i] << ": " << profile.getDistance(profiles[i]) << std::endl;

    // Print name of the file and identifier that takes min|max distance to the first one
    int pos_min_max = PosMinMax (profile, profiles, n_input_profiles, (min ? Minimum : Maximum));
    std::cout << (min ? "Nearest" : "Farthest") << " profile file: " << argv[pos_min_max + pos_first_file + 1] << std::endl
              << "Identifier of the "<< (min ? "nearest" : "farthest") << " profile: " << profiles[pos_min_max].getProfileId() << std::endl;
    
    // Deallocate the dynamic array of Profile
    delete[] profiles;
    
    return 0; 
}
