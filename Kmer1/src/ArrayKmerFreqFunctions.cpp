/*
 * Metodología de la Programación: Kmer1
 * Curso 2023/2024
 */

/** 
 * @file ArrayKmerFreqFunctions.cpp
 * @author Silvia Acid Carrillo <acid@decsai.ugr.es>
 * @author Andrés Cano Utrera <acu@decsai.ugr.es>
 * @author Luis Castillo Vidal <L.Castillo@decsai.ugr.es>
 * @author Javier Martínez Baena <jbaena@ugr.es>
 * 
 * Created on 27 October 2023, 12:00
 */


#include "ArrayKmerFreqFunctions.h"

void ReadArrayKmerFreq(KmerFreq array[], int dim, int& nElements) {
    
    std::string genome;
    int frequency;
    
    std::cin >> nElements;
    
    if (nElements > dim) {
        
        nElements = dim;
        
    }
    
    for (int i = 0; i < nElements; i++) {
        
        std::cin >> genome;
        std::cin >> frequency;
        
        array[i].setKmer(genome);
        array[i].setFrequency(frequency);
    }
}

void PrintArrayKmerFreq(KmerFreq array[], int nElements) {
    
    if (nElements > 0) {
        
        std::cout << nElements << std::endl;
    }
    
    for (int i = 0; i < nElements; i++) {
        
        std::cout << array[i].toString() << std::endl;
    }
}

void SwapElementsArrayKmerFreq(KmerFreq array[], int nElements, int first, int second) {
    
    if (first > nElements || second > nElements) {
        
        throw std::out_of_range(std::string("void SwapElementsArrayKmerFreq(KmerFreq array[], int nElements, int first, int second):") +
                                            " first or second is greater than nElements");
    }
    
    else {
        
        KmerFreq tmp = array[first];
        array[first] = array[second];
        array[second] = tmp;
    }
}

int FindKmerInArrayKmerFreq(KmerFreq array[], Kmer kmer, int initialPos, int finalPos) {
    
    bool found = false;
    int foundPos = std::string::npos;
    
    while (initialPos <= finalPos && !found) {
               
        if (array[initialPos].getKmer().toString() == kmer.toString()) {
            
            found = true;
            foundPos = initialPos;
        }
        
        else initialPos++;
    }
    
    return foundPos;
}

void SortArrayKmerFreq(KmerFreq array[], int nElements) {
    
    KmerFreq insert;
    
    for (int left = 0; left < nElements; left++) {
        
        insert = array[left];
        
        int i = left;
        
        while((i > 0) && ((insert.getFrequency() > array[i - 1].getFrequency()) || 
             ((insert.getFrequency() == array[i - 1].getFrequency()) && (insert.toString() < array[i - 1].toString())))) {
            
            array[i] = array[i - 1];
            
            i--;  
        }
        
        array[i] = insert; 
    }
    
}

void NormalizeArrayKmerFreq(KmerFreq array[], int& nElements, 
        const std::string& validNucleotides) {
    
    // Loop to traverse and normalize each one of the kmers in array
          // Normalize kmer i
    
    Kmer kmer;
    
    for (int i = 0; i < nElements; i++) {
        
        kmer = array[i].getKmer();
        kmer.normalize(validNucleotides);
        array[i].setKmer(kmer); 
        
    }
    
    // Loop to traverse the kmers in array from position 1 to position nElements-1
          // index = Position of array[i].getKmer() in the subarray that begins
          //         at position 0 and ends at position i-1
          // If array[i].getKmer() was found in the the subarray from 0 to i-1 
               // Accumulate the frequencies of the kmers at positions 
               //    index and i in the kmer at position index
               // Delete from the array, the kmer at position i 
    
    int index;
    int i = 1;
    
    while (i < nElements) {
       
        index = FindKmerInArrayKmerFreq(array, array[i].getKmer(), 0, i - 1);
        
       if (index != std::string::npos) {
           
           array[index].setFrequency(array[index].getFrequency() + array[i].getFrequency());
           DeletePosArrayKmerFreq(array, nElements, i);
       }
       
       else i++;
   } 
}

void DeletePosArrayKmerFreq(KmerFreq array[], int& nElements, int pos) {
    
    if (pos < 0 || pos > (nElements - 1)) {
       
        throw std::out_of_range(std::string("void DeletePosArrayKmerFreq(KmerFreq array[], int& nElements, int pos):") +
                                            "pos is less than 0 or greater than nElements - 1");
    }
    
    else {
        
        for (int i = pos; i < nElements - 1; i++) {
            
            array[i] = array[i + 1];
        }
        
        nElements--;
    }
}

void ZipArrayKmerFreq(KmerFreq array[], int& nElements, bool deleteMissing, int lowerBound) {
    
    int i = 0;
  
    while (i < nElements) {
        
        if ((deleteMissing && (array[i].getKmer().toString().find(Kmer().MISSING_NUCLEOTIDE) != std::string::npos))||
            (array[i].getFrequency() <= lowerBound)) {

            DeletePosArrayKmerFreq(array, nElements, i);
        }
        
        else i++;
    }   
}    