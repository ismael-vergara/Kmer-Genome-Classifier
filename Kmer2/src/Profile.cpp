/*
 * Metodología de la Programación: Kmer2
 * Curso 2023/2024
 */

/** 
 * @file Profile.cpp
 * @author Silvia Acid Carrillo <acid@decsai.ugr.es>
 * @author Andrés Cano Utrera <acu@decsai.ugr.es>
 * @author Luis Castillo Vidal <L.Castillo@decsai.ugr.es>
 * @author Javier Martínez Baena <jbaena@ugr.es>
 * 
 * Created on 29 January 2023, 11:00
 */


#include <fstream>

#include "Profile.h"

const std::string Profile::MAGIC_STRING_T="MP-KMER-T-1.0";

Profile::Profile() : _profileId("unknown"), _size(0) {}

Profile::Profile(int size) : _profileId("unknown"), _size(size) {
    
    if (size < 0) {
        
        throw std::out_of_range("Profile::Profile(int size): size is less than 0");
    }
    
    else if (size > DIM_VECTOR_KMER_FREQ) {
        
        throw std::out_of_range("Profile::Profile(int size): size is greater than DIM_VECTOR_KMER_FREQ");
    }
}

const std::string& Profile::getProfileId() const {
    
    return _profileId;
}

void Profile::setProfileId(const std::string& id) {
    
    _profileId = id;
}

const KmerFreq& Profile::at(int index) const {
    
    if (index < 0) {
        
        throw std::out_of_range("const KmerFreq& Profile::at (int index) const: index is less than 0");
    }
    
    else if (index >= _size) {
        
        throw std::out_of_range("const KmerFreq& Profile::at (int index) const: index is  greater than _size");
    }
    
    else {
        
        return _vectorKmerFreq[index];
    }
}

KmerFreq& Profile::at(int index) {
    
    if (index < 0) {
        
        throw std::out_of_range("KmerFreq& Profile::at (int index): index is less than 0");
    }
    
    else if (index >= _size) {
        
        throw std::out_of_range("KmerFreq& Profile::at (int index): index is  greater than _size");
    }
    
    else {
        
        return _vectorKmerFreq[index];
    }
}

int Profile::getSize() const {
    
    return _size;
}

int Profile::getCapacity() const {
    
    return DIM_VECTOR_KMER_FREQ;
}

int Profile::findKmer(const Kmer& kmer, int initialPos, int finalPos) const {
    
    bool found = false;
    size_t foundPos = std::string::npos;
    
    if (finalPos >= _size) {
        
        finalPos = _size - 1;
    }
    
    for (int i = initialPos; i <= finalPos && !found; i++) {  
        
        if (_vectorKmerFreq[i].getKmer().toString() == kmer.toString()) {
            
            found = true;
            foundPos = i;   
        }
    }
    
    return foundPos;  
}

int Profile::findKmer(const Kmer& kmer) const {
    
    bool found = false;
    size_t foundPos = std::string::npos;
    
    for (int i = 0; i < _size && !found; i++) {
        
        if (_vectorKmerFreq[i].getKmer().toString() == kmer.toString()) {
            
            found = true;
            foundPos = i;   
        }
    }
    
    return foundPos;
}

std::string Profile::toString() const {
    
    std::string string = _profileId + '\n' + std::to_string(_size) + '\n';
    
    for (int i = 0; i < _size; i++) {
        
        string += _vectorKmerFreq[i].toString() + '\n';
    }
    
    return string; 
}

void Profile::sort() {
    
    KmerFreq insert;
    
    for (int left = 1; left < _size; left++) {
        
        insert = _vectorKmerFreq[left];
        
        int i = left;
        
        while((i > 0) && ((insert.getFrequency() > _vectorKmerFreq[i - 1].getFrequency()) || 
              ((insert.getFrequency() == _vectorKmerFreq[i - 1].getFrequency()) &&
              (insert.toString() < _vectorKmerFreq[i - 1].toString())))) {
            
            _vectorKmerFreq[i] = _vectorKmerFreq[i - 1];
            
            i--;  
        }
        
        _vectorKmerFreq[i] = insert; 
    }
}

void Profile::save(const char fileName[]) const {
    
    std::ofstream output;
    
    output.open(fileName);
    
    if (output) {
        
        output << MAGIC_STRING_T << std::endl;
        output << toString();
        
        if (!output) {
            
            throw std::ios_base::failure("void Profile::save(const char fileName[]) const: an error ocurred while writing to the file");
        }
    }
    
    else throw std::ios_base::failure("void Profile::save(const char fileName[]) const: the given file cannot be opened");
    
    output.close();
}


void Profile::load(const char fileName[]) {
    
    std::ifstream input;
    
    input.open(fileName);
    
    if (!input) {
        
        throw std::ios_base::failure("void Profile::load(const char fileName[]): the given file cannot be opened");
    }
    
    else {
        
        std::string magic_string;
        getline(input,magic_string);
        
        if (magic_string != MAGIC_STRING_T) {
            
            throw std::invalid_argument("void Profile::load(const char fileName[]): an invalid magic string is found in the given file");
        }
        
        else {
            
            std::string profile;
            getline(input, profile);
            _profileId = profile;
            
            int nkmer_freqs;
            input >> nkmer_freqs;
            
            if (nkmer_freqs > DIM_VECTOR_KMER_FREQ) {
                
                throw std::out_of_range("void Profile::load(const char fileName[]): the number of kmers in the given file exceeds the maximum capacity");
            }
            
            else if (nkmer_freqs < 0) {
                
                throw std::out_of_range("void Profile::load(const char fileName[]): the number of kmers in the given file is negative");
            }
            
            else {
                
                _size = 0;
                
                for (int i = 0; i < nkmer_freqs; i++) {
                    
                    if (!input) {
                        
                        throw std::ios_base::failure("void Profile::load(const char fileName[]): an error ocurred while reading from the file");
                    }
                    
                    else {
                        
                        std::string kmer;
                        int frequency;
                        size_t found;
                        
                        input >> kmer;
                        input >> frequency;
                        
                        found = findKmer(Kmer(kmer));
                            
                        _vectorKmerFreq[i].setFrequency(frequency);
                        _vectorKmerFreq[i].setKmer(Kmer(kmer));
                            
                        _size++;
                    }
                }
            }
        }          
    }  

    input.close();
}

void Profile::append(const KmerFreq& kmerFreq) {
    
    size_t foundPos = findKmer(kmerFreq.getKmer());
    
    if (foundPos == std::string::npos) {
        
        if (_size == DIM_VECTOR_KMER_FREQ) {
            
            throw std::out_of_range("void Profile::append(const KmerFreq& kmerFreq): the array is full");
        }
        
        else {
            
            _vectorKmerFreq[_size].setKmer(kmerFreq.getKmer());
            _vectorKmerFreq[_size].setFrequency(kmerFreq.getFrequency());
            
            _size++;
        }      
    }
    
    else {
        
        _vectorKmerFreq[foundPos].setFrequency(_vectorKmerFreq[foundPos].getFrequency() + kmerFreq.getFrequency());      
    }
}

void Profile::normalize(const std::string& validNucleotides) {
    
    Kmer kmer;
    
    for (int i = 0; i < _size; i++) {
    
        kmer = _vectorKmerFreq[i].getKmer();
        kmer.normalize(validNucleotides);
        _vectorKmerFreq[i].setKmer(kmer);
    }
    
    size_t index;
    int i = 1;
    
    while (i < _size) {
       
       index = findKmer(_vectorKmerFreq[i].getKmer(), 0, i - 1);
        
       if (index != std::string::npos) {
           
           _vectorKmerFreq[index].setFrequency(_vectorKmerFreq[index].getFrequency() + _vectorKmerFreq[i].getFrequency());
           deletePos(i);
       }
       
       else i++;
   } 
}

void Profile::deletePos(int pos) {
    
    if (pos < 0) {
     
        throw std::out_of_range("void Profile::deletePos(int pos): pos is less than 0");
    }
    
    else if (pos >= _size) {
     
        throw std::out_of_range("void Profile::deletePos(int pos): pos is greater than _size");
    }
    
    else {
        
        for (int i = pos; i < _size - 1; i++) {
            
            _vectorKmerFreq[i] = _vectorKmerFreq[i + 1];
        }
    }
    
    _size--;
}

void Profile::zip(bool deleteMissing, int lowerBound) {
    
    int i = 0;
  
    while (i < _size) {
        
        if ((deleteMissing && (_vectorKmerFreq[i].getKmer().toString().find(Kmer().MISSING_NUCLEOTIDE) != std::string::npos))||
            (_vectorKmerFreq[i].getFrequency() <= lowerBound)) {

            deletePos(i);
        }
        
        else i++;
    }   
} 

void Profile::join(const Profile& profile) {
    
    int size = profile.getSize();
    
    if (_profileId == profile.getProfileId()) {
        
        for (int i = 0; i < size; i++) {
        
            append(profile.at(i));
        }
    }
}