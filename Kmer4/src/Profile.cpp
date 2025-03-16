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

Profile::Profile(): _profileId("unknown") {

    allocate(INITIAL_CAPACITY);
}

Profile::Profile(int size): _profileId("unknown") {
    
    if (size < 0) 
        throw std::out_of_range("Profile::Profile(int size): size is less than 0");
    
    allocate(size);
    _size = size;
}

Profile::Profile(const Profile& orig) {
    
    copy(orig);
}

Profile::~Profile() {
    
    deallocate();
}

Profile& Profile::operator=(const Profile& orig) {
  
    if (&orig != this)
        copy(orig);
    
    return *this;
}


const std::string& Profile::getProfileId() const {
    
    return _profileId;
}

void Profile::setProfileId(const std::string& id) {
    
    _profileId = id;
}

const KmerFreq& Profile::at(int index) const {
    
    if (index < 0) 
        throw std::out_of_range("const KmerFreq& Profile::at (int index) const: index is less than 0");
    
    else if (index >= _size) 
        throw std::out_of_range("const KmerFreq& Profile::at (int index) const: index is  greater than _size");
    
    return _vectorKmerFreq[index];
}

KmerFreq& Profile::at(int index) {
    
    if (index < 0) 
        throw std::out_of_range("KmerFreq& Profile::at (int index): index is less than 0");
    
    else if (index >= _size) 
        throw std::out_of_range("KmerFreq& Profile::at (int index): index is  greater than _size");
    
    return _vectorKmerFreq[index];
}

int Profile::getSize() const {
    
    return _size;
}

int Profile::getCapacity() const {
    
    return _capacity;
}

double Profile::getDistance(const Profile& otherProfile) const {
    
    int size_1 = getSize(), size_2 = otherProfile.getSize();
    if (size_1 == 0 || size_2 == 0) 
        throw std::invalid_argument("double Profile::getDistance(const Profile& otherProfile) const: the implicit object or the argument Profile object are empty, that is, they do not have any kmer");
    
    double sum = 0.0;
    for (int i = 0; i < size_1; i++) {
        size_t pos_found = otherProfile.findKmer(_vectorKmerFreq[i].getKmer());
        if (pos_found == std::string::npos)
            sum += abs(i - size_2);
 
        else sum += abs(i - pos_found);
    }
    
    return (sum / (size_1 * size_2));
}

int Profile::findKmer(const Kmer& kmer, int initialPos, int finalPos) const {
    
    bool found = false;
    size_t foundPos = std::string::npos;
    if (finalPos >= _size) 
        finalPos = _size - 1;
    
    for (int i = initialPos; i <= finalPos && !found; i++) {  
        if (_vectorKmerFreq[i].getKmer().toString() == kmer.toString()) {
            found = true;
            foundPos = i;   
        }
    }
    
    return foundPos;  
}

int Profile::findKmer(const Kmer& kmer) const {
    
    return findKmer(kmer, 0, _size - 1);
}

std::string Profile::toString() const {
    
    std::string string = getProfileId() + '\n' + std::to_string(_size) + '\n';
    for (int i = 0; i < _size; i++) 
        string += _vectorKmerFreq[i].toString() + '\n';
    
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
    if (!output) {
        output.close();
        throw std::ios_base::failure("void Profile::save(const char fileName[]) const: the given file cannot be opened");
    }
    
    output << MAGIC_STRING_T << std::endl;
    output << toString();
    if (!output) {
        output.close();
        throw std::ios_base::failure("void Profile::save(const char fileName[]) const: an error ocurred while writing to the file");   
    }
    
    output.close();
}


void Profile::load(const char fileName[]) {
    
    deallocate();
    std::ifstream input;
    input.open(fileName);
    if (!input) {
        input.close();
        throw std::ios_base::failure("void Profile::load(const char fileName[]): the given file cannot be opened");
    }
    
    std::string magic_string;
    getline(input,magic_string);
    if (magic_string != MAGIC_STRING_T) {
        input.close();
        throw std::invalid_argument("void Profile::load(const char fileName[]): an invalid magic string is found in the given file");
    }
    
    std::string profile;
    getline(input, profile);
    setProfileId(profile);     
    int nkmer_freqs;
    input >> nkmer_freqs;    
    if (nkmer_freqs < 0) {    
        input.close();
        throw std::out_of_range("void Profile::load(const char fileName[]): the number of kmers in the given file is negative");
    }
                
    for (int i = 0; i < nkmer_freqs; i++) {
        if (input) {
            std::string kmer;
            int frequency;
            input >> kmer >> frequency;
            KmerFreq kmerfreq;
            kmerfreq.setFrequency(frequency);
            kmerfreq.setKmer(kmer);
            append(kmerfreq); 
        } 
                    
        else {
            input.close();
            throw std::ios_base::failure("void Profile::load(const char fileName[]): an error ocurred while reading from the file");
        }                  
    }  

    input.close();
}

void Profile::append(const KmerFreq& kmerFreq) {
    
    size_t foundPos = findKmer(kmerFreq.getKmer());
    if (foundPos == std::string::npos) {
        if (_size == _capacity) 
            reallocate(_size + 1);
        
        _vectorKmerFreq[_size] = kmerFreq;
        _size++;
    }
    
    else _vectorKmerFreq[foundPos].setFrequency(_vectorKmerFreq[foundPos].getFrequency() + kmerFreq.getFrequency());      
}

void Profile::normalize(const std::string& validNucleotides) {
       
    for (int i = 0; i < _size; i++) {
        Kmer kmer = _vectorKmerFreq[i].getKmer();
        kmer.normalize(validNucleotides);
        _vectorKmerFreq[i].setKmer(kmer);
    }
    
    int i = 1;
    while (i < _size) {
       size_t index = findKmer(_vectorKmerFreq[i].getKmer(), 0, i - 1);
       if (index != std::string::npos) {
           _vectorKmerFreq[index].setFrequency(_vectorKmerFreq[index].getFrequency() + _vectorKmerFreq[i].getFrequency());
           deletePos(i);
       }
       
       else i++;
   } 
}

void Profile::deletePos(int pos) {
    
    if (pos < 0) 
        throw std::out_of_range("void Profile::deletePos(int pos): pos is less than 0");

    else if (pos >= _size) 
        throw std::out_of_range("void Profile::deletePos(int pos): pos is greater than _size");
    
    for (int i = pos; i < _size - 1; i++) 
        _vectorKmerFreq[i] = _vectorKmerFreq[i + 1];  

    _size--;
}

void Profile::zip(bool deleteMissing, int lowerBound) {
  
    int i = 0;
    while (i < _size) {
        if ((deleteMissing && (_vectorKmerFreq[i].getKmer().toString().find(Kmer().MISSING_NUCLEOTIDE) != std::string::npos))||
            (_vectorKmerFreq[i].getFrequency() <= lowerBound)) 
            deletePos(i);
        
        else i++;
    }   
} 

void Profile::join(const Profile& profile) {
      
    int size = profile.getSize();
    for (int i = 0; i < size; i++) 
        append(profile.at(i));
}

void Profile::allocate(int capacity) {
    
    if (capacity < 0) 
        throw std::out_of_range("void Profile::allocate(int capacity): capacity is less than 0");
    
    _vectorKmerFreq = new KmerFreq[capacity];
    _size = 0;
    _capacity = capacity;
}

void Profile::deallocate() {
    
    delete[] _vectorKmerFreq;
    _vectorKmerFreq = nullptr;
    _size = 0;
    _capacity = 0;
}

void Profile::reallocate(int capacity) {
    
    if (capacity < 0)
        throw std::out_of_range("void Profile::reallocate(int capacity): capacity is less than 0");
    
    else if (capacity < _capacity)
        throw std::out_of_range("void Profile::reallocate(int capacity): capacity is less than current _capacity");
    
    while (_capacity < capacity)
        _capacity += BLOCK_SIZE;
            
    KmerFreq* kmerfreq = new KmerFreq[_capacity];
    for (int i = 0; i < _size; i++)
        kmerfreq[i] = _vectorKmerFreq[i];
    
    delete[] _vectorKmerFreq;
    _vectorKmerFreq = kmerfreq;
}

void Profile::copy(const Profile& profile) {
    
    setProfileId(profile.getProfileId());
    reallocate(profile.getCapacity());
    for(int i = 0; i < profile.getSize(); i++)
        _vectorKmerFreq[i] = profile.at(i);

    _size = profile.getSize();
}