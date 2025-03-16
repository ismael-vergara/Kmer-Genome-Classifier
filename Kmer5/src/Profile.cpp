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
#include <limits>

#include "Profile.h"

const std::string Profile::MAGIC_STRING_T="MP-KMER-T-1.0";
const std::string Profile::MAGIC_STRING_B="MP-KMER-B-1.0";

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

    allocate(orig.getCapacity());
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
    
    else if (index >= getSize()) 
        throw std::out_of_range("const KmerFreq& Profile::at (int index) const: index is  greater than _size");
    
    return _vectorKmerFreq[index];
}

KmerFreq& Profile::at(int index) {
    
    if (index < 0) 
        throw std::out_of_range("KmerFreq& Profile::at (int index): index is less than 0");
    
    else if (index >= getSize()) 
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
        size_t pos_found = otherProfile.findKmer(at(i).getKmer());
        if (pos_found == std::string::npos)
            sum += abs(i - size_2);
 
        else sum += abs(i - pos_found);
    }
    
    return (sum / (size_1 * size_2));
}

int Profile::findKmer(const Kmer& kmer, int initialPos, int finalPos) const {
    
    bool found = false;
    size_t foundPos = std::string::npos;
    if (finalPos >= getSize()) 
        finalPos = getSize() - 1;
    
    for (int i = initialPos; i <= finalPos && !found; i++) {  
        if (at(i).getKmer() == kmer) {
            found = true;
            foundPos = i;   
        }
    }
    
    return foundPos;  
}

int Profile::findKmer(const Kmer& kmer) const {
    
    return findKmer(kmer, 0, getSize() - 1);
}

std::string Profile::toString() const {
    
    int size = getSize();
    std::string string = getProfileId() + '\n' + std::to_string(size) + '\n';
    for (int i = 0; i < size; i++) 
        string += at(i).toString() + '\n';
    
    return string; 
}

void Profile::sort() {
    
    KmerFreq insert;
    int size = getSize();
    for (int left = 1; left < size; left++) {
        insert = at(left);
        int i = left;
        while(i > 0 && insert > at(i - 1)) {
            at(i) = at(i - 1);
            i--;  
        }
        
        at(i) = insert; 
    }
}

void Profile::save(const char fileName[], char mode) const {
        
    if (mode != 't' && mode != 'b')
        throw std::invalid_argument("void Profile::save(const char fileName[], char mode) const: the given mode is not valid ('t' or 'b')");

    else {
        std::ofstream output(fileName);
        if (!output) {
            output.close();
            throw std::ios_base::failure("void Profile::save(const char fileName[]) const: the given file cannot be opened");
        }
        
        if (mode == 't')
            output << MAGIC_STRING_T << std::endl << *this;
        
        else {
            int size = getSize();
            output << MAGIC_STRING_B << std::endl << getProfileId() << std::endl << size << std::endl;
            for (int i = 0; i < size; i++) {
                at(i).write(output);
            }    
        }
        
        if (!output) {
            output.close();
            throw std::ios_base::failure("void Profile::save(const char fileName[]) const: an error ocurred while writing to the file");
        }
    }
}


void Profile::load(const char fileName[]) {
    
    deallocate();
    std::ifstream input(fileName);
    if (!input) {
        input.close();
        throw std::ios_base::failure("void Profile::load(const char fileName[]): the given file cannot be opened");
    }

    std::string magic_string;
    getline(input,magic_string);
    if (magic_string == MAGIC_STRING_T)
        input >> *this;
    
    else if (magic_string == MAGIC_STRING_B) {
        std::string profile;
        getline(input, profile);
        setProfileId(profile);
        int nkmer_freqs;
        input >> nkmer_freqs;
        input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        if (nkmer_freqs < 0) {    
            input.close();
            throw std::out_of_range("void Profile::load(const char fileName[]): the number of kmers in the given file is negative");
        }
        
        for (int i = 0; i < nkmer_freqs; i++) {
            if (input) {
                KmerFreq kmerFreq;
                kmerFreq.read(input);
                this->append(kmerFreq);
            }
            
            else {
                input.close();
                throw std::ios_base::failure("void Profile::load(const char fileName[]): an error ocurred while reading from the file");
            }
        }
    }
    
    else {
        input.close();
        throw std::invalid_argument("void Profile::load(const char fileName[]): an invalid magic string is found in the given file");
    }

    input.close();
}

void Profile::append(const KmerFreq& kmerFreq) {
    
    size_t foundPos = findKmer(kmerFreq.getKmer());
    if (foundPos == std::string::npos) {
        if (getSize() == getCapacity()) 
            reallocate(getSize() + 1);
        
        this->operator[](getSize()) = kmerFreq;
        _size++;
    }
    
    else 
        at(foundPos).setFrequency(at(foundPos).getFrequency() + kmerFreq.getFrequency());      
}

void Profile::normalize(const std::string& validNucleotides) {
    
    int size = getSize();   
    for (int i = 0; i < size; i++) {
        Kmer kmer = at(i).getKmer();
        kmer.normalize(validNucleotides);
        at(i).setKmer(kmer);
    }
      
    int i = 1;
    while (i < getSize()) {
        size_t index = findKmer(at(i).getKmer(), 0, i - 1);
        if (index != std::string::npos) {
            at(index).setFrequency(at(index).getFrequency() + at(i).getFrequency());
            deletePos(i);
        }
       
        else i++;
    } 
}

void Profile::deletePos(int pos) {
    
    int size = getSize();
    if (pos < 0) 
        throw std::out_of_range("void Profile::deletePos(int pos): pos is less than 0");

    else if (pos >= size) 
        throw std::out_of_range("void Profile::deletePos(int pos): pos is greater than _size");
    
    for (int i = pos; i < size - 1; i++) 
        at(i) = at(i + 1);  

    _size--;
}

void Profile::zip(bool deleteMissing, int lowerBound) {
  
    int i = 0;
    while (i < getSize()) {
        if ((deleteMissing && (at(i).getKmer().toString().find(Kmer::MISSING_NUCLEOTIDE) != std::string::npos))||
            (at(i).getFrequency() <= lowerBound)) 
            deletePos(i);
        
        else i++;
    }   
} 

void Profile::join(const Profile& profile) {
      
    int size = profile.getSize();
    for (int i = 0; i < size; i++) 
        append(profile.at(i));
}

const KmerFreq& Profile::operator[](int index) const {
    
    return _vectorKmerFreq[index];
}

KmerFreq& Profile::operator[](int index) {
    
    return _vectorKmerFreq[index];
}

Profile& Profile::operator+=(const KmerFreq& kmerFreq) {
    
    append(kmerFreq); // Como el método append no está como deprecated, lo usamos para no repetir código    
    return *this;
}

Profile& Profile::operator+=(const Profile& profile) {
    
    int size = profile.getSize();
    for (int i = 0; i < size; i++) 
        append(profile.at(i));
    
    return *this;
}

void Profile::allocate(int capacity) {
    
    if (capacity < 0) 
        throw std::out_of_range("void Profile::allocate(int capacity): capacity is less than 0");
    
    if (capacity == 0)
        _vectorKmerFreq = nullptr;
    
    else
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

    int size = getSize();
    if (capacity == 0) {
        _vectorKmerFreq = nullptr;
        _capacity = 0;
        _size = 0;
    }
    
    else if (capacity < 0)
        throw std::out_of_range("void Profile::reallocate(int capacity): capacity is less than 0");
    
    else if (capacity < _size)
        throw std::out_of_range("void Profile::reallocate(int capacity): capacity is less than current _size");
    
    else {
        while (getCapacity() < capacity)
            _capacity += BLOCK_SIZE;

        KmerFreq* kmerfreq = new KmerFreq[getCapacity()];
        for (int i = 0; i < size; i++)
            kmerfreq[i] = at(i);

        delete[] _vectorKmerFreq;
        _vectorKmerFreq = kmerfreq;
    }
}

void Profile::copy(const Profile& profile) {
    
    setProfileId(profile.getProfileId());
    for(int i = 0; i < profile.getSize(); i++)
        _vectorKmerFreq[i] = profile[i];

    _size = profile.getSize();
}

std::ostream& operator<<(std::ostream& os, const Profile& profile) {
    
    int size = profile.getSize();
    os << profile.getProfileId() << std::endl << size;
    for (int i = 0; i < size; i++)
        os << std::endl << profile.at(i);

    return os;
}

std::istream& operator>>(std::istream& is, Profile& profile) {
    
    profile.deallocate();
    std::string profile_string;
    getline(is, profile_string);
    profile.setProfileId(profile_string);     
    int nkmer_freqs;
    is >> nkmer_freqs;
    if (nkmer_freqs < 0)    
        throw std::out_of_range("void Profile::load(const char fileName[]): the number of kmers in the given file is negative");
    
    profile.allocate(nkmer_freqs);         
    for (int i = 0; i < nkmer_freqs; i++) {
        KmerFreq kmerFreq;
        Kmer kmer;
        int frequency;
        is >> kmer >> frequency;
        kmerFreq.setKmer(kmer);
        kmerFreq.setFrequency(frequency);
        profile.append(kmerFreq); 
    }              
    
    return is;
}
