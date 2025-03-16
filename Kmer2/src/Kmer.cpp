/*
 * Metodología de la Programación: Kmer0
 * Curso 2023/2024
 */

/** 
 * @file Kmer.cpp
 * @author Silvia Acid Carrillo <acid@decsai.ugr.es>
 * @author Andrés Cano Utrera <acu@decsai.ugr.es>
 * @author Luis Castillo Vidal <L.Castillo@decsai.ugr.es>
 * @author Javier Martínez Baena <jbaena@ugr.es>
 * 
 * Created on 24 October 2023, 14:00
 */

#include "Kmer.h"

using namespace std;

Kmer::Kmer(int k) {
    
    if (k < 1) {
            
        throw std::invalid_argument("Kmer::Kmer(int k): k is less than 1");
    }
    
    else {
        
        _text = string(k, MISSING_NUCLEOTIDE);
    }
}

Kmer::Kmer(const std::string& text) {
        
    if (text.empty()) {
            
        throw std::invalid_argument(std::string("Kmer(const std::string& text): ") +
                                                "text is an empty string");
    }
    
    else {
        
        _text = text;
    }    
}

int Kmer::getK() const {
    
    return _text.size();
}

int Kmer::size() const {
    
    return _text.size();
}

std::string Kmer::toString() const {
    
    return _text;
}

const char& Kmer::at(int index) const {
    
    if ((index < 0)||(index >= _text.size())) {
            
        throw std::out_of_range(std::string("const char& Kmer::at(int index) const: ") +
                                            "invalid position " + std::to_string(index));
    }
    
    else {
        
        return _text.at(index);
    }      
}
    
char& Kmer::at(int index) {
    
    if ((index < 0)||(index >= _text.size())) {
            
        throw std::out_of_range("char& Kmer::at(int index): invalid position " +
                                std::to_string(index));
    }
    
    else {
        
        return _text.at(index);
    }
}

void Kmer::toLower() {
    
    int length = _text.size();        
            
    for (int i = 0; i < length; i++) {
        
        at(i) = tolower(at(i));
    }
}

void Kmer::toUpper() {
    
    const int length = size();
    
    for (int i = 0; i < length; i++) {
        
        _text.at(i) = toupper(_text.at(i));
    }
}

void Kmer::normalize(const std::string& validNucleotides) {
    
    toUpper();
    
    int length = _text.size();
  
    for (int i = 0; i < length; i++) {
    
        if (!IsValidNucleotide(_text.at(i), validNucleotides)) {
            
            _text.at(i) = MISSING_NUCLEOTIDE;
        }
    }   
}

Kmer Kmer::complementary(const std::string& nucleotides, 
                         const std::string& complementaryNucleotides) const {
    
    if (nucleotides.size() != complementaryNucleotides.size()) {

        throw std::invalid_argument(std::string("Kmer Kmer::complementary(const std::string& nucleotides, ") +
                                    "const std::string& complementaryNucleotides) const:" +
                                    "The sizes of nucleotides and complementaryNucleotides are not the same");     
    }
    
    else {
        
        std::string complementary_kmer_string;
    
        int length = _text.size();
    
        for (int i = 0; i < length; i++) {
        
            if (IsValidNucleotide(_text.at(i), nucleotides))  {
            
                complementary_kmer_string.push_back(complementaryNucleotides.
                                          at(nucleotides.find(_text.at(i))));
            } 
        
            else {
            
                complementary_kmer_string.push_back(_text.at(i));
            }
        }   
    
        return Kmer(complementary_kmer_string);
    }
}

bool IsValidNucleotide (char nucleotide, const std::string& validNucleotides) {
    
    return (validNucleotides.find(nucleotide) != std::string::npos);
}

void ToLower(Kmer& kmer) {
    
    int length = kmer.size();
    
    for (int i = 0; i < length; i++) {
     
        kmer.at(i) = tolower(kmer.at(i));
    }
}

void ToUpper(Kmer& kmer) {
    
    int length = kmer.size();
    
    for (int i = 0; i < length; i++) {
        
        kmer.at(i) = toupper(kmer.at(i));
    }
}