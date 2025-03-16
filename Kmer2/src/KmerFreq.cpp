/*
 * Metodología de la Programación: Kmer1
 * Curso 2023/2024
 */

/** 
 * @file KmerFreq.cpp
 * @author Silvia Acid Carrillo <acid@decsai.ugr.es>
 * @author Andrés Cano Utrera <acu@decsai.ugr.es>
 * @author Luis Castillo Vidal <L.Castillo@decsai.ugr.es>
 * @author Javier Martínez Baena <jbaena@ugr.es>
 * 
 * Created on 27 de octubre de 2023, 11:03
 */

#include "KmerFreq.h"

KmerFreq::KmerFreq() : _kmer (Kmer()), _frequency(0) 

{}

const Kmer& KmerFreq::getKmer() const {
    
    return _kmer;
}

int KmerFreq::getFrequency() const {
    
    return _frequency;
}

void KmerFreq::setKmer(const Kmer& kmer) {
    
    _kmer = kmer;
}

void KmerFreq::setFrequency(int frequency) {
    
    if (frequency < 0) {
        
        throw std::out_of_range("frequency is negative");
    }
    
    else {
    
        _frequency = frequency;
    }
}

std::string KmerFreq::toString() const {
    
    return _kmer.toString() + " " + std::to_string(_frequency);
}

