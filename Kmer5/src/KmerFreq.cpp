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

KmerFreq::KmerFreq() : _kmer(Kmer()), _frequency(0) {}

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

    if (frequency < 0)
        throw std::out_of_range("void KmerFreq::setFrequency(int frequency): frequency is negative");

    _frequency = frequency;
}

std::string KmerFreq::toString() const {

    return getKmer().toString() + " " + std::to_string(getFrequency());
}

void KmerFreq::write(std::ostream& outputStream) const {

    getKmer().write(outputStream);
    outputStream.write(reinterpret_cast<const char*>(&_frequency), sizeof(_frequency)); 
}

void KmerFreq::read(std::istream& inputStream) {
    
    _kmer.read(inputStream);
    inputStream.read(reinterpret_cast<char*>(&_frequency), sizeof(_frequency));
}

std::ostream& operator<<(std::ostream& os, const KmerFreq& kmerFreq) {

    os << kmerFreq.getKmer() << " " << kmerFreq.getFrequency();

    return os;
}

std::istream& operator>>(std::istream& is, KmerFreq& kmerFreq) {
    
    Kmer kmer;
    int frequency;
    is >> kmer >> frequency;
    kmerFreq.setKmer(kmer);
    kmerFreq.setFrequency(frequency);

    return is;
}

bool operator>(const KmerFreq& kmerFreq1, const KmerFreq& kmerFreq2) {

    return kmerFreq1.getFrequency() > kmerFreq2.getFrequency() ||
           (kmerFreq1.getFrequency() == kmerFreq2.getFrequency() &&
            kmerFreq1.getKmer() > kmerFreq2.getKmer());
}

bool operator<(const KmerFreq& kmerFreq1, const KmerFreq& kmerFreq2) {

    return kmerFreq2 > kmerFreq1;
}

bool operator==(const KmerFreq& kmerFreq1, const KmerFreq& kmerFreq2) {

    return kmerFreq1.getKmer() == kmerFreq2.getKmer() &&
           kmerFreq1.getFrequency() == kmerFreq2.getFrequency();
}

bool operator!=(const KmerFreq& kmerFreq1, const KmerFreq& kmerFreq2) {

    return !(kmerFreq1 == kmerFreq2);
}

bool operator<=(const KmerFreq& kmerFreq1, const KmerFreq& kmerFreq2) {

    return kmerFreq2 > kmerFreq1 || kmerFreq1 == kmerFreq2;
}

bool operator>=(const KmerFreq& kmerFreq1, const KmerFreq& kmerFreq2) {

    return kmerFreq1 > kmerFreq2 || kmerFreq1 == kmerFreq2;
}