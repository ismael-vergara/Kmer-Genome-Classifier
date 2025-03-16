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

Kmer::Kmer(int k) {

    if (k < 1)
        throw std::invalid_argument("Kmer::Kmer(int k): k is less than 1");

    _text = std::string(k, MISSING_NUCLEOTIDE);
}

Kmer::Kmer(const std::string& text) {

    if (text.empty())
        throw std::invalid_argument("Kmer(const std::string& text): text is an empty string");

    _text = text;
}

int Kmer::getK() const {

    return _text.size();
}

int Kmer::size() const {

    return getK();
}

std::string Kmer::toString() const {

    return _text;
}

const char& Kmer::at(int index) const {

    if (index < 0 || index >= size())
        throw std::out_of_range("const char& Kmer::at(int index) const: invalid position " + std::to_string(index));

    return _text.at(index);
}

char& Kmer::at(int index) {

    if (index < 0 || index >= size())
        throw std::out_of_range("char& Kmer::at(int index): invalid position " + std::to_string(index));

    return _text.at(index);
}

void Kmer::toLower() {

    int size = getK();
    for (int i = 0; i < size; i++)
        at(i) = tolower(at(i));
}

void Kmer::toUpper() {

    int size = getK();
    for (int i = 0; i < size; i++)
        at(i) = toupper(at(i));
}

void Kmer::normalize(const std::string& validNucleotides) {

    toUpper();
    int size = getK();
    for (int i = 0; i < size; i++)
        if (!IsValidNucleotide(at(i), validNucleotides))
            at(i) = MISSING_NUCLEOTIDE;
}

Kmer Kmer::complementary(const std::string& nucleotides, const std::string& complementaryNucleotides) const { // for?

    if (nucleotides.size() != complementaryNucleotides.size())
        throw std::invalid_argument(std::string("Kmer Kmer::complementary(const std::string& nucleotides, ") +
            "const std::string& complementaryNucleotides) const:" +
            "The sizes of nucleotides and complementaryNucleotides are not the same");

    std::string complementary_kmer_string;
    int size = getK();
    for (int i = 0; i < size; i++) {
        if (IsValidNucleotide(at(i), nucleotides))
            complementary_kmer_string.push_back(complementaryNucleotides.at(nucleotides.find(at(i))));

        else
            complementary_kmer_string.push_back(at(i));
    }

    return Kmer(complementary_kmer_string);
}

void Kmer::write(std::ostream& outputStream) const {

    outputStream.write(_text.c_str(), sizeof(char) * getK() + 1);
}

void Kmer::read(std::istream& inputStream) {

    std::string kmer;
    char nucleotide;

    while (inputStream.get(nucleotide) && nucleotide != '\0')
        kmer += nucleotide;

    _text = kmer;
}

const char& Kmer::operator[](int index) const {
    
    return _text.at(index);
}

char& Kmer::operator[](int index) {
    
    return _text.at(index);
}
bool IsValidNucleotide(char nucleotide, const std::string& validNucleotides) {

    return (validNucleotides.find(nucleotide) != std::string::npos);
}

void ToLower(Kmer& kmer) {

    int size = kmer.getK();
    for (int i = 0; i < size; i++)
        kmer.at(i) = tolower(kmer.at(i));
}

void ToUpper(Kmer& kmer) {

    int size = kmer.getK();
    for (int i = 0; i < size; i++)
        kmer.at(i) = toupper(kmer.at(i));
}

std::ostream& operator<<(std::ostream& os, const Kmer& kmer) {

    os << kmer.toString();

    return os;
}

std::istream& operator>>(std::istream& is, Kmer& kmer) {

    std::string kmer_string;
    is >> kmer_string;
    kmer = Kmer(kmer_string);

    return is;
}

bool operator>(const Kmer& kmer1, const Kmer& kmer2) {

    return kmer1.toString() < kmer2.toString();
}

bool operator==(const Kmer& kmer1, const Kmer& kmer2) {

    return kmer1.toString() == kmer2.toString();
}