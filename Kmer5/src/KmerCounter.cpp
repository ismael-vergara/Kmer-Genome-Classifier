/*
 * Metodología de la Programación: Kmer5
 * Curso 2023/2024
 */

/** 
 * @file KmerCounter.cpp
 * @author Silvia Acid Carrillo <acid@decsai.ugr.es>
 * @author Andrés Cano Utrera <acu@decsai.ugr.es>
 * @author Luis Castillo Vidal <L.Castillo@decsai.ugr.es>
 * @author Javier Martínez Baena <jbaena@ugr.es>
 * 
 * Created on 22 December 2023, 10:00
 */

#include <complex>
#include <fstream>

#include "KmerCounter.h"

/**
 * DEFAULT_VALID_NUCLEOTIDES is a c-string that contains the set of characters
 * that will be considered as valid nucleotides. 

 * The constructor of the class KmerCounter uses this c-string as a 
 * default parameter. It is possible to use a different c-string if that
 * constructor is used with a different c-string
 */
const char* const KmerCounter::DEFAULT_VALID_NUCLEOTIDES = "ACGT";

KmerCounter::KmerCounter(int k, const std::string& validNucleotides) : _k(k), _validNucleotides(validNucleotides),
_allNucleotides(Kmer::MISSING_NUCLEOTIDE + validNucleotides) {

    allocate(getNumRows(), getNumCols());
    initFrequencies();
}

KmerCounter::KmerCounter(const KmerCounter& orig) {

    allocate(orig.getNumRows(), orig.getNumCols());
    copy(orig);
}

KmerCounter::~KmerCounter() {

    deallocate();
}

int KmerCounter::getNumNucleotides() const {

    return _allNucleotides.size();
}

int KmerCounter::getK() const {

    return _k;
}

int KmerCounter::getNumKmers() const {

    return std::pow(getNumNucleotides(), getK());
}

int KmerCounter::getNumberActiveKmers() const {

    int actives = 0, rows = getNumRows(), cols = getNumCols();
    for (int r = 0; r < rows; r++)
        for (int c = 0; c < cols; c++)
            if (this->operator()(r, c) > 0)
                actives++;

    return actives;
}

std::string KmerCounter::toString() const {
    std::string outputString = _allNucleotides + " " + std::to_string(_k) + "\n";

    for (int row = 0; row<this->getNumRows(); row++) {
        for (int col = 0; col<this->getNumCols(); col++) {
            outputString += std::to_string((*this)(row, col)) + " ";
        }
        outputString += "\n";
    }

    return outputString;
}

void KmerCounter::increaseFrequency(const Kmer& kmer, int frequency) {
    
    if (kmer.toString().find_first_not_of(_allNucleotides) != std::string::npos)
        throw std::invalid_argument("void KmerCounter::increaseFrequency(const Kmer& kmer, int frequency): the given kmer contains an invalid nucleotide.");

    int row, column;
    getRowColumn(kmer, row, column);
    if (row != -1 && column != -1)
        this->operator()(row, column) += frequency;
}

KmerCounter& KmerCounter::operator=(const KmerCounter& orig) {

    if (&orig != this)
        copy(orig);

    return *this;
}

KmerCounter& KmerCounter::operator+=(const KmerCounter& kc) {

    if (kc.getK() != getK())
        throw std::invalid_argument("KmerCounter& KmerCounter::operator+=(const KmerCounter& kc): kc has a different k");

    if (kc._allNucleotides != _allNucleotides)
        throw std::invalid_argument("KmerCounter& KmerCounter::operator+=(const KmerCounter& kc): kc has a different set of nucleotides");

    int row = getNumRows(), cols = getNumCols();
    for (int r = 0; r < row; r++)
        for (int c = 0; c < cols; c++)
            this->operator()(r, c) += kc(r, c);

    return *this;
}

void KmerCounter::calculateFrequencies(const char* fileName) {

    std::ifstream input(fileName);
    if (!input) {
        input.close();
        throw std::ios_base::failure("void KmerCounter::calculateFrequencies(const char* fileName): fileName can not be opened");
    }
    
    initFrequencies();
    std::string inputString;
    input >> inputString;
    int k = getK();
    int kmers = inputString.size() - k + 1;
    for (int i = 0; i < kmers; i++) {
        Kmer kmer(inputString.substr(i, k));
        kmer.normalize(_validNucleotides);
        increaseFrequency(kmer);
    } 
}

Profile KmerCounter::toProfile() const {

    Profile profile;
    KmerFreq kmerFreq;
    int rows = getNumRows(), cols = getNumCols();
    for (int r = 0; r < rows; r++)
        for (int c = 0; c < cols; c++) {
            if (this->operator()(r, c) > 0) {
                kmerFreq.setFrequency(this->operator()(r, c));
                kmerFreq.setKmer(getKmer(r, c));
                profile.append(kmerFreq);
            }
        }

    return profile;
}

int KmerCounter::getNumRows() const {

    return std::pow(getNumNucleotides(), (getK() + 1) / 2);
}

int KmerCounter::getNumCols() const {

    return std::pow(getNumNucleotides(), getK() / 2);
}

int KmerCounter::getIndex(const std::string& kmer) const {
    int index = 0;
    int base = 1;

    for (size_t i = 0; i < kmer.size(); i++) {
        size_t pos = _allNucleotides.find(kmer[kmer.size() - i - 1]);
        if (pos == std::string::npos)
            return -1;
        index += pos * base;
        base *= _allNucleotides.size();
    }
    return index;
}

std::string KmerCounter::getInvertedIndex(int index, int nCharacters) const {
    std::string result(nCharacters, Kmer::MISSING_NUCLEOTIDE);

    for (int i = result.size(); i > 0; i--) {
        result[i - 1] = _allNucleotides[index % _allNucleotides.size()];
        index = index / _allNucleotides.size();
    }
    return result;
}

void KmerCounter::getRowColumn(const Kmer& kmer, int& row, int& column) const {

    if (kmer.toString().find_first_not_of(_allNucleotides) != std::string::npos) {
        row = -1;
        column = -1;
    }
    else {
        int k = getK();
        std::string kmer_row = kmer.toString().substr(0, (k + 1) / 2);
        std::string kmer_col = kmer.toString().substr((k + 1) / 2);
        row = getIndex(kmer_row);
        column = getIndex(kmer_col);
    }
}

Kmer KmerCounter::getKmer(int row, int column) const {
    
    if (row < 0 || row >= getNumRows()) {
        
        throw std::invalid_argument("Kmer KmerCounter::getKmer(int row, int column) const: the given row is out of the correct bounds");
    }
    
    if (column < 0 || column >= getNumCols()) {
        
        throw std::invalid_argument("Kmer KmerCounter::getKmer(int row, int column) const: the given column is out of the correct bounds");
    }
    
    return Kmer(getInvertedIndex(row, (getK() + 1) / 2) + getInvertedIndex(column, getK() / 2));
}

void KmerCounter::initFrequencies() {

    int rows = getNumRows(), cols = getNumCols();
    for (int r = 0; r < rows; r++)
        for (int c = 0; c < cols; c++)
            this->operator()(r, c) = 0;
}

const int& KmerCounter::operator()(int row, int column) const {

    return _frequency[row][column];
}

int& KmerCounter::operator()(int row, int column) {

    return _frequency[row][column];
}

void KmerCounter::allocate(int rows, int colums) {

    _frequency = new int* [rows];
    _frequency[0] = new int[rows * colums];
    for (int i = 1; i < rows; i++)
        _frequency[i] = _frequency[i - 1] + colums;
}

void KmerCounter::deallocate() {

    delete[] _frequency[0];
    delete[] _frequency;
    _frequency = nullptr;
    _k = 0;
    _validNucleotides = "";
    _allNucleotides = "";
}

void KmerCounter::copy(const KmerCounter& kmerCounter) {

    int rows = kmerCounter.getNumRows(), cols = kmerCounter.getNumCols();
    int** frequency = new int*[rows];
    frequency[0] = new int [rows * cols];
    for (int i = 1; i < rows; i++)
        frequency[i] = frequency[i - 1] + cols;

    for (int r = 0; r < rows; r++)
        for (int c = 0; c < cols; c++)
            frequency[r][c] = kmerCounter(r, c);

    delete[] _frequency[0];
    delete[] _frequency;
    _frequency = frequency;
    _k = kmerCounter.getK();
    _validNucleotides = kmerCounter._validNucleotides;
    _allNucleotides = kmerCounter._allNucleotides;
}