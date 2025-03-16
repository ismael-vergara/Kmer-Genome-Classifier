# 🧬 Kmer - Genome Classifier

<p align="center">
  <img src="https://img.shields.io/badge/C%2B%2B-11-blue?style=for-the-badge&logo=c%2B%2B&logoColor=white">
  <img src="https://img.shields.io/badge/Linux/Mac/Windows-Compatible-green?style=for-the-badge">
  <img src="https://img.shields.io/badge/Status-Completed-brightgreen?style=for-the-badge">
</p>

## 📋 Project Overview
This project was developed as part of the **"Programming Methodology"** subject of the **Computer Science and Mathematics** dual degree in the University of Granada.

The main goal is to create a program capable of **analyzing and comparing genomes** from different species to identify similarities and differences between them.

The project is based on the analysis of **K-mers**, which are nucleotide sequences of length *k* extracted from DNA or RNA sequences. By analyzing the frequency and distribution of K-mers, this tool generates species profiles that can be used to classify unknown genomic sequences.

---

## 🚀 Main Features
- ✅ **K-mer Extraction:** Processes genomic sequences to extract K-mers of variable length.  
- ✅ **Species Profile Generation:** Builds species profiles (*Profile*) based on K-mer frequencies using the `LEARN` executable.  
- ✅ **Classification of Unknown Genomes:** Compares unknown genome sequences to known species profiles using the `CLASSIFY` executable.  
- ✅ **Binary and Text Input/Output:** Supports both text and binary file formats for profile management.  
- ✅ **Operator Overloading and Comparison:** Implements comparison operators for K-mers and K-mer frequencies to facilitate sorting and analysis.  

---

## 📦 Modules and Classes
### 📌 `Kmer`
Represents K-mers and handles their **manipulation, normalization, and I/O operations** (binary and text).

### 📌 `KmerFreq`
Associates a K-mer with its **frequency of appearance**. Includes comparison operators for sorting.

### 📌 `Profile`
Represents a species' **genomic profile** — a set of K-mers and their frequencies. Supports merging and manipulation of profiles.

### 📌 `KmerCounter`
Counts K-mers in DNA/RNA sequences and **builds profiles** from multiple individuals.

---

## 💻 Executables
### 🔬 `LEARN`
Creates a **genomic profile** for a species from a set of genome files.

#### **Usage:**
```
LEARN [-t|-b] [-p profileId] [-k kValue] [-n nucleotidesSet] [-o outputFile] input1.dna [input2.dna ...]
```

#### **Options:**
- `-t` → Output in **text mode**.
- `-b` → Output in **binary mode**.
- `-p` → Profile **ID** (species name, default: "unknown").
- `-k` → Length of **K-mers** (default: 5).
- `-n` → Valid **nucleotides** (default: "ACGT").
- `-o` → Output **file name** (default: "output.prf").

### 🔍 `CLASSIFY`
Classifies an **unknown genome sequence** by comparing it to one or more known species profiles.

#### **Usage:**
```sh
CLASSIFY unknown.dna profile1.prf [profile2.prf ...]
```

---

## 📂 File Structure
```bash
📦 Project Root
├── 📂 src/            # Source code files
│   ├── 📜 Kmer.h / Kmer.cpp
│   ├── 📜 KmerFreq.h / KmerFreq.cpp
│   ├── 📜 Profile.h / Profile.cpp
│   ├── 📜 KmerCounter.h / KmerCounter.cpp
│   ├── 📜 LEARN.cpp
│   ├── 📜 CLASSIFY.cpp
├── 📂 data/          # Example genome files
├── 📂 output/        # Generated profiles
└── 📜 README.md
```
---

## ⚙️ Compilation and Execution

### 🔧 **Compilation with g++**
```sh
g++ Kmer.cpp KmerFreq.cpp Profile.cpp KmerCounter.cpp LEARN.cpp -o learn
g++ Kmer.cpp KmerFreq.cpp Profile.cpp KmerCounter.cpp CLASSIFY.cpp -o classify
```

### 🚀 **Example Usage**
#### ✅ Generate a species profile:
```sh
./learn -p Human -k 5 -o human.prf human1.dna human2.dna
```
#### ✅ Classify an unknown genome:
```sh
./classify unknown.dna human.prf chimp.prf virus.prf
```

---

## ✅ Requirements
🔹 **C++11 or higher**  
🔹 **NetBeans** (optional, for project management)  
🔹 **Compatible with Linux, MacOS, and Windows** (with minor adjustments)  

---

## 🌐 Example Use Case
Given genome sequences of Homo sapiens, Chimpanzee, and SARS-CoV-2 (Covid-19), this tool can generate profiles for each species. Later, an unknown genome can be classified by comparing its profile to the existing profiles, identifying which species it is most similar to.
This technique is useful in fields like bioinformatics, genomics research, and virology for species identification and genetic analysis.

---

## 🎓 Educational Context
This project was designed for academic purposes within the "Programming Methodology" subject, helping me practice:

🎯 Object-Oriented Programming (OOP)  
🎯 Dynamic memory management  
🎯 Text and binary file handling  
🎯 Operator overloading in C++  
🎯 Modular and scalable software design for scientific computing  
