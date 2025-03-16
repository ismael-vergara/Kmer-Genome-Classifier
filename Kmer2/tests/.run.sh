touch tests//.timeout
CMD="valgrind --leak-check=full /home/ismael/Escritorio/MIMP/NetBeansProjects/Kmer2/dist/Debug/GNU-Linux/kmer2  tests/output/4pairs_5pairsDNA.prf ../Genomes/4pairsDNA.prf ../Genomes/5pairsDNA.prf 1> tests//.out12 2>&1"
eval $CMD
rm tests//.timeout
