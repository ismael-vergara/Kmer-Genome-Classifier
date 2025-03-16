touch tests//.timeout
CMD="valgrind --leak-check=full /home/ismael/Escritorio/MIMP/NetBeansProjects/Kmer3/dist/Debug/GNU-Linux/kmer3  ../Genomes/human1_k1.prf ../Genomes/human1.prf ../Genomes/human2.prf 1> tests//.out11 2>&1"
eval $CMD
rm tests//.timeout
