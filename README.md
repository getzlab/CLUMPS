This is the code repository for CLUMPS.  Please cite this paper if using this for your research:


Atanas Kamburov, Michael S. Lawrence, Paz Polak, Ignaty Leshchiner, Kasper Lage, Todd R. Golub, Eric S. Lander, and Gad Getz. "Comprehensive assessment of cancer missense mutation clustering in protein structures" PNAS October 6, 2015 112 (40) E5486-E5495; first published September 21, 2015 https://doi.org/10.1073/pnas.1516373112 

Contents




Steps for running CLUMPS

1.) Generate a genome-proteome mapped MAF file

python GPmapper.py ${tumor type} ${MAF file} ${GPmapped file(output)}

2.) Calculate Mutation Contexts for each mutation in the MAF (Output will be named ${MAF file}.mutSpectra.txt)

python calcMutationContexts.py ${MAF file}

3.) Calculate Mutation Frequencies (Output will be named ${MAF file}.sampleMutFreq.txt)

python calcSampleMutationFrequencies.py ${MAF file}

4.) Start server to host mutations

python ../ana/loadMutations.py ${GPmapped file} 9014 &

5.) Start server to host structures

python ../ana/loadMaps.py

6.) Make set file for your cohort (see set directory for examples)


7.) For each structure, run CLUMPS:

for i in $(seq 0 102)
do 
    for j in $(seq 0 1000)
    do
	python clumps2.py ../set/${set file} 2 1 ../res/huniprot2pdb.run18.filt.split/huniprot2pdb.run18.filt.txt_${i}.gz ${j} ${ttype} ${sampler} ${sampleMutFreq file}  ${maf.mutSpectra file}  ${Fixed with binary describing coverage (necessary if using CoverageSampler or MutspecCoverageSampler)}
    done
done

8.) Once CLUMPS is run, aggregate individual results

python clumps_postprocess.py ../set/${set file} ../res/{Output directory, should be in set file} ../dat/allCancerGenes.txt ../res/huniprot2pdb.run18.filt.txt.gz ${ttype}
