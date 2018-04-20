task clumps_prep {

    File huniprot2pdb
    Int scatterWidth
    File inMaf

    command <<<

        set -x

        ## EXTRACT THE MISSENSE MUTATIONS FROM MAF AND MAP THEM TO PROTEOME
        # get reference human proteome and genome
        wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640_9606.fasta.gz
        wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit
        python /sw/src/GPmapper.py ${inMaf}
        tar cfz mutfilesSplitByProtein.tar.gz splitByProtein/

        ## CALCULATE SAMPLE MUTATIONAL FREQUENCIES AND MUTATIONAL SPECTRA
        python /sw/src/calcSampleMutationFrequencies.py ${inMaf}
        python /sw/src/calcMutationContexts.py ${inMaf}

        ## SPLIT THE HUNIPROT2PDB MAP FILE
        huniprot2pdb_ungz=`echo ${huniprot2pdb} | sed -r "s/\.gz$//g"`
        zcat ${huniprot2pdb} > $huniprot2pdb_ungz
        ls -alh $huniprot2pdb_ungz
        split -d --number=l/${scatterWidth} -a 5 $huniprot2pdb_ungz huniprot2pdb_chunk_
        find huniprot2pdb_chunk_* -exec gzip {} \;

    >>>

    runtime {
        docker : "atanask/clumps_emprint@sha256:d0acd8c60d4f7bc1bea142f26ca8f5a370dfefcd47d259eca9c54d9a38e96a24"
        #docker : "clumps"
        disks: "local-disk 10 HDD"
        memory: "3 GB"  ## could be micro
        cpu: "1"
        preemptible: "3"
    }

    output {
        File mutations = "mutfilesSplitByProtein.tar.gz"
        File sampleMutFreq = "sampleMutFreq.txt"
        File sampleMutSpectra = "sampleMutSpectra.txt"
        Array[File] huniprot2pdbChunks = glob("huniprot2pdb_chunk_*")
    }
}




task clumps_test  {

    Int timeout
    Int nthreads
    File huniprot2pdbChunk
    Int lineId
    File mutationsTarball
    String ttype      ## tumor type. Use "PanCan" for a PanCancer analysis, otherwise specify one tumor type abbreviation that is present in your MAF
    String sampler    ## 'UniformSampler', 'CoverageSampler', 'MutspecCoverageSampler' ; this will be the simulator that generates the null model
    File sampleMutFreq
    File sampleMutSpectra
    File coveragetrack  ## whole-exome/genome coverage stats based on Panel of Normals (this file should be in some open-access bucket) 
    File coveragetrackidx  ## index file of the above 

    String cpu
    String memory
    String disks
    String preemptible


    command <<<

        set -x

        ## making sure we're writing to local disk (not boot disk)
        mkdir sw
        cp -vr /sw/* sw
        rm -rf /sw
        ln -vs $PWD/sw /sw
        find /sw
        
        ## process mutations tarball made by prep task
        tar xzvf ${mutationsTarball} && mv splitByProtein /sw/dat/

        ## get genome and preoteome references (needed for some of the samplers)
        wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit -P /sw/src/
        wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640_9606.fasta.gz -P /sw/src/

        ## download PDB files
        for struct in `zcat ${huniprot2pdbChunk} | cut -f 3 | cut -c 1-4 | sort | uniq` ; \
          do dirname=`echo $struct | cut -c 2-3` ; \
          wget -qr ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/$dirname/pdb$struct.ent.gz -P /sw/dat/ ; \
        done

        ## run clumps scan through all protein-structure pairs in this file
        if [ ${lineId} = -1 ]; then
            nlines=`zcat ${huniprot2pdbChunk} | wc -l`
            lines=$(seq 1 $nlines)
        else
            lines=${lineId}
        fi

        #echo 'sampleMutFreq', ${sampleMutFreq}
        #find ${sampleMutFreq}

        #echo "python clumps2.py ../set/14k.run3 ${timeout} ${nthreads} ${huniprot2pdbChunk} [LINE] ${ttype} ${sampler} ${sampleMutFreq} ${sampleMutSpectra} ${coveragetrack}"

        echo 'FIND COVERAGETRACK:'
        find ${coveragetrack}
        echo '#'
        find /cromwell_root/firecloud-tcga-open-access/tutorial/reference/PoNs/ReadCoveragePoNs

        START_DIR=$PWD
        cd /sw/src
        for line in $lines; \
          do python clumps2.py ../set/14k.run3 ${timeout} ${nthreads} ${huniprot2pdbChunk} $line ${ttype} ${sampler} ${sampleMutFreq} ${sampleMutSpectra} ${coveragetrack}; \
        done


        ## tar the results
        tar czf clumpsOut.tar.gz /sw/res
        mv -v clumpsOut.tar.gz $START_DIR
              
    >>>

    runtime {
        docker : "atanask/clumps_emprint@sha256:d0acd8c60d4f7bc1bea142f26ca8f5a370dfefcd47d259eca9c54d9a38e96a24"
        #docker : "clumps"
        disks: "${disks}"
        memory: "${memory}"
        cpu: "${cpu}"
        preemptible: "${preemptible}"
    }

    output {
        File clumpsOut = "clumpsOut.tar.gz"
    }
}


task clumps_nominate_candidates {
    
    Array[File] clumpsOut
    Array[String] chunkPaths

    command <<<

        set -x

        ## 
        for tar in ${sep=' ' clumpsOut}; do tar xzf $tar -C /; done
        tar czf clumpsOut.tar.gz /sw/res
        python /sw/src/nominateCandidates.py /sw/res/clumps ${sep=',' chunkPaths}

    >>>

    runtime {
        docker : "atanask/clumps_emprint@sha256:d0acd8c60d4f7bc1bea142f26ca8f5a370dfefcd47d259eca9c54d9a38e96a24"
        #docker : "clumps"
        disks: "local-disk 10 HDD"
        memory: "3 GB"
        cpu: "1"
        preemptible: "3"
    }

    output {
        File clumpsScanOut = "clumpsOut.tar.gz"
        Array[String] candidates_chunkpaths = read_lines("candidates_chunkpaths.txt")
        Array[Int] candidates_lineids = read_lines("candidates_lineids.txt")
        Array[Int] candidates_idx = read_lines("candidates_idx.txt")
    }
}


task clumps_make_report {
    File cancerGeneList
    String ttype

    File clumpsScanOut
    Array[File] clumpsCandidatesOut
    File huniprot2pdb
    File mutationsSplitByProtein

    command <<<

        set -x

        mkdir /scan
        tar xzvf ${clumpsScanOut} -C /scan
        mv /scan/sw/res/clumps /sw/res/
        
        find /sw/res/

        mkdir /cand
        for tar in ${sep=' ' clumpsCandidatesOut}; \
            do tar xzf $tar -C /cand; \
        done
        for fn in `ls -1 /cand/sw/res/clumps` ; \
            do cat /cand/sw/res/clumps/$fn >> /sw/res/clumps/$fn ; \
        done

        tar cfz clumpsConsolidatedOut.tar.gz /sw/res/clumps
        echo 'mutationsSplitByProtein', ${mutationsSplitByProtein}
        find ${mutationsSplitByProtein}
        tar xzvf ${mutationsSplitByProtein} && mv splitByProtein /sw/dat/
        python /sw/src/clumps_postprocess.py /sw/set/14k.run3 /sw/res/clumps ${cancerGeneList} ${huniprot2pdb} ${ttype}

    >>>

    runtime {
        docker: "atanask/clumps_emprint@sha256:d0acd8c60d4f7bc1bea142f26ca8f5a370dfefcd47d259eca9c54d9a38e96a24"
        #docker: "clumps"
        memory: "3 GB"
        cpu: "1"
        disks: "local-disk 10 HDD"
        preemptible: "3"
    }

    output {
        File clumpsConsolidatedOut = "clumpsConsolidatedOut.tar.gz"
        File resultsSummary = "clumps.summary.txt"
    }
}



workflow clumps_wf {
    call clumps_prep

    scatter (chunk in clumps_prep.huniprot2pdbChunks) {
        call clumps_test as clumps_test_scan {
            input:
                timeout = 15,
                nthreads = 1,
                huniprot2pdbChunk = chunk,
                lineId = -1,
                mutationsTarball = clumps_prep.mutations,
                sampleMutFreq = clumps_prep.sampleMutFreq,
                sampleMutSpectra = clumps_prep.sampleMutSpectra,
                cpu = "1",
                memory = "3 GB",
                disks = "local-disk 10 HDD",
                preemptible = "3"

        }
    }
    
    call clumps_nominate_candidates {
        input:
            chunkPaths = clumps_prep.huniprot2pdbChunks,
            clumpsOut = clumps_test_scan.clumpsOut

    }

    scatter (idx in clumps_nominate_candidates.candidates_idx) {
        call clumps_test as clumps_test_candidate {
            input:
                timeout = 0,
                nthreads = 16,
                huniprot2pdbChunk = clumps_nominate_candidates.candidates_chunkpaths[idx], 
                lineId = clumps_nominate_candidates.candidates_lineids[idx],
                mutationsTarball = clumps_prep.mutations,
                sampleMutFreq = clumps_prep.sampleMutFreq,
                sampleMutSpectra = clumps_prep.sampleMutSpectra,
                cpu = "16",
                memory = "5 GB",
                disks = "local-disk 10 HDD",
                preemptible = "3"

        }
    }

    call clumps_make_report {
        input:
            clumpsScanOut = clumps_nominate_candidates.clumpsScanOut,
            clumpsCandidatesOut = clumps_test_candidate.clumpsOut,
            mutationsSplitByProtein = clumps_prep.mutations
    }
}
