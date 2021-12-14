#-*- coding: utf-8

# This workflow details our scripts for skin MAG generation
# This workflow is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# MAG Snakemake workflow is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with MAG Snakemake workflow.  If not, see <https://www.gnu.org/licenses/>.

#'''
#   This is a basic framework for recovery and basic quality control of MAGs
#   To visualize the pipeline: snakemake --dag | dot -Tpng > dag.png
#'''

__maintainer__ = "Sara Kashaf"
__email__ = "sskashaf@ebi.ac.uk"



import os
from os.path import join
import sys
import glob
import pandas as pd
import csv


outfiles_all=[]

#collect all the eukaryotic bins using glob
EUKS,=glob_wildcards("eukcc_output/eukcc_{euks}/eukcc.tsv")
EUKS_filt,=glob_wildcards("QC_fasta/fasta_comp_cont/{euks}")

#QC euks
outfiles_all.append(expand(join("QC_fasta/all_tsv_parsed/{bins}_eukcc_parsed.csv"),bins=EUKS))

#dereplicate
outfiles_all.append("dRep/data_tables/Sdb.csv")

#dnadiff
outfiles_all.append("MAG_genbank/dnadiff_summary.tsv")


rule all:
   input: outfiles_all


rule parse_tsv:
    input:
        tsv="eukcc_output/eukcc_{euks}/eukcc.tsv",
    output:
        tsv="QC_fasta/all_tsv/{euks}_eukcc.tsv",
    priority: 1000
    shell:"""
       scp {input.tsv} {output.tsv}
"""


rule make_csv:
    input: "QC_fasta/all_tsv/{euks}_eukcc.tsv"
    output: "QC_fasta/all_tsv_parsed/{euks}_eukcc_parsed.csv"
    params:
        run="{euks}",
        qc="QC_fasta/eukcc_summ.tsv",
        f="fastas/{euks}",
        odir="QC_fasta/fasta_comp_cont"
    run:
        os.system("mkdir -p " + str(params.odir))
        infile=str(input)
        outfile=str(output)
        qc=str(params.qc)
        with open(infile) as inp:
             LoL=[x.strip().split('\t') for x in inp]
             row1=LoL[1]
             completeness=LoL[1][0]
             contamination=LoL[1][1]
             run=str(params.run)
             g = open(outfile, "w")
             if float(completeness)>=50 and float(contamination)<5:
                  g.write(run+','+completeness+','+contamination+'\n')
                  g.close()
                  h = open(params.qc, "a")
                  h.write(run+','+completeness+','+contamination+'\n')
                  h.close()
                  f = str(params.f)
                  odir=str(params.odir)
                  os.system("scp "+f+ " "+str(odir))
             else:
                  g.write("")
                  g.close()



rule make_genome_info:
   input: "QC_fasta/eukcc_summ.tsv"
   output: "QC_fasta/eukcc_metrics.csv"
   shell: "echo -e 'genome,completeness,contamination' | cat - {input} > {output}"



rule dereplicate:
   input: 
      parsed=expand("QC_fasta/all_tsv_parsed/{euks}_eukcc_parsed.csv",euks=EUKS),
      genome_info="QC_fasta/eukcc_metrics.csv"
   output: "dRep/data_tables/Sdb.csv"
   params:
      infolder="QC_fasta/fasta_comp_cont",
      outfolder="dRep"
   threads: 20
   shell:"""
        dRep dereplicate -p {threads} \
            {params.outfolder} -g {params.infolder}/*.fa \
            -pa 0.9 -sa 0.95 -nc 0.3 \
            -cm larger \
            --genomeInfo {input.genome_info} \
            -comp 50 -con 5 
"""



rule mash_dist:
   input: 
      bins="QC_fasta/fasta_comp_cont/{euks}",
      db="latest/mash/genbank-fungi.msh"
   output:
      join("MAG_genbank/mashdist/{euks}.tab"),
   threads: 1
   shell:
        "mash dist -p {threads} {input.db} {input.bins} > {output}"


rule best_mash:
    input:
        mashdist=join(
            "MAG_genbank/mashdist/{euks}.tab"
        ),
    output: "MAG_genbank/best_mash/{euks}.tab"
    threads: 1
    shell:
        """
        sort -gk3 {input.mashdist}|sed -n 1p >{output}
        """


rule dnadiff:
    input:
        bestmash="MAG_genbank/best_mash/{euks}.tab"
    output: "MAG_genbank/dnadiff/{euks}.report"
    params:
        outdir=join("MAG_genbank/dnadiff"),
        bins="{euks}",
        genomesdir="MAG_genbank/dnadiff/genomes"
    shell:
        """
        mkdir -p {params.genomesdir}
        while read col1 col2 rem
          do
            f="$(basename -- ${{col1}})"
            echo {params.genomesdir}/${{f%%.gz}}
            if [ ! -f {params.genomesdir}/${{f%%.gz}} ]; then
               echo "copying file over"
               scp ${{col1}} {params.genomesdir}
               gunzip {params.genomesdir}/${{f}}
            fi  
            echo 'dnadiff ${{col1}} ${{col2}} -p ${{col1%%.fasta}}_${{col2%%.fa}}_'
            dnadiff {params.genomesdir}/${{f%%.gz}} ${{col2}} -p {params.outdir}/{params.bins}
          done < {input.bestmash}
        """
 

rule parse_dnadiff:
    input:
        dnadiff=join("MAG_genbank/dnadiff/{euks}.report")
    output:
        join(
          "MAG_genbank/dnadiff_parsed/{euks}.tsv"
        ),
    run:
        outfile = str(output)
        f = open(input.dnadiff)
        data = f.read()
        first_line = data.split("\n", 1)[0]
        a = first_line.split(" ")
        ref = a[0]
        quer = a[1]
        with open(outfile, "w") as outf:
            path_dna = input.dnadiff
            base = os.path.basename(path_dna)
            base = base.split(".report")[0]
            with open(path_dna) as f:
                for line in f:
                    if "TotalBases" in line:
                        cols = line.split()
                        lenref = int(cols[1])
                        lenquer = int(cols[2])
                    if "AlignedBases" in line:
                        cols = line.split()
                        aliref = cols[1].split("(")[-1].split("%")[0]
                        alique = cols[2].split("(")[-1].split("%")[0]
                    if "AvgIdentity" in line:
                        cols = line.split()
                        ident = float(cols[1])
            line = "%s\t%s\t%i\t%.2f\t%i\t%.2f\t%.2f" % (ref, quer, lenref, float(aliref), lenquer, float(alique), float(ident))
            outf.writelines(line + "\n")



rule aggregate_dnadiff:
    input: expand("MAG_genbank/dnadiff_parsed/{euks}.tsv",euks=EUKS_filt)
    output: "MAG_genbank/dnadiff_summary.tsv"
    shell:
        """
        cat {input}>{output}
        """
