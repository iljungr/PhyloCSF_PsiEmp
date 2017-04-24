This project contains python 2.7 software and data files for computing PhyloCSF-PsiEmp on Anopheles and Drosophila alignments. This is a comparative genomics method for determining whether a given multispecies genomic alignment is from a protein-coding genomic region.

PhyloCSf-PsiEmp is described in "Evolutionary dynamics of abundant stop codon readthrough"  by Irwin Jungreis, Clara S. Chan, Robert M. Waterhouse, Gabriel Fields, Michael F. Lin, Manolis Kellis.  Molecular Biology and Evolution (2016). doi:10.1093/molbev/msw189 https://academic.oup.com/mbe/article/33/12/3108/2450093/Evolutionary-Dynamics-of-Abundant-Stop-Codon

It is a high-specificity variant of PhyloCSF-Psi, which is described here: Lin, M. F., Jungreis, I., & Kellis, M. (2011). PhyloCSF: a comparative genomics method to distinguish protein coding and non-coding regions. Bioinformatics (Oxford, England), 27(13), i275–82. doi:10.1093/bioinformatics/btr209 https://academic.oup.com/bioinformatics/article/27/13/i275/178183/PhyloCSF-a-comparative-genomics-method-to

Computing PhyloCSf-PsiEmp for one or more regions of the AgamP3 Anopheles gambiae or dm3  Drosophila melanogaster genome assemblies requires extracting alignments, computing raw PhyloCSF scores, and computing PhyloCSf-PsiEmp from these.

Alignments can be obtained using the “Fasta Out” option in CodAlignView after setting "Ancestor" to "None". CodAlignView is available here: https://data.broadinstitute.org/compbio1/cav.php. Use Alignment Set “AgamP3” for 21 Anopheles genomes, “AgamP3_19” for 19 Anopheles genomes (all except the two most distantly-related Anopheles species, A. darlingi and A. albimanus), “dm3” for 12 Drosophila genomes, or “dm3_20” for 20 Drosophila genomes. The Drosophila alignments are also available from other sources.

To compute raw PhyloCSF scores, download and install PhyloCSF as described here: https://github.com/mlin/PhyloCSF/wiki. Use the 12 flies, 20flies, or 21mosquitoes parameters as appropriate (the last can also be used for 21 or 19 mosquitoes).

To compute PhyloCSf-PsiEmp in python 2.7, create an instance of the class PsiEmp.PsiEmpEvaluator by invoking the static "load" method with one of the supplied .cp files. This instance can then be repeatedly called with the raw PhyloCSF score and the length of the region, measured in codons.

Example:
The 37-codon region between the stop codon of A. gambiae transcript AGAP006474-RA and the next in-frame stop codon (not including either stop codon) has coordinates chr2L:32674673-32674695+chr2L:32674763-32674850 on the "+" strand. Its alignment can be seen in CodAlignView here: https://data.broadinstitute.org/compbio1/cav.php?controlsState=show&intervals=chr2L:32674673-32674695+chr2L:32674763-32674850&alnset=AgamP3&ancestor=None. The alignment can be extracted by pressing the Fasta Out button and saving the resulting output to a fasta file.

Feeding that alignment into PhyloCSF yields a raw score of 187.9951.

To compute PhyloCSf-PsiEmp in python 2.7:
>>> from PsiEmp import PsiEmpEvaluator
>>> agam21PsiEmpEval = PsiEmpEvaluator.load('PsiEmp.AgamP3.cp')
>>> agam21PsiEmpEval(187.9951, 37)
29.27089486447831
