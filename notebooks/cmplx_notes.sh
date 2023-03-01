Complexes:
- Spliceosome --> U2-type spliceosome, LSM spliceosome, Small ribosomal subunit processome, Spliceosomal commitment complex, mRNA decay complex (UPF1, UPF2, UPF3B, DCP2, XRN1, XRN2, EXOSC2, EXOSC4, EXOSC10, PARN), C complex spliceosome, U2 small nuclear ribonucleoprotein complex
- Transcription factors --> Transcription factor complex TFIIA, Transcription factor complex TFIID, Transcription factor complex TFIIE, Transcription factor complex TFIIF, Transcription factor complex TFIIH
- Polycomb repressor complex --> PcG protein complex
- Chromatin remodeling, trithorax group proteins --> TrxG protein complex
- Fe-S clusters, MMS19 complex, CIA targeting complex --> CIA complex (cytosolic iron-sulfur cluster assembly)
- Voltage-gated ion channels (VGICs) --> VGIC complex
- Golgi to ER Traffic complex --> GET complex
- Adaptor complexes --> Adaptor complex AP-1, Adaptor complex AP-2
- tRNA multisynthetase complex --> 

***TODO: Sort out spliceosome and TrxG/chromatin remodeling naming schema; fix IFT-B complex, adaptor complex and RNP complex characterization labels

**NOTE: gamma secretase complex essential to Notch-like signaling, erlin proteins are ER lipid rafts, retromer complex is golgi

Characterization:
- "Known" = protein/protein family is known to interact with another protein in its cluster
- "Novel member" =  protein/protein family is a novel interactor with known protein complexes in its cluster
- "Uncharacterized" = interactions within this cluster of proteins are completely unknown or minimally characterized
- "Singleton" (?) = protein/protein families existing between clusters?


Complex labels:
c("60S cytosolic large ribosomal subunit", "40S cytosolic small ribosomal subunit", "Nuclear pore complex", "Mitochondrial respiratory chain complex I", "IFT-B complex", "IFT-A complex", "19S Proteasome", "26S Proteasome complex", "CCT complex", "ESCRT-III complex", "Exosome complex", "COP9 signalosome complex", "TRAPP complex", "Arp2/3 protein complex", "COPI vesicle coat complex")


### LRRC40 complex
KOG0472|LRRC40
 - 1 in human
 - 3 in tetts
 - 4 in chlre
KOG0046|Plastin/calponin
- 3 in human
- 1 in tettrs
- 1 in chlre
KOG1880|PP1R8
- 1 in human
- 2 in tetts
- 2 in chlre

# human
KOG0472|Q9H9A6|LRRC40
KOG1880|Q12972_HUMAN|PP1R8
KOG0046|P13796_HUMAN|LCP1
KOG0046|P13797|PLS3
KOG0046|Q14651|PLS1

for x in P13796 P13797 Q14651; do grep -A 15 $x human.fasta; done

scp -r lrrc40 rmcox@hfogcomp01.ccbb.utexas.edu:/stor/work/Marcotte/project/rmcox/leca/alphafold

# look for diffpep scripts
grep -rnw leca/ -e 'RAB35-IFT27_diffpep'
grep -rnw LECA_archive/ -e 'RAB35-IFT27_diffpep'