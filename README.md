
### Quick install and start ###
```
git clone --recurse-submodules https://github.com/zhangrengang/Centromics
cd Centromics

# install
conda env create -f Centromics.yaml
conda activate RepCent
./install.sh

# start
cd example_data
# long reads
centromics -l hifi.fq.gz -g ref.fa

# long reads + HiC data + ChIP data
centromics -l hifi.fq.gz -g ref.fa -pre hifi -chip chip.bam -hic merged_nodups.txt.gz
centromics -l ont*.fq.gz -g ref.fa -pre ont  -chip chip.bam -hic merged_nodups.txt.gz
```
