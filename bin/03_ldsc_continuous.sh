cfg=${1:-config.yml}
REF=$(yq '.paths.ref' $cfg)
BFILE=${REF}/1000G_EUR_Phase3_plink/1000G.EUR.QC
HM3=${REF}/listHM3.txt
LDSC=external/ldsc/ldsc.py
# loop over metrics_bed/*/*.bed …
