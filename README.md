# Measuring NOSH violation using IV-exposure variance effects

## Simulations

- sim7
    - NOSH violation
- sim8
    - Monotonicity violation
- sim9 
    - Power simulation for continuous outcome
- sim10
    - Multi-SNP MR simulation with constant IV-exposure effect
- pleiotropy
    - Explore the effect of horizonal pleiotropy on IV-exposure variance effects (which does not cause T1E inflation)
- sim11
    - Power simulation for binary outcome
- sim12
    - Multi-SNP MR simulation with varying IV-exposure effect

## Phenotypes

Copy pheno file to ramdisk

```sh
salloc --nodes=1 --cpus-per-task=21 --mem=80G --time=06:00:00 --partition=mrcieu
d=$(mktemp -d)
echo "copying pheno file to $d"
cp /mnt/storage/private/mrcieu/data/ukbiobank/phenotypic/applications/15825/2019-05-02/data/data.33352.csv "$d"/
```

Extract phenotypes

```sh
module load languages/r/3.6.0
Rscript pheno.R
```

## UK Biobank analysis

- dag
    - Causal diagram for lung function positive control
- smok
    - Smoking positive control
- mr
    - MR effects using instruments stratified by IV-exp effects
- pleiotropy
    - Effect of pleiotropy and LD with causal variant on IV-exp variance effects
- urate-gout
    - sex-stratified effect of urate on gout