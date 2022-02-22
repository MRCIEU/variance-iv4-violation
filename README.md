# Measuring violation of MR IV4 assumptions using variance effects

## Simulations

- Sim7
    - NOSH violation
- Sim8
    - Monotonicity violation
- Sim9 
    - Power simulation for continuous outcome
- Sim10
    - Power simulation for binary outcome

## UK Biobank analysis

### Phenotypes

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