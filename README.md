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

## Urate-by-sex effect on gout

Copy phenotype file

```
d=$(mktemp -d)
echo "copying pheno file to $d"
cp /mnt/storage/private/mrcieu/data/ukbiobank/phenotypic/applications/15825/2019-05-02/data/data.33352.csv "$d"/
```

Run analysis

```
Rscript urate-gout.R
```