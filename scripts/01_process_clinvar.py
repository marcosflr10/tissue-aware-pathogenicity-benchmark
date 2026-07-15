"""
01_process_clinvar.py

Author: Marcos Flores Ruipérez

Description:
Processes ClinVar variants to build the balanced reference dataset used
for benchmarking pathogenicity predictors.

Workflow:
1. Load ClinVar database.
2. Filter GRCh38 single nucleotide variants.
3. Keep only Benign and Pathogenic variants.
4. Assign tissue using GTEx tissue-specific gene lists.
5. Validate tissue assignment using phenotype annotations.
6. Balance classes by random downsampling.
7. Export the final reference dataset.

Input:
- clinvar.txt.gz
- genes_brain_tissue_specific.txt
- genes_heart_tissue_specific.txt
- genes_liver_tissue_specific.txt

Output:
- clinvar_balanced.csv
"""

import pandas as pd

# Load ClinVar dataset

clinvar = pd.read_csv("clinvar.txt.gz", sep="\t", low_memory=False)


# Filter relevant variants

clinvar_filtered = clinvar[
    (clinvar["Assembly"] == "GRCh38") &
    (clinvar["ClinicalSignificance"].isin(["Pathogenic", "Benign"])) &
    (clinvar["Type"] == "single nucleotide variant")
].copy()

print(f"Variants after filtering: {len(clinvar_filtered)}")


# Load tissue-specific gene lists

brain_genes = set(open("genes_brain_tissue_specific.txt").read().splitlines())
heart_genes = set(open("genes_heart_tissue_specific.txt").read().splitlines())
liver_genes = set(open("genes_liver_tissue_specific.txt").read().splitlines())

def assign_tissue(gene):
    if gene in brain_genes:
        return "brain"
    elif gene in heart_genes:
        return "heart"
    elif gene in liver_genes:
        return "liver"
    else:
        return None
      
def assign_tissue_multi(gene_field):
    genes = str(gene_field).split(";")
    for g in genes:
        if g in brain_genes:
            return "brain"
        elif g in heart_genes:
            return "heart"
        elif g in liver_genes:
            return "liver"
    return None


# Assign tissue using gene lists

clinvar_filtered["tissue"] = clinvar_filtered["GeneSymbol"].apply(assign_tissue_multi)

clinvar_final = clinvar_filtered.dropna(subset=["tissue"])

print(f"Variants after GTEx assignment: {len(clinvar_final)}")


# Validate tissue assignment using phenotype

phenotype_to_tissue = {
    "cardio": "heart",
    "cardiomyopathy": "heart",
    "arrhythmia": "heart",
    
    "neuro": "brain",
    "epilepsy": "brain",
    "alzheimer": "brain",
    
    "liver": "liver",
    "hepat": "liver",
    "cirrhosis": "liver"
}

def map_phenotype_to_tissue(pheno):
    if pd.isna(pheno):
        return None
    
    pheno = pheno.lower()
    
    for key, tissue in phenotype_to_tissue.items():
        if key in pheno:
            return tissue
    
    return None

clinvar_final["tissue_pheno"] = clinvar_final["PhenotypeList"].apply(map_phenotype_to_tissue)

clinvar_final = clinvar_final[(clinvar_final["tissue"] == clinvar_final["tissue_pheno"])]

print(f"Variants after phenotype validation: {len(clinvar_final)}")


# Balance classes
pathogenic = clinvar_final[
    clinvar_final["ClinicalSignificance"] == "Pathogenic"
]

benign = clinvar_final[
    clinvar_final["ClinicalSignificance"] == "Benign"
]

print(f"Pathogenic variants: {len(pathogenic)}")
print(f"Benign variants: {len(benign)}")

n_samples = min(len(pathogenic), len(benign))

pathogenic = pathogenic.sample(
    n=n_samples,
    random_state=42
)

benign = benign.sample(
    n=n_samples,
    random_state=42
)

clinvar_balanced = pd.concat(
    [pathogenic, benign],
    ignore_index=True
)

print(f"Balanced dataset: {len(clinvar_balanced)}")

print("\nFinal distribution:")
print(clinvar_balanced["ClinicalSignificance"].value_counts())


# Save dataset
output_file = "clinvar_balanced.csv"

clinvar_balanced.to_csv(
    output_file,
    index=False
)

print(f"\nDataset successfully saved as '{output_file}'")
print("Process completed successfully.")


clinvar_balanced.to_csv("clinvar_balanced.csv", index=False)
print("Archivo guardado: clinvar_balanced.csv")
