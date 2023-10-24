"""
Read genotype information from plink (bim, fam, bed) files.
"""

import numpy as np
import pandas as pd
from pandas_plink import read_plink1_bin, read_rel


def snp_query(genotypes, alleles, samples, selected_snps):
    r"""

    use to annotate samples for a small genomic region. example for CASP3:

    `
    from pandas_plink import read_plink
    from blood_response import plink
    alleles, samples, genotypes = read_plink(f'{gwas_dir}/imputed_merged')

    casp3_alleles_idx = (alleles.chrom == '4') & (alleles.pos >= 185560000) & (alleles.pos <= 185580000)
    genotypes_casp3, alleles_casp3, samples_casp3 = plink.snp_query(genotypes, alleles, samples, casp3_alleles_idx)
    `

    Parameters
    ----------
    genotypes : (`n_snps`, `n_inds`) array
        Genetic data
    alleles : pandas.DataFrame
        Variant annotation
    samples: pandas.DataFrame
        Sample annotation
    selected_snps : bool array
        Variant filter
    Returns
    -------
    genotypes_out : (`n_snps`, `n_inds`) array
        filtered genetic data
    alleles_out : dataframe
        filtered variant annotation
    samples_out : dataframe
        sample annotation with added genotype data
    """
    alleles_out = alleles[selected_snps].reset_index(drop=True)
    genotypes_out = genotypes[alleles_out.i.values].compute()
    alleles_out.i = pd.Series(np.arange(alleles_out.shape[0]), index=alleles_out.index)
    samples_out = samples.copy()

    # iterate over snps, map to allele values and add to sample annotation
    for snp_i in range(genotypes_out.shape[0]):
        snp = alleles_out.iloc[snp_i]
        geno_dict = {
            0: f'{snp.a0}{snp.a0}',
            1: f'{snp.a0}{snp.a1}',
            2: f'{snp.a1}{snp.a1}',
        }
        samples_out[snp.snp] = [geno_dict[i] if np.isfinite(i) else np.nan for i in genotypes_out[snp_i, :]]
    return genotypes_out, alleles_out, samples_out
