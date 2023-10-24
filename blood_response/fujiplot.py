"""
Input preparation for custom Fujiplots.
"""

import pandas as pd
import seaborn as sns
from matplotlib.colors import ListedColormap
import matplotlib


def prepare_fuji_input(plot_clumps, category, filename, analyses_root,
                       empty_ranges=False, split_genes=True, rename_genes=None,
                       sort_traits='count',
                       process_category=True):
    """
    prepare fujiplot input files
    category could be 'cell_type' or 'ptb_name'
    filename can be something like 'clumped_quantile_celltype'
    sort_traits can be 'count' or 'pvalue' or 'beta'
    """

    plot_clumps = plot_clumps.copy(deep=True)

    plot_clumps['cell_type'] = plot_clumps.trait.str.split('_', expand=True)[[0]]
    plot_clumps['RANGES'] = plot_clumps.RANGES.str.strip('[]')
    if not empty_ranges:
        plot_clumps = plot_clumps.loc[plot_clumps.RANGES != '']

    else:
        plot_clumps.loc[plot_clumps.RANGES == '', 'RANGES'] = 'NA'

    fuji_input = \
        plot_clumps.rename(columns={
            category: 'CATEGORY',
            'traitname': 'TRAIT',
            'chr': 'CHR',
            'pos': 'BP',
            'ID': 'MARKER',
            'RANGES': 'GENE'
        })[['CATEGORY', 'TRAIT', 'CHR', 'BP', 'MARKER', 'GENE', 'P', 'BETA', 'OBS_CT']]

    if rename_genes:
        replace_idx = fuji_input.GENE.isin(rename_genes.keys())
        fuji_input.loc[replace_idx, 'GENE'] = fuji_input.loc[replace_idx].GENE.replace(rename_genes)

    if split_genes:
        # explicitly replace the split_genes value if it's a boolean
        if split_genes is True:
            split_genes = ','
        # split gene names
        fuji_input = fuji_input. \
            set_index(['CATEGORY', 'TRAIT', 'CHR', 'BP', 'MARKER', 'P', 'BETA', 'OBS_CT']). \
            GENE.str.split(split_genes). \
            explode(). \
            reset_index()

    trait_counts = pd.DataFrame(fuji_input.TRAIT.value_counts()).reset_index().rename(
        columns={'index': 'TRAIT', 'TRAIT': 'trait_ct'},
    )

    fuji_input = fuji_input.merge(trait_counts)

    if sort_traits == 'trait_count':
        # sort by the number of traits
        fuji_input = fuji_input.sort_values('trait_ct', ascending=False)
    elif sort_traits == 'pvalue':
        fuji_input = fuji_input.sort_values('P', ascending=True)
    elif sort_traits == 'beta':
        fuji_input = fuji_input.sort_values('BETA', ascending=False)
    fuji_input['trait_str'] = fuji_input.TRAIT

    # sorting conditions, then only keeping up to 3 conditions per gene
    # fuji_input = fuji_input.groupby('GENE').head(n=10)
    fuji_input = fuji_input.groupby(['MARKER', 'GENE', 'TRAIT']).head(n=3)
    fuji_input['CATEGORY'] = fuji_input.CATEGORY.str.replace(' ', '_')
    fuji_input['TRAIT'] = fuji_input.TRAIT.str.replace(' ', '_')
    fuji_input['TRAIT'] = fuji_input.TRAIT.str.replace('|', '_', regex=False)

    # TODO use single locus id if the gene is the same
    markers = list(fuji_input.MARKER.drop_duplicates())

    locus_ids = pd.DataFrame. \
        from_dict({'MARKER': markers}). \
        merge(fuji_input[['MARKER', 'GENE']])

    locus_ids = locus_ids.drop_duplicates()
    genes = [g for g in locus_ids.GENE.unique() if g != 'NA']
    gene_locus_ids = pd.DataFrame.from_records(enumerate(sorted(genes))).rename(columns={0: 'LOCUS_ID', 1: 'GENE'})
    unmapped_locus_ids = pd.DataFrame.from_records(
        enumerate(sorted(locus_ids.loc[locus_ids.GENE == 'NA'].MARKER.unique()))).rename(
        columns={0: 'UNMAPPED_LOCUS_ID', 1: 'MARKER'})
    gene_locus_ids['LOCUS_ID'] = gene_locus_ids.LOCUS_ID + 1
    # add the gene ids to the index
    unmapped_locus_ids['UNMAPPED_LOCUS_ID'] = unmapped_locus_ids.UNMAPPED_LOCUS_ID + gene_locus_ids.shape[0]

    locus_ids = locus_ids.merge(gene_locus_ids, how='left')
    locus_ids = locus_ids.merge(unmapped_locus_ids, how='left')
    locus_ids.loc[locus_ids.GENE == 'NA', 'LOCUS_ID'] = locus_ids.loc[locus_ids.GENE == 'NA', 'UNMAPPED_LOCUS_ID']
    locus_ids['LOCUS_ID'] = locus_ids.LOCUS_ID.astype(int)

    fuji_input = fuji_input.merge(locus_ids[['MARKER', 'GENE', 'LOCUS_ID']])

    traitlist = fuji_input[['CATEGORY', 'TRAIT', 'GENE']].drop_duplicates()
    traitlist = traitlist.sort_values('CATEGORY')
    traitlist['CATEGORY2'] = traitlist.CATEGORY
    if process_category:
        # TODO write this so it applies to conditions and cell types
        traitlist['CATEGORY2'] = traitlist.CATEGORY2.str.replace('Baseline', 'Baseline 0h')
        traitlist['CATEGORY2'] = traitlist.CATEGORY2.str.replace('[0-9]h', '', regex=True)
        traitlist['CATEGORY2'] = traitlist.CATEGORY2.str.replace('[0-9]', '', regex=True)
        traitlist['CATEGORY2'] = traitlist.CATEGORY2.str.replace('j-q', '')
        traitlist['CATEGORY2'] = traitlist.CATEGORY2.str.replace('s-q', '')
        traitlist['CATEGORY2'] = traitlist.CATEGORY2.str.replace('pca', 'joint')
        traitlist['CATEGORY2'] = traitlist.CATEGORY2.str.replace('pica', 'joint')
        traitlist['CATEGORY2'] = traitlist.CATEGORY2.str.replace('umap', 'joint')
        traitlist['CATEGORY2'] = traitlist.CATEGORY2.str.replace('__overnight', '', regex=False)
        traitlist['CATEGORY2'] = traitlist.CATEGORY2.str.replace('[_\\-\\.]+$', '', regex=True)
        traitlist['CATEGORY'] = traitlist.CATEGORY2

    category_colors = pd.DataFrame(traitlist.CATEGORY2.value_counts()).reset_index().rename(
        columns={'index': 'CATEGORY2', 'CATEGORY2': 'ct'},
    )

    #cmap = cm.get_cmap('tab20', category_colors.shape[0])  # PiYG
    sns_pal = sns.color_palette("colorblind", category_colors.shape[0])
    cmap = ListedColormap(sns_pal.as_hex())
    colors = []
    for i in range(cmap.N):
        rgba = cmap(i)
        # rgb2hex accepts rgb or rgba
        colors.append(matplotlib.colors.rgb2hex(rgba))
    category_colors['COLOR'] = colors

    traitlist = traitlist.merge(category_colors)
    traitlist = traitlist.sort_values(['ct', 'CATEGORY2'], ascending=True)
    traitlist.to_csv(f'{analyses_root}/circos/traitlist_{category}_{filename}.txt', sep='\t')
    fuji_input.loc[fuji_input.GENE != ''].to_csv(f'{analyses_root}/circos/input_{category}_{filename}.txt', sep='\t')

    print(f'Unique genes: {len(fuji_input.GENE.unique())}')
    print(f'Unique traits: {len(fuji_input.TRAIT.unique())}')

    traitlist_styled = traitlist.style.apply(lambda row: [f'background-color: {row.COLOR}'] * traitlist.shape[1],
                                             axis=1)
    return traitlist, fuji_input, traitlist_styled
