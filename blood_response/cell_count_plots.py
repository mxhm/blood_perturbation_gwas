"""
Plots of cytometry data (2d histograms) for a single sample, multiple samples or comparisons between samples.
"""

import pandas as pd
from plotnine import *
from tqdm.auto import tqdm
from scipy.sparse import csc_matrix, coo_matrix
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt

from .plink import snp_query
from .gwas_preprocessing import standardize_study_id


def plot_scatter(x, y, c, name, marker=None, size=10, alpha=0.5):
    plt.scatter(x, y, c=c, marker=marker, cmap='Spectral', alpha=alpha, s=size)
    plt.gca().set_aspect('equal', 'datalim')
    plt.colorbar()
    plt.title(name, fontsize=12)
    # plt.show()


def generate_boxplots(
        snps_df,
        metadata,
        pheno_column,
        genotypes, alleles, samples
        ):

    """

    use like this:
    picas_lead_snps = all_lead_snps. \
      loc[all_lead_snps.projection == 'picas_all']. \
      drop_duplicates(['perturbation', 'projection', 'uniqID']). \
      sort_values('P')

    save_as_pdf_pages(
        generate_boxplots(picas_lead_snps, 'picas_all'),
        filename=f'{analysis_path}/boxplots_picas_plink_minp.pdf')

    :param snps_df:
    :param alleles:
    :param samples:
    :param metadata:
    :return:
    """
    gplts = []
    pheno_snp_list = []
    pheno_plot_list = []
    for idx, lead_snp in tqdm(list(snps_df.iterrows())):
        sel_snps = ((alleles.chrom == str(lead_snp.chr)) & (alleles.pos == lead_snp.pos))
        genotypes_, alleles_, samples_ = snp_query(genotypes, alleles, samples, sel_snps)
        # renaming the SNP column from the SNP id to 'genotype'
        samples_.columns = [c for c in samples_.columns[:-1]] + ['genotype']
        subject_snp_data = samples_[['iid', 'genotype']]
        if subject_snp_data.shape[0] > 0:
            meta_snp = metadata.merge(subject_snp_data, left_on='IID', right_on='iid')
            pheno_snp = meta_snp[['genotype', pheno_column]]
            pheno_snp.loc[pheno_snp.genotype.isna(), 'genotype'] = 'NA'
            pheno_snp['genotype'] = pheno_snp.genotype.str.replace('/', '')
            pheno_snp_m = pd.melt(pheno_snp, id_vars='genotype', var_name='phenotype')
            title = f'{lead_snp.locus}'
            pheno_plot = pheno_snp_m.dropna()
            snp_counts = pheno_plot.genotype.value_counts()
            snp_counts = snp_counts[snp_counts >= 3]
            snp_counts_dict = {c[0]: c[1] for c in snp_counts.items() if c[0] != 'NA'}
            if 'NA' in snp_counts.index:
                snp_counts_dict['NA'] = snp_counts['NA']
            snp_legend = [f'{i} ({c})' for i, c in snp_counts_dict.items()]

            pheno_plot = pheno_plot.loc[pheno_plot.genotype.isin(snp_counts_dict)]
            pheno_plot.genotype = pheno_plot.genotype.astype('category')
            pheno_plot.genotype = pheno_plot.genotype.cat.set_categories(snp_counts_dict.keys())
            pheno_plot['locus'] = lead_snp.locus
            #pheno_plot['chr'] = lead_snp.chr
            #pheno_plot['pos'] = lead_snp.pos

            # add all information from the snps_df
            pheno_plot = pheno_plot.merge(snps_df)

            gplt = (ggplot(pheno_plot, aes('genotype', 'value', fill='genotype', color='genotype')) +
                    geom_sina(size=2, shape='.', alpha=1) +
                    geom_violin(draw_quantiles=[0.25, 0.5, 0.75], color='black', fill=None) +
                    scale_fill_hue(labels=snp_legend) +
                    scale_color_hue(labels=snp_legend) +
                    labs(title=title) +
                    theme_minimal())

            gplts.append(gplt)
            pheno_snp_list.append(pheno_snp)
            pheno_plot_list.append(pheno_plot)
    return gplts, pheno_snp_list, pheno_plot_list


# # plot raw and transformed values for all dimensions
# def generate_boxplots_projected(snps_df,
#                       g, bim, fam,
#                       metadata,
#                       projection,
#                       platform='XN10'):
#     """
#
#     use like this:
#     picas_white_lead_snps = all_lead_snps. \
#       loc[all_lead_snps.projection == 'picas_all']. \
#       drop_duplicates(['perturbation', 'projection', 'uniqID']). \
#       sort_values('P')
#
#     save_as_pdf_pages(
#         generate_boxplots(picas_lead_snps, 'picas_all'),
#         filename=f'{analysis_path}/boxplots_picas_plink_minp.pdf')
#
#     :param snps_df:
#     :param projection:
#     :param bim:
#     :param fam:
#     :param metadata:
#     :return:
#     """
#     for idx, lead_snp in tqdm(list(snps_df.iterrows())):
#
#         sel_snps = ((bim.chrom == str(lead_snp.chr)) & (bim.pos == lead_snp.pos))
#         g_, bim_, fam_ = snp_query(g, bim, fam, sel_snps)
#         subject_snp_data = fam_[['iid', 'genotype']]
#         if subject_snp_data.shape[0] > 0:
#             meta_snp = metadata.merge(subject_snp_data, left_on='IID', right_on='iid')
#             pheno_columns = list(meta_snp.columns[
#                                      meta_snp.columns.str.startswith(f'ptb{lead_snp.perturbation_id}-{platform}')])
#
#             pheno_snp = meta_snp[['genotype'] + pheno_columns]
#             pheno_snp.loc[pheno_snp.genotype.isna(), 'genotype'] = 'NA'
#             pheno_snp['genotype'] = pheno_snp.genotype.str.replace('/', '')
#             pheno_snp_m = pd.melt(pheno_snp, id_vars='genotype', var_name='phenotype')
#             pheno_snp_m[['perturbation', 'platform', 'dimension', 'transformation']] = \
#                 pheno_snp_m.phenotype.str.split('-', expand=True)
#
#             title = f'{lead_snp.perturbation}|{lead_snp.algorithm}{lead_snp.lat_dim}|{lead_snp.uniqID}|{lead_snp.p:.1}|{lead_snp.re}'
#             pheno_plot = pheno_snp_m.loc[pheno_snp_m.dimension.str.startswith(projection)].dropna()
#             snp_counts = pheno_plot.genotype.value_counts()
#             snp_counts = (snp_counts /
#                           (len(pheno_plot.dimension.unique()) *
#                            len(pheno_plot.transformation.unique()))).astype(int)
#             snp_counts = snp_counts[snp_counts >= 3]
#             snp_counts_dict = {c[0]: c[1] for c in snp_counts.items() if c[0] != 'NA'}
#             if 'NA' in snp_counts.index:
#                 snp_counts_dict['NA'] = snp_counts['NA']
#             snp_legend = [f'{i} ({c})' for i, c in snp_counts_dict.items()]
#
#             pheno_plot = pheno_plot.loc[pheno_plot.genotype.isin(snp_counts_dict)]
#             pheno_plot.genotype = pheno_plot.genotype.astype('category')
#             pheno_plot.genotype = pheno_plot.genotype.cat.set_categories(snp_counts_dict.keys())
#
#             gplt = (ggplot(pheno_plot, aes('genotype', 'value', fill='genotype', color='genotype')) +
#                     geom_sina(size=2, shape='.', alpha=1) +
#                     geom_violin(draw_quantiles=[0.25, 0.5, 0.75], color='black', fill=None) +
#                     scale_fill_hue(labels=snp_legend) +
#                     scale_color_hue(labels=snp_legend) +
#                     facet_grid('transformation~dimension', scales='free_y') +
#                     labs(title=title) +
#                     theme_minimal())
#
#             yield gplt


# @lru_cache(maxsize=None)
def numpy_to_df(indices,
                sysmex_raw,
                channels=['WNR', 'WDF', 'RET', 'PLT-F'],
                views=[['Side Fluorescence Signal', 'Forward Scatter Signal'],
                       ['Side Fluorescence Signal', 'Side Scatter Signal'],
                       ['Side Scatter Signal', 'Forward Scatter Signal']],
                pseudocount=None, gaussian_sigma=None):
    """
    converts matrix format input to a long dataframe (across channels and views)
    """
    plot_data = sysmex_raw[indices, :, :, :].sum(axis=0)
    all_counts_df = pd.DataFrame()
    i = 0
    for channel in channels:
        for view in views:
            count_coo = csc_matrix(plot_data[:, :, i]).tocoo(copy=False)
            if pseudocount is not None:
                count_coo = coo_matrix(count_coo.todense() + pseudocount)
            if gaussian_sigma is not None:
                count_coo = coo_matrix(gaussian_filter(count_coo.todense(), sigma=gaussian_sigma))

            # flip two of the views to better match the sysmex plots
            if view in [['Side Fluorescence Signal', 'Forward Scatter Signal'],
                        ['Side Scatter Signal', 'Forward Scatter Signal']]:
                view = view[::-1]
                count_df = pd.DataFrame(
                    {'x': count_coo.row,
                     'y': count_coo.col,
                     'count': count_coo.data})[['x', 'y', 'count']]. \
                    sort_values(['x', 'y']). \
                    reset_index(drop=True)
            else:
                count_df = pd.DataFrame(
                    {'x': count_coo.col,
                     'y': count_coo.row,
                     'count': count_coo.data})[['x', 'y', 'count']]. \
                    sort_values(['x', 'y']). \
                    reset_index(drop=True)

            # add normalized count
            count_df['count_norm'] = count_df['count'] / count_df['count'].sum()
            count_df['channel'] = channel
            count_df['view'] = 'x: ' + view[1].replace('Signal', '') + '\n' + 'y: ' + view[0].replace('Signal', '')
            count_df['xlab'] = view[1]
            count_df['ylab'] = view[0]
            all_counts_df = all_counts_df.append(count_df)
            i += 1
    return all_counts_df


def raw_plot(sysmex_raw,
             metadata,
             plot_type,
             group,
             genotypes, alleles, samples,
             chrom=None,
             pos=None,
             gaussian_sigma=3,
             pseudocount=None,
             title=""
             ):
    if chrom is not None:
        title = f'{title}{chrom}:{pos}'
        sel_snps = ((alleles.chrom == str(chrom)) & (alleles.pos == pos))
        genotypes_, alleles_, samples_ = snp_query(genotypes, alleles, samples, sel_snps)
        # TODO handling only single SNP for plots (the last one)
        samples_.columns = [c for c in samples_.columns[:-1]] + ['genotype']
        subject_snp_data = samples_[['iid', 'genotype']]
        # the patient_sample_id is the cleaned up version
        subject_snp_data['standard_study_id'] = subject_snp_data['iid'].apply(standardize_study_id)
        metadata = metadata.merge(subject_snp_data, left_on='patient_sample_id', right_on='standard_study_id')

    # returning empty objects if only plotting totals
    count_df_dict = {}
    wide_count_df = None

    def wide_counts_merge(count_df_dict, keys):
        assert len(keys) == 2
        wide_count_df = count_df_dict[keys[0]]. \
            merge(count_df_dict[keys[1]], on=['x', 'y', 'channel', 'view'])
        wide_count_df['count_norm_ratio'] = wide_count_df.count_norm_x / wide_count_df.count_norm_y
        wide_count_df['count_norm_diff'] = wide_count_df.count_norm_x - wide_count_df.count_norm_y
        wide_count_df['count_norm_diff_abs'] = wide_count_df.count_norm_diff.abs()
        wide_count_df = wide_count_df.merge(
            wide_count_df.groupby(['channel', 'view']).count_norm_diff_abs.max().reset_index(name='max_diff')
        )
        wide_count_df['count_norm_diff_scaled'] = wide_count_df.count_norm_diff / wide_count_df.max_diff
        return wide_count_df

    if metadata.shape[0] > 0:
        if plot_type == 'total':
            plot_df = numpy_to_df(tuple(metadata.idx), sysmex_raw,
                                  pseudocount=pseudocount, gaussian_sigma=gaussian_sigma)
            gplt = (ggplot(plot_df, aes('x', 'y', fill='count')) +
                    geom_tile() +
                    scale_fill_continuous(trans='log10') +
                    facet_grid('channel~view') +
                    theme_minimal(base_size=8))
            # theme_gray(base_size=7))

        elif plot_type == 'ratio' or plot_type == 'difference':
            groups = metadata[group].unique()
            # keep only homozygous genotypes for now
            # could also try some sort of additive model
            if group == 'genotype':
                groups = set(groups).intersection({'AA', 'CC', 'GG', 'TT'})
            for grp in groups:
                metadata_ptb = metadata.loc[metadata[group] == grp]
                plot_df = numpy_to_df(tuple(metadata_ptb.idx),
                                      pseudocount=pseudocount, gaussian_sigma=gaussian_sigma)
                plot_df['group'] = grp
                count_df_dict[grp] = plot_df
            groups = tuple(count_df_dict.keys())
            wide_count_df = wide_counts_merge(count_df_dict, groups)

        if group is not None:
            assert len(groups) == 2, 'Need exactly two groups for comparison plots.'
            if plot_type == 'ratio':
                legend_name = f'ratio' \
                              f'\n{wide_count_df.group_x[0]} / ' \
                              f'{wide_count_df.group_y[0]}\n'
                gplt = (ggplot(wide_count_df, aes('x', 'y', fill='count_norm_ratio')) +
                        geom_tile() +
                        scale_fill_gradient2(trans='log10', name=legend_name) +
                        facet_grid('channel~view') +
                        labs(title=title) +
                        theme_gray(base_size=7))

            elif plot_type == 'difference':
                legend_name = f'scaled difference' \
                              f'\n{wide_count_df.group_x[0]} - ' \
                              f'{wide_count_df.group_y[0]}\n'
                gplt = (ggplot(wide_count_df, aes('x', 'y', fill='count_norm_diff_scaled')) +
                        geom_tile() +
                        scale_fill_gradient2(trans='pseudo_log', name=legend_name) +
                        facet_grid('channel~view') +
                        labs(title=title) +
                        theme_gray(base_size=7))
            else:
                raise NotImplementedError

    return count_df_dict, plot_df, wide_count_df, gplt


def filter_plot_df(plot_df, channel_views):
    return pd.concat([
        plot_df.loc[(plot_df.channel == channel) & (plot_df.xlab == view[0]) & (plot_df.ylab == view[1])]
        for channel, view in channel_views.items()
    ])


def raw_plot_subset(sysmex_raw,
                    metadata,
                    plot_type,
                    group,
                    genotypes, alleles, samples,
                    chrom=None,
                    pos=None,
                    channel_views={
                        'WNR': ['Side Fluorescence Signal', 'Forward Scatter Signal'],
                        'WDF': ['Side Scatter Signal', 'Side Fluorescence Signal'],
                        'RET': ['Side Fluorescence Signal', 'Forward Scatter Signal'],
                        'PLT-F': ['Side Fluorescence Signal', 'Forward Scatter Signal']},
                    gaussian_sigma=2,
                    pseudocount=None,
                    title="",
                    legend_position=None
                    ):
    if chrom is not None:
        title = f'{title}{chrom}:{pos}'
        sel_snps = ((alleles.chrom == str(chrom)) & (alleles.pos == pos))
        genotypes_, alleles_, samples_ = snp_query(genotypes, alleles, samples, sel_snps)
        # TODO handling only single SNP for plots (the last one)
        samples_.columns = [c for c in samples_.columns[:-1]] + ['genotype']
        subject_snp_data = samples_[['iid', 'genotype']]
        # the patient_sample_id is the cleaned up version
        subject_snp_data['standard_study_id'] = subject_snp_data['iid'].apply(standardize_study_id)
        metadata = metadata.merge(subject_snp_data, left_on='patient_sample_id', right_on='standard_study_id')

    # returning empty objects if only plotting totals
    count_df_dict = {}
    wide_count_df = None

    def wide_counts_merge(count_df_dict, group_count_dict):
        assert len(group_count_dict.keys()) == 2
        # this sorts the dict by value, with the less frequent group appearing first - this is the effect/treatment group
        group_count_tuples = sorted(group_count_dict.items(), key=lambda x: x[1])
        wide_count_df = count_df_dict[group_count_tuples[0][0]]. \
            merge(count_df_dict[group_count_tuples[1][0]], on=['x', 'y', 'channel', 'view'])
        wide_count_df['group_count_x'] = group_count_tuples[0][1]
        wide_count_df['group_count_y'] = group_count_tuples[1][1]
        wide_count_df['count_norm_ratio'] = wide_count_df.count_norm_x / wide_count_df.count_norm_y
        wide_count_df['count_norm_diff'] = wide_count_df.count_norm_x - wide_count_df.count_norm_y
        wide_count_df['count_norm_diff_abs'] = wide_count_df.count_norm_diff.abs()
        wide_count_df = wide_count_df.merge(
            wide_count_df.groupby(['channel', 'view']).count_norm_diff_abs.max().reset_index(name='max_diff')
        )
        wide_count_df['count_norm_diff_scaled'] = wide_count_df.count_norm_diff / wide_count_df.max_diff
        return wide_count_df

    if metadata.shape[0] > 0:
        if plot_type == 'total':
            plot_df = numpy_to_df(tuple(metadata.idx), sysmex_raw, pseudocount=pseudocount,
                                  gaussian_sigma=gaussian_sigma)
            plot_df = filter_plot_df(plot_df, channel_views)
            gplt = (ggplot(plot_df, aes('x', 'y', fill='count')) +
                    geom_tile() +
                    scale_fill_continuous(trans='log10') +
                    facet_wrap('~ channel') +
                    ylab('') +
                    xlab('') +
                    theme_minimal(base_size=8))

        elif plot_type == 'ratio' or plot_type == 'difference':
            group_count_dict = dict(metadata.value_counts(group))
            # keep only homozygous genotypes for now
            # could also try some sort of additive model
            if group == 'genotype':
                group_count_dict = {k: v for k, v in group_count_dict.items() if k in ['AA', 'CC', 'GG', 'TT']}
            for grp, grp_count in group_count_dict.items():
                metadata_ptb = metadata.loc[metadata[group] == grp]
                plot_df = numpy_to_df(tuple(metadata_ptb.idx), sysmex_raw, pseudocount=pseudocount,
                                      gaussian_sigma=gaussian_sigma)
                plot_df = filter_plot_df(plot_df, channel_views)
                plot_df['group'] = grp
                plot_df['group_count'] = grp_count
                count_df_dict[grp] = plot_df
            wide_count_df = wide_counts_merge(count_df_dict, group_count_dict)

        if group is not None:
            assert len(group_count_dict.keys()) == 2, 'Need exactly two groups for comparison plots.'
            if plot_type == 'ratio':
                legend_name = f'ratio\n' \
                              f'{wide_count_df.group_x[0]} ({wide_count_df.group_count_x[0]}) \n' \
                              f'\ {wide_count_df.group_y[0]} ({wide_count_df.group_count_y[0]})\n'
                gplt = (ggplot(wide_count_df, aes('x', 'y', fill='count_norm_ratio')) +
                        geom_tile() +
                        scale_fill_gradient2(trans='log10', name=legend_name) +
                        facet_wrap('~ channel') +
                        ylab('') +
                        xlab('') +
                        labs(title=title) +
                        theme_minimal(base_size=8, ) +
                        theme(legend_position=legend_position)
                        )

            elif plot_type == 'difference':
                legend_name = f'scaled difference\n' \
                              f'{wide_count_df.group_x[0]} ({wide_count_df.group_count_x[0]})\n' \
                              f'- {wide_count_df.group_y[0]} ({wide_count_df.group_count_y[0]})\n'
                gplt = (ggplot(wide_count_df, aes('x', 'y', fill='count_norm_diff_scaled')) +
                        geom_tile() +
                        scale_fill_gradient2(trans='pseudo_log', name=legend_name) +
                        facet_wrap('~ channel') +
                        ylab('') +
                        xlab('') +
                        labs(title=title) +
                        theme_minimal(base_size=8, ) +
                        theme(legend_position=legend_position)
                        )
            else:
                raise NotImplementedError

    return count_df_dict, plot_df, wide_count_df, gplt
