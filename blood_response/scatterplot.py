"""
Scatterplots for phenotype QC.
"""

from plotnine import *
import matplotlib.pyplot as plt


def plot_scatter(x, y, c, name, marker=None, size=10, alpha=0.5):
    plt.scatter(x, y, c=c, marker=marker, cmap='Spectral', alpha=alpha, s=size)
    plt.gca().set_aspect('equal', 'datalim')
    plt.colorbar()
    plt.title(name, fontsize=12)
    # plt.show()


def plot_scatter_2d(
        meta,
        prefix,
        basepath,
        plot_vars_discrete=['CYT', 'factor(keep_pred)', 'outlier',
                            'Sex', 'Race', 'Ethnicity', 'Tobacco',
                            'Transplant', 'Glucocorticoid', 'Transplant'],
        plot_vars_cont=['Perturbation', 'perturbation_str',
                        'analysis_date', 'analysis_time',
                        'true_score', 'Age', 'BMI']
):
    highlight_samples = meta. \
        drop_duplicates(subset=['PatientID', 'Perturbation']). \
        groupby('Perturbation'). \
        tail(n=1)['sample']

    # # TODO continuous vs discrete scale
    # plot_vars_discrete = ['CYT', 'factor(keep_pred)', 'Sex', 'Race', 'Ethnicity', 'Tobacco', 'Transplant']
    # # too many classes in perturbation_str - so using the cont plotting version
    # plot_vars_cont = ['' 'perturbation_str', 'analysis_date', 'analysis_time', 'true_score', 'Age']

    for plot_var in set(plot_vars_discrete).intersection(meta.columns):
        (ggplot(meta, aes(f'{prefix}20', f'{prefix}21', color=plot_var)) +
         geom_point(alpha=0.05) +
         geom_rug(alpha=0.25) +
         geom_text(aes(label='Perturbation'), data=meta.loc[meta['sample'].isin(highlight_samples)], size=10,
                   color='black',
                   adjust_text=dict(arrowprops=dict(arrowstyle="-", color='k', lw=0.5))) +
         theme_minimal() +
         scale_color_brewer(type='qual', palette='Set1')). \
            save(f'{basepath}/scatter2d_{prefix}_{plot_var}.pdf', width=10, height=10)

    for plot_var in set(plot_vars_cont).intersection(meta.columns):
        (ggplot(meta, aes(f'{prefix}20', f'{prefix}21', color=plot_var)) +
         geom_point(alpha=0.05) +
         geom_rug(alpha=0.25) +
         geom_text(aes(label='Perturbation'), data=meta.loc[meta['sample'].isin(highlight_samples)], size=10,
                   color='black',
                   adjust_text=dict(arrowprops=dict(arrowstyle="-", color='k', lw=0.5))) +
         theme_minimal()
         ). \
            save(f'{basepath}/scatter2d_{prefix}_{plot_var}.pdf', width=10, height=10)


    # faceted plots by perturbation
    for plot_var in set(plot_vars_discrete).intersection(meta.columns):
        (ggplot(meta, aes(f'{prefix}20', f'{prefix}21', color=plot_var)) +
         geom_point(alpha=0.15) +
         facet_wrap('~Perturbation') +
         geom_rug(alpha=0.25) +
         theme_minimal() +
         scale_color_brewer(type='qual', palette='Set1')). \
            save(f'{basepath}/scatter2d_{prefix}_facet_{plot_var}.pdf', width=20, height=20)

    for plot_var in set(plot_vars_cont).intersection(meta.columns):
        (ggplot(meta, aes(f'{prefix}20', f'{prefix}21', color=plot_var)) +
         geom_point(alpha=0.15) +
         facet_wrap('~Perturbation') +
         geom_rug(alpha=0.25) +
         theme_minimal()). \
            save(f'{basepath}/scatter2d_{prefix}_facet_{plot_var}.pdf', width=20, height=20)
