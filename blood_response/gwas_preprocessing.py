"""
GWAS preprocessing and phenotype QC.
"""

import pandas as pd
import matplotlib.pyplot as plt
import dill
from tqdm.auto import tqdm
from sklearn.preprocessing import quantile_transform

from sklearn.pipeline import make_pipeline
from sklearn.decomposition import PCA, FastICA, NMF
from sklearn.preprocessing import scale, robust_scale, minmax_scale, quantile_transform, power_transform
from scipy import stats
from joblib import delayed, dump

import umap
from .utils import ParallelExecutor
from .scatterplot import plot_scatter_2d


def drop_char(number_string):
    return ''.join(filter(str.isdigit, number_string))


def standardize_study_id(study_id, fixed_ids=None, drop_x=True):
    # standardize study ids to have this format:
    # PET-0001
    # JC-00001
    # FH-0001
    # X
    standardized_id = pd.NA

    if fixed_ids is not None and study_id in fixed_ids.index:
        study_id = fixed_ids.loc[study_id, 'standard_study_id']

    try:
        if study_id.startswith(('J', 'PC', 'W', 'w')) or (not drop_x and study_id.startswith('X')):
            clinic, number = study_id.split('-')
            standardized_id = f'{clinic.upper()}-{int(drop_char(number)):05}'
        elif study_id.startswith(('PET', 'AZ')):
            clinic, number = study_id.split('-')
            standardized_id = f'{clinic}-{int(drop_char(number)):04}'
        elif study_id.startswith('FH'):
            split_ids = study_id.split('-')
            # standardized_id = f'{clinic}-{int(number):04}-{family_id}'
            standardized_id = f'{split_ids[0]}-{int(drop_char(split_ids[1])):05}'
    except:
        pass

    return standardized_id


def pheno_quantile_transform(df, pheno_columns):
    return pd.DataFrame(
        quantile_transform(df[pheno_columns], output_distribution='normal', n_quantiles=5000),
        index=df.index,
        columns=[c + '-q' for c in pheno_columns])


def write_gcta_phenotypes(pheno_data, outpath, name):
    # Input file format
    # test.phen (no header line; columns are family ID, individual ID and phenotypes)
    for col in pheno_data.columns:
        if col not in ['#FID', 'IID']:
            pheno_data[['#FID', 'IID', col]].to_csv(f'{outpath}/{name}_{col}.phen',
                                                    sep='\t', header=False, index=False, na_rep='NA')


def prepare_phenotype_data(
        pheno_data,
        gwas_path,
        out_path,
        id_type='PatientID',
        covariate_cols=['BATCH', 'Sex', 'Age', 'Race'] + ['PC' + str(i+1) for i in range(10)],
        bedfile='imputed_merged',
        name=None,
        covariates=None,
        quantile_transform=True,
        plot_phenotype=True,
        table_phenotype=True
):
    """
    expects a df with patientid/mrn as index and plink coded phenotypes:
        categorical: 0 missing, 1 control, 2 case
        continuous: missing values are NaN
    """

    # filter any phenotypes that are constant
    # TODO dropping even if constant only in the white cohort
    def drop_constant_pheno(pheno_df, return_pheno_df=True):
        pheno_std = pheno_df.loc[pheno_df.Race == 'White', pheno_columns].std().reset_index(name='sd')
        constant_cols = []
        if any(pheno_std.sd == 0):
            constant_cols = list(pheno_std.loc[pheno_std.sd == 0, 'index'])
            print(f'Removing constant traits: {constant_cols}')
        if return_pheno_df:
            return pheno_df.copy(deep=True).drop(constant_cols, axis=1)
        else:
            return constant_cols

    if name is None:
        name = bedfile

    # make output dir for this data
    from pathlib import Path
    basepath = Path(out_path)
    if not basepath.exists():
        basepath.mkdir()

    # keep track of the trait columns
    pheno_columns = list(set(pheno_data.columns).difference(covariate_cols))

    # read in patient ids
    patient_ids = pd.read_csv(f'{gwas_path}/patient_ids.tsv', sep='\t', dtype=str)
    print(f'Subjects with GWAS data: \n'
          f'{len(set(pheno_data.index).intersection(patient_ids[id_type]))} out of {len(set(pheno_data.index))}')
    pheno_data = pheno_data.merge(patient_ids, left_index=True, right_on=id_type)

    # read in covariates including pca
    if covariates is None:
        covariates = pd.read_csv(f'{gwas_path}/covariate.txt', sep='\t', dtype=str)
    pca = pd.read_csv(f'{gwas_path}/{bedfile}_pca.eigenvec', sep='\t', dtype=str)

    # CHANGED: only add covariates that are not already in the pheno_data df
    pheno_data = pheno_data.merge(
        covariates[['IID'] + [c for c in covariates.columns if c not in pheno_data.columns]], on='IID'
    )
    pheno_data = pheno_data.merge(pca.drop('#FID', axis=1), on='IID')
    covariate_pca_header = ['#FID', 'IID'] + covariate_cols

    # make sure all traits are numeric
    if not all(pheno_data[pheno_columns] == 'float64'):
        print(f'WARNING: processing phenotype columns, converting all str cols to numeric:'
              f'{pheno_data[pheno_columns].dtypes}')
        pheno_data[pheno_columns] = pheno_data[pheno_columns].apply(pd.to_numeric)

    # remove any constant trait columns
    constant_cols = drop_constant_pheno(pheno_data, return_pheno_df=False)
    pheno_columns = [c for c in pheno_columns if c not in constant_cols]
    pheno_data = pheno_data.drop(constant_cols, axis=1)

    # add quantile transformed version of traits
    if quantile_transform:
        quantile_data = pheno_quantile_transform(pheno_data[pheno_columns], pheno_columns)
        pheno_data = pheno_data.merge(quantile_data, left_index=True, right_index=True)
        pheno_columns = pheno_columns + list(quantile_data.columns)

    # drop any duplicated subjects
    if any(pheno_data.PatientID.duplicated()):
        print(f'Warning: removing {pheno_data.PatientID.duplicated().sum()} duplicated PatientIDs.')
        pheno_data = pheno_data.drop_duplicates(subset='PatientID')
    if any(pheno_data.IID.duplicated()):
        print(f'Warning: removing {pheno_data.IID.duplicated().sum()} duplicated IIDs.')
        pheno_data = pheno_data.drop_duplicates(subset='IID')

    # drop anyone with missing covariates
    if any(pheno_data[covariate_pca_header].isna()):
        print(f'Warning: removing subjects with NA entries in covariates: \n '
              f'{pheno_data[covariate_pca_header].isna().sum(axis=0)}')
        pheno_data = pheno_data.dropna(subset=covariate_pca_header, axis=0)

    # change categorical batch covariates to start with string
    if 'BATCH' in covariate_cols:
        pheno_data['BATCH'] = 'b' + pheno_data['BATCH'].astype(str)

    # write out the covariates
    pheno_data[covariate_pca_header]. \
        to_csv(f'{out_path}/{name}_covariates.tsv', index=False, header=True, sep='\t')

    # write out the covariates for gcta
    # quantitative covariates
    qcov_cols = sorted(list(set(covariate_cols).intersection(['Age', 'bmi', 'hr'] + ['PC' + str(i+1) for i in range(10)])))
    pheno_data[['#FID', 'IID'] + qcov_cols]. \
        to_csv(f'{out_path}/{name}_qcov.tsv', index=False, header=False, sep='\t')

    # categorical covariates
    ccov_cols = sorted(list(set(covariate_cols).intersection(['BATCH', 'Sex', 'Race', 'Model'])))
    pheno_data[['#FID', 'IID'] + ccov_cols]. \
        to_csv(f'{out_path}/{name}_ccov.tsv', index=False, header=False, sep='\t')

    pheno_data.loc[pheno_data.Race == 'White'][covariate_pca_header]. \
        to_csv(f'{out_path}/{name}_covariates_white.tsv', index=False, header=True, sep='\t')

    # write out the phenotypes
    pheno_data[['#FID', 'IID'] + pheno_columns]. \
        to_csv(f'{out_path}/{name}_pheno.tsv', index=False, header=True, na_rep='NA', sep='\t')

    pheno_data.loc[pheno_data.Race == 'White'][['#FID', 'IID'] + pheno_columns]. \
        to_csv(f'{out_path}/{name}_pheno_white.tsv', index=False, header=True, na_rep='NA', sep='\t')

    # write_gcta_phenotypes(pheno_data[['#FID', 'IID'] + pheno_columns], out_path, name=name)
    # write_gcta_phenotypes(pheno_data.loc[pheno_data.Race == 'White'][['#FID', 'IID'] + pheno_columns],
    #                       out_path, name=name + '_white')

    pheno_data.to_csv(f'{out_path}/{name}_pheno_cov_data.tsv', index=False, header=True, na_rep='NA', sep='\t')

    # Make histogram that shows distribution of numeric phenotypes:
    if plot_phenotype:
        plot_pheno(pheno_data, pheno_columns, f'{out_path}/{name}_phenotype_histogram.pdf')
        # num_pheno = int(len(pheno_columns))
        # if num_pheno % 2 == 0:
        #     #fig_pheno, ax = plt.subplots(int(num_pheno / 2), 2, figsize=(10, num_pheno * 2.5))
        #     fig_pheno, ax = plt.subplots(2, int(num_pheno / 2), figsize=(num_pheno * 2, 8))
        # else:
        #     fig_pheno, ax = plt.subplots(num_pheno, 1, figsize=(6, num_pheno * 3))
        # pheno_data[pheno_columns].hist(ax=ax)
        # fig_pheno.savefig(f'{out_path}/{name}_phenotype_histogram.pdf')
        # # close the figure so it doesn't get displayed in interactive sessions
        # plt.close(fig_pheno)

    # Make table that contains basic statistics of numeric phenotypes:
    if table_phenotype:
        tb_pheno = pheno_data[pheno_columns].describe()
        tb_pheno.to_csv(f'{out_path}/{name}_pheno_stat.tsv', index=True, header=True, na_rep='NA', sep='\t')


def plot_pheno(pheno_data, pheno_columns, filename):
    # Make histogram that shows distribution of numeric phenotypes:
    num_pheno = int(len(pheno_columns))
    if num_pheno % 2 == 0:
        #fig_pheno, ax = plt.subplots(int(num_pheno / 2), 2, figsize=(10, num_pheno * 2.5))
        fig_pheno, ax = plt.subplots(2, int(num_pheno / 2), figsize=(num_pheno * 2, 8))
    else:
        fig_pheno, ax = plt.subplots(num_pheno, 1, figsize=(6, num_pheno * 3))
    pheno_data[pheno_columns].hist(ax=ax)
    fig_pheno.savefig(filename)
    # close the figure so it doesn't get displayed in interactive sessions
    plt.close(fig_pheno)


def filter_outliers_mad(df, feature, outlier_cutoff=4, return_stats=False):
    median = df[feature].median()
    estimated_sd = stats.median_abs_deviation(df[feature], scale="normal", nan_policy='omit')
    upper = median + outlier_cutoff * estimated_sd
    lower = median - outlier_cutoff * estimated_sd
    if return_stats:
        return median, estimated_sd, upper, lower
    else:
        return (df[feature] <= upper) & (df[feature] >= lower)


# added 21/02/26: plotting for manually gated data
def project_phenotypes(
        pheno_df,
        pheno_cols,
        output_dir,
        impute=False,
        load_embeddings=False,
        reducers=['pcaj', 'pcas', 'picas', 'picaj', 'umapj', 'umaps'],
        n_components=[2],
        outlier_prefix='picas',
        outlier_sd_cutoff=4,
        make_plots=True,
        plot_outliers=True,
        n_jobs=20,
        random_state=42,
):

    """
    Using for GWAS preprocessing like this:

        # optionally impute values if only a few are missing
        imputer = KNNImputer(n_neighbors=5, weights="uniform")
        transformed = imputer.fit_transform(pheno_df[pheno_cols])
        pheno_df[list(pheno_cols)] = transformed

        # project with PCA, ICA and UMAP, flag outliers and make plots with various covariates
        pheno_df_projected = project_phenotypes(
        pheno_df,
        pheno_cols,
        output_dir=analysis_path + f'/projections/{channel}/',
        impute=False, load_embeddings=False,
        make_plots=True, plot_outliers=True
    )
    """

    # make output dir for this data
    from pathlib import Path
    basepath = Path(f'{output_dir}')
    if not basepath.exists():
        basepath.mkdir(parents=True, exist_ok=True)

    if impute:
        from sklearn.impute import KNNImputer
        imputer = KNNImputer(n_neighbors=5, weights="uniform")
        pheno_df[pheno_cols] = imputer.fit_transform(pheno_df[pheno_cols])

    if 'DATE' in pheno_df:
        pheno_df['analysis_date'] = pd.to_datetime(pheno_df.DATE).dt.strftime("%Y%m")
        pheno_df['analysis_month'] = (pd.to_datetime(pheno_df.DATE) - pd.to_datetime('2019')).dt.days / 30  # month since 2019
        pheno_df['analysis_time'] = pd.to_numeric(pd.to_datetime(pheno_df.ETIM).dt.strftime("%H%M"))

    reducers_ = {
        'umapj': umap.UMAP,
        'umaps': umap.UMAP,
        'pcaj': lambda *args, **kwargs: PCA(*args, **kwargs, svd_solver='randomized'),
        'pcas': lambda *args, **kwargs: PCA(*args, **kwargs, svd_solver='randomized'),
        # CHANGED doing ICA directly instead of PCA to 100 components and then ICA
        # 'picaj': lambda *args, **kwargs: make_pipeline(PCA(n_components=100, random_state=random_state),
        #                                                FastICA(*args, **kwargs)),
        # 'picas': lambda *args, **kwargs: make_pipeline(PCA(n_components=100, random_state=random_state),
        #                                                FastICA(*args, **kwargs)),
        'picaj': lambda *args, **kwargs: FastICA(*args, **kwargs),
        'picas': lambda *args, **kwargs: FastICA(*args, **kwargs),
        'nmfj': lambda *args, **kwargs: NMF(*args, **kwargs),
        'nmfs': lambda *args, **kwargs: NMF(*args, **kwargs),
    }

    parallel_executor = ParallelExecutor(n_jobs=n_jobs)

    # TODO ignoring any arguments that are being passed to wrapper
    def reduce_separately(reducer, latent_df):
        # define reducer task
        def reduce_(ptb):
            latent_ptb = latent_df.loc[latent_df.Perturbation == ptb]
            embedding = reducer.fit_transform(latent_ptb[pheno_cols])
            return pd.DataFrame(embedding, index=latent_ptb.index)

        # process each perturbation separately
        perturbations = latent_df.Perturbation.unique()
        embeddings = parallel_executor(total=len(perturbations))(delayed(reduce_)(ptb) for ptb in perturbations)
        embeddings = pd.concat(embeddings)
        # match the order of the latent df and return a plain numpy array
        return embeddings.loc[latent_df.index].values

    lat_dims = []
    for reducer_name in tqdm(reducers):
        for n_comp in n_components:
            if load_embeddings:
                with open(basepath / f'{reducer_name}_embedding_{n_comp}.pic', 'rb') as f:
                    embedding = dill.load(f)
            else:
                reducer = reducers_[reducer_name](random_state=random_state, n_components=n_comp)
                if reducer_name.endswith('j'):
                    embedding = reducer.fit_transform(pheno_df[pheno_cols])
                elif reducer_name.endswith('s'):
                    embedding = reduce_separately(reducer, pheno_df)
                with open(basepath / f'{reducer_name}_embedding_{n_comp}.pic', 'wb') as f:
                    dill.dump(embedding, f)

            for col in range(embedding.shape[1]):
                lat_dim = f'{reducer_name}{n_comp}{col}'
                lat_dims.append(lat_dim)
                pheno_df[lat_dim] = embedding[:, col]

    pheno_df = pheno_df.rename(columns={'sample_id': 'sample'})

    # flag global outliers based on the joint PCA projection
    pheno_df['outlier'] = ~(
            filter_outliers_mad(pheno_df, f'{outlier_prefix}20', outlier_cutoff=outlier_sd_cutoff) &
            filter_outliers_mad(pheno_df, f'{outlier_prefix}21', outlier_cutoff=outlier_sd_cutoff)
    )

    if make_plots:
        if plot_outliers:
            (basepath / 'with_outliers/').mkdir(exist_ok=True)
            parallel_executor(total=len(reducers))(
                delayed(plot_scatter_2d)(
                    meta=pheno_df[['PatientID', 'sample', 'Perturbation', 'outlier', f'{reducer}20', f'{reducer}21']],
                    basepath=basepath / 'with_outliers/', prefix=reducer
                ) for reducer in reducers
            )

        parallel_executor(total=len(reducers))(
            delayed(plot_scatter_2d)(
                meta=pheno_df.loc[~pheno_df.outlier], basepath=basepath, prefix=reducer
            ) for reducer in reducers
        )

    return pheno_df
