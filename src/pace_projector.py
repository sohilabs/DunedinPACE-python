import pandas as pd
import numpy as np
from rpy2.robjects.packages import importr
import rpy2.robjects as ro

preprocessCore = importr("preprocessCore")
base = importr("base")


def pace_projector(
    betas, mPACE_Models, proportion_of_probes_required=0.8
) -> pd.DataFrame:
    """
    Calculate PACE scores for a set of samples using a set of mPACE models
    """
    # Loop through models
    model_results = {}
    for model_name in mPACE_Models["model_names"]:
        # Make sure it has been converted to a DataFrame
        if not isinstance(betas, pd.DataFrame) or not all(
            betas.applymap(np.isreal).all()
        ):
            raise ValueError("betas DataFrame is not numeric!")

        probe_overlap = len(
            set(betas.index).intersection(mPACE_Models["model_probes"][model_name])
        ) / len(mPACE_Models["model_probes"][model_name])
        probe_overlap_background = len(
            set(betas.index).intersection(
                mPACE_Models["gold_standard_probes"][model_name]
            )
        ) / len(mPACE_Models["gold_standard_probes"][model_name])

        # Make sure enough of the probes are present in the data file
        if (
            probe_overlap < proportion_of_probes_required
            or probe_overlap_background < proportion_of_probes_required
        ):
            result = pd.Series([np.nan] * betas.shape[1], index=betas.columns)
        else:
            # Work with a numeric matrix of betas
            betas_mat = betas.loc[
                betas.index.intersection(
                    mPACE_Models["gold_standard_probes"][model_name]
                ),
                :,
            ].copy()

            # If probes don't exist, we'll add them as rows of values based on their mean in the gold standard dataset
            probes_not_in_matrix = set(
                mPACE_Models["gold_standard_probes"][model_name]
            ).difference(betas_mat.index)
            for probe in probes_not_in_matrix:
                tmp_mat = pd.DataFrame(0, index=[probe], columns=betas_mat.columns)
                tmp_mat.loc[probe, :] = mPACE_Models["gold_standard_means"][model_name][
                    probe
                ]
                betas_mat = pd.concat([betas_mat, tmp_mat])

            # Identify samples with too many missing probes and remove them from the matrix
            samples_to_remove = betas_mat.columns[
                (betas_mat.count() / betas_mat.shape[0]) < proportion_of_probes_required
            ]
            if len(samples_to_remove) > 0:
                betas_mat.drop(columns=samples_to_remove, inplace=True)

            if betas_mat.shape[1] > 0:
                # Identify missingness on a probe level
                pct_values_present = betas_mat.count() / betas_mat.shape[1]

                # If they're missing values, but less than the proportion required, we impute to the cohort mean
                probes_to_adjust = pct_values_present.index[
                    (pct_values_present < 1)
                    & (pct_values_present >= proportion_of_probes_required)
                ]
                betas_mat.loc[probes_to_adjust, :] = betas_mat.loc[
                    probes_to_adjust, :
                ].apply(lambda x: x.fillna(x.mean()), axis=1)

                # If they're missing too many values, everyone's value gets replaced with the mean from the Dunedin cohort
                probes_to_replace_with_mean = pct_values_present.index[
                    pct_values_present < proportion_of_probes_required
                ]
                for probe in probes_to_replace_with_mean:
                    betas_mat.loc[probe, :] = mPACE_Models["model_means"][model_name][
                        probe
                    ]

                # Normalize the matrix to the gold standard dataset
                # We convert betas and gold_standard_means to R matrices and vectors because that's what the R function expects
                betas_norm = preprocessCore.normalize_quantiles_use_target(
                    base.as_matrix(ro.conversion.py2rpy(betas_mat)),
                    ro.vectors.FloatVector(
                        mPACE_Models["gold_standard_means"][model_name]
                    ),
                )  # Assuming normalize_quantiles_use_target is defined elsewhere
                betas_norm = pd.DataFrame(betas_norm)
                betas_norm.index = betas_mat.index
                betas_norm.columns = betas_mat.columns

                # Calculate score:
                score = mPACE_Models["model_intercept"][model_name] + np.sum(
                    betas_norm.loc[mPACE_Models["model_probes"][model_name], :].T
                    * mPACE_Models["model_weights"][model_name],
                    axis=1,
                )
                score.index = betas_norm.columns

                if len(samples_to_remove) > 0:
                    score_tmp = pd.Series(
                        [np.nan] * len(samples_to_remove), index=samples_to_remove
                    )
                    score = pd.concat([score, score_tmp])
                score = score.loc[betas.columns]
                result = score
            else:
                result = pd.Series(
                    [np.nan] * betas_mat.shape[1], index=betas_mat.columns
                )
        model_results[model_name] = result
    return pd.DataFrame(model_results)
