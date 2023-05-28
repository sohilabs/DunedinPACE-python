import rdata
from rpy2.robjects import pandas2ri
import pandas as pd
from pace_projector import pace_projector

# Activate automatic conversion of R objects into pandas data frames
pandas2ri.activate()

# Parse and convert the file containing the beta values
beta_matrix = rdata.conversion.convert(
    rdata.parser.parse_file("data/example_betas.rda")
)["example_betas"]
# Convert the beta values to a pandas DataFrame
betas_df = pd.DataFrame(
    beta_matrix.assays.data.listData["counts"],
    index=beta_matrix.NAMES,
    columns=beta_matrix.colData.rownames,
)

# Parse and convert the file containing the mPACE models
mPACE_Models = rdata.conversion.convert(rdata.parser.parse_file("data/sysdata.rda"))[
    "mPACE_Models"
]

# Print the first 5 rows of the betas_df DataFrame
print(f"Betas Matrix samples:\n{betas_df.head()}\n")

# Print the first 5 values of the mPACE_Models regression weights
print(
    f"mPACE_Model example regression weights:\n{list(mPACE_Models['model_weights'].values())[0][:5]}\n"
)

# Compute pace_values
print("Computing pace_values...\n")
pace_values = pace_projector(betas_df, mPACE_Models)

# Print pace_values
print(f"Computed Pace values:\n{pace_values}\n")

# Validate matches R output
R_output = "7786915023_R02C02 7786915135_R04C02 7471147149_R06C01 7786915035_R05C01 7786923035_R01C01\n1.0972737         1.1191063         0.9958676         1.3322271     1.1275386"
r_output_df = pd.DataFrame(
    list(map(float, R_output.split("\n")[1].split())),
    index=R_output.split("\n")[0].split(),
    columns=["DunedinPACE"],
)
print(f"R output to match:\n{r_output_df}\n")
# Compare for each row value, the abs difference between the R and Python output is less than 1e-3
print(
    f"Absolute diff tolerance below 1e-3:\n{(pace_values - r_output_df).abs() < 1e-3}"
)
