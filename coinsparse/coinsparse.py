import click
import fnmatch
import numpy as np
import os
import os.path as op
import pandas as pd
import re
import warnings

from glob import glob
from tqdm import tqdm


def _get_coins_path(coins_basename, coins_dir):
    coins_pattern = "_".join(
        [
            "[0-9][0-9][0-9][0-9]",
            coins_basename,
            "[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9].csv",
        ]
    )
    coins_regex = re.compile(fnmatch.translate(coins_pattern), re.IGNORECASE)
    return [
        op.join(coins_dir, fn) for fn in os.listdir(coins_dir) if coins_regex.match(fn)
    ][0]


def parse_single_xlsx(
    data_dict_filename,
    coins_dir,
    derivative_subjects=None,
    derivative_name=None,
    count_subs=True,
):
    """
    Parameters
    ----------
    data_dict_filename : str
        Filename for an individual data dictionary excel file

    coins_dir : str
        Directory containing all of the COINS phenotypical data

    derivative_subjects : set, optional
        Set containing subject identifiers of subjects with derivative data

    derivative_name : str, optional
        The name of the derivative data with which to intersect the phenotypical
        responses

    count_subs : bool, default=True
        If True, count the number of non-null responses for each variable and,
        optionally, the number of those subjects with derivative data.
    """
    assessment = pd.read_excel(data_dict_filename, engine="openpyxl", header=0).columns[
        0
    ]
    header = 2 if "ESWAN" in data_dict_filename else 1

    df = pd.read_excel(data_dict_filename, engine="openpyxl", header=header).dropna(
        axis="columns", how="all"
    )

    df.rename(
        columns={
            "Question ": "Question",
            "Variable": "Variable Name",
            "Value": "Values",
            "Value Label": "Value Labels",
            "Subtest": "Question",
            "Item": "Question",
            "Scores": "Question",
        },
        inplace=True,
    )

    df.dropna(how="all", subset=df.columns.drop("Question").tolist(), inplace=True)

    if "SWAN" in data_dict_filename:
        df["Value Labels"] = df["Value Labels"].values[0]

    df.drop(
        columns=[
            "Run 1",
            "Run 2",
            "Run 3",
            "Run 4",
            "Run 4 ",
            "Run 5",
            "Run 6",
            "Description",
            "Value Scale",
            "5-Trial Adjusting Delay (Temporal Discounting) Task Protocol",
            "Unnamed: 1",
            "Unnamed: 2",
            "Unnamed: 3",
            "Unnamed: 4",
            "Unnamed: 5",
            "Unnamed: 6",
            "Unnamed: 7",
            "Response Values",
        ],
        inplace=True,
        errors="ignore",
    )

    df["Assessment"] = assessment

    coins_filename = (
        op.splitext(op.basename(data_dict_filename))[0]
        .replace("DigitSpan", "Digit_Span")
        .replace("CDI", "CDI2")
        .replace("CTOPP", "CTOPP_2")
        .replace("KSADS_C", "KSADS_CMI")
        .replace("TRF_P", "TRF_Pre")
        .replace("ACE_Spatial_Span", "ACE_Spatial")
    )

    df["COINS Filename"] = coins_filename

    if coins_filename == "CELF":
        df["COINS Filename"] = "CELF_5_Screen"
    elif coins_filename == "CELF_Meta":
        df["COINS Filename"] = "CELF5_Meta"
    elif coins_filename == "Diagnosis_ClinicianConsensus":
        df["COINS Filename"] = "ConsensusDx"
    elif coins_filename == "Diagnosis_KSADS":
        df["COINS Filename"] = "KSADS"
    elif coins_filename == "FGA":
        df["COINS Filename"] = "FitnessGram_Adult"
    elif coins_filename == "FGC":
        df["COINS Filename"] = "FitnessGram"
    elif coins_filename == "PhenX_School":
        df["COINS Filename"] = "PhenX_SchoolRisk"
    elif coins_filename == "WAIS_abb":
        df["COINS Filename"] = "WAIS_Abbreviated"
    elif coins_filename in ["KSADS_P", "NIDA", "PreInt_FamHx", "Quotient_Ext"]:
        df["COINS Filename"] = np.nan
        if count_subs:
            df["non-na count"] = np.nan
            if derivative_subjects is not None:
                df["with " + derivative_name] = np.nan
        return df

    coins_filename = pd.unique(df["COINS Filename"])[0]
    coins_path = _get_coins_path(coins_filename, coins_dir)

    if count_subs:
        try:
            response_df = (
                pd.read_csv(coins_path, low_memory=False)
                .drop(0)
                .drop_duplicates(subset="EID", ignore_index=True)
            )
        except KeyError:
            response_df = (
                pd.read_csv(coins_path, low_memory=False)
                .drop(0)
                .drop_duplicates(subset="EID #1", ignore_index=True)
            )

        non_na_count = len(response_df) - pd.DataFrame(response_df.isna().sum())
        non_na_count.columns = ["non-na count"]

        if derivative_subjects is not None:
            num_derivative_subs = [
                len(
                    derivative_subjects.intersection(
                        set(response_df.dropna(subset=[var])["EID"])
                    )
                )
                if var in response_df.columns
                else np.nan
                for var in df["Variable Name"]
            ]
            df["with " + derivative_name] = num_derivative_subs

        df = df.merge(
            non_na_count,
            how="left",
            left_on="Variable Name",
            right_index=True,
        )

    return df


def create_hbn_summary(
    data_dict_dir,
    coins_dir,
    assessments=None,
    derivatives_file=None,
    derivative_name=None,
    count_subs=True,
):
    """
    Parameters
    ----------
    data_dict_dir : str
        Directory containing the HBN data dictionaries

    coins_dir : str
        Directory containing the COINS data

    assessments : str, optional
        If provided, limit the output to this assessments only. This string must
        match the name of the data dictionary file.

    derivatives_file : str, optional
        Filename containing a derivatives dataset with a "subjectID" column

    derivative_name : str, optional
        The name of the derivative data with which to intersect the phenotypical responses

    count_subs : bool, default=True
        If True, count the number of non-null responses for each variable and,
        optionally, the number of those subjects with derivative data.
    """
    data_dict_files = glob(op.join(data_dict_dir, "*.xlsx"))
    data_dict_files = [fn for fn in data_dict_files if "~$" not in fn and " " not in fn]

    if assessments is not None:
        if isinstance(assessments, str):
            assessment_list = [assessments]
        elif all(isinstance(el, str) for el in assessments):
            assessment_list = list(assessments)
        else:
            raise TypeError("assessments must be a string or a sequence of strings.")

        data_dict_files = [
            fn
            for fn in data_dict_files
            if op.splitext(op.basename(fn))[0] in assessment_list
        ]

        if not len(data_dict_files) > 0:
            raise ValueError(
                "The assessment(s) that you provided does not exist in the "
            )

    deriv_name = derivative_name

    if derivatives_file is not None:
        derivatives_subs = set(
            pd.unique(pd.read_csv(derivatives_file, usecols=["subjectID"])["subjectID"])
        )
        if derivative_name is None:
            deriv_name = "derivative"
    else:
        derivatives_subs = None
        deriv_name = None

    df_summary = pd.concat(
        [
            parse_single_xlsx(
                data_dict_filename=fn,
                coins_dir=coins_dir,
                derivative_subjects=derivatives_subs,
                derivative_name=deriv_name,
                count_subs=count_subs,
            )
            for fn in tqdm(data_dict_files)
        ]
    ).reset_index(drop=True)

    df_summary.sort_values(
        by="COINS Filename",
        inplace=True,
    )

    df_summary.reset_index(drop=True, inplace=True)

    return df_summary


def _get_variable_info(variable, data_dict_dir, coins_dir, return_df=False):
    # Create a dataframe with all assessments and variable names
    # Don't bother to count non-na subjects at this point
    with warnings.catch_warnings():
        # this will suppress all warnings in this block
        warnings.simplefilter("ignore", category=UserWarning)
        df_data_dict = create_hbn_summary(
            data_dict_dir=data_dict_dir, coins_dir=coins_dir, count_subs=False
        )

    # Confirm that the requested variables actually exist in COINS
    allowed_vars = df_data_dict["Variable Name"].to_list()
    for v in variable:
        if not v in allowed_vars:
            raise click.ClickException(
                f"{v} is not a valid COINS variable. Please consult the data "
                "dictionary to find a valid variable name."
            )

    # Keep only the variables requested by the user
    var_mask = df_data_dict["Variable Name"].isin(variable)
    df_data_dict = df_data_dict[var_mask]

    # We will store sets of subject IDs in a list
    # Each element will be a set of subject IDs that are present for a variable
    # or derivative. At the end, we will take the intersection of all of the sets.
    subject_sets = []

    # If return_df is True, we will return a dataframe with all of the merged
    # phenotypic data. Store each assessment's dataframe in a list and merge later
    var_dataframes = []

    # Iterate through the COINS files
    for coins_file in pd.unique(df_data_dict["COINS Filename"]):
        coins_path = _get_coins_path(coins_file, coins_dir)
        subject_col = "EID"
        try:
            response_df = (
                pd.read_csv(coins_path, low_memory=False)
                .drop(0)
                .drop_duplicates(subset="EID", ignore_index=True)
            )
        except KeyError:
            response_df = (
                pd.read_csv(coins_path, low_memory=False)
                .drop(0)
                .drop_duplicates(subset="EID #1", ignore_index=True)
            )
            subject_col = "EID #1"

        df_assessment = df_data_dict[df_data_dict["COINS Filename"] == coins_file]

        if return_df:
            var_columns = [subject_col] + [
                var for var in df_assessment["Variable Name"]
            ]
            var_dataframes.append(
                response_df[var_columns].set_index(subject_col, drop=True)
            )

        # Now iterate through each requested variable in this assessment
        for var in df_assessment["Variable Name"]:
            # Drop na responses and collect the subject IDs
            subject_sets.append(
                set(pd.unique(response_df.dropna(subset=[var])[subject_col]))
            )

    subjects = set.intersection(*subject_sets)
    if return_df:
        return subjects, pd.concat(var_dataframes, axis="columns", join="inner")
    else:
        return subjects


@click.group()
def cli():
    """Summarize phenotypic data from the Healthy Brain Network

    This software assumes that the data dictionaries and phenotypic files are in the
    format provided by the COINS data exchange.

    If you don't currently have access to HBN phenotypic data, follow the instructions at
    http://fcon_1000.projects.nitrc.org/indi/cmi_healthy_brain_network/Pheno_Access.html
    to gain access.
    """
    pass


@cli.command()
@click.option("--output", type=str, default="output.csv", help="Output file path.")
@click.option(
    "-a",
    "--assessment",
    multiple=True,
    default=[],
    help=(
        "Assessments with which to limit the output. These assessments must "
        "correspond to filenames (without the xlsx extension) in the data "
        "dictionary directory. To specify multiple assessments, provide the "
        "option multiple times. e.g '-a FSQ -a WIAT'."
    ),
)
@click.option(
    "--deriv-file",
    help=(
        "Path to csv file containing derivatives subjects. If supplied, "
        "this must contain a column titled 'subjectID'. If supplied, the output "
        "csv will compute the number of subjects in the intersection between the "
        "derivatives file and each assessment."
    ),
    type=str,
    default=None,
)
@click.option(
    "--deriv-name",
    help=(
        "Name of the derivative represented by the deriv-file argument. If "
        "not supplied, we will attempt to infer it from the filename."
    ),
    type=str,
    default=None,
)
@click.argument("data_dict_dir", type=str)
@click.argument("coins_dir", type=str)
def summarize(output, assessment, deriv_file, deriv_name, data_dict_dir, coins_dir):
    """Create a csv file containing summary information on HBN phenotypical assessments.

    This software assumes that the data dictionaries and phenotypic files are in the
    format provided by the COINS data exchange.

    If you don't currently have access to HBN phenotypic data, follow the instructions at
    http://fcon_1000.projects.nitrc.org/indi/cmi_healthy_brain_network/Pheno_Access.html
    to gain access.

    \b
    Arguments:\b
      DATA_DICT_DIR       Path to the data dictionary directory\b
      COINS_DIR           Path to the COINS data directory\b
    """
    if not assessment:
        assessment_arg = None
    else:
        assessment_arg = assessment

    with warnings.catch_warnings():
        # this will suppress all warnings in this block
        warnings.simplefilter("ignore", category=UserWarning)
        df = create_hbn_summary(
            data_dict_dir=data_dict_dir,
            coins_dir=coins_dir,
            assessments=assessment_arg,
            derivatives_file=deriv_file,
            derivative_name=deriv_name,
        )

    df.to_csv(output)
    click.echo("Output saved to " + output)


@cli.command("list-subjects")
@click.option(
    "-d",
    "--deriv-file",
    multiple=True,
    default=[],
    help=(
        "Path to csv file containing derivatives subjects. If supplied, "
        "this must contain a column titled 'subjectID'. If supplied, the output "
        "subject list will contain subjects with BOTH assessment responses and "
        "derivative data. To specify multiple derivatives, provide the option"
        "multiple times. e.g. '-d file1.csv -d file2.csv'."
    ),
)
@click.option(
    "-v",
    "--variable",
    multiple=True,
    default=[],
    help=(
        "Return subjects who have non-null responses to these variables. These variable names must "
        "correspond to variable names in the data dictionaries."
        "To specify multiple variables, provide the "
        "option multiple times. e.g '-v FSQ_08 -v WIAT_RC_Stnd'."
    ),
)
@click.option(
    "-o",
    "--output",
    type=click.File("w"),
    default="-",
    help="Redirect output to this file.",
)
@click.argument("data_dict_dir", type=str)
@click.argument("coins_dir", type=str)
def list_subjects(deriv_file, variable, output, data_dict_dir, coins_dir):
    """Return the list of subjects who have BOTH valid responses to certain variables AND derivative data.

    This software assumes that the data dictionaries and phenotypic files are in the
    format provided by the COINS data exchange.

    If you don't currently have access to HBN phenotypic data, follow the instructions at
    http://fcon_1000.projects.nitrc.org/indi/cmi_healthy_brain_network/Pheno_Access.html
    to gain access.

    \b
    Arguments:\b
      DATA_DICT_DIR       Path to the data dictionary directory\b
      COINS_DIR           Path to the COINS data directory\b
    """
    if not (deriv_file or variable):
        raise click.ClickException(
            "You must specify either derivatives (using the -d option) or "
            "variables (using the -v option)."
        )

    if variable:
        pheno_subs = _get_variable_info(
            variable, data_dict_dir, coins_dir, return_df=False
        )

    # We will store sets of subject IDs in a list
    # Each element will be a set of subject IDs that are present for a derivative
    # At the end, we will take the intersection of all of the sets.
    deriv_subs = []
    for dfile in deriv_file:
        deriv_subs.append(
            set(pd.unique(pd.read_csv(dfile, usecols=["subjectID"])["subjectID"]))
        )

    if deriv_subs:
        deriv_subs = set.intersection(*deriv_subs)
        if variable:
            subjects = list(pheno_subs.intersection(deriv_subs))
        else:
            subjects = list(deriv_subs)
    else:
        subjects = list(pheno_subs)

    click.echo("\n".join(subjects), file=output)
    return subjects


@cli.command()
@click.option(
    "--outdir", type=str, default="merged_data", help="Output directory path."
)
@click.option("-s", "--suffix", type=str, default="_merged", help="Output file path.")
@click.option(
    "-d",
    "--deriv-file",
    multiple=True,
    default=[],
    help=(
        "Path to csv file containing derivatives subjects. If supplied, "
        "this must contain a column titled 'subjectID'. If supplied, the output "
        "subject list will contain subjects with BOTH assessment responses and "
        "derivative data. To specify multiple derivatives, provide the option"
        "multiple times. e.g. '-d file1.csv -d file2.csv'."
    ),
)
@click.option(
    "-v",
    "--variable",
    multiple=True,
    default=[],
    help=(
        "Return subjects who have non-null responses to these variables. These variable names must "
        "correspond to variable names in the data dictionaries."
        "To specify multiple variables, provide the "
        "option multiple times. e.g '-v FSQ_08 -v WIAT_RC_Stnd'."
    ),
)
@click.argument("data_dict_dir", type=str)
@click.argument("coins_dir", type=str)
def merge(outdir, suffix, deriv_file, variable, data_dict_dir, coins_dir):
    """Create subset versions of phenotypic and derivative data files.

    This utility will find all subjects that have BOTH phenotypic data for the
    requested variables and derivative data in the supplied derivative files.
    It will then create subset datafile for each of these inputs.

    This software assumes that the data dictionaries and phenotypic files are in the
    format provided by the COINS data exchange.

    If you don't currently have access to HBN phenotypic data, follow the instructions at
    http://fcon_1000.projects.nitrc.org/indi/cmi_healthy_brain_network/Pheno_Access.html
    to gain access.

    \b
    Arguments:\b
      DATA_DICT_DIR       Path to the data dictionary directory\b
      COINS_DIR           Path to the COINS data directory\b
    """
    if not (deriv_file or variable):
        raise click.ClickException(
            "You must specify either derivatives (using the -d option) or "
            "variables (using the -v option)."
        )

    if variable:
        pheno_subs, df_pheno = _get_variable_info(
            variable, data_dict_dir, coins_dir, return_df=True
        )

    # We will store sets of subject IDs in a list
    # Each element will be a set of subject IDs that are present for a derivative
    # At the end, we will take the intersection of all of the sets.
    deriv_subs = []
    deriv_dfs = {}
    for dfile in deriv_file:
        deriv_dfs[dfile] = pd.read_csv(dfile)
        deriv_subs.append(set(pd.unique(deriv_dfs[dfile]["subjectID"])))

    if deriv_subs:
        deriv_subs = set.intersection(*deriv_subs)
        if variable:
            subjects = list(pheno_subs.intersection(deriv_subs))
        else:
            subjects = list(deriv_subs)
    else:
        subjects = list(pheno_subs)

    os.makedirs(outdir, exist_ok=True)

    for dfile in deriv_file:
        mask = deriv_dfs[dfile]["subjectID"].isin(subjects)
        deriv_dfs[dfile][mask].to_csv(
            op.join(outdir, suffix.join(op.splitext(op.basename(dfile))))
        )

    df_pheno = df_pheno.filter(items=subjects, axis="index")
    df_pheno.index.rename("subjectID", inplace=True)
    df_pheno.to_csv(op.join(outdir, "pheno" + suffix + ".csv"))


if __name__ == "__main__":
    cli()
