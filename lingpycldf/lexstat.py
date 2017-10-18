#!/usr/bin/env python

"""Analyze pycldf word list using lingpy"""

import json

import argparse

import pycldf
from pycldf.util import Path
from pycldf.dataset import Dataset
from pycldf.cli import _get_dataset

import lingpy
from lingpy.align.sca import Alignments
from lingpy.compare.lexstat import LexStat


def get_dataset(fname):
    """Load a CLDF dataset.

    Load the file as `json` CLDF metadata description file, or as metadata-free
    dataset contained in a single csv file.

    The distinction is made depending on the file extension: `.json` files are
    loaded as metadata descriptions, all other files are matched against the
    CLDF module specifications. Directories are checked for the presence of
    any CLDF datasets in undefined order of the dataset types.

    Parameters
    ----------
    fname : str or Path
        Path to a CLDF dataset

    Returns
    -------
    Dataset
    """
    fname = Path(fname)
    if not fname.exists():
        raise FileNotFoundError(
            '{:} does not exist'.format(fname))
    if fname.suffix == '.json':
        return Dataset.from_metadata(fname)
    return Dataset.from_data(fname)


def to_lingpy(wordlist, replace_tab=" ", replace_newline=" "):
    """Write a CLDF wordlist in LingPy format.

    Write all rows from a CLDF wordlist to a LingPy-readable file. The LingPy
    file parser is extremely naïve and cannot understand quoted cells, so every
    "\\t" (field separator) and "\\n" (row separator) inside the wordlist cells
    will be replaced by `replace_tab` and `replace_newline` respectively.

    NOTE: Currently, this function can only convert the following properties,
    naming the columns as follows: id→REFERENCE, languageReference→DOCULECT,
    parameterReference→CONCEPT, soundSequence→TOKENS.

    The ID values of the original wordlist will be stored in a REFERENCE column
    for LingPy, because LingPy expects its IDs to be subsequent integers.

    Parameters
    ----------
    wordlist : Wordlist
        The CLDF wordlist to be converted
    replace_tab : str
        String to replace tabs inside cells with (Default value = " ")
    replace_newline : str
        String to replace newlines inside cells with (Default value = " ")

    Returns
    -------
    tempfile.NamedTemporaryFile
        A named temporary file containing the LingPy-formatted wordlist
    """
    lpwl = {}

    def lingpy_write(entries):
        lpwl[len(lpwl)] = [
            entry.replace("\t", replace_tab
            ).replace("\n", replace_newline)
            if type(entry) == str
            else entry
            for entry in entries]

    lingpy_write(["ID", "REFERENCE", "DOCULECT", "CONCEPT", "IPA", "TOKENS"])
    reference = wordlist[("FormTable", "id")].name
    doculect = wordlist[("FormTable", "languageReference")].name
    concept = wordlist[("FormTable", "parameterReference")].name
    tokens = wordlist[("FormTable", "soundSequence")].name
    for r, row in enumerate(
            wordlist[wordlist.primary_table].iterdicts()):
        if not row[tokens]:
            continue
        lingpy_row = [
            r+1, row[reference], row[doculect], row[concept], ''.join(row[tokens]), [x for x in row[tokens] if x not in "(,_.-;)"]
        ]
        lingpy_write(lingpy_row)
    return lpwl


def cognatetable_from_lingpy(lingpy, column="cogid"):
    """Generate a CLDF cognate table of data read from LingPy

    Parameters
    ----------
    lingpy : lingpy.LexStat
        The LingPy object containing cognate codes

    column : str
        The name of the LingPy column containing the cognate classes (Default
        value = "cogid")

    Returns
    -------
    cognates : [{`column`: `value`}]
        A list of rows, which can be written to a standard CLDF CognateTable

    """
    cognates = []
    for r, row in enumerate(lexstat._data.values()):
        row = {k: row[v] for k, v in lexstat.header.items()}
        cognates.append({
            "ID": r,
            "Form_ID": row["reference"],
            "Cognateset_ID": row[column],
            "Alignment": row["alignment"],
            "Source": ["LexStat"]})
    return cognates

def find_bad_tokens(wordlist):
    """Collect which bad symbols appear in which forms."""
    bad_tokens = {}
    for k, segments, form_id in wordlist.iter_rows('tokens', "reference"):
        classes = lingpy.tokens2class(segments, 'dolgo')
        for token, cls in zip(segments, classes):
            if cls == "0":
                bad_tokens.setdefault(token, []).append(form_id)
    return bad_tokens

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0] + """

    Generate a CognatesTable, using LingPy's LexStat algorithm.

    TODO: Expose more cluster() arguments; deal with TableSet descriptions;
    test""")
    parser.add_argument("wordlist", type=get_dataset, default=None, nargs="?")
    parser.add_argument(
        "--method", choices={'sca', 'lexstat', 'edit-dist', 'turchin'}, default="sca",
        help="The method that shall be used for the calculation of distances")
    parser.add_argument(
        "--cluster-method", choices={'upgma', 'single', 'complete', 'mcl', 'infomap'},
        default='upgma',
        help="The method used to identify clusters")
    parser.add_argument(
        "--threshold", type=float, default=False,
        help="Use this threshold for the cluster algorithm")
    parser.add_argument(
        "--overwrite", action="store_true", default=False,
        help="Overwrite an existing CognateTable if one exists.")
    parser.add_argument(
        "--bad-tokens-log", type=argparse.FileType("w"), default=None,
        help="File to write a list of bad tokens to")
    args = parser.parse_args()

    # Load the word list into a LingPy compatible format
    if args.wordlist is None:
        try:
            wordlist = get_dataset("Wordlist-metadata.json")
        except FileNotFoundError:
            wordlist = get_dataset("forms.csv")
    else:
        wordlist = args.wordlist
    try:
        wordlist.add_component("CognateTable")
    except ValueError:
        if args.overwrite:
            pass
        else:
            print("DataSet already has a CognateTable. To drop existing cognate data, use `--overwrite`.")
            sys.exit(2)
    lpwl = to_lingpy(wordlist)

    # Use LingPy functionality
    lexstat = LexStat(lpwl, check=False, segments="tokens")
    if args.bad_tokens_log:
        json.dump(find_bad_tokens(lexstat), args.bad_tokens_log)

    # Prepare analysis
    if args.method != 'sca':
        lexstat.get_scorer(preprocessing=False, runs=10000, ratio=(2,1), vscale=1.0)
    lexstat.cluster(method=args.method, cluster_method=args.cluster_method, ref="cogid",
                    threshold=args.threshold)
    lexstat = Alignments(lexstat, segments="tokens")
    lexstat.align(model="sca")
    lexstat.output("tsv", filename="with_lexstat_and_alignment")

    # Create new CognateTable and write it to there
    cognate_table = wordlist["CognateTable"]
    cognate_table.write(cognatetable_from_lingpy(lexstat))

