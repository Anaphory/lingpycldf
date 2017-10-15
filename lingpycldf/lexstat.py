#!/usr/bin/env python

"""Analyze pycldf word list using lingpy"""

import tempfile

import argparse

import pycldf
from lingpy.compare.lexstat import LexStat
from pycldf.util import Path
from pycldf.dataset import Dataset
from pycldf.cli import _get_dataset
from clldutils.clilib import ParserError

import column_names


def get_dataset(fname):
    fname = Path(fname)
    if not fname.exists() or not fname.is_file():
        raise ParserError(
            '{:} is not an existing directory'.format(
                fname))
    if fname.suffix == '.json':
        return Dataset.from_metadata(fname)
    return Dataset.from_data(fname)


def lingpy_write(entries):
    lpwl.write(
        "\t".join(map(str, entries)).encode('utf-8'))
    lpwl.write(b"\n")


parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
parser.add_argument("wordlist", type=_get_dataset, default=None, nargs="?")
"""Write a CognatesTable. In later programming, deal with TableSet
descriptions/overwriting data/existing codes etc., but first try to get the
basics done."""
args = parser.parse_args([])

if args.wordlist is None:
    try:
        wordlist = get_dataset("Wordlist-metadata.json")
    except ParserError:
        wordlist = get_dataset("forms.csv")
else:
    wordlist = args.wordlist

lpwl = tempfile.NamedTemporaryFile()

lingpy_write(["ID", "REFERENCE", "DOCULECT", "CONCEPT", "TOKENS"])
reference = wordlist[("FormTable", "id")].name
doculect = wordlist[("FormTable", "languageReference")].name
concept = wordlist[("FormTable", "parameterReference")].name
tokens = wordlist[("FormTable", "soundSequence")].name
for r, row in enumerate(
        wordlist[wordlist.primary_table].iterdicts()):
    if not row[tokens]:
        continue
    lingpy_row = [
        r, row[reference], row[doculect], row[concept], ' '.join(row[tokens])
    ]
    lingpy_write(lingpy_row)

wordlist.add_component("CognateTable")
cognate_table = wordlist["CognateTable"]

lpwl.seek(0, 0)

lexstat = LexStat(lpwl.name)
lexstat.cluster(method="sca", cluster_method="upgma")
lexstat.cluster(method="lexstat", cluster_method="upgma")

cognates = []
for r, row in enumerate(lexstat._data.values()):
    row = {k: row[v] for k, v in lexstat.header.items()}
    cognates.append({
        "ID": r,
        "Form_ID": row["reference"],
        "Cognateset_ID": row["scaid"],
        "Cognateset_source": ["LexStat"]})
cognate_table.write(cognates)

lpwl.close()
