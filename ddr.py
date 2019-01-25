# -*- coding: utf-8 -*-

"""Download, convert, and summarize disease-disease relationships with BEL graph."""

import os
from urllib.request import urlretrieve
from zipfile import ZipFile

import pandas as pd
from pybel import BELGraph
from pybel.dsl import Pathology

data_url = 'http://science.sciencemag.org/highwire/filestream/628238/field_highwire_adjunct_files/1/Datasets_S1-S4.zip'
path = os.path.join(os.path.expanduser('~'), 'Desktop', 'data.zip')
columns = [
    'disease_A',
    'disease_B',
    's_AB (observed)',
    'd_AB (observed)',
    'z (full rand)',
    'p (full rand)',
    'q (full rand)',
    'z (MeSH rand)',
    'p (MeSH rand)',
    'q (MeSH rand)',
]


def download_data() -> None:
    """Download the data from Menche, J., *et al.* 2015."""
    urlretrieve(data_url, path)


def extract_data() -> pd.DataFrame:
    """Load the data."""
    with ZipFile(path) as myzip:
        with myzip.open('data/DataS4_disease_pairs.tsv') as file:
            return pd.read_csv(file, sep="\t", skiprows=33, names=columns)


def make_graph(df: pd.DataFrame) -> BELGraph:
    """Convert the data to a BEL graph."""
    graph = BELGraph()

    for _, (disease_a, disease_b) in df[['disease_A', 'disease_B']].iterrows():
        if not disease_a or not disease_b:
            continue
        graph.add_association(
            Pathology("MeSH", disease_a),
            Pathology("MeSH", disease_b),
            citation="25700523",
            evidence="from ddr",
        )

    return graph


def main():
    """Download, convert, and summarize disease-disease relationships with BEL graph."""
    if not os.path.exists(path):
        download_data()
    df = extract_data()
    graph = make_graph(df)
    graph.summarize()


if __name__ == '__main__':
    main()
