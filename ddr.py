# -*- coding: utf-8 -*-

"""Download, convert, and summarize disease-disease relationships with BEL graph."""

import os
import sys
from urllib.request import urlretrieve
from zipfile import ZipFile

import click
import pandas as pd
from pybel import BELGraph
from pybel.dsl import Pathology
from tqdm import tqdm

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


def make_graph():
    """Convert the data to a BEL graph."""
    if not os.path.exists(path):
        download_data()
    df = extract_data()
    return _make_graph(df)


def _make_graph(df: pd.DataFrame) -> BELGraph:
    graph = BELGraph(
        name="disease-disease relationships",
        version="1.0.0",
    )
    '''
    From lit: In summary, we conclude that the proposed separation measure sAB offers a robust quantification
    of the network-based relationship between diseases. As expected, we find that most disease pairs
    are clearly separated (sAB > 0), however, we also identified a considerable number of disease pairs
    with statistically significant overlap (sAB < 0).
    '''
    for _, (disease_a, disease_b, netwrk_sep) in tqdm(df[['disease_A', 'disease_B', 's_AB (observed)']].iterrows(),
                                                      total=len(df.index)):
        if not disease_a or not disease_b or netwrk_sep > 0:
            continue
        graph.add_association(
            Pathology("MeSH", disease_a),
            Pathology("MeSH", disease_b),
            citation="25700523",
            evidence="from ddr",
        )

    return graph


@click.command()
@click.option("-o", "--file", type=click.File("w"), default=sys.stdout, help="file to output")
@click.option("--format", default="nodelink")
def main(file, format):
    """Download, convert, and summarize disease-disease relationships with BEL graph."""
    graph = make_graph()
    graph.summarize()
    graph.serialize(format, file)


if __name__ == '__main__':
    main()
