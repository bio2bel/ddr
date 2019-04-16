# -*- coding: utf-8 -*-

"""Download, convert, and summarize disease-disease relationships with BEL graph."""

import os
from typing import Mapping, Optional
from zipfile import ZipFile

import pandas as pd
from tqdm import tqdm

from bio2bel.downloading import make_downloader
from bio2bel.manager.bel_manager import BELManagerMixin
from pybel import BELGraph
from pybel.dsl import Pathology
from .constants import DATA_DIR, MODULE_NAME

DATA_URL = 'http://science.sciencemag.org/highwire/filestream/628238/field_highwire_adjunct_files/1/Datasets_S1-S4.zip'
DATA_PATH = os.path.join(DATA_DIR, 'data.zip')
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

download_data = make_downloader(DATA_URL, DATA_PATH)
"""Download the data from Menche, J., *et al.* 2015."""


def extract_data(path: Optional[str] = None) -> pd.DataFrame:
    """Load the data."""
    with ZipFile(path or DATA_PATH) as zip_file:
        with zip_file.open('data/DataS4_disease_pairs.tsv') as file:
            return pd.read_csv(file, sep="\t", skiprows=33, names=columns)


def make_graph():
    """Convert the data to a BEL graph."""
    if not os.path.exists(DATA_PATH):
        download_data()
    df = extract_data()
    return _make_graph(df)


def _make_graph(df: pd.DataFrame,
                use_tqdm: bool = False,
                min_network_separation: float = 0,
                ) -> BELGraph:
    graph = BELGraph(
        name="Disease-disease relationships",
        version="1.0.0",
    )
    '''
    From lit: In summary, we conclude that the proposed separation measure sAB offers a robust quantification
    of the network-based relationship between diseases. As expected, we find that most disease pairs
    are clearly separated (sAB > 0), however, we also identified a considerable number of disease pairs
    with statistically significant overlap (sAB < 0).
    '''
    it = df[['disease_A', 'disease_B', 's_AB (observed)']].iterrows()
    if use_tqdm:
        it = tqdm(it, total=len(df.index), desc='generating BEL')
    for _, (disease_a, disease_b, network_separation) in it:
        if not disease_a or not disease_b or network_separation > min_network_separation:
            continue
        graph.add_association(
            Pathology("mesh", disease_a),
            Pathology("mesh", disease_b),
            citation="25700523",
            evidence="from ddr",
            annotations={
                'bio2bel': MODULE_NAME,
                's_AB': network_separation,
            }
        )

    return graph


class Manager(BELManagerMixin):
    """Disease-disease relationships."""

    def __init__(self, *args, **kwargs):
        self.graph = make_graph()

    @classmethod
    def _get_connection(self):
        pass

    def summarize(self) -> Mapping[str, int]:
        """Summarize the database."""
        return dict(
            associations=self.count_relations(),
            diseases=self.count_diseases(),
        )

    @staticmethod
    def is_populated() -> bool:
        """Return if the database is populated."""
        return True

    def count_diseases(self) -> int:
        """Count the number of diseases in the database."""
        return self.graph.number_of_nodes()

    def count_relations(self) -> int:
        """Count the number of disease-disease relationships."""
        return self.graph.number_of_edges()

    def to_bel(self) -> BELGraph:
        """Export as a BEL graph"""
        return self.graph


main = Manager.get_cli()

if __name__ == '__main__':
    main()
