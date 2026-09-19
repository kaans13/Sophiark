from __future__ import annotations

import csv
import json
from pathlib import Path
import pickle
import socket
import sqlite3
from tempfile import TemporaryDirectory

import pytest

from src.research.providers.local import (
    ComplexMembership,
    FamilyMembership,
    HUMAN_TAXON_ID,
    LocalComplexProvider,
    LocalEvidenceAdapters,
    LocalHGNCFamilyProvider,
    LocalMyGeneProvider,
    LocalProviderStatus,
    LocalRegulatoryProvider,
    LocalUniProtProvider,
    MOUSE_TAXON_ID,
)


def _pickle(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as handle:
        pickle.dump(value, handle)


def _tsv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_construction_is_io_free_and_missing_sources_never_use_network(monkeypatch) -> None:
    with TemporaryDirectory() as directory:
        root = Path(directory)
        before = tuple(root.rglob("*"))
        adapters = LocalEvidenceAdapters.from_project_root(root)
        assert tuple(root.rglob("*")) == before

        def forbidden_network(*_args, **_kwargs):
            raise AssertionError("local providers must never access the network")

        monkeypatch.setattr(socket, "create_connection", forbidden_network)
        results = (
            adapters.mygene.get_annotations(["ENSP1"], taxon_id=HUMAN_TAXON_ID),
            adapters.families.get_families(["RAMP1"], taxon_id=HUMAN_TAXON_ID),
            adapters.complexes.get_memberships(["RAMP1"], taxon_id=HUMAN_TAXON_ID),
            adapters.uniprot.get_context(["ENSP1"], taxon_id=HUMAN_TAXON_ID),
            adapters.regulatory.get_relations(["RAMP1"], taxon_id=HUMAN_TAXON_ID),
        )
        assert all(result.status is LocalProviderStatus.UNAVAILABLE for result in results)
        assert all(result.missing_sources for result in results)
        assert tuple(root.rglob("*")) == before


def test_mygene_cache_is_exact_and_species_scoped() -> None:
    with TemporaryDirectory() as directory:
        root = Path(directory)
        human = root / "human.pkl"
        mouse = root / "mouse.pkl"
        _pickle(
            human,
            {
                "ENSP-H": {
                    "symbol": "RAMP1",
                    "name": "human receptor activity modifying protein 1",
                    "go_bp": [{"id": "GO:1", "term": "human process"}],
                    "go_cc": ["plasma membrane"],
                    "go_mf": [],
                }
            },
        )
        _pickle(
            mouse,
            {
                "ENSMUSP-M": {
                    "symbol": "Ramp1",
                    "name": "mouse receptor activity modifying protein 1",
                    "go_bp": ["mouse process"],
                    "go_cc": [],
                    "go_mf": [],
                }
            },
        )
        provider = LocalMyGeneProvider(human_cache_path=human, mouse_cache_path=mouse)

        human_hit = provider.get_annotations(["ENSP-H"], taxon_id=HUMAN_TAXON_ID)
        mouse_hit = provider.get_annotations(["ENSMUSP-M"], taxon_id=MOUSE_TAXON_ID)
        assert human_hit.status is LocalProviderStatus.AVAILABLE
        assert human_hit.records[0].symbol == "RAMP1"
        assert human_hit.records[0].provenance.details["taxon_id"] == HUMAN_TAXON_ID
        assert mouse_hit.records[0].symbol == "Ramp1"
        assert mouse_hit.records[0].provenance.details["taxon_id"] == MOUSE_TAXON_ID

        assert provider.get_annotations(["ENSMUSP-M"], taxon_id=HUMAN_TAXON_ID).status is LocalProviderStatus.NO_MATCH
        assert provider.get_annotations(["ENSP-H"], taxon_id=MOUSE_TAXON_ID).status is LocalProviderStatus.NO_MATCH

        go_record = human_hit.records[0].go_bp[0]
        with pytest.raises(TypeError):
            go_record["term"] = "mutated"


def test_hgnc_is_authoritative_and_complex_membership_is_not_family() -> None:
    with TemporaryDirectory() as directory:
        root = Path(directory)
        hgnc = root / "hgnc.tsv"
        _tsv(
            hgnc,
            [
                "HGNC ID",
                "Approved symbol",
                "Approved name",
                "Previous symbols",
                "Alias symbols",
                "NCBI Gene ID",
                "Ensembl gene ID",
                "Group ID",
                "Group name",
            ],
            [
                {
                    "HGNC ID": "HGNC:1",
                    "Approved symbol": "RAMP1",
                    "Approved name": "receptor activity modifying protein 1",
                    "Previous symbols": "OLD_RAMP1",
                    "Alias symbols": "ALIAS_RAMP1",
                    "NCBI Gene ID": "10267",
                    "Ensembl gene ID": "ENSG1",
                    "Group ID": "500",
                    "Group name": "Receptor activity modifying proteins",
                }
            ],
        )
        family = LocalHGNCFamilyProvider(hgnc)
        hit = family.get_families(["OLD_RAMP1"], taxon_id=HUMAN_TAXON_ID)
        assert hit.status is LocalProviderStatus.AVAILABLE
        assert isinstance(hit.records[0], FamilyMembership)
        assert hit.records[0].matched_via == "previous_symbol"
        assert hit.records[0].group_id == "500"

        # Similar prefixes never synthesize an authoritative family.
        no_guess = family.get_families(["RAMP9"], taxon_id=HUMAN_TAXON_ID)
        assert no_guess.status is LocalProviderStatus.NO_MATCH
        assert "no symbol-prefix heuristic" in no_guess.message
        assert family.get_families(["Ramp1"], taxon_id=MOUSE_TAXON_ID).status is LocalProviderStatus.UNSUPPORTED_TAXON

        processed = root / "complex.pkl"
        portal = root / "complex_portal.tsv"
        _pickle(processed, {"RAMP1": ["RAMP1-CALCRL complex"]})
        _tsv(
            portal,
            [
                "#Complex ac",
                "Recommended name",
                "Taxonomy identifier",
                "Identifiers (and stoichiometry) of molecules in complex",
                "Expanded participant list",
            ],
            [
                {
                    "#Complex ac": "CPX-1",
                    "Recommended name": "RAMP1-CALCRL complex",
                    "Taxonomy identifier": "9606",
                    "Identifiers (and stoichiometry) of molecules in complex": "P12345(1)|Q99999(1)",
                    "Expanded participant list": "P12345(1)|Q99999(1)",
                }
            ],
        )
        complexes = LocalComplexProvider(
            human_processed_path=processed,
            complex_portal_path=portal,
        ).get_memberships(["RAMP1", "P12345"], taxon_id=HUMAN_TAXON_ID)
        assert {record.identifier_type for record in complexes.records} == {"gene_symbol", "uniprot"}
        assert all(isinstance(record, ComplexMembership) for record in complexes.records)
        assert all(record.evidence_type == "COMPLEX_MEMBERSHIP" for record in complexes.records)
        assert all(not isinstance(record, FamilyMembership) for record in complexes.records)
        assert all(record.provenance.details["not_family_annotation"] is True for record in complexes.records)


def test_partial_uniprot_uses_read_only_sqlite_and_declares_missing_fields() -> None:
    with TemporaryDirectory() as directory:
        root = Path(directory)
        database = root / "structural_binding.db"
        connection = sqlite3.connect(database)
        connection.executescript(
            """
            CREATE TABLE ensp_uniprot_map(ensp TEXT PRIMARY KEY, uniprot TEXT);
            CREATE TABLE structural_descriptors(
                ensp TEXT PRIMARY KEY, uniprot TEXT, length INTEGER,
                mean_plddt REAL, sasa REAL, largest_pocket_volume REAL,
                hydrophobicity REAL, radius_gyration REAL
            );
            CREATE TABLE bindingdb_summary(
                uniprot TEXT PRIMARY KEY, num_ligands INTEGER,
                avg_ki REAL, min_ki REAL, avg_ic50 REAL, min_ic50 REAL
            );
            INSERT INTO ensp_uniprot_map VALUES ('ENSP1', 'P12345');
            INSERT INTO structural_descriptors VALUES
                ('ENSP1', 'P12345', 420, 91.5, 123.4, 55.0, -0.2, 20.1);
            INSERT INTO bindingdb_summary VALUES
                ('P12345', 7, 12.0, 1.0, 30.0, 2.0);
            """
        )
        connection.commit()
        connection.close()
        before = (database.stat().st_size, database.stat().st_mtime_ns)

        provider = LocalUniProtProvider(database)
        result = provider.get_context(["ENSP1", "P12345"], taxon_id=HUMAN_TAXON_ID)
        after = (database.stat().st_size, database.stat().st_mtime_ns)
        assert result.status is LocalProviderStatus.PARTIAL
        assert len(result.records) == 1
        assert result.records[0].uniprot_accession == "P12345"
        assert result.records[0].structure["mean_plddt"] == 91.5
        assert result.records[0].ligand_summary["num_ligands"] == 7
        assert "protein_family" in result.records[0].unavailable_fields
        assert "function" in result.records[0].unavailable_fields
        assert before == after
        assert not Path(str(database) + "-wal").exists()

        assert provider.get_context(["P12345"], taxon_id=MOUSE_TAXON_ID).status is LocalProviderStatus.UNSUPPORTED_TAXON
        missing = root / "missing.db"
        assert LocalUniProtProvider(missing).get_context(["P1"], taxon_id=HUMAN_TAXON_ID).status is LocalProviderStatus.UNAVAILABLE
        assert not missing.exists()


def test_regulatory_sources_filter_taxon_and_preserve_source_types() -> None:
    with TemporaryDirectory() as directory:
        root = Path(directory)
        human_trrust = root / "human_trrust.pkl"
        mouse_trrust = root / "mouse_trrust.pkl"
        _pickle(human_trrust, {"SHARED": [("HUMAN_TF", "Activation")]})
        _pickle(mouse_trrust, {"Shared": [("MouseTf", "Activation")]})

        signor = root / "signor.tsv"
        signor_fields = [
            "ENTITYA", "IDA", "ENTITYB", "IDB", "EFFECT", "MECHANISM",
            "TAX_ID", "PMID", "DIRECT", "SCORE", "SIGNOR_ID",
        ]
        _tsv(
            signor,
            signor_fields,
            [
                {
                    "ENTITYA": "HUMAN_SIGNOR",
                    "IDA": "P1",
                    "ENTITYB": "SHARED",
                    "IDB": "P2",
                    "EFFECT": "up-regulates activity",
                    "MECHANISM": "binding",
                    "TAX_ID": "9606",
                    "PMID": "11",
                    "DIRECT": "YES",
                    "SCORE": "0.9",
                    "SIGNOR_ID": "SIGNOR-H",
                },
                {
                    "ENTITYA": "MOUSE_SIGNOR",
                    "IDA": "Q1",
                    "ENTITYB": "Shared",
                    "IDB": "Q2",
                    "EFFECT": "down-regulates activity",
                    "MECHANISM": "binding",
                    "TAX_ID": "10090",
                    "PMID": "22",
                    "DIRECT": "YES",
                    "SCORE": "0.8",
                    "SIGNOR_ID": "SIGNOR-M",
                },
            ],
        )

        omni_fields = [
            "source", "target", "source_genesymbol", "target_genesymbol",
            "is_directed", "is_stimulation", "is_inhibition",
            "consensus_direction", "consensus_stimulation",
            "consensus_inhibition", "sources", "references",
        ]
        human_omni, mouse_omni = root / "omni_h.tsv", root / "omni_m.tsv"
        _tsv(
            human_omni,
            omni_fields,
            [{
                "source": "P3", "target": "P4", "source_genesymbol": "HUMAN_OMNI",
                "target_genesymbol": "SHARED", "is_directed": "1", "is_stimulation": "1",
                "is_inhibition": "0", "consensus_direction": "1",
                "consensus_stimulation": "1", "consensus_inhibition": "0",
                "sources": "SIGNOR;Reactome", "references": "SIGNOR:33",
            }],
        )
        _tsv(
            mouse_omni,
            omni_fields,
            [{
                "source": "Q3", "target": "Q4", "source_genesymbol": "MOUSE_OMNI",
                "target_genesymbol": "Shared", "is_directed": "1", "is_stimulation": "0",
                "is_inhibition": "1", "consensus_direction": "1",
                "consensus_stimulation": "0", "consensus_inhibition": "1",
                "sources": "SIGNOR", "references": "SIGNOR:44",
            }],
        )
        human_meta, mouse_meta = root / "omni_h.json", root / "omni_m.json"
        human_meta.write_text(json.dumps({"organism_taxid": 9606, "release": "fixture-h"}), encoding="utf-8")
        mouse_meta.write_text(json.dumps({"organism_taxid": 10090, "release": "fixture-m"}), encoding="utf-8")

        provider = LocalRegulatoryProvider(
            human_trrust_path=human_trrust,
            mouse_trrust_path=mouse_trrust,
            signor_raw_path=signor,
            human_omnipath_path=human_omni,
            mouse_omnipath_path=mouse_omni,
            human_omnipath_metadata_path=human_meta,
            mouse_omnipath_metadata_path=mouse_meta,
        )
        human = provider.get_relations(["SHARED"], taxon_id=HUMAN_TAXON_ID)
        mouse = provider.get_relations(["Shared"], taxon_id=MOUSE_TAXON_ID)

        assert human.status is LocalProviderStatus.AVAILABLE
        assert mouse.status is LocalProviderStatus.AVAILABLE
        assert {record.database for record in human.records} == {"TRRUST", "SIGNOR", "OmniPath"}
        assert {record.database for record in mouse.records} == {"TRRUST", "SIGNOR", "OmniPath"}
        assert all(record.taxon_id == HUMAN_TAXON_ID for record in human.records)
        assert all(record.taxon_id == MOUSE_TAXON_ID for record in mouse.records)
        assert all("MOUSE" not in record.source_symbol.upper() for record in human.records)
        assert all("HUMAN" not in record.source_symbol.upper() for record in mouse.records)
        assert all(record.evidence_type == "DIRECTED_REGULATION" for record in human.records + mouse.records)
        assert all(record.provenance.locator for record in human.records + mouse.records)
