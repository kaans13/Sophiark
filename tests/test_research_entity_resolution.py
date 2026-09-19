from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock

from src.research.config import DEFAULT_RESEARCH_CONFIG, ResearchConfig
from src.research.entity_resolution import (
    EntityResolver,
    ExternalCandidateSet,
    external_cache_key,
    load_human_project_index,
    load_mouse_project_index,
)
from src.research.models import (
    AliasType,
    Entity,
    EntityAlias,
    EntityType,
    ResolutionStatus,
)


HUMAN_HEADER = (
    "HGNC ID\tApproved symbol\tApproved name\tStatus\tLocus type\tPrevious symbols\t"
    "Alias symbols\tChromosome\tNCBI Gene ID\tEnsembl gene ID\n"
)


class EntityResolverTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        symbols = self.root / "symbols.csv"
        symbols.write_text(
            "gene,Symbol\n"
            "ENSP_A,SHARED\n"
            "ENSP_B,AMB\n"
            "ENSP_C,AMB\n"
            "ENSP_CANON,CANON_SYMBOL\n"
            "ENSP_OTHER,ENSP_CANON\n",
            encoding="utf-8",
        )
        ensg = self.root / "ensp_to_ensg.csv"
        ensg.write_text(
            "ENSP,ENSG\n"
            "ENSP_A,ENSG_A\n"
            "ENSP_B,ENSG_B\n"
            "ENSP_C,ENSG_C\n"
            "ENSP_CANON,ENSG_CANON\n"
            "ENSP_OTHER,ENSG_OTHER\n",
            encoding="utf-8",
        )
        hgnc = self.root / "hgnc.tsv"
        hgnc.write_text(
            HUMAN_HEADER
            + "HGNC:1\tSHARED\tShared human protein\tApproved\tgene with protein product\tOLD_SHARED\tA_ALIAS\t1\t101\tENSG_A\n"
            + "HGNC:2\tAMB\tAmbiguous B\tApproved\tgene with protein product\tOLD_AMB\tDUP_ALIAS\t2\t102\tENSG_B\n"
            + "HGNC:3\tAMB\tAmbiguous C\tApproved\tgene with protein product\t\tDUP_ALIAS\t3\t103\tENSG_C\n",
            encoding="utf-8",
        )
        mouse = self.root / "mouse_info.tsv"
        mouse.write_text(
            "#string_protein_id\tpreferred_name\n"
            "10090.ENSMUSP_A\tSHARED\n"
            "10090.ENSMUSP_B\tMouseOnly\n",
            encoding="utf-8",
        )

        self.human_index = load_human_project_index(
            self.root,
            symbol_path=symbols,
            ensp_to_ensg_path=ensg,
            hgnc_path=hgnc,
        )
        self.mouse_index = load_mouse_project_index(self.root, string_info_path=mouse)
        self.resolver = EntityResolver((self.human_index, self.mouse_index))

    def tearDown(self) -> None:
        self.temp.cleanup()

    def test_exact_canonical_precedes_a_same_text_symbol(self) -> None:
        result = self.resolver.resolve("ENSP_CANON", taxon_id=9606)
        self.assertEqual(result.status, ResolutionStatus.EXACT)
        self.assertEqual(result.entity.canonical_id, "ENSP_CANON")
        self.assertEqual(result.source, "project_canonical")

    def test_existing_project_mapping_and_symbol_are_exact(self) -> None:
        by_ensg = self.resolver.resolve("ENSG_A", taxon_id=9606)
        by_symbol = self.resolver.resolve("shared", taxon_id=9606)
        by_ncbi = self.resolver.resolve("101", taxon_id=9606)
        self.assertEqual(by_ensg.status, ResolutionStatus.EXACT)
        self.assertEqual(by_symbol.status, ResolutionStatus.EXACT)
        self.assertEqual(by_ncbi.status, ResolutionStatus.EXACT)
        self.assertEqual(by_ensg.entity.canonical_id, "ENSP_A")
        self.assertEqual(by_symbol.entity.canonical_id, "ENSP_A")
        self.assertEqual(by_ncbi.entity.canonical_id, "ENSP_A")

    def test_known_alias_has_alias_status_and_provenance(self) -> None:
        result = self.resolver.resolve("old_shared", taxon_id=9606)
        self.assertEqual(result.status, ResolutionStatus.ALIAS)
        self.assertEqual(result.entity.canonical_id, "ENSP_A")
        self.assertEqual(result.matched_alias.alias_type, AliasType.SYNONYM)
        self.assertIn("HGNC", result.matched_alias.source)

    def test_ambiguous_symbol_and_alias_never_select_a_candidate(self) -> None:
        for query in ("AMB", "DUP_ALIAS"):
            with self.subTest(query=query):
                result = self.resolver.resolve(query, taxon_id=9606)
                self.assertEqual(result.status, ResolutionStatus.AMBIGUOUS)
                self.assertIsNone(result.entity)
                self.assertEqual(
                    {candidate.canonical_id for candidate in result.candidates},
                    {"ENSP_B", "ENSP_C"},
                )

    def test_ambiguity_stops_before_external_fallback(self) -> None:
        fallback = Mock(
            return_value=ExternalCandidateSet(
                candidates=(
                    Entity(9606, EntityType.PROTEIN, "EXTERNAL_CHOICE", symbol="AMB"),
                ),
                source="test_provider",
            )
        )
        resolver = EntityResolver((self.human_index,), external_fallback=fallback)
        result = resolver.resolve("AMB", taxon_id=9606, allow_external=True)
        self.assertEqual(result.status, ResolutionStatus.AMBIGUOUS)
        fallback.assert_not_called()

    def test_same_symbol_is_resolved_only_within_requested_species(self) -> None:
        human = self.resolver.resolve("SHARED", taxon_id=9606)
        mouse = self.resolver.resolve("SHARED", taxon_id=10090)
        self.assertEqual(human.entity.canonical_id, "ENSP_A")
        self.assertEqual(mouse.entity.canonical_id, "ENSMUSP_A")
        self.assertEqual(human.entity.logical_key, "9606:protein:ENSP_A")
        self.assertEqual(mouse.entity.logical_key, "10090:protein:ENSMUSP_A")
        prefixed = self.resolver.resolve("10090.ENSMUSP_A", taxon_id=10090)
        self.assertEqual(prefixed.status, ResolutionStatus.ALIAS)

    def test_unresolved_does_not_use_symbol_prefix_as_family_inference(self) -> None:
        result = self.resolver.resolve("SHARED2", taxon_id=9606)
        self.assertEqual(result.status, ResolutionStatus.UNRESOLVED)
        self.assertIsNone(result.entity)
        self.assertEqual(result.candidates, ())

    def test_local_external_cache_is_used_before_an_optional_fallback(self) -> None:
        cached_entity = Entity(9606, EntityType.PROTEIN, "ENSP_CACHE", symbol="CACHE_ME")
        cache_alias = EntityAlias(AliasType.UNIPROT, "P12345", "local UniProt cache")
        cache = {
            external_cache_key("P12345", 9606): ExternalCandidateSet(
                candidates=(cached_entity,),
                source="local UniProt cache",
                matched_alias=cache_alias,
            )
        }
        fallback = Mock()
        resolver = EntityResolver(
            (self.human_index,), local_external_cache=cache, external_fallback=fallback,
        )
        result = resolver.resolve("P12345", taxon_id=9606, allow_external=True)
        self.assertEqual(result.status, ResolutionStatus.EXTERNAL_RESOLVED)
        self.assertEqual(result.entity, cached_entity)
        self.assertEqual(result.source, "local UniProt cache")
        fallback.assert_not_called()

    def test_negative_cache_and_external_opt_in_prevent_implicit_requests(self) -> None:
        key = external_cache_key("NO_HIT", 9606)
        fallback = Mock(
            return_value=ExternalCandidateSet(
                candidates=(Entity(9606, EntityType.PROTEIN, "ENSP_EXT"),),
                source="test_provider",
            )
        )
        cached = EntityResolver(
            (self.human_index,),
            local_external_cache={key: ExternalCandidateSet((), "negative cache")},
            external_fallback=fallback,
        )
        self.assertEqual(
            cached.resolve("NO_HIT", taxon_id=9606, allow_external=True).status,
            ResolutionStatus.UNRESOLVED,
        )
        fallback.assert_not_called()

        uncached = EntityResolver((self.human_index,), external_fallback=fallback)
        self.assertEqual(
            uncached.resolve("NO_HIT", taxon_id=9606).status,
            ResolutionStatus.UNRESOLVED,
        )
        fallback.assert_not_called()
        resolved = uncached.resolve("NO_HIT", taxon_id=9606, allow_external=True)
        self.assertEqual(resolved.status, ResolutionStatus.EXTERNAL_RESOLVED)
        fallback.assert_called_once()

    def test_external_candidates_must_preserve_species_and_explicit_contract(self) -> None:
        cross_species = Mock(
            return_value=ExternalCandidateSet(
                candidates=(Entity(10090, EntityType.PROTEIN, "ENSMUSP_WRONG"),),
                source="bad_provider",
            )
        )
        result = EntityResolver(external_fallback=cross_species).resolve(
            "WRONG", taxon_id=9606, allow_external=True,
        )
        self.assertEqual(result.status, ResolutionStatus.UNRESOLVED)
        self.assertIn("taxon", result.message)

        invalid_contract = Mock(return_value=Entity(9606, EntityType.PROTEIN, "ENSP_GUESS"))
        result = EntityResolver(external_fallback=invalid_contract).resolve(
            "GUESS", taxon_id=9606, allow_external=True,
        )
        self.assertEqual(result.status, ResolutionStatus.UNRESOLVED)
        self.assertIn("invalid candidate contract", result.message)

    def test_empty_query_and_missing_files_are_safe(self) -> None:
        result = self.resolver.resolve("  ", taxon_id=9606)
        self.assertEqual(result.status, ResolutionStatus.UNRESOLVED)
        missing = load_human_project_index(self.root / "not-present")
        self.assertEqual(missing.entity_count, 0)


class ResearchConfigTests(unittest.TestCase):
    def test_external_context_and_all_providers_are_default_off(self) -> None:
        config = DEFAULT_RESEARCH_CONFIG
        self.assertFalse(config.research_context_enabled)
        self.assertFalse(config.external_context_enabled)
        for provider in ("uniprot", "interpro", "quickgo", "reactome", "europepmc", "unknown"):
            self.assertFalse(config.provider_enabled(provider))

    def test_provider_requires_both_global_and_provider_flags(self) -> None:
        provider_only = ResearchConfig(provider_uniprot_enabled=True)
        global_and_provider = ResearchConfig(
            external_context_enabled=True,
            provider_uniprot_enabled=True,
        )
        self.assertFalse(provider_only.provider_enabled("UniProt"))
        self.assertTrue(global_and_provider.provider_enabled("Uni_Prot"))
        self.assertTrue(global_and_provider.presentation_limits.family_threshold_presentation_only)


if __name__ == "__main__":
    unittest.main()
