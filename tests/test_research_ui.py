"""Offline, injected-UI tests for the Research Explorer projections."""

from __future__ import annotations

from dataclasses import FrozenInstanceError, replace
import io
from types import SimpleNamespace
import unittest
import zipfile

import pandas as pd

from src.research.config import PresentationLimits, ResearchConfig
from src.research.contracts import CONTRACTS
from src.research.facade import (
    ExplorerSection,
    ResearchExportFormat,
    ResearchExplorerView,
    build_research_export_tables,
    build_research_explorer_view,
    research_explorer_download_payload,
    render_research_explorer,
)
from src.research.local_context import (
    LocalAnnotation,
    LocalEntityContext,
    LocalResearchContext,
    LocalSourceStatus,
)
from src.research.models import (
    Entity,
    EntityResolution,
    EntityType,
    EvidenceAssertion,
    Observation,
    ObservationScope,
    ObservationType,
    ProvenanceKind,
    ProvenanceRecord,
    RelationshipType,
    ResolutionStatus,
)
from src.research.observations import DetectedObservation, ObservationDetectionResult
from src.research.providers import (
    ComplexMembership,
    FamilyMembership,
    LocalProviderStatus,
)
from src.research.relationships import RelationshipBuildResult
from src.research.semantic_adapter import (
    FunctionalContextKind,
    FunctionalContextRow,
    ResponseDirection,
    ResultRole,
    ResultSource,
    SemanticRecord,
    SemanticSimulationResult,
)


class FakeUI:
    def __init__(self, section: str = "Overview", selected_id: str | None = None) -> None:
        self.section = section
        self.selected_id = selected_id
        self.calls: list[tuple[str, tuple, dict]] = []

    def _record(self, name, *args, **kwargs):
        self.calls.append((name, args, kwargs))

    def segmented_control(self, *args, **kwargs):
        self._record("segmented_control", *args, **kwargs)
        return self.section

    def selectbox(self, _label, options, **kwargs):
        self._record("selectbox", _label, options, **kwargs)
        return self.selected_id if self.selected_id in options else options[0]

    def __getattr__(self, name):
        def recorder(*args, **kwargs):
            self._record(name, *args, **kwargs)
        return recorder

    @property
    def methods(self) -> tuple[str, ...]:
        return tuple(call[0] for call in self.calls)

    @property
    def text(self) -> str:
        return " ".join(str(value) for _, args, kwargs in self.calls for value in (*args, *kwargs.values()))


def _prov(source: str, kind: ProvenanceKind = ProvenanceKind.LOCAL_ANNOTATION) -> ProvenanceRecord:
    return ProvenanceRecord(source, kind, version="fixture-v1", locator=f"fixture://{source}")


def _resolution(entity: Entity) -> EntityResolution:
    return EntityResolution(
        entity.canonical_id, entity.taxon_id, entity.entity_type,
        ResolutionStatus.EXACT, entity=entity, source="fixture",
    )


def _bundle(*, observation_count: int = 8, empty: bool = False):
    config = ResearchConfig(research_context_enabled=True)
    snapshot_prov = _prov("snapshot", ProvenanceKind.SOPHIARK_COMPUTED)
    hgnc_prov = _prov("HGNC", ProvenanceKind.CURATED_DATABASE)
    ontology_prov = _prov("GO", ProvenanceKind.ONTOLOGY)
    reactome_prov = _prov("Reactome", ProvenanceKind.PATHWAY)
    literature_prov = _prov("Europe PMC", ProvenanceKind.LITERATURE)
    scope = ObservationScope("top_n", 20, 2, None, 20, "Pancreas", "Homo sapiens", 9606)

    entities = (
        Entity(9606, EntityType.PROTEIN, "ENSP_T", "TARGET"),
        Entity(9606, EntityType.PROTEIN, "ENSP_1", "R1"),
        Entity(9606, EntityType.PROTEIN, "ENSP_2", "R2"),
    )
    roles = (
        (ResultRole.PERTURBATION_TARGET,),
        (ResultRole.POSITIVE_REDISTRIBUTION, ResultRole.FDR_SUPPORTED),
        (ResultRole.NETWORK_LOSS,),
    )
    directions = (
        ResponseDirection.UNAVAILABLE,
        ResponseDirection.GAIN,
        ResponseDirection.LOSS,
    )
    records = tuple(
        SemanticRecord(
            ResultSource.REPORT, index, entity.canonical_id, entity.symbol,
            roles[index], directions[index],
            values={"Gümrük_Kapisi": index == 1},
        )
        for index, entity in enumerate(entities)
    )
    family_records = (
        FamilyMembership(9606, "R1", "R1", "HGNC:7", "Fixture family", "approved_symbol", (), None, None, hgnc_prov),
        FamilyMembership(9606, "R2", "R2", "HGNC:7", "Fixture family", "approved_symbol", (), None, None, hgnc_prov),
    )
    annotations = (
        (
            LocalAnnotation("go_biological_process", {"id": "GO:1", "term": "Shared process"}, "GO", provenance=(ontology_prov,)),
            LocalAnnotation("reactome_pathway", {"id": "R-HSA-1", "name": "Reactome fixture"}, "Reactome", provenance=(reactome_prov,)),
        ),
        (
            LocalAnnotation("authoritative_family_membership", family_records[0], "HGNC", provenance=(hgnc_prov,)),
            LocalAnnotation("complex_membership", ComplexMembership(9606, "ENSP_1", "ensembl", "CPX-1", "Not a family", hgnc_prov), "Complex Portal", provenance=(hgnc_prov,)),
            LocalAnnotation("gene_name", "Fixture response protein", "local_mygene", provenance=(hgnc_prov,)),
            LocalAnnotation("go_biological_process", {"id": "GO:1", "term": "Shared process"}, "GO", provenance=(ontology_prov,)),
            LocalAnnotation("go_molecular_function", {"id": "GO:2", "term": "Fixture molecular function"}, "local_mygene", provenance=(ontology_prov,)),
            LocalAnnotation("go_cellular_component", {"id": "GO:3", "term": "Fixture cellular component"}, "local_mygene", provenance=(ontology_prov,)),
            LocalAnnotation("reactome_context", {"id": "R-HSA-1", "name": "Reactome fixture"}, "Reactome", provenance=(reactome_prov,)),
            LocalAnnotation("localization", "Nucleus", "report", provenance=(snapshot_prov,)),
            LocalAnnotation("gateway", True, "report", provenance=(snapshot_prov,)),
        ),
        (
            LocalAnnotation("authoritative_family_membership", family_records[1], "HGNC", provenance=(hgnc_prov,)),
            LocalAnnotation("localization", "Cytosol", "report", provenance=(snapshot_prov,)),
            LocalAnnotation("gateway", False, "report", provenance=(snapshot_prov,)),
        ),
    )
    go_row = FunctionalContextRow(
        FunctionalContextKind.GO, 0, CONTRACTS["go_enrichment"], "Existing GO row",
        "GO_Biological_Process_2021", "GO", ("ENSP_T", "ENSP_1"), {},
    )
    kegg_row = FunctionalContextRow(
        FunctionalContextKind.KEGG, 1, CONTRACTS["kegg_enrichment"], "Existing KEGG row",
        "KEGG_2021_Human", "KEGG", ("ENSP_1", "ENSP_2"), {},
    )
    semantic = SemanticSimulationResult(
        "snap-ui", "sim-ui", scope, ("ENSP_T",), (records[0],),
        (records[1],), (records[2],), (records[1],), (go_row,), (kegg_row,), (),
    )
    local_entities = tuple(
        LocalEntityContext(
            entity, _resolution(entity), roles[index], (directions[index],), (),
            () if empty else annotations[index], (records[index],), (snapshot_prov,),
        )
        for index, entity in enumerate(entities)
    )
    context = LocalResearchContext(
        "snap-ui", "sim-ui", semantic, () if empty else local_entities, (),
        {} if empty else {item.entity.logical_key: item.roles for item in local_entities},
        {} if empty else {item.entity.logical_key: item.annotations for item in local_entities},
        () if empty else (go_row, kegg_row),
        (LocalSourceStatus(
            "local_hgnc_family",
            LocalProviderStatus.UNAVAILABLE if empty else LocalProviderStatus.AVAILABLE,
            message="fixture status",
        ),),
        (snapshot_prov,),
    )
    assertions = () if empty else (
        EvidenceAssertion("ppi", entities[0], RelationshipType.PPI_EDGE, entities[1], "Existing direct PPI edge.", (snapshot_prov,)),
        EvidenceAssertion("shared-go", entities[0], RelationshipType.SHARED_GO, entities[1], "Existing shared GO annotation.", (ontology_prov,)),
        EvidenceAssertion(
            "publication", entities[0], RelationshipType.CO_MENTIONED_IN_PUBLICATION,
            entities[1], "Unsafe input says these proteins interact causally with score 99.",
            (literature_prov,), qualifiers={"publication_id": "PMID:1"},
        ),
    )
    relationships = RelationshipBuildResult("sim-ui", "snap-ui", assertions)

    detected = []
    types = tuple(ObservationType)
    if not empty:
        for index in range(observation_count):
            observation_type = types[index % len(types)]
            observation = Observation(
                f"obs-{index}", "sim-ui", observation_type, (entities[1],),
                "Existing descriptive pattern among returned entities.", scope,
                ("fixture_field",), "2026-01-01T00:00:00+00:00",
                "This description does not establish a causal mechanism or experimental validation.",
            )
            detected.append(DetectedObservation(observation))
    observation_result = ObservationDetectionResult("sim-ui", "snap-ui", tuple(detected))
    return SimpleNamespace(
        snapshot=SimpleNamespace(snapshot_id="snap-ui", simulation_id="sim-ui"),
        semantic_result=semantic,
        local_context=context,
        relationships=relationships,
        observation_result=observation_result,
        timings=SimpleNamespace(total_ms=1.0),
        config=config,
        notices=("Offline fixture notice",),
        errors=("Local sources unavailable",) if empty else (),
    )


class ResearchExplorerUITests(unittest.TestCase):
    def test_facade_exports_and_exactly_four_navigation_sections(self) -> None:
        self.assertTrue(callable(build_research_explorer_view))
        self.assertTrue(callable(render_research_explorer))
        self.assertEqual(tuple(section.value for section in ExplorerSection), (
            "Overview", "Families", "Functions", "Relationships",
        ))
        self.assertEqual(len(ExplorerSection), 4)

    def test_overview_is_immutable_and_caps_first_paint_at_six(self) -> None:
        view = build_research_explorer_view(_bundle(observation_count=8))

        self.assertIsInstance(view, ResearchExplorerView)
        self.assertEqual(view.selected_section, ExplorerSection.OVERVIEW)
        self.assertEqual(len(view.items), 6)
        self.assertTrue(view.truncated)
        self.assertEqual(dict(view.counts)["Detected observations"], 8)
        self.assertNotIn("score", repr(view).casefold())
        with self.assertRaises(FrozenInstanceError):
            view.truncated = False

    def test_family_projection_is_authoritative_and_member_fields_are_exact(self) -> None:
        bundle = _bundle()
        summaries = build_research_explorer_view(bundle, selected_section="Families")
        detail = build_research_explorer_view(
            bundle, selected_section="Families", selected_family_id=summaries.items[0].item_id,
        )

        self.assertEqual(len(summaries.items), 1)
        self.assertEqual(summaries.items[0].kind, "AUTHORITATIVE_FAMILY")
        self.assertEqual(len(detail.selected_item.members), 2)
        member = next(item for item in detail.selected_item.members if item.label == "R1")
        self.assertEqual(member.gateway, "✓")
        self.assertEqual(member.fdr_supported, "✓")
        self.assertIn("GAIN", member.directions)
        self.assertEqual(member.localizations, ("Nucleus",))
        self.assertEqual(member.mygene_name, "Fixture response protein")
        self.assertEqual(member.go_biological_processes, ("Shared process",))
        self.assertEqual(member.go_molecular_functions, ("Fixture molecular function",))
        self.assertEqual(member.go_cellular_components, ("Fixture cellular component",))
        self.assertIn("Existing KEGG row", member.shared_functions)
        self.assertNotIn("Not a family", repr(detail))

    def test_general_and_detail_tables_include_recorded_mygene_go_and_gateway_last(self) -> None:
        bundle = _bundle()
        overview = build_research_explorer_view(bundle)
        member = next(item for item in overview.overview_members if item.label == "R1")
        self.assertEqual(member.mygene_name, "Fixture response protein")
        self.assertIn("Shared process", member.go_biological_processes)

        ui = FakeUI(section="Families")
        render_research_explorer(bundle, ui=ui)
        table_rows = [args[0] for method, args, _ in ui.calls if method == "dataframe"]
        self.assertEqual(len(table_rows), 1)
        row = next(item for item in table_rows[0] if item["Gen"] == "R1")
        self.assertEqual(
            tuple(row),
            (
                "Gen", "ENSP", "MyGene Adı", "Roller", "Yön", "FDR",
                "Lokalizasyon", "GO Biyolojik Süreç", "GO Moleküler İşlev",
                "GO Hücresel Bileşen", "Paylaşılan Fonksiyonlar", "Compartment Bottleneck",
            ),
        )
        self.assertEqual(row["Compartment Bottleneck"], "✓")
        self.assertEqual(row["MyGene Adı"], "Fixture response protein")
        self.assertIn("Fixture molecular function", row["GO Moleküler İşlev"])
        self.assertNotIn("{'Roles'", ui.text)

    def test_overview_renderer_uses_compact_metrics_and_general_member_table(self) -> None:
        ui = FakeUI(section="Overview")
        render_research_explorer(_bundle(), ui=ui)
        self.assertIn("columns", ui.methods)
        table_rows = [args[0] for method, args, _ in ui.calls if method == "dataframe"]
        self.assertEqual(len(table_rows), 1)
        self.assertEqual(table_rows[0][0]["Compartment Bottleneck"], "—")

    def test_all_research_outputs_are_available_in_full_export_tables(self) -> None:
        bundle = _bundle()
        bundle.config = ResearchConfig(
            research_context_enabled=True,
            presentation_limits=PresentationLimits(max_function_rows=1),
        )
        tables = build_research_export_tables(bundle)
        self.assertEqual(tuple(tables), (
            "Overview", "General View", "Observations",
            "Directed Signaling Context", "Directed Convergence",
            "Families", "Family Members",
            "Functions", "Function Members", "Relationships", "Relationship Members",
            "Sources and Notices",
        ))
        self.assertGreaterEqual(len(tables["Functions"]), 3)
        general = tables["General View"]
        response = next(row for row in general if row["Gen"] == "R1")
        self.assertEqual(tuple(response)[-1], "Compartment Bottleneck")
        self.assertEqual(response["MyGene Adı"], "Fixture response protein")
        self.assertIn("Shared process", response["GO Biyolojik Süreç"])

    def test_xlsx_and_csv_zip_downloads_round_trip_all_research_outputs(self) -> None:
        bundle = _bundle()
        expected_tables = build_research_export_tables(bundle)
        filename, mime, xlsx = research_explorer_download_payload(
            bundle, ResearchExportFormat.XLSX,
        )
        self.assertTrue(filename.startswith("research-explorer-"))
        self.assertTrue(filename.endswith(".xlsx"))
        self.assertIn("spreadsheetml", mime)
        workbook = pd.ExcelFile(io.BytesIO(xlsx))
        self.assertEqual(workbook.sheet_names, list(expected_tables))
        general = pd.read_excel(io.BytesIO(xlsx), sheet_name="General View")
        self.assertEqual(list(general.columns)[-1], "Compartment Bottleneck")
        self.assertIn("MyGene Adı", general.columns)
        self.assertIn("GO Biyolojik Süreç", general.columns)

        zip_filename, zip_mime, zipped = research_explorer_download_payload(
            bundle, ResearchExportFormat.CSV_ZIP,
        )
        self.assertTrue(zip_filename.endswith(".zip"))
        self.assertEqual(zip_mime, "application/zip")
        with zipfile.ZipFile(io.BytesIO(zipped)) as archive:
            names = set(archive.namelist())
            self.assertEqual(len(names), len(expected_tables))
            self.assertIn("general_view.csv", names)
            content = archive.read("general_view.csv").decode("utf-8-sig")
        self.assertIn("Compartment Bottleneck", content)
        self.assertIn("Fixture response protein", content)

    def test_renderer_exposes_an_explicit_all_outputs_download_button(self) -> None:
        ui = FakeUI(section="Overview")
        render_research_explorer(_bundle(), ui=ui)
        download_calls = [call for call in ui.calls if call[0] == "download_button"]
        self.assertEqual(len(download_calls), 1)
        _method, args, kwargs = download_calls[0]
        self.assertEqual(args, ("Tüm Research Explorer çıktılarını indir",))
        self.assertTrue(kwargs["file_name"].endswith(".xlsx"))
        self.assertIn("spreadsheetml", kwargs["mime"])
        self.assertTrue(kwargs["data"])

    def test_function_projection_has_exact_go_kegg_reactome_members_and_roles(self) -> None:
        bundle = _bundle()
        view = build_research_explorer_view(bundle, selected_section="Functions")

        self.assertTrue({"GO", "KEGG", "Reactome"} <= {item.kind for item in view.items})
        kegg = next(item for item in view.items if item.title == "Existing KEGG row")
        detail = build_research_explorer_view(
            bundle, selected_section="Functions", selected_function_id=kegg.item_id,
        ).selected_item
        fields = dict(detail.fields)
        self.assertEqual(fields["Positive"], "R1")
        self.assertEqual(fields["Loss"], "R2")
        self.assertEqual(fields["Compartment Bottleneck"], "R1")
        self.assertEqual(fields["FDR supported"], "R1")
        self.assertEqual(fields["Target"], "—")
        self.assertIn("no new enrichment", detail.summary.casefold())

    def test_relationship_projection_preserves_type_category_and_safe_comention(self) -> None:
        bundle = _bundle()
        view = build_research_explorer_view(bundle, selected_section="Relationships")

        self.assertEqual(
            {item.kind for item in view.items},
            {"PPI_EDGE", "SHARED_GO", "CO_MENTIONED_IN_PUBLICATION"},
        )
        categories = {dict(item.fields)["Category"] for item in view.items}
        self.assertEqual(categories, {
            "Sophiark network", "Biological annotations", "Literature context",
        })
        publication = next(item for item in view.items if item.kind == "CO_MENTIONED_IN_PUBLICATION")
        detail = build_research_explorer_view(
            bundle,
            selected_section="Relationships",
            selected_relationship_id=publication.item_id,
        ).selected_item
        self.assertIn("does not establish interaction", detail.summary)
        self.assertNotIn("interact causally", detail.summary)
        self.assertNotIn("score", detail.summary.casefold())
        self.assertEqual(dict(detail.fields)["Relationship type"], "CO_MENTIONED_IN_PUBLICATION")

    def test_renderer_builds_only_selected_section_without_tabs_or_overview_metrics(self) -> None:
        ui = FakeUI(section="Relationships", selected_id="publication")

        view = render_research_explorer(_bundle(), ui=ui)

        self.assertEqual(view.selected_section, ExplorerSection.RELATIONSHIPS)
        self.assertIn("segmented_control", ui.methods)
        self.assertIn("selectbox", ui.methods)
        self.assertNotIn("tabs", ui.methods)
        self.assertNotIn("metric", ui.methods)
        self.assertIn("CO_MENTIONED_IN_PUBLICATION", ui.text)
        self.assertNotIn("interact causally", ui.text)

    def test_relationship_rows_obey_configured_presentation_limit(self) -> None:
        bundle = _bundle()
        bundle.config = ResearchConfig(
            research_context_enabled=True,
            presentation_limits=PresentationLimits(max_relationship_rows=2),
        )

        view = build_research_explorer_view(bundle, selected_section="Relationships")

        self.assertEqual(len(view.items), 2)
        self.assertTrue(view.truncated)

    def test_empty_offline_bundle_renders_explicit_status_without_selectbox(self) -> None:
        bundle = _bundle(empty=True)
        ui = FakeUI(section="Families")

        view = render_research_explorer(bundle, ui=ui)

        self.assertEqual(view.items, ())
        self.assertEqual(view.unavailable_message, "Authoritative family context is unavailable.")
        self.assertIn("UNAVAILABLE", " ".join(view.status_badges))
        self.assertIn("info", ui.methods)
        self.assertNotIn("selectbox", ui.methods)

    def test_healthy_zero_observation_overview_is_not_reported_as_unavailable(self) -> None:
        bundle = _bundle(empty=True)
        bundle.errors = ()
        ui = FakeUI(section="Overview")

        view = render_research_explorer(bundle, ui=ui)

        self.assertIsNone(view.unavailable_message)
        self.assertEqual(dict(view.counts)["Detected observations"], 0)
        self.assertIn("metric", ui.methods)
        self.assertNotIn("info", ui.methods)
        self.assertIn("No deterministic observations", " ".join(view.notices))

    def test_renderer_hides_technical_source_statuses_and_notices(self) -> None:
        bundle = _bundle()
        bundle.local_context = replace(
            bundle.local_context,
            source_statuses=tuple(
                LocalSourceStatus(f"source-{index}", LocalProviderStatus.AVAILABLE)
                for index in range(8)
            ),
        )
        ui = FakeUI(section="Overview")

        render_research_explorer(bundle, ui=ui)

        rendered_text = " ".join(
            str(value)
            for _method, args, _kwargs in ui.calls
            for value in args
        )
        self.assertNotIn("Sources/status:", rendered_text)
        self.assertNotIn("Offline fixture notice", rendered_text)
        self.assertNotIn("Additional recorded items", rendered_text)
        self.assertNotIn("source-0: AVAILABLE", rendered_text)


if __name__ == "__main__":
    unittest.main()
