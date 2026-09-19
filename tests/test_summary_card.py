from src.ui.summary_card import analysis_summary_html, summary_cell_html


def test_summary_markup_is_compact_html_not_indented_markdown_code():
    cell = summary_cell_html("Directed aday", "20")
    markup = analysis_summary_html(
        title="Analiz Özeti", target_html="ANO2", context_label="Lung · DIRECTED",
        cells=(cell,),
    )
    assert markup.startswith('<div class="sk-readout-wrap">')
    assert "\n" not in markup
    assert '<div class="sk-readout-grid"><div class="sk-readout-cell">' in markup
    assert "Directed aday" in markup and ">20<" in markup


def test_directed_beta_labels_are_visible_in_engine_selector_and_result_heading():
    source = open("app.py", encoding="utf-8").read()
    assert '"Directed · BETA" if mode == "Directed"' in source
    assert 'ui_t(\'directed_calculation_engine\')' in source
