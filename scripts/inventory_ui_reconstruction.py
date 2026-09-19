"""Record source-level UI coverage and immutable pre-redesign source evidence."""
from pathlib import Path
import ast
import hashlib
import json
import shutil

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/ui_reconstruction"
KINDS = set("title header subheader markdown caption write text info warning error success exception metric dataframe table json code tabs expander container columns selectbox multiselect radio segmented_control slider select_slider number_input text_input text_area checkbox toggle button download_button file_uploader form form_submit_button plotly_chart pyplot altair_chart bar_chart line_chart image iframe status spinner progress".split())


def main():
    OUT.mkdir(exist_ok=True, parents=True)
    paths = [ROOT / "app.py", *sorted((ROOT / "src").rglob("*.py"))]
    hashes = {}
    items = []
    for path in paths:
        relative = path.relative_to(ROOT).as_posix()
        hashes[relative] = hashlib.sha256(path.read_bytes()).hexdigest()
        if not (relative == "app.py" or relative.startswith("src/ui/") or relative.endswith("/ui.py")):
            continue
        source = path.read_text(encoding="utf-8-sig")
        tree = ast.parse(source)
        for node in ast.walk(tree):
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and node.func.attr in KINDS:
                receiver = ast.unparse(node.func.value)
                if receiver in {"pd", "log", "logging"} or receiver.endswith("session_state"):
                    continue
                expression = ast.get_source_segment(source, node) or ast.unparse(node)
                items.append(dict(id=f"UI-{len(items)+1:04}", old_file=relative, old_line=node.lineno,
                    kind=node.func.attr, old_element=expression, current_purpose=node.func.attr,
                    new_location=relative, preserved="PENDING", functionally_tested="UNVERIFIED",
                    notes="Internal page (not registered by app.py)" if '/pages/0' in relative else "Includes conditional branches; source inventory is not runtime acceptance"))
    # Never replace the initial baseline when the script is rerun.
    for name, value in [("source-baseline.json", hashes), ("coverage-manifest.json", items)]:
        target = OUT / name
        if not target.exists():
            target.write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding="utf-8")
    for path in [ROOT / "app.py", ROOT / ".streamlit/config.toml", *sorted((ROOT / "src/ui").rglob("*.py")), ROOT / "src/state.py"]:
        target = OUT / "before" / path.relative_to(ROOT)
        target.parent.mkdir(exist_ok=True, parents=True)
        if not target.exists():
            shutil.copy2(path, target)
    print(f"Inventoried {len(items)} UI calls; recorded {len(hashes)} source hashes.")


if __name__ == "__main__":
    main()
