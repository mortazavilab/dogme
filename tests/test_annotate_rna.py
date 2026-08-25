import importlib.util
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "annotateRNA.py"
module_spec = importlib.util.spec_from_file_location("annotateRNA", SCRIPT)
annotate_rna = importlib.util.module_from_spec(module_spec)
module_spec.loader.exec_module(annotate_rna)


def test_filter_novel_models_excludes_antisense_transcripts_only_when_disabled():
    novel_models = {
        "known": {"id": "known_tx", "strand": "+"},
        "antisense": {"id": "antisense_tx", "strand": "-"},
    }
    abundance_counter = Counter({
        ("gene1", "known_tx", "Known"): 3,
        ("gene2", "antisense_tx", "Antisense"): 2,
    })

    assert annotate_rna.filter_novel_models_for_export(
        novel_models, abundance_counter, include_antisense=False
    ) == {"known": novel_models["known"]}
    assert annotate_rna.filter_novel_models_for_export(
        novel_models, abundance_counter, include_antisense=True
    ) == novel_models