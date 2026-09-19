"""Local process boundary for the unchanged, CPU-intensive dose engine.

Only application-created local payloads are loaded. The interpreter, working
directory, graph, scores and engine arguments are the current local runtime's.
"""
from __future__ import annotations

import importlib
import json
import os
from pathlib import Path
import pickle
import random
import subprocess
import sys
import time
import traceback
import uuid

import numpy as np

ROOT = Path(__file__).resolve().parents[2]


def start_dose_job(*, motor_name, graph, scores, output_dir, targets, hill_n, paralog_boost, damping):
    if motor_name not in {"src.main", "src.mouse.main"}:
        raise ValueError("Unsupported dose engine")
    folder = ROOT / "outputs" / "ui_jobs" / uuid.uuid4().hex
    folder.mkdir(parents=True)
    payload = dict(motor_name=motor_name, graph=graph, scores=scores,
        output_dir=output_dir, targets=targets, hill_n=hill_n,
        paralog_boost=paralog_boost, damping=damping,
        random_state=random.getstate(), numpy_state=np.random.get_state())
    with (folder / "input.pkl").open("wb") as stream:
        pickle.dump(payload, stream, protocol=pickle.HIGHEST_PROTOCOL)
    with (folder / "worker.log").open("wb") as log:
        process = subprocess.Popen([sys.executable, "-m", "src.ui.dose_job", str(folder)],
            cwd=ROOT, stdout=log, stderr=subprocess.STDOUT,
            creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0)
    return dict(folder=str(folder), started=time.time(), process=process)


def _work(folder):
    folder = Path(folder).resolve()
    if folder.parent != (ROOT / "outputs" / "ui_jobs").resolve():
        raise ValueError("Job outside local UI job directory")
    try:
        with (folder / "input.pkl").open("rb") as stream:
            payload = pickle.load(stream)
        motor = importlib.import_module(payload["motor_name"])
        random.setstate(payload["random_state"])
        np.random.set_state(payload["numpy_state"])
        result = motor.run_pharmacological_dose_response(
            payload["graph"], payload["scores"], payload["output_dir"],
            spesifik_hedefler=payload["targets"], hill_n=payload["hill_n"],
            paralog_boost=payload["paralog_boost"], damping=payload["damping"])
        with (folder / "result.tmp").open("wb") as stream:
            pickle.dump(result, stream, protocol=pickle.HIGHEST_PROTOCOL)
        (folder / "result.tmp").replace(folder / "result.pkl")
    except Exception:
        (folder / "error.json").write_text(json.dumps({"error": traceback.format_exc()}), encoding="utf-8")
        raise


if __name__ == "__main__":
    _work(sys.argv[1])
