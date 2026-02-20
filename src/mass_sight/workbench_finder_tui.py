"""Textual TUI for finding Metabolomics Workbench studies/analyses by disease.

This is intentionally lightweight:
- Start from the Workbench disease index (`/study_id/ST/disease`).
- Let users drill down disease -> study -> analyses.
- Export selected analysis IDs to a reproducible selection manifest.
"""

from __future__ import annotations

import asyncio
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from textual import on
from textual.app import App, ComposeResult
from textual.screen import Screen
from textual.widgets import DataTable, Footer, Header, Input, Static

from .workbench import (
    fetch_workbench_all_study_diseases,
    fetch_workbench_all_study_summaries,
    fetch_workbench_study_analyses,
)


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _norm(s: object) -> str:
    return str(s or "").strip()


def _canonical_study_id(raw: object) -> str:
    s = _norm(raw).upper()
    if not s:
        return ""
    if s.startswith("ST"):
        return s
    digits = "".join(ch for ch in s if ch.isdigit())
    if not digits:
        return s
    try:
        return f"ST{int(digits):06d}"
    except Exception:
        return s


def _canonical_analysis_id(raw: object) -> str:
    s = _norm(raw).upper()
    if not s:
        return ""
    if s.startswith("AN"):
        return s
    digits = "".join(ch for ch in s if ch.isdigit())
    if not digits:
        return s
    try:
        return f"AN{int(digits):06d}"
    except Exception:
        return s


def _row_key_value(key: object) -> str:
    """Extract the underlying user-provided key from Textual DataTable RowKey objects."""
    if key is None:
        return ""
    v = getattr(key, "value", None)
    if v is not None:
        return str(v)
    return str(key)


def _datatable_cursor_row_key(table: DataTable) -> object | None:
    """Return the RowKey for the current cursor position (Textual v8 compatible)."""
    try:
        if table.row_count <= 0:
            return None
        if table.cursor_row is None:
            return None
        cell_key = table.coordinate_to_cell_key(table.cursor_coordinate)
        return cell_key.row_key
    except Exception:
        return None


def _parse_int(raw: object) -> Optional[int]:
    s = _norm(raw)
    if not s:
        return None
    try:
        return int(float(s))
    except Exception:
        return None


def _samples_bucket(n_samples_raw: object) -> str:
    """Bucket sample counts to avoid a high-cardinality facet."""
    n = _parse_int(n_samples_raw)
    if n is None or n < 0:
        return "<missing>"
    if n <= 5:
        return "1-5"
    if n <= 10:
        return "6-10"
    if n <= 20:
        return "11-20"
    if n <= 50:
        return "21-50"
    if n <= 100:
        return "51-100"
    return "101+"


def _canonical_ion_mode(raw: object) -> str:
    s = _norm(raw).lower()
    if not s:
        return ""
    if "pos" in s:
        return "positive"
    if "neg" in s:
        return "negative"
    return s


def _canonical_chromatography(raw: object) -> str:
    s = _norm(raw).lower()
    if not s:
        return ""
    if "hilic" in s:
        return "hilic"
    if ("reverse" in s) or ("reversed" in s) or ("rp" in s) or ("c18" in s):
        return "reversed_phase"
    if "gc" in s:
        return "gc"
    return s


def _canonical_instrument_type(raw: object) -> str:
    s = _norm(raw)
    if not s:
        return ""
    return s.upper()


FACET_ANALYSIS_TYPE = "analysis_type"
FACET_SPECIES = "species"
FACET_INSTITUTE = "institute"
FACET_N_SAMPLES_BUCKET = "n_samples_bucket"
FACET_ION_MODE = "ion_mode"
FACET_CHROMATOGRAPHY = "chromatography"
FACET_MS_INSTRUMENT_TYPE = "ms_instrument_type"


_FACETS_ORDER: List[str] = [
    FACET_ANALYSIS_TYPE,
    FACET_SPECIES,
    FACET_N_SAMPLES_BUCKET,
    FACET_ION_MODE,
    FACET_CHROMATOGRAPHY,
    FACET_MS_INSTRUMENT_TYPE,
    FACET_INSTITUTE,
]


def _facet_label(facet: str) -> str:
    return {
        FACET_ANALYSIS_TYPE: "analysis_type (study summary)",
        FACET_SPECIES: "species (study summary)",
        FACET_N_SAMPLES_BUCKET: "n_samples (bucketed)",
        FACET_ION_MODE: "ion_mode (per-study analyses)",
        FACET_CHROMATOGRAPHY: "chromatography (per-study analyses)",
        FACET_MS_INSTRUMENT_TYPE: "ms_instrument_type (per-study analyses)",
        FACET_INSTITUTE: "institute (study summary)",
    }.get(str(facet), str(facet))


def _is_advanced_facet(facet: str) -> bool:
    return facet in {FACET_ION_MODE, FACET_CHROMATOGRAPHY, FACET_MS_INSTRUMENT_TYPE}


@dataclass(frozen=True)
class FinderConfig:
    cache_dir: Path
    refresh: bool
    timeout_s: float
    retries: int
    delay_s: float
    out_path: Path


class DiseaseScreen(Screen):
    BINDINGS = [
        ("q", "app.quit", "Quit"),
        ("ctrl+q", "app.quit", "Quit"),
        ("enter", "select", "Select"),
        ("r", "refresh", "Refresh"),
        ("/", "focus_filter", "Filter"),
        ("escape", "focus_table", "Table"),
    ]

    def __init__(self) -> None:
        super().__init__()
        self._rows: List[Tuple[str, int]] = []
        self._filter: str = ""

    def compose(self) -> ComposeResult:
        yield Header()
        yield Static("Select disease (counts are number of studies).", id="title")
        yield Input(placeholder="Filter diseases...", id="filter")
        yield DataTable(id="table")
        yield Footer()

    async def on_mount(self) -> None:
        table = self.query_one("#table", DataTable)
        table.add_columns("Disease", "n_studies")
        await self._load()
        table.focus()

    def action_focus_filter(self) -> None:
        self.query_one("#filter", Input).focus()

    def action_focus_table(self) -> None:
        self.query_one("#table", DataTable).focus()

    async def _load(self) -> None:
        await self.app.ensure_indices_loaded()  # type: ignore[attr-defined]
        mapping: Dict[str, set[str]] = self.app.disease_to_studies  # type: ignore[attr-defined]
        rows = [(d, len(sts)) for d, sts in mapping.items()]
        rows.sort(key=lambda x: (-x[1], x[0].lower()))
        self._rows = rows
        self._render_table()

    def _render_table(self) -> None:
        table = self.query_one("#table", DataTable)
        table.clear(columns=False)
        f = self._filter.lower().strip()
        for disease, n in self._rows:
            if f and f not in disease.lower():
                continue
            table.add_row(disease, str(int(n)), key=disease)

    @on(Input.Changed, "#filter")
    def _on_filter_changed(self, event: Input.Changed) -> None:
        self._filter = str(event.value or "")
        self._render_table()

    @on(DataTable.CellSelected, "#table")
    def _on_cell_selected(self, event: DataTable.CellSelected) -> None:
        disease = _row_key_value(event.cell_key.row_key)
        if not disease:
            return
        self.app.push_screen(StudyScreen(disease))  # type: ignore[attr-defined]

    def action_select(self) -> None:
        table = self.query_one("#table", DataTable)
        row_key = _datatable_cursor_row_key(table)
        if row_key is None:
            return
        self.app.push_screen(StudyScreen(_row_key_value(row_key)))  # type: ignore[attr-defined]

    async def action_refresh(self) -> None:
        await self.app.ensure_indices_loaded(force_refresh=True)  # type: ignore[attr-defined]
        await self._load()


class StudyScreen(Screen):
    BINDINGS = [
        ("q", "app.quit", "Quit"),
        ("ctrl+q", "app.quit", "Quit"),
        ("b", "app.pop_screen", "Back"),
        ("g", "facets", "Facets"),
        ("x", "clear_facets", "Clear facets"),
        ("enter", "select", "Select"),
        ("r", "refresh", "Refresh"),
        ("/", "focus_filter", "Filter"),
        ("escape", "focus_table", "Table"),
    ]

    def __init__(self, disease: str) -> None:
        super().__init__()
        self.disease = str(disease)
        self._rows: List[Dict[str, str]] = []
        self._filter: str = ""

    def compose(self) -> ComposeResult:
        yield Header()
        yield Static(f"Disease: {self.disease}", id="title")
        yield Input(placeholder="Filter studies (id/title/species/institute/type/samples)...", id="filter")
        yield DataTable(id="table")
        yield Static("", id="status")
        yield Footer()

    async def on_mount(self) -> None:
        table = self.query_one("#table", DataTable)
        table.add_columns("study_id", "analysis_type", "n_samples", "study_title", "species", "institute")
        await self._load()
        table.focus()

    def action_focus_filter(self) -> None:
        self.query_one("#filter", Input).focus()

    def action_focus_table(self) -> None:
        self.query_one("#table", DataTable).focus()

    async def _load(self) -> None:
        await self.app.ensure_indices_loaded()  # type: ignore[attr-defined]
        self._rows = self.app.get_disease_study_rows(self.disease)  # type: ignore[attr-defined]
        self._render_table()

    def _render_table(self) -> None:
        table = self.query_one("#table", DataTable)
        table.clear(columns=False)
        status = self.query_one("#status", Static)
        f = self._filter.lower().strip()
        filters = self.app.get_facet_filters(self.disease)  # type: ignore[attr-defined]
        n_total = 0
        n_shown = 0
        for r in self._rows:
            n_total += 1
            if not self.app.study_row_passes_filters(self.disease, r):  # type: ignore[attr-defined]
                continue
            blob = " ".join(
                [
                    r.get("study_id", ""),
                    r.get("study_title", ""),
                    r.get("analysis_type", ""),
                    r.get("n_samples", ""),
                    r.get("species", ""),
                    r.get("institute", ""),
                ]
            ).lower()
            if f and f not in blob:
                continue
            sid = r.get("study_id", "")
            table.add_row(
                r.get("study_id", ""),
                r.get("analysis_type", ""),
                r.get("n_samples", ""),
                r.get("study_title", ""),
                r.get("species", ""),
                r.get("institute", ""),
                key=sid,
            )
            n_shown += 1

        facet_summary = self.app.format_facet_filters_summary(self.disease, filters=filters)  # type: ignore[attr-defined]
        status.update(f"Showing {n_shown}/{n_total} studies. {facet_summary}")

    @on(Input.Changed, "#filter")
    def _on_filter_changed(self, event: Input.Changed) -> None:
        self._filter = str(event.value or "")
        self._render_table()

    @on(DataTable.CellSelected, "#table")
    def _on_cell_selected(self, event: DataTable.CellSelected) -> None:
        sid = _canonical_study_id(_row_key_value(event.cell_key.row_key))
        if not sid:
            return
        self.app.push_screen(AnalysisScreen(study_id=sid, disease=self.disease))  # type: ignore[attr-defined]

    def action_select(self) -> None:
        table = self.query_one("#table", DataTable)
        row_key = _datatable_cursor_row_key(table)
        if row_key is None:
            return
        sid = _canonical_study_id(_row_key_value(row_key))
        self.app.push_screen(AnalysisScreen(study_id=sid, disease=self.disease))  # type: ignore[attr-defined]

    async def action_refresh(self) -> None:
        await self.app.ensure_indices_loaded(force_refresh=True)  # type: ignore[attr-defined]
        await self._load()

    def action_facets(self) -> None:
        # Re-render on return, since facet selections live in app state.
        self.app.push_screen(DiseaseFacetsScreen(self.disease), callback=lambda _: self._render_table())  # type: ignore[attr-defined]

    def action_clear_facets(self) -> None:
        self.app.clear_facet_filters(self.disease)  # type: ignore[attr-defined]
        self._render_table()


class DiseaseFacetsScreen(Screen):
    """Facet filter UI scoped to a single disease."""

    BINDINGS = [
        ("q", "app.quit", "Quit"),
        ("ctrl+q", "app.quit", "Quit"),
        ("b", "app.pop_screen", "Back"),
        ("tab", "next_facet", "Next facet"),
        ("shift+tab", "prev_facet", "Prev facet"),
        ("l", "load_advanced", "Load advanced"),
        ("space", "toggle", "Toggle"),
        ("enter", "toggle", "Toggle"),
        ("c", "clear", "Clear facet"),
        ("x", "clear_all", "Clear all"),
        ("/", "focus_filter", "Filter"),
        ("escape", "focus_table", "Table"),
    ]

    def __init__(self, disease: str) -> None:
        super().__init__()
        self.disease = str(disease)
        self._facet_idx = 0
        self._filter: str = ""
        self._rows: List[Tuple[str, int]] = []
        self._advanced_missing: int = 0
        self._advanced_task: Optional[asyncio.Task] = None

    @property
    def facet(self) -> str:
        try:
            return _FACETS_ORDER[int(self._facet_idx) % len(_FACETS_ORDER)]
        except Exception:
            return FACET_ANALYSIS_TYPE

    def compose(self) -> ComposeResult:
        yield Header()
        yield Static("", id="title")
        yield Input(placeholder="Filter facet values...", id="filter")
        yield DataTable(id="table")
        yield Static("", id="status")
        yield Footer()

    async def on_mount(self) -> None:
        table = self.query_one("#table", DataTable)
        table.add_columns("sel", "value", "n_studies")
        await self._load()
        table.focus()

    def on_unmount(self) -> None:
        if self._advanced_task is not None and not self._advanced_task.done():
            try:
                self._advanced_task.cancel()
            except Exception:
                pass

    def action_focus_filter(self) -> None:
        self.query_one("#filter", Input).focus()

    def action_focus_table(self) -> None:
        self.query_one("#table", DataTable).focus()

    def _update_title(self) -> None:
        title = self.query_one("#title", Static)
        title.update(
            f"Disease: {self.disease} | Facet: {_facet_label(self.facet)}  "
            f"(tab/shift+tab to change, space/enter to toggle)"
        )

    async def _load(self) -> None:
        self._update_title()
        status = self.query_one("#status", Static)
        status.update("")

        await self.app.ensure_indices_loaded()  # type: ignore[attr-defined]

        if _is_advanced_facet(self.facet):
            # Load advanced facet metadata only for studies that pass all other filters.
            target_ids = [
                _canonical_study_id(r.get("study_id"))
                for r in (self.app.get_disease_study_rows(self.disease) or [])  # type: ignore[attr-defined]
                if self.app.study_row_passes_filters(self.disease, r, exclude_facet=self.facet)  # type: ignore[attr-defined]
            ]
            target_ids = [x for x in target_ids if x]
            missing = [sid for sid in target_ids if sid not in self.app.study_analysis_facets_by_id]  # type: ignore[attr-defined]
            self._advanced_missing = int(len(missing))

            if self._advanced_missing > 0:
                table = self.query_one("#table", DataTable)
                table.clear(columns=False)
                table.add_row("", "<press 'l' to load per-study analyses>", str(self._advanced_missing), key="__load_required__")
                status.update(
                    f"Advanced facet requires loading per-study analyses for {self._advanced_missing} studies "
                    f"(current filter). Press 'l' to load. Tip: narrow via analysis_type/species/n_samples first."
                )
                return

        self._advanced_missing = 0
        self._rows = self._compute_facet_counts()
        self._render_table()

        filters = self.app.get_facet_filters(self.disease)  # type: ignore[attr-defined]
        status.update(self.app.format_facet_filters_summary(self.disease, filters=filters))  # type: ignore[attr-defined]

    def _facet_values_for_study(self, *, study_id: str, row: Dict[str, str]) -> List[str]:
        facet = self.facet
        if facet == FACET_ANALYSIS_TYPE:
            return [_norm(row.get("analysis_type")) or "<missing>"]
        if facet == FACET_SPECIES:
            return [_norm(row.get("species")) or "<missing>"]
        if facet == FACET_INSTITUTE:
            return [_norm(row.get("institute")) or "<missing>"]
        if facet == FACET_N_SAMPLES_BUCKET:
            return [_samples_bucket(row.get("n_samples"))]
        if facet in {FACET_ION_MODE, FACET_CHROMATOGRAPHY, FACET_MS_INSTRUMENT_TYPE}:
            facets_by_study = self.app.study_analysis_facets_by_id.get(study_id) or {}  # type: ignore[attr-defined]
            values = facets_by_study.get(facet) or set()
            if not values:
                return ["<missing>"]
            return sorted(str(v) for v in values if str(v))
        return ["<missing>"]

    def _compute_facet_counts(self) -> List[Tuple[str, int]]:
        app = self.app  # type: ignore[attr-defined]
        rows = app.get_disease_study_rows(self.disease)
        counts: Dict[str, int] = {}

        # Facet counts are computed on the subset passing all other facet filters.
        for r in rows:
            if not app.study_row_passes_filters(self.disease, r, exclude_facet=self.facet):
                continue
            sid = _canonical_study_id(r.get("study_id"))
            for v in self._facet_values_for_study(study_id=sid, row=r):
                counts[v] = counts.get(v, 0) + 1

        items = list(counts.items())
        items.sort(key=lambda x: (-x[1], x[0].lower()))
        return items

    def _render_table(self) -> None:
        table = self.query_one("#table", DataTable)
        table.clear(columns=False)
        selected = self.app.get_facet_filters(self.disease).get(self.facet, set())  # type: ignore[attr-defined]
        f = self._filter.lower().strip()
        for value, n in self._rows:
            if f and f not in value.lower():
                continue
            mark = "x" if value in selected else ""
            table.add_row(mark, value, str(int(n)), key=value)

    @on(Input.Changed, "#filter")
    def _on_filter_changed(self, event: Input.Changed) -> None:
        self._filter = str(event.value or "")
        self._render_table()

    @on(DataTable.CellSelected, "#table")
    def _on_cell_selected(self, event: DataTable.CellSelected) -> None:
        # DataTable consumes Enter; treat cell selection as a toggle.
        if self._advanced_missing > 0:
            return
        value = _row_key_value(event.cell_key.row_key)
        if not value:
            return
        if value == "__load_required__":
            return
        self.app.toggle_facet_filter(self.disease, self.facet, value)  # type: ignore[attr-defined]
        self._rows = self._compute_facet_counts()
        self._render_table()
        status = self.query_one("#status", Static)
        filters = self.app.get_facet_filters(self.disease)  # type: ignore[attr-defined]
        status.update(self.app.format_facet_filters_summary(self.disease, filters=filters))  # type: ignore[attr-defined]

    async def action_next_facet(self) -> None:
        self._facet_idx = int(self._facet_idx) + 1
        await self._load()

    async def action_prev_facet(self) -> None:
        self._facet_idx = int(self._facet_idx) - 1
        await self._load()

    def action_toggle(self) -> None:
        if self._advanced_missing > 0:
            return
        table = self.query_one("#table", DataTable)
        row_key = _datatable_cursor_row_key(table)
        if row_key is None:
            return
        value = _row_key_value(row_key)
        if not value:
            return
        if value == "__load_required__":
            return
        self.app.toggle_facet_filter(self.disease, self.facet, value)  # type: ignore[attr-defined]
        self._rows = self._compute_facet_counts()
        self._render_table()
        status = self.query_one("#status", Static)
        filters = self.app.get_facet_filters(self.disease)  # type: ignore[attr-defined]
        status.update(self.app.format_facet_filters_summary(self.disease, filters=filters))  # type: ignore[attr-defined]

    def action_clear(self) -> None:
        self.app.clear_facet_filters(self.disease, facet=self.facet)  # type: ignore[attr-defined]
        self._rows = self._compute_facet_counts()
        self._render_table()
        status = self.query_one("#status", Static)
        filters = self.app.get_facet_filters(self.disease)  # type: ignore[attr-defined]
        status.update(self.app.format_facet_filters_summary(self.disease, filters=filters))  # type: ignore[attr-defined]

    def action_clear_all(self) -> None:
        self.app.clear_facet_filters(self.disease)  # type: ignore[attr-defined]
        self._rows = self._compute_facet_counts()
        self._render_table()
        status = self.query_one("#status", Static)
        filters = self.app.get_facet_filters(self.disease)  # type: ignore[attr-defined]
        status.update(self.app.format_facet_filters_summary(self.disease, filters=filters))  # type: ignore[attr-defined]

    def _safe_status_update(self, message: str) -> None:
        try:
            self.query_one("#status", Static).update(str(message))
        except Exception:
            return

    async def _run_advanced_load(self, target_ids: List[str]) -> None:
        try:
            await self.app.ensure_disease_advanced_facets_loaded(  # type: ignore[attr-defined]
                self.disease,
                study_ids=target_ids,
                progress_cb=lambda msg: self._safe_status_update(msg),
            )
        except asyncio.CancelledError:
            return
        except Exception as exc:
            self._safe_status_update(f"Advanced facet load failed: {exc}")
            return
        try:
            await self._load()
        except Exception:
            return

    def action_load_advanced(self) -> None:
        if not _is_advanced_facet(self.facet):
            self._safe_status_update("Advanced loading is only needed for ion_mode/chrom/instrument facets.")
            return

        target_ids = [
            _canonical_study_id(r.get("study_id"))
            for r in (self.app.get_disease_study_rows(self.disease) or [])  # type: ignore[attr-defined]
            if self.app.study_row_passes_filters(self.disease, r, exclude_facet=self.facet)  # type: ignore[attr-defined]
        ]
        target_ids = [x for x in target_ids if x]
        missing = [sid for sid in target_ids if sid not in self.app.study_analysis_facets_by_id]  # type: ignore[attr-defined]
        if not missing:
            self._safe_status_update("Advanced facet metadata already loaded for the current filter.")
            self._advanced_missing = 0
            return

        if self._advanced_task is not None and not self._advanced_task.done():
            return

        self._safe_status_update(
            f"Loading per-study analyses for {len(missing)} studies (cached). This may take a while..."
        )
        self._advanced_task = asyncio.create_task(self._run_advanced_load(target_ids))


class AnalysisScreen(Screen):
    BINDINGS = [
        ("q", "app.quit", "Quit"),
        ("ctrl+q", "app.quit", "Quit"),
        ("b", "app.pop_screen", "Back"),
        ("space", "toggle", "Toggle"),
        ("enter", "toggle", "Toggle"),
        ("s", "save", "Save"),
        ("r", "refresh", "Refresh"),
        ("/", "focus_filter", "Filter"),
        ("escape", "focus_table", "Table"),
    ]

    def __init__(self, *, study_id: str, disease: str) -> None:
        super().__init__()
        self.study_id = _canonical_study_id(study_id)
        self.disease = str(disease)
        self._rows: List[Dict[str, str]] = []
        self._filter: str = ""

    def compose(self) -> ComposeResult:
        yield Header()
        yield Static(f"Study: {self.study_id} (Disease: {self.disease})", id="title")
        yield Input(placeholder="Filter analyses (id/type/ion_mode/chrom/instrument/summary)...", id="filter")
        yield DataTable(id="table")
        yield Static("", id="status")
        yield Footer()

    async def on_mount(self) -> None:
        table = self.query_one("#table", DataTable)
        table.add_columns(
            "sel",
            "analysis_id",
            "analysis_type",
            "ion_mode",
            "chromatography_type",
            "ms_instrument_type",
            "analysis_summary",
        )
        await self._load()
        table.focus()

    def action_focus_filter(self) -> None:
        self.query_one("#filter", Input).focus()

    def action_focus_table(self) -> None:
        self.query_one("#table", DataTable).focus()

    async def _load(self) -> None:
        cfg: FinderConfig = self.app.cfg  # type: ignore[attr-defined]
        status = self.query_one("#status", Static)
        status.update("Loading analyses...")
        obj = await asyncio.to_thread(
            fetch_workbench_study_analyses,
            self.study_id,
            cache_root=cfg.cache_dir,
            timeout_s=float(cfg.timeout_s),
            retries=int(cfg.retries),
            delay_s=float(cfg.delay_s),
            refresh=bool(cfg.refresh),
        )
        if not isinstance(obj, dict):
            status.update("Failed to load analyses for this study.")
            self._rows = []
            self._render_table()
            return

        rows: List[Dict[str, str]] = []
        for rec in obj.values():
            if not isinstance(rec, dict):
                continue
            aid = _canonical_analysis_id(rec.get("analysis_id"))
            if not aid:
                continue
            rows.append(
                {
                    "analysis_id": aid,
                    "analysis_summary": _norm(rec.get("analysis_summary") or ""),
                    "analysis_type": _norm(rec.get("analysis_type") or ""),
                    "ion_mode": _norm(rec.get("ion_mode") or ""),
                    "chromatography_type": _norm(rec.get("chromatography_type") or rec.get("chromatography system") or ""),
                    "ms_instrument_type": _norm(rec.get("ms_instrument_type") or ""),
                }
            )

        rows.sort(key=lambda r: r["analysis_id"])
        self._rows = rows

        # Update analysis_id -> study_id mapping (used when exporting selection).
        self.app.analysis_to_study.update({r["analysis_id"]: self.study_id for r in rows})  # type: ignore[attr-defined]

        status.update(f"{len(rows)} analyses loaded. Selected: {len(self.app.selected_analyses)}")  # type: ignore[attr-defined]
        self._render_table()

    def _render_table(self) -> None:
        table = self.query_one("#table", DataTable)
        table.clear(columns=False)
        selected: set[str] = self.app.selected_analyses  # type: ignore[attr-defined]
        f = self._filter.lower().strip()
        for r in self._rows:
            blob = " ".join(
                [
                    r.get("analysis_id", ""),
                    r.get("analysis_type", ""),
                    r.get("ion_mode", ""),
                    r.get("chromatography_type", ""),
                    r.get("ms_instrument_type", ""),
                    r.get("analysis_summary", ""),
                ]
            ).lower()
            if f and f not in blob:
                continue
            aid = r.get("analysis_id", "")
            mark = "x" if aid in selected else ""
            table.add_row(
                mark,
                aid,
                r.get("analysis_type", ""),
                r.get("ion_mode", ""),
                r.get("chromatography_type", ""),
                r.get("ms_instrument_type", ""),
                r.get("analysis_summary", ""),
                key=aid,
            )

    @on(Input.Changed, "#filter")
    def _on_filter_changed(self, event: Input.Changed) -> None:
        self._filter = str(event.value or "")
        self._render_table()

    @on(DataTable.CellSelected, "#table")
    def _on_cell_selected(self, event: DataTable.CellSelected) -> None:
        # DataTable consumes Enter; treat cell selection as a toggle event.
        aid = _canonical_analysis_id(_row_key_value(event.cell_key.row_key))
        if not aid:
            return
        selected: set[str] = self.app.selected_analyses  # type: ignore[attr-defined]
        if aid in selected:
            selected.remove(aid)
        else:
            selected.add(aid)
        status = self.query_one("#status", Static)
        status.update(f"{len(self._rows)} analyses loaded. Selected: {len(selected)}")
        self._render_table()

    def action_toggle(self) -> None:
        table = self.query_one("#table", DataTable)
        row_key = _datatable_cursor_row_key(table)
        if row_key is None:
            return
        aid = _canonical_analysis_id(_row_key_value(row_key))
        if not aid:
            return
        selected: set[str] = self.app.selected_analyses  # type: ignore[attr-defined]
        if aid in selected:
            selected.remove(aid)
        else:
            selected.add(aid)
        status = self.query_one("#status", Static)
        status.update(f"{len(self._rows)} analyses loaded. Selected: {len(selected)}")
        self._render_table()

    async def action_refresh(self) -> None:
        self.app.cfg = FinderConfig(  # type: ignore[attr-defined]
            cache_dir=self.app.cfg.cache_dir,  # type: ignore[attr-defined]
            refresh=True,
            timeout_s=self.app.cfg.timeout_s,  # type: ignore[attr-defined]
            retries=self.app.cfg.retries,  # type: ignore[attr-defined]
            delay_s=self.app.cfg.delay_s,  # type: ignore[attr-defined]
            out_path=self.app.cfg.out_path,  # type: ignore[attr-defined]
        )
        await self._load()

    def action_save(self) -> None:
        self.app.save_selection(disease_hint=self.disease)  # type: ignore[attr-defined]
        self.app.exit()  # type: ignore[attr-defined]


class WorkbenchFinderApp(App):
    CSS = """
    #title {
        padding: 1 1;
    }
    #status {
        padding: 0 1;
        color: $text-muted;
    }
    """

    def __init__(self, cfg: FinderConfig) -> None:
        super().__init__()
        self.cfg = cfg
        self._indices_loaded = False
        self.disease_to_studies: Dict[str, set[str]] = {}
        self.study_summary_by_id: Dict[str, Dict[str, Any]] = {}
        # Per-disease facet filters used on the studies screen.
        self.disease_facet_filters: Dict[str, Dict[str, set[str]]] = {}
        # Advanced per-study facets derived from `/study_id/<ST>/analysis` (loaded on-demand).
        self.study_analysis_facets_by_id: Dict[str, Dict[str, set[str]]] = {}
        self._advanced_facets_loaded_diseases: set[str] = set()
        self.selected_analyses: set[str] = set()
        self.analysis_to_study: Dict[str, str] = {}
        self.saved_selection_path: Optional[Path] = None

    async def on_mount(self) -> None:
        self.push_screen(DiseaseScreen())

    async def ensure_indices_loaded(self, *, force_refresh: bool = False) -> None:
        if self._indices_loaded and not force_refresh:
            return
        cfg = self.cfg

        diseases_task = asyncio.to_thread(
            fetch_workbench_all_study_diseases,
            cache_root=cfg.cache_dir,
            timeout_s=float(cfg.timeout_s),
            retries=int(cfg.retries),
            delay_s=float(cfg.delay_s),
            refresh=bool(cfg.refresh or force_refresh),
        )
        summaries_task = asyncio.to_thread(
            fetch_workbench_all_study_summaries,
            cache_root=cfg.cache_dir,
            timeout_s=float(cfg.timeout_s),
            retries=int(cfg.retries),
            delay_s=float(cfg.delay_s),
            refresh=bool(cfg.refresh or force_refresh),
        )
        diseases_obj, summaries_obj = await asyncio.gather(diseases_task, summaries_task)
        if not isinstance(diseases_obj, dict) or not isinstance(summaries_obj, dict):
            raise RuntimeError("Failed to fetch Workbench indices (disease and/or summary).")

        # disease index rows are like: {"Study ID": "ST000010", "Disease": "Cancer"}
        mapping: Dict[str, set[str]] = {}
        for rec in diseases_obj.values():
            if not isinstance(rec, dict):
                continue
            sid = _canonical_study_id(rec.get("Study ID") or rec.get("study_id") or rec.get("study") or "")
            dis = _norm(rec.get("Disease") or rec.get("disease") or "")
            if not sid or not dis:
                continue
            mapping.setdefault(dis, set()).add(sid)
        self.disease_to_studies = mapping

        # summary index rows include: study_id, study_title, species, institute, etc.
        by_id: Dict[str, Dict[str, Any]] = {}
        for rec in summaries_obj.values():
            if not isinstance(rec, dict):
                continue
            sid = _canonical_study_id(rec.get("study_id") or rec.get("Study ID") or "")
            if not sid:
                continue
            by_id[sid] = dict(rec)
        self.study_summary_by_id = by_id

        self._indices_loaded = True

    def get_disease_study_rows(self, disease: str) -> List[Dict[str, str]]:
        """Build rows for the studies table from the cached disease + summary indices."""
        disease = str(disease)
        study_ids = sorted(self.disease_to_studies.get(disease, set()))
        rows: List[Dict[str, str]] = []
        for sid in study_ids:
            srec = self.study_summary_by_id.get(sid) or {}
            rows.append(
                {
                    "study_id": sid,
                    "study_title": _norm(srec.get("study_title") or srec.get("title") or ""),
                    "analysis_type": _norm(srec.get("analysis_type") or ""),
                    "n_samples": _norm(srec.get("number_of_samples") or ""),
                    "species": _norm(srec.get("species") or srec.get("organism") or ""),
                    "institute": _norm(srec.get("institute") or ""),
                }
            )
        rows.sort(key=lambda r: (r.get("study_id", "")))
        return rows

    def get_facet_filters(self, disease: str) -> Dict[str, set[str]]:
        disease = str(disease)
        filters = self.disease_facet_filters.setdefault(disease, {})
        for facet in _FACETS_ORDER:
            filters.setdefault(facet, set())
        return filters

    def clear_facet_filters(self, disease: str, *, facet: Optional[str] = None) -> None:
        filters = self.get_facet_filters(disease)
        if facet is None:
            for f in list(filters.keys()):
                filters[f].clear()
            return
        filters.setdefault(str(facet), set()).clear()

    def toggle_facet_filter(self, disease: str, facet: str, value: str) -> None:
        filters = self.get_facet_filters(disease)
        facet = str(facet)
        value = str(value)
        selected = filters.setdefault(facet, set())
        if value in selected:
            selected.remove(value)
        else:
            selected.add(value)

    def format_facet_filters_summary(self, disease: str, *, filters: Optional[Dict[str, set[str]]] = None) -> str:
        filters = filters if filters is not None else self.get_facet_filters(disease)
        parts: List[str] = []
        for facet in _FACETS_ORDER:
            selected = sorted(filters.get(facet) or [])
            if not selected:
                continue
            label = {
                FACET_ANALYSIS_TYPE: "analysis_type",
                FACET_SPECIES: "species",
                FACET_INSTITUTE: "institute",
                FACET_N_SAMPLES_BUCKET: "n_samples",
                FACET_ION_MODE: "ion_mode",
                FACET_CHROMATOGRAPHY: "chrom",
                FACET_MS_INSTRUMENT_TYPE: "instrument",
            }.get(facet, facet)
            if len(selected) <= 2:
                parts.append(f"{label}={','.join(selected)}")
            else:
                parts.append(f"{label}={len(selected)} selected")
        if not parts:
            return "Facets: none (press 'g' to add filters)"
        return "Facets: " + "; ".join(parts)

    def is_advanced_facets_loaded(self, disease: str) -> bool:
        return str(disease) in self._advanced_facets_loaded_diseases

    def study_row_passes_filters(
        self,
        disease: str,
        row: Dict[str, str],
        *,
        exclude_facet: Optional[str] = None,
    ) -> bool:
        """Return True if a study row passes all active facet filters (except exclude_facet)."""
        disease = str(disease)
        exclude_facet = str(exclude_facet) if exclude_facet else None
        filters = self.get_facet_filters(disease)
        study_id = _canonical_study_id(row.get("study_id"))

        def _selected(facet: str) -> List[str]:
            if exclude_facet and facet == exclude_facet:
                return []
            return sorted(filters.get(facet) or [])

        sel = _selected(FACET_ANALYSIS_TYPE)
        if sel:
            v = _norm(row.get("analysis_type")) or "<missing>"
            if v not in sel:
                return False

        sel = _selected(FACET_SPECIES)
        if sel:
            v = _norm(row.get("species")) or "<missing>"
            if v not in sel:
                return False

        sel = _selected(FACET_INSTITUTE)
        if sel:
            v = _norm(row.get("institute")) or "<missing>"
            if v not in sel:
                return False

        sel = _selected(FACET_N_SAMPLES_BUCKET)
        if sel:
            v = _samples_bucket(row.get("n_samples"))
            if v not in sel:
                return False

        if _selected(FACET_ION_MODE) or _selected(FACET_CHROMATOGRAPHY) or _selected(FACET_MS_INSTRUMENT_TYPE):
            facets_by_study = self.study_analysis_facets_by_id.get(study_id) or {}

            sel = _selected(FACET_ION_MODE)
            if sel:
                vals = set(facets_by_study.get(FACET_ION_MODE) or [])
                if not vals or not (vals & set(sel)):
                    return False

            sel = _selected(FACET_CHROMATOGRAPHY)
            if sel:
                vals = set(facets_by_study.get(FACET_CHROMATOGRAPHY) or [])
                if not vals or not (vals & set(sel)):
                    return False

            sel = _selected(FACET_MS_INSTRUMENT_TYPE)
            if sel:
                vals = set(facets_by_study.get(FACET_MS_INSTRUMENT_TYPE) or [])
                if not vals or not (vals & set(sel)):
                    return False

        return True

    async def ensure_disease_advanced_facets_loaded(
        self,
        disease: str,
        *,
        study_ids: Optional[Sequence[str]] = None,
        progress_cb=None,
        max_concurrency: int = 3,
    ) -> None:
        """Load per-study analysis listings for a disease to enable advanced facets."""
        disease = str(disease)
        await self.ensure_indices_loaded()

        cfg = self.cfg

        if study_ids is None:
            target_ids = sorted(self.disease_to_studies.get(disease, set()))
        else:
            target_ids = sorted({_canonical_study_id(x) for x in study_ids if _canonical_study_id(x)})

        to_load = [sid for sid in target_ids if sid not in self.study_analysis_facets_by_id]
        total = int(len(to_load))
        if total <= 0:
            # Mark the disease "fully loaded" if we now have all its studies' facets.
            all_ids = self.disease_to_studies.get(disease, set())
            if all_ids and all(sid in self.study_analysis_facets_by_id for sid in all_ids):
                self._advanced_facets_loaded_diseases.add(disease)
            return

        sem = asyncio.Semaphore(max(1, int(max_concurrency)))
        done = 0

        async def _one(sid: str) -> None:
            nonlocal done
            async with sem:
                obj = await asyncio.to_thread(
                    fetch_workbench_study_analyses,
                    sid,
                    cache_root=cfg.cache_dir,
                    timeout_s=float(cfg.timeout_s),
                    retries=int(cfg.retries),
                    delay_s=float(cfg.delay_s),
                    refresh=bool(cfg.refresh),
                )
                ion: set[str] = set()
                chrom: set[str] = set()
                inst: set[str] = set()
                if isinstance(obj, dict):
                    for rec in obj.values():
                        if not isinstance(rec, dict):
                            continue
                        v = _canonical_ion_mode(rec.get("ion_mode"))
                        if v:
                            ion.add(v)
                        v = _canonical_chromatography(rec.get("chromatography_type") or rec.get("chromatography system") or "")
                        if v:
                            chrom.add(v)
                        v = _canonical_instrument_type(rec.get("ms_instrument_type"))
                        if v:
                            inst.add(v)
                self.study_analysis_facets_by_id[sid] = {
                    FACET_ION_MODE: ion,
                    FACET_CHROMATOGRAPHY: chrom,
                    FACET_MS_INSTRUMENT_TYPE: inst,
                }
                done += 1
                if progress_cb:
                    try:
                        progress_cb(f"Loaded advanced facets: {done}/{total} studies")
                    except Exception:
                        pass

        await asyncio.gather(*[_one(sid) for sid in to_load])

        all_ids = self.disease_to_studies.get(disease, set())
        if all_ids and all(sid in self.study_analysis_facets_by_id for sid in all_ids):
            self._advanced_facets_loaded_diseases.add(disease)

    def save_selection(self, *, disease_hint: str = "") -> Path:
        cfg = self.cfg
        analysis_ids = sorted({_canonical_analysis_id(x) for x in self.selected_analyses if _canonical_analysis_id(x)})
        study_ids = sorted({_canonical_study_id(self.analysis_to_study.get(a, "")) for a in analysis_ids if self.analysis_to_study.get(a)})
        payload = {
            "schema_version": "mass_sight.workbench_selection.v1",
            "created_at_utc": _utc_now_iso(),
            "disease": str(disease_hint or ""),
            "study_ids": [str(x) for x in study_ids if str(x)],
            "analysis_ids": [str(x) for x in analysis_ids if str(x)],
            "notes": "Generated by mass_sight find (Textual TUI).",
        }
        out_path = Path(cfg.out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        self.saved_selection_path = out_path
        return out_path


def run_workbench_finder_tui(
    *,
    out_path: Path,
    cache_dir: Path,
    refresh: bool = False,
    timeout_s: float = 60.0,
    retries: int = 3,
    delay_s: float = 0.25,
) -> Optional[Path]:
    cfg = FinderConfig(
        cache_dir=Path(cache_dir),
        refresh=bool(refresh),
        timeout_s=float(timeout_s),
        retries=int(retries),
        delay_s=float(delay_s),
        out_path=Path(out_path),
    )
    app = WorkbenchFinderApp(cfg)
    app.run()
    return app.saved_selection_path
