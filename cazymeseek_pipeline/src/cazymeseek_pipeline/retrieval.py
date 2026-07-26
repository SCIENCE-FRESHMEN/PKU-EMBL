"""CAZymeSeek's four-store exact-plus-semantic knowledge retrieval.

The EB05 project document lists CAZypedia, CAZy, CGC and substrate resources.
The supplied vector database README specifies ChromaDB, local all-MiniLM-L6-v2,
cosine distance, family_id and cazyme_families_base join keys.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import chromadb
from sentence_transformers import SentenceTransformer


@dataclass(frozen=True)
class Store:
    directory: str
    collection: str
    family_field: str


STORES = {
    "cazypedia": Store("01_cazypedia", "cazypedia", "family_id"),
    "cazy": Store("02_cazy_database", "cazy_db", "family_id"),
    "cgc": Store("03_cgc_database", "cgc_db", "cazyme_families_base"),
    "reaction": Store("04_substrate_specificity", "substrate_specificity", "family_id"),
}


class CAZymeSeekRetriever:
    """Retrieve evidence records; MiniLM never predicts an enzyme or substrate."""

    def __init__(self, root: str | Path):
        self.root = Path(root)
        self.model = SentenceTransformer(str(self.root / "models" / "all-MiniLM-L6-v2"))
        self.collections = {}
        self.supplementary: dict[str, dict[str, Any]] = {}
        for name, spec in STORES.items():
            client = chromadb.PersistentClient(path=str(self.root / spec.directory))
            self.collections[name] = client.get_collection(spec.collection)
            with (self.root / spec.directory / "supplementary.json").open(encoding="utf-8") as handle:
                self.supplementary[name] = json.load(handle)

    def search(self, family_id: str | None, query: str, top_k: int = 3) -> list[dict[str, Any]]:
        """Fuse exact project-key links (1.0) with stored cosine semantic similarity.

        Similarity is a retrieval-rank value only. The documents specify no biological
        validation threshold, so this class deliberately does not filter by one.
        """
        hits: dict[tuple[str, str], dict[str, Any]] = {}
        embedding = self.model.encode([query], normalize_embeddings=True).tolist()
        for name, spec in STORES.items():
            collection = self.collections[name]
            semantic = collection.query(query_embeddings=embedding, n_results=top_k,
                                        include=["documents", "metadatas", "distances"])
            self._add_semantic(hits, name, semantic)
            if family_id:
                where = {spec.family_field: {"$contains": family_id}} if name == "cgc" else {"family_id": {"$eq": family_id}}
                exact = collection.get(where=where, include=["documents", "metadatas"])
                self._add_exact(hits, name, exact)
        return sorted(hits.values(), key=lambda row: row["confidence"], reverse=True)

    def _payload(self, name: str, item_id: str, document: str, metadata: dict[str, Any]) -> dict[str, Any]:
        full = self.supplementary[name].get(item_id, {})
        return {"store": name, "id": item_id, "text": full.get("text", document),
                "metadata": full if full else metadata}

    def _add_semantic(self, hits, name, result):
        for item_id, doc, meta, distance in zip(result["ids"][0], result["documents"][0], result["metadatas"][0], result["distances"][0]):
            item = self._payload(name, item_id, doc, meta)
            # The project indexes use cosine distance; ChromaDB returns distance = 1-cosine.
            item.update(match="semantic", confidence=round(max(0.0, 1.0 - float(distance)), 4))
            hits[(name, item_id)] = item

    def _add_exact(self, hits, name, result):
        for item_id, doc, meta in zip(result["ids"], result["documents"], result["metadatas"]):
            item = self._payload(name, item_id, doc, meta)
            item.update(match="exact_family", confidence=1.0)
            hits[(name, item_id)] = item
