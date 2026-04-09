# Wiki Usage

- Before diving into code you haven't read, check `wiki/index.md` for relevant articles
- The wiki has two sections:
  - **Codebase wiki** (`wiki/pipeline/`, `wiki/concepts/`, `wiki/reference/`, `wiki/guides/`) — covers pipeline stages, ABCD method, isolation, shower shapes, unfolding, systematics, config schema, constants sync, MC samples, file inventory, and step-by-step guides
  - **Physics wiki** (`wiki/physics/`) — domain knowledge from ~70 research papers, organized by concept (measurements, theory, techniques, detector). Use `wiki/physics/index.md` as entry point
- For **physics domain questions** (how other experiments handle isolation, NLO theory, detector performance, what measurements exist at different sqrt(s)), read the relevant `wiki/physics/` article FIRST
- For **paper-specific details**, read the raw summary at `wiki/raw/papers/{category}/{arxiv_id}.md`
- For **methodology comparisons** (e.g., "how does ATLAS do ABCD vs PPG12?"), the `wiki/physics/techniques/` articles synthesize across experiments
- If you discover something not covered in the wiki, flag it to the user
- Do NOT blindly trust wiki claims — verify against current code when acting on them (the wiki was built at a point in time and code may have changed since)
- When modifying constants (cross-section weights, pT bins, luminosity), check `wiki/reference/constants-sync.md` for all locations that must stay in sync
- When adding a systematic variation, follow `wiki/guides/adding-a-systematic.md`
- The wiki's `_reports/` directory contains raw exploration reports with more detail than the articles
- The physics wiki's `wiki/raw/manifest.md` tracks all ingested papers and their availability tiers
