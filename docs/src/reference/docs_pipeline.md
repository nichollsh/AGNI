# Documentation build pipeline

This page briefly describes how the AGNI documentation site is built. For the palette,
fonts, logo usage, and badge conventions, see [Visual language](@ref "Visual language").

## Build pipeline

The site is built by [`docs/make.jl`](https://github.com/nichollsh/AGNI/blob/main/docs/make.jl)
using [Documenter.jl](https://documenter.juliadocs.org/), with two plugins:

* [`DocumenterCitations`](https://github.com/JuliaDocs/DocumenterCitations.jl) — expands
  `@bibliography` blocks from `docs/src/assets/refs.bib` (see
  [Bibliography](@ref "Bibliography")).
* [`DocumenterPages`](https://github.com/asinghvi17/DocumenterPages.jl) — provides the
  `PageNode` construct used to build the nested sidebar (How-to guides / Tutorials /
  Explanation / Reference), following the [Diátaxis](https://diataxis.fr/) framework.

Before `makedocs` runs, `make.jl` compiles the custom theme from SCSS to CSS using
`DocumenterTools.Themes.compile`. The file `docs/src/assets/style.scss` (palette, fonts, sidebar styling) is concatenated with `docs/src/assets/lightdefs.scss` (Documenter theme variable overrides) into `light.scss`. Then, `light.scss` is copied to `dark.scss`. AGNI currently uses one visual theme for both light and dark mode, rather than a separate dark palette.


`makedocs` then renders the page tree defined in `make.jl`, pulling in local assets
(`assets/style.css`, `assets/logo.ico`). All three site fonts are bundled as local
variable TTF files under `docs/src/assets/fonts/` and loaded via `@font-face` rules in
`style.scss` — there are no remote font requests. Docstrings are pulled automatically
from the source code into `reference/api.md`.

Deployment (`deploydocs`) pushes the built site to GitHub Pages via the
`documentation.yml` GitHub Actions workflow, with `push_preview=true` so that pull
requests get a preview build.

## Keeping this page in sync

Anyone changing `docs/make.jl` (page structure, plugins, `deploydocs` config) should update this page in the same change. Changes to the palette, fonts, logo, or badge conventions belong in [Visual language](@ref "Visual language") instead.

