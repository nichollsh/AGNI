# Visual language

This page is the precise reference for AGNI documentation site's palette, fonts, logo usage, and badge conventions. For how the site is built, see [Documentation build pipeline](@ref "Documentation build pipeline").

All prose should adopt British English spellings and grammar, rather than American English.

## Palette

Currently defined in `docs/src/assets/style.scss`:

| Role | Colour |
|:-----|:------|
| Main colour (`$maincolor`) — sidebar background, bold/heading text | `#1c2b4b` (dark navy) |
| Secondary colour (`$secondcolor`) — links and accents | `#1B6FA8` |
| Page background (`$mainwhite`) | `#F2F5F7` (off-white) |
| Body text (`$mainblack`) | `#10151B` |

AGNI currently uses one visual theme for both light and dark mode.

### Highlight colours

The brand palette also defines a set of highlight colours beyond what is currently in `style.scss`:

| Role | Colour | Incorporated into `style.scss`? |
|:-----|:------|:--------------------------|
| White (background) | `#F2F5F7` | Yes — `$mainwhite` |
| Text | `#10151B` | Yes — `$mainblack` |
| Links / actions | `#1B6FA8` | Yes — `$secondcolor` |
| Highlight | `#C2362B` | No |
| Highlight | `#C8860F` | No |
| Highlight | `#57A05C` | No |


## Fonts

Three families, bundled locally as variable-weight TTF files under `docs/src/assets/fonts/` (SIL OFL-1.1 licensed; license text ships alongside each font as `OFL-<family>.txt`) and declared via `@font-face` rules at the top of
`docs/src/assets/style.scss`. Fonts are self-hosted rather than loaded from Google Fonts.

| Family | Axis range | Used for |
|:-------|:-----------|:---------|
| **Sora** | `wght` 100–800 | Titles, headings (`h1`–`h6`), and the wordmark (`.docs-package-name`). |
| **Instrument Sans** (`$family-sans-serif`) | `wght` 400–700, `wdth` 75–100 | Prose, documentation body text, and interface chrome (sidebar, search box, buttons). |
| **Spline Sans Mono** (`$family-monospace`) | `wght` 300–700 | Code (inline and blocks), the version selector, and small metadata text (captions, timestamps). |

When updating a font file, re-check its `fvar` axis table (e.g. with `fontTools.ttLib`) before assuming the `font-weight`/`font-stretch` ranges in `style.scss` still match.

## Logo

`docs/src/assets/logo_title_{light,dark}.svg` — the light/dark variants are selected in Markdown via the `display-light-only` / `display-dark-only` CSS classes (see the top of `docs/src/index.md`). The favicon uses `assets/logo.ico`.

## Badges

Status badges (coverage, CI, test counts) are embedded as raw HTML `<img>` tags with an explicit `height:20px;width:auto` inline style, so that badges from different sources (Shields.io, Codecov, a self-hosted gist SVG) render at a consistent height without distorting their aspect ratio. See [Testing suite](@ref "Testing suite") for
the current example.

## Keeping this page in sync

Anyone changing the SCSS files under `docs/src/assets/`, the logo/font assets, or the
badge-embedding convention should update this page in the same change.

