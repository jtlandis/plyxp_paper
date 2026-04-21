// biorxiv.typ — Typst template inspired by bioRxiv.cls (Royle lab / PNAS-style)

// -----------------------------------------------------------------------------
// USER-FACING INTERFACE
// -----------------------------------------------------------------------------

// Document modes: "submit" (one column, line numbers, no logo) vs "preprint"
#let biorxiv-mode = "preprint" // or "submit"

// Page geometry: A4, total={170mm,247mm}, left=20mm, top=25mm
// -> Right margin = 210 - 20 - 170 = 20mm; bottom margin = 297 - 25 - 247 = 25mm.
#set page(
  paper: "a4",
  margin: (
    top: 25mm,
    bottom: 25mm,
    left: 20mm,
    right: 20mm,
  ),
  columns: if biorxiv-mode == "submit" { 1 } else { 2 },
  column-gap: 7mm,
)

// Basic typography: you can swap fonts if you want closer LaTeX matching.
#define main-sans "Linux Biolinum"  // roughly like lmodern sans
#define main-serif "Linux Libertine" // roughly like lmodern/Times combo
#set text(font: main-serif, size: 10pt, leading: 12pt)
#set par(justify: true, first-line-indent: 0pt)

// Colors
#define black50 = rgb(0.5, 0.5, 0.5)
#define color0 = rgb(0, 0, 0)
#define color1 = rgb(0.16, 0.16, 0.16)

// --- Line numbers (submit mode only) -----------------------------------------
#show: if biorxiv-mode == "submit" {
  show par: it => {
    numbering(
      it,
      position: "left",
      text: n => text(size: 7pt, font: main-sans, fill: black50, str(n)),
      step: 1,
      // adjust gap as desired
      gap: 1em,
    )
  }
} else {
  // No automatic line numbering
  _
}

// -----------------------------------------------------------------------------
// META FIELDS AND GLOBAL STATE
// -----------------------------------------------------------------------------

// We'll hold meta state in a dictionary captured by closures
#let _meta = (
  title: none,
  authors: (),
  affiliations: (),
  short-title: none,
  lead-author: none,
)

// Short title (for footer)
#let short-title(title: none) = {
  // returns nothing, just updates _meta via shadowing
  #let _meta = _meta + (short-title: title)
}

// Lead author (for footer)
#let lead-author(name: none) = {
  #let _meta = _meta + (lead-author: name)
}

// Author & affiliation helpers
#let author(name, affil: none, orcid: none, corresponding: false) = {
  (
    name: name,
    affil: affil,
    orcid: orcid,
    corresponding: corresponding,
  )
}

#let affiliation(id, text) = (id: id, text: text)

// We'll define a container that stores the doc meta and prints the title page.
#let biorxiv(
  title,
  authors: (),
  affiliations: (),
  short-title: none,
  lead-author: none,
  content,
) = {
  // Flatten authors/affiliations into a dictionary we can reuse
  #let _lead =
    if lead-author != none { lead-author }
    else if len(authors) > 0 { authors.at(0).name }
    else { none }

  // Cache meta; in Typst we just keep them as locals and pass to header/footer
  #let meta = (
    title: title,
    authors: authors,
    affiliations: affiliations,
    short-title:
      if short-title != none { short-title }
      else { title },
    lead-author: _lead,
  )

  // Set up page headers/footers via `show page`
  #show: doc => {
    show page.with(
      header: page-header(meta),
      footer: page-footer(meta),
    ): doc
  }

  // Render the front-matter title block, then main content
  #biorxiv-title-block(meta)
  #content
}

// -----------------------------------------------------------------------------
// TITLE BLOCK (similar to \maketitle / @maketitle)
// -----------------------------------------------------------------------------

#let biorxiv-title-block(meta) = {
  // Title font: large sans, bold; mode-dependent size
  #let title-size =
    if biorxiv-mode == "submit" { 18pt } else { 18pt }
  #let title-leading =
    if biorxiv-mode == "submit" { 20pt } else { 20pt }

  // Author font / affil font (mode-dependent)
  #let author-size = 9pt
  #let author-leading = if biorxiv-mode == "submit" { 10pt } else { 11pt }
  #let affil-size = if biorxiv-mode == "submit" { 9pt } else { 7pt }
  #let affil-leading = if biorxiv-mode == "submit" { 10pt } else { 8pt }

  // Prepare author line: "Name^1,2, Name^2, and Name^3"
  #let author-text = {
    // group authors with indices to allow "et al." logic elsewhere if needed
    let parts = ()
    for a in meta.authors {
      let sup = if a.affil == none {
        ""
      } else if type(a.affil) == "array" {
        // {1,2} -> "^{1,2}"
        "^(" + str(a.affil.join(",")) + ")"
      } else {
        "^(" + str(a.affil) + ")"
      }

      let label = a.name + sup
      // append ORCID as small icon-like in brackets if present
      if a.orcid != none {
        label = label + " [" + a.orcid + "]"
      }
      parts.push(label)
    }
    if len(parts) == 0 {
      ""
    } else if len(parts) == 1 {
      parts.at(0)
    } else if len(parts) == 2 {
      parts.at(0) + " and " + parts.at(1)
    } else {
      // Oxford-comma style listing; roughly matches `Authands{, and }`
      parts.slice(0, -1).join(", ") + ", and " + parts.last()
    }
  }

  // Map affiliation id -> text for printing
  #let affil-map = meta.affiliations

  // Typeset title + authors + affiliations flush-left
  align(left)[
    // Title
    #block(
      spacing: 1em,
      text(
        font: main-sans,
        weight: "bold",
        size: title-size,
        leading: title-leading,
        meta.title,
      ),
    )

    // Author line
    #block(
      spacing: 0.4em,
      text(
        font: main-sans,
        weight: "bold",
        size: author-size,
        leading: author-leading,
        fill: color1,
        author-text,
      ),
    )

    // Affiliations in single paragraph separated by "; "
    #if len(affil-map) > 0 [
      #block(
        spacing: 1.4em,
        text(
          font: main-sans,
          size: affil-size,
          leading: affil-leading,
          fill: color1,
          {
            for (idx, aff) in affil-map.enumerate() {
              let sup = (idx + 1)
              // expect aff.type == (id, text) or just text
              let content =
                if type(aff) == "dictionary" and "text" in aff {
                  aff.text
                } else {
                  aff
                }
              if idx > 0 { "; " }
              raw("$^{" + str(sup) + "}$ " + content)
            }
          },
        ),
      )
    ]
  ]

  // Vertical space below title block
  v(1.4em)
}

// -----------------------------------------------------------------------------
// ABSTRACT / KEYWORDS / CORRESPONDENCE / ACKS / CONTRIBUTIONS / INTERESTS
// -----------------------------------------------------------------------------

#let abstract(body) = {
  // submit: normalsize bold; else: small or bold
  let size = if biorxiv-mode == "submit" { 10pt } else { 8pt }
  let leading = if biorxiv-mode == "submit" { 12pt } else { 11pt }

  // boldface or not: we mimic \bfabstracttrue default
  let is-bold = true

  block(spacing: 1em, inset: (top: 1em))[
    #text(
      font: main-sans,
      size: size,
      leading: leading,
      weight: if is-bold { "bold" } else { "regular" },
      body,
    )
  ]
}

#let keywords(body) = {
  let size = if biorxiv-mode == "submit" { 10pt } else { 7pt }
  let leading = if biorxiv-mode == "submit" { 11pt } else { 10pt }
  block(spacing: 0.5em)[
    #text(
      font: if biorxiv-mode == "submit" { main-serif } else { main-serif },
      size: size,
      leading: leading,
      weight: if biorxiv-mode == "submit" { "regular" } else { "bold" },
      "Keywords: " + body,
    )
  ]
}

#let correspondence(body) = {
  let size = if biorxiv-mode == "submit" { 10pt } else { 9pt }
  let leading = if biorxiv-mode == "submit" { 11pt } else { 10pt }
  block(spacing: 0.5em)[
    #text(
      font: main-sans,
      size: size,
      leading: leading,
      weight: if biorxiv-mode == "submit" { "regular" } else { "bold" },
      "Correspondence: ",
    )
    #text(
      font: main-sans,
      size: size,
      leading: leading,
      style: "italic",
      body,
    )
  ]
}

#let heading-block(label, body) = {
  let size = if biorxiv-mode == "submit" { 10pt } else { 9pt }
  let leading = if biorxiv-mode == "submit" { 11pt } else { 10pt }
  block(spacing: 0.5em, inset: (top: 1em, bottom: 0.3em))[
    #text(
      font: main-sans,
      size: size,
      leading: leading,
      weight: "bold",
      label,
    )
    #text(
      font: main-sans,
      size: size,
      leading: leading,
      weight: "regular",
      body,
    )
  ]
}

#let acknowledgements(body) = heading-block("ACKNOWLEDGEMENTS", body)
#let contributions(body) = heading-block("AUTHOR CONTRIBUTIONS", body)
#let interests(body) = heading-block("COMPETING FINANCIAL INTERESTS", body)

// -----------------------------------------------------------------------------
// SECTIONING: "Supplementary Note N: Title", etc.
// -----------------------------------------------------------------------------

// Section numbering like: "Supplementary Note 1: Title"
#show heading.where(level: 1): it => {
  let num = it.numbering("1")
  // if unnumbered, just show the title
  if num == none {
    heading(
      level: 1,
      outlined: it.outlined,
      supplement: none,
      body: [
        #text(font: main-sans, size: 14pt, weight: "bold", it.body)
      ],
    )
  } else {
    heading(
      level: 1,
      outlined: it.outlined,
      supplement: none,
      body: [
        #text(
          font: main-sans,
          size: 14pt,
          weight: "bold",
          "Supplementary Note " + num + ": " + it.body,
        )
      ],
    )
  }
}

// Subsection: "A. Title"
#show heading.where(level: 2): it => {
  let num = it.numbering("A")
  heading(
    level: 2,
    outlined: it.outlined,
    supplement: none,
    body: [
      #text(font: main-sans, size: 10pt, weight: "bold", num + ". " + it.body)
    ],
  )
}

// Subsubsection: small sans, italic, run-in style
#show heading.where(level: 3): it => {
  let num = it.numbering("1")
  // Run-in effect
  block(spacing: 0.2em)[
    #text(
      font: main-sans,
      size: 9pt,
      weight: "bold",
      style: "italic",
      num + ". " + it.body + ". ",
    )
  ]
}

// Paragraph heading: small bold sans, run-in
#let para-heading(body) = {
  text(font: main-sans, size: 9pt, weight: "bold", body + " ")
}

// -----------------------------------------------------------------------------
// FOOTERS / HEADERS (fancyhdr-like)
// -----------------------------------------------------------------------------

// Equivalent to \aBioRXiv, \bBioRXiv, \cBioRXiv combined with dAuthor.
// We'll keep it simple: left footer: lead author (+ et al.) | short title
// right footer: page X–LastPage
#let page-footer(meta) = {
  // lead author + "et al." if > 2
  let n-auth = len(meta.authors)
  let lead = meta.lead-author
  let left-auth = {
    if lead == none { "" }
    else if n-auth > 2 { lead + " et al." }
    else { lead }
  }


  let short = meta.short-title
  let left-part = {
    if short == none or short == "" {
      left-auth
    } else if left-auth == "" {
      short
    } else {
      left-auth + " | " + short
    }
  }


  // Page X–Y style
  let total-pages = counter(page).final()
  let right-part = str(counter(page).get()) + "–" + str(total-pages)

  // submit vs preprint: do not show logo in submit mode
  let logo = {
    if biorxiv-mode == "preprint" {
      // If you have a logo image, enable:
      // image("bioRxiv_logo.png", width: 1cm)
      none
    } else { none }
  }


  // Compose footer with two columns
  box(
    width: 100%,
    inset: (top: 0.5em),
  )[
    #align(left)[
      #text(font: main-sans, size: 7pt, fill: color1, left-part)
    ]
    #align(right)[
      #if logo != none [
        #logo
        #text(" | ")
      ]
      #text(font: main-sans, size: 7pt, fill: color1, right-part)
    ]
  ]
}

#let page-header(meta) = {
  // Keep headers empty (like \headrulewidth=0), but you could add section title
  none
}

// -----------------------------------------------------------------------------
// FIGURE CAPTIONS
// -----------------------------------------------------------------------------

// Approximate small caption format: small sans, bold "Figure" name.
#show figure: it => {
  show it: set text(font: main-sans, size: if biorxiv-mode == "submit" { 9pt } else { 7pt }, leading: if biorxiv-mode == "submit" { 12pt } else { 9pt })
  it
}

// -----------------------------------------------------------------------------
// SIMPLE BIBLIOGRAPHY STYLE (natbib-like feel)
// -----------------------------------------------------------------------------

// You can plug any citation style you like. Here we just tweak text size.
#show bibliography: it => {
  // Small sans-serif bibliography
  show it: set text(font: main-sans, size: if biorxiv-mode == "submit" { 9pt } else { 6pt }, leading: if biorxiv-mode == "submit" { 10pt } else { 8pt })
  it
}
