// BioRxiv-style Typst template adapted from Quarto's default article template
// and the Royle-lab bioRxiv.cls.

// Main entrypoint used by typst-show.typ
#let article(
  title: none,
  subtitle: none,
  authors: none,
  date: none,
  abstract: none,
  abstract-title: none,
  cols: 1,
  margin: (x: 1.25in, y: 1.25in),
  paper: "us-letter",
  lang: "en",
  region: "US",
  font: "libertinus serif",
  fontsize: 11pt,
  title-size: 1.5em,
  subtitle-size: 1.25em,
  heading-family: "libertinus serif",
  heading-weight: "bold",
  heading-style: "normal",
  heading-color: black,
  heading-line-height: 0.65em,
  sectionnumbering: none,
  pagenumbering: "1",
  toc: false,
  toc_title: none,
  toc_depth: none,
  toc_indent: 1.5em,
  doc,
) = {
  // ---------------------------------------------------------------------------
  // BioRxiv-style page geometry (overrides some arguments for consistency)
  // - A4
  // - geometry: total={170mm,247mm}, left=20mm, top=25mm
  //   => margins: left/right 20mm, top/bottom 25mm
  // - columns controlled by `cols` (1 or 2)
  // ---------------------------------------------------------------------------
  let paper-size = "a4"
  let biorxiv-margin = (
    top: 25mm,
    bottom: 25mm,
    left: 20mm,
    right: 20mm,
  )

  set page(
    paper: paper-size,
    margin: biorxiv-margin,
    numbering: pagenumbering,
    columns: cols,
    //column-gap: 7mm,
  )

  // ---------------------------------------------------------------------------
  // Base text and paragraph settings
  // ---------------------------------------------------------------------------
  set par(justify: true, first-line-indent: 0pt)
  // Approximate 10pt text with modest leading
  set par(leading: 12pt)
  set text(
    lang: lang,
    region: region,
    font: "libertinus serif",
    size: 10pt,
  )

  // Section numbering behavior is still controlled by Quarto
  set heading(numbering: sectionnumbering)

  // ---------------------------------------------------------------------------
  // TITLE BLOCK (BioRxiv-style)
  //   - Flush-left title in large bold sans
  //   - Authors line below, then affiliations
  //   - Reuses Quarto's `authors` array (name, affiliation, email)
  // ---------------------------------------------------------------------------

  if title != none {
    // Use sans for the title
    let title-font = font
    let title-size = 18pt
    let title-leading = 20pt

    // Author font / affil font (single style here; could branch for submit mode)
    let author-font = font
    let author-size = 9pt
    let author-leading = 11pt
    let affil-font = font
    let affil-size = 7pt
    let affil-leading = 8pt

    // Build author line from Quarto's author metadata
    let author-line = if authors == none or authors.len() == 0 {
      ""
    } else {
      let parts = ()
      for a in authors {
        // If Quarto has `affiliation` as text, we just show the name (no superscripts).
        // If you move to numeric affiliations, we can change this to add ^1,2, etc.
        let label = a.name
        parts.push(label)
      }
      if parts.len() == 1 {
        parts.at(0)
      } else if parts.len() == 2 {
        parts.at(0) + " and " + parts.at(1)
      } else {
        parts.slice(0, -1).join(", ") + ", and " + parts.last()
      }
    }

    // Collect unique affiliations from authors
    let affils = if authors == none { () } else {
      let seen = ()
      for a in authors {
        if "affiliation" in a and a.affiliation != none {
          if not seen.contains(a.affiliation) {
            seen.push(a.affiliation)
          }
        }
      }
      seen
    }

    // Title + authors + affiliations aligned left (BioRxiv style)
    align(left)[
      // Title
      #block(spacing: 1em)[
        #set par(leading: title-leading) // Extra leading for title block
        #text(
          font: title-font,
          size: title-size,
          weight: "bold",
        )[ #title ]
      ]

      // Authors
      #if author-line != "" [
        #block(spacing: 0.4em)[
          #set par(leading: author-leading) // Extra leading for author block
          #text(
            font: author-font,
            size: author-size,
            weight: "bold",
          )[ #author-line ]
        ]
      ]

      // Affiliations (single paragraph, separated by "; ")
      #if affils.len() > 0 [
        #set par(leading: affil-leading) // Extra leading for affil block
        #block(spacing: 1.4em)[
          #text(
            font: affil-font,
            size: affil-size,
          )[
            // Join with "; "
            #for (idx, aff) in affils.enumerate() {
              if idx > 0 { "; " }
              aff
            }
          ]
        ]
      ]
    ]
  }

  // ---------------------------------------------------------------------------
  // DATE (if provided by Quarto)
  // ---------------------------------------------------------------------------
  if date != none {
    // BioRxiv doesn't usually print a big centered date in the title block,
    // but we allow a small date below if requested.
    align(left)[
      #block(inset: (top: 0.5em, bottom: 1em))[
        #set par(leading: 8pt) // Extra leading for date block
        #text(
          font: font,
          size: 7pt,
        )[ #date ]
      ]
    ]
  }

  // ---------------------------------------------------------------------------
  // ABSTRACT (BioRxiv-like: small sans, bold "Abstract" heading)
  // Quarto passes both abstract-title and abstract as strings/blocks.
  // ---------------------------------------------------------------------------
  if abstract != none {
    // BioRxiv uses sans bold for the abstract in preprint mode.
    let abs-font = font
    let abs-size = 8pt
    let abs-leading = 11pt

    block(inset: (top: 1em, bottom: 1em))[
      #set par(leading: abs-leading) // Extra leading for abstract block
      #text(
        font: abs-font,
        size: abs-size,
        weight: "bold",
      )[
        // Use provided abstract-title or default to "Abstract"
        #if abstract-title != none {
          abstract-title
        }
      ]
      #h(0.75em)
      #text(
        font: abs-font,
        size: abs-size,
      )[ #abstract ]
    ]
  }

  // ---------------------------------------------------------------------------
  // TABLE OF CONTENTS (if requested)
  // ---------------------------------------------------------------------------
  if toc {
    let t = if toc_title == none { auto } else { toc_title }
    block(above: 0em, below: 2em)[
      #outline(
        title: t,
        depth: toc_depth,
        indent: toc_indent,
      );
    ]
  }

  // ---------------------------------------------------------------------------
  // HEADINGS: approximate bioRxiv typography (sans, bold, small).
  // This is a *light* mapping; if you want Supplementary Note style,
  // we can specialize heading(level: 1/2/3) further.
  // ---------------------------------------------------------------------------
  show heading.where(level: 1): it => {
    set heading(supplement: none)
    set text(font: font, size: 12pt, weight: "bold")
    it
  }

  show heading.where(level: 2): it => {
    set heading(supplement: none)
    set text(font: font, size: 10pt, weight: "bold")
    it
  }

  show heading.where(level: 3): it => {
    set heading(supplement: none)
    set text(font: font, size: 9pt, weight: "bold", style: "italic")
    it
  }

  // ---------------------------------------------------------------------------
  // TABLES: keep Quarto's simple styling
  // ---------------------------------------------------------------------------
  set table(
    inset: 6pt,
    stroke: none,
  )

  // ---------------------------------------------------------------------------
  // MAIN DOCUMENT CONTENT
  // columns(): we already set page.columns via `cols`, so just emit doc.
  // ---------------------------------------------------------------------------
  doc
}
