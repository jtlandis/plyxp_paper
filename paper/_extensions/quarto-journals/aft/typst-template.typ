
// This is an example typst template (based on the default template that ships
// with Quarto). It defines a typst function named 'article' which provides
// various customization options. This function is called from the
// 'typst-show.typ' file (which maps Pandoc metadata function arguments)
//
// If you are creating or packaging a custom typst template you will likely
// want to replace this file and 'typst-show.typ' entirely. You can find
// documentation on creating typst templates and some examples here:
//   - https://typst.app/docs/tutorial/making-a-template/
//   - https://github.com/typst/templates


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
  set page(
    paper: paper,
    margin: margin,
    numbering: pagenumbering,
  )
  set par(justify: true)
  set text(lang: lang, region: region, font: font, size: fontsize)
  set heading(numbering: sectionnumbering)
  if title != none {
    align(left)[#block[
      #set par(leading: heading-line-height)
      #if (
        heading-family != none
          or heading-weight != "bold"
          or heading-style != "normal"
          or heading-color != black
          or heading-decoration == "underline"
          or heading-background-color != none
      ) {
        set text(font: heading-family, weight: heading-weight, style: heading-style, fill: heading-color)
        text(size: title-size)[#title]
        if subtitle != none {
          parbreak()
          text(size: subtitle-size)[#subtitle]
        }
      } else {
        text(weight: "bold", size: title-size)[#title]
        if subtitle != none {
          parbreak()
          text(weight: "bold", size: subtitle-size)[#subtitle]
        }
      }
    ]]
  }
  // build out authors
  if authors != none {
    let count = authors.len()
    let affils = (:)
    let n = 0
    for author in authors {
      text(weight: "semibold")[#author.name]
      if author.affiliation != none {
        // typst-show uses ||| delemiter between affiliations
        let auth_aff = author.affiliation.at("text").split("|||")
        let aff_num = none
        let first = true
        for aff in auth_aff {
          let aff = aff.trim()
          let num = affils.at(aff, default: 0)
          // affiliation not found
          if num == 0 {
            n += 1
            num = n
            // insert it
            affils.insert(aff, n)
          }
          if not first {
            aff_num = aff_num + [,]
          }
          first = false
          aff_num = aff_num + [#num]
        }
        if aff_num != none {
          super[#aff_num]
        }
        [   ]
      }
    }
    [\ ]
    n = 0
    let first = true
    for key in affils.keys() {
      set text(size: fontsize - 2pt)
      n += 1
      if not first {
        [;   ]
      }
      super[#n] + [#key]
      first = false
    }
  }

  if date != none {
    align(center)[#block(inset: 1em)[
      #date
    ]]
  }

  if abstract != none {
    if abstract-title == "Abstract" {
      // likely not set by default
      abstract-title = none
    } else {
      abstract-title = [#text(weight: "semibold")[#abstract-title] #h(1em)]
    }
    abstract = block[
      #abstract-title #text(weight: "semibold")[#abstract]
    ]
  }

  if toc {
    let title = if toc_title == none {
      auto
    } else {
      toc_title
    }
    block(above: 0em, below: 2em)[
      #outline(
        title: toc_title,
        depth: toc_depth,
        indent: toc_indent,
      );
    ]
  }

  doc = abstract + doc

  if cols == 1 {
    doc
  } else {
    columns(cols, doc)
  }
}

#show figure: it => {
  place(
    top,
    scope: "parent",
    float: true,
    clearance: 2em,
  )[#it]
}

#set table(
  inset: 6pt,
  stroke: none,
)


// #show bibliography: it => { none }
