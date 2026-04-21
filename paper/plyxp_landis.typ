// Some definitions presupposed by pandoc's typst output.
#let blockquote(body) = [
  #set text( size: 0.92em )
  #block(inset: (left: 1.5em, top: 0.2em, bottom: 0.2em))[#body]
]

#let horizontalrule = line(start: (25%,0%), end: (75%,0%))

#let endnote(num, contents) = [
  #stack(dir: ltr, spacing: 3pt, super[#num], contents)
]

#show terms: it => {
  it.children
    .map(child => [
      #strong[#child.term]
      #block(inset: (left: 1.5em, top: -0.4em))[#child.description]
      ])
    .join()
}

// Some quarto-specific definitions.

#show raw.where(block: true): set block(
    fill: luma(230),
    width: 100%,
    inset: 8pt,
    radius: 2pt
  )

#let block_with_new_content(old_block, new_content) = {
  let d = (:)
  let fields = old_block.fields()
  fields.remove("body")
  if fields.at("below", default: none) != none {
    // TODO: this is a hack because below is a "synthesized element"
    // according to the experts in the typst discord...
    fields.below = fields.below.abs
  }
  return block.with(..fields)(new_content)
}

#let empty(v) = {
  if type(v) == str {
    // two dollar signs here because we're technically inside
    // a Pandoc template :grimace:
    v.matches(regex("^\\s*$")).at(0, default: none) != none
  } else if type(v) == content {
    if v.at("text", default: none) != none {
      return empty(v.text)
    }
    for child in v.at("children", default: ()) {
      if not empty(child) {
        return false
      }
    }
    return true
  }

}

// Subfloats
// This is a technique that we adapted from https://github.com/tingerrr/subpar/
#let quartosubfloatcounter = counter("quartosubfloatcounter")

#let quarto_super(
  kind: str,
  caption: none,
  label: none,
  supplement: str,
  position: none,
  subrefnumbering: "1a",
  subcapnumbering: "(a)",
  body,
) = {
  context {
    let figcounter = counter(figure.where(kind: kind))
    let n-super = figcounter.get().first() + 1
    set figure.caption(position: position)
    [#figure(
      kind: kind,
      supplement: supplement,
      caption: caption,
      {
        show figure.where(kind: kind): set figure(numbering: _ => numbering(subrefnumbering, n-super, quartosubfloatcounter.get().first() + 1))
        show figure.where(kind: kind): set figure.caption(position: position)

        show figure: it => {
          let num = numbering(subcapnumbering, n-super, quartosubfloatcounter.get().first() + 1)
          show figure.caption: it => {
            num.slice(2) // I don't understand why the numbering contains output that it really shouldn't, but this fixes it shrug?
            [ ]
            it.body
          }

          quartosubfloatcounter.step()
          it
          counter(figure.where(kind: it.kind)).update(n => n - 1)
        }

        quartosubfloatcounter.update(0)
        body
      }
    )#label]
  }
}

// callout rendering
// this is a figure show rule because callouts are crossreferenceable
#show figure: it => {
  if type(it.kind) != str {
    return it
  }
  let kind_match = it.kind.matches(regex("^quarto-callout-(.*)")).at(0, default: none)
  if kind_match == none {
    return it
  }
  let kind = kind_match.captures.at(0, default: "other")
  kind = upper(kind.first()) + kind.slice(1)
  // now we pull apart the callout and reassemble it with the crossref name and counter

  // when we cleanup pandoc's emitted code to avoid spaces this will have to change
  let old_callout = it.body.children.at(1).body.children.at(1)
  let old_title_block = old_callout.body.children.at(0)
  let old_title = old_title_block.body.body.children.at(2)

  // TODO use custom separator if available
  let new_title = if empty(old_title) {
    [#kind #it.counter.display()]
  } else {
    [#kind #it.counter.display(): #old_title]
  }

  let new_title_block = block_with_new_content(
    old_title_block, 
    block_with_new_content(
      old_title_block.body, 
      old_title_block.body.body.children.at(0) +
      old_title_block.body.body.children.at(1) +
      new_title))

  block_with_new_content(old_callout,
    block(below: 0pt, new_title_block) +
    old_callout.body.children.at(1))
}

// 2023-10-09: #fa-icon("fa-info") is not working, so we'll eval "#fa-info()" instead
#let callout(body: [], title: "Callout", background_color: rgb("#dddddd"), icon: none, icon_color: black, body_background_color: white) = {
  block(
    breakable: false, 
    fill: background_color, 
    stroke: (paint: icon_color, thickness: 0.5pt, cap: "round"), 
    width: 100%, 
    radius: 2pt,
    block(
      inset: 1pt,
      width: 100%, 
      below: 0pt, 
      block(
        fill: background_color, 
        width: 100%, 
        inset: 8pt)[#text(icon_color, weight: 900)[#icon] #title]) +
      if(body != []){
        block(
          inset: 1pt, 
          width: 100%, 
          block(fill: body_background_color, width: 100%, inset: 8pt, body))
      }
    )
}


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
      [#author.name]
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
      n += 1
      if not first {
        [,   ]
      }
      super[#n] + [#key]
      first = false
    }
    // let ncols = calc.min(count, 3)
    // grid(
    //   columns: (1fr,) * ncols,
    //   row-gutter: 1.5em,
    //   ..authors.map(author => align(center)[
    //     #author.name \
    //     #author.affiliation \
    //     #author.email
    //   ])
    // )
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


// Typst custom formats typically consist of a 'typst-template.typ' (which is
// the source code for a typst template) and a 'typst-show.typ' which calls the
// template's function (forwarding Pandoc metadata values as required)
//
// This is an example 'typst-show.typ' file (based on the default template
// that ships with Quarto). It calls the typst function named 'article' which
// is defined in the 'typst-template.typ' file.
//
// If you are creating or packaging a custom typst template you will likely
// want to replace this file and 'typst-template.typ' entirely. You can find
// documentation on creating typst templates here and some examples here:
//   - https://typst.app/docs/tutorial/making-a-template/
//   - https://github.com/typst/templates
#show: doc => article(
  title: [Efficient and Tidy Manipulation of Annotated Matrix Data with plyxp],
  authors: (
    ( name: [Justin T. Landis],
      affiliation: [Bioinformatics and Computation Biology Program||| Department of Genetics],
      email: [] ),
    ( name: [Michael I. Love],
      affiliation: [Department of Genetics||| Department of Biostatistics],
      email: [] ),
    ),
  abstract: [In order to manipulate high-dimensional omics data, such as bulk or single cell gene expression counts matrices, a bioinformatics analyst must typically learn domain-specific functions and syntax. These matrix-centric functions and syntax may not be as straightforward as working with tidy data analytic principles using tools such as dplyr applied to tabular data. We propose an expressive grammar for manipulating annotated matrix data, with syntax to compute on and update matrix data and tabular row and column metadata, including row-wise or column-wise grouped operations. This grammar works by defining multiple contexts, and providing pronouns for specific recall and assignment in these contexts. The plyxp package is implementation of this grammar for the R/Bioconductor ecosystem, with efficient abstractions for the SummarizedExperiment class. We demonstrate plyxp's efficiency compared to alternative approaches on data manipulation tasks requiring computation across contexts.

],
  abstract-title: "Abstract",
  pagenumbering: "1",
  toc_title: [Table of contents],
  toc_depth: 3,
  cols: 2,
  doc,
)


#show raw.where(block: true): it => {
  show raw.line: it => {
    text(fill: gray)[#it.number]
    h(1em)
    it.body
  }
  set text(size: 6pt);
  it
}

= Introduction
<introduction>
High-throughput omics data are often summarized into one or more matrices representing counts or relative abundance of molecules in a set of experiments, as summarized across various genomic features (rows), and across different experiments or samples (columns). During subsequent analysis, it is critical to preserve links between rows and columns of the matrix data with relevant feature and sample metadata tables. For example, one may wish to divide a matrix of RNA-seq transcript counts by the transcript length (metadata associated with the transcript features along the rows) or by the library size (metadata associated with samples along the columns). This need for robust and efficient handling of annotated matrix data motivated the development of various data classes, including but not limited to the #emph[ExpressionSet] @Gentleman2004 and the #emph[SummarizedExperiment] @Huber2015 in R/Bioconductor, and #emph[AnnData] @Virshup2024 in Python. Such annotated matrix data classes and associated methods ensure that structural consistency is maintained between data and metadata during common operations such as subsetting, re-ordering, and row or column binding.

Operating on these data objects requires learning domain-specific functions and syntax. Performing complex operations such as grouping over sets of rows and columns and computing further data summarization may require writing additional control structures such as for loops or using less well-known functional programming helpers, which can obfuscate the purpose of analysis code. The #emph[tidyomics] project was created within R/Bioconductor to help address this issue @Hutchison2024, providing methods and coding paradigms for working with common R data classes representing genomic ranges @Lee2019@Lee2020, count matrices @Mangiola2021bulk, and Seurat objects @Mangiola2021Seurat. The project enables intuitive manipulation and visualization without requiring users to learn domain-specific syntax. #emph[tidyomics] packages and workflows enable chained data-manipulative operations that resemble the style and framework of tidy data analysis @Wickham2014 as implemented in #emph[dplyr] @dplyr and other #emph[tidyverse] @tidyverse packages, or chained operations for python #emph[DataFrames] using #emph[pandas] @McKinney2022. Due to the popularity of the #emph[tidyverse] among newer R users, #emph[tidyomics] increases the accessibility of omic data analysis.

Here, we propose a new expressive grammar for access and manipulation of annotated matrix data and a framework for efficiently performing the operations implied by the syntax. Data within the annotated matrix object are lazily bound to a series of environments that are evaluated when the user forces their symbols. Thus, a user has more freedom in how they choose to work with their data. We introduce the #emph[plyxp] R package that implements this grammar, allowing the use of #emph[dplyr] syntax with the Bioconductor #emph[SummarizedExperiment] (SE) data class.

The #emph[tidySummarizedExperiment] package, which was released with Bioconductor 3.12 in 2020, also provides #emph[dplyr];-like access to SE objects within the #emph[tidyomics] project. The #emph[tidySummarizedExperiment] package allows users to approach the SE as an abstracted data frame, conceptually unwrapping the compact matrix representation into a tidy, long-form version. #emph[tidySummarizedExperiment] facilitates data exploration and visualization, allowing datasets to be directly piped into #emph[ggplot] plotting functions.

However some complex computations involving simultaneous manipulation of matrix values and annotation metadata were slower than necessary, a limitation addressable through data-masking strategies.

#emph[plyxp];'s use of data masking from the #emph[rlang] @rlang package provides gains in efficiency and is currently used to speed up existing functionality within #emph[tidySummarizedExperiment];.

Further, #emph[plyxp] reinforces the structure of the annotated matrix through additional syntax compared to #emph[tidySummarizedExperiment];, often resulting in more verbose, but more specific evaluation.

= Syntax
<syntax>
We describe a grammar in which users can access, modify, and append data across different contexts of a dataset, where here we mean "context" to refer to a collection of objects of the same length (for vectors) or dimension (for matrices), which collectively describe either features, samples, or a matrix of values over features and samples. For annotated matrix data with $n$ rows and $p$ columns, we define three contexts: the `assays` consist of one or more matrices of size $n times p$; the `rows` consist of row metadata (a collection of vectors of length $n$ describing the features); and the `cols` consist of column metadata (a collection of vectors of length $p$ describing the samples/cells). For example, an object representing an RNA-seq experiment might have an `assays` context consisting of a counts matrix and a transcripts-per-million (TPM) matrix across $n$ genes and $p$ samples, a `rows` context with feature metadata such as alternative gene identifiers and gene length, and a `cols` context with sample metadata such as library size and sample condition group.

By default, expressions are evaluated against the named symbols within the `assays` context. The user can specify to evaluate against named symbols in the row or column metadata by decorating their expressions with sentinel functions, #raw(lang:"r", "rows(...)") or #raw(lang:"r", "cols(...)");. Beyond communicating where expressions are evaluated, this syntax also enforces where the resulting vector or matrix will be stored. For example, with the #raw(lang:"r", "mutate()") command, the new objects are bound to their respective context within the object, but for the #raw(lang:"r", "filter()") command, instead of adding data, these results inform how the rows or columns should be subset.

```r
data |>
  mutate(
    log10_counts = log10(counts + 1),
    cols(treated = (dex == "trt"))
  )
```

```r
data |>
  filter(
    rows(gene_biotype == "lincRNA"),
    cols(avgLength >= 126)
  )
```

The above example demonstrates this grammar as implemented in #emph[plyxp];.

The user computes the logarithm of the counts matrix plus a pseudocount, and creates a new logical variable in the column metadata, each of which are stored in their respective contexts. In the next step of the pipe, rows and columns are filtered based on available metadata about the genes and the samples. Here the base R `|>` symbol refers to piping the results of a function into the first place argument of the following function. As in #emph[dplyr];, the `mutate()` verb implies an addition of an object (matrix or vector) to the corresponding evaluation context, whereas `filter()` subsets the data (described more later).

The #emph[plyxp] syntax is designed to be convenient and flexible, providing users with multiple tools to access or manipulate their data across contexts. For example, consider the following:

```r
data |>
  mutate(
    norm_cts = counts/.cols$size_factor,
    cts_bp = counts/.rows$gene_length
  )
```

Using pronouns `.cols` and `.rows`, the user accesses vectors from another context. A "pronoun" is a special object of the #emph[rlang] package that provides an unambiguous symbol lookup into named lists or environments. In the above expression, an element-wise division of the counts matrix via data within the columns or rows is done by replicating the data to the appropriate size automatically.

Alternatively, in the computation of `gene_ave` in the following code chunk, the user employs the efficient row-wise means operator `rowMeans`. Alternatively, for `gene_ave2` the expression `.assays$counts` within the row context is resolved to a list where each element is a row of the counts matrix. A list output for the `.assays` pronoun offers flexibility as lists can be supplied to the apply family of functions or using #emph[purrr] @purrr as shown below.

```r
data |>
  mutate(
    rows(
      gene_ave = rowMeans(.assays_asis$counts),
      gene_ave2 = map_dbl(.assays$counts, mean)
    )
  )
```

In the example code above, the `data` object can be instantiated from any #emph[SummarizedExperiment] using the #raw(lang:"r", "new_plyxp()") function (see Supplementary Note for more details).

= Semantics
<semantics>
#emph[plyxp] is an implementation of this grammar for annotated matrix data stored in R/Bioconductor's #emph[SummarizedExperiment] class, allowing the use of #emph[dplyr] verbs with the following expected outcomes and behaviors. #emph[plyxp] enforces conventions that eliminate ambiguity while providing flexibility in terms of data shape and computation:

+ Expressions are scoped explicitly to an evaluation context.
+ Evaluation contexts are connected via pronouns, which provide both an `asis` and size-stable version of the requested data.
+ #emph[dplyr] verbs (#raw(lang:"r", "mutate()");, #raw(lang:"r", "filter()");, etc.) only support expressions that would result in the same structure of object as input (SE object).

Evaluation contexts and their connection via pronouns are diagrammed in Figure #ref(<fig-plyxp-diagram>, supplement: [Figure]). We note that #emph[plyxp] relies on Bioconductor-native calls, e.g.~subsetting with square brackets `[,]`, and is thus never faster than using the SE interface directly. There will always be some overhead in generating the necessary abstraction to enable implicit references into the object as required by #emph[dplyr];. However, we aim that #emph[plyxp] always be as close to Bioconductor-native speed as possible, via line profiling of the function definitions using #emph[profvis];.

== Pronouns and Data Masks
<pronouns-and-data-masks>
#emph[plyxp] enables expressive evaluation by leveraging #emph[rlang] data masks, a series of basic R environments that are evaluated with the #emph[rlang] tidy evaluation engine @rlang[#cite(<Wickham2019>, form: "prose");]. #emph[plyxp] exploits several utilities from the #emph[rlang] package, such as lazy evaluation and #emph[rlang] pronouns, to manage these environments. Here, "lazy evaluation" refers to expressions whose computation is delayed until the value is accessed. A developer may choose any number of environments for their data mask that are necessary for their use case. With multi-context scenarios that require multiple data masks, as with annotated matrix data, we suggest the final parent ancestor environment of the mask is reserved for binding #emph[rlang] pronouns. This suggestion is motivated by R's lexical scoping, which allows variable bindings to be searched through environment parents during evaluation @rbase. Placing pronouns in the final parent of the mask, #emph[i.g.] the final environment prior to the calling environment, will ensure the user's data within the context are prioritized over potential access to other contexts. The #emph[rlang] pronouns can be conceptually thought of as a pointer to another environment. #emph[rlang] pronouns are utilized within the #emph[dplyr] package to disambiguate where to find a symbol, in either the calling environment or user's data with `.env`and`.data` pronouns respectively.

In #emph[plyxp];, these pronouns are used to connect each context through two pronouns, an `asis` variant and a "size-stable" variant. The `asis`variant is to provide the user with access to the actual data stored on the object, as there may be efficient ways to operate on its native shape #emph[e.g.] the row-wise and column-wise matrix operations in base R and #emph[matrixStats];. The size-stable variant provides the requested data reshaped to match the size of the current context. The size-stable variant is intended to reduce friction from pulling data from other contexts, as there is no guarantee the target data`asis` will match the expected output size of the evaluation context.

With annotated matrix data, there are predictable transformations that may be used to guarantee size stability. For example, consider requesting row or column data from within the `assays` evaluation context. To make the row data size stable, we need only to replicate the target vector as many times as the number of columns in our matrix context, and then concatenate these vectors. For column data within the `assays` evaluation context, we replicate each element of the target vector as many times as the number of rows in our matrix context. Requesting assay data from within the `rows` or `cols` evaluation context, we may reshape the matrix by returning a list where each element is a row or column of the matrix, respectively.

#figure([
#box(image("./application_note_F1_005.png"))
], caption: figure.caption(
position: bottom, 
[
Diagram of the data masking structure in #emph[plyxp];. Grey rectangles represent environments. Parent environments are linked to child environments by dotted lines where the circle indicates the parent and the arrow indicates the child environment. Full environment hierarchies are not explicitly shown. Example variable bindings are shown within environments. Data mask contexts are grouped by color where the `assays`, `rows`, and `cols` contexts are represented by the green, purple, and yellow blocks respectively. Red lines show connections from pronouns of one context to environments that are children of other contexts. For brevity, only extended environments from the `assays` context into the `rows` and `cols` contexts are shown. Entry points for each context is indicated by the blue star. Through #emph[rlang];'s tidy evaluation engine, symbols within the data mask are evaluated before searching the calling environment.
]), 
kind: "quarto-float-fig", 
supplement: "Figure", 
)
<fig-plyxp-diagram>


== Grouping
<grouping>
In #emph[dplyr];, specifying groups with #raw(lang:"r", "group_by()") augments the behavior of further #emph[dplyr] verbs. #emph[dplyr] grouping operations traditionally include additional metadata on the target object describing how the data should be partitioned. This strategy works well for data structures that are rectangular, but may lead to ambiguities for data structures that deviate from this structure, such as annotated matrix data. For tabular data, there is only one dimension (along the rows) in which those data can be partitioned, and no matter how the rows are partitioned, a new table can be constructed from the results. This is not true for matrices. One such example would be to attempt to group by matrix values \> 100. Here, there is no guarantee the partitioning of arbitrary elements of a matrix would create valid sub-matrices. It is for this reason that #emph[plyxp] follows strict heuristics for grouping and only allows groupings that partition rows or columns of the matrix.

This restriction coincides with the natural extension of applying grouping to tabular data from #emph[dplyr];. We may conceptually regard the row data and the column data as independent and grouped as normal tabular structures. The row data and the column data are related to one another through the matrix, and thus the matrix observes all grouping partitions imposed on both dimensions. Since the row data and the column data groupings are independent, the `rows` and `cols` contexts only consider their own respective grouping during evaluation. In other words, the `assays` context will always respect groupings about the `rows` and `cols` contexts, but the `rows` context drops any `cols` grouping and `cols` context drops any `rows` grouping.

= Performance Comparison
<performance-comparison>
All evaluation within #emph[plyxp] is intended to take place within one of the specified contexts mentioned previously. Due to the particular shape of annotated matrix data, each context has a unique notion of size and a unique set of symbols within. Additional syntax is used to declare where an expression evaluates and there is no additional cost to accessing data `asis`. An alternative to this approach would be to force all variables within the same context, such as an instantiated or abstracted #emph[data.frame] where all contexts are joined together.

A single context would require variables within the assays to be unwound into one dimensional vectors, and data within the row and column metadata tables to be replicated appropriately. Additional performance considerations are discussed in more detail in the supplemental note. When comparing between the two approaches, namely #emph[plyxp] compared to annotated matrix data unwound into a #emph[data.frame];, #emph[plyxp] runs $tilde.op 14$ times faster when scaling the counts matrix by gene length using the #emph[airway] dataset @Himes2014.

= Conclusion
<conclusion>
#emph[plyxp] provides a framework for operating on annotated matrix data. As data masking was possible via the #emph[rlang] package, we were able to introduce novel syntax for operating on SE objects within Bioconductor that is ergonomically inspired by #emph[dplyr] style and unambiguously expressive. #emph[tidySummarizedExperiment] has incorporated #emph[plyxp] methods internally for some operations as of Bioconductor 3.22 release in Fall 2025, such that users of both #emph[plyxp] and #emph[tidySummarizedExperiment] benefit from the performance gains described above, as compared to operations that would require creation of a long form tabular representation of the annotated matrix data. We imagine that the grammar described here could be applied to manipulate #emph[AnnData] with #emph[pandas] methods. We plan to release the internals of #emph[plyxp] that construct and manage the data masks so that developers can choose to create their own #emph[dplyr] style for additional R classes.

= Code Availability
<code-availability>
The #emph[plyxp] software package is available from Bioconductor, installable with `BiocManager::install()`, and at #link("https://github.com/jtlandis/plyxp");. All code used in this paper, including comparisons of run time to alternative implementations and reproducible code chunks appearing above, is available at #link("https://github.com/jtlandis/plyxp_paper")

#bibliography("references.bib")

