// Minimalist report template.
//
//   #import "template.typ": report
//   #show: report.with(title: "...", authors: ((name: "...", id: "..."),))

#let accent = rgb("#1f3a5f")
#let muted  = luma(110)
#let rule   = luma(200)

#let report(
  title: none,
  subtitle: none,
  authors: (),
  course: none,
  institution: none,
  date: datetime.today(),
  abstract: none,
  body,
) = {
  set document(title: title, author: authors.map(a => a.name))

  set page(
    paper: "a4",
    margin: (x: 2.0cm, top: 2.2cm, bottom: 2.2cm),
    header: context {
      if counter(page).get().first() > 1 {
        set text(8pt, fill: muted)
        smallcaps(title)
        h(1fr)
        course
        v(-0.4em)
        line(length: 100%, stroke: 0.4pt + rule)
      }
    },
    footer: context {
      set align(center)
      set text(8pt, fill: muted)
      counter(page).display("1")
    },
  )

  set text(font: "New Computer Modern", size: 10.5pt, lang: "en")
  set par(justify: true, leading: 0.62em, spacing: 1.1em)

  // Headings: numbered, compact, accent colour on the number only.
  set heading(numbering: "1.1")
  show heading: it => {
    set text(weight: "bold")
    block(above: 1.6em, below: 0.9em, {
      if it.numbering != none {
        text(fill: accent, counter(heading).display(it.numbering))
        h(0.6em)
      }
      it.body
    })
  }
  show heading.where(level: 1): set text(13pt)
  show heading.where(level: 2): set text(11pt)
  show heading.where(level: 3): set text(10.5pt, style: "italic")

  // Code: monospace, light left rule for blocks.
  show raw: set text(font: "DejaVu Sans Mono", size: 8.8pt)
  show raw.where(block: true): it => block(
    width: 100%,
    breakable: false,
    inset: (left: 0.9em, y: 0.5em),
    stroke: (left: 1.5pt + rule),
    it,
  )

  // Tables: booktabs style, horizontal rules only.
  set table(
    stroke: (_, y) => (
      top: if y == 0 { 0.8pt } else if y == 1 { 0.4pt } else { 0pt },
    ),
    inset: (x: 0.6em, y: 0.45em),
  )
  show table: it => block(stroke: (bottom: 0.8pt), it)
  show table: set par(justify: false)
  show table.cell.where(y: 0): set text(weight: "bold")

  // Figures: small caption, label in accent colour.
  set figure(gap: 0.8em)
  show figure.caption: it => {
    set text(9pt)
    text(fill: accent, weight: "bold", it.supplement + " " + context it.counter.display(it.numbering) + ".")
    h(0.3em)
    it.body
  }
  show figure.where(kind: table): set figure.caption(position: top)

  set math.equation(numbering: "(1)")
  show link: set text(fill: accent)
  show ref: set text(fill: accent)
  set list(indent: 0.6em)
  set enum(indent: 0.6em)

  // Title block.
  v(0.3cm)
  {
    set align(left)
    text(20pt, weight: "bold", title)
    if subtitle != none { v(0.1em); align(text(12pt, fill: muted, subtitle), center); v(0.2em) }
    if institution != none { align(text(9pt, fill: muted, institution), center); v(0.2em)}
    v(0.4em)
    line(length: 100%, stroke: 0.6pt + accent)
    v(0.6em)
    grid(
      columns: (1fr, auto),
      align: (left + top, right + top),
      authors.map(a => [#a.name #if "id" in a { text(fill: muted, "(" + a.id + ")") }]).join(linebreak()),
      text(fill: muted)[#course \ #date.display("[day] [month repr:long] [year]")],
    )
  }

  if abstract != none {
    v(1.2em)
    block(inset: (x: 1.5em), {
      set text(9.5pt)
      text(weight: "bold", "Abstract.")
      h(0.4em)
      abstract
    })
  }

  v(1em)
  body
}

#let panel = rgb("#eef2f7")
#let warn-bg = rgb("#fff3cd")
#let warn-border = rgb("#d4a017")

// Callout with a coloured left rule and a bold label line.
#let note-box(
  body,
  label: "Note",
  bg: warn-bg,
  border: warn-border,
) = block(
  width: 100%,
  fill: bg,
  stroke: (left: 3pt + border),
  inset: (x: 10pt, y: 8pt),
  radius: 2pt,
  {
    text(weight: "bold", fill: border.darken(30%), label)
    parbreak()
    body
  },
)

// Numbered stage of a process: a filled circle with the number, then a
// bold title and its body.
#let step-block(
  num,
  title,
  body,
) = block(
  width: 100%,
  fill: panel,
  stroke: (left: 3pt + accent),
  inset: (x: 10pt, y: 8pt),
  radius: 2pt,
  breakable: false,
  grid(
    columns: (22pt, 1fr),
    column-gutter: 10pt,
    align: (center + top, left + top),
    box(
      width: 22pt,
      height: 22pt,
      radius: 11pt,
      fill: accent,
      align(center + horizon, text(fill: white, weight: "bold", str(num))),
    ),
    {
      text(weight: "bold", fill: accent, title)
      parbreak()
      body
    },
  ),
)

// Placeholder for content still to be written; visible in the PDF.
#let todo(body) = text(fill: rgb("#b03a2e"), [\[TODO: #body\]])
