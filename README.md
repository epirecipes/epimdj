# epimdj

**Work in progress!**

This repository contains Julia ports of the R code from Ottar Bjornstad's [Epidemics: Models and Data Using R, Second Edition](https://link.springer.com/book/10.1007/978-3-031-12056-5).

Each chapter is in a separate folder, as a single [Quarto notebook](https://quarto.org/). These can be rendered using the Quarto command line interface. For example, to render `chapter2.qmd` to Markdown, use the following. 

```bash
quarto render chapter2.qmd --to markdown
```

Style-wise, the code is intended to match the R syntax as closely as possible, and so is not idiomatic Julia. Where there is a clear advantage to using Julia-specific packages and APIs, then this will be given in addition.

## Pluto.jl port of Shiny apps

Ports of the accompanying Shiny apps at [https://github.com/objornstad/ecomodelmarkdowns](https://github.com/objornstad/ecomodelmarkdowns) using [Pluto.jl](https://plutojl.org/) can be found in the pluto directory.
