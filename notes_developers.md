## Notes for developers

Recommended steps using Rstudio (requires to set up `roxygen2` & RStudio)

- Build documentation (`command + shift + D`)
- Run the tests (`command + shift + T`)
- Check everything (`command + shift + E`), which also run tests

To load the package (`command + shift + D`)

The R CMD check report should look like this:

```
── R CMD check results ───────────────────────────────── MotrpacBicQC 0.1.0 ────
Duration: 41.6s

0 errors ✓ | 0 warnings ✓ | 0 notes ✓

R CMD check succeeded
```
