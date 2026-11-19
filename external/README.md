# Third-party software bundled with PACT

PACT depends on the work of other authors. We are grateful to them and
redistribute their code unmodified, under their original license terms.
Please cite their work as well as PACT when you publish results.

## CircStat — Circular Statistics Toolbox for MATLAB

* **Folder:** `external/CircStat/`
* **Author:** Philipp Berens (Berens Lab, University of Tübingen)
* **Contributors:** Marc J. Velasco, Tal Krasovsky, and the contributors to
  the GitHub repository.
* **Source:** <https://github.com/circstat/circstat-matlab>
  (the authors' official repository; also distributed on MATLAB Central
  File Exchange, ID 10676).
* **Version bundled:** `master` at commit
  `302acc161d33bdb5c432243b946b178d88c361b8` (24 May 2019), the most recent
  code published by the authors as of November 2026. The repository has no
  tagged releases.
* **License:** BSD (see `external/CircStat/license.txt`, Copyright (c) 2011,
  Philipp Berens). The upstream README states "Distributed under Open Source
  BSD License"; the license text is retained from the authors' File Exchange
  distribution because the GitHub repository does not carry a separate
  license file.
* **Modifications:** none. The files are an exact copy of the upstream
  commit above. (The `examples/` folder and its data files from the earlier
  File Exchange 2012a package are not included; PACT does not use them.)

**Please cite:**

> P. Berens, *CircStat: A MATLAB Toolbox for Circular Statistics*,
> Journal of Statistical Software, 31(10), 2009.
> doi:10.18637/jss.v031.i10

### How PACT uses CircStat

PACT calls only `circ_r`, `circ_mean`, `circ_rad2ang`, `circ_rtest` (z
statistic only, in the HFO-frequency scan) and `circ_plot` (polar angular
histogram). The regression test `code/tests/test_pact_circstat.m` checks
that these functions return the values PACT expects.

### History

Earlier PACT releases shipped CircStat 2012a in a top-level `CircStat2012a/`
folder without acknowledging its authors in the PACT documentation. We
apologize for that omission. In November 2026 the toolbox was moved here,
updated to the authors' latest published code, and credited properly.
