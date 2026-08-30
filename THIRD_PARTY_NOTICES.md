# Third-party notices

## Thermochimica SGTE magnetic Gibbs formulation

The SGTE magnetic Gibbs-energy evaluators in
`src/corrosion/base/MSTDBTCData.C` and
`validation/corrosion/calibration/src/msr_corrosion_bv/mstdb.py` are adapted
from `src/gem/CompGibbsMagnetic.f90` in
[ORNL-CEES/Thermochimica](https://github.com/ORNL-CEES/thermochimica), pinned
for this work at commit
[`0c35c8d7d1cf2084b4e2ca5d6608f7dcdf60adad`](https://github.com/ORNL-CEES/thermochimica/commit/0c35c8d7d1cf2084b4e2ca5d6608f7dcdf60adad).
Thermochimica is licensed under the BSD 3-Clause License reproduced below.
Thermochimica itself is not vendored by OpenPronghorn.

```text
BSD 3-Clause License

Copyright 2017 UT-Battelle, LLC
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.
* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from this
  software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

## MSTDB-TC data

No authentic MSTDB-TC database file or raw coefficient record is included in
OpenPronghorn. Tiny test-only ChemSage-format fixtures use chemical species
identifiers with invented thermodynamic coefficients; their filenames and/or
first `System` line mark them as synthetic. Checked calibration reports and an
environment-gated integration test record a finite set of scalar outputs from
an authorized V3.1 run as regression anchors; they contain no database records.
Authorized users supply MSTDB-TC database files at runtime, and the
implementation records and verifies their edition and SHA-256 digests.
