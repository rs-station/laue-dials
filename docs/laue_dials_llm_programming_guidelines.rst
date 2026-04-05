.. _laue_dials_llm_programming_guidelines:

======================================
laue-dials LLM Programming Guidelines
======================================

Process
=======

Before Writing Code
-------------------

- Pose clarifying questions to grasp the objective completely. If
  ambiguous, request specifics rather than making assumptions.
- Study existing code before making changes. Understand the function's
  purpose, its callers, and its test coverage.
- Create a plan addressing:

  - Problem statement
  - Overall strategy
  - Detailed step-by-step execution
  - Unit tests for all new and modified code
  - Potentially affected existing tests

- Present the plan for feedback and begin implementation only after
  agreement.

While Writing Code
------------------

- Divide implementation into distinct steps. Present the current state
  after each step for review.
- Run ``tox`` (or ``tox -e default``) after each step. Resolve failures
  before continuing.
- **Perform self-review before presenting each step.** Ask yourself:

  - Did I address edge cases? (None, empty arrays, unindexed reflections,
    missing flex columns, wrong types)
  - Does this work on macOS and Linux? (The supported platforms per
    ``setup.cfg``.)
  - Are there scoping problems?
  - Does this match surrounding code style? (Indentation, naming, PHIL
    conventions, logging style, error handling)
  - Would I suggest further changes if asked?

- Upon completion, perform comprehensive verification:

  - Parse-check all modified files (``python -m py_compile <file>``)
  - Run ``tox`` to execute the full test suite and linters
  - Check for line-length violations (``tox -e lint``)
  - Confirm imports work correctly

After Writing Code
------------------

- Update docstrings and ``docs/cli/`` pages for any user-visible changes.
- If a new command-line script was added, verify its entry point is
  registered in ``setup.cfg`` and that ``tox`` installs and discovers it
  correctly.

Parameter System (PHIL)
========================

All command-line parameters are defined using ``libtbx.phil``. This is
non-negotiable for any ``command_line/`` script; it is what allows users
to pass ``param=value`` on the command line, write override files, and
get diff-phil logging of non-default values. Never use ``argparse``
directly in a ``command_line/`` script.

Anatomy of a Command-Line Script
---------------------------------

Every ``command_line/`` file follows this pattern:

.. code-block:: python

    import libtbx.phil
    from dials.util import show_mail_handle_errors
    from dials.util.options import (
        ArgumentParser,
        reflections_and_experiments_from_files,
    )

    phil_scope = libtbx.phil.parse("""
    output {
      experiments = 'output.expt'
        .type = str
        .help = "Output experiment list filename."
      log = 'laue.myscript.log'
        .type = str
        .help = "Log filename."
    }
    """)

    @show_mail_handle_errors()
    def run(args=None, *, phil=phil_scope):
        parser = ArgumentParser(
            usage="laue.myscript [options] input.expt input.refl",
            phil=phil,
            read_reflections=True,
            read_experiments=True,
            check_format=False,   # True only when pixel data is needed
            epilog=help_message,
        )
        params, options = parser.parse_args(args=args, show_diff_phil=True)
        ...

    if __name__ == "__main__":
        run()

``show_mail_handle_errors()`` must wrap every ``run()`` function; it
provides consistent DIALS error reporting.

Do Not Invent Parameters
------------------------

LLMs frequently hallucinate PHIL parameter names or DIALS API keyword
arguments that do not exist.

- Before using a PHIL parameter, verify it exists in the relevant
  program's ``phil_scope`` definition or in the DIALS source. Do not
  guess names from context.
- Before calling a DIALS or dxtbx Python API, check the actual function
  signature in the source. Do not assume keyword arguments from the
  function name.
- When uncertain, ask: "Does ``lam_min`` exist in this PHIL scope?"
  rather than assuming.

PHIL naming conventions:

- Use ``output.experiments``, ``output.reflections``, and ``output.log``
  for output filenames.
- Use ``include scope <other_phil>`` to compose PHIL scopes from DIALS
  programs rather than duplicating parameter definitions.
- Use ``*`` only for ``choice`` types to mark the default; for
  ``str``/``int``/``float`` just supply the default value directly.
- Always provide a ``.help`` string for every parameter. Multi-line help
  strings use ``\`` for continuation — do not use triple-quoted strings:

  .. code-block:: none

      lam_min = 1.0
        .type = float
        .help = "Minimum wavelength of the beam bandpass in Ångströms. \
                 Reflections assigned wavelengths below this value are \
                 rejected."

- Valid ``.type`` values: ``int``, ``float``, ``str``, ``bool``,
  ``path``, ``choice``, ``strings``, ``ints``, ``floats``.
- To override defaults from an included scope, construct a second PHIL
  string and fetch it into the working scope. See the
  ``indexer_phil`` / ``refiner_phil`` pattern in ``index.py``.

Code Style
==========

Formatting
----------

- Line length: **88 characters** (Black-compatible; enforced by flake8
  in ``setup.cfg``).
- Indentation: **4 spaces**.
- Formatting is enforced by **Black** and **isort** via pre-commit. Run
  ``tox -e lint`` before submitting changes.

Naming
------

- Function and class names should be descriptive: ``assign_wavelengths()``
  not ``aw()``.
- Private helper functions use a leading underscore: ``_compute_residual()``.
- Constants use ``UPPER_SNAKE_CASE``: ``MAX_SPOTS_PER_IMAGE = 10000``.
- Boolean variables should read naturally as assertions: ``is_harmonic``,
  ``has_wavelength``, ``use_scan_varying``.
- Follow established crystallographic naming conventions throughout:
  use ``d_min`` not ``resolution``, ``lam_min``/``lam_max`` not
  ``wavelength_min``/``wavelength_max``, ``s0``/``s1`` for beam
  direction vectors, ``miller_index`` or ``hkl`` not ``reflection_index``.

Cleanup
-------

- Remove trailing whitespace from all lines (enforced by pre-commit).
- Eliminate unused imports entirely. If an import is flagged by the
  linter but is genuinely needed (e.g. re-exported or used in a string
  context), suppress the warning with ``# noqa`` and a brief comment
  explaining why.

Docstrings
----------

All public functions and classes must have docstrings. Use the
**NumPy/Google hybrid** style already present in the codebase:

.. code-block:: python

    def my_function(expts, refls, dmin):
        """
        One-line summary.

        Longer description if needed.

        Args:
            expts (dxtbx.model.ExperimentList): Description.
            refls (dials.array_family.flex.reflection_table): Description.
            dmin (float): Maximum resolution in Ångströms.

        Returns:
            new_refls (dials.array_family.flex.reflection_table): Description.
        """

Always include the full type path in ``Args`` (e.g.
``dials.array_family.flex.reflection_table``, not just ``table``).

For methods that intentionally swallow exceptions, document this with a
``Never raises:`` note explaining the rationale:

.. code-block:: python

    def attempt_wavelength_assignment(refl):
        """
        Attempt to assign a wavelength to a reflection.

        ...

        Never raises:
            Unindexed reflections (id == -1) have no wavelength and are
            silently skipped; callers handle missing assignments downstream.
        """

Type Hints
----------

Type hints are optional and their use varies across the codebase. Match
the convention of the file being edited: if it uses hints, add them to
new functions; if not, do not introduce them. When adding hints, import
from ``typing`` for broad Python compatibility (``Optional``, ``List``,
``Tuple``). Annotate function signatures where helpful, but avoid
annotating every local variable.

Logging
-------

- Use the standard ``logging`` module. Obtain a logger at module level:

  .. code-block:: python

      import logging
      logger = logging.getLogger("laue-dials.command_line.myscript")

- In command-line scripts, attach both a ``FileHandler`` (writing to
  ``params.output.log``) and a ``StreamHandler`` (stdout). Mirror log
  handlers to ``dials``, ``dxtbx``, and ``xfel`` loggers so DIALS output
  is captured. See ``index.py`` for the canonical logging setup block.
- Log section headers with ``"*" * 80`` banners and timing with
  ``time.time()`` differences.
- Use ``logger.info(msg, *args)`` (lazy formatting), not
  ``logger.info(f"... {val}")`` (eager formatting).

Error Handling
--------------

- Use ``libtbx.utils.Sorry`` for user-facing errors (bad input, missing
  files, invalid parameters). Use standard Python exceptions for
  programming errors:

  .. code-block:: python

      from libtbx.utils import Sorry

      # User error — shown without traceback
      if lam_min >= lam_max:
          raise Sorry(
              "lam_min must be less than lam_max "
              f"(got {lam_min} >= {lam_max})"
          )

      # Programming error — should never happen
      assert refls is not None, "refls must be set before integration"

- Never use bare ``except:``; catch at minimum ``except Exception``.
- Do not raise exceptions inside ``run()`` for expected user errors
  (missing files, bad parameters); let the PHIL/ArgumentParser machinery
  handle these, or log a clear message and return early.

Imports
-------

- Standard library imports first, then third-party, then local
  (``laue_dials.*``). isort enforces this.
- Defer heavy DIALS internal imports to inside functions when they are
  only needed in a specific code path. This reduces startup time and
  avoids circular imports. See ``gen_beam_models`` in ``laue.py`` for an
  example.

Testing
=======

Framework
---------

Tests use **pytest**. Test files live under ``tests/algorithms/`` and
mirror the structure of ``src/laue_dials/algorithms/``.

Writing Tests
-------------

- Name test functions ``test_<function_name>`` or
  ``test_<class_name>_<method>``.
- Use ``np.allclose`` for floating-point comparisons; never ``==``.
- For crystallographic correctness checks, assert the mathematical
  property explicitly rather than comparing to a hard-coded numeric
  result.

  .. code-block:: python

      # Good
      assert np.allclose(U @ U.T, np.eye(3))

      # Avoid
      assert np.allclose(U[0, 0], 0.57735)

- Construct minimal synthetic inputs (small numpy arrays, toy unit cells)
  rather than depending on real diffraction data files where possible.
- When tests do require data files, place them under ``tests/data/`` and
  load them with paths relative to ``__file__``.

Running Tests
-------------

Use ``tox`` to run tests in the same isolated environment that CI uses:

.. code-block:: bash

    tox                                          # full test suite (default env)
    tox -- tests/algorithms/test_diffgeo.py      # single file via posargs
    tox -- -k test_hkl2ray                       # single test by name

``tox`` installs DIALS via conda-forge and the package under test before
running pytest, so results match CI exactly. Direct ``pytest`` invocations
work only if DIALS is already on the active ``PATH``.

Adding a New Command-Line Script
=================================

1. Create ``src/laue_dials/command_line/myscript.py`` using the anatomy
   described above.
2. Register the entry point in ``setup.cfg`` under
   ``[options.entry_points]``:

   .. code-block:: ini

       laue.myscript = laue_dials.command_line.myscript:run

3. Keep algorithm logic in ``src/laue_dials/algorithms/``, not inline in
   the command-line script. The command-line file should only handle
   argument parsing, logging setup, file I/O, and timing.
4. Add unit tests to ``tests/algorithms/`` for the new algorithm module.
5. Add a documentation page to ``docs/cli/``.

Adding a New Algorithm
======================

1. Place the implementation in ``src/laue_dials/algorithms/``.
2. Pure numerical functions should not import ``libtbx`` or ``dials``;
   reserve those imports for the wrapping functions that connect to DIALS
   data types.
3. Vectorise over reflections using ``numpy`` wherever possible. Avoid
   ``for refl in refls.rows()`` loops at Python level; these are very
   slow for large datasets.
4. Accept both ``ExperimentList``/``reflection_table`` and plain numpy
   arrays at the algorithm level where it simplifies testing.
5. Document units (e.g. Ångströms for wavelengths and resolution,
   degrees for angles) explicitly in every docstring parameter that
   carries a physical unit.

Common Pitfalls
===============

1. **Modifying a reflection table in-place accidentally.** ``flex``
   arrays are reference types. Always call ``refls.copy()`` before
   modifying columns if the original must be preserved.

2. **Forgetting** ``check_format=True`` **when pixel data is needed.**
   ``ArgumentParser(..., check_format=False)`` is faster but prevents
   reading raw image pixels. Scripts that only manipulate ``.expt``/
   ``.refl`` files should use ``check_format=False``.

3. **Redefining PHIL parameters already in an included scope.** When
   using ``include scope``, all parameters from the included scope are
   present. Override defaults by constructing a second PHIL string and
   fetching it into the working scope; do not redefine the parameter.

4. **Using scalar values where arrays are expected.** Many algorithm
   functions take ``numpy`` arrays; passing a scalar will silently
   broadcast and may produce wrong shapes.

5. **Path separators.** Use ``os.path.join()`` for constructing file
   paths; never hardcode ``/`` or ``\``.

6. **None-safety on** ``dict.get()``. ``d.get("key", "")`` returns the
   default ``""`` only when the key is *absent*; it returns ``None`` when
   the key *exists* with value ``None``. Calling ``.lower()`` or
   ``.strip()`` on the result then crashes:

   .. code-block:: python

       # Wrong — crashes when value is None
       name = params.get("output_file", "").lower()

       # Right — handles both absent and None
       name = (params.get("output_file") or "").lower()

   This pattern arises wherever PHIL scope extracts or dicts from JSON
   are processed.

7. **Floating-point comparison with** ``==``. Never use ``==`` to compare
   floating-point values from crystallographic calculations, including
   wavelengths, resolution limits, and scattering vector components:

   .. code-block:: python

       # Wrong
       assert lam == 1.0

       # Right
       assert np.isclose(lam, 1.0)

   In tests, prefer ``np.allclose`` for arrays and ``np.isclose`` for
   scalars. This applies equally in algorithm code, not just tests.

8. **flex boolean indexing and reduction.** ``flex`` arrays do not support
   NumPy-style boolean indexing or method-style reductions:

   .. code-block:: python

       # Wrong — flex does not support [] with a boolean array
       subset = refls[mask]

       # Right
       subset = refls.select(flex.bool(mask))

       # Wrong — flex arrays have no .sum() method
       total = arr.sum()

       # Right
       total = flex.sum(arr)

   When converting between flex and NumPy, use ``.as_numpy_array()`` and
   the appropriate flex constructor (``flex.double``, ``flex.int``,
   ``flex.bool``).

9. **Space group comparison.** Never compare ``gemmi.SpaceGroup`` objects
   as strings — different notations represent the same group:

   .. code-block:: python

       # Wrong — "P 21 21 21" and "P212121" are the same group
       if sg.xhm() == "P 21 21 21":
           ...

       # Right
       if sg == gemmi.SpaceGroup("P 21 21 21"):
           ...

   Use ``rs.utils.is_absent(hkl, spacegroup)`` from
   ``reciprocalspaceship`` for systematic absence checks; do not
   reimplement them by hand.

10. **Crystallographic unit conventions.** Wrong units produce silently
    incorrect results with no error. The conventions used throughout
    this codebase are:

    - Wavelengths: **Ångströms** (not nm or keV)
    - Resolution / ``d_min``: **Ångströms**
    - Angles: **degrees** (PHIL inputs) or **radians** (internal
      rotation matrices — check the context)
    - Scattering vectors ``s0``, ``s1``: **unit vectors** (normalise
      before use)
    - ``1/d²``: used for KDE training in ``outliers.py`` — do not
      confuse with ``d`` itself

Serialization
=============

When writing classes that serialize to or from dicts (e.g. for JSON
persistence or inter-process communication):

- Implement ``to_dict()`` / ``from_dict()`` and verify they survive a
  round-trip.
- ``from_dict()`` must tolerate missing keys so that data serialized by
  an older version of the code can still be loaded:

  .. code-block:: python

      @classmethod
      def from_dict(cls, d):
          obj = cls()
          # Use .get() with defaults, never d["key"]
          obj.lam_min = d.get("lam_min", 1.0)
          obj.lam_max = d.get("lam_max", 1.5)
          return obj

- Test the round-trip explicitly:

  .. code-block:: python

      def test_serialization_round_trip():
          original = MyClass(lam_min=1.0, lam_max=1.5)
          restored = MyClass.from_dict(original.to_dict())
          assert np.isclose(restored.lam_min, original.lam_min)
          assert np.isclose(restored.lam_max, original.lam_max)

Code Presentation
=================

When presenting code changes, avoid placeholder comments such as
``# ... rest of function unchanged ...`` inside code blocks, as output
is often pasted directly. For small changes, show the entire modified
function. For large files, clearly indicate the exact old and new text
with surrounding context. If showing only part of a file, mark the
boundaries explicitly with line numbers or ellipsis comments.

Checklist Before Presenting Code
=================================

- [ ] All changed files parse cleanly (``python -m py_compile <file>``)
- [ ] ``tox`` passes (full test suite)
- [ ] ``tox -e lint`` passes (Black, isort, flake8)
- [ ] No bare ``except:`` blocks added
- [ ] Any intentionally swallowed exceptions have a ``Never raises:``
      docstring contract explaining the rationale
- [ ] No hardcoded path separators — use ``os.path.join()``
- [ ] None-safety: ``(x.get("key") or "")`` not ``x.get("key", "")``
- [ ] No ``==`` comparisons on floating-point crystallographic values —
      use ``np.isclose`` / ``np.allclose``
- [ ] flex boolean indexing uses ``.select(flex.bool(...))`` not ``[]``
- [ ] All PHIL parameters verified against the relevant scope definition —
      none invented from context
- [ ] Physical units (Å, degrees) documented in docstrings for all new
      parameters that carry a unit
- [ ] New public functions and classes have docstrings with full type
      paths in ``Args``
- [ ] If a new command-line script was added, its entry point is
      registered in ``setup.cfg``
- [ ] If a new command-line script was added, a ``docs/cli/`` page exists

Pre-Commit and CI
=================

The repository uses pre-commit hooks for formatting and linting. Run
them via tox before pushing:

.. code-block:: bash

    tox -e lint           # runs pre-commit (Black, isort, flake8, …)
    tox                   # runs the full test suite

Both environments are required to pass before a PR is mergeable. CI
(GitHub Actions) runs both on every pull request. Do not submit PRs that
break existing tests or linting checks.
