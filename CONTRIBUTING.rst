============
Contributing
============

Welcome to ``laue-dials`` developer's guide.

This document focuses on outlining development processes, but `other kinds of contributions`_ are also
appreciated.

Please note that all users and contributors are expected to be **open,
considerate, reasonable, honest, and respectful** when contributing
ideas, issues, documentation, or code.

Issue Reports
=============

If you experience bugs or general issues with ``laue-dials``, please have a look
on the `issue tracker`_. If you don't see anything useful there, please feel
free to fire an issue report.

New issue reports should include information about your programming environment
(e.g., operating system, Python version) and steps to reproduce the problem.
Please try also to simplify the reproduction steps to a very minimal example
that still illustrates the problem you are facing. By removing other factors,
you help us to identify the root cause of the issue.


Documentation Improvements
==========================

You can help improve ``laue-dials`` docs by making them more readable and coherent, or
by adding missing information and correcting mistakes.

``laue-dials`` documentation uses Sphinx_ as its main documentation compiler.
This means that the docs are kept in the same repository as the project code, and
that any documentation update is done in the same way was a code contribution.
When working on documentation changes in your local machine, you can
compile them using ``tox``

    tox -e docs

and use Python's built-in web server for a preview in your web browser
(``http://localhost:8000``)::

    python3 -m http.server --directory 'docs/_build/html'


Code Contributions
==================

To the extent possible, this software is to be written such that it functions
as an extension of `DIALS`_ for the user. Software should be written in a modular
form with a command-line interface for common crystallographic tasks needed during
data analysis. The scope of this package is limited to Laue
(i.e. wide spectral-bandwidth) crystallographic experiments, although multiple
types of experiments under this umbrella can be supported.


Submit an issue
---------------

Before you work on any non-trivial code contribution it's best to first create
a report in the `issue tracker`_ to start a discussion on the subject.
This often provides additional considerations and avoids unnecessary work.
Please add an issue even when working on new features that don't involve changes
to existing code, to prevent overlaps in work with others.

Create an environment
---------------------

Before you start coding, we recommend creating an isolated `virtual
environment`_ to avoid any problems with your installed Python packages.
This can easily be done via Miniconda_::

    conda create -n laue-dials python=3 six virtualenv pytest pytest-cov
    conda activate laue-dials

Clone the repository
--------------------

#. Create an user account on |the repository service| if you do not already have one.
#. Fork the project repository_: click on the *Fork* button near the top of the
   page. This creates a copy of the code under your account on |the repository service|.
#. Clone this copy to your local disk::

    git clone git@github.com:YourLogin/laue-dials.git
    cd laue-dials

#. You should run::

    pip install -U pip setuptools -e .

   to be able to import the package under development in the Python REPL.

#. Install |pre-commit|_::

    pip install pre-commit
    pre-commit install

   ``laue-dials`` comes with a lot of hooks configured to automatically help the
   developer to check the code being written. Compliance with all hooks is
   necessary to contribute code to maintain code quality.

Implement your changes
----------------------

#. Create a branch to hold your changes::

    git checkout -b my-feature

   and start making changes. Never work on the main branch!

#. Start your work on this branch. Don't forget to add docstrings_ to new
   functions, modules and classes. Sphinx will automatically build
   documentation for your work.

#. Add yourself to the list of contributors in ``AUTHORS.rst`` if you are
   not already listed.

#. When you’re done editing, do::

    git add <MODIFIED FILES>
    git commit

   to record your changes in git_.

   Please make sure to see the validation messages from |pre-commit|_ and fix
   any eventual issues.
   This should automatically use flake8_/black_ to check/fix the code style
   in a way that is compatible with the project. Any unfixed issues will
   result in a rejected contribution.

   Don't forget to add unit tests and documentation in case your
   contribution adds an additional feature and is not just a bugfix.

   Moreover, writing a `descriptive commit message`_ is mandatory.
   In case of doubt, you can check the commit history with::

      git log --graph --decorate --pretty=oneline --abbrev-commit --all

   to look for recurring communication patterns.

#. Please check that your changes don't break any unit tests with::

    tox

   (after having installed ``tox-conda`` with ``pip install tox-conda`` or ``pipx``).

   You can also use ``tox-conda`` to run several other pre-configured tasks in the
   repository. Try ``tox -av`` to see a list of the available checks.

Submit your contribution
------------------------

#. If everything works fine, push your local branch to |the repository service| with::

    git push -u origin my-feature

#. Go to the web page of your fork and click |contribute button|
   to send your changes for review. Find more detailed information in
   `creating a PR`_.


Troubleshooting
---------------

The following tips can be used when facing problems to build or test the
package:

#. Make sure to fetch all the tags from the upstream repository_.
   The command ``git describe --abbrev=0 --tags`` should return the version you
   are expecting. If you are trying to run CI scripts in a fork repository,
   make sure to push all the tags.
   You can also try to remove all the egg files or the complete egg folder, i.e.,
   ``.eggs``, as well as the ``*.egg-info`` folders in the ``src`` folder or
   potentially in the root of your project.

#. Sometimes ``tox-conda`` misses out when new dependencies are added, especially to
   ``setup.cfg`` and ``docs/requirements.txt``. If you find any problems with
   missing dependencies when running a command with ``tox-conda``, try to recreate the
   ``tox`` environment using the ``-r`` flag. For example, instead of::

    tox -e docs

   Try running::

    tox -r -e docs

#. Make sure to have a reliable ``tox-conda`` installation that uses the correct
   Python version (e.g., 3.7+). When in doubt you can run::

    tox --version
    # OR
    which tox

#. `Pytest can drop you`_ in an interactive session in the case an error occurs.
   In order to do that you need to pass a ``--pdb`` option (for example by
   running ``tox -- -k <NAME OF THE FAILING TEST> --pdb``).
   You can also setup breakpoints manually instead of using the ``--pdb`` option.


Maintainer tasks
================

Releases
--------

If you are part of the group of maintainers and have correct user permissions
on PyPI_, the following steps can be used to release a new version for
``laue-dials``:

#. Make sure all unit tests are successful.
#. Update ``CHANGELOG.rst`` with new features and changes for the new release.
#. Run ``git pull``, resolve any merge conflicts, and then ``git push`` the source code.
#. Tag the current commit on the main branch with a release tag, e.g., ``git tag -a v0.1 -m 'Version message'``.
#. Push the new tag to the upstream repository_, e.g., ``git push origin v0.1``
#. Navigate to ``https://github.com/rs-station/laue-dials/releases/new``.
#. Select the appropriate tag and write a description for the release.
#. Set as a pre-release if necessary, and then publish the release on Github.
#. After Github Actions workflows have executed, check PyPI to ensure they worked correctly.

.. <-- Documentation variables -->
.. _repository: https://github.com/rs-station/laue-dials
.. _issue tracker: https://github.com/rs-station/laue-dials/issues

.. |the repository service| replace:: GitHub
.. |contribute button| replace:: "Create pull request"
.. |virtualenv| replace:: ``virtualenv``
.. |pre-commit| replace:: ``pre-commit``


.. _black: https://pypi.org/project/black/
.. _CommonMark: https://commonmark.org/
.. _contribution-guide.org: https://www.contribution-guide.org/
.. _creating a PR: https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request
.. _descriptive commit message: https://chris.beams.io/posts/git-commit
.. _DIALS: https://dials.github.io/index.html
.. _docstrings: https://www.sphinx-doc.org/en/master/usage/extensions/napoleon.html
.. _first-contributions tutorial: https://github.com/firstcontributions/first-contributions
.. _flake8: https://flake8.pycqa.org/en/stable/
.. _git: https://git-scm.com
.. _Miniconda: https://docs.conda.io/en/latest/miniconda.html
.. _MyST: https://myst-parser.readthedocs.io/en/latest/syntax/syntax.html
.. _other kinds of contributions: https://opensource.guide/how-to-contribute
.. _pre-commit: https://pre-commit.com/
.. _PyPI: https://pypi.org/
.. _PyScaffold's contributor's guide: https://pyscaffold.org/en/stable/contributing.html
.. _Pytest can drop you: https://docs.pytest.org/en/stable/how-to/failures.html#using-python-library-pdb-with-pytest
.. _Python Software Foundation's Code of Conduct: https://www.python.org/psf/conduct/
.. _reStructuredText: https://www.sphinx-doc.org/en/master/usage/restructuredtext/
.. _Sphinx: https://www.sphinx-doc.org/en/master/
.. _TestPyPI: https://test.pypi.org
.. _virtual environment: https://realpython.com/python-virtual-environments-a-primer/
.. _virtualenv: https://virtualenv.pypa.io/en/stable/

.. _GitHub's fork and pull request workflow: https://guides.github.com/activities/forking/
.. _GitHub web interface: https://docs.github.com/en/repositories/working-with-files/managing-files/editing-files
.. _GitHub's code editor: https://docs.github.com/en/repositories/working-with-files/managing-files/editing-files
