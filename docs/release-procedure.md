# Release procedure

1. Update version number to, e.g. 0.5.0, in `pyunfold/__version__.py`

2. Make sure `docs/changelog.rst` has up to date release notes

3. Commit change and push to `master`

       git add <changed-files>
       git commit -m "Version 0.5.0 release"
       git push origin master

4. Run `python setup.py upload`. This will

    - Build source and wheel (universal) distributions
    - Upload PyCondor package to PyPI via twine
    - Push git tag for release to GitHub

5. Update version to, e.g. `0.6.0.dev0`, in `pyunfold/__version__.py`

6. Add a new  section for the unreleased version to `docs/source/changelog.rst`

7. Commit change and push to `master`

        git add <changed-files>
        git commit -m "Bump version to 0.6.0.dev0 release"
        git push origin master

8. Update version and hash in the [recipe at conda-forge/pyunfold-feedstock](https://github.com/conda-forge/pyunfold-feedstock/blob/master/recipe/meta.yaml) & submit a pull request.
    
    Note: there is a conda-forge bot that may take care of opening this PR automatically.
